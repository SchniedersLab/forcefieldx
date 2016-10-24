package ffx.potential.extended;

import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;
import java.util.logging.Logger;
import java.util.stream.Stream;

import static java.lang.String.format;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwaldCart;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;

/**
 *
 * @author slucore
 */
public class ExtendedSystem {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());

    // ESV variables
    private int numESVs;
    private int nAtoms;
    private List<ExtendedVariable> esvList;
    private StringBuilder esvLogger;
    private final boolean ESV_DEBUG = false;
    private boolean esvTerm, phTerm, vdwTerm, mpoleTerm;
    // Potential Objects
    private final MolecularAssembly mola;
    private final ForceField ff;
    private final VanDerWaals vdw;
    private final ParticleMeshEwaldQI pme;  // must not be global frame
    private boolean usePME = false;
    private boolean useVDW = true;
    private double[] esvGrad;
    // Application-specific variables
    private final double pHconst;
    private double currentTemperature = 298.15;

    public ExtendedSystem() {
        logger.warning("Invoked ExtendedSystem null/temporary constuctor.");
        mola = null;
        nAtoms = 0;
        ff = null;
        vdw = null;
        pme = null;
        pHconst = 7.4;
        esvTerm = false;
        phTerm = false;
        vdwTerm = false;
        mpoleTerm = false;
        esvList = new ArrayList<>();
    }

    public ExtendedSystem(MolecularAssembly mola) {
        this(mola, 7.4);
    }

    public ExtendedSystem(MolecularAssembly mola, double pH) {
        this.mola = mola;
        if (mola == null) {
            logger.warning("Null mola in ESV...");
            throw new IllegalArgumentException();
        } else {
            logger.info("What mola looks like to ES: " + mola.toString());
        }
        if (mola.getAtomArray() != null) {
            nAtoms = mola.getAtomArray().length;
        }
        Potential potential = mola.getPotentialEnergy();
        if (!(potential instanceof ForceFieldEnergy)) {
            logger.warning("ExtendedSystem currently supported only for ForceFieldEnergy potentials.");
        }
        this.ff = mola.getForceField();
        
        esvTerm = ff.getBoolean(ForceField.ForceFieldBoolean.ESVTERM, false);
        phTerm = ff.getBoolean(ForceField.ForceFieldBoolean.PHTERM, false);
        vdwTerm = ff.getBoolean(ForceField.ForceFieldBoolean.VDWTERM, true);
        mpoleTerm = ff.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, false);
        
        if (!esvTerm) {
            logger.severe("Extended system created while !esvTerm.");
        }
        
        ForceFieldEnergy ffe = null;
        if (potential instanceof ForceFieldEnergy) {
            ffe = (ForceFieldEnergy) potential;
        } else {
            logger.warning("Extended system only suppoted by force field potentials.");
        }

        if (!vdwTerm) {
            vdw = null;
        } else {
            vdw = ffe.getVdwNode();
            if (vdw == null) {
                logger.warning("Extended system found null vanderWaals object.");
            }
        }
        
        if (!mpoleTerm) {
            pme = null;
        } else {
            ParticleMeshEwald pmeNode = ffe.getPmeNode();
            if (pmeNode instanceof ParticleMeshEwaldQI) {
                this.pme = (ParticleMeshEwaldQI) pmeNode;
            } else if (pmeNode instanceof ParticleMeshEwaldCart) {
                logger.warning("Extended system cannot operate with global-frame ParticleMeshEwald.");
                this.pme = null;
            } else {
                logger.warning(format("Extended system constructed with null or invalid PME: %s", pmeNode.toString()));
                this.pme = null;
            }
        }

        this.pHconst = pH;
        esvList = new ArrayList<>();
    }

    /**
     * Read-only; modifying the returned list has no effect. TODO LOW modify
     * remaining calls to this function to instead use the improved and
     * naturally parallelizable stream().
     */
    public List<ExtendedVariable> getESVList() {
        return esvList;
    }
    
    public void setTemperature(double temp) {
        currentTemperature = temp;
    }

    /**
     * Stream the extended variable list.
     */
    public Stream<ExtendedVariable> stream() {
        return esvList.stream();
    }

    public int num() {
        if (numESVs != esvList.size()) {
            logger.severe("Programming error: ExtendedSystem.num()");
        }
        return numESVs;
    }

    /**
     * Get the zero/unity bias and the pH energy term. PME and vdW not included
     * here as their ESV dependence isn't decomposable.
     */
    public double biases() {
        if (!esvTerm) {
            logger.warning("Called ExtendedSystem energy while !esvTerm.");
        }
        if (esvList == null || esvList.isEmpty()) {
            logger.warning("Called for extended energy with null/empty esvList.");
        }
        double biasEnergySum = 0.0;
        double phEnergySum = 0.0;

        StringBuilder sb = new StringBuilder();
        for (ExtendedVariable esv : esvList) {
            sb.append(format("%.2f ", esv.getLambda()));
            biasEnergySum += esv.getBiasEnergy();
            if (phTerm && esv instanceof TitrationESV) {
                phEnergySum += ((TitrationESV) esv).getPhEnergy(pHconst, currentTemperature);
            }
        }

        esvLogger.append(format(" [ldh: %s, bias: %g, pH: %g] ", sb.toString(), biasEnergySum, phEnergySum));
        logger.info(esvLogger.toString());
        esvLogger = new StringBuilder();

        return biasEnergySum + phEnergySum;
    }

    /**
     * Update the position of all ESV particles.
     */
    public void propagateESVs(double temperature, double dt, Double currentTimePs) {
        currentTemperature = temperature;
        double[] dEdLdh = getdEdLdh(false);
        if (esvList != null && !esvList.isEmpty()) {
            for (ExtendedVariable esv : esvList) {
                double was = esv.getLambda();
                esv.propagateLamedh(dEdLdh[esv.index], temperature, dt);
                if (currentTimePs != null) {
                    logger.info(format(" Propagating ESV[%d]: %g --> %g @ temp,psec: %g %g",
                            esv.index, was, esv.getLambda(), temperature, currentTimePs));
                } else {
                    logger.info(format(" Propagating ESV[%d]: %g --> %g @ temp: %g",
                            esv.index, was, esv.getLambda(), temperature));
                }
            }
            if (usePME) {
                pme.sourceLamedh();
            }
        }
    }

    /**
     * Returns *the difference* between full bonded energy and lamedh-scaled
     * bonded energy. Lamedh_Bonded = esvBondedEnergy() = FFE.totalBondedEnergy
     * - esvBondedCorrection();
     */
    public double bonded(boolean[] termFlags, BondedTerm[][] termArrays,
            boolean gradient, boolean lambdaBondedTerms) {
        return calcBondedTerms(true, termFlags, termArrays, gradient, lambdaBondedTerms);
    }

    /**
     * Returns the lamedh-scaled bonded energy. Lamedh_Bonded =
     * esvBondedEnergy() = FFE.totalBondedEnergy - esvBondedCorrection();
     */
    public double totalBonded(boolean[] termFlags, BondedTerm[][] termArrays,
            boolean gradient, boolean lambdaBondedTerms) {
        return calcBondedTerms(false, termFlags, termArrays, gradient, lambdaBondedTerms);
    }

    /**
     * FFE passes in: List<Object[]> termArrays = new ArrayList<>(Arrays.asList(
     * bonds, angles, stretchBends, ureyBradleys, outOfPlaneBends, torsions,
     * piOrbitalTorsions, torsionTorsions, improperTorsions, restraintBonds));
     */
    private double calcBondedTerms(boolean correctionOnly,
            boolean[] termFlags, BondedTerm[][] termArrays,
            boolean gradient, boolean lambdaBondedTerms) {
        double bondEnergy = 0.0, angleEnergy = 0.0, stretchBendEnergy = 0.0,
                ureyBradleyEnergy = 0.0, outOfPlaneBendEnergy = 0.0, torsionEnergy = 0.0,
                piOrbitalTorsionEnergy = 0.0, torsionTorsionEnergy = 0.0,
                improperTorsionEnergy = 0.0, restraintBondEnergy = 0.0;

        if (termFlags[0]) {
            Bond[] bonds = (Bond[]) termArrays[0];
            for (int i = 0; i < bonds.length; i++) {
                Bond b = bonds[i];
                if (lambdaBondedTerms && !b.applyLambda()) {
                    continue;
                }
                double be = b.energy(gradient);
                bondEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(b)) * be
                        : (be * lamedhScaling(b));
            }
        }

        if (termFlags[1]) {
            Angle[] angles = (Angle[]) termArrays[1];
            for (int i = 0; i < angles.length; i++) {
                Angle a = angles[i];
                if (lambdaBondedTerms && !a.applyLambda()) {
                    continue;
                }
                double ae = a.energy(gradient);
                angleEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(a)) * ae
                        : (ae * lamedhScaling(a));
            }
        }

        if (termFlags[2]) {
            StretchBend[] stretchBends = (StretchBend[]) termArrays[2];
            for (int i = 0; i < stretchBends.length; i++) {
                StretchBend stretchBend = stretchBends[i];
                if (lambdaBondedTerms && !stretchBend.applyLambda()) {
                    continue;
                }
                double sbe = stretchBend.energy(gradient);
                stretchBendEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(stretchBend)) * sbe
                        : (sbe * lamedhScaling(stretchBend));
            }
        }

        if (termFlags[3]) {
            UreyBradley[] ureyBradleys = (UreyBradley[]) termArrays[3];
            for (int i = 0; i < ureyBradleys.length; i++) {
                UreyBradley ureyBradley = ureyBradleys[i];
                if (lambdaBondedTerms && !ureyBradley.applyLambda()) {
                    continue;
                }
                double ube = ureyBradley.energy(gradient);
                ureyBradleyEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(ureyBradley)) * ube
                        : (ube * lamedhScaling(ureyBradley));
            }
        }

        if (termFlags[4]) {
            OutOfPlaneBend[] outOfPlaneBends = (OutOfPlaneBend[]) termArrays[4];
            for (int i = 0; i < outOfPlaneBends.length; i++) {
                OutOfPlaneBend outOfPlaneBend = outOfPlaneBends[i];
                if (lambdaBondedTerms && !outOfPlaneBend.applyLambda()) {
                    continue;
                }
                double oope = outOfPlaneBend.energy(gradient);
                outOfPlaneBendEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(outOfPlaneBend)) * oope
                        : (outOfPlaneBend.energy(gradient) * lamedhScaling(outOfPlaneBend));
            }
        }

        if (termFlags[5]) {
            Torsion[] torsions = (Torsion[]) termArrays[5];
            for (int i = 0; i < torsions.length; i++) {
                Torsion torsion = torsions[i];
                if (lambdaBondedTerms && !torsion.applyLambda()) {
                    continue;
                }
                double te = torsion.energy(gradient);
                torsionEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(torsion)) * te
                        : (torsion.energy(gradient) * lamedhScaling(torsion));
            }
        }

        if (termFlags[6]) {
            PiOrbitalTorsion[] piOrbitalTorsions = (PiOrbitalTorsion[]) termArrays[6];
            for (int i = 0; i < piOrbitalTorsions.length; i++) {
                PiOrbitalTorsion piOrbitalTorsion = piOrbitalTorsions[i];
                if (lambdaBondedTerms && !piOrbitalTorsion.applyLambda()) {
                    continue;
                }
                double pote = piOrbitalTorsion.energy(gradient);
                piOrbitalTorsionEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(piOrbitalTorsion)) * pote
                        : (piOrbitalTorsion.energy(gradient) * lamedhScaling(piOrbitalTorsion));
            }
        }

        if (termFlags[7]) {
            TorsionTorsion[] torsionTorsions = (TorsionTorsion[]) termArrays[7];
            for (int i = 0; i < torsionTorsions.length; i++) {
                TorsionTorsion torsionTorsion = torsionTorsions[i];
                if (lambdaBondedTerms && !torsionTorsion.applyLambda()) {
                    continue;
                }
                double tte = torsionTorsion.energy(gradient);
                torsionTorsionEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(torsionTorsion)) * tte
                        : (torsionTorsion.energy(gradient) * lamedhScaling(torsionTorsion));
            }
        }

        if (termFlags[8]) {
            ImproperTorsion[] improperTorsions = (ImproperTorsion[]) termArrays[8];
            for (int i = 0; i < improperTorsions.length; i++) {
                ImproperTorsion improperTorsion = improperTorsions[i];
                if (lambdaBondedTerms && !improperTorsion.applyLambda()) {
                    continue;
                }
                double ite = improperTorsion.energy(gradient);
                improperTorsionEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(improperTorsion)) * ite
                        : (ite * lamedhScaling(improperTorsion));
            }
        }

        if (termFlags[9]) {
            RestraintBond[] restraintBonds = (RestraintBond[]) termArrays[9];
            for (int i = 0; i < restraintBonds.length; i++) {
                RestraintBond rb = restraintBonds[i];
                if (lambdaBondedTerms && !rb.applyLambda()) {
                    continue;
                }
                double rbe = rb.energy(gradient);
                restraintBondEnergy += (correctionOnly)
                        ? (1.0 - lamedhScaling(rb)) * rbe
                        : (rbe * lamedhScaling(rb));
            }
        }

        double esvBonded = bondEnergy + angleEnergy + stretchBendEnergy
                + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy
                + piOrbitalTorsionEnergy + torsionTorsionEnergy
                + improperTorsionEnergy + restraintBondEnergy;
        if (correctionOnly) {
            esvBonded = -esvBonded;     // Such that E_corrected = E_normal `+` esvBonded.
        }
        return esvBonded;
    }

    /**
     * Get amount by which this term should be scaled d/t any ESVs.
     *
     * @param rols
     * @return
     */
    private double lamedhScaling(ROLS rols) {
        double totalScale = 1.0;
        for (ExtendedVariable esv : esvList) {
            OptionalDouble scale = esv.getROLSScaling(rols);
            if (scale.isPresent()) {
                totalScale *= esv.getROLSScaling(rols).getAsDouble();
                if (ESV_DEBUG) {
                    esvLogger.append(String.format(" Scaling by ESV[%d] (%4.2f) : ROLS %s (%4.2f)\n",
                            esv.index, scale.getAsDouble(), rols.toString(), totalScale));
                }
            }
        }
        return totalScale;
    }

    /**
     *
     * @param esv
     */
    public void addVariable(ExtendedVariable esv) {
        if (esvList == null) {
            esvList = new ArrayList<>();
        }
        esvList.add(esv);
        numESVs = esvList.size();
        esvTerm = true;
        esvGrad = new double[numESVs];
        
        if (esv instanceof TitrationESV) {
            phTerm = true;
        }

        logger.info(String.format(" ExtendedSystem acquired ESV: %s\n", esv));
        if (esvLogger == null) {
            esvLogger = new StringBuilder(String.format(" ESV Scaling: \n"));
        }
    }

    /**
     * @return [numESVs] gradient w.r.t. each lamedh
     */
    public double[] getdEdLdh(boolean lambdaBondedTerms) {
        StringBuilder sb = new StringBuilder(format(" ESV derivative components: \n"));
        List<double[]> terms = new ArrayList<>();
        double[] vdwGrad = new double[numESVs];
        double[] pmeGrad = new double[numESVs];
        double[] biasGrad = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            biasGrad[i] = esvList.get(i).getdBiasdLdh();
        }
        terms.add(biasGrad);
        if (!lambdaBondedTerms) {
            if (useVDW) {
                vdw.getdEdLdh(esvGrad);
            }
            if (usePME) {
//                pme.getdEdLdh(esvGrad);
            }
//            if (restraintBondTerm) {
//                for (int i = 0; i < nRestraintBonds; i++) {
//                    // TODO terms.add(restraintBonds[i].getdEdLdh());
//                }
//            }
//            if (ncsTerm && ncsRestraint != null) {
//                // TODO terms.add(ncsRestraint.getdEdLdh());
//            }
//            if (restrainTerm && !coordRestraints.isEmpty()) {
//                for (CoordRestraint restraint : coordRestraints) {
//                    // TODO terms.add(restraint.getdEdLdh());
//                }
//            }
//            if (comTerm && comRestraint != null) {
//                // TODO terms.add(comRestraint.getdEdLdh());
//            }
        }
        for (int i = 0; i < numESVs; i++) {
            sb.append(format("  Bias %d: %g\n", i, biasGrad[i]));
            if (useVDW) {
                sb.append(format("  vdW  %d: %g\n", i, vdwGrad[i]));
            }
            if (usePME) {
                sb.append(format("  PME  %d: %g\n", i, pmeGrad[i]));
            }
        }
        logger.config(sb.toString());
        if (terms.isEmpty()) {
            return new double[numESVs]; // zeroes
        }
        return esvGrad;
    }

    /**
     * @return [numESVs][nAtoms]
     */
    public double[][] getdEdXdLdh() {
        List<double[][]> terms = new ArrayList<>();
        if (vdwTerm) {
            //terms.add(vdw.getdEdXdLdh());
        }
        if (mpoleTerm) {
            // TODO terms.add(particleMeshEwald.getdEdXdLdh());
        }
//        if (restraintBondTerm) {
//            for (int i = 0; i < nRestraintBonds; i++) {
//                // TODO terms.add(restraintBonds[i].getdEdXdLdh());
//            }
//        }
//        if (ncsTerm && ncsRestraint != null) {
//            // TODO terms.add(ncsRestraint.getdEdXdLdh());
//        }
//        if (restrainTerm && !coordRestraints.isEmpty()) {
//            for (CoordRestraint restraint : coordRestraints) {
//                // TODO terms.add(restraint.getdEdXdLdh());
//            }
//        }
//        if (comTerm && comRestraint != null) {
//            // TODO terms.add(comRestraint.getdEdXdLdh());
//        }
        if (terms.isEmpty()) {
            return new double[numESVs][3 * nAtoms];
        }
        return eleSum2DArrays(terms, numESVs, 3 * nAtoms);
    }

    public double[] getd2EdLdh2(boolean lambdaBondedTerms) {
        List<double[]> terms = new ArrayList<>();
        double[] bias = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            bias[i] = esvList.get(i).getBiasEnergy();
        }
        terms.add(bias);
        if (!lambdaBondedTerms) {
            if (vdwTerm) {
                // terms.add(vdw.getd2EdLdh2());
            }
            if (mpoleTerm) {
//                 TODO terms.add(particleMeshEwald.getd2EdLdh2());
            }
//            if (restraintBondTerm) {
//                for (int i = 0; i < nRestraintBonds; i++) {
//                    // TODO terms.add(restraintBonds[i].getd2EdLdh2());
//                }
//            }
//            if (ncsTerm && ncsRestraint != null) {
//                // TODO terms.add(ncsRestraint.getd2EdLdh2());
//            }
//            if (restrainTerm && !coordRestraints.isEmpty()) {
//                for (CoordRestraint restraint : coordRestraints) {
//                    // TODO terms.add(restraint.getd2EdLdh2());
//                }
//            }
//            if (comTerm && comRestraint != null) {
//                // TODO terms.add(comRestraint.getd2EdLdh2());
//            }
        }
        if (terms.isEmpty()) {
            return new double[numESVs]; // zeroes
        }
        return eleSum1DArrays(terms, numESVs);
    }

    /**
     * Element-wise sum over a list of 1D double arrays. Benchmarks faster than
     * Java8's stream API.
     */
    private static double[] eleSum1DArrays(List<double[]> terms, int numESVs) {
        double[] termSum = new double[terms.size()];
        for (int iTerm = 0; iTerm < terms.size(); iTerm++) {
            double[] currentTerm = terms.get(iTerm);
            if (currentTerm.length != numESVs) {
                logger.warning(format("iTerm %d length: %d, numESVs: %d", iTerm, terms.get(iTerm).length, numESVs));
                throw new IndexOutOfBoundsException();
            }
            for (int iESV = 0; iESV < numESVs; iESV++) {
                termSum[iESV] += currentTerm[iESV];
            }
        }
        return termSum;
    }

    /**
     * Element-wise sum over a list of 2D double arrays.
     */
    private static double[][] eleSum2DArrays(List<double[][]> terms, int numESVs, int nVars) {
        if (terms == null || terms.isEmpty()) {
            throw new NullPointerException("Summing an empty or null terms list.");
        }
        double[][] termSum = new double[numESVs][nVars];
        for (int iTerm = 0; iTerm < terms.size(); iTerm++) {
            double[][] currentTerm = terms.get(iTerm);
            if (currentTerm.length != numESVs) {
                throw new IndexOutOfBoundsException();
            }
            for (int iESV = 0; iESV < numESVs; iESV++) {
                if (currentTerm[iESV].length != nVars) {
                    throw new IndexOutOfBoundsException();
                }
                for (int iAtom = 0; iAtom < nVars; iAtom++) {
                    termSum[iESV][iAtom] += currentTerm[iESV][iAtom];
                }
            }
        }
        return termSum;
    }
}
