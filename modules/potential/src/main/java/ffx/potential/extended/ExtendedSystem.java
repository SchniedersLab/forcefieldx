package ffx.potential.extended;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwaldCart;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;

import static ffx.potential.extended.ThermoConstants.FAIL_SEVERE;
import static ffx.potential.extended.ThermoConstants.prop;

/**
 *
 * @author slucore
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());
    
    // Properties
    private static final boolean esvBiasTerm = prop("esv-biasTerm", true);
    private static final boolean esvTempOverride = prop("esv-tempOverride", false);

    // ESV variables
    private int numESVs;
    private int nAtoms;
    private List<ExtendedVariable> esvList;
    private ExtendedVariable[] esvByAtom;
    private StringBuilder esvLogger;
    private final boolean ESV_DEBUG = false;
    private boolean esvTerm, phTerm, vdwTerm, mpoleTerm;
    // Potential Objects
    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final VanDerWaals vdw;
    private final ParticleMeshEwaldQI pme;  // must not be global frame
    private Double temperature = null;
    private Double propagationTemperature = ThermoConstants.roomTemperature;
    private boolean usePME = false;
    private boolean useVDW = true;
    private double latestDiscrBias, latestPhBias;

    public ExtendedSystem(MolecularAssembly mola) {
        this.mola = mola;
        if (mola == null) {
            logger.warning("Null mola in ESV...");
            throw new IllegalArgumentException();
        }
        if (mola.getAtomArray() != null) {
            nAtoms = mola.getAtomArray().length;
        }
        Potential potential = mola.getPotentialEnergy();
        if (!(potential instanceof ForceFieldEnergy)) {
            logger.warning("ExtendedSystem currently supported only for ForceFieldEnergy potentials.");
        }
        ffe = (ForceFieldEnergy) potential;
        ForceField ff = mola.getForceField();
        
        esvTerm = ff.getBoolean(ForceField.ForceFieldBoolean.ESVTERM, false);
        phTerm = ff.getBoolean(ForceField.ForceFieldBoolean.PHTERM, false);
        vdwTerm = ff.getBoolean(ForceField.ForceFieldBoolean.VDWTERM, true);
        mpoleTerm = ff.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, false);
        
        ForceFieldEnergy ffe = null;
        if (potential instanceof ForceFieldEnergy) {
            ffe = (ForceFieldEnergy) potential;
        } else {
            logger.warning("Extended system only suppoted by force field potentials.");
        }

        if (!vdwTerm) {
            vdw = null;
        } else {
            useVDW = true;
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

        esvList = new ArrayList<>();
    }
    
    public synchronized boolean hasAtom (int i) {
        return esvByAtom[i] != null;
    }
    
    public synchronized ExtendedVariable atomEsv(int i) {
        return esvByAtom[i];
    }
    
    public synchronized int atomEsvId(int i) {
        if (esvByAtom[i] != null && esvByAtom[i].index != 0) {
            logger.severe("Failure to index.");   // TODO REMOVE
        }
        return (esvByAtom[i] != null) ? esvByAtom[i].index : -1;
    }
    
    public synchronized double atomLambda(int i) {
        return (esvByAtom[i] != null) ? esvByAtom[i].getLambda() : 1.0;
    }
    
    /**
     * For testing and debugging.
     * The ExtendedSystem notifies VdW and PME of lambda changes.
     * ExtendedVariables are responsible for updating their own bonded terms.
     */
    public void setAllLambdas(double lambda) {
        for (ExtendedVariable esv : this) {
            esv.setLambda(lambda);
        }
        vdw.updateEsvLambda();
        // pme.updateEsvLambda();
    }

    public int num() {
        if (numESVs != esvList.size()) {
            logger.severe("Programming error: ExtendedSystem.num()");
        }
        return numESVs;
    }

    /**
     * Get the zero/unity bias and the pH energy term.
     */
    public double biasEnergy(Double temperature) {
        if (!esvBiasTerm) {
            return 0.0;
        }
        if (!esvTerm) {
            logger.warning("Called ExtendedSystem energy while !esvTerm.");
            return 0.0;
        }
        if (esvList == null || esvList.isEmpty()) {
            logger.warning("Called for extended energy with null/empty esvList.");
            return 0.0;
        }
        if (esvTempOverride) {
            temperature = ThermoConstants.roomTemperature;
        } else if (temperature == null) {
//            logger.warning("Called for extended bias energy without temperature.");
            temperature = propagationTemperature;
        }
        double biasEnergySum = 0.0;
        double phEnergySum = 0.0;

        StringBuilder sb = new StringBuilder();
        for (ExtendedVariable esv : this) {
            sb.append(format("%.2f ", esv.getLambda()));
            biasEnergySum += esv.discretizationBiasEnergy();
            latestDiscrBias = biasEnergySum;
            if (phTerm && esv instanceof TitrationESV) {
                phEnergySum += ((TitrationESV) esv).phBiasEnergy(temperature);
                latestPhBias = phEnergySum;
            }
        }
        return biasEnergySum + phEnergySum;
    }

    /**
     * Update the position of all ESV particles.
     */
    public void propagateESVs(Double temperature, double dt, double currentTimePs) {
        if (esvTempOverride) {
            temperature = ThermoConstants.roomTemperature;
        } else if (temperature == null) {
            temperature = propagationTemperature;
        } else {
            propagationTemperature = temperature;
        }
        double[] dEdLdh = getdEdL(false, temperature);
        if (esvList != null && !esvList.isEmpty()) {
            for (ExtendedVariable esv : esvList) {
                double was = esv.getLambda();
                esv.propagate(dEdLdh[esv.index], dt, temperature);
                if (esv instanceof TitrationESV) {
                    logger.info(format(" Propagating ESV[%d]: %g --> %g @ psec,temp,lbd,zuBias,phBias: %g %g %.2f %.2f %.2f",
                            esv.index, was, esv.getLambda(), currentTimePs, temperature, 
                            esv.getLambda(), esv.discretizationBiasEnergy(),
                            ((TitrationESV) esv).phBiasEnergy(temperature)));
                } else {
                    logger.info(format(" Propagating ESV[%d]: %g --> %g @ psec,temp,lbd,zuBias: %g %g %.2f %.2f",
                            esv.index, was, esv.getLambda(), currentTimePs, temperature, 
                            esv.getLambda(), esv.discretizationBiasEnergy()));
                }
            }
            if (usePME) {
                pme.sourceLamedh();
            }
        }
    }
    
    /**
     *
     * @param esv
     */
    public void addVariable(ExtendedVariable esv) {
        if (esvList == null) {
            esvList = new ArrayList<>();
        }
        if (!esv.isReady()) {
            esv.readyup();
        }
        esvList.add(esv);
        numESVs = esvList.size();
        esvTerm = true;
        
        if (esv instanceof TitrationESV) {
            phTerm = true;
        }

        logger.info(String.format(" ExtendedSystem acquired ESV: %s\n", esv));
        if (esvLogger == null) {
            esvLogger = new StringBuilder(String.format(" ESV Scaling: \n"));
        }
        
        Atom[] allAtoms = mola.getAtomArray();
        nAtoms = allAtoms.length;
        if (esvByAtom == null || esvByAtom.length < nAtoms) {
            esvByAtom = new ExtendedVariable[nAtoms];
        }
        for (int i = 0; i < nAtoms; i++) {
            final Atom ai = allAtoms[i];
            for (ExtendedVariable var : this) {
                if (var.containsAtom(ai)) {
                    esvByAtom[i] = var;
                    if (i != ai.xyzIndex - 1) {
                        logger.warning(format("mismatch i,xyzIndexLessOne: %d %d", i, ai.xyzIndex - 1));
                    }
                }
            }
        }
        for (int i = 0; i < nAtoms; i++) {
            final Atom ai = allAtoms[i];
            ExtendedVariable varA = esvByAtom[i];
            ExtendedVariable varB = ai.getESV();
            if ((varA == null && varB != null) || (varA != null && varB == null)) {
                FAIL_SEVERE("What1.");
            } else if (varA != varB) {
                FAIL_SEVERE("What2.");
            }
        }
    }

    /**
     * @return [numESVs] gradient w.r.t. each lamedh
     */
    public double[] getdEdL(boolean lambdaBondedTerms, Double temperature) {
        if (esvTempOverride) {
            temperature = ThermoConstants.roomTemperature;
        } else if (temperature == null) {
//            logger.warning("Accurate pH bias requires temperature information.");
            temperature = propagationTemperature;
        }
        double[] vdwGrad = null, pmeGrad = null;
        double[] discrBiasGrad = new double[numESVs];
        double[] phBiasGrad = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            if (!esvBiasTerm) {
                discrBiasGrad[i] = 0.0;
                phBiasGrad[i] = 0.0;
                continue;
            }
            ExtendedVariable esv = esvList.get(i);
            discrBiasGrad[i] = esv.discretizationBiasDeriv();
            if (esv instanceof TitrationESV) {
                phBiasGrad[i] += ((TitrationESV) esv).phBiasDeriv(temperature);
            }
        }
//        if (!lambdaBondedTerms) {
            if (useVDW) {
                vdwGrad = vdw.getdEdEsv();
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
//        }
        StringBuilder sb = new StringBuilder(format(" ESV derivative components: \n"));
        double esvGrad[] = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            esvGrad[i] = 0.0;
            sb.append(format("  Discr %d: %g\n", i, discrBiasGrad[i]));
            esvGrad[i] += discrBiasGrad[i];
            if (esvList.get(i) instanceof TitrationESV) {
                sb.append(format("  pH    %d: %g\n", i, phBiasGrad[i]));
                esvGrad[i] += phBiasGrad[i];
            }
            if (useVDW) {
                sb.append(format("  vdW   %d: %g\n", i, vdwGrad[i]));
                esvGrad[i] += vdwGrad[i];
            }
            if (usePME) {
                sb.append(format("  PME   %d: %g\n", i, pmeGrad[i]));
                esvGrad[i] += pmeGrad[i];
            }
        }
        logger.config(sb.toString());
        return esvGrad;
    }

    public double getLatestDiscrBias() {
        return latestDiscrBias;
    }
    public double getLatestPhBias() {
        return latestPhBias;
    }
    
    /**
     * @deprecated Until second derivatives w.r.t. esvLambda become tractable.
     */
    @Deprecated
    public double[][] getdEdXdLdh() {
        List<double[][]> terms = new ArrayList<>();
        if (vdwTerm) {
            //terms.add(vdw.getdEdXdLdh());
        }
        if (mpoleTerm) {
            // TODO terms.add(particleMeshEwald.getdEdXdLdh());
        }
        if (terms.isEmpty()) {
            return new double[numESVs][3 * nAtoms];
        }
        return eleSum2DArrays(terms, numESVs, 3 * nAtoms);
    }

    /**
     * @deprecated Until second derivatives w.r.t. esvLambda become tractable.
     */
    @Deprecated
    public double[] getd2EdLdh2(boolean lambdaBondedTerms) {
        List<double[]> terms = new ArrayList<>();
        double[] bias = new double[numESVs];
        for (ExtendedVariable esv : esvList) {
            bias[esv.index] = esvList.get(esv.index).discretizationBiasEnergy();
        }
        terms.add(bias);
        if (!lambdaBondedTerms) {
            if (vdwTerm) {
                // TODO terms.add(vdw.getd2EdLdh2());
            }
            if (mpoleTerm) {
                // TODO terms.add(particleMeshEwald.getd2EdLdh2());
            }
        }
        if (terms.isEmpty()) {
            return new double[numESVs];
        }
        return eleSum1DArrays(terms, numESVs);
    }

    public String getLambdaList() {
        StringBuilder sb = new StringBuilder();
        sb.append("L: ");
        for (int i = 0; i < numESVs; i++) {
            sb.append(format(" %5.2f", esvList.get(i).getLambda()));
        }
        return sb.toString();
    }
    
    /**
     * Element-wise sum over a list of 1D double arrays. 
     * This implementation benchmarks faster than the equivalent Java8 stream() API.
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

    /**
     * Allows simple iteration over ESV via "for (ExtendedVariable : ExtendedSystem)".
     */
    @Override
    public Iterator<ExtendedVariable> iterator() {
        return esvList.iterator();
    }
}
