package ffx.potential.extended;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Stream;

import static java.lang.String.format;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
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
    private final ForceFieldEnergy ffe;
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
        ffe = null;
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
//            logger.warning("Called ExtendedSystem energy while !esvTerm.");
            return 0.0;
        }
        if (esvList == null || esvList.isEmpty()) {
//            logger.warning("Called for extended energy with null/empty esvList.");
            return 0.0;
        }
        double biasEnergySum = 0.0;
        double phEnergySum = 0.0;

        StringBuilder sb = new StringBuilder();
        for (ExtendedVariable esv : esvList) {
            sb.append(format("%.2f ", esv.getLambda()));
            biasEnergySum += esv.getDiscretizationBiasEnergy();
            if (phTerm && esv instanceof TitrationESV) {
                phEnergySum += ((TitrationESV) esv).getPhBiasEnergy(currentTemperature);
            }
        }

//        logger.info(esvLogger.toString());
//        esvLogger = new StringBuilder();

        return biasEnergySum + phEnergySum;
    }

    /**
     * Update the position of all ESV particles.
     */
    public void propagateESVs(double temperature, double dt, double currentTimePs) {
        currentTemperature = temperature;
        double[] dEdLdh = getdEdL(false);
        if (esvList != null && !esvList.isEmpty()) {
            for (ExtendedVariable esv : esvList) {
                double was = esv.getLambda();
                esv.propagate(dEdLdh[esv.index], temperature, dt);
                if (esv instanceof TitrationESV) {
                    logger.info(format(" Propagating ESV[%d]: %g --> %g @ psec,temp,lbd,zuBias,phBias: %g %g %.2f %.2f %.2f",
                            esv.index, was, esv.getLambda(), currentTimePs, temperature, 
                            esv.getLambda(), esv.getDiscretizationBiasEnergy(),
                            ((TitrationESV) esv).getPhBiasEnergy(temperature)));
                } else {
                    logger.info(format(" Propagating ESV[%d]: %g --> %g @ psec,temp,lbd,zuBias: %g %g %.2f %.2f",
                            esv.index, was, esv.getLambda(), currentTimePs, temperature, 
                            esv.getLambda(), esv.getDiscretizationBiasEnergy()));
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
    public double[] getdEdL(boolean lambdaBondedTerms) {
        StringBuilder sb = new StringBuilder(format(" ESV derivative components: \n"));
        List<double[]> terms = new ArrayList<>();
        double[] vdwGrad = new double[numESVs];
        double[] pmeGrad = new double[numESVs];
        double[] biasGrad = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            biasGrad[i] = esvList.get(i).getdDiscretizationBiasdL();
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
            bias[i] = esvList.get(i).getDiscretizationBiasEnergy();
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
