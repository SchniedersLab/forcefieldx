package ffx.potential.extended;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtUtils.SB;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;

import static ffx.potential.extended.ExtUtils.DebugHandler.buglog;
import static ffx.potential.extended.ExtUtils.prop;

/**
 *
 * @author slucore
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());
    
    // Properties
    private static final boolean esvBiasTerm = prop("esv-biasTerm", true);
    private static final boolean esvTempOverride = prop("esv-tempOverride", true);
    private static final boolean esvInterpLambda = prop("esv-interpLambda", false);
    private Double propagationTemperature = ExtConstants.roomTemperature;
    
    // ESV variables
    private int numESVs;
    private int nAtomsExt;
    private List<ExtendedVariable> esvList;
    private ExtendedVariable[] esvByAtom;
    private boolean esvTerm, phTerm, vdwTerm, mpoleTerm;
    private Atom[] atomsExt;
    private Atom[] atomsExtVdw;
    private int[] moleculeExt;     // Corresponds to extendedAtoms[]
    private List<Atom> atomListOne;
    private List<Atom> atomListZro;
    
    // Potential Objects
    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final VanDerWaals vdw;
    private final ParticleMeshEwaldQI pme;  // must not be global frame
//    private final NeighborList extendedNeighborList;
    
    // TODO set default true
    private final boolean esvUsePme = prop(new String[]{"esv-usePME", "esv-usePme", "esv-usepme"}, false);
    private final boolean esvUseVdw = prop(new String[]{"esv-useVDW", "esv-useVdw", "esv-useVdW", "esv-usevdw"}, true);
    private double latestDiscrBias, latestPhBias;

    public ExtendedSystem(MolecularAssembly molecularAssembly) {
        if (molecularAssembly == null) {
            throw new IllegalArgumentException();
        }
        mola = molecularAssembly;
        if (mola.getAtomArray() != null) {
            nAtomsExt = mola.getAtomArray().length;
        }
        Potential potential = mola.getPotentialEnergy();
        if (!(potential instanceof ForceFieldEnergy)) {
            logger.severe("ExtendedSystem supported only for ForceFieldEnergy potentials.");
        }
        ffe = (ForceFieldEnergy) potential;
        
        ForceField ff = mola.getForceField();
        vdwTerm = ff.getBoolean(ForceField.ForceFieldBoolean.VDWTERM, true);
        mpoleTerm = ff.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, false);   // TODO set default true

        vdw = (vdwTerm && esvUseVdw) ? ffe.getVdwNode() : null;
//        extendedNeighborList = (vdw != null) ? vdw.getExtendedNeighborList() : null;
        
        if (mpoleTerm && esvUsePme) {
            ParticleMeshEwald pmeNode = ffe.getPmeNode();
            if (pmeNode instanceof ParticleMeshEwaldQI) {
                pme = (ParticleMeshEwaldQI) pmeNode;
            } else {
                logger.severe("Extended system supported only for quasi-internal ParticleMeshEwald.");
                pme = null;
            }
        } else {
            pme = null;
        }

        esvList = new ArrayList<>();
        atomListOne = new ArrayList<>();
        atomListZro = new ArrayList<>();
        atomListOne.addAll(Arrays.asList(mola.getAtomArray()));
        atomsExt = mola.getAtomArray();
        moleculeExt = mola.getMoleculeNumbers();
    }
    
    public boolean hasAtom (int i) {
        return esvByAtom[i] != null;
    }
    
    public ExtendedVariable atomEsv(int i) {
        return esvByAtom[i];
    }
    
    public int atomEsvId(int i) {
        return (esvByAtom[i] != null) ? esvByAtom[i].index : -1;
    }
    
    public double atomLambda(int i) {
        return (esvByAtom[i] != null) ? esvByAtom[i].getLambda() : 1.0;
    }
        
    public Atom[] getExtendedAtomArray() {
        return atomsExt;
    }
    
    public int[] getExtendedMoleculeArray() {
        return moleculeExt;
    }
    
    public boolean useEsvInterpLambda() {
        return esvInterpLambda;
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
        if (esvUseVdw) {
            vdw.updateEsvLambda();
        }
        if (esvUsePme) {
            pme.updateEsvLambda();
        }
    }

    public int count() {
        return esvList.size();
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
            temperature = ExtConstants.roomTemperature;
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
            temperature = ExtConstants.roomTemperature;
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
                logger.config(format(" Propagating ESV[%d]: %g --> %g @ psec,temp,L,bias: %g %g %.2f %.2f %.2f",
                        esv.index, was, esv.getLambda(), currentTimePs, temperature, 
                        esv.getLambda(), esv.totalBiasEnergy(temperature)));
            }
            if (esvUseVdw) {
                vdw.updateEsvLambda();
            }
            if (esvUsePme) {
                pme.updateEsvLambda();
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
        } else {
            logger.warning("Debug?");
        }
        
        buglog("\n ExtendedSystem acquired ESV: %s\n", esv);
        atomListOne = new ArrayList<>();
        atomListOne.addAll(Arrays.asList(mola.getAtomArray()));
        atomListZro.addAll(esv.getAtomListZro());
        int nAtomsFG = atomListOne.size();
        nAtomsExt = atomListOne.size() + atomListZro.size();
        atomsExt = new Atom[nAtomsExt];
        moleculeExt = new int[nAtomsExt];
        
        SB.logfn("atoms,extAtoms,total: %d %d %d", atomListOne.size(), atomListZro.size(), nAtomsExt);
        for (int i = 0; i < nAtomsExt; i++) {
            if (i < nAtomsFG) {
                atomsExt[i] = atomListOne.get(i);
            } else {
                if (i == nAtomsFG) {
                    SB.logfn(" -- AtomsZro -- ");
                }
                atomsExt[i] = atomListZro.get(i - nAtomsFG);
            }
            moleculeExt[i] = atomsExt[i].getMoleculeNumber();
            SB.logfn(" %s", atomsExt[i]);
        }
//        SB.print();
        SB.clear();
        
        if (esvByAtom == null || esvByAtom.length < nAtomsExt) {
            esvByAtom = new ExtendedVariable[nAtomsExt];
        }
        for (int i = 0; i < nAtomsExt; i++) {
            final Atom ai = atomsExt[i];
            esvByAtom[i] = ai.getESV();
        }
        if (vdwTerm) {
            vdw.updateEsvLambda();
        }
        if (mpoleTerm) {
            pme.updateEsvLambda();
        }
    }

    /**
     * @return [numESVs] gradient w.r.t. each ESV
     */
    public double[] getdEdL(boolean lambdaBondedTerms, Double temperature) {
        if (esvTempOverride) {
            temperature = ExtConstants.roomTemperature;
        } else if (temperature == null) {
            temperature = propagationTemperature;
        }
        SB.logfn(" ESV derivative components: ");
        
        double[] vdwDeriv = (esvUseVdw) ? vdw.getdEdEsv() : null;
        double[] pmeDeriv = (esvUsePme) ? pme.getdEdEsv() : null;
        double esvDeriv[] = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            esvDeriv[i] = 0.0;
            if (esvBiasTerm) {
                esvDeriv[i] += esvList.get(i).totalBiasDeriv(temperature);
            }
            if (esvUseVdw) {
                SB.logfn("  vdW   %d: %g\n", i, vdwDeriv[i]);
                esvDeriv[i] += vdwDeriv[i];
            }
            if (esvUsePme) {
                SB.logfn("  PME   %d: %g\n", i, pmeDeriv[i]);
                esvDeriv[i] += pmeDeriv[i];
            }
        }
        SB.print();
        return esvDeriv;
    }

    public String getDecomposition() {
        return  format("    %16s %16.8f\n", "Discretization", latestDiscrBias) +
                format("    %16s %16.8f\n", "pH Electrostat", latestPhBias);
    }
    
    /**
     * @deprecated Until second derivatives w.r.t. esvLambda become tractable.
     */
    @Deprecated
    public double[][] getdEdXdL() {
        List<double[][]> terms = new ArrayList<>();
        if (vdwTerm) {
            //terms.add(vdw.getdEdXdEsv());
        }
        if (mpoleTerm) {
            //terms.add(particleMeshEwald.getdEdXdEsv());
        }
        if (terms.isEmpty()) {
            return new double[numESVs][3 * nAtomsExt];
        }
        return eleSum2DArrays(terms, numESVs, 3 * nAtomsExt);
    }

    /**
     * @deprecated Until second derivatives w.r.t. esvLambda become tractable.
     */
    @Deprecated
    public double[] getd2EdL2(boolean lambdaBondedTerms) {
        List<double[]> terms = new ArrayList<>();
        double[] bias = new double[numESVs];
        for (ExtendedVariable esv : esvList) {
            bias[esv.index] = esvList.get(esv.index).discretizationBiasEnergy();
        }
        terms.add(bias);
        if (!lambdaBondedTerms) {
            if (vdwTerm) {
                //terms.add(vdw.getd2EdEsv2());
            }
            if (mpoleTerm) {
                //terms.add(particleMeshEwald.getd2EdEsv2());
            }
        }
        if (terms.isEmpty()) {
            return new double[numESVs];
        }
        return eleSum1DArrays(terms, numESVs);
    }

    public String getLambdaList() {
        StringBuilder sb = new StringBuilder();
        sb.append("L:");
        for (int i = 0; i < numESVs; i++) {
            sb.append(format(" %6.4f", esvList.get(i).getLambda()));
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
