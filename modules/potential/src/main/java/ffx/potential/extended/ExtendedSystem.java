/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.extended;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;

import edu.rit.pj.reduction.SharedDouble;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.extended.ExtUtils.SB;
import ffx.potential.nonbonded.Multipole;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;

import static ffx.potential.extended.ExtUtils.logf;
import static ffx.potential.extended.ExtUtils.prop;

/**
 *
 * @author slucore
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());
    public static boolean esvSystemActive = false;
    private static boolean printedConfig = false;

    // Properties: Static Final
    public static final boolean esvTempOverride = prop("esv-tempOverride", true);
    public static final boolean esvUsePme = prop("esv-usePme", false);
    public static final boolean esvUseVdw = prop("esv-useVdw", true);
    public static final boolean ffePrintOverride = prop("ffe-printOverride", false);
    public static final boolean esvScaleBonded = prop("esv-scaleBonded", true);
    public static final boolean esvDecomposeBonded = prop("esv-decomposeBonded", true);
    // Properties: Mutable
    private boolean esvBiasTerm = prop("esv-biasTerm", true);

    // Atom Lists
    /**
     * Foreground/background atom lists and arrays based on the
     * "Titratable Hydrogens Only" inclusion criterion.
     * Heavy atoms of deprotonated forms are not present. Used by VdW.
     */
    private List<Atom> atomList;
    private Atom[] atomsExt;
    private int[] moleculeExt;
    private int nAtomsExt;
    private ExtendedVariable[] esvForShared;
    private ExtendedVariable[] esvForUnshared;
    private int[] fg2bgIdx;
    /**
     * Foreground/background atoms lists and arrays based on the
     * "All Atoms" inclusion criterion. Used by bonded terms.
     */
//    private Atom[] atomsExtPlusBackground;
//    private int[] moleculeExtPlusBackground;

    // ESV variables
    private int numESVs;
    private List<ExtendedVariable> esvList;
    private boolean esvTerm, phTerm, vdwTerm, mpoleTerm;
    private Double propagationTemperature;

    // Potential Objects
    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final VanDerWaals vdw;
    private final ParticleMeshEwaldQI pme;  // must not be global frame

    public ExtendedSystem(MolecularAssembly molecularAssembly) {
        if (molecularAssembly == null) {
            throw new IllegalArgumentException();
        }
        mola = molecularAssembly;
        Potential potential = mola.getPotentialEnergy();
        if (!(potential instanceof ForceFieldEnergy)) {
            logger.severe("ExtendedSystem supported only for ForceFieldEnergy potentials.");
        }
        ffe = (ForceFieldEnergy) potential;
        ffe.setPrintOverride(ffePrintOverride);

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
        propagationTemperature = ExtConstants.roomTemperature;

        // Initialize atom arrays with the existing assembly.
        atomsExt = mola.getAtomArray();
        nAtomsExt = atomsExt.length;
        moleculeExt = mola.getMoleculeNumbers();
        atomList = Arrays.asList(mola.getAtomArray());        
        esvForUnshared = new ExtendedVariable[nAtomsExt];
        esvForShared = new ExtendedVariable[nAtomsExt];
    }

    /**
     * Getters for all extended atom properties needed  by other cl
     */
    public Atom[] getExtendedAtoms() {
        return atomsExt;
    }
    public int[] getExtendedMolecule() {
        return moleculeExt;
    }
    public boolean isExtended(int i) {
        return isShared(i) || isUnshared(i);
    }
    public boolean isShared(int i) {
        return esvForShared[i] != null;
    }
    public boolean isUnshared(int i) {
        return esvForUnshared[i] != null;
    }
    public ExtendedVariable getEsvForAtom(int i) {
        return (isShared(i)) ? esvForShared[i]
                : (isUnshared(i)) ? esvForUnshared[i]
                : Defaults.esv;
    }    
    public int getEsvIdForAtom(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).index : Defaults.esvId;
    }
    public double getLambda(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getLambda() : Defaults.lambda;
    }
    public double getLambdaSwitch(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getLambdaSwitch() : Defaults.lambdaSwitch;
    }    
    public double getSwitchDeriv(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getSwitchDeriv() : Defaults.switchDeriv;
    }

    public void setLambda(int esvID, double lambda) {
        esvList.get(esvID).setLambda(lambda);
        updateListeners();
    }

    public void updateListeners() {
        if (esvScaleBonded) {
            updateBondedEsvLambda();
        }
        if (esvUseVdw) {
            vdw.updateEsvLambda();
        }
        if (esvUsePme) {
            pme.updateEsvLambda();
        }
    }
    
    public void initializeBackgroundMultipoles(List<Atom> atomsBackground) {
        for (int i = 0; i < atomsBackground.size(); i++) {
            Atom bg = atomsBackground.get(i);
            Multipole multipole = Multipole.buildMultipole(bg, mola.getForceField());
            if (multipole == null) {
                logger.severe(format("No multipole could be assigned to atom %s of type %s.",
                        bg.toString(), bg.getAtomType()));
            }
        }
    }

    public int n() {
        return esvList.size();
    }

    public void setEsvBiasTerm(boolean set) {
        esvBiasTerm = set;
    }

    /**
     * Get ESV biases such as discretization, pH, etc.
     * This method public and final for error-checking; new ESVs should override biasEnergy().
     */
    public final double getBiasEnergy(Double temperature) {
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
        return biasEnergy(temperature);
    }

    /**
     * Get the
     * @param temperature
     * @return
     */
    private double biasEnergy(double temperature) {
        double biasEnergySum = 0.0;
        for (ExtendedVariable esv : this) {
            biasEnergySum += esv.getTotalBias(temperature, false);
        }
        return biasEnergySum;
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
        if (esvList == null || esvList.isEmpty()) {
            return;
        }
        double[] dedl = getdEdL(temperature, false, false);
        for (ExtendedVariable esv : esvList) {
            double oldLambda = esv.getLambda();
            esv.propagate(dedl[esv.index], dt, temperature);
            double newLambda = esv.getLambda();
            logger.info(format(" Propagating ESV[%d]: %g --> %g @ psec,temp,bias: %g %g %.2f",
                    esv.index, oldLambda, newLambda,
                    currentTimePs, temperature, esv.getTotalBias(temperature, false)));
        }
        updateListeners();
    }

    /**
     *
     * @param esv
     */
    public void addVariable(ExtendedVariable esv) {
        if (!printedConfig) {
            ExtUtils.printExtConfig();
            printedConfig = true;
        }
        logf(" ExtendedSystem acquired ESV: %s", esv);
        esvSystemActive = true;
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
        
        for (int i = 0; i < nAtomsExt; i++) {
            if (esv.viewUnsharedAtoms().contains(atomsExt[i])) {
                esvForUnshared[i] = esv;
            } else {
                esvForShared[i] = esv;
            }
        }
        
        /* Background atoms don't get automatically typed by PME since they're
         * disconnected from the molecular assembly; must be done manually. */
        initializeBackgroundMultipoles(esv.viewBackgroundAtoms());
        
        fg2bgIdx = new int[nAtomsExt];
        for (int i = 0; i < nAtomsExt; i++) {
            Atom atom = atomsExt[i];
            if (atom.getEsv() != null) {
                Atom back = atom.getEsv().getBackground(atom);
                if (atom == null) {
    //                fg2bgIdx[i] = null;
                } else if (back != null) {
                    fg2bgIdx[i] = (back != null) ? back.getIndex() : -1;
                }
            }
        }
        
        updateListeners();
    }
    
    public Atom[] getExtendedAndBackgroundAtoms() {
        Atom[] extended = getExtendedAtoms();
        List<Atom> background = new ArrayList<>();
        for (ExtendedVariable esv : this) {
            background.addAll(esv.viewBackgroundAtoms());
        }
        List<Atom> mega = new ArrayList<>();
        mega.addAll(Arrays.asList(extended));
        mega.addAll(background);
        return mega.toArray(new Atom[0]);
    }
    
    public int[] getExtendedAndBackgroundMolecule() {
        Atom[] mega = getExtendedAndBackgroundAtoms();
        int[] molecule = new int[mega.length];
        for (int i = 0; i < molecule.length; i++) {
            molecule[i] = mega[i].getMoleculeNumber();
        }
        return molecule;
    }
    
    public final void crashDump() {
        SB.nlogfn("*************");
        SB.nlogfn(" Crash Dump:");
        SB.logfn("   All Atoms:");
        for (Atom atom : mola.getAtomArray()) {
            SB.logfn("     %s", atom.toString());
        }
        SB.print();
        for (ExtendedVariable esv : esvList) {
            esv.describe();
        }
        SB.nlogfn("*************");
        SB.print();
    }

    /**
     * @return [numESVs] gradient w.r.t. each ESV
     */
    public double[] getdEdL(Double temperature, boolean lambdaBondedTerms, boolean print) {
        double esvDeriv[] = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            esvDeriv[i] = getdEdL(i, temperature, lambdaBondedTerms, print);
        }
        return esvDeriv;
    }

    public double[] getdEdL() {
        return getdEdL(null, false, false);
    }

    public double getdEdL(int esvID) {
        return getdEdL(esvID, null, false, false);
    }

    public double getdEdL(int esvID, boolean print) {
        return getdEdL(esvID, null, false, print);
    }

    public double getdVdwdL(int esvID) {
        return (esvUseVdw) ? vdw.getdEdEsv(esvID) : 0.0;
    }

    public double getdEdL(int esvID, Double temperature, boolean lambdaBondedTerms, boolean print) {
        if (esvTempOverride) {
            temperature = ExtConstants.roomTemperature;
        } else if (temperature == null) {
            temperature = propagationTemperature;
        }
        ExtendedVariable esv = esvList.get(esvID);
        double esvDeriv = 0.0;
        SB.logfn(" %s derivative components: ", esv.getName());
        if (esvBiasTerm) {
            esvDeriv += esv.getTotalBiasDeriv(temperature, print);
        }
        if (esvUseVdw) {
            double dVdw = vdw.getdEdEsv(esvID);
            SB.logfn("  vdW    %d: %g", esvID, dVdw);
            esvDeriv += dVdw;
        }
        if (esvUsePme) {
            double dPme = pme.getdEdEsv(esvID);
            SB.logfn("  PME    %d: %g", esvID, dPme);
            esvDeriv += dPme;
        }
        if (esvScaleBonded) {
            double dBonded = esv.getBondedDeriv();
            esvDeriv += dBonded;
            if (print && esvDecomposeBonded) {
                // Decompose the bonded derivative into term types.
                SB.logfn("  Bonded %d: %g", esvID, dBonded);
                // Foreground portion:
                double fgSum = 0.0;
                HashMap<Class<? extends BondedTerm>,SharedDouble> fgMap =
                        esv.getBondedDerivDecomp();
                for (SharedDouble dub : fgMap.values()) {
                    fgSum += dub.get();
                }
                SB.logf("    Foreground: %9.6f" , fgSum);
                for (Class<? extends BondedTerm> clas : fgMap.keySet()) {
                    SB.logf("\n      %-16s %9.6f",
                            clas.getName().replaceAll("ffx.potential.bonded.", "") + ":",
                            fgMap.get(clas).get());
                }
                // Background portion:
                double bgSum = 0.0;
                HashMap<Class<? extends BondedTerm>,SharedDouble> bgMap =
                        esv.getBackgroundBondedDerivDecomp();
                for (SharedDouble dub : bgMap.values()) {
                    bgSum += dub.get();
                }
                SB.logf("\n    Background: %9.6f" , bgSum);
                for (Class<? extends BondedTerm> clas : bgMap.keySet()) {
                    SB.logf("\n      %-16s %9.6f",
                            clas.getName().replaceAll("ffx.potential.bonded.", "") + ":",
                            bgMap.get(clas).get());
                }
            }
        }
        SB.printIf(print);
        return esvDeriv;
    }

    public String getBiasDecomposition() {
        if (!esvBiasTerm) {
            return "";
        }
        double discrBias = 0.0;
        double phBias = 0.0;
        for (ExtendedVariable esv : esvList) {
            discrBias += esv.getDiscrBias();
            if (esv instanceof TitrationESV) {
                phBias += ((TitrationESV) esv).getPhBias(propagationTemperature);
            }
        }
        return  format("    %16s %16.8f\n", "Discretization", discrBias) +
                format("    %16s %16.8f\n", "pH Electrostat", phBias);
    }

    public String getLambdaList() {
        StringBuilder sb = new StringBuilder();
        sb.append(format("Lambdas (%d):", esvList.size()));
        for (int i = 0; i < numESVs; i++) {
            sb.append(format(" %4.2f->%4.2f",
                    esvList.get(i).getLambda(), esvList.get(i).getLambdaSwitch()));
        }
        return sb.toString();
    }

    /**
     * Parallel Java constructs in ForceFieldEnergy, VanDerWaals, and ParticleMeshEwald
     * loop over atom indices. This maps the range {0,nAtoms} to {0,numESVs} so
     * that parallelization of ESVs can be done in the same constructs.
     */
    public int parallelRangeConversion(int index) {
        return (int) Math.floor((index / (double) nAtomsExt) * numESVs);
    }

    public void updateBondedEsvLambda() {
        for (ExtendedVariable esv : this) {
            esv.updateBondedLambdas();
        }
    }
    
    public void describe() {
        SB.nlogf(" AtomsExtH: ");
        for (int i = 0; i < atomsExt.length; i++) {
            SB.nlogf(" %3d  %s", i, atomsExt[i].toString());
        }
        SB.print();
    }

    /**
     * Allows simple iteration over ESV via "for (ExtendedVariable : ExtendedSystem)".
     */
    @Override
    public Iterator<ExtendedVariable> iterator() {
        return esvList.iterator();
    }

    public double getLargestLambda() {
        double ret = 0.0;
        for (ExtendedVariable esv : this) {
            ret = (esv.getLambda() > ret) ? esv.getLambda() : ret;
        }
        return ret;
    }
    
    /**
     * These populate the order-n preloaded lambda parameter arrays in VdW and
     * PME in the absence of an attached ESV.
     */
    public static final class Defaults {
        private Defaults() {}   // value singleton
        public static final ExtendedVariable esv = null;
        public static final int esvId = -1;
        public static final double lambda = 1.0;
        public static final double lambdaSwitch = 1.0;
        public static final double switchDeriv = 0.0;
    }
}
