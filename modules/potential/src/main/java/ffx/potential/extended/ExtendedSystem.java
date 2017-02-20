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
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;

import static ffx.potential.extended.ExtUtils.prop;

/**
 *
 * @author slucore
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());
    public static boolean esvSystemActive = false;
    private static boolean printedConfig = false;

    private final EsvConfiguration config = new EsvConfiguration();
    public class EsvConfiguration {
        // Properties: Static Final
        public final boolean esvVdwFlag         = prop("esv-useVdw",true);
        public final boolean esvPmeFlag         = prop("esv-usePme", true);
        public final boolean esvTempOverride    = prop("esv-tempOverride", true);
        public final boolean esvPrintOverride   = prop("esv-printOverride", false);
        public final boolean esvScaleBonded     = prop("esv-scaleBonded", true);
        public final boolean backgroundBondedHookup = prop("esv-backgroundBonded", true);
        public final boolean esvScaleUnshared   = prop("esv-scaleUnshared", true);
        public final boolean esvDecomposeBonded = prop("esv-decomposeBonded", true);
        public final boolean esvPropagation     = prop("esv-propagation", false);
        public final Double biasOverride        = prop("esv-biasOverride", Double.NaN);
        public final double thetaMass           = prop("esv-thetaMass", 1.0e-18);            // from OSRW, reasonably 100 a.m.u.
        public final double thetaFriction       = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps

        // Properties: Mutable
        public boolean esvBiasTerm = prop("esv-biasTerm", true);
    }
    
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
    private Double lastKnownTemperature;

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
        ffe.setPrintOverride(config.esvPrintOverride);

        ForceField ff = mola.getForceField();
        vdwTerm = ff.getBoolean(ForceField.ForceFieldBoolean.VDWTERM, true);
        mpoleTerm = ff.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, false);   // TODO set default true

        vdw = (vdwTerm && config.esvVdwFlag) ? ffe.getVdwNode() : null;
        if (config.esvVdwFlag && !vdwTerm) {
            logger.severe("Error: esvVdw without vdwTerm.");
        }
//        extendedNeighborList = (vdw != null) ? vdw.getExtendedNeighborList() : null;

        if (config.esvPmeFlag) {
            if (!mpoleTerm) {
                logger.severe("Error: esvPme without mpoleTerm.");
            }
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
        lastKnownTemperature = ExtConstants.roomTemperature;

        // Initialize atom arrays with the existing assembly.
        atomsExt = mola.getAtomArray();
        nAtomsExt = atomsExt.length;
        moleculeExt = mola.getMoleculeNumbers();
        atomList = Arrays.asList(mola.getAtomArray());        
        esvForUnshared = new ExtendedVariable[nAtomsExt];
        esvForShared = new ExtendedVariable[nAtomsExt];
    }
    
    public EsvConfiguration getConfig() {
        return config;
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

    public void setLambda(int esvId, double lambda) {
        esvList.get(esvId).setLambda(lambda);
        updateListeners();
    }
    
    public void setLambda(String esvIdStr, double lambda) {
        int esvId = esvIdStr.toUpperCase().charAt(0) - 'A';
        setLambda(esvId, lambda);
    }
    
    public void setLambda(char esvIdChar, double lambda) {
        setLambda(String.valueOf(esvIdChar), lambda);
    }
    
    /**
     * DEBUG: to be called after rotating multipoles in PME.
     */
    public void updateMultipoles() {
        for (ExtendedVariable esv : this) {
            esv.updateMultipoleTypes();
        }
    }

    protected void updateListeners() {
        if (config.esvVdwFlag) {
            vdw.updateEsvLambda();
        }
        if (config.esvPmeFlag) {
            pme.updateEsvLambda();
        }
    }

    public int size() {
        return esvList.size();
    }

    protected void setEsvBiasTerm(boolean set) {
        config.esvBiasTerm = set;
    }
    
    public final double getBiasEnergy() {
        return getBiasEnergy(lastKnownTemperature);
    }

    /**
     * Get ESV biases such as discretization, pH, etc.
     * This method public and final for error-checking; new ESVs should override biasEnergy().
     */
    public final double getBiasEnergy(double temperature) {
        if (!config.esvBiasTerm) {
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
        if (config.esvTempOverride) {
            return biasEnergy(ExtConstants.roomTemperature);
        } else {
            return biasEnergy(temperature);
        }
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

    public void propagateESVs(double dt, double currentTimePs) {
        propagateESVs(lastKnownTemperature, dt, currentTimePs);
    }
    
    /**
     * Update the position of all ESV particles via langevin dynamics.
     * Temperature controls propagation speed and also affects (pH-) bias energy.
     */
    public void propagateESVs(double temperature, double dt, double currentTimePs) {
        if (config.esvTempOverride) {
            temperature = ExtConstants.roomTemperature;
        } else {
            lastKnownTemperature = temperature;
        }
        if (esvList == null || esvList.isEmpty()) {
            return;
        }
        double[] dedl = ExtendedSystem.this.getEsvLambdaDerivatives(temperature, false);
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
    
    public void setTemperature(double set) {
        lastKnownTemperature = set;
    }

    /**
     *
     * @param esv
     */
    public void addVariable(ExtendedVariable esv) {
        if (!printedConfig) {
            ExtUtils.printEsvConfig();
            printedConfig = true;
        }
        SB.logfp(" ExtendedSystem acquired ESV: %s", esv);
        esvSystemActive = true;
        if (esvList == null) {
            esvList = new ArrayList<>();
        }
        esv.readyup(mola.getForceField());
        esvList.add(esv);
        numESVs = esvList.size();
        esvTerm = true;
        if (esv instanceof TitrationESV) {
            phTerm = true;
        }
        
        for (int i = 0; i < nAtomsExt; i++) {
            if (esv.viewUnsharedAtoms().contains(atomsExt[i])) {
                esvForUnshared[i] = esv;
            } else if (esv.viewSharedAtoms().contains(atomsExt[i])) {
                esvForShared[i] = esv;
            }
        }
        
        fg2bgIdx = new int[nAtomsExt];
        for (int i = 0; i < nAtomsExt; i++) {
            Atom atom = atomsExt[i];
            if (atom.getEsv() != null) {
                Atom bg = atom.getEsv().getBackgroundForAtom(atom);
                fg2bgIdx[i] = (bg != null) ? bg.getIndex() : -1;
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
    
    protected final void crashDump() {
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

    public double getTotalEsvLambdaDerivative() {
        return Arrays.stream(getEsvLambdaDerivatives(false)).sum();
    }
    
    /**
     * [numESVs] gradient with respect to each ESV
     */
    public double[] getEsvLambdaDerivatives(double temperature, boolean print) {
        double esvDeriv[] = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            esvDeriv[i] = getdEdL(i, temperature, print);
        }
        return esvDeriv;
    }

    public double[] getEsvLambdaDerivatives(boolean print) {
        return getEsvLambdaDerivatives(lastKnownTemperature, print);
    }

    public double getdEdL(int esvID) {
        return getdEdL(esvID, lastKnownTemperature, false);
    }

    public double getdEdL(int esvID, boolean print) {
        return getdEdL(esvID, lastKnownTemperature, print);
    }

    /**
     * Called reflectively from Groovy by script testESVs.
     */
    public double getdVdwdL(int esvID) {
        return (config.esvVdwFlag) ? vdw.getdEdEsv(esvID) : 0.0;
    }
    
    /**
     * Called reflectively from Groovy by script testESVs.
     */
    public double getdPermRealdL(int esvID) {
        return (config.esvPmeFlag) ? pme.getdEdEsv(esvID) : 0.0;
    }

    private double getdEdL(int esvID, double temperature, boolean print) {
        if (config.esvTempOverride) {
            temperature = ExtConstants.roomTemperature;
        }
        ExtendedVariable esv = esvList.get(esvID);
        double esvDeriv = 0.0;
        SB.logfn(" %s derivative components: ", esv.getName());
        if (config.esvBiasTerm) {
            esvDeriv += esv.getTotalBiasDeriv(temperature, print);
        }
        if (config.esvVdwFlag) {
            double dVdw = vdw.getdEdEsv(esvID);
            SB.logfn("  vdW    %d: %g", esvID, dVdw);
            esvDeriv += dVdw;
        }
        if (config.esvPmeFlag) {
            double dPme = pme.getdEdEsv(esvID);
            SB.logfn("  PME    %d: %g", esvID, dPme);
            esvDeriv += dPme;
        }
        if (config.esvScaleBonded) {
            double dBonded = esv.getBondedDeriv();
            esvDeriv += dBonded;
            if (print && config.esvDecomposeBonded) {
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
        if (!config.esvBiasTerm) {
            return "";
        }
        double discrBias = 0.0;
        double phBias = 0.0;
        for (ExtendedVariable esv : esvList) {
            discrBias += esv.getDiscrBias();
            if (esv instanceof TitrationESV) {
                phBias += ((TitrationESV) esv).getPhBias(lastKnownTemperature);
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
    private static final class Defaults {
        private Defaults() {}   // value singleton
        public static final ExtendedVariable esv = null;
        public static final int esvId = -1;
        public static final double lambda = 1.0;
        public static final double lambdaSwitch = 1.0;
        public static final double switchDeriv = 0.0;
    }
}
