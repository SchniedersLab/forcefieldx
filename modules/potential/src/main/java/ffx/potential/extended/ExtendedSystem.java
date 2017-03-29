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
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import edu.rit.pj.reduction.SharedDouble;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.SBLogger.SB;

/**
 *
 * @author slucore
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());
    public static boolean esvSystemActive = false;
    private int indexer = 0;

    /**
     * Stores configuration of system properties at instantiation of ExtendedSystem.
     */
    public final EsvConfig config;
    
    public static class EsvConfig {
        // Properties
        public final boolean esvVanDerWaals         = prop("esv-vdw", true);
        public final boolean esvElectrostatics      = prop("esv-pme", true);
        public final boolean esvBonded              = prop("esv-bonded", true);
        public final boolean esvScaleUnshared       = prop("esv-scaleUnshared", true);
        public final boolean biasTerm               = prop("esv-biasTerm", true);
        public final double discretizeBiasMagnitude = prop("esv-biasMagnitude", 1.0);
        public final boolean forceRoomTemperature   = prop("esv-forceRoomTemp", false);
        public final boolean esvDecomposeBonded     = prop("esv-decomposeBonded", true);
        public final boolean esvDecomposePme        = prop("esv-decomposePme", true);
        public final boolean esvAllowPropagation    = prop("esv-propagation", true);
        public final double thetaMass               = prop("esv-thetaMass", 1.0e-18);        // from OSRW, reasonably 100 a.m.u.
        public final double thetaFriction           = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps
        /* Nonlinear (switched) multipoles needs to be investigated before use. */
        public final boolean nonlinearMultipoles    = prop("esv-nonlinearMultipoles", false);
        public final boolean backgroundAtomsInFFE   = prop("esv-backgroundAtomsInFFE", true);
        public final boolean cloneXyzIndices        = prop("esv-cloneXyzIndices", true);
        
        public void print() {
			ExtUtils.printConfigSet("ExtendedSystem Config:", System.getProperties(), "esv");
        }
    }    
    
    // Atom Lists
    private Atom[] extendedAtoms;
    private int[] extendedMolecule;
    private int nAtomsExt;
    private ExtendedVariable[] esvForShared;
    private ExtendedVariable[] esvForUnshared;
    private int[] fg2bgIdx;

    // ESV variables
    public final Double constantSystemPh;
    private int numESVs;
    private List<ExtendedVariable> esvList;
    private boolean phTerm, vdwTerm, mpoleTerm;
    private Double currentTemperature;

    // Potential Objects
    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final VanDerWaals vdw;
    private final ParticleMeshEwaldQI pme;

    /**
     * @param constPh may be passed null if no TitrationESVs are added
     */
    public ExtendedSystem(MolecularAssembly mola, Double constPh) {
        if (mola == null) {
            throw new IllegalArgumentException();
        }
        this.config = new EsvConfig();
        this.mola = mola;
        this.constantSystemPh = constPh;
        Potential potential = mola.getPotentialEnergy();
        if (potential == null) {
            logger.severe("No potential energy found?");
        }
        if (!(potential instanceof ForceFieldEnergy)) {
            logger.severe("ExtendedSystem supported only for ForceFieldEnergy potentials.");
        }
        ffe = (ForceFieldEnergy) potential;

        ForceField ff = mola.getForceField();
        vdwTerm = ff.getBoolean(ForceField.ForceFieldBoolean.VDWTERM, true);
        mpoleTerm = ff.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, false);   // TODO set default true

        vdw = (vdwTerm && config.esvVanDerWaals) ? ffe.getVdwNode() : null;
        if (config.esvVanDerWaals && !vdwTerm) {
            logger.severe("Conflict: esvVdw without vdwTerm.");
        }
//        extendedNeighborList = (vdw != null) ? vdw.getExtendedNeighborList() : null;

        if (config.esvElectrostatics) {
            if (!mpoleTerm) {
                logger.severe("Conflict: esvPme without mpoleTerm.");
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
        currentTemperature = ExtConstants.roomTemperature;

        // Initialize atom arrays with the existing assembly.
        extendedAtoms = mola.getAtomArray();
        extendedMolecule = mola.getMoleculeNumbers();
        nAtomsExt = extendedAtoms.length;
        esvForUnshared = new ExtendedVariable[nAtomsExt];
        esvForShared = new ExtendedVariable[nAtomsExt];
        config.print();
    }
	    
    int requestIndexing() {
        return indexer++;
    }
    
    public void initializeBackgroundMultipoles(List<Atom> atomsBackground) {
        ExtUtils.initializeBackgroundMultipoles(atomsBackground, mola.getForceField());
    }
    
    public void populate(List<String> residueIDs) {
        // Locate the Residue identified by the given resid.
        Polymer[] polymers = mola.getChains();
        for (String token : residueIDs) {
            char chainID = token.charAt(0);
            int resNum = Integer.parseInt(token.substring(1));
            Residue target = null;
            for (Polymer p : polymers) {
                char pid = p.getChainID().charValue();
                if (pid == chainID) {
                    for (Residue res : p.getResidues()) {
                        if (res.getResidueNumber() == resNum) {
                            target = res;
                            break;
                        }
                    }
                    if (target != null) {
                        break;
                    }
                }
            }
            if (target == null) {
                SB.crash("Couldn't find target residue " + token);
            }

            MultiResidue titrating = TitrationUtils.titratingMultiresidueFactory(mola, target);
            TitrationESV esv = new TitrationESV(this, titrating);
            this.addVariable(esv);
        }
    }
    
    public void populate(String[] residueIDs) {
        populate(Arrays.asList(residueIDs));
    }
    
    public void populate(String residueIDs) {
        String[] tokens =
                  (residueIDs.split(".").length > 1) ? residueIDs.split(".")
                : (residueIDs.split(",").length > 1) ? residueIDs.split(",")
                : new String[]{residueIDs};
        populate(tokens);
    }
    
    public void setLambdas(String[] lambdas) {
        if (lambdas.length != numESVs) {
            throw new IllegalArgumentException();
        }
        for (int i = 0; i < lambdas.length; i++) {
            double lambda = Double.parseDouble(lambdas[i]);
            setLambda(i, lambda);
        }
    }    
    public void setLambdas(String lambdas) {
        String[] tokens =
                  (lambdas.split(",").length > 1) ? lambdas.split(",")
                : new String[]{lambdas};
        setLambdas(tokens);
    }
    
    public EsvConfig getConfig() {
        return config;
    }
    
    /**
     * All atoms of the fully-protonated system (not just those affected by this system).
     */
    public Atom[] getExtendedAtoms() {
        return extendedAtoms;
    }
    /**
     * Companion to getExtendedAtoms() for vdw::setAtoms and pme::setAtoms.
     */
    public int[] getExtendedMolecule() {
        return extendedMolecule;
    }
    /**
     * Whether the Atom at extendedAtoms[i] is affected by an ESV of this system.
     */
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
        return (isExtended(i)) ? getEsvForAtom(i).esvIndex : Defaults.esvId;
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
     * Allows PME to request updated scaling after multipole rotation.
     */
    public void updateMultipoles() {
        for (ExtendedVariable esv : this) {
            esv.updateMultipoleTypes();
        }
    }

    protected void updateListeners() {
        if (config.esvVanDerWaals) {
            vdw.updateEsvLambda();
        }
        if (config.esvElectrostatics) {
            pme.updateEsvLambda();
        }
    }

    public int size() {
        return esvList.size();
    }
    
    public final double getBiasEnergy() {
        return getBiasEnergy(currentTemperature);
    }

    /**
     * Get ESV biases such as discretization, pH, etc.
     * This method public and final for error-checking; new ESVs should override biasEnergy().
     */
    public final double getBiasEnergy(double temperature) {
        if (!config.biasTerm) {
            return 0.0;
        }
        if (esvList == null || esvList.isEmpty()) {
            logger.warning("Requested energy from empty/null esvSystem.");
            return 0.0;
        }
        if (config.forceRoomTemperature) {
            return biasEnergy(ExtConstants.roomTemperature);
        } else {
            return biasEnergy(temperature);
        }
    }

    private double biasEnergy(double temperature) {
        double biasEnergySum = 0.0;
        for (ExtendedVariable esv : this) {
            biasEnergySum += esv.getTotalBias(temperature, false);
        }
        return biasEnergySum;
    }

    public void propagateESVs(double dt, double currentTimePs) {
        propagateESVs(currentTemperature, dt, currentTimePs);
    }
    
    /**
     * Update the position of all ESV particles via langevin dynamics.
     * Temperature controls propagation speed and also affects (pH-) bias energy.
     */
    public void propagateESVs(double temperature, double dt, double currentTimePs) {
        if (config.forceRoomTemperature) {
            temperature = ExtConstants.roomTemperature;
        } else {
            currentTemperature = temperature;
        }
        if (esvList == null || esvList.isEmpty()) {
            return;
        }
        double[] dedl = ExtendedSystem.this.getEsvLambdaDerivatives(temperature, false);
        for (ExtendedVariable esv : esvList) {
            double oldLambda = esv.getLambda();
            esv.propagate(dedl[esv.esvIndex], dt, temperature);
            double newLambda = esv.getLambda();
            if (logger.isLoggable(Level.FINEST)) {
                logger.log(Level.FINEST, format(" Propagating ESV[%d]: %g --> %g @ psec,temp,bias: %g %g %.2f",
                        esv.esvIndex, oldLambda, newLambda,
                        currentTimePs, temperature, esv.getTotalBias(temperature, false)));
            }
        }
        updateListeners();
    }
    
    public void setTemperature(double set) {
        currentTemperature = set;
    }

    /**
     * Prefer ExtendedSystem::populate to manual ESV creation.
     */
    public void addVariable(ExtendedVariable esv) {
        SB.logfp(" ExtendedSystem acquired ESV: %s", esv);
        esvSystemActive = true;
        if (esvList == null) {
            esvList = new ArrayList<>();
        }
        if (esvList.contains(esv)) {
            SB.warning("Attempted to add duplicate variable %s to system.", esv.toString());
            return;
        }
        esvList.add(esv);
        
        numESVs = esvList.size();
        if (esv instanceof TitrationESV) {
            if (constantSystemPh == null) {
                SB.crash("Added TitrationESV to an ExtendedSystem lacking pH.");
            }
            phTerm = true;
        }
        
        for (int i = 0; i < nAtomsExt; i++) {
            if (esv.viewUnsharedAtoms().contains(extendedAtoms[i])) {
                esvForUnshared[i] = esv;
            } else if (esv.viewSharedAtoms().contains(extendedAtoms[i])) {
                esvForShared[i] = esv;
            }
        }
        
        fg2bgIdx = new int[nAtomsExt];
        for (int i = 0; i < nAtomsExt; i++) {
            Atom atom = extendedAtoms[i];
            if (atom.getEsv() != null) {
                Atom bg = atom.getEsv().getBackgroundForAtom(atom);
                fg2bgIdx[i] = (bg != null) ? bg.getIndex() : -1;
            }
        }
        
        updateListeners();
    }
    
    /**
     * Used only by ForceFieldEnergy and only once; we'd prefer to be rid of this altogether.
     * Background atoms are not true degrees of freedom.
     */
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
        return getEsvLambdaDerivatives(currentTemperature, print);
    }

    public double getdEdL(int esvID) {
        return getdEdL(esvID, currentTemperature, false);
    }

    public double getdEdL(int esvID, boolean print) {
        return getdEdL(esvID, currentTemperature, print);
    }

    /**
     * Called reflectively from Groovy by script testESVs.
     */
    public double getdVdwdL(int esvID) {
        if (!config.esvVanDerWaals) {
            logger.warning("Called for VdW ESV deriv while !esvVdw.");
            return 0.0;
        }
        return vdw.getdEdEsv(esvID);
    }
    
    /**
     * Called reflectively from Groovy by script testESVs.
     */
    public double getdPermRealdL(int esvID) {
        if (!config.esvElectrostatics) {
            logger.warning("Called for PermReal ESV deriv while !esvPme.");
            return 0.0;
        }
        return pme.getdPermRealdEsv(esvID);
    }
    
    public double getdPermRecipSelf_dL(int esvID) {
        if (!config.esvElectrostatics) {
            logger.warning("Called for Reciprocal ESV deriv while !esvPme.");
            return 0.0;
        }
        return pme.getdPermRecipSelf_dEsv(esvID);
    }
    
    public double getdPermRecipMpole_dL(int esvID) {
        if (!config.esvElectrostatics) {
            logger.warning("Called for Reciprocal ESV deriv while !esvPme.");
            return 0.0;
        }
        return pme.getdPermRecipMpole_dEsv(esvID);
    }
    
    public void printdEdL() {
        for (int i = 0; i < numESVs; i++) {
            getdEdL(i, true);
        }
    }

    private double getdEdL(int esvID, double temperature, boolean print) {
        if (config.forceRoomTemperature) {
            temperature = ExtConstants.roomTemperature;
        }
        ExtendedVariable esv = esvList.get(esvID);
        double esvDeriv = 0.0;
        if (config.biasTerm) {
            final double dBias = esv.getTotalBiasDeriv(temperature, print);
            SB.logfn("   %-18.18s %9.4f", "Biases:", dBias);
            final double dDiscr = esv.getDiscrBiasDeriv();
            SB.logfn("     %-16.16s %9.4f", "Discretizer:", dDiscr);
            if (esv instanceof TitrationESV) {
                final double dPh = ((TitrationESV) esv).getPhBiasDeriv(temperature);
                SB.logfn("     %-16.16s %9.4f", "Acidostat:", dPh);
            }
            esvDeriv += dBias;
        }
        if (config.esvVanDerWaals) {
            final double dVdw = vdw.getdEdEsv(esvID);
            SB.logfn("   %-18.18s %9.4f", "van der Waals:", dVdw);
            esvDeriv += dVdw;
        }
        if (config.esvElectrostatics) {
            final double dPme = pme.getdEdEsv(esvID);
            SB.logfn("   %-18s %9.4f", "Electrostatics:", dPme);
            if (config.esvDecomposePme) {
                double permReal = pme.getdPermRealdEsv(esvID);
                double permRecipSelf = pme.getdPermRecipSelf_dEsv(esvID);
                double permRecipMpole = pme.getdPermRecipMpole_dEsv(esvID);
                SB.logfn("     %-16s %9.4f", "PermReal:", permReal);
                SB.logfn("     %-16s %9.4f", "PermRcpSelf:", permRecipSelf);
                SB.logfn("     %-16s %9.4f", "PermRecipMpole:", permRecipMpole);
            }
            esvDeriv += dPme;
        }
        if (config.esvBonded) {
            final double dBonded = esv.getBondedDeriv();
            SB.logfn("   %-18s %9.4f", "Bonded:", dBonded);
            esvDeriv += dBonded;
            /* If desired, decompose bonded contribution into component types from foreground and background. */
            if (config.esvDecomposeBonded) {
                // Foreground portion:
                double fgSum = 0.0;
                HashMap<Class<? extends BondedTerm>,SharedDouble> fgMap =
                        esv.getBondedDerivDecomp();
                for (SharedDouble dub : fgMap.values()) {
                    fgSum += dub.get();
                }
                SB.logf("     %-16s %9.4f", "Foreground:" , fgSum);
                for (Class<? extends BondedTerm> clas : fgMap.keySet()) {
                    SB.nlogf("       %-14s %9.4f",
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
                SB.nlogf("     %-16s %9.4f", "Background:" , bgSum);
                for (Class<? extends BondedTerm> clas : bgMap.keySet()) {
                    SB.nlogf("       %-14s %9.4f",
                            clas.getName().replaceAll("ffx.potential.bonded.", "") + ":",
                            bgMap.get(clas).get());
                }
            }
        }
        if (Double.isNaN(esvDeriv) || !Double.isFinite(esvDeriv)) {
            SB.warning("NaN/Inf lambda derivative: %s", this);
        }
        SB.headern(" %-20s %9.4f", format("dUd%s:", esv.getName()), esvDeriv);
        SB.printIf(print);
        return esvDeriv;
    }
    
    public ExtendedVariable getEsv(int esvID) {
        return esvList.get(esvID);
    }

    public String getBiasDecomposition() {
        if (!config.biasTerm) {
            return "";
        }
        double discrBias = 0.0;
        double phBias = 0.0;
        for (ExtendedVariable esv : esvList) {
            discrBias += esv.getDiscrBias();
            if (esv instanceof TitrationESV) {
                phBias += ((TitrationESV) esv).getPhBias(currentTemperature);
            }
        }
        return  format("    %-16s %16.8f\n", "Discretizer", discrBias) +
                format("    %-16s %16.8f\n", "Acidostat", phBias);
    }

    public String getLambdaList() {
        if (numESVs < 1) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        sb.append(format("ExtendedVariables: "));
        for (int i = 0; i < numESVs; i++) {
            if (i > 0) {
                sb.append(", ");
            }
            sb.append(format("%6.4f", esvList.get(i).getLambda()));
        }
        return sb.toString();
    }
    
    public double[] getLambdas() {
        double[] lambdas = new double[numESVs];
        for (int i = 0; i < lambdas.length; i++) {
            lambdas[i] = esvList.get(i).getLambda();
        }
        return lambdas;
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
        public static final double switchDeriv = 1.0;
    }
}
