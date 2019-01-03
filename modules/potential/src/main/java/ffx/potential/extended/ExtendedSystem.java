/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
import ffx.potential.PotentialComponent;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;
import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;

/**
 * <p>ExtendedSystem class.</p>
 *
 * @author Stephen LuCore
 *
 * @since 1.0
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());
    /**
     * Constant <code>esvSystemActive=false</code>
     */
    public static boolean esvSystemActive = false;
    private int indexer = 0;
    private static final StringBuilder sb = new StringBuilder();

    /**
     * Stores configuration of system properties at instantiation of
     * ExtendedSystem.
     */
    public final ExtendedSystemConfig config;

    /**
     * Stores all the default values listed in prop() calls below.
     */
    public static final ExtendedSystemConfig DefaultConfig;

    static {
        /*
         * During static initialization, clear System properties of "esv." keys
         * to load the default ExtendedSystemConfig, then put them back.
         */
        synchronized (System.getProperties()) {
            HashMap<String, String> bak = new HashMap<>();
            System.getProperties().forEach((Object k, Object v) -> {
                String key = (String) k;
                if (key.startsWith("esv.")) {
                    bak.put(key, (String) v);
                }
            });
            for (String k : bak.keySet()) {
                System.clearProperty(k);
            }
            DefaultConfig = new ExtendedSystemConfig();
            System.getProperties().putAll(bak);
        }
    }

    public static class ExtendedSystemConfig {

        /**
         * Use polarizability scaling by default.
         */
        public final boolean bonded = prop("esv.bonded", true);
        public final boolean vanDerWaals = prop("esv.vanDerWaals", true);
        public final boolean electrostatics = prop("esv.electrostatics", true);
        public final boolean polarization = prop("esv.polarization", true);
        public final boolean biasTerm = prop("esv.biasTerm", true);
        public final boolean verbose = prop("esv.verbose", false);
        public final boolean decomposeBonded = prop("esv.decomposeBonded", false);
        public final boolean decomposeDeriv = prop("esv.decomposeDeriv", false);

        /**
         * Note that without the Lambda-Switch, the derivative dPol/dEsv is
         * incorrect at L=0.0 and L=1.0
         */
        public final boolean allowLambdaSwitch = prop("esv.allowLambdaSwitch", true);

        public final boolean nonlinearMultipoles = prop("esv.nonlinearMultipoles", false);    // sigmoid lambda Mpole switch
        public final double discrBias = prop("esv.biasMagnitude", 1.0);
        public final boolean forceRoomTemp = prop("esv.forceRoomTemp", false);
        public final boolean propagation = prop("esv.propagation", true);
        public final double thetaMass = prop("esv.thetaMass", 1.0e-18);    // from OSRW, reasonably 100 a.m.u.
        public final double thetaFriction = prop("esv.thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps

        // Options below are untested and/or dangerous if changed.
        public final boolean cloneXyzIndices = prop("esv.cloneXyzIndices", true);    // set bg_idx = fg_idx

        public static void print(ExtendedSystemConfig config) {
            List<Field> fields = Arrays.asList(ExtendedSystemConfig.class.getDeclaredFields());
            Collections.sort(fields, (Field t, Field t1)
                    -> String.CASE_INSENSITIVE_ORDER.compare(t.getName(), t1.getName()));
            for (int i = 0, col = 0; i < fields.size(); i++) {
                if (++col > 3) {
                    sb.append(format("\n"));
                    col = 1;
                }
                String key = fields.get(i).getName() + ":";
                try {
                    Object obj = fields.get(i).get(config);
                    sb.append(format(" %-30s %7.7s          ", key, obj));
                } catch (IllegalAccessException ignored) {
                }
            }
            sb.append(format(" %-30s %7.7s          %-30s %7.7s          %-30s %7.7s"
                            + "\n %-30s %7.7s          %-30s %7.7s          %-30s %7.7s"
                            + "\n %-30s %7.7s",
                    "polarization", System.getProperty("polarization"),
                    "scf-algorithm", System.getProperty("scf-algorithm"),
                    "polar-eps", System.getProperty("polar-eps"),
                    "use-charges", System.getProperty("use-charges"),
                    "use-dipoles", System.getProperty("use-dipoles"),
                    "use-quadrupoles", System.getProperty("use-quadrupoles"),
                    "grid-method", System.getProperty("grid-method")));
            sb.append(format("\n"));
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
    private Double constantSystemPh;
    private int numESVs;
    private List<ExtendedVariable> esvList;
    private boolean vdwTerm, mpoleTerm;
    private Double currentTemperature;

    // Potential Objects
    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final VanDerWaals vdw;
    private final ParticleMeshEwaldQI pme;

    /**
     * Construct extended system with a default configuration and/or system
     * properties.
     *
     * @param mola a {@link ffx.potential.MolecularAssembly} object.
     */
    public ExtendedSystem(MolecularAssembly mola) {
        this(mola, new ExtendedSystemConfig());
    }

    /**
     * Construct extended system with the provided configuration.
     *
     * @param mola   a {@link ffx.potential.MolecularAssembly} object.
     * @param config a {@link ffx.potential.extended.ExtendedSystem.ExtendedSystemConfig} object.
     */
    public ExtendedSystem(MolecularAssembly mola, ExtendedSystemConfig config) {
        if (mola == null) {
            throw new IllegalArgumentException();
        }
        this.config = config;
        this.mola = mola;
        Potential potential = mola.getPotentialEnergy();
        if (potential == null) {
            logger.severe("No potential energy found?");
        }
        if (!(potential instanceof ForceFieldEnergy)) {
            logger.warning("ExtendedSystem supported only for ForceFieldEnergy potentials.");
            throw new ClassCastException();
        }
        ffe = (ForceFieldEnergy) potential;

        ForceField ff = mola.getForceField();
        vdwTerm = ff.getBoolean(ForceField.ForceFieldBoolean.VDWTERM, true);
        mpoleTerm = ff.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, true);

        VanDerWaals vdwNode = null;
        ParticleMeshEwaldQI pmeNode = null;
        try {
            Method getVdw = ffe.getClass().getDeclaredMethod("getVdwNode");
            Method getPme = ffe.getClass().getDeclaredMethod("getPmeNode");
            getVdw.setAccessible(true);
            getPme.setAccessible(true);
            vdwNode = (VanDerWaals) getVdw.invoke(ffe);
            pmeNode = (ParticleMeshEwaldQI) getPme.invoke(ffe);
        } catch (ReflectiveOperationException | IllegalArgumentException | SecurityException ex) {
            /* To run ExtendedSystem under a SecurityManager, expose as public the methods
             * ForceFieldEnergy::getVdwNode and ForceFieldEnergy::getPmeNode.			*/
            logger.log(Level.SEVERE, "Extended system unable to reflectively obtain VdW, PME nodes.", ex);
        } catch (ClassCastException ex) {
            logger.warning("Extended system supported only for quasi-internal ParticleMeshEwald.");
            throw ex;
        }
        vdw = (vdwTerm && config.vanDerWaals) ? vdwNode : null;
        pme = (mpoleTerm && config.electrostatics) ? pmeNode : null;
        if (config.vanDerWaals && !vdwTerm) {
            logger.severe("Conflict: esvVdw without vdwTerm.");
        }
        if (config.electrostatics && !mpoleTerm) {
            logger.severe("Conflict: esvElectrostatics without mpoleTerm.");
        }

        esvList = new ArrayList<>();
        currentTemperature = ExtConstants.roomTemperature;

        // Initialize atom arrays with the existing assembly.
        extendedAtoms = mola.getAtomArray();
        extendedMolecule = mola.getMoleculeNumbers();
        nAtomsExt = extendedAtoms.length;
        esvForUnshared = new ExtendedVariable[nAtomsExt];
        esvForShared = new ExtendedVariable[nAtomsExt];
        if (config.verbose) {
            ExtendedSystemConfig.print(config);
        }
    }

    int requestIndexing() {
        return indexer++;
    }

    /**
     * <p>initializeBackgroundMultipoles.</p>
     *
     * @param atomsBackground a {@link java.util.List} object.
     */
    public void initializeBackgroundMultipoles(List<Atom> atomsBackground) {
        ExtUtils.initializeBackgroundMultipoles(atomsBackground, mola.getForceField());
    }

    /**
     * <p>populate.</p>
     *
     * @param residueIDs a {@link java.util.List} object.
     */
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
                logger.severe("Couldn't find target residue " + token);
            }

            MultiResidue titrating = TitrationUtils.titratingMultiresidueFactory(mola, target);
            TitrationESV esv = new TitrationESV(this, titrating);
            this.addVariable(esv);
        }
    }

    /**
     * <p>populate.</p>
     *
     * @param residueIDs an array of {@link java.lang.String} objects.
     */
    public void populate(String[] residueIDs) {
        populate(Arrays.asList(residueIDs));
    }

    /**
     * <p>populate.</p>
     *
     * @param residueIDs a {@link java.lang.String} object.
     */
    public void populate(String residueIDs) {
        String[] tokens
                = (residueIDs.split(".").length > 1) ? residueIDs.split(".")
                : (residueIDs.split(",").length > 1) ? residueIDs.split(",")
                : new String[]{residueIDs};
        populate(tokens);
    }

    /**
     * <p>Getter for the field <code>config</code>.</p>
     *
     * @return a {@link ffx.potential.extended.ExtendedSystem.ExtendedSystemConfig} object.
     */
    public ExtendedSystemConfig getConfig() {
        return config;
    }

    /**
     * All atoms of the fully-protonated system (not just those affected by this
     * system).
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getExtendedAtoms() {
        return extendedAtoms;
    }

    /**
     * Companion to getExtendedAtoms() for vdw::setAtoms and pme::setAtoms.
     *
     * @return an array of {@link int} objects.
     */
    public int[] getExtendedMolecule() {
        return extendedMolecule;
    }

    /**
     * Whether the Atom at extendedAtoms[i] is affected by an ESV of this
     * system.
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isExtended(int i) {
        return isShared(i) || isUnshared(i);
    }

    /**
     * <p>isAlphaScaled.</p>
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isAlphaScaled(int i) {
        return isExtended(i) && isTitratableHydrogen(extendedAtoms[i]);
    }

    /**
     * <p>isShared.</p>
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isShared(int i) {
        return esvForShared[i] != null;
    }

    /**
     * <p>isUnshared.</p>
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isUnshared(int i) {
        return esvForUnshared[i] != null;
    }

    /**
     * <p>getEsvForAtom.</p>
     *
     * @param i a int.
     * @return a {@link ffx.potential.extended.ExtendedVariable} object.
     */
    public ExtendedVariable getEsvForAtom(int i) {
        return (isShared(i)) ? esvForShared[i]
                : (isUnshared(i)) ? esvForUnshared[i]
                : null;
    }

    /**
     * <p>getEsvIndex.</p>
     *
     * @param i a int.
     * @return a {@link java.lang.Integer} object.
     */
    public Integer getEsvIndex(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).esvIndex : null;
    }

    /**
     * <p>getLambda.</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getLambda(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getLambda() : Defaults.lambda;
    }

    /**
     * <p>getLambdaSwitch.</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getLambdaSwitch(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getLambdaSwitch() : Defaults.lambdaSwitch;
    }

    /**
     * <p>getSwitchDeriv.</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getSwitchDeriv(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getSwitchDeriv() : Defaults.switchDeriv;
    }

    /**
     * <p>setLambda.</p>
     *
     * @param esvId  a int.
     * @param lambda a double.
     */
    public void setLambda(int esvId, double lambda) {
        if (esvId >= numESVs) {
            logger.warning("Requested an invalid ESV id.");
        }
        getEsv(esvId).setLambda(lambda);
        updateListeners();
    }

    /**
     * <p>setLambda.</p>
     *
     * @param esvIdChar a char.
     * @param lambda    a double.
     */
    public void setLambda(char esvIdChar, double lambda) {
        setLambda(esvIdChar - 'A', lambda);
    }

    /**
     * <p>getLambda.</p>
     *
     * @param esvIdChar a char.
     * @return a double.
     */
    public double getLambda(char esvIdChar) {
        return getEsv(esvIdChar - 'A').getLambda();
    }

    /**
     * <p>updateListeners.</p>
     */
    protected void updateListeners() {
        if (config.vanDerWaals) {
            vdw.updateEsvLambda();
        }
        if (config.electrostatics) {
            pme.updateEsvLambda();
        }
    }

    /**
     * <p>size.</p>
     *
     * @return a int.
     */
    public int size() {
        return esvList.size();
    }

    /**
     * <p>getBiasEnergy.</p>
     *
     * @return a double.
     */
    public final double getBiasEnergy() {
        return getBiasEnergy(currentTemperature);
    }

    /**
     * Get ESV biases such as discretization, pH, etc. This method public and
     * final for error-checking; new ESVs should override biasEnergy().
     *
     * @param temperature a double.
     * @return a double.
     */
    public final double getBiasEnergy(double temperature) {
        if (!config.biasTerm) {
            return 0.0;
        }
        if (esvList == null || esvList.isEmpty()) {
            logger.warning("Requested energy from empty/null esvSystem.");
            return 0.0;
        }
        if (config.forceRoomTemp) {
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

    /**
     * Update the position of all ESV particles via langevin dynamics.
     * Temperature controls propagation speed and also affects (pH-) bias
     * energy.
     *
     * @param temperature   a double.
     * @param dt            a double.
     * @param currentTimePs a double.
     */
    public void propagateESVs(double temperature, double dt, double currentTimePs) {
        if (config.forceRoomTemp) {
            temperature = ExtConstants.roomTemperature;
        } else {
            currentTemperature = temperature;
        }
        if (esvList == null || esvList.isEmpty()) {
            return;
        }
        double[] dedl = ExtendedSystem.this.getDerivatives();
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

    /**
     * <p>setTemperature.</p>
     *
     * @param set a double.
     */
    public void setTemperature(double set) {
        currentTemperature = set;
    }

    /**
     * Prefer ExtendedSystem::populate to manual ESV creation.
     *
     * @param esv a {@link ffx.potential.extended.ExtendedVariable} object.
     */
    public void addVariable(ExtendedVariable esv) {
        esvSystemActive = true;
        if (esvList == null) {
            esvList = new ArrayList<>();
        }
        if (esvList.contains(esv)) {
            logger.warning(format("Attempted to add duplicate variable %s to system.", esv.toString()));
            return;
        }
        esvList.add(esv);

        numESVs = esvList.size();
        if (esv instanceof TitrationESV) {
            if (constantSystemPh == null) {
                logger.severe("Set ExtendedSystem (constant) pH before adding TitrationESVs.");
            }
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
            if (isExtended(i)) {
                Atom bg = atom.getEsv().getBackgroundForAtom(atom);
                fg2bgIdx[i] = (bg != null) ? bg.getIndex() : -1;
            }
        }

        updateListeners();
    }

    /**
     * <p>setConstantPh.</p>
     *
     * @param pH a double.
     */
    public void setConstantPh(double pH) {
        if (constantSystemPh != null) {
            logger.severe("Attempted to modify an existing constant pH value.");
        }
        constantSystemPh = pH;
    }

    /**
     * <p>getConstantPh.</p>
     *
     * @return a double.
     */
    public double getConstantPh() {
        if (constantSystemPh == null) {
            logger.severe("Requested an unset system pH value.");
        }
        return constantSystemPh.doubleValue();
    }

    /**
     * Used only by ForceFieldEnergy and only once; we'd prefer to be rid of
     * this altogether. Background atoms are not true degrees of freedom.
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
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

    /**
     * <p>getExtendedAndBackgroundMolecule.</p>
     *
     * @return an array of {@link int} objects.
     */
    public int[] getExtendedAndBackgroundMolecule() {
        Atom[] mega = getExtendedAndBackgroundAtoms();
        int[] molecule = new int[mega.length];
        for (int i = 0; i < molecule.length; i++) {
            molecule[i] = mega[i].getMoleculeNumber();
        }
        return molecule;
    }

    /**
     * <p>crashDump.</p>
     */
    protected final void crashDump() {
        sb.append("*************");
        sb.append(" Crash Dump:");
        sb.append("   All Atoms:");
        for (Atom atom : mola.getAtomArray()) {
            sb.append(format("     %s", atom.toString()));
        }
        logger.info(sb.toString());
        for (ExtendedVariable esv : esvList) {
            esv.describe();
        }
        sb.append("*************");
        logger.info(sb.toString());
    }

    /**
     * Potential gradient with respect to each ESV; used to propagate langevin
     * dynamics.
     *
     * @return an array of {@link double} objects.
     */
    public double[] getDerivatives() {
        double esvDeriv[] = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            esvDeriv[i] = getDerivative(i);
        }
        return esvDeriv;
    }

    /**
     * <p>getDerivativeComponent.</p>
     *
     * @param dd    a {@link ffx.potential.PotentialComponent} object.
     * @param esvID a int.
     * @return a double.
     */
    public double getDerivativeComponent(PotentialComponent dd, int esvID) {
        switch (dd) {
            case Topology:
            case ForceFieldEnergy:
                return getDerivative(esvID);
            case VanDerWaals:
                return vdw.getEsvDerivative(esvID);
            case Bonded:
                return esvList.get(esvID).getBondedDeriv();
            case Bias:
            case pHMD:
                return esvList.get(esvID).getTotalBiasDeriv(currentTemperature, false);
            case Discretizer:
                return esvList.get(esvID).getDiscrBiasDeriv();
            case Acidostat:
                return ((TitrationESV) esvList.get(esvID)).getPhBiasDeriv(currentTemperature);
            case Multipoles:
                return pme.getEsvDerivative(esvID);
            case Permanent:
                return pme.getEsvDeriv_Permanent(esvID);
            case PermanentRealSpace:
                return pme.getEsvDeriv_PermReal(esvID);
            case PermanentSelf:
                return pme.getEsvDeriv_PermSelf(esvID);
            case PermanentReciprocal:
                return pme.getEsvDeriv_PermRecip(esvID);
            case Induced:
                return pme.getEsvDeriv_Induced(esvID);
            case InducedRealSpace:
                return pme.getEsvDeriv_IndReal(esvID);
            case InducedSelf:
                return pme.getEsvDeriv_IndSelf(esvID);
            case InducedReciprocal:
                return pme.getEsvDeriv_IndRecip(esvID);
            default:
                throw new AssertionError(dd.name());
        }
    }

    private double getDerivative(int esvID) {
        final double temperature = (config.forceRoomTemp)
                ? ExtConstants.roomTemperature : currentTemperature;
        final boolean p = config.decomposeDeriv;
        ExtendedVariable esv = esvList.get(esvID);
        double esvDeriv = 0.0;
        final String format = " %-20.20s %2.2s %9.4f";
        if (config.biasTerm) {
            final double dBias = esv.getTotalBiasDeriv(temperature, true);
            if (p) {
                sb.append(format("  Biases:", "", dBias));
            }
            final double dDiscr = esv.getDiscrBiasDeriv();
            if (p) {
                sb.append(format("    Discretizer:", ">", dDiscr));
            }
            if (esv instanceof TitrationESV) {
                final double dPh = ((TitrationESV) esv).getPhBiasDeriv(temperature);
                if (p) {
                    sb.append(format("    Acidostat:", ">", dPh));
                }
            }
            esvDeriv += dBias;
        }
        if (config.vanDerWaals) {
            final double dVdw = vdw.getEsvDerivative(esvID);
            if (p) {
                sb.append(format("  VanDerWaals:", "", dVdw));
            }
            esvDeriv += dVdw;
        }
        if (config.electrostatics) {
            final double permanent = pme.getEsvDeriv_Permanent(esvID);
            esvDeriv += permanent;
            if (p) {
                sb.append(format("  PermanentElec:", "", permanent));
            }
            double permReal = pme.getEsvDeriv_PermReal(esvID);
            double permSelf = pme.getEsvDeriv_PermSelf(esvID);
            double permRecip = pme.getEsvDeriv_PermRecip(esvID);
            if (p) {
                sb.append(format("    PermReal:", ">", permReal));
            }
            if (p) {
                sb.append(format("    PermRcpSelf:", ">", permSelf));
            }
            if (p) {
                sb.append(format("    PermRecipMpole:", ">", permRecip));
            }
            if (config.polarization) {
                final double induced = pme.getEsvDeriv_Induced(esvID);
                esvDeriv += induced;
                if (p) {
                    sb.append(format("  Polarization:", "", induced));
                }
                double indReal = pme.getEsvDeriv_IndReal(esvID);
                double indSelf = pme.getEsvDeriv_IndSelf(esvID);
                double indRecip = pme.getEsvDeriv_IndRecip(esvID);
                if (p) {
                    sb.append(format("    IndReal:", ">", indReal));
                }
                if (p) {
                    sb.append(format("    IndSelf:", ">", indSelf));
                }
                if (p) {
                    sb.append(format("    IndRecip:", ">", indRecip));
                }
            }
        }
        if (config.bonded) {
            final double dBonded = esv.getBondedDeriv();
            if (p) {
                sb.append(format("  Bonded:", "", dBonded));
            }
            esvDeriv += dBonded;
            /* If desired, decompose bonded contribution into component types from foreground and background. */
            if (config.decomposeBonded) {
                // Foreground portion:
                double fgSum = 0.0;
                HashMap<Class<? extends BondedTerm>, SharedDouble> fgMap
                        = esv.getBondedDerivDecomp();
                for (SharedDouble dub : fgMap.values()) {
                    fgSum += dub.get();
                }
                if (p) {
                    sb.append(format("    Foreground:", ">", fgSum));
                }
                for (Class<? extends BondedTerm> clas : fgMap.keySet()) {
                    if (p) {
                        sb.append(format("      " + clas.getName().replaceAll("ffx.potential.bonded.", "") + ":",
                                ">>", fgMap.get(clas).get()));
                    }
                }
                // Background portion:
                double bgSum = 0.0;
                HashMap<Class<? extends BondedTerm>, SharedDouble> bgMap
                        = esv.getBackgroundBondedDerivDecomp();
                for (SharedDouble dub : bgMap.values()) {
                    bgSum += dub.get();
                }
                if (p) {
                    sb.append(format("    Background:", ">", bgSum));
                }
                for (Class<? extends BondedTerm> clas : bgMap.keySet()) {
                    if (p) {
                        sb.append(format("      " + clas.getName().replaceAll("ffx.potential.bonded.", "") + ":",
                                ">>", bgMap.get(clas).get()));
                    }
                }
            }
        }
        if (Double.isNaN(esvDeriv) || !Double.isFinite(esvDeriv)) {
            logger.warning(format("NaN/Inf lambda derivative: %s", this));
        }
        if (p) {
            sb.insert(0, format(" %-21.21s %-2.2s %9.4f", format("dUd%s:", esv.getName()), "", esvDeriv));
        }
        if (p) {
            logger.info(sb.toString());
        }
        return esvDeriv;
    }

    /**
     * <p>getEnergyComponent.</p>
     *
     * @param component a {@link ffx.potential.PotentialComponent} object.
     * @return a double.
     */
    public double getEnergyComponent(PotentialComponent component) {
        double uComp = 0.0;
        switch (component) {
            case Bias:
            case pHMD:
                return getBiasEnergy();
            case Discretizer:
                for (int i = 0; i < numESVs; i++) {
                    uComp += esvList.get(i).getDiscrBias();
                }
                return uComp;
            case Acidostat:
                for (int i = 0; i < numESVs; i++) {
                    uComp += ((TitrationESV) esvList.get(i)).getPhBias(currentTemperature);
                }
                return uComp;
            default:
                throw new AssertionError(component.name());
        }
    }

    /**
     * <p>getEsv.</p>
     *
     * @param esvID a int.
     * @return a {@link ffx.potential.extended.ExtendedVariable} object.
     */
    public ExtendedVariable getEsv(int esvID) {
        return esvList.get(esvID);
    }

    /**
     * <p>getBiasDecomposition.</p>
     *
     * @return a {@link java.lang.String} object.
     */
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
        return format("    %-16s %16.8f\n", "Discretizer", discrBias)
                + format("    %-16s %16.8f\n", "Acidostat", phBias);
    }

    /**
     * <p>getLambdaList.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getLambdaList() {
        if (numESVs < 1) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numESVs; i++) {
            if (i > 0) {
                sb.append(", ");
            }
            sb.append(format("%6.4f", esvList.get(i).getLambda()));
        }
        return sb.toString();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Allows simple iteration over ESV via "for (ExtendedVariable :
     * ExtendedSystem)".
     */
    @Override
    public Iterator<ExtendedVariable> iterator() {
        return esvList.iterator();
    }

    /**
     * These populate the order-n preloaded lambda parameter arrays in VdW and
     * PME in the absence of an attached ESV.
     */
    private static final class Defaults {

        private Defaults() {
        }   // value singleton

        public static final ExtendedVariable esv = null;
        public static final Integer esvId = null;
        public static final double lambda = 1.0;
        public static final double lambdaSwitch = 1.0;
        public static final double switchDeriv = 1.0;
    }
}
