// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.potential.extended;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;
import static java.lang.String.format;

import edu.rit.pj.reduction.SharedDouble;
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

import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;

import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.lang.reflect.Field;
import java.util.*;
import java.util.logging.Logger;

/**
 * ExtendedSystem class.
 *
 * @author Stephen LuCore
 * @since 1.0
 */
public class ExtendedSystem implements Iterable<ExtendedVariable> {

    /**
     * Stores all the default values listed in prop() calls below.
     */
    public static final ExtendedSystemConfig DefaultConfig;

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());

    public static final double THETA_MASS = 10.0; //Atomic Mass Units

    public static final double THETA_FRICTION = 5.0; // 1/psec

    /**
     * Constant <code>esvSystemActive=false</code>
     */
    public static boolean esvSystemActive = false;

    static {
        /*
         * During static initialization, clear System properties of "esv." keys
         * to load the default ExtendedSystemConfig, then put them back.
         */
        synchronized (System.getProperties()) {
            HashMap<String, String> bak = new HashMap<>();
            System.getProperties()
                    .forEach(
                            (Object k, Object v) -> {
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

    /**
     * Stores configuration of system properties at instantiation of ExtendedSystem.
     */
    public final ExtendedSystemConfig config;
    /**
     * MolecularAssembly instance.
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * VanDerWaals instance.
     */
    private final VanDerWaals vanDerWaals;
    /**
     * PME instance.
     */
    private final ParticleMeshEwaldQI particleMeshEwaldQI;
    /**
     * Count of ESV variables. TODO: Repalce wtih esvList.size()?
     */
    private int indexer = 0;
    /**
     * List of extended atoms.
     */
    private final Atom[] extendedAtoms;
    /**
     * List of extended molecules.
     */
    private final int[] extendedMolecule;
    /**
     * Number of extended atoms.
     */
    private final int nAtomsExt;
    /**
     * ExtendedVariable list for shared atoms.
     */
    private final ExtendedVariable[] esvForShared;
    /**
     * Extended Variable list for unshared atoms.
     */
    private final ExtendedVariable[] esvForUnshared;
    /**
     * System PH.
     */
    private Double constantSystemPh;
    /**
     * Number of ESV.
     */
    private int numESVs;
    /**
     * List of ESV instances.
     */
    private List<ExtendedVariable> esvList;
    /**
     * Target system temperature.
     */
    private Double currentTemperature;
    /**
     * Current value of theta for each ESV.
     */
    public double[] theta_position;
    /**
     * Current theta velocity for each ESV.
     */
    public double[] theta_velocity;
    /**
     * Current theta acceleration for each ESV.
     */
    public double[] theta_accel;
    /**
     * Mass of each theta particle.
     */
    public double[] theta_mass;
    /**
     * Friction for the ESV system
     */
    public final double thetaFriction;

    /**
     * Construct extended system with a default configuration and/or system properties.
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
        this.config = config;
        this.molecularAssembly = mola;

        ForceFieldEnergy forceFieldEnergy = mola.getPotentialEnergy();
        if (forceFieldEnergy == null) {
            logger.severe("No potential energy found?");
        }

        CompositeConfiguration properties = molecularAssembly.getProperties();
        thetaFriction = properties.getDouble("esv.friction",ExtendedSystem.THETA_FRICTION);

        ForceField ff = mola.getForceField();
        boolean vdwTerm = ff.getBoolean("VDWTERM", true);
        boolean mpoleTerm = ff.getBoolean("MPOLETERM", true);

        VanDerWaals vdwNode = forceFieldEnergy.getVdwNode();
        ParticleMeshEwaldQI pmeNode = forceFieldEnergy.getPmeQiNode();
        vanDerWaals = (vdwTerm && config.vanDerWaals) ? vdwNode : null;
        particleMeshEwaldQI = (mpoleTerm && config.electrostatics) ? pmeNode : null;
        if (config.vanDerWaals && !vdwTerm) {
            logger.severe("Conflict: esvVdw without vdwTerm.");
        }
        if (config.electrostatics && !mpoleTerm) {
            logger.severe("Conflict: esvElectrostatics without mpoleTerm.");
        }

        esvList = new ArrayList<>();
        currentTemperature = Constants.ROOM_TEMPERATURE;

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

        updateListeners();
    }

    /**
     * initializeBackgroundMultipoles.
     *
     * @param atomsBackground a {@link java.util.List} object.
     */
    public void initializeBackgroundMultipoles(List<Atom> atomsBackground) {
        ExtUtils.initializeBackgroundMultipoles(atomsBackground, molecularAssembly.getForceField());
    }

    public void createMDThetaArrays() {
        theta_position = new double[numESVs];
        theta_velocity = new double[numESVs];
        theta_accel = new double[numESVs];
        theta_mass = new double[numESVs];

        //Theta masses should always be the same for each ESV
        double mass = getEsv(0).getThetaMass();
        Arrays.fill(theta_mass, mass);
        collectThetaValues();
    }

    public MolecularAssembly getMolecularAssembly() {
        return molecularAssembly;
    }

    /**
     * Iterate over all Extended Variables in Extended System and collect thetas, velocities, and accelerations into arrays.
     */
    public void collectThetaValues() {
        for (ExtendedVariable esv : esvList) {
            theta_position[esv.getEsvIndex()] = esv.getTheta();
            theta_velocity[esv.getEsvIndex()] = esv.getThetaVelocity();
            theta_accel[esv.getEsvIndex()] = esv.getThetaAccel();
      /*logger.info(format("ESV: %d Theta: %g, Theta_v: %g, Theta_a: %g",
              esv.getEsvIndex(),esv.getTheta(),esv.getTheta_velocity(),esv.getTheta_accel()));*/
        }
    }

    /**
     * Send all theta information stored in Extended System arrays back to respective Extended Variables.
     */
    public void setThetaValues() {
        for (ExtendedVariable esv : esvList) {
            esv.setTheta(theta_position[esv.getEsvIndex()]);
            esv.setThetaVelocity(theta_velocity[esv.getEsvIndex()]);
            esv.setThetaAccel(theta_accel[esv.getEsvIndex()]);
        }
    }

    /**
     * getBiasDecomposition.
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
     * getBiasEnergy.
     *
     * @return a double.
     */
    public final double getBiasEnergy() {
        return getBiasEnergy(currentTemperature);
    }

    /**
     * Get ESV biases such as discretization, pH, etc. This method public and final for
     * error-checking; new ESVs should override biasEnergy().
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
            return biasEnergy(Constants.ROOM_TEMPERATURE);
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
     * Getter for the field <code>config</code>.
     *
     * @return a {@link ffx.potential.extended.ExtendedSystem.ExtendedSystemConfig} object.
     */
    public ExtendedSystemConfig getConfig() {
        return config;
    }

    /**
     * getConstantPh.
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
     * setConstantPh.
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
     * getDerivativeComponent.
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
                return vanDerWaals.getEsvDerivative(esvID);
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
                return particleMeshEwaldQI.getEsvDerivative(esvID);
            case Permanent:
                return particleMeshEwaldQI.getEsvDeriv_Permanent(esvID);
            case PermanentRealSpace:
                return particleMeshEwaldQI.getEsvDeriv_PermReal(esvID);
            case PermanentSelf:
                return particleMeshEwaldQI.getEsvDeriv_PermSelf(esvID);
            case PermanentReciprocal:
                return particleMeshEwaldQI.getEsvDeriv_PermRecip(esvID);
            case Induced:
                return particleMeshEwaldQI.getEsvDeriv_Induced(esvID);
            case InducedRealSpace:
                return particleMeshEwaldQI.getEsvDeriv_IndReal(esvID);
            case InducedSelf:
                return particleMeshEwaldQI.getEsvDeriv_IndSelf(esvID);
            case InducedReciprocal:
                return particleMeshEwaldQI.getEsvDeriv_IndRecip(esvID);
            default:
                throw new AssertionError(dd.name());
        }
    }

    /**
     * Potential gradient with respect to each ESV; used to propagate langevin dynamics.
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

    private double getDerivative(int esvID) {
        StringBuilder sb = new StringBuilder();
        final double temperature =
                (config.forceRoomTemp) ? Constants.ROOM_TEMPERATURE : currentTemperature;
        final boolean p = config.decomposeDeriv;
        ExtendedVariable esv = esvList.get(esvID);
        double esvDeriv = 0.0;
        final String format = " %-20.20s %2.2s %9.4f";
        if (config.biasTerm) {
            final double dBias = esv.getTotalBiasDeriv(temperature, false);
            if (p) {
                sb.append(format("  Biases: %9.4f", dBias));
            }
            final double dDiscr = esv.getDiscrBiasDeriv();
            if (p) {
                sb.append(format("    Discretizer: %9.4f", dDiscr));
            }
            if (esv instanceof TitrationESV) {
                final double dPh = ((TitrationESV) esv).getPhBiasDeriv(temperature);
                if (p) {
                    sb.append(format("    Acidostat: %9.4f", dPh));
                }
            }
            esvDeriv += dBias;
        }
        if (config.vanDerWaals) {
            final double dVdw = vanDerWaals.getEsvDerivative(esvID);
            if (p) {
                sb.append(format("  VanDerWaals: %9.4f", dVdw));
            }
            esvDeriv += dVdw;
        }
        if (config.electrostatics) {
            final double permanent = particleMeshEwaldQI.getEsvDeriv_Permanent(esvID);
            esvDeriv += permanent;
            if (p) {
                sb.append(format("  PermanentElec: %9.4f", permanent));
            }
            double permReal = particleMeshEwaldQI.getEsvDeriv_PermReal(esvID);
            double permSelf = particleMeshEwaldQI.getEsvDeriv_PermSelf(esvID);
            double permRecip = particleMeshEwaldQI.getEsvDeriv_PermRecip(esvID);
            if (p) {
                sb.append(format("    PermReal: %9.4f", permReal));
            }
            if (p) {
                sb.append(format("    PermRcpSelf: %9.4f", permSelf));
            }
            if (p) {
                sb.append(format("    PermRecipMpole: %9.4f", permRecip));
            }
            if (config.polarization) {
                final double induced = particleMeshEwaldQI.getEsvDeriv_Induced(esvID);
                esvDeriv += induced;
                if (p) {
                    sb.append(format("  Polarization: %9.4f", induced));
                }
                double indReal = particleMeshEwaldQI.getEsvDeriv_IndReal(esvID);
                double indSelf = particleMeshEwaldQI.getEsvDeriv_IndSelf(esvID);
                double indRecip = particleMeshEwaldQI.getEsvDeriv_IndRecip(esvID);
                if (p) {
                    sb.append(format("    IndReal: %9.4f", indReal));
                }
                if (p) {
                    sb.append(format("    IndSelf: %9.4f", indSelf));
                }
                if (p) {
                    sb.append(format("    IndRecip: %9.4f", indRecip));
                }
            }
        }
        if (config.bonded) {
            final double dBonded = esv.getBondedDeriv();
            if (p) {
                sb.append(format("  Bonded: %9.4f", dBonded));
            }
            esvDeriv += dBonded;
            /* If desired, decompose bonded contribution into component types from foreground and background. */
            if (config.decomposeBonded) {
                // Foreground portion:
                double fgSum = 0.0;
                HashMap<Class<? extends BondedTerm>, SharedDouble> fgMap = esv.getBondedDerivDecomp();
                for (SharedDouble dub : fgMap.values()) {
                    fgSum += dub.get();
                }
                if (p) {
                    sb.append(format("    Foreground: %9.4f", fgSum));
                }
                for (Class<? extends BondedTerm> clas : fgMap.keySet()) {
                    if (p) {
                        sb.append(
                                format(
                                        "      " + clas.getName().replaceAll("ffx.potential.bonded.", "") + ": " +
                                                fgMap.get(clas).get()));
                    }
                }
                // Background portion:
                double bgSum = 0.0;
                HashMap<Class<? extends BondedTerm>, SharedDouble> bgMap =
                        esv.getBackgroundBondedDerivDecomp();
                for (SharedDouble dub : bgMap.values()) {
                    bgSum += dub.get();
                }
                if (p) {
                    sb.append(format("    Background: %9.4f", bgSum));
                }
                for (Class<? extends BondedTerm> clas : bgMap.keySet()) {
                    if (p) {
                        sb.append(
                                format(
                                        "      " + clas.getName().replaceAll("ffx.potential.bonded.", "") + ": " +
                                                bgMap.get(clas).get()));
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
     * getEnergyComponent.
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
     * getEsv.
     *
     * @param esvID a int.
     * @return a {@link ffx.potential.extended.ExtendedVariable} object.
     */
    public ExtendedVariable getEsv(int esvID) {
        return esvList.get(esvID);
    }

    /**
     * getEsvForAtom.
     *
     * @param i a int.
     * @return a {@link ffx.potential.extended.ExtendedVariable} object.
     */
    public ExtendedVariable getEsvForAtom(int i) {
        return (isShared(i)) ? esvForShared[i] : (isUnshared(i)) ? esvForUnshared[i] : null;
    }

    /**
     * getEsvIndex.
     *
     * @param i a int.
     * @return a {@link java.lang.Integer} object.
     */
    public Integer getEsvIndex(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).esvIndex : null;
    }

    /**
     * Used only by ForceFieldEnergy and only once; we'd prefer to be rid of this altogether.
     * Background atoms are not true degrees of freedom.
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
     * getExtendedAndBackgroundMolecule.
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
     * All atoms of the fully-protonated system (not just those affected by this system).
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
     * getLambda.
     *
     * @param i a int.
     * @return a double.
     */
    public double getLambda(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getLambda() : Defaults.lambda;
    }

    /**
     * getLambda.
     *
     * @param esvIdChar a char.
     * @return a double.
     */
    public double getLambda(char esvIdChar) {
        return getEsv(esvIdChar - 'A').getLambda();
    }

    /**
     * getLambdaList.
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
     * getLambdaSwitch.
     *
     * @param i a int.
     * @return a double.
     */
    public double getLambdaSwitch(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getLambdaSwitch() : Defaults.lambdaSwitch;
    }

    /**
     * getNumESVs
     *
     * @return a int num of ESVs
     */
    public int getNumESVs() {
        return numESVs;
    }

    /**
     * getSwitchDeriv.
     *
     * @param i a int.
     * @return a double.
     */
    public double getSwitchDeriv(int i) {
        return (isExtended(i)) ? getEsvForAtom(i).getSwitchDeriv() : Defaults.switchDeriv;
    }

    /**
     * isAlphaScaled.
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isAlphaScaled(int i) {
        return isExtended(i) && isTitratableHydrogen(extendedAtoms[i]);
    }

    /**
     * Whether the Atom at extendedAtoms[i] is affected by an ESV of this system.
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isExtended(int i) {
        return isShared(i) || isUnshared(i);
    }

    /**
     * isShared.
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isShared(int i) {
        return esvForShared[i] != null;
    }

    /**
     * isUnshared.
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isUnshared(int i) {
        return esvForUnshared[i] != null;
    }

    /**
     * Processes lambda values based on propagation of theta value from Stochastic integrator in Molecular dynamics
     */
    public void preForce() {
        for (ExtendedVariable esv : esvList) {
            double sinTheta = sin(theta_position[esv.getEsvIndex()]);
            double oldLambda = esv.getLambda();
            esv.updateLambda(sinTheta * sinTheta, true);
            updateListeners();
            double newLambda = esv.getLambda();
            logger.info(format(" Propagating ESV[%d]: %g --> %g ", esv.esvIndex, oldLambda, newLambda));
        }
    }

    /**
     * Applies a chain rule term to the derivative to account for taking a derivative of lambda = sin(theta)^2
     *
     * @return
     */
    public double[] postForce() {
        double[] dEdL = ExtendedSystem.this.getDerivatives();
        for (ExtendedVariable esv : esvList) {
            //logger.info(format("dEdL: %g", dEdL[esv.getEsvIndex()]));
            dEdL[esv.getEsvIndex()] = dEdL[esv.getEsvIndex()] * sin(2 * theta_position[esv.getEsvIndex()]);
            //logger.info(format("dEdL*sin(2x): %g", dEdL[esv.getEsvIndex()]));
        }
        return dEdL;
    }

    /**
     * setLambda.
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
     * setLambda.
     *
     * @param esvIdChar a char.
     * @param lambda    a double.
     */
    public void setLambda(char esvIdChar, double lambda) {
        setLambda(esvIdChar - 'A', lambda);
    }

    /**
     * updateListeners.
     */
    private void updateListeners() {
        if (config.vanDerWaals) {
            vanDerWaals.updateEsvLambda();
        }
        if (config.electrostatics) {
            particleMeshEwaldQI.updateEsvLambda();
        }
    }

    /**
     * {@inheritDoc}
     *
     * <p>Allows simple iteration over ESV via "for (ExtendedVariable : ExtendedSystem)".
     */
    @Override
    public Iterator<ExtendedVariable> iterator() {
        return esvList.iterator();
    }

    /**
     * size.
     *
     * @return a int.
     */
    public int size() {
        return esvList.size();
    }

    int requestIndexing() {
        return indexer++;
    }

    /**
     * populate.
     *
     * @param residueIDs a {@link java.util.List} object.
     */
    public void populate(List<String> residueIDs) {
        // Locate the Residue identified by the given resid.
        Polymer[] polymers = molecularAssembly.getChains();
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

            MultiResidue titrating = TitrationUtils.titratingMultiresidueFactory(molecularAssembly, target);
            TitrationESV esv = new TitrationESV(this, titrating);
            this.addVariable(esv);
        }
    }

    /**
     * populate.
     *
     * @param residueIDs an array of {@link java.lang.String} objects.
     */
    public void populate(String[] residueIDs) {
        populate(Arrays.asList(residueIDs));
    }

    /**
     * populate.
     *
     * @param residueIDs a {@link java.lang.String} object.
     */
    public void populate(String residueIDs) {
        String[] tokens =
                (residueIDs.split(".").length > 1)
                        ? residueIDs.split(".")
                        : (residueIDs.split(",").length > 1)
                        ? residueIDs.split(",")
                        : new String[]{residueIDs};
        populate(tokens);
    }

    /**
     * setTemperature.
     *
     * @param set a double.
     */
    public void setTemperature(double set) {
        currentTemperature = set;
    }

    public static class ExtendedSystemConfig {
        public final boolean tautomer = prop("esv.tautomer", false);
        public final boolean bonded = prop("esv.bonded", false);
        public final boolean vanDerWaals = prop("esv.vanDerWaals", true);
        public final boolean electrostatics = prop("esv.electrostatics", true);
        public final boolean polarization = prop("esv.polarization", true);
        public final boolean biasTerm = prop("esv.biasTerm", true);
        public final boolean verbose = prop("esv.verbose", false);
        public final boolean decomposeBonded = prop("esv.decomposeBonded", false);
        public final boolean decomposeDeriv = prop("esv.decomposeDeriv", false);

        /**
         * Note that without the Lambda-Switch, the derivative dPol/dEsv is incorrect at L=0.0 and L=1.0
         */
        public final boolean allowLambdaSwitch = prop("esv.allowLambdaSwitch", true);
        public final boolean nonlinearMultipoles = prop("esv.nonlinearMultipoles", false); // sigmoid lambda Mpole switch
        public final double discrBias = prop("esv.biasMagnitude", 0.0);
        public final boolean forceRoomTemp = prop("esv.forceRoomTemp", false);
        public final boolean propagation = prop("esv.propagation", true);

        // Options below are untested and/or dangerous if changed.
        public final boolean cloneXyzIndices = prop("esv.cloneXyzIndices", true); // set bg_idx = fg_idx

        public static void print(ExtendedSystemConfig config) {
            List<Field> fields = Arrays.asList(ExtendedSystemConfig.class.getDeclaredFields());
            Collections.sort(
                    fields,
                    (Field t, Field t1) -> String.CASE_INSENSITIVE_ORDER.compare(t.getName(), t1.getName()));
            StringBuilder sb = new StringBuilder();
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
            sb.append(
                    format(
                            " %-30s %7.7s          %-30s %7.7s          %-30s %7.7s"
                                    + "\n %-30s %7.7s          %-30s %7.7s          %-30s %7.7s"
                                    + "\n %-30s %7.7s",
                            "polarization",
                            System.getProperty("polarization"),
                            "scf-algorithm",
                            System.getProperty("scf-algorithm"),
                            "polar-eps",
                            System.getProperty("polar-eps"),
                            "use-charges",
                            System.getProperty("use-charges"),
                            "use-dipoles",
                            System.getProperty("use-dipoles"),
                            "use-quadrupoles",
                            System.getProperty("use-quadrupoles"),
                            "grid-method",
                            System.getProperty("grid-method")));
            sb.append(format("\n"));
        }
    }

    /**
     * These populate the order-n preloaded lambda parameter arrays in VdW and PME in the absence of
     * an attached ESV.
     */
    private static final class Defaults {

        public static final ExtendedVariable esv = null;
        public static final Integer esvId = null;
        public static final double lambda = 1.0;
        public static final double lambdaSwitch = 1.0;
        public static final double switchDeriv = 1.0;

        private Defaults() {
        } // value singleton
    }
}
