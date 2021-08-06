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

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.PotentialComponent;
import ffx.potential.bonded.*;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.lang.reflect.Field;
import java.util.*;
import java.util.logging.Logger;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sin;

/**
 * ExtendedSystem class.
 *
 * @author Stephen LuCore
 * @since 1.0
 */
public class NewExtendedSystem implements Iterable<NewExtendedVariable> {

    /**
     * Stores all the default values listed in prop() calls below.
     */
    public static final NewExtendedSystemConfig DefaultConfig;

    private static final Logger logger = Logger.getLogger(NewExtendedSystem.class.getName());

    public static final double THETA_MASS = 1.0; //Atomic Mass Units

    public static final double THETA_FRICTION = 5.0; // psec^-1

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
            DefaultConfig = new NewExtendedSystemConfig();
            System.getProperties().putAll(bak);
        }
    }

    /**
     * Stores configuration of system properties at instantiation of ExtendedSystem.
     */
    public final NewExtendedSystemConfig config;
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
    //TODO: Create data structure to map titration and tautomer ESVs back to an atom.
    /**
     * Extended Variable list for tirating and tautomerizing atoms.
     * 1st dimension refers to Atom that has the ESV
     * 2nd dimension refers to ESV type: [0] = titrationESV, [1]=tautomerESV
     */
    private final NewExtendedVariable[][] esvsForAtom;
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
    private List<NewExtendedVariable> esvList;
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
     * @param mola a {@link MolecularAssembly} object.
     */
    public NewExtendedSystem(MolecularAssembly mola) {
        this(mola, new NewExtendedSystemConfig());
    }

    /**
     * Construct extended system with the provided configuration.
     *
     * @param mola   a {@link MolecularAssembly} object.
     * @param config a {@link NewExtendedSystem.NewExtendedSystemConfig} object.
     */
    public NewExtendedSystem(MolecularAssembly mola, NewExtendedSystemConfig config) {
        this.config = config;
        this.molecularAssembly = mola;

        ForceFieldEnergy forceFieldEnergy = mola.getPotentialEnergy();
        if (forceFieldEnergy == null) {
            logger.severe("No potential energy found?");
        }

        CompositeConfiguration properties = molecularAssembly.getProperties();
        thetaFriction = properties.getDouble("esv.friction", NewExtendedSystem.THETA_FRICTION);

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
        esvsForAtom = new NewExtendedVariable[2][nAtomsExt];
        if (config.verbose) {
            NewExtendedSystemConfig.print(config);
        }
    }

    /**
     * Prefer ExtendedSystem::populate to manual ESV creation.
     *
     * @param esv a {@link NewExtendedVariable} object.
     */
    public void addVariable(NewExtendedVariable esv) {
        esvSystemActive = true;
        if (esvList == null) {
            esvList = new ArrayList<>();
        }
        if (esvList.contains(esv)) {
            logger.warning(format("Attempted to add duplicate variable %s to system.", esv));
            return;
        }
        esvList.add(esv);

        numESVs = esvList.size();
        if (esv instanceof NewTitrationESV) {
            if (constantSystemPh == null) {
                logger.severe("Set ExtendedSystem (constant) pH before adding TitrationESVs.");
            }
        }
        for (int i = 0; i < nAtomsExt; i++) {
            if(esv instanceof NewTitrationESV){
                if(esv.viewExtendedAtoms().contains(extendedAtoms[i])){
                    esvsForAtom[i][0] = esv;
                }
            }
            if(esv instanceof NewTautomerESV){
                if(esv.viewExtendedAtoms().contains(extendedAtoms[i])){
                    esvsForAtom[i][1] = esv;
                }
            }
        }
        updateListeners();
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
        for (NewExtendedVariable esv : esvList) {
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
        for (NewExtendedVariable esv : esvList) {
            esv.setTheta(theta_position[esv.getEsvIndex()]);
            esv.setThetaVelocity(theta_velocity[esv.getEsvIndex()]);
            esv.setThetaAccel(theta_accel[esv.getEsvIndex()]);
        }
    }

    /**
     * getBiasDecomposition.
     *
     * @return a {@link String} object.
     */
    public String getBiasDecomposition() {
        if (!config.biasTerm) {
            return "";
        }
        double discrBias = 0.0;
        double phBias = 0.0;
        for (NewExtendedVariable esv : esvList) {
            discrBias += esv.getDiscrBias();
            if (esv instanceof NewTitrationESV) {
                phBias += ((NewTitrationESV) esv).getPhBias(currentTemperature);
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
        for (NewExtendedVariable esv : this) {
            biasEnergySum += esv.getTotalBias(temperature, false);
        }
        return biasEnergySum;
    }

    /**
     * Getter for the field <code>config</code>.
     *
     * @return a {@link NewExtendedSystem.NewExtendedSystemConfig} object.
     */
    public NewExtendedSystemConfig getConfig() {
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
        return constantSystemPh;
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
     * @param dd    a {@link PotentialComponent} object.
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
            case Bias:
            case pHMD:
                return esvList.get(esvID).getTotalBiasDeriv(currentTemperature, false);
            case Discretizer:
                return esvList.get(esvID).getDiscrBiasDeriv();
            case Acidostat:
                return ((NewTitrationESV) esvList.get(esvID)).getPhBiasDeriv(currentTemperature);
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
        double[] esvDeriv = new double[numESVs];
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
        NewExtendedVariable esv = esvList.get(esvID);
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
            if (esv instanceof NewTitrationESV) {
                final double dPh = ((NewTitrationESV) esv).getPhBiasDeriv(temperature);
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
     * @param component a {@link PotentialComponent} object.
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
                    uComp += ((NewTitrationESV) esvList.get(i)).getPhBias(currentTemperature);
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
     * @return a {@link NewExtendedVariable} object.
     */
    public NewExtendedVariable getEsv(int esvID) {
        return esvList.get(esvID);
    }

    /**
     * getEsvForAtom.
     *
     * @param i a int.
     * @return a {@link NewExtendedVariable} object.
     */
    public NewExtendedVariable[] getEsvForAtom(int i) {
        NewExtendedVariable[] newExtendedVariable = new NewExtendedVariable[2];
        newExtendedVariable[0] = esvsForAtom[i][0];
        newExtendedVariable[1] = esvsForAtom[i][1];

        return newExtendedVariable;
    }

    /**
     * getEsvIndex.
     *
     * @param i a int.
     * @return a {@link Integer} object.
     */
    public Integer[] getEsvIndex(int i) {
        Integer[] indexes = new Integer[2];
        for(int j=0; j < 2; j++){
            if(esvsForAtom[i][j] != null){
                indexes[j] = esvsForAtom[i][j].esvIndex;
            }
            else{
                indexes[j] = null;
            }
        }
        return indexes;
    }

    /**
     * All atoms of the fully-protonated system (not just those affected by this system).
     *
     * @return an array of {@link Atom} objects.
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
    public double[] getLambda(int i) {
        double[] lambdas = new double[2];
        for(int j=0; j < 2; j++){
            if(esvsForAtom[i][j] != null){
                lambdas[j] = esvsForAtom[i][j].lambda;
            }
            else{
                lambdas[j] = Defaults.lambda;
            }
        }
        return lambdas;
    }

    /**
     * getLambdaList.
     *
     * @return a {@link String} object.
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
    public double[] getLambdaSwitch(int i) {
        double[] lambdaSwitch = new double[2];
        for(int j=0; j < 2; j++){
            if(esvsForAtom[i][j] != null){
                lambdaSwitch[j] = esvsForAtom[i][j].getLambdaSwitch();
            }
            else{
                lambdaSwitch[j] = Defaults.lambdaSwitch;
            }
        }
        return lambdaSwitch;
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
    public double[] getSwitchDeriv(int i) {
        double[] lambdaSwitchDeriv = new double[2];
        for(int j=0; j < 2; j++){
            if(esvsForAtom[i][j] != null){
                lambdaSwitchDeriv[j] = esvsForAtom[i][j].getSwitchDeriv();
            }
            else{
                lambdaSwitchDeriv[j] = Defaults.switchDeriv;
            }
        }
        return lambdaSwitchDeriv;
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
        return (esvsForAtom[i][0] != null || esvsForAtom[i][1] != null);
    }

    /**
     * Processes lambda values based on propagation of theta value from Stochastic integrator in Molecular dynamics
     */
    public void preForce() {
        for (NewExtendedVariable esv : esvList) {
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
     * @return dE/dL a double[]
     */
    public double[] postForce() {
        double[] dEdL = NewExtendedSystem.this.getDerivatives();
        for (NewExtendedVariable esv : esvList) {
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
    public Iterator<NewExtendedVariable> iterator() {
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
     * setTemperature.
     *
     * @param set a double.
     */
    public void setTemperature(double set) {
        currentTemperature = set;
    }

    public static class NewExtendedSystemConfig {
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

        public static void print(NewExtendedSystemConfig config) {
            List<Field> fields = Arrays.asList(NewExtendedSystemConfig.class.getDeclaredFields());
            fields.sort((Field t, Field t1) -> String.CASE_INSENSITIVE_ORDER.compare(t.getName(), t1.getName()));
            StringBuilder sb = new StringBuilder();
            for (int i = 0, col = 0; i < fields.size(); i++) {
                if (++col > 3) {
                    sb.append("\n");
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
            sb.append("\n");
            logger.info(sb.toString());
        }
    }

    /**
     * These populate the order-n preloaded lambda parameter arrays in VdW and PME in the absence of
     * an attached ESV.
     */
    private static final class Defaults {

        public static final NewExtendedVariable esv = null;
        public static final Integer esvId = null;
        public static final double lambda = 1.0;
        public static final double lambdaSwitch = 1.0;
        public static final double switchDeriv = 1.0;

        private Defaults() {
        } // value singleton
    }
}
