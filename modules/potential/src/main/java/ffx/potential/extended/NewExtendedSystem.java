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
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.lang.reflect.Field;
import java.util.*;
import java.util.logging.Logger;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;
import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

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

    public static final double DISCR_BIAS = 1.0; // kcal/mol

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
     * Array of booleans that is initialized to match the number of atoms in the molecular assembly
     * noting whether the atom is extended.
     */
    private final boolean[] isExtended;
    /**
     * System PH.
     */
    private double constantSystemPh = 7.4;
    /**
     * Number of ESV.
     */
    private int numESVs;

    /**
     * Target system temperature.
     */
    private double currentTemperature = Constants.ROOM_TEMPERATURE;

    private List<Residue> titratingResidueList;
    private List<Residue> tautomerizingResidueList;
    private List<Residue> extendedResidueList;

    private double initialTitrationLambda = 1.0;
    private double initialTautomerLambda = 1.0;
    private double[] lambdaArray;

    /**
     * Current value of theta for each ESV.
     */
    public double[] thetaPosition;
    /**
     * Current theta velocity for each ESV.
     */
    public double[] thetaVelocity;
    /**
     * Current theta acceleration for each ESV.
     */
    public double[] thetaAccel;
    /**
     * Mass of each theta particle. Different theta mass for each particle are not supported.
     */
    public double[] thetaMassArray;
    /**
     * The system defined theta mass of the fictional particle. Used to fill theta mass array.
     */
    private double thetaMass;
    /**
     * Friction for the ESV system
     */
    public final double thetaFriction;

    private final double discrBiasMag;

    //TODO: Loop over protonation ESVs to sum the discr, acidostat, and Fmod energy terms.
    //TODO: Collect partial derivs for each term. Keep all in one place.

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
        thetaMass = properties.getDouble("esv.mass", NewExtendedSystem.THETA_MASS);
        discrBiasMag = properties.getDouble("discretize.bias", DISCR_BIAS);
        vanDerWaals = forceFieldEnergy.getVdwNode();
        particleMeshEwaldQI = forceFieldEnergy.getPmeQiNode();

        titratingResidueList = new ArrayList<>();
        tautomerizingResidueList = new ArrayList<>();
        extendedResidueList = new ArrayList<>();
        // Initialize atom arrays with the existing assembly.
        Atom[] atoms = mola.getAtomArray();
        isExtended = new boolean[atoms.length];
    }

    /**
     * Prefer ExtendedSystem::populate to manual ESV creation.
     *
     * @param residue a {@link Residue} object.
     */
    public void addVariable(Residue residue) {
        esvSystemActive = true;
        if(isTitrable(residue)){
            titratingResidueList.add(residue);
            if(isTautomer(residue)) {
                tautomerizingResidueList.add(residue);
            }
            extendAtoms(residue);
            updateListeners();
        }
        else{
            logger.warning(format("The residue %s is not supported for titration", residue.getAminoAcid3().name()));
            return;
        }


    }

    private void extendAtoms(Residue residue){
        List<Atom> currentAtomList = residue.getAtomList();
        Atom[] atoms = molecularAssembly.getAtomArray();
        for (int i = 0; i < atoms.length; i++) {
            if (currentAtomList.contains(atoms[i])) {
                this.isExtended[i] = true;
            }
        }
    }

    public void createMDThetaArrays() {
        extendedResidueList = titratingResidueList;
        extendedResidueList.addAll(tautomerizingResidueList);
        int size = extendedResidueList.size();

        lambdaArray = new double[size];
        thetaPosition = new double[size];
        thetaVelocity = new double[size];
        thetaAccel = new double[size];
        thetaMassArray = new double[size];

        //Theta masses should always be the same for each ESV
        Arrays.fill(thetaMassArray, thetaMass);

        for (int i=0; i< extendedResidueList.size(); i++){
            if(i < titratingResidueList.size()){
                lambdaArray[i] = initialTitrationLambda;
                thetaPosition[i] = Math.asin(Math.sqrt(initialTitrationLambda));
                Random random = new Random();
                thetaVelocity[i] = random.nextGaussian() * sqrt(kB * 298.15 / thetaMass);
                double dUdTheta = getDerivative(i)* sin(2*thetaPosition[i]);
                thetaAccel[i] = -Constants.KCAL_TO_GRAM_ANG2_PER_PS2 * dUdTheta / thetaMass;
            }
            else{
                lambdaArray[i] = initialTautomerLambda;
                thetaPosition[i] = Math.asin(Math.sqrt(initialTautomerLambda));
                Random random = new Random();
                thetaVelocity[i] = random.nextGaussian() * sqrt(kB * 298.15 / thetaMass);
                double dUdTheta = getDerivative(i)* sin(2*thetaPosition[i]);
                thetaAccel[i] = -Constants.KCAL_TO_GRAM_ANG2_PER_PS2 * dUdTheta / thetaMass;
            }
        }
    }

    public boolean isExtended(int i) {
        return isExtended[i];
    }

    public boolean isExtended(Residue residue){
        return extendedResidueList.contains(residue);
    }

    private Residue getResidueFromAtom(int atomIndex){
        Atom atom = molecularAssembly.getAtomList().get(atomIndex);
        int residueIndex = atom.getResidueNumber();
        return molecularAssembly.getResidueList().get(residueIndex);
    }

    private boolean isTitrable(Residue residue){
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        if(AA3 == AminoAcidUtils.AminoAcid3.ASP || AA3 == AminoAcidUtils.AminoAcid3.ASD || AA3 == AminoAcidUtils.AminoAcid3.ASH ||
                AA3 == AminoAcidUtils.AminoAcid3.GLU || AA3 == AminoAcidUtils.AminoAcid3.GLD || AA3 == AminoAcidUtils.AminoAcid3.GLH ||
                AA3 == AminoAcidUtils.AminoAcid3.HIS || AA3 == AminoAcidUtils.AminoAcid3.HID || AA3 == AminoAcidUtils.AminoAcid3.HIE ||
                AA3 == AminoAcidUtils.AminoAcid3.LYS || AA3 == AminoAcidUtils.AminoAcid3.LYD){
            return true;
        }
        else{
            return false;
        }
    }

    private boolean isTautomer(Residue residue){
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        if(AA3 == AminoAcidUtils.AminoAcid3.ASP || AA3 == AminoAcidUtils.AminoAcid3.ASD || AA3 == AminoAcidUtils.AminoAcid3.ASH ||
            AA3 == AminoAcidUtils.AminoAcid3.GLU || AA3 == AminoAcidUtils.AminoAcid3.GLD || AA3 == AminoAcidUtils.AminoAcid3.GLH ||
            AA3 == AminoAcidUtils.AminoAcid3.HIS || AA3 == AminoAcidUtils.AminoAcid3.HID || AA3 == AminoAcidUtils.AminoAcid3.HIE){
            return true;
        }
        else{
            return false;
        }
    }

    public void setTitrationLambda(int resIndex, double lambda){
        Residue residue = molecularAssembly.getResidueList().get(resIndex);
        if(titratingResidueList.contains(residue)){
            int index = titratingResidueList.indexOf(residue);
            thetaPosition[index] = Math.asin(Math.sqrt(lambda));
            lambdaArray[index] = lambda;
            updateListeners();
        }
        else{
            logger.warning(format("This residue %s is not titrating.", residue.getName()));
        }
    }

    public void setTautomerLambda(int resIndex, double lambda){
        Residue residue = molecularAssembly.getResidueList().get(resIndex);
        if(tautomerizingResidueList.contains(residue)){
            // The correct index in the theta arrays for tautomer coordinates is after the titration list.
            // So titrationList.size() + tautomerIndex should match with appropriate spot in thetaPosition, etc.
            int index = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
            thetaPosition[index] = Math.asin(Math.sqrt(lambda));
            lambdaArray[index] = lambda;
            updateListeners();
        }
        else{
            logger.warning(format("This residue %s does not have any titrating tautomers.", residue.getName()));
        }
    }

    private void updateLambdas(){
        for(int i = 0; i < extendedResidueList.size(); i++){
            double sinTheta = Math.sin(thetaPosition[i]);
            lambdaArray[i] = sinTheta * sinTheta;
        }
        updateListeners();
    }

    public double getTitrationLambda(int i) {
        if(isExtended(i)){
            Residue residue = getResidueFromAtom(i);
            int index = titratingResidueList.indexOf(residue);
            return lambdaArray[index];
        }
        else{
            return Defaults.lambda;
        }
    }

    public double getTautomerLambda(int i){
        if(isExtended(i)){
            Residue residue = getResidueFromAtom(i);
            if(tautomerizingResidueList.contains(residue)){
                int index = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
                return lambdaArray[index];
            }
            else{
                return Defaults.lambda;
            }
        }
        else{
            return Defaults.lambda;
        }
    }

    private double getDiscrBias(){
        double discrBias = 0.0;
        for(int i=0; i <extendedResidueList.size();i++){
            double lambda = lambdaArray[i];
            double bias = - 4.0 * discrBiasMag * (lambda - 0.5) * (lambda - 0.5);
            discrBias += bias;
        }
        return discrBias;
    }

    private double[] getDiscrBiasDeriv(){
        double[] discrBiasDeriv = new double[extendedResidueList.size()];
        for(int i=0; i <extendedResidueList.size(); i++){
            double lambda = lambdaArray[i];
            discrBiasDeriv[i] = -8.0 * discrBiasMag * (lambda - 0.5);
        }
        return discrBiasDeriv;
    }

    private double getPhBias(){
        double pHBias = 0.0;
        for(Residue residue : titratingResidueList){
            int titrationResIndex = titratingResidueList.indexOf(residue);
            double titrationLambda = lambdaArray[titrationResIndex];
            double pKa1 = ConstantPhUtils.
            if(isTautomer(residue)){
                int tautomerResIndex = tautomerizingResidueList.indexOf(residue);
                double tautomerLambda = lambdaArray[tautomerResIndex+titratingResidueList.size()];
            }

        }
    }

    /**
     * getBiasDecomposition.
     *
     * @return a {@link String} object.
     */
    /*public String getBiasDecomposition() {
        if (!config.biasTerm) {
            return "";
        }
        double discrBias = 0.0;
        double phBias = 0.0;
        for (Residue residue : extendedResidueList) {
            discrBias += esv.getDiscrBias();
            if (esv instanceof NewTitrationESV) {
                phBias += ((NewTitrationESV) esv).getPhBias(currentTemperature);
            }
        }
        return format("    %-16s %16.8f\n", "Discretizer", discrBias)
                + format("    %-16s %16.8f\n", "Acidostat", phBias);
    }*/

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
        return constantSystemPh;
    }

    /**
     * setConstantPh.
     *
     * @param pH a double.
     */
    public void setConstantPh(double pH) {
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
        return molecularAssembly.getAtomArray();
    }

    /**
     * Companion to getExtendedAtoms() for vdw::setAtoms and pme::setAtoms.
     *
     * @return an array of {@link int} objects.
     */
    public int[] getExtendedMolecule() {
        return molecularAssembly.getMoleculeNumbers();
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
        return isExtended(i) && isTitratableHydrogen(molecularAssembly.getAtomArray()[i]);
    }

    /**
     * Processes lambda values based on propagation of theta value from Stochastic integrator in Molecular dynamics
     */
    public void preForce() {
        for (NewExtendedVariable esv : esvList) {
            double sinTheta = sin(thetaPosition[esv.getEsvIndex()]);
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
        double[] dEdTheta = new double[dEdL.length];
        for (NewExtendedVariable esv : esvList) {
            int esvIndex = esv.getEsvIndex();
            dEdTheta[esvIndex] = dEdL[esvIndex] * sin(2 * thetaPosition[esvIndex]);
        }
        return dEdTheta;
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
