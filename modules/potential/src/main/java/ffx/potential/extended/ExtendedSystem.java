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
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TitrationUtils;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.*;
import java.util.logging.Logger;

import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.*;

import edu.rit.pj.reduction.SharedDouble;

/**
 * ExtendedSystem class.
 *
 * @author Andrew Thiel
 * @since 1.0
 */
public class ExtendedSystem {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());

    private static final double THETA_MASS = 1.0; //Atomic Mass Units

    private static final double THETA_FRICTION = 5.0; // psec^-1

    private static final double DISCR_BIAS = 1.0; // kcal/mol

    private static final double LOG10 = log(10.0);
    private static final int discrBiasIndex = 0;
    private static final int pHBiasIndex = 1;
    private static final int modelBiasIndex = 2;
    private static final int dDiscr_dTitrIndex = 3;
    private static final int dPh_dTitrIndex = 4;
    private static final int dModel_dTitrIndex = 5;
    private static final int dDiscr_dTautIndex = 6;
    private static final int dPh_dTautIndex = 7;
    private static final int dModel_dTautIndex = 8;
    /**
     * Array of ints that is initialized to match the number of atoms in the molecular assembly.
     * 1 indicates that the tautomer lambda direction is normal.
     * -1 indicates that the tautomer lambda direction is reversed (1-x).
     * 0 indicates that the atom is not a tautomerizing atom.
     */
    public final int[] tautomerDirections;
    /**
     * Array of doubles that is initialized to match the number of atoms in the molecular assembly.
     * Only elements that match a titrating atom will have their lambda updated.
     */
    private final double[] titrationLambdas;
    /**
     * Array of doubles that is initialized to match the number of atoms in the molecular assembly.
     * Only elements that match a tautomerizing atom will have their lambda updated.
     */
    private final double[] tautomerLambdas;
    /**
     * Shared double that is initialized to match the number of ESVs in the system.
     * Once reduced, will equal either dU_Titr/dLambda or dU_Taut/dLambda for specific ESV
     */
    private final SharedDouble[] esvVdwDerivs;

    private final SharedDouble[] esvPermElecDerivs;
    private final SharedDouble[] esvIndElecDerivs;
    /**
     * Array of AminoAcid3 initialized  to match the number of atoms in the system.
     * Used to know how to apply vdW or electrostatic ESV terms for the atom.
     */
    public final AminoAcid3[] residueNames;
    /**
     * MolecularAssembly instance.
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * Titration Utils instance. This instance is the master copy that will be distributed to electrostatics classes
     * when an Extended System is attached.
     */
    private final TitrationUtils titrationUtils;
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
     * noting whether the atom is titrating. Note that any side chain atom that belongs to a titrating residue
     * will be flagged as titrating for purposes of scaling electrostatic parameters.
     */
    private final boolean[] isTitrating;
    /**
     * Array of booleans that is initialized to match the number of atoms in the assembly to note whether an atom is
     * specifically a titrating hydrogen.
     */
    private final boolean[] isTitratingHydrogen;
    /**
     * Array of booleans that is initialized to match the number of atoms in the molecular assembly
     * noting whether the atom is tautomerizing. Note that any side chain atom that belongs to a tautomerizing residue
     * will be flagged as tautomerizing for purposes of scaling electrostatic parameters.
     */
    private final boolean[] isTautomerizing;
    /**
     * Array of ints that is initialized to match the number of atoms in the molecular assembly.
     * Elements correspond to residue index in the titratingResidueList. Only set for titrating residues, -1 otherwise.
     */
    private final int[] titrationIndexMap;
    /**
     * Array of ints that is initialized to match the number of atoms in the molecular assembly.
     * Elements correspond to residue index in the tautomerizingResidueList. Only set for tautomerizing residues, -1 otherwise.
     */
    private final int[] tautomerIndexMap;
    /**
     * List of titrating residues.
     */
    private final List<Residue> titratingResidueList;
    /**
     * List of tautomerizing residues.
     */
    private final List<Residue> tautomerizingResidueList;
    /**
     * Concatenated list of titrating residues + tautomerizing residues.
     */
    private final List<Residue> extendedResidueList;
    /**
     * The system defined theta mass of the fictional particle. Used to fill theta mass array.
     */
    private final double thetaMass;
    /**
     * Friction for the ESV system
     */
    private final double thetaFriction;
    /**
     * Array of lambda values that matches residues in extendedResidueList
     */
    private final double[] extendedLambdas;
    /**
     * Descritizer Bias Magnitude. Default is 1 kcal/mol.
     */
    private final double discrBiasMag;
    // Controls for turning of certain terms for testing.
    private final boolean doVDW;
    private final boolean doElectrostatics;
    private final boolean doBias;
    private final boolean doPolarization;
    /**
     * System PH.
     */
    private double constantSystemPh = 7.4;
    /**
     * Target system temperature.
     */
    private double currentTemperature = Constants.ROOM_TEMPERATURE;
    /**
     * Current value of theta for each ESV.
     */
    private double[] thetaPosition;
    /**
     * Current theta velocity for each ESV.
     */
    private double[] thetaVelocity;
    /**
     * Current theta acceleration for each ESV.
     */
    private double[] thetaAccel;
    /**
     * Mass of each theta particle. Different theta mass for each particle are not supported.
     */
    private double[] thetaMassArray;

    /**
     * Construct extended system with the provided configuration.
     *
     * @param mola a {@link MolecularAssembly} object.
     */
    public ExtendedSystem(MolecularAssembly mola) {
        this.molecularAssembly = mola;

        ForceField forceField = mola.getForceField();
        ForceFieldEnergy forceFieldEnergy = mola.getPotentialEnergy();
        if (forceFieldEnergy == null) {
            logger.severe("No potential energy found?");
        }

        CompositeConfiguration properties = molecularAssembly.getProperties();
        titrationUtils = new TitrationUtils(forceField);
        thetaFriction = properties.getDouble("esv.friction", ExtendedSystem.THETA_FRICTION);
        thetaMass = properties.getDouble("esv.mass", ExtendedSystem.THETA_MASS);
        discrBiasMag = properties.getDouble("discretize.bias.magnitude", DISCR_BIAS);
        double initialTitrationLambda = properties.getDouble("lambda.titration.initial", 1.0);
        double initialTautomerLambda = properties.getDouble("lambda.tautomer.initial", 1.0);
//        boolean bonded = properties.getBoolean("esv.bonded", false);
        doVDW = properties.getBoolean("esv.vdW", true);
        doElectrostatics = properties.getBoolean("esv.elec", true);
        doBias = properties.getBoolean("esv.bias", true);
        doPolarization = properties.getBoolean("esv.polarization", true);
//        boolean verbose = properties.getBoolean("esv.verbose", false);
//        boolean decomposeBonded = properties.getBoolean("esv.decomposeBonded", false);
//        boolean decomposeDeriv = properties.getBoolean("esv.decomposeDeriv", false);
//        boolean nonlinearMultipoles = properties.getBoolean("esv.nonlinearMultipoles", false); // sigmoid lambda Mpole switch
//        boolean forceRoomTemp = properties.getBoolean("esv.forceRoomTemp", false);
//        boolean propagation = properties.getBoolean("esv.propagation", true);

        vanDerWaals = forceFieldEnergy.getVdwNode();
        particleMeshEwaldQI = forceFieldEnergy.getPmeQiNode();

        titratingResidueList = new ArrayList<>();
        tautomerizingResidueList = new ArrayList<>();
        extendedResidueList = new ArrayList<>();
        // Initialize atom arrays with the existing assembly.
        Atom[] atoms = mola.getAtomArray();

        isTitrating = new boolean[atoms.length];
        isTitratingHydrogen = new boolean[atoms.length];
        isTautomerizing = new boolean[atoms.length];
        titrationLambdas = new double[atoms.length];
        tautomerLambdas = new double[atoms.length];
        titrationIndexMap = new int[atoms.length];
        tautomerIndexMap = new int[atoms.length];
        tautomerDirections = new int[atoms.length];
        residueNames = new AminoAcid3[atoms.length];

        Arrays.fill(isTitrating, false);
        Arrays.fill(isTitratingHydrogen, false);
        Arrays.fill(isTautomerizing, false);
        Arrays.fill(titrationLambdas, 1.0);
        Arrays.fill(tautomerLambdas, 0.0);
        Arrays.fill(titrationIndexMap, -1);
        Arrays.fill(tautomerIndexMap, -1);
        Arrays.fill(tautomerDirections, 0);
        Arrays.fill(residueNames, AminoAcid3.UNK);

        // Cycle through each residue to determine if it is titratable or tautomerizing.
        // If a residue is one of these, add to titration or tautomer lists.
        // Next, loop through all atoms and check to see if the atom belongs to this residue.
        // If the atom does belong to this residue, set all corresponding variables in the respective titration or tautomer array (size = numAtoms).
        // Store the index of the residue in the respective list into a map array (size = numAtoms).
        List<Residue> residueList = molecularAssembly.getResidueList();
        for (Residue residue : residueList) {
            if (isTitrable(residue)) {
                titratingResidueList.add(residue);
                List<Atom> atomList = residue.getSideChainAtoms();
                for (Atom atom : atomList) {
                    int atomIndex = atom.getIndex();
                    residueNames[atomIndex] = residue.getAminoAcid3();
                    isTitrating[atomIndex] = true;
                    titrationLambdas[atomIndex] = initialTitrationLambda;
                    int titrationIndex = titratingResidueList.indexOf(residue);
                    titrationIndexMap[atomIndex] = titrationIndex;
                    isTitratingHydrogen[atomIndex] = TitrationUtils.isTitratingHydrogen(residue.getAminoAcid3(), atom);
                }
                // If is a tautomer, it must also be titrating.
                if (isTautomer(residue)) {
                    tautomerizingResidueList.add(residue);
                    for (Atom atom : atomList) {
                        int atomIndex = atom.getIndex();
                        isTautomerizing[atomIndex] = true;
                        tautomerLambdas[atomIndex] = initialTautomerLambda;
                        int tautomerIndex = tautomerizingResidueList.indexOf(residue);
                        tautomerIndexMap[atomIndex] = tautomerIndex;
                        tautomerDirections[atomIndex] = TitrationUtils.getTitratingHydrogenDirection(residue.getAminoAcid3(), atom);
                    }
                }
            }
        }

        //Concatenate titratingResidueList and tautomerizingResidueList
        extendedResidueList.addAll(titratingResidueList);
        extendedResidueList.addAll(tautomerizingResidueList);
        //Arrays that are sent to integrator are based on extendedResidueList size
        int size = extendedResidueList.size();
        extendedLambdas = new double[size];
        thetaPosition = new double[size];
        thetaVelocity = new double[size];
        thetaAccel = new double[size];
        thetaMassArray = new double[size];
        esvVdwDerivs = new SharedDouble[size];
        esvPermElecDerivs = new SharedDouble[size];
        esvIndElecDerivs =  new SharedDouble[size];
        for(int i=0; i < size; i++){
            esvVdwDerivs[i] = new SharedDouble(0.0);
            esvPermElecDerivs[i] = new SharedDouble(0.0);
            esvIndElecDerivs[i] = new SharedDouble(0.0);
        }

        //Theta masses should always be the same for each ESV
        Arrays.fill(thetaMassArray, thetaMass);

        for (int i = 0; i < extendedResidueList.size(); i++) {
            if (i < titratingResidueList.size()) {
                initializeThetaArrays(i, initialTitrationLambda);
            } else {
                initializeThetaArrays(i, initialTautomerLambda);
            }
        }
    }

    /**
     * During constructor, initialize arrays that will hold theta positions, velocities, and accelerations.
     * Positions determined from starting lambda.
     * Velocities randomly set according to Maxwell Boltzmann Distribution based on temperature.
     * Accelerations determined from initial forces.
     * @param index index of ExtendedResidueList for which to set these values.
     * @param lambda starting lambda value for each ESV.
     */
    private void initializeThetaArrays(int index, double lambda) {
        extendedLambdas[index] = lambda;
        thetaPosition[index] = Math.asin(Math.sqrt(lambda));
        Random random = new Random();
        thetaVelocity[index] = random.nextGaussian() * sqrt(kB * 298.15 / thetaMass);
        double dUdL = 0.0;//getDerivatives()[index];
        double dUdTheta = dUdL * sin(2 * thetaPosition[index]);
        thetaAccel[index] = -Constants.KCAL_TO_GRAM_ANG2_PER_PS2 * dUdTheta / thetaMass;
    }

    public void initEsvVdw(){
        for (int i = 0; i < extendedResidueList.size(); i++) {
            esvVdwDerivs[i].set(0.0);
        }
    }

    public void initEsvPermElec(){
        for (int i = 0; i < extendedResidueList.size(); i++) {
            esvPermElecDerivs[i].set(0.0);
        }
    }

    public void initEsvIndElec(){
        for (int i = 0; i < extendedResidueList.size(); i++) {
            esvIndElecDerivs[i].set(0.0);
        }
    }

    public boolean isTitrating(int atomIndex) {
        return isTitrating[atomIndex];
    }

    public boolean isTitratingHydrogen(int atomIndex) {
        return isTitratingHydrogen[atomIndex];
    }

    public boolean isTautomerizing(int atomIndex) {
        return isTautomerizing[atomIndex];
    }

    public boolean isExtended(Residue residue) {
        return extendedResidueList.contains(residue);
    }

    public boolean isTitrable(Residue residue) {
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        return AA3.isConstantPhTitratable;
    }

    public boolean isTautomer(Residue residue) {
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        return AA3.isConstantPhTautomer;
    }

    public double getTitrationLambda(Residue residue) {
        if (titratingResidueList.contains(residue)) {
            int resIndex = titratingResidueList.indexOf(residue);
            return extendedLambdas[resIndex];
        } else {
            return 1.0;
        }
    }

    public double getTitrationLambda(int atomIndex){
        return titrationLambdas[atomIndex];
    }

    public int getTitrationESVIndex(int i){
        return titrationIndexMap[i];
    }

    public double getTautomerLambda(Residue residue) {
        if (tautomerizingResidueList.contains(residue)) {
            int resIndex = tautomerizingResidueList.indexOf(residue);
            return extendedLambdas[titratingResidueList.size() + resIndex];
        } else {
            return 1.0;
        }
    }

    public double getTautomerLambda(int atomIndex){
        return tautomerLambdas[atomIndex];
    }

    public int getTautomerESVIndex(int i){
        return tautomerIndexMap[i];
    }

    public void setTitrationLambda(Residue residue, double lambda) {
        if (titratingResidueList.contains(residue)) {
            int index = titratingResidueList.indexOf(residue);
            extendedLambdas[index] = lambda;
            thetaPosition[index] = Math.asin(Math.sqrt(lambda));
            List<Atom> currentAtomList = residue.getSideChainAtoms();
            for (Atom atom : currentAtomList) {
                int atomIndex = atom.getIndex();
                titrationLambdas[atomIndex] = lambda;
            }
        } else {
            logger.warning(format("This residue %s is not titrating.", residue.getName()));
        }
    }

    public void setTautomerLambda(Residue residue, double lambda) {
        if (tautomerizingResidueList.contains(residue)) {
            // The correct index in the theta arrays for tautomer coordinates is after the titration list.
            // So titrationList.size() + tautomerIndex should match with appropriate spot in thetaPosition, etc.
            int index = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
            extendedLambdas[index] = lambda;
            thetaPosition[index] = Math.asin(Math.sqrt(lambda));
            List<Atom> currentAtomList = residue.getSideChainAtoms();
            for (Atom atom : currentAtomList) {
                int atomIndex = atom.getIndex();
                tautomerLambdas[atomIndex] = lambda;
            }
        } else {
            logger.warning(format("This residue %s does not have any titrating tautomers.", residue.getName()));
        }
    }

    /**
     * Update all theta (lambda) postions after each move from the Stochastic integrator
     */
    private void updateLambdas() {
        //This will prevent recalculating multiple sinTheta*sinTheta that are the same number.
        for (int i = 0; i < extendedResidueList.size(); i++) {
            double sinTheta = Math.sin(thetaPosition[i]);
            double oldLambda = extendedLambdas[i];
            extendedLambdas[i] = sinTheta * sinTheta;
        }
        for (int i = 0; i < molecularAssembly.getAtomArray().length; i++) {
            int mappedTitrationIndex = titrationIndexMap[i];
            int mappedTautomerIndex = tautomerIndexMap[i] + titratingResidueList.size();
            if (isTitrating(i) && mappedTitrationIndex != -1) {
                titrationLambdas[i] = extendedLambdas[mappedTitrationIndex];
            }
            if (isTautomerizing(i) && mappedTautomerIndex >= titratingResidueList.size()) {
                tautomerLambdas[i] = extendedLambdas[mappedTautomerIndex];
            }
        }
    }

    public List<Residue> getTitratingResidueList() {
        return titratingResidueList;
    }

    public List<Residue> getTautomerizingResidueList() {
        return tautomerizingResidueList;
    }

    public List<Residue> getExtendedResidueList() {
        return extendedResidueList;
    }

    public double[] getThetaPosition() {
        return thetaPosition;
    }

    public double[] getThetaVelocity() {
        return thetaVelocity;
    }

    public double[] getThetaAccel() {
        return thetaAccel;
    }

    public double[] getThetaMassArray() {
        return thetaMassArray;
    }

    public double getThetaFriction() {
        return thetaFriction;
    }

    public double[] getExtendedLambdas() {
        return extendedLambdas;
    }

    /**
     * getLambdaList.
     *
     * @return a {@link String} object.
     */
    public String getLambdaList() {
        if (extendedResidueList.size() < 1) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < extendedResidueList.size(); i++) {
            if (i == 0) {
                sb.append("Titration Lambdas: ");
            }
            if (i > 0) {
                sb.append(", ");
            }
            if (i == titratingResidueList.size()) {
                sb.append("\n Tautomer Lambdas: ");
            }
            sb.append(format("%6.4f", extendedLambdas[i]));
        }
        return sb.toString();
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
     * get array of dU/dL for each titrating residue
     * @return esvDeriv a double[]
     */
    public double[] getDerivatives() {
        int numESVs = extendedResidueList.size();
        double[] esvDeriv = new double[numESVs];
        double[] biasDerivComponents = new double[9];

        for (Residue residue : titratingResidueList) {
            int resTitrIndex = titratingResidueList.indexOf(residue);
            //Bias Terms
            if (doBias){
                getBiasTerms(residue, biasDerivComponents);
                //Sum up titration bias derivs
                esvDeriv[resTitrIndex] += biasDerivComponents[dDiscr_dTitrIndex] + biasDerivComponents[dPh_dTitrIndex] - biasDerivComponents[dModel_dTitrIndex];
                //Sum up tautomer bias derivs
                if (isTautomer(residue)) {
                    int resTautIndex = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
                    esvDeriv[resTautIndex] += biasDerivComponents[dDiscr_dTautIndex] + biasDerivComponents[dPh_dTautIndex] - biasDerivComponents[dModel_dTautIndex];
                }
            }

            if (doVDW) {
                //TODO: Add vdW
                esvDeriv[resTitrIndex] += getVdwDeriv(resTitrIndex);
                //Sum up tautomer bias derivs
                if (isTautomer(residue)) {
                    int resTautIndex = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
                    esvDeriv[resTautIndex] += getVdwDeriv(resTautIndex);
                }
            }

            if (doElectrostatics) {
                esvDeriv[resTitrIndex] += getPermElecDeriv(resTitrIndex);
                //Sum up tautomer bias derivs
                if (isTautomer(residue)) {
                    int resTautIndex = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
                    esvDeriv[resTautIndex] += getPermElecDeriv(resTautIndex);
                }
                if(doPolarization){
                    esvDeriv[resTitrIndex] += getIndElecDeriv(resTitrIndex);
                    if(isTautomer(residue)){
                        int resTautIndex = tautomerizingResidueList.indexOf(residue) + titratingResidueList.size();
                        esvDeriv[resTautIndex] += getIndElecDeriv(resTautIndex);
                    }
                }
            }

        }


        return esvDeriv;
    }

    /**
     * Sum up total bias (Ubias = UpH + Udiscr - Umod)
     * @return totalBiasEnergy
     */
    public double getBiasEnergy() {
        double totalBiasEnergy = 0.0;
        double[] biasEnergyComponents = new double[9];
        //Do not double count residues in tautomer list.
        for (Residue residue : titratingResidueList) {
            getBiasTerms(residue, biasEnergyComponents);
            double biasEnergy = biasEnergyComponents[discrBiasIndex] + biasEnergyComponents[pHBiasIndex] - biasEnergyComponents[modelBiasIndex];
            totalBiasEnergy += biasEnergy;
        }
        return totalBiasEnergy;
    }

    /**
     * Collect respective pH, model, and discr bias terms and their derivatives for each titrating residue.
     * @param residue
     * @param biasEnergyAndDerivs
     */
    public void getBiasTerms(Residue residue, double[] biasEnergyAndDerivs) {
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        double titrationLambda = getTitrationLambda(residue);
        double discrBias;
        double pHBias;
        double modelBias;
        double dDiscr_dTitr;
        double dDiscr_dTaut;
        double dPh_dTitr;
        double dPh_dTaut;
        double dMod_dTitr;
        double dMod_dTaut;
        //If bias terms shouldn't be computed, set AA3 to UNK so that default case executes and all terms are set to zero.
        if(!doBias){
            AA3 = AminoAcid3.UNK;
        }
        switch (AA3) {
            case ASD:
            case ASH:
            case ASP:
                // Discr Bias & Derivs
                double tautomerLambda = getTautomerLambda(residue);
                discrBias = -4.0 * discrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                discrBias += -4.0 * discrBiasMag * (tautomerLambda - 0.5) * (tautomerLambda - 0.5);
                dDiscr_dTitr = -8.0 * discrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = -8.0 * discrBiasMag * (tautomerLambda - 0.5);

                // pH Bias & Derivs
                double pKa1 = TitrationUtils.Titration.ASHtoASP.pKa;
                double pKa2 = pKa1;
                pHBias = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda)
                        * (tautomerLambda * (pKa1 - constantSystemPh) + (1.0 - tautomerLambda) * (pKa2 - constantSystemPh));
                dPh_dTitr = LOG10 * Constants.R * currentTemperature * -1.0
                        * (tautomerLambda * (pKa1 - constantSystemPh) + (1.0 - tautomerLambda) * (pKa2 - constantSystemPh));
                dPh_dTaut = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda)
                        * ((pKa1 - constantSystemPh) - (pKa2 - constantSystemPh));

                // Model Bias & Derivs
                double refEnergy = TitrationUtils.Titration.ASHtoASP.refEnergy;
                double lambdaIntercept = TitrationUtils.Titration.ASHtoASP.lambdaIntercept;
                modelBias = refEnergy * ((1 - titrationLambda) - lambdaIntercept) * ((1 - titrationLambda) - lambdaIntercept);
                dMod_dTitr = -2.0 * refEnergy * ((1 - titrationLambda) - lambdaIntercept);
                dMod_dTaut = 0.0;
                break;
            case GLD:
            case GLH:
            case GLU:
                // Discr Bias & Derivs
                tautomerLambda = getTautomerLambda(residue);
                discrBias = -4.0 * discrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                discrBias += -4.0 * discrBiasMag * (tautomerLambda - 0.5) * (tautomerLambda - 0.5);
                dDiscr_dTitr = -8.0 * discrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = -8.0 * discrBiasMag * (tautomerLambda - 0.5);

                // pH Bias & Derivs
                pKa1 = TitrationUtils.Titration.GLHtoGLU.pKa;
                pKa2 = pKa1;
                pHBias = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda)
                        * (tautomerLambda * (pKa1 - constantSystemPh) + (1.0 - tautomerLambda) * (pKa2 - constantSystemPh));
                dPh_dTitr = LOG10 * Constants.R * currentTemperature * -1.0
                        * (tautomerLambda * (pKa1 - constantSystemPh) + (1.0 - tautomerLambda) * (pKa2 - constantSystemPh));
                dPh_dTaut = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda)
                        * ((pKa1 - constantSystemPh) - (pKa2 - constantSystemPh));

                // Model Bias & Derivs
                refEnergy = TitrationUtils.Titration.GLHtoGLU.refEnergy;
                lambdaIntercept = TitrationUtils.Titration.GLHtoGLU.lambdaIntercept;
                modelBias = refEnergy * ((1 - titrationLambda) - lambdaIntercept) * ((1 - titrationLambda) - lambdaIntercept);
                dMod_dTitr = -2.0 * refEnergy * ((1 - titrationLambda) - lambdaIntercept);
                dMod_dTaut = 0.0;
                break;
            case HIS:
            case HID:
            case HIE:
                // Discr Bias & Derivs
                tautomerLambda = getTautomerLambda(residue);
                double tautomerLambdaSquared = tautomerLambda * tautomerLambda;
                discrBias = -4.0 * discrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                discrBias += -4.0 * discrBiasMag * (tautomerLambda - 0.5) * (tautomerLambda - 0.5);
                dDiscr_dTitr = -8.0 * discrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = -8.0 * discrBiasMag * (tautomerLambda - 0.5);

                // pH Bias & Derivs
                // At tautomerLambda=1 HIE is fully on.
                pKa1 = TitrationUtils.Titration.HIStoHIE.pKa;
                // At tautomerLambda=0 HID is fully on.
                pKa2 = TitrationUtils.Titration.HIStoHID.pKa;
                pHBias = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda)
                        * (tautomerLambda * (pKa1 - constantSystemPh) + (1.0 - tautomerLambda) * (pKa2 - constantSystemPh));
                dPh_dTitr = LOG10 * Constants.R * currentTemperature * -1.0
                        * (tautomerLambda * (pKa1 - constantSystemPh) + (1.0 - tautomerLambda) * (pKa2 - constantSystemPh));
                dPh_dTaut = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda)
                        * ((pKa1 - constantSystemPh) - (pKa2 - constantSystemPh));

                // Model Bias & Derivs
                double refEnergyHID = TitrationUtils.Titration.HIStoHID.refEnergy;
                double lambdaInterceptHID = TitrationUtils.Titration.HIStoHID.lambdaIntercept;
                double refEnergyHIE = TitrationUtils.Titration.HIStoHIE.refEnergy;
                double lambdaInterceptHIE = TitrationUtils.Titration.HIStoHIE.lambdaIntercept;
                double refEnergyHIDtoHIE = TitrationUtils.Titration.HIDtoHIE.refEnergy;
                double lambdaInterceptHIDtoHIE = TitrationUtils.Titration.HIDtoHIE.lambdaIntercept;

                double coeff4 = -2 * refEnergyHID * lambdaInterceptHID;
                double coeff3 = -2 * refEnergyHIE * lambdaInterceptHIE - coeff4;
                double coeff2 = refEnergyHID;
                double coeff1 = -2 * refEnergyHIDtoHIE * lambdaInterceptHIDtoHIE - coeff3;
                double coeff0 = refEnergyHIDtoHIE;
                double oneMinusTitrationLambda = (1 - titrationLambda);
                double oneMinusTitrationLambdaSquared = (1 - titrationLambda) * (1 - titrationLambda);
                modelBias = oneMinusTitrationLambdaSquared * (coeff0 * tautomerLambdaSquared + coeff1 * tautomerLambda + coeff2)
                        + oneMinusTitrationLambda * (coeff3 * tautomerLambda + coeff4);
                dMod_dTitr = -2 * oneMinusTitrationLambda * (coeff0 * tautomerLambdaSquared + coeff1 * tautomerLambda + coeff2)
                        - (coeff3 * tautomerLambda + coeff4);
                dMod_dTaut = 2 * coeff0 * tautomerLambda * oneMinusTitrationLambdaSquared + coeff1 * oneMinusTitrationLambdaSquared + coeff3 * oneMinusTitrationLambda;
                break;
            case LYS:
            case LYD:
                // Discr Bias & Derivs
                discrBias = -4.0 * discrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                dDiscr_dTitr = -8.0 * discrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = 0.0;

                // pH Bias & Derivs
                pKa1 = TitrationUtils.Titration.LYStoLYD.pKa;
                pHBias = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda) * (pKa1 - constantSystemPh);
                dPh_dTitr = LOG10 * Constants.R * currentTemperature * -1.0 * (pKa1 - constantSystemPh);
                dPh_dTaut = 0.0;

                // Model Bias & Derivs
                refEnergy = TitrationUtils.Titration.LYStoLYD.refEnergy;
                lambdaIntercept = TitrationUtils.Titration.LYStoLYD.lambdaIntercept;
                modelBias = refEnergy * ((1 - titrationLambda) - lambdaIntercept) * ((1 - titrationLambda) - lambdaIntercept);
                dMod_dTitr = -2.0 * refEnergy * ((1 - titrationLambda) - lambdaIntercept);
                dMod_dTaut = 0.0;
                break;
            default:
                discrBias = 0.0;
                pHBias = 0.0;
                modelBias = 0.0;
                dDiscr_dTitr = 0.0;
                dDiscr_dTaut = 0.0;
                dPh_dTitr = 0.0;
                dPh_dTaut = 0.0;
                dMod_dTitr = 0.0;
                dMod_dTaut = 0.0;
                break;
        }
        biasEnergyAndDerivs[0] = discrBias;
        biasEnergyAndDerivs[1] = pHBias;
        biasEnergyAndDerivs[2] = modelBias;
        biasEnergyAndDerivs[3] = dDiscr_dTitr;
        biasEnergyAndDerivs[4] = dPh_dTitr;
        biasEnergyAndDerivs[5] = dMod_dTitr;
        biasEnergyAndDerivs[6] = dDiscr_dTaut;
        biasEnergyAndDerivs[7] = dPh_dTaut;
        biasEnergyAndDerivs[8] = dMod_dTaut;
    }

    /**
     * getBiasDecomposition.
     *
     * @return a {@link String} object.
     */
    public String getBiasDecomposition() {
        double discrBias = 0.0;
        double phBias = 0.0;
        double modelBias = 0.0;
        double[] biasEnergyAndDerivs = new double[9];
        for (Residue residue : extendedResidueList) {
            getBiasTerms(residue, biasEnergyAndDerivs);
            discrBias += biasEnergyAndDerivs[discrBiasIndex];
            phBias += biasEnergyAndDerivs[pHBiasIndex];
            //Reminder: Ubias = UpH + Udiscr - Umod
            modelBias -= biasEnergyAndDerivs[modelBiasIndex];
        }
        return format("    %-16s %16.8f\n", "Discretizer", discrBias)
                + format("    %-16s %16.8f\n", "Acidostat", phBias)
                + format("    %-16s %16.8f\n", "Fmod", modelBias);
    }

    /**
     * getEnergyComponent for use by ForceFieldEnergy
     *
     * @param component a {@link ffx.potential.PotentialComponent} object.
     * @return a double.
     */
    public double getEnergyComponent(PotentialComponent component) {
        double uComp = 0.0;
        double[] biasEnergyAndDerivs = new double[9];
        switch (component) {
            case Bias:
            case pHMD:
                return getBiasEnergy();
            case DiscretizeBias:
                for (Residue residue : extendedResidueList) {
                    getBiasTerms(residue, biasEnergyAndDerivs);
                    uComp += biasEnergyAndDerivs[discrBiasIndex];
                }
                return uComp;
            case ModelBias:
                for (Residue residue : extendedResidueList) {
                    getBiasTerms(residue, biasEnergyAndDerivs);
                    //Reminder: Ubias = UpH + Udiscr - Umod
                    uComp -= biasEnergyAndDerivs[modelBiasIndex];
                }
                return uComp;
            case pHBias:
                for (Residue residue : extendedResidueList) {
                    getBiasTerms(residue, biasEnergyAndDerivs);
                    uComp += biasEnergyAndDerivs[pHBiasIndex];
                }
                return uComp;
            default:
                throw new AssertionError(component.name());
        }
    }

    /**
     * Calculate prefactor for scaling the van der Waals based on titration/tautomer state if titrating proton
     * @param atomIndex
     * @param vdwPrefactorAndDerivs
     */
    public void getVdwPrefactor(int atomIndex, double[] vdwPrefactorAndDerivs) {
        double prefactor = 1.0;
        double titrationDeriv = 0.0;
        double tautomerDeriv = 0.0;
        if (!isTitratingHydrogen(atomIndex) || !doVDW) {
            vdwPrefactorAndDerivs[0] = prefactor;
            vdwPrefactorAndDerivs[1] = titrationDeriv;
            vdwPrefactorAndDerivs[2] = tautomerDeriv;
            return;
        }
        AminoAcid3 AA3 = residueNames[atomIndex];
        switch (AA3) {
            case ASD:
            case GLD:
                if (tautomerDirections[atomIndex] == 1) {
                    prefactor = titrationLambdas[atomIndex] * tautomerLambdas[atomIndex];
                    titrationDeriv = tautomerLambdas[atomIndex];
                    tautomerDeriv = titrationLambdas[atomIndex];
                } else if (tautomerDirections[atomIndex] == -1) {
                    prefactor = titrationLambdas[atomIndex] * (1.0 - tautomerLambdas[atomIndex]);
                    titrationDeriv = (1.0 - tautomerLambdas[atomIndex]);
                    tautomerDeriv = -titrationLambdas[atomIndex];
                }
                break;
            case HIS:
                if (tautomerDirections[atomIndex] == 1) {
                    prefactor = (1.0 - titrationLambdas[atomIndex]) * tautomerLambdas[atomIndex] + titrationLambdas[atomIndex];
                    titrationDeriv = -tautomerLambdas[atomIndex] + 1.0;
                    tautomerDeriv = (1 - titrationLambdas[atomIndex]);
                } else if (tautomerDirections[atomIndex] == -1) {
                    prefactor = (1.0 - titrationLambdas[atomIndex]) * (1.0 - tautomerLambdas[atomIndex]) + titrationLambdas[atomIndex];
                    titrationDeriv = tautomerLambdas[atomIndex];
                    tautomerDeriv = -(1.0 - titrationLambdas[atomIndex]);
                }
                break;
            case LYS:
                prefactor = titrationLambdas[atomIndex];
                titrationDeriv = 1.0;
                tautomerDeriv = 0.0;
                break;
        }
        vdwPrefactorAndDerivs[0] = prefactor;
        vdwPrefactorAndDerivs[1] = titrationDeriv;
        vdwPrefactorAndDerivs[2] = tautomerDeriv;
    }

    /**
     * Add van der Waals deriv to appropriate dU/dL term.
     * @param atomI
     * @param vdwEnergy van der Waals energy calculated with no titration/tautomer scaling
     * @param vdwPrefactorAndDerivI
     * @param vdwPrefactorJ
     */
    public void addVdwDeriv(int atomI, double vdwEnergy, double[] vdwPrefactorAndDerivI, double vdwPrefactorJ) {
        if (!isTitratingHydrogen(atomI)) {
            return;
        }

        //Sum up dU/dL for titration ESV if atom i is titrating hydrogen
        //Sum up dU/dL for tautomer ESV if atom i is titrating hydrogen
        int titrationEsvIndex = titrationIndexMap[atomI];
        int tautomerEsvIndex = tautomerIndexMap[atomI] + titratingResidueList.size();
        double dTitr_dLambda;
        double dTaut_dLambda;

        dTitr_dLambda = vdwPrefactorAndDerivI[1] * vdwPrefactorJ * vdwEnergy;
        dTaut_dLambda = vdwPrefactorAndDerivI[2] * vdwPrefactorJ * vdwEnergy;

        esvVdwDerivs[titrationEsvIndex].addAndGet(dTitr_dLambda);
        if(tautomerEsvIndex >= titratingResidueList.size()){
            esvVdwDerivs[tautomerEsvIndex].addAndGet(dTaut_dLambda);
        }
    }

    public void addPermElecDeriv(int atomI, double titrationEnergy, double tautomerEnergy){
        int titrationEsvIndex = titrationIndexMap[atomI];
        int tautomerEsvIndex = tautomerIndexMap[atomI] + titratingResidueList.size();
        esvPermElecDerivs[titrationEsvIndex].addAndGet(titrationEnergy);
        if(tautomerEsvIndex >= titratingResidueList.size()){
            esvPermElecDerivs[tautomerEsvIndex].addAndGet(tautomerEnergy);
        }
    }

    public void addIndElecDeriv(int atomI, double titrationEnergy, double tautomerEnergy){
        int titrationEsvIndex = titrationIndexMap[atomI];
        int tautomerEsvIndex = tautomerIndexMap[atomI] + titratingResidueList.size();
        esvIndElecDerivs[titrationEsvIndex].addAndGet(titrationEnergy);
        if(tautomerEsvIndex >= titratingResidueList.size()){
            esvIndElecDerivs[tautomerEsvIndex].addAndGet(tautomerEnergy);
        }
    }

    /**
     * get total dUvdw/dL for the selected extended system variable
     * @param esvID
     * @return
     */
    public double getVdwDeriv(int esvID) {
        return esvVdwDerivs[esvID].get();
    }

    public double getPermElecDeriv(int esvID){return esvPermElecDerivs[esvID].get();}

    public double getIndElecDeriv(int esvID){return esvIndElecDerivs[esvID].get();}

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
     *
     * @return titrationUtils master copy from ExtendedSystem.
     */
    public TitrationUtils getTitrationUtils(){
        return titrationUtils;
    }

    /**
     * setTemperature.
     *
     * @param set a double.
     */
    public void setTemperature(double set) {
        currentTemperature = set;
    }

    /**
     * Processes lambda values based on propagation of theta value from Stochastic integrator in Molecular dynamics
     */
    public void preForce() {
        updateLambdas();
    }

    /**
     * Applies a chain rule term to the derivative to account for taking a derivative of lambda = sin(theta)^2
     *
     * @return dE/dL a double[]
     */
    public double[] postForce() {
        double[] dEdL = ExtendedSystem.this.getDerivatives();
        double[] dEdTheta = new double[dEdL.length];
        for (int i = 0; i < extendedResidueList.size(); i++) {
            dEdTheta[i] = dEdL[i] * sin(2 * thetaPosition[i]);
        }
        return dEdTheta;
    }

    /*private void updateMultipoleTypes() {
        for (Residue residue : extendedResidueList){
            AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
            List<Atom> atomList = residue.get;
            for(Atom atom : atomList){
                int atomIndex = atom.getArrayIndex();
                MultipoleType atomMultipoleType = atom.getMultipoleType();
                double[] esvMultipole = TitrationUtils.getMultipole(AA3, atomIndex, titrationLambdas[atomIndex],
                        tautomerLambdas[atomIndex], atomMultipoleType.getMultipole());
                double[] esvTitrationDotMultipole = TitrationUtils.getMultipoleTitrationDeriv(AA3, atomIndex,
                        titrationLambdas[atomIndex], tautomerLambdas[atomIndex], atomMultipoleType.getMultipole());
                double[] esvTautomerDotMultipole = TitrationUtils.getMultipoleTautomerDeriv(AA3, atomIndex,
                        titrationLambdas[atomIndex], tautomerLambdas[atomIndex], atomMultipoleType.getMultipole());
                MultipoleType esvType = new MultipoleType(esvMultipole,atomMultipoleType.frameAtomTypes,
                        atomMultipoleType.frameDefinition, false);
                MultipoleType esvTitrationDotType = new MultipoleType(esvTitrationDotMultipole, atomMultipoleType.frameAtomTypes,
                        atomMultipoleType.frameDefinition, false);
                MultipoleType esvTautomerDotType =  new MultipoleType(esvTautomerDotMultipole, atomMultipoleType.frameAtomTypes,
                        atomMultipoleType.frameDefinition, false);
                //TODO: Detect hydrogen and scale alpha.
                atom.setEsv( esvType, esvTitrationDotType);
            }
        }
    }*/
}
