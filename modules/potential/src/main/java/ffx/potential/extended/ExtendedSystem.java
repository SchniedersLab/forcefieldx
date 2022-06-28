// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2022.
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

import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.reduction.SharedDouble;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.PotentialComponent;
import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TitrationUtils;
import ffx.potential.parsers.ESVFilter;
import ffx.utilities.Constants;
import ffx.utilities.FileUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * ExtendedSystem class.
 *
 * @author Andrew Thiel
 * @since 1.0
 */
public class ExtendedSystem implements Potential {

    private static final Logger logger = Logger.getLogger(ExtendedSystem.class.getName());

    private static final double THETA_MASS = 5.0; //Atomic Mass Units

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
     * Array of AminoAcid3 initialized  to match the number of atoms in the system.
     * Used to know how to apply vdW or electrostatic ESV terms for the atom.
     */
    public final AminoAcid3[] residueNames;
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
     * Extended Atoms
     */
    private final Atom[] extendedAtoms;
    /**
     * Extended Molecules
     */
    private final int[] extendedMolecules;
    /**
     * Titration Utils instance. This instance is the master copy that will be distributed to electrostatics classes
     * when an Extended System is attached.
     */
    private final TitrationUtils titrationUtils;

    /**
     * Number of atoms in the molecular assembly. Since all protons are instantiated at start, this int will not change.
     */
    private final int nAtoms;
    /**
     * Number of ESVs attached to the molecular assembly. This number is the sum of tautomerizing + titrating ESVs.
     */
    private final int nESVs;
    /**
     * Number of titrating ESVs attached to the molecular assembly.
     */
    private final int nTitr;
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
    private final boolean[] isTitratingSulfur;
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

    // Controls for turning of certain terms for testing.
    private final boolean doVDW;
    private final boolean doElectrostatics;
    private final boolean doBias;
    private final boolean doPolarization;
    private final boolean fixTitrationState;
    private final boolean fixTautomerState;
    private final boolean lockStates;
    /**
     * Current value of theta for each ESV.
     */
    final private double[] thetaPosition;
    /**
     * Current theta velocity for each ESV.
     */
    final private double[] thetaVelocity;
    /**
     * Current theta acceleration for each ESV.
     */
    final private double[] thetaAccel;
    /**
     * Mass of each theta particle. Different theta mass for each particle are not supported.
     */
    final private double[] thetaMassArray;
    /**
     * 3D array to store the titration and tautomer population states for each ESV
     */
    final private int[][][] esvHistogram;
    /**
     * Descritizer Bias Magnitude. Default is 1 kcal/mol.
     */
    private final double titrBiasMag;
    private final double tautBiasMag;
    private final double HIStitrBiasMag;
    private final double HIStautBiasMag;
    /**
     * Coefficients that define the per residue type Fmod polynomial
     * refEnergy * titrationLambda^2 + lambdaIntercept * titrationLambda
     */
    final private double ASHrefEnergy;
    final private double ASHlambdaIntercept;
    final private double GLHrefEnergy;
    final private double GLHlambdaIntercept;
    final private double LYSrefEnergy;
    final private double LYSlambdaIntercept;
    final private double CYSrefEnergy;
    final private double CYSlambdaIntercept;
    final private double HIDrefEnergy;
    final private double HIDlambdaIntercept;
    final private double HIErefEnergy;
    final private double HIElambdaIntercept;
    /**
     * Coefficients that define the tautomer component of the bivariate Histidine Fmod
     * refEnergy * tautomerLambda^2 + lambdaIntercept * tautomerLambda
     */
    final private double HIDtoHIErefEnergy;
    final private double HIDtoHIElambdaIntercept;
    /**
     * Dynamics restart file.
     */
    File restartFile = null;
    /**
     * Filter to parse the dynamics restart file.
     */
    ESVFilter esvFilter = null;
    /**
     * System PH.
     */
    private double constantSystemPh = 7.4;
    /**
     * Target system temperature.
     */
    private double currentTemperature = Constants.ROOM_TEMPERATURE;

    /**
     * Construct extended system with the provided configuration.
     *
     * @param mola a {@link MolecularAssembly} object.
     */
    public ExtendedSystem(MolecularAssembly mola, final File esvFile) {
        extendedAtoms = mola.getActiveAtomArray();
        extendedMolecules = mola.getMoleculeNumbers();

        ForceField forceField = mola.getForceField();
        ForceFieldEnergy forceFieldEnergy = mola.getPotentialEnergy();
        if (forceFieldEnergy == null) {
            logger.severe("No potential energy found?");
        }
        CompositeConfiguration properties = mola.getProperties();
        titrationUtils = new TitrationUtils(forceField);

        doVDW = properties.getBoolean("esv.vdW", true);
        doElectrostatics = properties.getBoolean("esv.elec", true);
        doBias = properties.getBoolean("esv.bias", true);
        doPolarization = properties.getBoolean("esv.polarization", true);
//        boolean bonded = properties.getBoolean("esv.bonded", false);
//        boolean verbose = properties.getBoolean("esv.verbose", false);
//        boolean decomposeBonded = properties.getBoolean("esv.decomposeBonded", false);
//        boolean decomposeDeriv = properties.getBoolean("esv.decomposeDeriv", false);
//        boolean nonlinearMultipoles = properties.getBoolean("esv.nonlinearMultipoles", false); // sigmoid lambda Mpole switch
//        boolean forceRoomTemp = properties.getBoolean("esv.forceRoomTemp", false);
//        boolean propagation = properties.getBoolean("esv.propagation", true);
        thetaFriction = properties.getDouble("esv.friction", ExtendedSystem.THETA_FRICTION);
        thetaMass = properties.getDouble("esv.mass", ExtendedSystem.THETA_MASS);
        titrBiasMag = properties.getDouble("titration.bias.magnitude", DISCR_BIAS);
        tautBiasMag = properties.getDouble("tautomer.bias.magnitude", DISCR_BIAS);
        HIStitrBiasMag = properties.getDouble("HIS.titration.bias.magnitude", DISCR_BIAS);
        HIStautBiasMag = properties.getDouble("HIS.tautomer.bias.magnitude", DISCR_BIAS);
        lockStates = properties.getBoolean("lock.esv.states", false); // Prevents setTitrationLambda/setTautomerLambda
        double initialTitrationLambda = properties.getDouble("lambda.titration.initial", 0.5);
        double initialTautomerLambda = properties.getDouble("lambda.tautomer.initial", 0.5);
        boolean guessTitrState = properties.getBoolean("guess.titration.state", false);
        fixTitrationState = properties.getBoolean("fix.titration.lambda", false);
        fixTautomerState = properties.getBoolean("fix.tautomer.lambda", false);
        ASHrefEnergy = properties.getDouble("ASH.ref.energy", TitrationUtils.Titration.ASHtoASP.refEnergy);
        ASHlambdaIntercept = properties.getDouble("ASH.lambda.intercept", TitrationUtils.Titration.ASHtoASP.lambdaIntercept);
        GLHrefEnergy = properties.getDouble("GLH.ref.energy", TitrationUtils.Titration.GLHtoGLU.refEnergy);
        GLHlambdaIntercept = properties.getDouble("GLH.lambda.intercept", TitrationUtils.Titration.GLHtoGLU.lambdaIntercept);
        LYSrefEnergy = properties.getDouble("LYS.ref.energy", TitrationUtils.Titration.LYStoLYD.refEnergy);
        LYSlambdaIntercept = properties.getDouble("LYS.lambda.intercept", TitrationUtils.Titration.LYStoLYD.lambdaIntercept);
        CYSrefEnergy = properties.getDouble("CYS.ref.energy", TitrationUtils.Titration.CYStoCYD.refEnergy);
        CYSlambdaIntercept = properties.getDouble("CYS.lambda.intercept", TitrationUtils.Titration.CYStoCYD.lambdaIntercept);
        HIDrefEnergy = properties.getDouble("HID.ref.energy", TitrationUtils.Titration.HIStoHID.refEnergy);
        HIDlambdaIntercept = properties.getDouble("HID.lambda.intercept", TitrationUtils.Titration.HIStoHID.lambdaIntercept);
        HIErefEnergy = properties.getDouble("HIE.ref.energy", TitrationUtils.Titration.HIStoHIE.refEnergy);
        HIElambdaIntercept = properties.getDouble("HIE.lambda.intercept", TitrationUtils.Titration.HIStoHIE.lambdaIntercept);
        HIDtoHIErefEnergy = properties.getDouble("HIDtoHIE.ref.energy", TitrationUtils.Titration.HIDtoHIE.refEnergy);
        HIDtoHIElambdaIntercept = properties.getDouble("HIDtoHIE.lambda.intercept", TitrationUtils.Titration.HIDtoHIE.lambdaIntercept);

        titratingResidueList = new ArrayList<>();
        tautomerizingResidueList = new ArrayList<>();
        extendedResidueList = new ArrayList<>();
        // Initialize atom arrays with the existing assembly.
        Atom[] atoms = mola.getAtomArray();
        nAtoms = atoms.length;
        isTitrating = new boolean[nAtoms];
        isTitratingHydrogen = new boolean[nAtoms];
        isTitratingSulfur = new boolean[nAtoms];
        isTautomerizing = new boolean[nAtoms];
        titrationLambdas = new double[nAtoms];
        tautomerLambdas = new double[nAtoms];
        titrationIndexMap = new int[nAtoms];
        tautomerIndexMap = new int[nAtoms];
        tautomerDirections = new int[nAtoms];
        residueNames = new AminoAcid3[nAtoms];

        Arrays.fill(isTitrating, false);
        Arrays.fill(isTitratingHydrogen, false);
        Arrays.fill(isTitratingSulfur, false);
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
        List<Residue> residueList = mola.getResidueList();
        // Use only a list that contains AminoAcid residues so remove Nucleic Acid residues
        residueList.removeIf(residue -> (residue.getResidueType() == Residue.ResidueType.NA));
        for (Residue residue : residueList) {
            if (isTitratable(residue)) {
                titratingResidueList.add(residue);
                List<Atom> atomList = residue.getSideChainAtoms();
                for (Atom atom : atomList) {
                    int atomIndex = atom.getArrayIndex();
                    residueNames[atomIndex] = residue.getAminoAcid3();
                    isTitrating[atomIndex] = true;
                    if(guessTitrState){
                        double guessTitrLambda = initialTitrationState(residue, initialTitrationLambda);
                        titrationLambdas[atomIndex] = guessTitrLambda;
                    } else{
                        titrationLambdas[atomIndex] = initialTitrationLambda;
                    }
                    int titrationIndex = titratingResidueList.indexOf(residue);
                    titrationIndexMap[atomIndex] = titrationIndex;
                    isTitratingHydrogen[atomIndex] = TitrationUtils.isTitratingHydrogen(residue.getAminoAcid3(), atom);
                    isTitratingSulfur[atomIndex] = TitrationUtils.isTitratingSulfur(residue.getAminoAcid3(), atom);
                }
                // If is a tautomer, it must also be titrating.
                if (isTautomer(residue)) {
                    tautomerizingResidueList.add(residue);
                    for (Atom atom : atomList) {
                        int atomIndex = atom.getArrayIndex();
                        isTautomerizing[atomIndex] = true;
                        tautomerLambdas[atomIndex] = initialTautomerLambda;
                        int tautomerIndex = tautomerizingResidueList.indexOf(residue);
                        tautomerIndexMap[atomIndex] = tautomerIndex;
                        tautomerDirections[atomIndex] = TitrationUtils.getTitratingHydrogenDirection(residue.getAminoAcid3(), atom);
                    }
                }
                // Test the multipole frames during unit testing.
                assert (titrationUtils.testResidueTypes(residue));
            }
        }

        //Concatenate titratingResidueList and tautomerizingResidueList
        extendedResidueList.addAll(titratingResidueList);
        extendedResidueList.addAll(tautomerizingResidueList);
        //Arrays that are sent to integrator are based on extendedResidueList size
        nESVs = extendedResidueList.size();
        nTitr = titratingResidueList.size();
        extendedLambdas = new double[nESVs];
        thetaPosition = new double[nESVs];
        thetaVelocity = new double[nESVs];
        thetaAccel = new double[nESVs];
        thetaMassArray = new double[nESVs];
        esvHistogram = new int[nTitr][10][10];
        esvVdwDerivs = new SharedDouble[nESVs];
        esvPermElecDerivs = new SharedDouble[nESVs];
        esvIndElecDerivs = new SharedDouble[nESVs];
        for (int i = 0; i < nESVs; i++) {
            esvVdwDerivs[i] = new SharedDouble(0.0);
            esvPermElecDerivs[i] = new SharedDouble(0.0);
            esvIndElecDerivs[i] = new SharedDouble(0.0);
        }

        //Theta masses should always be the same for each ESV
        Arrays.fill(thetaMassArray, thetaMass);

        for (int i = 0; i < nESVs; i++) {
            if (i < nTitr) {
                if(guessTitrState){
                    Residue residue = extendedResidueList.get(i);
                    double initialTitrLambda = initialTitrationState(residue, initialTitrationLambda);
                    initializeThetaArrays(i, initialTitrLambda);
                } else{
                    initializeThetaArrays(i, initialTitrationLambda);
                }
            } else {
                initializeThetaArrays(i, initialTautomerLambda);
            }
        }
        if (esvFilter == null) {
            esvFilter = new ESVFilter(mola.getName());
        }
        if (esvFile == null) {
            String firstFileName = FilenameUtils.removeExtension(mola.getFile().getAbsolutePath());
            restartFile = new File(firstFileName + ".esv");
        } else {
            if (!esvFilter.readESV(esvFile, thetaPosition, thetaVelocity, thetaAccel)) {
                String message = " Could not load the restart file - dynamics terminated.";
                logger.log(Level.WARNING, message);
                throw new IllegalStateException(message);
            } else {
                restartFile = esvFile;
                updateLambdas();
            }
        }
    }

    /**
     * During constructor, initialize arrays that will hold theta positions, velocities, and accelerations.
     * Positions determined from starting lambda.
     * Velocities randomly set according to Maxwell Boltzmann Distribution based on temperature.
     * Accelerations determined from initial forces.
     *
     * @param index  index of ExtendedResidueList for which to set these values.
     * @param lambda starting lambda value for each ESV.
     */
    private void initializeThetaArrays(int index, double lambda) {
        extendedLambdas[index] = lambda;
        thetaPosition[index] = Math.asin(Math.sqrt(lambda));
        Random random = new Random();
        thetaVelocity[index] = random.nextGaussian() * sqrt(kB * 298.15 / thetaMass);
        double dUdL = getDerivatives()[index];
        double dUdTheta = dUdL * sin(2 * thetaPosition[index]);
        thetaAccel[index] = -Constants.KCAL_TO_GRAM_ANG2_PER_PS2 * dUdTheta / thetaMass;
    }

    /**
     * get array of dU/dL for each titrating residue
     *
     * @return esvDeriv a double[]
     */
    public double[] getDerivatives() {
        double[] esvDeriv = new double[nESVs];
        double[] biasDerivComponents = new double[9];

        for (Residue residue : titratingResidueList) {
            int resTitrIndex = titratingResidueList.indexOf(residue);
            //Bias Terms
            if (doBias) {
                getBiasTerms(residue, biasDerivComponents);
                //Sum up titration bias derivs
                esvDeriv[resTitrIndex] += biasDerivComponents[dDiscr_dTitrIndex] + biasDerivComponents[dPh_dTitrIndex] + biasDerivComponents[dModel_dTitrIndex];
                //Sum up tautomer bias derivs
                if (isTautomer(residue)) {
                    int resTautIndex = tautomerizingResidueList.indexOf(residue) + nTitr;
                    esvDeriv[resTautIndex] += biasDerivComponents[dDiscr_dTautIndex] + biasDerivComponents[dPh_dTautIndex] + biasDerivComponents[dModel_dTautIndex];
                }
            }

            if (doVDW) {
                esvDeriv[resTitrIndex] += getVdwDeriv(resTitrIndex);
                //Sum up tautomer bias derivs
                if (isTautomer(residue)) {
                    int resTautIndex = tautomerizingResidueList.indexOf(residue) + nTitr;
                    esvDeriv[resTautIndex] += getVdwDeriv(resTautIndex);
                }
            }

            if (doElectrostatics) {
                esvDeriv[resTitrIndex] += getPermElecDeriv(resTitrIndex);
                //Sum up tautomer bias derivs
                if (isTautomer(residue)) {
                    int resTautIndex = tautomerizingResidueList.indexOf(residue) + nTitr;
                    esvDeriv[resTautIndex] += getPermElecDeriv(resTautIndex);
                }
                if (doPolarization) {
                    esvDeriv[resTitrIndex] += getIndElecDeriv(resTitrIndex);
                    if (isTautomer(residue)) {
                        int resTautIndex = tautomerizingResidueList.indexOf(residue) + nTitr;
                        esvDeriv[resTautIndex] += getIndElecDeriv(resTautIndex);
                    }
                }
            }
        }
        return esvDeriv;
    }

    /**
     * Collect respective pH, model, and discr bias terms and their derivatives for each titrating residue.
     *
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
        if (!doBias) {
            AA3 = AminoAcid3.UNK;
        }
        switch (AA3) {
            case ASD:
            case ASH:
            case ASP:
                // Discr Bias & Derivs
                double tautomerLambda = getTautomerLambda(residue);
                discrBias = -4.0 * titrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                discrBias += -4.0 * tautBiasMag * (tautomerLambda - 0.5) * (tautomerLambda - 0.5);
                dDiscr_dTitr = -8.0 * titrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = -8.0 * tautBiasMag * (tautomerLambda - 0.5);

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
                double refEnergy = ASHrefEnergy;
                double lambdaIntercept = ASHlambdaIntercept;
                modelBias = refEnergy * titrationLambda*titrationLambda + lambdaIntercept * titrationLambda;
                dMod_dTitr = 2 * refEnergy * titrationLambda + lambdaIntercept;
                dMod_dTaut = 0.0;
                break;
            case GLD:
            case GLH:
            case GLU:
                // Discr Bias & Derivs
                tautomerLambda = getTautomerLambda(residue);
                discrBias = -4.0 * titrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                discrBias += -4.0 * tautBiasMag * (tautomerLambda - 0.5) * (tautomerLambda - 0.5);
                dDiscr_dTitr = -8.0 * titrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = -8.0 * tautBiasMag * (tautomerLambda - 0.5);

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
                refEnergy = GLHrefEnergy;
                lambdaIntercept = GLHlambdaIntercept;
                modelBias = refEnergy * titrationLambda*titrationLambda + lambdaIntercept * titrationLambda;
                dMod_dTitr = 2 * refEnergy * titrationLambda + lambdaIntercept;
                dMod_dTaut = 0.0;
                break;
            case HIS:
            case HID:
            case HIE:
                // Discr Bias & Derivs
                tautomerLambda = getTautomerLambda(residue);

                discrBias = -4.0 * HIStitrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                discrBias += -4.0 * HIStautBiasMag * (tautomerLambda - 0.5) * (tautomerLambda - 0.5);
                dDiscr_dTitr = -8.0 * HIStitrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = -8.0 * HIStautBiasMag * (tautomerLambda - 0.5);

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

                double coeffA0 = HIDtoHIErefEnergy;
                double coeffA1 = HIDtoHIElambdaIntercept;
                double coeffB0 = HIErefEnergy;
                double coeffB1 = HIElambdaIntercept;
                double coeffC0 = HIDrefEnergy;
                double coeffC1 = HIDlambdaIntercept;
                double titrationLambdaSquared = titrationLambda * titrationLambda;
                double tautomerLambdaSquared = tautomerLambda * tautomerLambda;
                double oneMinusTitrationLambda = (1.0 - titrationLambda);
                double oneMinusTautomerLambda = (1.0 - tautomerLambda);
                double coeffBSum = coeffB0 + coeffB1;
                double coeffCSum = coeffC0 + coeffC1;
                modelBias = oneMinusTitrationLambda * (coeffA0 * tautomerLambdaSquared + coeffA1 * tautomerLambda) +
                            tautomerLambda * (coeffB0 * titrationLambdaSquared + coeffB1 * titrationLambda) +
                            oneMinusTautomerLambda * (coeffC0 * titrationLambdaSquared + coeffC1 * titrationLambda) +
                            //Enforce that HIS(titration==1) state is equal energy no matter tautomer value
                            titrationLambda * (coeffCSum - coeffBSum) * tautomerLambda;
                dMod_dTitr = -(coeffA0 * tautomerLambdaSquared + coeffA1 * tautomerLambda) +
                             tautomerLambda * (2.0 * coeffB0 * titrationLambda + coeffB1) +
                             oneMinusTautomerLambda * (2.0 * coeffC0 * titrationLambda + coeffC1) +
                             tautomerLambda * (coeffCSum - coeffBSum);
                dMod_dTaut = oneMinusTitrationLambda * (2.0 * coeffA0 * tautomerLambda + coeffA1) +
                             (coeffB0 * titrationLambdaSquared + coeffB1 * titrationLambda) +
                             -(coeffC0 * titrationLambdaSquared + coeffC1 * titrationLambda) +
                             titrationLambda * (coeffCSum - coeffBSum);

                break;
            case LYS:
            case LYD:
                // Discr Bias & Derivs
                discrBias = -4.0 * titrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                dDiscr_dTitr = -8.0 * titrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = 0.0;

                // pH Bias & Derivs
                pKa1 = TitrationUtils.Titration.LYStoLYD.pKa;
                pHBias = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda) * (pKa1 - constantSystemPh);
                dPh_dTitr = LOG10 * Constants.R * currentTemperature * -1.0 * (pKa1 - constantSystemPh);
                dPh_dTaut = 0.0;

                // Model Bias & Derivs
                refEnergy = LYSrefEnergy;
                lambdaIntercept = LYSlambdaIntercept;
                modelBias = refEnergy * titrationLambda * titrationLambda + lambdaIntercept * titrationLambda;
                dMod_dTitr = 2 * refEnergy * titrationLambda + lambdaIntercept;
                dMod_dTaut = 0.0;
                break;
            case CYS:
            case CYD:
                // Discr Bias & Derivs
                discrBias = -4.0 * titrBiasMag * (titrationLambda - 0.5) * (titrationLambda - 0.5);
                dDiscr_dTitr = -8.0 * titrBiasMag * (titrationLambda - 0.5);
                dDiscr_dTaut = 0.0;

                // pH Bias & Derivs
                pKa1 = TitrationUtils.Titration.CYStoCYD.pKa;
                pHBias = LOG10 * Constants.R * currentTemperature * (1.0 - titrationLambda) * (pKa1 - constantSystemPh);
                dPh_dTitr = LOG10 * Constants.R * currentTemperature * -1.0 * (pKa1 - constantSystemPh);
                dPh_dTaut = 0.0;

                // Model Bias & Derivs
                refEnergy = CYSrefEnergy;
                lambdaIntercept = CYSlambdaIntercept;
                modelBias = refEnergy * titrationLambda * titrationLambda + lambdaIntercept * titrationLambda;
                dMod_dTitr = 2 * refEnergy * titrationLambda + lambdaIntercept;
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
        biasEnergyAndDerivs[2] = -modelBias;
        biasEnergyAndDerivs[3] = dDiscr_dTitr;
        biasEnergyAndDerivs[4] = dPh_dTitr;
        biasEnergyAndDerivs[5] = -dMod_dTitr;
        biasEnergyAndDerivs[6] = dDiscr_dTaut;
        biasEnergyAndDerivs[7] = dPh_dTaut;
        biasEnergyAndDerivs[8] = -dMod_dTaut;
    }

    public double getTitrationLambda(Residue residue) {
        if (titratingResidueList.contains(residue)) {
            int resIndex = titratingResidueList.indexOf(residue);
            return extendedLambdas[resIndex];
        } else {
            return 1.0;
        }
    }

    public double getTautomerLambda(Residue residue) {
        if (tautomerizingResidueList.contains(residue)) {
            int resIndex = tautomerizingResidueList.indexOf(residue);
            return extendedLambdas[nTitr + resIndex];
        } else {
            return 1.0;
        }
    }

    /**
     * get total dUvdw/dL for the selected extended system variable
     *
     * @param esvID
     * @return
     */
    public double getVdwDeriv(int esvID) {
        return esvVdwDerivs[esvID].get();
    }

    public double getPermElecDeriv(int esvID) {
        return esvPermElecDerivs[esvID].get();
    }

    public double getIndElecDeriv(int esvID) {
        return esvIndElecDerivs[esvID].get();
    }

    public boolean isTitratable(Residue residue) {
        if(residue.getResidueType() == Residue.ResidueType.NA){
            return false;
        }
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        return AA3.isConstantPhTitratable;
    }

    public boolean isTautomer(Residue residue) {
        if(residue.getResidueType() == Residue.ResidueType.NA){
            return false;
        }
        AminoAcidUtils.AminoAcid3 AA3 = residue.getAminoAcid3();
        return AA3.isConstantPhTautomer;
    }

    /**
     * Update all theta (lambda) postions after each move from the Stochastic integrator
     */
    private void updateLambdas() {
        //This will prevent recalculating multiple sinTheta*sinTheta that are the same number.
        for (int i = 0; i < nESVs; i++) {
            //Check to see if titration/tautomer lambdas are to be fixed
            if ((!fixTitrationState && i < nTitr) || (!fixTautomerState && i >= nTitr)) {
                double sinTheta = Math.sin(thetaPosition[i]);
                extendedLambdas[i] = sinTheta * sinTheta;
            }
        }
        for (int i = 0; i < nAtoms; i++) {
            int mappedTitrationIndex = titrationIndexMap[i];
            int mappedTautomerIndex = tautomerIndexMap[i] + nTitr;
            if (isTitrating(i) && mappedTitrationIndex != -1) {
                titrationLambdas[i] = extendedLambdas[mappedTitrationIndex];
            }
            if (isTautomerizing(i) && mappedTautomerIndex >= nTitr) {
                tautomerLambdas[i] = extendedLambdas[mappedTautomerIndex];
            }
        }
        setESVHistogram();
    }

    public boolean isTitrating(int atomIndex) {
        return isTitrating[atomIndex];
    }

    public boolean isTautomerizing(int atomIndex) {
        return isTautomerizing[atomIndex];
    }

    private void setESVHistogram() {
        for (Residue residue : titratingResidueList) {
            int index = titratingResidueList.indexOf(residue);
            if (residue.getAminoAcid3().equals(AminoAcid3.LYS)) {
                double titrLambda = getTitrationLambda(residue);
                esvHistogram(index, titrLambda);
            } else {
                double titrLambda = getTitrationLambda(residue);
                double tautLambda = getTautomerLambda(residue);
                esvHistogram(index, titrLambda, tautLambda);
            }
        }
    }

    private void esvHistogram(int esv, double lambda) {
        int value = (int) (lambda * 10.0);
        //Cover the case where lambda could be exactly 1.0
        if (value == 10) {
            value = 9;
        }
        esvHistogram[esv][value][0]++;
    }

    private void esvHistogram(int esv, double titrLambda, double tautLambda) {
        int titrValue = (int) (titrLambda * 10.0);
        //Cover the case where lambda could be exactly 1.0
        if (titrValue == 10) {
            titrValue = 9;
        }
        int tautValue = (int) (tautLambda * 10.0);
        if (tautValue == 10) {
            tautValue = 9;
        }
        esvHistogram[esv][titrValue][tautValue]++;
    }

    public void initEsvVdw() {
        for (int i = 0; i < nESVs; i++) {
            esvVdwDerivs[i].set(0.0);
        }
    }

    public void initEsvPermElec() {
        for (int i = 0; i < nESVs; i++) {
            esvPermElecDerivs[i].set(0.0);
        }
    }

    public void initEsvIndElec() {
        for (int i = 0; i < nESVs; i++) {
            esvIndElecDerivs[i].set(0.0);
        }
    }

    public boolean isExtended(Residue residue) {
        return extendedResidueList.contains(residue);
    }

    public double getTitrationLambda(int atomIndex) {
        return titrationLambdas[atomIndex];
    }

    public int getTitrationESVIndex(int i) {
        return titrationIndexMap[i];
    }

    public double getTautomerLambda(int atomIndex) {
        return tautomerLambdas[atomIndex];
    }

    public int getTautomerESVIndex(int i) {
        return tautomerIndexMap[i];
    }

    public void setTitrationLambda(Residue residue, double lambda) {
        if (titratingResidueList.contains(residue) && !lockStates) {
            int index = titratingResidueList.indexOf(residue);
            extendedLambdas[index] = lambda;
            thetaPosition[index] = Math.asin(Math.sqrt(lambda));
            List<Atom> currentAtomList = residue.getSideChainAtoms();
            for (Atom atom : currentAtomList) {
                int atomIndex = atom.getArrayIndex();
                titrationLambdas[atomIndex] = lambda;
            }
        } /*else {
            logger.warning(format("This residue %s is not titrating or locked by user property.", residue.getName()));
        }*/
    }

    public void setTautomerLambda(Residue residue, double lambda) {
        if (tautomerizingResidueList.contains(residue) && !lockStates) {
            // The correct index in the theta arrays for tautomer coordinates is after the titration list.
            // So titrationList.size() + tautomerIndex should match with appropriate spot in thetaPosition, etc.
            int index = tautomerizingResidueList.indexOf(residue) + nTitr;
            extendedLambdas[index] = lambda;
            thetaPosition[index] = Math.asin(Math.sqrt(lambda));
            List<Atom> currentAtomList = residue.getSideChainAtoms();
            for (Atom atom : currentAtomList) {
                int atomIndex = atom.getArrayIndex();
                tautomerLambdas[atomIndex] = lambda;
            }
        } /*else {
            logger.warning(format("This residue %s does not have any titrating tautomers.", residue.getName()));
        }*/
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

    public double getThetaMass() {
        return thetaMass;
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
        if (nESVs < 1) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < nESVs; i++) {
            if (i == 0) {
                sb.append("\n  Titration Lambdas: ");
            }
            if (i > 0) {
                sb.append(", ");
            }
            if (i == nTitr) {
                sb.append("\n  Tautomer Lambdas: ");
            }
            sb.append(format("%6.4f", extendedLambdas[i]));
        }
        return sb.toString();
    }

    public double[] getLambdaArray(){
        double[] lambdaArray = new double[extendedResidueList.size()];
        System.arraycopy(extendedLambdas,0, lambdaArray, 0, lambdaArray.length);
        return lambdaArray;
    }
    public void setLambdaArray(double[] lambdaArray){
        System.arraycopy(lambdaArray, 0, extendedLambdas, 0, extendedLambdas.length);
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
        return extendedMolecules;
    }

    /**
     * Sum up total bias (Ubias = UpH + Udiscr - Umod)
     *
     * @return totalBiasEnergy
     */
    public double getBiasEnergy() {
        double totalBiasEnergy = 0.0;
        double[] biasEnergyComponents = new double[9];
        //Do not double count residues in tautomer list.
        for (Residue residue : titratingResidueList) {
            getBiasTerms(residue, biasEnergyComponents);
            double biasEnergy = biasEnergyComponents[discrBiasIndex] + biasEnergyComponents[pHBiasIndex] + biasEnergyComponents[modelBiasIndex];
            totalBiasEnergy += biasEnergy;
        }
        return totalBiasEnergy;
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
        for (Residue residue : titratingResidueList) {
            getBiasTerms(residue, biasEnergyAndDerivs);
            discrBias += biasEnergyAndDerivs[discrBiasIndex];
            phBias += biasEnergyAndDerivs[pHBiasIndex];
            //Reminder: Ubias = UpH + Udiscr - Umod
            modelBias += biasEnergyAndDerivs[modelBiasIndex];
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
                for (Residue residue : titratingResidueList) {
                    getBiasTerms(residue, biasEnergyAndDerivs);
                    uComp += biasEnergyAndDerivs[discrBiasIndex];
                }
                return uComp;
            case ModelBias:
                for (Residue residue : titratingResidueList) {
                    getBiasTerms(residue, biasEnergyAndDerivs);
                    //Reminder: Ubias = UpH + Udiscr - Umod
                    uComp += biasEnergyAndDerivs[modelBiasIndex];
                }
                return uComp;
            case pHBias:
                for (Residue residue : titratingResidueList) {
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
     *
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
            case CYS:
                prefactor = titrationLambdas[atomIndex];
                titrationDeriv = 1.0;
                tautomerDeriv = 0.0;
                break;
        }
        vdwPrefactorAndDerivs[0] = prefactor;
        vdwPrefactorAndDerivs[1] = titrationDeriv;
        vdwPrefactorAndDerivs[2] = tautomerDeriv;
    }

    public boolean isTitratingHydrogen(int atomIndex) {
        return isTitratingHydrogen[atomIndex];
    }

    public boolean isTitratingSulfur(int atomIndex){ return isTitratingSulfur[atomIndex];}

    /**
     * Add van der Waals deriv to appropriate dU/dL term.
     *
     * @param atomI
     * @param vdwEnergy             van der Waals energy calculated with no titration/tautomer scaling
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
        int tautomerEsvIndex = tautomerIndexMap[atomI] + nTitr;
        double dTitr_dLambda;
        double dTaut_dLambda;

        dTitr_dLambda = vdwPrefactorAndDerivI[1] * vdwPrefactorJ * vdwEnergy;
        dTaut_dLambda = vdwPrefactorAndDerivI[2] * vdwPrefactorJ * vdwEnergy;

        esvVdwDerivs[titrationEsvIndex].addAndGet(dTitr_dLambda);
        if (tautomerEsvIndex >= nTitr) {
            esvVdwDerivs[tautomerEsvIndex].addAndGet(dTaut_dLambda);
        }
    }

    public void addPermElecDeriv(int atomI, double titrationEnergy, double tautomerEnergy) {
        int titrationEsvIndex = titrationIndexMap[atomI];
        int tautomerEsvIndex = tautomerIndexMap[atomI] + nTitr;
        esvPermElecDerivs[titrationEsvIndex].addAndGet(titrationEnergy);
        if (tautomerEsvIndex >= nTitr) {
            esvPermElecDerivs[tautomerEsvIndex].addAndGet(tautomerEnergy);
        }
    }

    public void addIndElecDeriv(int atomI, double titrationEnergy, double tautomerEnergy) {
        int titrationEsvIndex = titrationIndexMap[atomI];
        int tautomerEsvIndex = tautomerIndexMap[atomI] + nTitr;
        esvIndElecDerivs[titrationEsvIndex].addAndGet(titrationEnergy);
        if (tautomerEsvIndex >= nTitr) {
            esvIndElecDerivs[tautomerEsvIndex].addAndGet(tautomerEnergy);
        }
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
     * @return titrationUtils master copy from ExtendedSystem.
     */
    public TitrationUtils getTitrationUtils() {
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
        for (int i = 0; i < nESVs; i++) {
            dEdTheta[i] = dEdL[i] * sin(2 * thetaPosition[i]);
            //logger.info("dEdL["+i+"]: "+dEdL[i]);
        }
        return dEdTheta;
    }

    public void writeRestart() {
        String esvName = FileUtils.relativePathTo(restartFile).toString();
        if (esvFilter.writeESV(restartFile, thetaPosition, thetaVelocity, thetaAccel)) {
            logger.info(" Wrote dynamics restart file to " + esvName);
        } else {
            logger.info(" Writing dynamics restart file to " + esvName + " failed");
        }
    }

    public void writeLambdaHistogram() {
        StringBuilder tautomerHeader = new StringBuilder("      X→ ");
        for (int k = 0; k < 10; k++) {
            double lb = (double) k / 10;
            double ub = (double) (k + 1) / 10;
            tautomerHeader.append(String.format("%1$10s", "[" + lb + "-" + ub + "]"));
        }
        tautomerHeader.append("\nλ↓");
        for (int i = 0; i < nTitr; i++) {
            logger.info(format("ESV: %s (%d) \n", titratingResidueList.get(i).getName(), i));
            logger.info(tautomerHeader.toString());
            for (int j = 0; j < 10; j++) {
                double lb = (double) j / 10;
                double ub = (j + 1.0) / 10;
                StringBuilder histogram = new StringBuilder();
                for (int k = 0; k < 10; k++) {
                    StringBuilder hisvalue = new StringBuilder();
                    hisvalue.append(String.format("%1$10s", esvHistogram[i][j][k]));
                    histogram.append(hisvalue);
                }
                logger.info("[" + lb + "-" + ub + "]" + histogram);
            }
            logger.info("\n");
        }
    }

    //Naive guess as to what the best starting state should be based purely on the acidostat term.
    private double initialTitrationState(Residue residue, double initialLambda){
        AminoAcid3 AA3 = residue.getAminoAcid3();
        double initialTitrationLambda;
        switch (AA3){
            case ASD:
                initialTitrationLambda = (constantSystemPh < TitrationUtils.Titration.ASHtoASP.pKa) ? 1.0 : 0.0;
                break;
            case GLD:
                initialTitrationLambda = (constantSystemPh < TitrationUtils.Titration.GLHtoGLU.pKa) ? 1.0 : 0.0;
                break;
            case HIS:
                initialTitrationLambda = (constantSystemPh < TitrationUtils.Titration.HIStoHID.pKa) ? 1.0 : 0.0;
                break;
            case LYS:
                initialTitrationLambda = (constantSystemPh < TitrationUtils.Titration.LYStoLYD.pKa) ? 1.0 : 0.0;
                break;
            case CYS:
                initialTitrationLambda = (constantSystemPh < TitrationUtils.Titration.CYStoCYD.pKa) ? 1.0 : 0.0;
                break;
            default:
                initialTitrationLambda = initialLambda;
                break;
        }
        return initialTitrationLambda;
    }

    @Override
    public double energy(double[] x) {
        return 0;
    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {
        return 0;
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        return new double[0];
    }

    @Override
    public double[] getCoordinates(double[] parameters) {
        return new double[0];
    }

    @Override
    public STATE getEnergyTermState() {
        return null;
    }

    @Override
    public void setEnergyTermState(STATE state) {

    }

    @Override
    public double[] getMass() {
        return new double[0];
    }

    @Override
    public int getNumberOfVariables() {
        return 0;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        return new double[0];
    }

    @Override
    public double[] getScaling() {
        return new double[0];
    }

    @Override
    public void setScaling(double[] scaling) {

    }

    @Override
    public double getTotalEnergy() {
        return 0;
    }

    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return new VARIABLE_TYPE[0];
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        return new double[0];
    }

    @Override
    public void setAcceleration(double[] acceleration) {

    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {

    }

    @Override
    public void setVelocity(double[] velocity) {

    }
}
