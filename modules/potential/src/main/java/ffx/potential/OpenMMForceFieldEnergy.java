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
package ffx.potential;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.lang.SystemUtils;
import static org.apache.commons.math3.util.FastMath.sqrt;

import simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod;
import simtk.openmm.OpenMMLibrary;
import simtk.openmm.OpenMMLibrary.OpenMM_Boolean;
import simtk.openmm.OpenMMLibrary.OpenMM_NonbondedForce_NonbondedMethod;
import simtk.openmm.OpenMM_Vec3;

import static simtk.openmm.OpenMMAmoebaLibrary.*;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent12;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent13;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent14;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent15;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent11;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_Bisector;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_NoAxisType;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ThreeFold;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZBisect;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZOnly;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZThenX;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_NoCutoff;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_PME;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Direct;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Extrapolated;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Mutual;
import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static simtk.openmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static simtk.openmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePair;
import static simtk.openmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePairNoExclusions;
import static simtk.openmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_SingleParticle;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.*;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.AngleType.AngleFunction;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.BondType.BondFunction;
import ffx.potential.parameters.ImproperTorsionType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiTorsionType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWType;
import ffx.potential.utils.EnergyException;

import static ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_SIZE.RADIUS;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_TYPE.R_MIN;
import static ffx.potential.nonbonded.VanDerWaalsForm.VDW_TYPE.LENNARD_JONES;
import ffx.potential.parameters.ForceField.ForceFieldDouble;

/**
 * Compute the potential energy and derivatives using OpenMM.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class OpenMMForceFieldEnergy extends ForceFieldEnergy {

    private static final Logger logger = Logger.getLogger(OpenMMForceFieldEnergy.class.getName());

    private static PointerByReference openMMPlatform;

    /**
     * OpenMM System.
     */
    private PointerByReference openMMSystem;
    /**
     * OpenMM Integrator.
     */
    private PointerByReference openMMIntegrator;
    /**
     * OpenMM Context.
     */
    private PointerByReference openMMContext;
    /**
     * OpenMM State.
     */
    private PointerByReference openMMState;

    private PointerByReference openMMForces;
    private PointerByReference openMMPositions;
    private PointerByReference openMMVelocities;
    private PointerByReference thermostat;
    /**
     * OpenMM center-of-mass motion remover.
     */
    private PointerByReference commRemover = null;

    /**
     * OpenMM AMOEBA Force References.
     */
    private PointerByReference amoebaBondForce = null;
    private PointerByReference amoebaAngleForce = null;
    private PointerByReference amoebaInPlaneAngleForce = null;
    private PointerByReference amoebaUreyBradleyForce = null;
    private PointerByReference amoebaOutOfPlaneBendForce = null;
    private PointerByReference amoebaStretchBendForce = null;
    private PointerByReference amoebaTorsionForce = null;
    private PointerByReference amoebaImproperTorsionForce = null;
    private PointerByReference amoebaPiTorsionForce = null;
    private PointerByReference amoebaTorsionTorsionForce = null;
    private PointerByReference amoebaVDWForce = null;
    private PointerByReference amoebaMultipoleForce = null;
    private PointerByReference amoebaGeneralizedKirkwoodForce = null;
    private PointerByReference amoebaWcaDispersionForce = null;
    /**
     * OpenMM Fixed Charge Force References.
     */
    private PointerByReference fixedChargeNonBondedForce = null;
    private PointerByReference customGBForce = null;

    private double lambda = 1.0;

    /**
     * Use flag for each atom.
     */
    private boolean use[];

    /**
     * Size of step to take in lambda for finite differences.
     */
    private double fdDLambda = 0.001;

    /**
     * OpenMMForceFieldEnergy constructor; offloads heavy-duty computation to an
     * OpenMM Platform while keeping track of information locally.
     *
     * @param molecularAssembly Assembly to contruct energy for.
     * @param platform OpenMM platform to be used.
     * @param restraints Harmonic coordinate restraints.
     * @param nThreads Number of threads to use in the underlying
     * ForceFieldEnergy.
     */
    protected OpenMMForceFieldEnergy(MolecularAssembly molecularAssembly, Platform platform, List<CoordRestraint> restraints, int nThreads) {
        super(molecularAssembly, restraints, nThreads);

        //super.energy(false, true);
        logger.info(" Initializing OpenMM\n");

        if (openMMPlatform == null) {
            loadPlatform(platform);
        }

        // Create the OpenMM System
        openMMSystem = OpenMM_System_create();
        logger.info(" Created OpenMM System");

        // Load atoms.
        addAtoms();

        // Add Bond Force.
        addBondForce();

        // Add Angle Force.
        addAngleForce();
        addInPlaneAngleForce();

        // Add Stretch-Bend Force.
        addStretchBendForce();

        // Add Urey-Bradley Force.
        addUreyBradleyForce();

        // Out-of Plane Bend Force.
        addOutOfPlaneBendForce();

        // Add Torsion Force.
        addTorsionForce();

        // Add Improper Torsion Force.
        addImproperTorsionForce();

        // Add Pi-Torsion Force.
        addPiTorsionForce();

        // Add Torsion-Torsion Force.
        addTorsionTorsionForce();

        // Add coordinate restraints.
        addHarmonicRestraintForce();

        VanDerWaals vdW = super.getVdwNode();
        if (vdW != null) {
            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            if (vdwForm.vdwType == LENNARD_JONES) {
                addFixedChargeNonBondedForce();
            } else {
                // Add vdW Force.
                addAmoebaVDWForce();

                // Add Multipole Force.
                addAmoebaMultipoleForce();
            }
        }

        // Set periodic box vectors.
        setDefaultPeriodicBoxVectors();

        openMMIntegrator = OpenMM_VerletIntegrator_create(0.001);

        // Create a openMMContext.
        openMMContext = OpenMM_Context_create_2(openMMSystem, openMMIntegrator, openMMPlatform);

        // Set initial positions.
        loadFFXPositionToOpenMM();

        int infoMask = OpenMM_State_Positions;
        infoMask += OpenMM_State_Forces;
        infoMask += OpenMM_State_Energy;

        openMMState = OpenMM_Context_getState(openMMContext, infoMask, 0);
        openMMForces = OpenMM_State_getForces(openMMState);
        double openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(openMMState) / OpenMM_KJPerKcal;

        logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", openMMPotentialEnergy));
        ForceField forceField = molecularAssembly.getForceField();
        fdDLambda = forceField.getDouble(ForceFieldDouble.FD_DLAMBDA, 0.00001);

        OpenMM_State_destroy(openMMState);
    }

    /**
     * Load an OpenMM Platform
     */
    private void loadPlatform(Platform platform) {
        // Print out the OpenMM Version.
        Pointer version = OpenMM_Platform_getOpenMMVersion();
        logger.log(Level.INFO, " OpenMM Version: {0}", version.getString(0));

        // Print out the OpenMM plugin directory.
        Pointer pluginDir = OpenMM_Platform_getDefaultPluginsDirectory();
        logger.log(Level.INFO, " OpenMM Plugin Dir: {0}", pluginDir.getString(0));

        /**
         * Load plugins and print out plugins.
         *
         * Call the method twice to avoid a bug in OpenMM where not all
         * platforms are list after the first call.
         */
        PointerByReference platforms = OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));
        OpenMM_StringArray_destroy(platforms);
        String plugDirString = pluginDir.getString(0);
        if (SystemUtils.IS_OS_WINDOWS) {
            plugDirString = plugDirString + "/plugins";
        }
        platforms = OpenMM_Platform_loadPluginsFromDirectory(plugDirString);

        int numPlatforms = OpenMM_Platform_getNumPlatforms();
        boolean cuda = false;
        logger.log(Level.INFO, " Number of OpenMM Plugins: {0}", numPlatforms);
        for (int i = 0; i < numPlatforms; i++) {
            String platformString = stringFromArray(platforms, i);
            logger.log(Level.INFO, " Plugin Library :{0}", platformString);
            if (platformString.toUpperCase().contains("AMOEBACUDA")) {
                cuda = true;
            }
        }
        OpenMM_StringArray_destroy(platforms);

        /**
         * Extra logging to print out plugins that failed to load.
         */
        if (logger.isLoggable(Level.FINE)) {
            PointerByReference pluginFailers = OpenMM_Platform_getPluginLoadFailures();
            int numFailures = OpenMM_StringArray_getSize(pluginFailers);
            for (int i = 0; i < numFailures; i++) {
                Pointer message = OpenMM_StringArray_get(pluginFailers, i);
                logger.log(Level.FINE, " Plugin load failure: {0}", message.getString(0));
            }
            OpenMM_StringArray_destroy(pluginFailers);
        }

        if (cuda && platform != Platform.OMM_REF) {
            openMMPlatform = OpenMM_Platform_getPlatformByName("CUDA");
            // OpenMM_Platform_setPropertyDefaultValue(platform, stringPtr("Precision"), stringPtr("mixed"));
            logger.info(" Created OpenMM AMOEBA CUDA Plaform");
        } else {
            openMMPlatform = OpenMM_Platform_getPlatformByName("Reference");
            logger.info(" Created OpenMM AMOEBA Reference Plaform");
        }

    }

    /**
     * Create a JNA Pointer to a String.
     *
     * @param string WARNING: assumes ascii-only string
     * @return pointer.
     */
    private Pointer pointerForString(String string) {
        Pointer pointer = new Memory(string.length() + 1);
        pointer.setString(0, string);
        return pointer;
    }

    /**
     * Returns the platform array as a String
     *
     * @param stringArray
     * @param i
     * @return String
     */
    private String stringFromArray(PointerByReference stringArray, int i) {
        Pointer platformPtr = OpenMM_StringArray_get(stringArray, i);
        if (platformPtr == null) {
            return null;
        }
        return platformPtr.getString(0);
    }

    @Override
    public void destroy() throws Exception {
        super.destroy();
        freeOpenMM();
        logger.info(" Destroyed the Context, Integrator, and OpenMMSystem.");
    }

    @Override
    public void finalize() throws Throwable {
        // Safer to leave super.finalize() in, even though right now that calls Object.finalize().
        logger.info(" OpenMMForceFieldEnergy instance is being finalized.");
        super.finalize();
        if (destroyed) {
            logger.info(String.format(" Finalize called on a destroyed OpenMM ForceFieldEnergy %s", this.toString()));
        } else {
            destroy();
        }
    }

    /**
     * Destroys pointer references to Context, Integrator and System to free up
     * memory.
     */
    private void freeOpenMM() {
        if (openMMContext != null) {
            OpenMM_Context_destroy(openMMContext);
            openMMContext = null;
        }
        if (openMMIntegrator != null) {
            OpenMM_Integrator_destroy(openMMIntegrator);
            openMMIntegrator = null;
        }
        if (openMMSystem != null) {
            OpenMM_System_destroy(openMMSystem);
            openMMSystem = null;
        }
    }

    /**
     * Adds atoms from the molecular assembly to the OpenMM System and reports
     * to the user the number of particles added.
     */
    private void addAtoms() {
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        use = new boolean[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            use[i] = atom.getUse();
            OpenMM_System_addParticle(openMMSystem, atom.getMass());
        }
        logger.log(Level.INFO, " Added particles ({0})", nAtoms);
    }

    /**
     * Adds a center-of-mass motion remover to the Potential. Not advised for
     * anything not running MD using the OpenMM library (i.e.
     * OpenMMMolecularDynamics). Has caused bugs with the FFX MD class.
     */
    public void addCOMMRemover() {
        addCOMMRemover(false);
    }

    /**
     * Adds a center-of-mass motion remover to the Potential. Not advised for
     * anything not running MD using the OpenMM library (i.e.
     * OpenMMMolecularDynamics). Has caused bugs with the FFX MD class.
     *
     * @param addIfDuplicate Add a CCOM remover even if it already exists
     */
    public void addCOMMRemover(boolean addIfDuplicate) {
        if (commRemover == null || addIfDuplicate) {
            if (commRemover != null) {
                logger.warning(" Adding a second center-of-mass remover; this is probably incorrect!");
            }
            int frequency = 100;
            commRemover = OpenMM_CMMotionRemover_create(frequency);
            OpenMM_System_addForce(openMMSystem, commRemover);
            logger.log(Level.INFO, " Added center of mass motion remover (frequency: {0})", frequency);
        } else {
            logger.warning(" Attempted to add a second center-of-mass motion remover when one already exists!");
        }
    }

    private void addBondForce() {
        Bond bonds[] = super.getBonds();
        if (bonds == null || bonds.length < 1) {
            return;
        }
        int nBonds = bonds.length;
        amoebaBondForce = OpenMM_AmoebaBondForce_create();
        double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

        for (int i = 0; i < nBonds; i++) {
            Bond bond = bonds[i];
            int i1 = bond.getAtom(0).getXyzIndex() - 1;
            int i2 = bond.getAtom(1).getXyzIndex() - 1;
            BondType bondType = bond.bondType;
            OpenMM_AmoebaBondForce_addBond(amoebaBondForce, i1, i2,
                    bond.bondType.distance * OpenMM_NmPerAngstrom,
                    kParameterConversion * bondType.forceConstant * BondType.units);

        }

        if (bonds[0].bondType.bondFunction == BondFunction.QUARTIC) {
            OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaBondForce,
                    BondType.cubic / OpenMM_NmPerAngstrom);
            OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaBondForce,
                    BondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
        }

        OpenMM_System_addForce(openMMSystem, amoebaBondForce);
        logger.log(Level.INFO, " Added bonds ({0})", nBonds);
    }

    private void addAngleForce() {
        Angle angles[] = super.getAngles();
        if (angles == null || angles.length < 1) {
            return;
        }
        int nAngles = angles.length;
        List<Angle> normalAngles = new ArrayList<>();
        // Sort all normal angles from in-plane angles
        for (int i = 0; i < nAngles; i++) {
            if (angles[i].getAngleMode() == Angle.AngleMode.NORMAL) {
                normalAngles.add(angles[i]);
            }
        }
        nAngles = normalAngles.size();
        if (nAngles < 1) {
            return;
        }
        amoebaAngleForce = OpenMM_AmoebaAngleForce_create();
        for (int i = 0; i < nAngles; i++) {
            Angle angle = normalAngles.get(i);
            int i1 = angle.getAtom(0).getXyzIndex() - 1;
            int i2 = angle.getAtom(1).getXyzIndex() - 1;
            int i3 = angle.getAtom(2).getXyzIndex() - 1;
            int nh = angle.nh;
            OpenMM_AmoebaAngleForce_addAngle(amoebaAngleForce, i1, i2, i3,
                    angle.angleType.angle[nh], OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant);
        }

        if (angles[0].angleType.angleFunction == AngleFunction.SEXTIC) {
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(amoebaAngleForce, AngleType.cubic);
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(amoebaAngleForce, AngleType.quartic);
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(amoebaAngleForce, AngleType.quintic);
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(amoebaAngleForce, AngleType.sextic);
        }

        OpenMM_System_addForce(openMMSystem, amoebaAngleForce);
        logger.log(Level.INFO, " Added angles ({0})", nAngles);
    }

    private void addInPlaneAngleForce() {
        Angle angles[] = super.getAngles();
        if (angles == null || angles.length < 1) {
            return;
        }
        int nAngles = angles.length;
        List<Angle> inPlaneAngles = new ArrayList<>();
        //Sort all in-plane angles from normal angles
        for (int i = 0; i < nAngles; i++) {
            if (angles[i].getAngleMode() == Angle.AngleMode.IN_PLANE) {
                inPlaneAngles.add(angles[i]);
            }
        }
        nAngles = inPlaneAngles.size();
        if (nAngles < 1) {
            return;
        }
        amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create();
        for (int i = 0; i < nAngles; i++) {
            Angle angle = inPlaneAngles.get(i);
            int i1 = angle.getAtom(0).getXyzIndex() - 1;
            int i2 = angle.getAtom(1).getXyzIndex() - 1;
            int i3 = angle.getAtom(2).getXyzIndex() - 1;
            int i4 = angle.getAtom4().getXyzIndex() - 1;
            int nh = angle.nh;
            OpenMM_AmoebaInPlaneAngleForce_addAngle(amoebaInPlaneAngleForce, i1, i2, i3, i4,
                    angle.angleType.angle[nh], OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant);
        }
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic(amoebaInPlaneAngleForce, AngleType.cubic);
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic(amoebaInPlaneAngleForce, AngleType.quartic);
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic(amoebaInPlaneAngleForce, AngleType.quintic);
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic(amoebaInPlaneAngleForce, AngleType.sextic);
        OpenMM_System_addForce(openMMSystem, amoebaInPlaneAngleForce);
        logger.log(Level.INFO, " Added in-plane angles ({0})", nAngles);
    }

    private void addUreyBradleyForce() {
        UreyBradley ureyBradleys[] = super.getUreyBradleys();
        if (ureyBradleys == null || ureyBradleys.length < 1) {
            return;
        }
        amoebaUreyBradleyForce = OpenMM_AmoebaBondForce_create();
        double kParameterConversion = UreyBradleyType.units * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
        int nUreys = ureyBradleys.length;
        for (int i = 0; i < nUreys; i++) {
            UreyBradley ureyBradley = ureyBradleys[i];
            int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
            int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
            UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
            OpenMM_AmoebaBondForce_addBond(amoebaUreyBradleyForce, i1, i2,
                    ureyBradleyType.distance * OpenMM_NmPerAngstrom,
                    ureyBradleyType.forceConstant * kParameterConversion);
        }

        OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaUreyBradleyForce,
                UreyBradleyType.cubic / OpenMM_NmPerAngstrom);
        OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaUreyBradleyForce,
                UreyBradleyType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));

        OpenMM_System_addForce(openMMSystem, amoebaUreyBradleyForce);
        logger.log(Level.INFO, " Added Urey-Bradleys ({0})", nUreys);
    }

    private void addOutOfPlaneBendForce() {
        OutOfPlaneBend outOfPlaneBends[] = super.getOutOfPlaneBends();
        if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
            return;
        }
        amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create();
        int nOutOfPlaneBends = outOfPlaneBends.length;
        for (int i = 0; i < nOutOfPlaneBends; i++) {
            OutOfPlaneBend outOfPlaneBend = outOfPlaneBends[i];
            int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
            int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
            int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
            int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
            OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;

            OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(amoebaOutOfPlaneBendForce, i1, i2, i3, i4,
                    OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * OutOfPlaneBendType.units);
        }
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.cubic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.quartic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.quintic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.sextic);
        OpenMM_System_addForce(openMMSystem, amoebaOutOfPlaneBendForce);
        logger.log(Level.INFO, " Added Out of Plane Bends ({0})", nOutOfPlaneBends);
    }

    private void addStretchBendForce() {
        StretchBend stretchBends[] = super.getStretchBends();
        if (stretchBends == null || stretchBends.length < 1) {
            return;
        }
        int nStretchBends = stretchBends.length;
        amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create();
        for (int i = 0; i < nStretchBends; i++) {
            StretchBend stretchBend = stretchBends[i];
            int i1 = stretchBend.getAtom(0).getXyzIndex() - 1;
            int i2 = stretchBend.getAtom(1).getXyzIndex() - 1;
            int i3 = stretchBend.getAtom(2).getXyzIndex() - 1;
            double angle = stretchBend.angleEq;
            double beq0 = stretchBend.bond0Eq;
            double beq1 = stretchBend.bond1Eq;
            double fc0 = stretchBend.force0;
            double fc1 = stretchBend.force1;
            OpenMM_AmoebaStretchBendForce_addStretchBend(amoebaStretchBendForce, i1, i2, i3,
                    beq0 * OpenMM_NmPerAngstrom, beq1 * OpenMM_NmPerAngstrom, OpenMM_RadiansPerDegree * angle,
                    (OpenMM_KJPerKcal / OpenMM_NmPerAngstrom) * fc0, (OpenMM_KJPerKcal / OpenMM_NmPerAngstrom) * fc1);

        }
        OpenMM_System_addForce(openMMSystem, amoebaStretchBendForce);
        logger.log(Level.INFO, " Added Stretch Bends ({0})", nStretchBends);
    }

    private void addTorsionForce() {
        Torsion torsions[] = super.getTorsions();
        if (torsions == null || torsions.length < 1) {
            return;
        }
        int nTorsions = torsions.length;
        amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create();
        for (int i = 0; i < nTorsions; i++) {
            Torsion torsion = torsions[i];
            int a1 = torsion.getAtom(0).getXyzIndex() - 1;
            int a2 = torsion.getAtom(1).getXyzIndex() - 1;
            int a3 = torsion.getAtom(2).getXyzIndex() - 1;
            int a4 = torsion.getAtom(3).getXyzIndex() - 1;
            TorsionType torsionType = torsion.torsionType;
            int nTerms = torsionType.phase.length;
            for (int j = 0; j < nTerms; j++) {
                OpenMM_PeriodicTorsionForce_addTorsion(amoebaTorsionForce,
                        a1, a2, a3, a4, j + 1,
                        torsionType.phase[j] * OpenMM_RadiansPerDegree,
                        OpenMM_KJPerKcal * torsion.units * torsionType.amplitude[j]);
            }
        }

        OpenMM_System_addForce(openMMSystem, amoebaTorsionForce);
        logger.log(Level.INFO, " Added Torsions ({0})", nTorsions);
    }

    private void addImproperTorsionForce() {
        ImproperTorsion impropers[] = super.getImproperTorsions();
        if (impropers == null || impropers.length < 1) {
            return;
        }
        int nImpropers = impropers.length;
        amoebaImproperTorsionForce = OpenMM_PeriodicTorsionForce_create();

        for (int i = 0; i < nImpropers; i++) {
            ImproperTorsion improperTorsion = impropers[i];
            int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
            int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
            int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
            int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
            ImproperTorsionType improperTorsionType = improperTorsion.improperType;
            OpenMM_PeriodicTorsionForce_addTorsion(amoebaImproperTorsionForce,
                    a1, a2, a3, a4, improperTorsionType.periodicity,
                    improperTorsionType.phase * OpenMM_RadiansPerDegree,
                    OpenMM_KJPerKcal * improperTorsion.units
                    * improperTorsion.scaleFactor * improperTorsionType.k);
        }
        OpenMM_System_addForce(openMMSystem, amoebaImproperTorsionForce);
        logger.log(Level.INFO, " Added improper torsions ({0})", nImpropers);
    }

    private void addPiTorsionForce() {
        PiOrbitalTorsion piOrbitalTorsions[] = super.getPiOrbitalTorsions();
        if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
            return;
        }
        int nPiOrbitalTorsions = piOrbitalTorsions.length;
        amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create();
        double units = PiTorsionType.units;
        for (int i = 0; i < nPiOrbitalTorsions; i++) {
            PiOrbitalTorsion piOrbitalTorsion = piOrbitalTorsions[i];
            int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
            int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
            int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
            int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
            int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
            int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
            PiTorsionType type = piOrbitalTorsion.piTorsionType;
            OpenMM_AmoebaPiTorsionForce_addPiTorsion(amoebaPiTorsionForce,
                    a1, a2, a3, a4, a5, a6,
                    OpenMM_KJPerKcal * type.forceConstant * units);
        }
        OpenMM_System_addForce(openMMSystem, amoebaPiTorsionForce);
        logger.log(Level.INFO, " Added Pi-Orbital Torsions ({0})", nPiOrbitalTorsions);
    }

    private void addTorsionTorsionForce() {
        TorsionTorsion torsionTorsions[] = super.getTorsionTorsions();
        if (torsionTorsions == null || torsionTorsions.length < 1) {
            return;
        }
        /**
         * Load the torsion-torsions.
         */

        int nTypes = 0;
        LinkedHashMap<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();

        int nTorsionTorsions = torsionTorsions.length;
        amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create();
        for (int i = 0; i < nTorsionTorsions; i++) {
            TorsionTorsion torsionTorsion = torsionTorsions[i];
            int ia = torsionTorsion.getAtom(0).getXyzIndex() - 1;
            int ib = torsionTorsion.getAtom(1).getXyzIndex() - 1;
            int ic = torsionTorsion.getAtom(2).getXyzIndex() - 1;
            int id = torsionTorsion.getAtom(3).getXyzIndex() - 1;
            int ie = torsionTorsion.getAtom(4).getXyzIndex() - 1;

            TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
            String key = torsionTorsionType.getKey();
            /**
             * Check if the TorTor parameters have already been added to the
             * Hash.
             */
            int gridIndex = 0;
            if (torTorTypes.containsKey(key)) {
                /**
                 * If the TorTor has been added, get its (ordered) index in the
                 * Hash.
                 */
                int index = 0;
                for (String entry : torTorTypes.keySet()) {
                    if (entry.equalsIgnoreCase(key)) {
                        gridIndex = index;
                        break;
                    } else {
                        index++;
                    }
                }
            } else {
                /**
                 * Add the new TorTor.
                 */
                torTorTypes.put(key, torsionTorsionType);
                gridIndex = nTypes;
                nTypes++;
            }

            Atom atom = torsionTorsion.getChiralAtom();
            int iChiral = -1;
            if (atom != null) {
                iChiral = atom.getXyzIndex() - 1;
            }
            OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(amoebaTorsionTorsionForce,
                    ia, ib, ic, id, ie, iChiral, gridIndex);
        }
        /**
         * Load the Torsion-Torsion parameters.
         */
        PointerByReference values = OpenMM_DoubleArray_create(6);
        int gridIndex = 0;
        for (String key : torTorTypes.keySet()) {
            TorsionTorsionType torTorType = torTorTypes.get(key);
            int nx = torTorType.nx;
            int ny = torTorType.ny;
            double tx[] = torTorType.tx;
            double ty[] = torTorType.ty;
            double f[] = torTorType.energy;
            double dx[] = torTorType.dx;
            double dy[] = torTorType.dy;
            double dxy[] = torTorType.dxy;
            /**
             * Create the 3D grid.
             */
            PointerByReference grid3D = OpenMM_3D_DoubleArray_create(nx, ny, 6);
            int xIndex = 0;
            int yIndex = 0;
            for (int j = 0; j < nx * ny; j++) {
                int addIndex = 0;
                OpenMM_DoubleArray_set(values, addIndex++, tx[xIndex]);
                OpenMM_DoubleArray_set(values, addIndex++, ty[yIndex]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * f[j]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dx[j]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dy[j]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dxy[j]);
                OpenMM_3D_DoubleArray_set(grid3D, yIndex, xIndex, values);
                xIndex++;
                if (xIndex == nx) {
                    xIndex = 0;
                    yIndex++;
                }
            }
            OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(amoebaTorsionTorsionForce, gridIndex++, grid3D);
            OpenMM_3D_DoubleArray_destroy(grid3D);
        }
        OpenMM_DoubleArray_destroy(values);
        OpenMM_System_addForce(openMMSystem, amoebaTorsionTorsionForce);
        logger.log(Level.INFO, " Added Torsion-Torsions ({0})", nTorsionTorsions);
    }

    /**
     * Uses arithmetic mean to define sigma and geometric mean for epsilon.
     */
    private void addFixedChargeNonBondedForce() {
        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return;
        }
        /**
         * Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
         * for epsilon is supported.
         */
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        if (vdwForm.vdwType != LENNARD_JONES
                || vdwForm.radiusRule != ARITHMETIC
                || vdwForm.epsilonRule != GEOMETRIC) {
            logger.log(Level.SEVERE, String.format(" Unsuppporterd van der Waals functional form."));
            return;
        }

        fixedChargeNonBondedForce = OpenMM_NonbondedForce_create();

        /**
         * OpenMM vdW force requires a diameter (i.e. not radius).
         */
        double radScale = 1.0;
        if (vdwForm.radiusSize == RADIUS) {
            radScale = 2.0;
        }
        /**
         * OpenMM vdw force requires atomic sigma values (i.e. not r-min).
         */
        if (vdwForm.radiusType == R_MIN) {
            radScale /= 1.122462048309372981;
        }

        /**
         * Add particles.
         */
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            double charge = 0.0;
            MultipoleType multipoleType = atom.getMultipoleType();
            if (multipoleType != null) {
                charge = multipoleType.charge;
            }
            VDWType vdwType = atom.getVDWType();
            double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
            double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
            OpenMM_NonbondedForce_addParticle(fixedChargeNonBondedForce, charge, sigma, eps);
        }
        /**
         * Define 1-4 scale factors.
         */
        double lj14Scale = vdwForm.getScale14();
        double coulomb14Scale = 1.0 / 1.2;

        ParticleMeshEwald pme = super.getPmeNode();
        Bond bonds[] = super.getBonds();
        if (bonds != null && bonds.length > 0) {
            int nBonds = bonds.length;
            PointerByReference bondArray;
            bondArray = OpenMM_BondArray_create(0);
            for (int i = 0; i < nBonds; i++) {
                Bond bond = bonds[i];
                int i1 = bond.getAtom(0).getXyzIndex() - 1;
                int i2 = bond.getAtom(1).getXyzIndex() - 1;
                OpenMM_BondArray_append(bondArray, i1, i2);
            }
            if (pme != null) {
                coulomb14Scale = pme.getScale14();
            }
            OpenMM_NonbondedForce_createExceptionsFromBonds(fixedChargeNonBondedForce, bondArray, coulomb14Scale, lj14Scale);
            OpenMM_BondArray_destroy(bondArray);
        }

        Crystal crystal = super.getCrystal();
        if (crystal.aperiodic()) {
            OpenMM_NonbondedForce_setNonbondedMethod(fixedChargeNonBondedForce,
                    OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_NoCutoff);
        } else {
            OpenMM_NonbondedForce_setNonbondedMethod(fixedChargeNonBondedForce,
                    OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_PME);
            if (pme != null) {
                // Units of the Ewald coefficient are A^-1; Multiply by AngstromsPerNM to convert to (Nm^-1).
                double aEwald = OpenMM_AngstromsPerNm * pme.getEwaldCoefficient();
                int nx = pme.getReciprocalSpace().getXDim();
                int ny = pme.getReciprocalSpace().getYDim();
                int nz = pme.getReciprocalSpace().getZDim();
                OpenMM_NonbondedForce_setPMEParameters(fixedChargeNonBondedForce, aEwald, nx, ny, nz);
            }
        }

        NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
        double off = nonbondedCutoff.off;
        double cut = nonbondedCutoff.cut;
        OpenMM_NonbondedForce_setCutoffDistance(fixedChargeNonBondedForce, OpenMM_NmPerAngstrom * off);

        OpenMM_NonbondedForce_setUseSwitchingFunction(fixedChargeNonBondedForce, OpenMM_True);
        if (cut == off) {
            logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
            if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
                logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
                cut = 1E40;
            } else {
                logger.info(String.format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut, off));
                cut *= 0.99;
            }
        }
        OpenMM_NonbondedForce_setSwitchingDistance(fixedChargeNonBondedForce, OpenMM_NmPerAngstrom * cut);

        OpenMM_NonbondedForce_setUseDispersionCorrection(fixedChargeNonBondedForce, OpenMM_False);

        OpenMM_Force_setForceGroup(fixedChargeNonBondedForce, 1);
        OpenMM_System_addForce(openMMSystem, fixedChargeNonBondedForce);
        logger.log(Level.INFO, String.format(" Added fixed charge non-bonded force."));

        GeneralizedKirkwood gk = super.getGK();
        if (gk != null) {
            addCustomGBForce();
        }
    }

    private void addCustomGBForce() {
        GeneralizedKirkwood gk = super.getGK();
        if (gk == null) {
            return;
        }

        customGBForce = OpenMM_CustomGBForce_create();
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "q");
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "radius");
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "scale");
        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "solventDielectric", 78.3);
        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "soluteDielectric", 1.0);
        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "dOffset", gk.getDielecOffset() * OpenMMLibrary.OpenMM_NmPerAngstrom); // Factor of 0.1 for Ang to nm.
        OpenMM_CustomGBForce_addComputedValue(customGBForce, "I",
                // "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                // "step(r+sr2-or1)*0.5*((1/L^3-1/U^3)/3+(1/U^4-1/L^4)/8*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);"
                "0.5*((1/L^3-1/U^3)/3.0+(1/U^4-1/L^4)/8.0*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);"
                + "U=r+sr2;"
                // + "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                + "C=2/3*(1/or1^3-1/L^3)*step(sr2-r-or1);"
                // + "L=step(or1-D)*or1 + (1-step(or1-D))*D;"
                // + "D=step(r-sr2)*(r-sr2) + (1-step(r-sr2))*(sr2-r);"
                + "L = step(sr2 - r1r)*sr2mr + (1 - step(sr2 - r1r))*L;"
                + "sr2mr = sr2 - r;"
                + "r1r = radius1 + r;"
                + "L = step(r1sr2 - r)*radius1 + (1 - step(r1sr2 - r))*L;"
                + "r1sr2 = radius1 + sr2;"
                + "L = r - sr2;"
                + "sr2 = scale2 * radius2;"
                + "or1 = radius1; or2 = radius2",
                OpenMM_CustomGBForce_ParticlePairNoExclusions);

        OpenMM_CustomGBForce_addComputedValue(customGBForce, "B",
                // "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                // "psi=I*or; or=radius-0.009"
                "step(BB-radius)*BB + (1 - step(BB-radius))*radius;"
                + "BB = 1 / ( (3.0*III)^(1.0/3.0) );"
                + "III = step(II)*II + (1 - step(II))*1.0e-9/3.0;"
                + "II = maxI - I;"
                + "maxI = 1/(3.0*radius^3)",
                OpenMM_CustomGBForce_SingleParticle);

        double sTens = gk.getSurfaceTension();
        logger.info(String.format(" FFX surface tension: %9.5g kcal/mol/Ang^2", sTens));
        sTens *= OpenMMLibrary.OpenMM_KJPerKcal;
        sTens *= 100.0; // 100 square Angstroms per square nanometer.
        logger.info(String.format(" OpenMM surface tension: %9.5g kJ/mol/nm^2", sTens));
        String surfaceTension = Double.toString(sTens);

        OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
                surfaceTension
                + "*(radius+0.14+dOffset)^2*((radius+dOffset)/B)^6/6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B",
                OpenMM_CustomGBForce_SingleParticle);

        /**
         * Particle pair term is the generalized Born cross term.
         */
        OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
                "-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                + "f=sqrt(r^2+B1*B2*exp(-r^2/(2.455*B1*B2)))",
                OpenMM_CustomGBForce_ParticlePair);

        double baseRadii[] = gk.getBaseRadii();
        double overlapScale[] = gk.getOverlapScale();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        PointerByReference doubleArray = OpenMM_DoubleArray_create(0);
        for (int i = 0; i < nAtoms; i++) {
            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_DoubleArray_append(doubleArray, multipoleType.charge);
            OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadii[i]);
            OpenMM_DoubleArray_append(doubleArray, overlapScale[i]);
            OpenMM_CustomGBForce_addParticle(customGBForce, doubleArray);
            OpenMM_DoubleArray_resize(doubleArray, 0);
        }
        OpenMM_DoubleArray_destroy(doubleArray);

        double cut = gk.getCutoff();
        OpenMM_CustomGBForce_setCutoffDistance(customGBForce, cut);
        OpenMM_Force_setForceGroup(customGBForce, 1);
        OpenMM_System_addForce(openMMSystem, customGBForce);

        logger.log(Level.INFO, " Added generalized Born force");
    }

    private void addAmoebaVDWForce() {
        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return;
        }

        amoebaVDWForce = OpenMM_AmoebaVdwForce_create();
        OpenMM_System_addForce(openMMSystem, amoebaVDWForce);
        OpenMM_Force_setForceGroup(amoebaVDWForce, 1);

        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
        Crystal crystal = super.getCrystal();

        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        /**
         * Note that the API says it wants a SIGMA value.
         */
        if (vdwForm.radiusType == VanDerWaalsForm.RADIUS_TYPE.R_MIN) {
            //radScale *= 1.122462048309372981;
        }

        int ired[] = vdW.getReductionIndex();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            OpenMM_AmoebaVdwForce_addParticle(amoebaVDWForce,
                    ired[i], OpenMM_NmPerAngstrom * vdwType.radius * radScale,
                    OpenMM_KJPerKcal * vdwType.wellDepth,
                    vdwType.reductionFactor);
        }

        // OpenMM_AmoebaVdwForce_setSigmaCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.radiusRule.name()));
        // OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.epsilonRule.name()));
        OpenMM_AmoebaVdwForce_setCutoffDistance(amoebaVDWForce, nonbondedCutoff.off * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_Boolean.OpenMM_False);

        if (crystal.aperiodic()) {
            OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
                    OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff);
        } else {
            OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
                    OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_CutoffPeriodic);
        }

        /**
         * Create exclusion lists.
         */
        PointerByReference exclusions = OpenMM_IntArray_create(0);
        double mask[] = new double[nAtoms];
        Arrays.fill(mask, 1.0);
        for (int i = 0; i < nAtoms; i++) {
            OpenMM_IntArray_append(exclusions, i);
            vdW.applyMask(mask, i);
            for (int j = 0; j < nAtoms; j++) {
                if (mask[j] == 0.0) {
                    OpenMM_IntArray_append(exclusions, j);
                }
            }
            vdW.removeMask(mask, i);
            OpenMM_AmoebaVdwForce_setParticleExclusions(amoebaVDWForce, i, exclusions);
            OpenMM_IntArray_resize(exclusions, 0);
        }
        OpenMM_IntArray_destroy(exclusions);
        logger.log(Level.INFO, " Added van der Waals force.");
    }

    private void addAmoebaMultipoleForce() {
        ParticleMeshEwald pme = super.getPmeNode();
        if (pme == null) {
            return;
        }
        int axisAtom[][] = pme.getAxisAtoms();
        double dipoleConversion = OpenMM_NmPerAngstrom;
        double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
        double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                * OpenMM_NmPerAngstrom;
        double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

        amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();

        double polarScale = 1.0;
        if (pme.getPolarizationType() != Polarization.MUTUAL) {
            OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Direct);
            if (pme.getPolarizationType() == Polarization.NONE) {
                polarScale = 0.0;
            }
        } else {
            ForceField forceField = molecularAssembly.getForceField();
            String algorithm = forceField.getString(ForceField.ForceFieldString.SCF_ALGORITHM, "CG");
            ParticleMeshEwald.SCFAlgorithm scfAlgorithm;

            try {
                algorithm = algorithm.replaceAll("-", "_").toUpperCase();
                scfAlgorithm = ParticleMeshEwald.SCFAlgorithm.valueOf(algorithm);
            } catch (Exception e) {
                scfAlgorithm = ParticleMeshEwald.SCFAlgorithm.CG;
            }

            switch (scfAlgorithm) {
                case EPT:
                    logger.info(" Using extrapolated perturbation theory approximation instead of full SCF calculations. Not supported in FFX reference implementation.");
                    OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Extrapolated);
                    PointerByReference exptCoefficients = OpenMM_DoubleArray_create(4);
                    OpenMM_DoubleArray_set(exptCoefficients, 0, -0.154);
                    OpenMM_DoubleArray_set(exptCoefficients, 1, 0.017);
                    OpenMM_DoubleArray_set(exptCoefficients, 2, 0.657);
                    OpenMM_DoubleArray_set(exptCoefficients, 3, 0.475);
                    OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(amoebaMultipoleForce,exptCoefficients);
                    OpenMM_DoubleArray_destroy(exptCoefficients);
                    break;
                case CG:
                case SOR:
                default:
                    OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Mutual);
                    break;
            }
        }

        PointerByReference dipoles = OpenMM_DoubleArray_create(3);
        PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            MultipoleType multipoleType = atom.getMultipoleType();
            PolarizeType polarType = atom.getPolarizeType();

            /**
             * Define the frame definition.
             */
            int axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
            switch (multipoleType.frameDefinition) {
                case ZONLY:
                    axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
                    break;
                case ZTHENX:
                    axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
                    break;
                case BISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_Bisector;
                    break;
                case ZTHENBISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
                    break;
                case TRISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                    break;
                default:
                    break;
            }

            /**
             * Load local multipole coefficients.
             */
            for (int j = 0; j < 3; j++) {
                OpenMM_DoubleArray_set(dipoles, j, multipoleType.dipole[j] * dipoleConversion);

            }
            int l = 0;
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    OpenMM_DoubleArray_set(quadrupoles, l++, multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0);
                }
            }

            int zaxis = 0;
            int xaxis = 0;
            int yaxis = 0;
            int refAtoms[] = axisAtom[i];
            if (refAtoms != null) {
                zaxis = refAtoms[0];
                if (refAtoms.length > 1) {
                    xaxis = refAtoms[1];
                    if (refAtoms.length > 2) {
                        yaxis = refAtoms[2];
                    }
                }
            }

            /**
             * Add the multipole.
             */
            OpenMM_AmoebaMultipoleForce_addMultipole(amoebaMultipoleForce,
                    multipoleType.charge, dipoles, quadrupoles,
                    axisType, zaxis, xaxis, yaxis,
                    polarType.thole,
                    polarType.pdamp * dampingFactorConversion,
                    polarType.polarizability * polarityConversion * polarScale);
        }
        OpenMM_DoubleArray_destroy(dipoles);
        OpenMM_DoubleArray_destroy(quadrupoles);

        Crystal crystal = super.getCrystal();
        if (!crystal.aperiodic()) {
            OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_PME);
            OpenMM_AmoebaMultipoleForce_setCutoffDistance(amoebaMultipoleForce,
                    pme.getEwaldCutoff() * OpenMM_NmPerAngstrom);
            OpenMM_AmoebaMultipoleForce_setAEwald(amoebaMultipoleForce,
                    pme.getEwaldCoefficient() / OpenMM_NmPerAngstrom);

            double ewaldTolerance = 1.0e-04;
            OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(amoebaMultipoleForce, ewaldTolerance);

            PointerByReference gridDimensions = OpenMM_IntArray_create(3);
            ReciprocalSpace recip = pme.getReciprocalSpace();
            OpenMM_IntArray_set(gridDimensions, 0, recip.getXDim());
            OpenMM_IntArray_set(gridDimensions, 1, recip.getYDim());
            OpenMM_IntArray_set(gridDimensions, 2, recip.getZDim());
            OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(amoebaMultipoleForce, gridDimensions);
            OpenMM_IntArray_destroy(gridDimensions);
        } else {
            OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_NoCutoff);
        }

        OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(amoebaMultipoleForce, 500);
        OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(amoebaMultipoleForce, pme.getPolarEps());

        int ip11[][] = pme.getPolarization11();
        int ip12[][] = pme.getPolarization12();
        int ip13[][] = pme.getPolarization13();

        ArrayList<Integer> list12 = new ArrayList<>();
        ArrayList<Integer> list13 = new ArrayList<>();
        ArrayList<Integer> list14 = new ArrayList<>();

        PointerByReference covalentMap = OpenMM_IntArray_create(0);
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            list12.clear();
            list13.clear();
            list14.clear();

            for (Bond bond : ai.getBonds()) {
                int index = bond.get1_2(ai).getIndex() - 1;
                OpenMM_IntArray_append(covalentMap, index);
                list12.add(index);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Angle angle : ai.getAngles()) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    int index = ak.getIndex() - 1;
                    if (!list12.contains(index)) {
                        list13.add(index);
                        OpenMM_IntArray_append(covalentMap, index);
                    }
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Torsion torsion : ai.getTorsions()) {
                Atom ak = torsion.get1_4(ai);
                if (ak != null) {
                    int index = ak.getIndex() - 1;
                    if (!list12.contains(index)
                            && !list13.contains(index)) {
                        list14.add(index);
                        OpenMM_IntArray_append(covalentMap, index);
                    }
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Atom ak : ai.get1_5s()) {
                int index = ak.getIndex() - 1;
                if (!list12.contains(index)
                        && !list13.contains(index)
                        && !list14.contains(index)) {
                    OpenMM_IntArray_append(covalentMap, index);
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent15, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (int j = 0; j < ip11[i].length; j++) {
                OpenMM_IntArray_append(covalentMap, ip11[i][j]);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent11, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

//            for (int j = 0; j < ip12[i].length; j++) {
//                OpenMM_IntArray_append(covalentMap, ip12[i][j]);
//            }
//            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
//                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent12, covalentMap);
//            OpenMM_IntArray_resize(covalentMap, 0);
//
//            for (int j = 0; j < ip13[i].length; j++) {
//                OpenMM_IntArray_append(covalentMap, ip13[i][j]);
//            }
//            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
//                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent13, covalentMap);
//            OpenMM_IntArray_resize(covalentMap, 0);
//
//            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
//                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent14, covalentMap);
        }

        OpenMM_IntArray_destroy(covalentMap);

        OpenMM_System_addForce(openMMSystem, amoebaMultipoleForce);
        OpenMM_Force_setForceGroup(amoebaMultipoleForce, 1);

        logger.log(Level.INFO, " Added polarizable multipole force.");

        GeneralizedKirkwood gk = super.getGK();
        if (gk != null) {
            addGKForce();
        }

    }

    private void addGKForce() {

        GeneralizedKirkwood gk = super.getGK();

        amoebaGeneralizedKirkwoodForce = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(amoebaGeneralizedKirkwoodForce, 78.3);
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(amoebaGeneralizedKirkwoodForce, 1.0);

        double overlapScale[] = gk.getOverlapScale();
        double baseRadii[] = gk.getBaseRadii();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(amoebaGeneralizedKirkwoodForce,
                    multipoleType.charge, OpenMM_NmPerAngstrom * baseRadii[i], overlapScale[i]);
        }

        OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(amoebaGeneralizedKirkwoodForce, 1.4 * OpenMM_NmPerAngstrom);

        NonPolar nonpolar = gk.getNonPolarModel();
        switch (nonpolar) {
            case BORN_SOLV:
            case BORN_CAV_DISP:
            default:
                // Configure a Born Radii based surface area term.
                double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal
                        * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm;
                OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(amoebaGeneralizedKirkwoodForce, OpenMM_True);
                OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(amoebaGeneralizedKirkwoodForce, -surfaceTension);
                break;
            case CAV:
            case CAV_DISP:
            case HYDROPHOBIC_PMF:
            case NONE:
                // This NonPolar model does not use a Born Radii based surface area term.
                OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(amoebaGeneralizedKirkwoodForce, OpenMM_False);
                break;
        }
        OpenMM_System_addForce(openMMSystem, amoebaGeneralizedKirkwoodForce);

        switch (nonpolar) {
            case CAV_DISP:
            case BORN_CAV_DISP:
                addWCAForce();
                break;
            case CAV:
            case HYDROPHOBIC_PMF:
            case BORN_SOLV:
            case NONE:
            default:
            // WCA force is not being used.
        }

        logger.log(Level.INFO, " Added generalized Kirkwood force.");
    }

    private void addWCAForce() {

        double epso = 0.1100;
        double epsh = 0.0135;
        double rmino = 1.7025;
        double rminh = 1.3275;
        double awater = 0.033428;
        double slevy = 1.0;
        double dispoff = 0.26;
        double shctd = 0.81;

        VanDerWaals vdW = super.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create();

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;

        for (int i = 0; i < nAtoms; i++) {
            // cdispTotal += nonpol__.cdisp[ii];
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            double radius = vdwType.radius;
            double eps = vdwType.wellDepth;
            OpenMM_AmoebaWcaDispersionForce_addParticle(amoebaWcaDispersionForce,
                    OpenMM_NmPerAngstrom * radius * radScale,
                    OpenMM_KJPerKcal * eps);
        }

        OpenMM_AmoebaWcaDispersionForce_setEpso(amoebaWcaDispersionForce, epso * OpenMM_KJPerKcal);
        OpenMM_AmoebaWcaDispersionForce_setEpsh(amoebaWcaDispersionForce, epsh * OpenMM_KJPerKcal);
        OpenMM_AmoebaWcaDispersionForce_setRmino(amoebaWcaDispersionForce, rmino * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaWcaDispersionForce_setRminh(amoebaWcaDispersionForce, rminh * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaWcaDispersionForce_setDispoff(amoebaWcaDispersionForce, dispoff * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaWcaDispersionForce_setAwater(amoebaWcaDispersionForce,
                awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
        OpenMM_AmoebaWcaDispersionForce_setSlevy(amoebaWcaDispersionForce, slevy);
        OpenMM_AmoebaWcaDispersionForce_setShctd(amoebaWcaDispersionForce, shctd);

        OpenMM_System_addForce(openMMSystem, amoebaWcaDispersionForce);
        logger.log(Level.INFO, " Added WCA dispersion force.");

    }

    /**
     * Adds harmonic restraints (CoordRestraint objects) to OpenMM as a custom
     * external force.
     */
    private void addHarmonicRestraintForce() {
        for (CoordRestraint restraint : super.getCoordRestraints()) {
            double forceConst = restraint.getForceConstant();
            forceConst *= OpenMMLibrary.OpenMM_KJPerKcal;
            forceConst *= (OpenMMLibrary.OpenMM_AngstromsPerNm * OpenMMLibrary.OpenMM_AngstromsPerNm);
            Atom[] restAtoms = restraint.getAtoms();
            int nRestAts = restraint.getNumAtoms();
            double[][] oCoords = restraint.getOriginalCoordinates();
            for (int i = 0; i < nRestAts; i++) {
                oCoords[i][0] *= OpenMMLibrary.OpenMM_NmPerAngstrom;
                oCoords[i][1] *= OpenMMLibrary.OpenMM_NmPerAngstrom;
                oCoords[i][2] *= OpenMMLibrary.OpenMM_NmPerAngstrom;
            }

            PointerByReference theRestraint = OpenMM_CustomExternalForce_create("k*periodicdistance(x,y,z,x0,y0,z0)^2");
            OpenMM_CustomExternalForce_addGlobalParameter(theRestraint, "k", forceConst);
            OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "x0");
            OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "y0");
            OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "z0");

            PointerByReference xyzOrigArray = OpenMM_DoubleArray_create(3);
            for (int i = 0; i < nRestAts; i++) {
                int ommIndex = restAtoms[i].getXyzIndex() - 1;
                for (int j = 0; j < 3; j++) {
                    OpenMM_DoubleArray_set(xyzOrigArray, j, oCoords[i][j]);
                }
                OpenMM_CustomExternalForce_addParticle(theRestraint, ommIndex, xyzOrigArray);
            }
            OpenMM_DoubleArray_destroy(xyzOrigArray);

            OpenMM_System_addForce(openMMSystem, theRestraint);
        }
    }

    /**
     * Update parameters if the Use flags changed.
     */
    private void updateParameters() {

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        List<Integer> changedIndices = new ArrayList<>();
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            boolean useI = atom.getUse();
            if (useI != use[i]) {
                changedIndices.add(i);
                use[i] = useI;
            }
        }

        if (changedIndices.isEmpty()) {
            return;
        }
        int[] useChanged = changedIndices.stream().mapToInt(Integer::intValue).toArray();

        // Update fixed charge non-bonded parameters.
        if (fixedChargeNonBondedForce != null) {
            updateFixedChargeNonBondedForce(atoms, useChanged);
        }

        // Update fixed charge GB parameters.
        if (customGBForce != null) {
            updateCustomGBForce(atoms, useChanged);
        }

        // Update AMOEBA vdW parameters.
        if (amoebaVDWForce != null) {
            updateAmoebaVDWForce(atoms, useChanged);
        }

        // Update AMOEBA polarizable multipole parameters.
        if (amoebaMultipoleForce != null) {
            updateAmoebaMultipoleForce(atoms, useChanged);
        }

        // Update GK force.
        if (amoebaGeneralizedKirkwoodForce != null) {
            updateAmoebaGeneralizedKirkwoodForce(atoms, useChanged);
        }

        // Update WCA Force.
        if (amoebaWcaDispersionForce != null) {
            updateWCAForce(atoms, useChanged);
        }
    }

    /**
     * Updates the AMOEBA van der Waals force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     * @param useChanged Indices of atoms with changed Use flags
     */
    private void updateAmoebaVDWForce(Atom[] atoms, int[] useChanged) {
        VanDerWaals vdW = super.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();

        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        /**
         * Note that the API says it wants a SIGMA value.
         */
        if (vdwForm.radiusType == VanDerWaalsForm.RADIUS_TYPE.R_MIN) {
            //radScale *= 1.122462048309372981;
        }

        int ired[] = vdW.getReductionIndex();
        /*Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {*/
        int nAtoms = useChanged.length;
        for (int j = 0; j < nAtoms; j++) {
            int i = useChanged[j];
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            double useFactor = 1.0;
            if (!use[i]) {
                useFactor = 0.0;
            }
            double eps = OpenMM_KJPerKcal * vdwType.wellDepth * useFactor;
            OpenMM_AmoebaVdwForce_setParticleParameters(amoebaVDWForce,
                    i, ired[i], OpenMM_NmPerAngstrom * vdwType.radius * radScale,
                    eps, vdwType.reductionFactor);
        }
        OpenMM_AmoebaVdwForce_updateParametersInContext(amoebaVDWForce, openMMContext);
    }

    /**
     * Updates the fixed-charge non-bonded force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     * @param useChanged Indices of atoms with changed Use flags
     */
    private void updateFixedChargeNonBondedForce(Atom[] atoms, int[] useChanged) {
        VanDerWaals vdW = super.getVdwNode();
        /**
         * Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
         * for epsilon is supported.
         */
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        if (vdwForm.vdwType != LENNARD_JONES
                || vdwForm.radiusRule != ARITHMETIC
                || vdwForm.epsilonRule != GEOMETRIC) {
            logger.log(Level.SEVERE, String.format(" Unsuppporterd van der Waals functional form."));
            return;
        }

        /**
         * OpenMM vdW force requires a diameter (i.e. not radius).
         */
        double radScale = 1.0;
        if (vdwForm.radiusSize == RADIUS) {
            radScale = 2.0;
        }
        /**
         * OpenMM vdw force requires atomic sigma values (i.e. not r-min).
         */
        if (vdwForm.radiusType == R_MIN) {
            radScale /= 1.122462048309372981;
        }

        /**
         * Update parameters.
         */
        /*Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {*/
        int nAtoms = useChanged.length;
        for (int j = 0; j < nAtoms; j++) {
            int i = useChanged[j];
            Atom atom = atoms[i];
            double useFactor = 1.0;
            if (!use[i]) {
                useFactor = 0.0;
            }
            double charge = 0.0;
            MultipoleType multipoleType = atom.getMultipoleType();
            if (multipoleType != null) {
                charge = multipoleType.charge * useFactor;
            }
            VDWType vdwType = atom.getVDWType();
            double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
            double eps = OpenMM_KJPerKcal * vdwType.wellDepth * useFactor;
            OpenMM_NonbondedForce_setParticleParameters(fixedChargeNonBondedForce, i, charge, sigma, eps);
        }
        OpenMM_NonbondedForce_updateParametersInContext(fixedChargeNonBondedForce, openMMContext);
    }

    /**
     * Updates the custom GB force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     * @param useChanged Indices of atoms with changed Use flags
     */
    private void updateCustomGBForce(Atom[] atoms, int[] useChanged) {
        GeneralizedKirkwood gk = super.getGK();
        double baseRadii[] = gk.getBaseRadii();
        double overlapScale[] = gk.getOverlapScale();
        PointerByReference doubleArray = OpenMM_DoubleArray_create(0);

        /*Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {*/
        int nAtoms = useChanged.length;
        for (int j = 0; j < nAtoms; j++) {
            int i = useChanged[j];
            Atom atom = atoms[i];
            double useFactor = 1.0;
            if (!use[i]) {
                useFactor = 0.0;
            }
            MultipoleType multipoleType = atom.getMultipoleType();
            double charge = multipoleType.charge * useFactor;
            double oScale = overlapScale[i] * useFactor;
            OpenMM_DoubleArray_append(doubleArray, charge);
            OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadii[i]);
            OpenMM_DoubleArray_append(doubleArray, oScale);
            OpenMM_CustomGBForce_setParticleParameters(customGBForce, i, doubleArray);
            OpenMM_DoubleArray_resize(doubleArray, 0);
        }
        OpenMM_DoubleArray_destroy(doubleArray);
        OpenMM_CustomGBForce_updateParametersInContext(customGBForce, openMMContext);
    }

    /**
     * Updates the Amoeba electrostatic multipolar force for change in Use
     * flags.
     *
     * @param atoms Array of all Atoms in the system
     * @param useChanged Indices of atoms with changed Use flags
     */
    private void updateAmoebaMultipoleForce(Atom[] atoms, int[] useChanged) {
        ParticleMeshEwald pme = super.getPmeNode();
        int axisAtom[][] = pme.getAxisAtoms();
        double dipoleConversion = OpenMM_NmPerAngstrom;
        double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
        double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                * OpenMM_NmPerAngstrom;
        double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

        double polarScale = 1.0;
        if (pme.getPolarizationType() == Polarization.NONE) {
            polarScale = 0.0;
        }

        PointerByReference dipoles = OpenMM_DoubleArray_create(3);
        PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

        /*Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {*/
        int nAtoms = useChanged.length;
        for (int iInd = 0; iInd < nAtoms; iInd++) {
            int i = useChanged[iInd];
            Atom atom = atoms[i];
            MultipoleType multipoleType = atom.getMultipoleType();
            PolarizeType polarType = atom.getPolarizeType();

            double useFactor = 1.0;
            if (!use[i]) {
                useFactor = 0.0;
            }

            /**
             * Define the frame definition.
             */
            int axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
            switch (multipoleType.frameDefinition) {
                case ZONLY:
                    axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
                    break;
                case ZTHENX:
                    axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
                    break;
                case BISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_Bisector;
                    break;
                case ZTHENBISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
                    break;
                case TRISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                    break;
                default:
                    break;
            }

            /**
             * Load local multipole coefficients.
             */
            for (int j = 0; j < 3; j++) {
                OpenMM_DoubleArray_set(dipoles, j, multipoleType.dipole[j] * dipoleConversion * useFactor);

            }
            int l = 0;
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    OpenMM_DoubleArray_set(quadrupoles, l++, multipoleType.quadrupole[j][k]
                            * quadrupoleConversion / 3.0 * useFactor);
                }
            }

            int zaxis = 0;
            int xaxis = 0;
            int yaxis = 0;
            int refAtoms[] = axisAtom[i];
            if (refAtoms != null) {
                zaxis = refAtoms[0];
                if (refAtoms.length > 1) {
                    xaxis = refAtoms[1];
                    if (refAtoms.length > 2) {
                        yaxis = refAtoms[2];
                    }
                }
            }

            /**
             * Add the multipole.
             */
            OpenMM_AmoebaMultipoleForce_setMultipoleParameters(amoebaMultipoleForce, i,
                    multipoleType.charge * useFactor, dipoles, quadrupoles,
                    axisType, zaxis, xaxis, yaxis,
                    polarType.thole,
                    polarType.pdamp * dampingFactorConversion,
                    polarType.polarizability * polarityConversion * polarScale * useFactor);
        }
        OpenMM_DoubleArray_destroy(dipoles);
        OpenMM_DoubleArray_destroy(quadrupoles);

        OpenMM_AmoebaMultipoleForce_updateParametersInContext(amoebaMultipoleForce, openMMContext);
    }

    @Override
    public void setLambda(double lambda) {
        if (lambda >= 0 && lambda <= 1) {
            this.lambda = lambda;
            super.setLambda(lambda);

            Atom[] atoms = molecularAssembly.getAtomArray();
            int nAtoms = atoms.length;
            List<Atom> lambdaList = new ArrayList<>();

            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                if (atom.applyLambda()) {
                    lambdaList.add(atom);
                }
            }

            Atom[] atomArray = new Atom[lambdaList.size()];
            for (int i = 0; i < lambdaList.size(); i++) {
                atomArray[i] = lambdaList.get(i);
            }

            if (amoebaMultipoleForce != null) {
                scaleAmoebaMultipoleForceByLambda(atomArray, lambda);
            }
        } else {
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }

    }

    private void scaleAmoebaMultipoleForceByLambda(Atom[] atoms, double lambda) {
        ParticleMeshEwald pme = super.getPmeNode();
        int axisAtom[][] = pme.getAxisAtoms();
        double dipoleConversion = OpenMM_NmPerAngstrom;
        double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
        double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                * OpenMM_NmPerAngstrom;
        double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

        double polarScale = 1.0;
        if (pme.getPolarizationType() == Polarization.NONE) {
            polarScale = 0.0;
        }

        PointerByReference dipoles = OpenMM_DoubleArray_create(3);
        PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            MultipoleType multipoleType = atom.getMultipoleType();
            PolarizeType polarType = atom.getPolarizeType();

            double lambdaFactor = lambda;
            if (!atom.applyLambda()) {
                lambdaFactor = 1.0;
            }

            /**
             * Define the frame definition.
             */
            int axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
            switch (multipoleType.frameDefinition) {
                case ZONLY:
                    axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
                    break;
                case ZTHENX:
                    axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
                    break;
                case BISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_Bisector;
                    break;
                case ZTHENBISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
                    break;
                case TRISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                    break;
                default:
                    break;
            }

            /**
             * Load local multipole coefficients.
             */
            for (int j = 0; j < 3; j++) {
                OpenMM_DoubleArray_set(dipoles, j, multipoleType.dipole[j] * dipoleConversion * lambdaFactor);

            }
            int l = 0;
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    OpenMM_DoubleArray_set(quadrupoles, l++, multipoleType.quadrupole[j][k]
                            * quadrupoleConversion / 3.0 * lambdaFactor);
                }
            }

            int zaxis = 0;
            int xaxis = 0;
            int yaxis = 0;
            int refAtoms[] = axisAtom[i];
            if (refAtoms != null) {
                zaxis = refAtoms[0];
                if (refAtoms.length > 1) {
                    xaxis = refAtoms[1];
                    if (refAtoms.length > 2) {
                        yaxis = refAtoms[2];
                    }
                }
            }

            /**
             * Add the multipole.
             */
            OpenMM_AmoebaMultipoleForce_setMultipoleParameters(amoebaMultipoleForce, i,
                    multipoleType.charge * lambdaFactor, dipoles, quadrupoles,
                    axisType, zaxis, xaxis, yaxis,
                    polarType.thole,
                    polarType.pdamp * dampingFactorConversion,
                    polarType.polarizability * polarityConversion * polarScale * lambdaFactor);
        }
        OpenMM_DoubleArray_destroy(dipoles);
        OpenMM_DoubleArray_destroy(quadrupoles);

        OpenMM_AmoebaMultipoleForce_updateParametersInContext(amoebaMultipoleForce, openMMContext);
    }

    /**
     * Updates the AMOEBA Generalized Kirkwood force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     * @param useChanged Indices of atoms with changed Use flags
     */
    private void updateAmoebaGeneralizedKirkwoodForce(Atom[] atoms, int[] useChanged) {
        GeneralizedKirkwood gk = super.getGK();
        double overlapScale[] = gk.getOverlapScale();
        double baseRadii[] = gk.getBaseRadii();
        /*Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {*/
        int nAtoms = useChanged.length;
        for (int j = 0; j < nAtoms; j++) {
            int i = useChanged[j];
            double useFactor = 1.0;
            if (!use[i]) {
                useFactor = 0.0;
            }
            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(amoebaGeneralizedKirkwoodForce, i,
                    multipoleType.charge * useFactor, OpenMM_NmPerAngstrom * baseRadii[i], overlapScale[i] * useFactor);
        }
        OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(amoebaGeneralizedKirkwoodForce, openMMContext);
    }

    /**
     * Updates the WCA force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     * @param useChanged Indices of atoms with changed Use flags
     */
    private void updateWCAForce(Atom[] atoms, int[] useChanged) {
        VanDerWaals vdW = super.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }
        /*Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {*/
        int nAtoms = useChanged.length;
        for (int j = 0; j < nAtoms; j++) {
            int i = useChanged[j];
            double useFactor = 1.0;
            if (!use[i]) {
                useFactor = 0.0;
            }

            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            double radius = vdwType.radius;
            double eps = vdwType.wellDepth;
            OpenMM_AmoebaWcaDispersionForce_setParticleParameters(amoebaWcaDispersionForce, i,
                    OpenMM_NmPerAngstrom * radius * radScale,
                    OpenMM_KJPerKcal * eps * useFactor);
        }

        OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(amoebaWcaDispersionForce, openMMContext);
    }

    /**
     * Evaluates energy both with OpenMM and reference potential, and returns
     * the difference FFX-OpenMM.
     *
     * @param x Coordinate array
     * @param verbose
     * @return Energy discrepancy
     */
    public double energyVsFFX(double[] x, boolean verbose) {
        double ffxE = super.energy(x, verbose);
        double thisE = energy(x, verbose);
        return ffxE - thisE;
    }

    /**
     * Evaluates energy and gradients both with OpenMM and reference potential,
     * and returns the difference FFX-OpenMM.
     *
     * @param x Coordinate array
     * @param gFFX Array for FFX gradients to be stored in
     * @param gOMM Array for OpenMM gradients to be stored in
     * @param verbose
     * @return Energy discrepancy
     */
    public double energyAndGradVsFFX(double[] x, double[] gFFX, double[] gOMM, boolean verbose) {
        double ffxE = super.energyAndGradient(x, gFFX, verbose);
        double thisE = energyAndGradient(x, gOMM, verbose);
        return ffxE - thisE;
    }

    /**
     * Returns the current energy. Preferred is to use the methods with explicit
     * coordinate/gradient arrays.
     *
     * @return Current energy.
     */
    @Override
    public double energy() {
        return energy(false, false);
    }

    @Override
    public double energy(double[] x) {
        return energy(x, false);
    }

    @Override
    public double energy(double[] x, boolean verbose) {

        if (lambdaBondedTerms) {
            return 0.0;
        }

        updateParameters();

        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }
        setCoordinates(x);
        loadFFXPositionToOpenMM();

        int infoMask = OpenMM_State_Energy;
        openMMState = OpenMM_Context_getState(openMMContext, infoMask, 0);
        double e = OpenMM_State_getPotentialEnergy(openMMState) / OpenMM_KJPerKcal;

        if (verbose) {
            logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", e));
        }

        /**
         * Rescale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
            }
        }

        OpenMM_State_destroy(openMMState);
        return e;
    }

    @Override
    public double energyAndGradient(double x[], double g[]) {
        return energyAndGradient(x, g, false);
    }

    @Override
    public double energyAndGradient(double x[], double g[], boolean verbose) {
        if (lambdaBondedTerms) {
            return 0.0;
        }

        /**
         * Un-scale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }
        setCoordinates(x);
        loadFFXPositionToOpenMM();

        int infoMask = OpenMM_State_Energy;
        infoMask += OpenMM_State_Forces;

        openMMState = OpenMM_Context_getState(openMMContext, infoMask, 0);
        double e = OpenMM_State_getPotentialEnergy(openMMState) / OpenMM_KJPerKcal;

        if (verbose) {
            logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", e));
        }

        openMMForces = OpenMM_State_getForces(openMMState);

        getGradients(g);
        /**
         * Scale the coordinates and gradients.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }

        OpenMM_State_destroy(openMMState);
        return e;
    }

    @Override
    public void setCrystal(Crystal crystal) {
        super.setCrystal(crystal);
        setDefaultPeriodicBoxVectors();
        loadFFXPositionToOpenMM();
    }

    /**
     * <p>
     * getGradients</p>
     *
     * @param g an array of double.
     * @return
     */
    @Override
    public double[] getGradients(double g[]) {
        assert (g != null);
        int n = getNumberOfVariables();
        if (g.length < n) {
            g = new double[n];
        }
        int index = 0;
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                OpenMM_Vec3 posInNm = OpenMM_Vec3Array_get(openMMForces, i);
                /**
                 * Convert OpenMM Forces in KJ/Nm into an FFX gradient in
                 * Kcal/A.
                 */
                double gx = -posInNm.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double gy = -posInNm.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double gz = -posInNm.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                if (Double.isNaN(gx) || Double.isInfinite(gx)
                        || Double.isNaN(gy) || Double.isInfinite(gy)
                        || Double.isNaN(gz) || Double.isInfinite(gz)) {
                    /*String message = format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                            a.toString(), gx, gy, gz);*/
                    StringBuilder sb = new StringBuilder(format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                            a.toString(), gx, gy, gz));
                    double[] vals = new double[3];
                    a.getVelocity(vals);
                    sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getAcceleration(vals);
                    sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getPreviousAcceleration(vals);
                    sb.append(format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));

                    //logger.severe(message);
                    throw new EnergyException(sb.toString());
                }
                a.setXYZGradient(gx, gy, gz);
                g[index++] = gx;
                g[index++] = gy;
                g[index++] = gz;
            }
        }
        return g;
    }

    /**
     * Loads positions into OpenMM from the FFX data structure.
     */
    public final void loadFFXPositionToOpenMM() {
        if (openMMPositions == null) {
            openMMPositions = OpenMM_Vec3Array_create(0);
        } else {
            OpenMM_Vec3Array_resize(openMMPositions, 0);
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        OpenMM_Vec3.ByValue posInNm = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            posInNm.x = atom.getX() * OpenMM_NmPerAngstrom;
            posInNm.y = atom.getY() * OpenMM_NmPerAngstrom;
            posInNm.z = atom.getZ() * OpenMM_NmPerAngstrom;
            OpenMM_Vec3Array_append(openMMPositions, posInNm);
        }

        // Load positions into the openMMContext.
        OpenMM_Context_setPositions(openMMContext, openMMPositions);
    }

    /**
     * setOpenMMPositions takes in an array of doubles generated by the DYN
     * reader method and appends these values to a Vec3Array. Finally this
     * method sets the created Vec3Array as the positions of the openMMContext.
     *
     * @param x
     * @param numberOfVariables
     */
    public void setOpenMMPositions(double x[], int numberOfVariables) {
        if (openMMPositions == null) {
            openMMPositions = OpenMM_Vec3Array_create(0);
        } else {
            OpenMM_Vec3Array_resize(openMMPositions, 0);
        }
        OpenMM_Vec3.ByValue pos = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < numberOfVariables; i = i + 3) {
            pos.x = x[i] * OpenMM_NmPerAngstrom;
            pos.y = x[i + 1] * OpenMM_NmPerAngstrom;
            pos.z = x[i + 2] * OpenMM_NmPerAngstrom;
            OpenMM_Vec3Array_append(openMMPositions, pos);
        }

        OpenMM_Context_setPositions(openMMContext, openMMPositions);
    }

    /**
     * setOpenMMVelocities takes in an array of doubles generated by the DYN
     * reader method and appends these values to a Vec3Array. Finally this
     * method sets the created Vec3Arrat as the velocities of the openMMContext.
     *
     * @param v
     * @param numberOfVariables
     */
    public void setOpenMMVelocities(double v[], int numberOfVariables) {
        if (openMMVelocities == null) {
            openMMVelocities = OpenMM_Vec3Array_create(0);
        } else {
            OpenMM_Vec3Array_resize(openMMVelocities, 0);
        }
        OpenMM_Vec3.ByValue vel = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < numberOfVariables; i = i + 3) {
            vel.x = v[i] * OpenMM_NmPerAngstrom;
            vel.y = v[i + 1] * OpenMM_NmPerAngstrom;
            vel.z = v[i + 2] * OpenMM_NmPerAngstrom;
            OpenMM_Vec3Array_append(openMMVelocities, vel);
        }
        OpenMM_Context_setVelocities(openMMContext, openMMVelocities);
    }

    /**
     * getOpenMMPositions takes in a PointerByReference containing the position
     * information of the openMMContext. This method creates a Vec3Array that
     * contains the three dimensional information of the positions of the atoms.
     * The method then adds these values to a new double array x and returns it
     * to the method call
     *
     * @param positions
     * @param numberOfVariables
     * @param x
     * @return x
     */
    public double[] getOpenMMPositions(PointerByReference positions, int numberOfVariables, double x[]) {
        if (x == null || x.length < numberOfVariables) {
            x = new double[numberOfVariables];
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            int offset = i * 3;
            OpenMM_Vec3 pos = OpenMM_Vec3Array_get(positions, i);
            x[offset] = pos.x * OpenMM_AngstromsPerNm;
            x[offset + 1] = pos.y * OpenMM_AngstromsPerNm;
            x[offset + 2] = pos.z * OpenMM_AngstromsPerNm;
            Atom atom = atoms[i];
            atom.moveTo(x[offset], x[offset + 1], x[offset + 2]);
        }
        return x;
    }

    /**
     * getOpenMMVelocities takes in a PointerByReference containing the velocity
     * information of the openMMContext. This method creates a Vec3Array that
     * contains the three dimensional information of the velocities of the
     * atoms. This method then adds these values to a new double array v and
     * returns it to the method call
     *
     * @param velocities
     * @param numberOfVariables
     * @return
     */
    public double[] getOpenMMVelocities(PointerByReference velocities, int numberOfVariables, double v[]) {
        if (v == null || v.length < numberOfVariables) {
            v = new double[numberOfVariables];
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            int offset = i * 3;
            OpenMM_Vec3 vel = OpenMM_Vec3Array_get(velocities, i);
            v[offset] = vel.x * OpenMM_AngstromsPerNm;
            v[offset + 1] = vel.y * OpenMM_AngstromsPerNm;
            v[offset + 2] = vel.z * OpenMM_AngstromsPerNm;
        }
        return v;
    }

    /**
     * getOpenMMAccelerations takes in a PointerByReference containing the force
     * information of the openMMContext. This method creates a Vec3Array that
     * contains the three dimensional information of the forces on the atoms.
     * This method then adds these values (divided by mass, effectively turning
     * them into accelerations) to a new double array a and returns it to the
     * method call
     *
     * @param accelerations
     * @param numberOfVariables
     * @param mass
     * @return
     */
    public double[] getOpenMMAccelerations(PointerByReference accelerations, int numberOfVariables,
            double[] mass, double[] a) {
        if (a == null || a.length < numberOfVariables) {
            a = new double[numberOfVariables];
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            int offset = i * 3;
            OpenMM_Vec3 acc = OpenMM_Vec3Array_get(accelerations, i);
            a[offset] = (acc.x * 10.0) / mass[i];
            a[offset + 1] = (acc.y * 10.0) / mass[i + 1];
            a[offset + 2] = (acc.z * 10.0) / mass[i + 2];
        }
        return a;
    }

    private void setDefaultPeriodicBoxVectors() {

        OpenMM_Vec3 a = new OpenMM_Vec3();
        OpenMM_Vec3 b = new OpenMM_Vec3();
        OpenMM_Vec3 c = new OpenMM_Vec3();

        Crystal crystal = super.getCrystal();

        if (!crystal.aperiodic()) {
            a.x = crystal.a * OpenMM_NmPerAngstrom;
            a.y = 0.0 * OpenMM_NmPerAngstrom;
            a.z = 0.0 * OpenMM_NmPerAngstrom;
            b.x = 0.0 * OpenMM_NmPerAngstrom;
            b.y = crystal.b * OpenMM_NmPerAngstrom;
            b.z = 0.0 * OpenMM_NmPerAngstrom;
            c.x = 0.0 * OpenMM_NmPerAngstrom;
            c.y = 0.0 * OpenMM_NmPerAngstrom;
            c.z = crystal.c * OpenMM_NmPerAngstrom;
            OpenMM_System_setDefaultPeriodicBoxVectors(openMMSystem, a, b, c);
        }
    }

    /**
     * getIntegrator returns the integrator used for the openMMContext
     *
     * @return openMMIntegrator
     */
    public PointerByReference getIntegrator() {
        return openMMIntegrator;
    }

    /**
     * setIntegrator takes in a parameters to determine which integrator the
     * user requested during the start up of the simulation. A switch statement
     * is used with Strings as the variable to determine between Lengevin,
     * Brownian, Custom, Compound and Verlet integrator
     *
     * @param integrator
     * @param timeStep
     * @param frictionCoeff
     * @param temperature
     * @param collisionFreq
     */
    public void setIntegrator(String integrator, double timeStep, double frictionCoeff, double temperature, double collisionFreq) {
        OpenMM_Context_destroy(openMMContext);
        double dt = timeStep * 1.0e-3;
        switch (integrator) {
            case "LANGEVIN":
                openMMIntegrator = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, dt);
                break;
            /*
            case "BROWNIAN":
                openMMIntegrator = OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, dt);
                break;
            case "CUSTOM":
                openMMIntegrator = OpenMM_CustomIntegrator_create(dt);
                break;
            case "COMPOUND":
                openMMIntegrator = OpenMM_CompoundIntegrator_create();
                break;
             */
            case "VERLET":
            default:
                openMMIntegrator = OpenMM_VerletIntegrator_create(dt);
            //thermostat = OpenMM_AndersenThermostat_create(temperature, collisionFreq);
            //OpenMM_System_addForce(openMMSystem, thermostat);
        }
        //logger.info(String.format(" Created %s OpenMM Integrator", integrator));

        // Create a openMMContext.
        openMMContext = OpenMM_Context_create_2(openMMSystem, openMMIntegrator, openMMPlatform);

        // Set initial positions.
        loadFFXPositionToOpenMM();

        int infoMask = OpenMM_State_Positions;
        infoMask += OpenMM_State_Forces;
        infoMask += OpenMM_State_Energy;

        openMMState = OpenMM_Context_getState(openMMContext, infoMask, 0);
        openMMForces = OpenMM_State_getForces(openMMState);
        double openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(openMMState) / OpenMM_KJPerKcal;

        //logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", openMMPotentialEnergy));
    }

    /**
     * Returns the openMMContext created for the OpenMMForceFieldEnergy object
     *
     * @return openMMContext
     */
    public PointerByReference getContext() {
        return openMMContext;
    }

    /**
     * Sets the finite-difference step size used for getdEdL.
     *
     * @param fdDLambda FD step size.
     */
    public void setFdDLambda(double fdDLambda) {
        this.fdDLambda = fdDLambda;
    }

    private boolean doOpenMMdEdL = false;
    private boolean doFFXdEdL = true;

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        if (!lambdaTerm) {
            return 0.0;
        }

        double currentLambda = lambda;
        double width = fdDLambda;
        double ePlus;
        double eMinus;
        double dEdL = 0.0;
        double openMMdEdL = 0.0;
        double ffxdEdL = 0.0;

        // Small optimization to only create the x array once.
        double[] x = new double[getNumberOfVariables()];
        getCoordinates(x);

        if (doOpenMMdEdL) {
            width = fdDLambda;
            if (currentLambda + fdDLambda > 1.0) {
                logger.fine(" Could not test the upper point, as current lambda + fdDL > 1");
                ePlus = energy(x);
                setLambda(currentLambda - fdDLambda);
                eMinus = energy(x);
            } else if (currentLambda - fdDLambda < 0.0) {
                logger.fine(" Could not test the lower point, as current lambda - fdDL < 1");
                eMinus = energy(x);
                setLambda(currentLambda + fdDLambda);
                ePlus = energy(x);
            } else {
                setLambda(currentLambda + fdDLambda);
                ePlus = energy(x);
                setLambda(currentLambda - fdDLambda);
                eMinus = energy(x);
                //logger.info(String.format(" OpenMM %12.8f %12.8f", ePlus, eMinus));
                width *= 2.0;
            }
            // Reset Lambda.
            setLambda(currentLambda);
            openMMdEdL = (ePlus - eMinus) / width;
            //logger.info(String.format(" Step: %16.10f OpenMM: %16.10f", fdDLambda, openMMdEdL));
            dEdL = openMMdEdL;
        }

        if (doFFXdEdL || !doOpenMMdEdL) {
            width = fdDLambda;
            // This section technically not robust to the case that fdDLambda > 0.5.
            // However, that should be an error case checked when fdDLambda is set.
            super.setLambda(1.0);
            if (currentLambda + fdDLambda > 1.0) {
                logger.fine(" Could not test the upper point, as current lambda + fdDL > 1");
                super.setLambdaMultipoleScale(currentLambda);
                ePlus = super.energy(x, false);
                // ePlus = super.getTotalElectrostaticEnergy();
                super.setLambdaMultipoleScale(currentLambda - fdDLambda);
                eMinus = super.energy(x, false);
                // eMinus = super.getTotalElectrostaticEnergy();
            } else if (currentLambda - fdDLambda < 0.0) {
                logger.fine(" Could not test the lower point, as current lambda - fdDL < 1");
                super.setLambdaMultipoleScale(currentLambda);
                eMinus = super.energy(x, false);
                // eMinus = super.getTotalElectrostaticEnergy();
                super.setLambdaMultipoleScale(currentLambda + fdDLambda);
                ePlus = super.energy(x, false);
                // ePlus = super.getTotalElectrostaticEnergy();
            } else {
                super.setLambdaMultipoleScale(currentLambda + fdDLambda);
                ePlus = super.energy(x, false);
                // ePlus = super.getTotalElectrostaticEnergy();
                super.setLambdaMultipoleScale(currentLambda - fdDLambda);
                eMinus = super.energy(x, false);
                // eMinus = super.getTotalElectrostaticEnergy();
                //logger.info(String.format(" FFX    %12.8f %12.8f", ePlus, eMinus));
                width *= 2.0;
            }
            super.setLambdaMultipoleScale(currentLambda);
            ffxdEdL = (ePlus - eMinus) / width;
            //logger.info(String.format(" Step: %16.10f FFX:    %16.10f", fdDLambda, ffxdEdL));
            dEdL = ffxdEdL;
        }

        return dEdL;
    }

    /**
     * {@inheritDoc}
     *
     * @param gradients
     */
    @Override
    public void getdEdXdL(double gradients[]) {
        // Note for OpenMMForceFieldEnergy this method is not implemented.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        // Note for OpenMMForceFieldEnergy this method is not implemented.
        return 0.0;
    }
}
