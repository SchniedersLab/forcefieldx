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

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import static org.apache.commons.math3.util.FastMath.sqrt;

import simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod;
import simtk.openmm.OpenMMLibrary.OpenMM_Boolean;
import simtk.openmm.OpenMM_Vec3;

import static simtk.openmm.OpenMMAmoebaLibrary.*;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent12;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent13;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent14;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent15;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent11;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent12;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent13;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent14;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_Bisector;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_NoAxisType;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ThreeFold;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZBisect;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZOnly;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZThenX;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_NoCutoff;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_PME;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Direct;
import static simtk.openmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Mutual;
import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static simtk.openmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.*;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PiTorsionType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWType;

import static ffx.potential.parameters.ForceField.toPropertyForm;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;

/**
 * Compute the potential energy and derivatives using OpenMM.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class OpenMMForceFieldEnergy extends ForceFieldEnergy {

    private static final Logger logger = Logger.getLogger(OpenMMForceFieldEnergy.class.getName());

    private ForceFieldEnergy ffxForceFieldEnergy;

    private PointerByReference openMMSystem;
    private PointerByReference openMMIntegrator;
    private PointerByReference platform;
    private PointerByReference context;
    private PointerByReference state;

    private PointerByReference initialPosInNm;
    private PointerByReference openMMForces;

    /**
     * OpenMMForceFieldEnergy constructor.
     *
     * @param molecularAssembly
     */
    public OpenMMForceFieldEnergy(MolecularAssembly molecularAssembly) {
        super(molecularAssembly);

        ffxForceFieldEnergy = molecularAssembly.getPotentialEnergy();

        logger.info("\n\n Initializing OpenMM");

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
        platforms = OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));

        platforms = OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));
        int numPlatforms = OpenMM_Platform_getNumPlatforms();
        boolean cuda = false;
        logger.log(Level.INFO, " Number of OpenMM Plugins: {0}", numPlatforms);
        for (int i = 0; i < numPlatforms; i++) {
            Pointer platformPtr = OpenMM_StringArray_get(platforms, i);
            String platform = platformPtr.getString(0);
            logger.log(Level.INFO, " Plugin Library :{0}", platform);
            if (platform.toUpperCase().contains("AMOEBACUDA")) {
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

        // Create the OpenMM System
        openMMSystem = OpenMM_System_create();
        logger.info(" Created OpenMM System");

        openMMIntegrator = OpenMM_VerletIntegrator_create(0.001);
        logger.info(" Created OpenMM Integrator");

        if (cuda) {
            platform = OpenMM_Platform_getPlatformByName("CUDA");
            logger.info(" Created OpenMM AMOEBA CUDA Plaform");
        } else {
            platform = OpenMM_Platform_getPlatformByName("Reference");
            logger.info(" Created OpenMM AMOEBA Reference Plaform");
        }

        // Load atoms.
        addAtoms();

        // CCOM remover.
        addCCOMRemover();

        // Add Bond Force.
        addBonds();

        // Add Angle Force.
        addAngles();
        addInPlaneAngles();

        // Add Urey-Bradley Force.
        addUreyBradleys();

        // Add Stretch-Bend Force.
        addStretchBendForce();

        // Add Torsion Force.
        addTorsions();

        // Add Pi-Torsion Force.
        addPiTorsions();

        // Add Torsion-Torsion Force.
        addTorsionTorsions();

        // Add vdW Force.
        addVDWForce();

        // Add Multipole Force.
        addMultipoleForce();

        // Set periodic box vectors.
        setDefaultPeriodicBoxVectors();

        // Set initial position.
        loadPositions();

        // Create a context.
        context = OpenMM_Context_create_2(openMMSystem, openMMIntegrator, platform);

        // Load positions into the context.
        OpenMM_Context_setPositions(context, initialPosInNm);

        int infoMask = OpenMM_State_Positions;
        infoMask += OpenMM_State_Forces;
        infoMask += OpenMM_State_Energy;

        state = OpenMM_Context_getState(context, infoMask, 0);

        openMMForces = OpenMM_State_getForces(state);
        double openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;

        logger.log(Level.INFO, " OpenMM Energy: {0}", openMMPotentialEnergy);

        // Free the OpenMM System.
        freeOpenMM();
        logger.info(" Destroyed the Context, Integrator, and OpenMMSystem.");
    }

    private void freeOpenMM() {
        OpenMM_Context_destroy(context);
        OpenMM_Integrator_destroy(openMMIntegrator);
        OpenMM_System_destroy(openMMSystem);
    }

    private void addAtoms() {
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            OpenMM_System_addParticle(openMMSystem, atom.getMass());
        }
        logger.log(Level.INFO, " Added particles ({0})", nAtoms);
    }

    private void addCCOMRemover() {
        int frequency = 100;
        PointerByReference cMMotionRemover = OpenMM_CMMotionRemover_create(frequency);
        OpenMM_System_addForce(openMMSystem, cMMotionRemover);
        logger.log(Level.INFO, " Added center of mass motion remover (frequency: {0})", frequency);
    }

    private void addBonds() {
        Bond bonds[] = ffxForceFieldEnergy.getBonds();
        int nBonds = bonds.length;

        if (nBonds < 1) {
            return;
        }

        PointerByReference amoebaBondForce = OpenMM_AmoebaBondForce_create();
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
        OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaBondForce,
                BondType.cubic / OpenMM_NmPerAngstrom);
        OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaBondForce,
                BondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));

        OpenMM_System_addForce(openMMSystem, amoebaBondForce);
        logger.log(Level.INFO, " Added bonds ({0})", nBonds);
    }

    private void addAngles() {
        Angle angles[] = ffxForceFieldEnergy.getAngles();
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
        PointerByReference amoebaAngleForce = OpenMM_AmoebaAngleForce_create();
        for (int i = 0; i < nAngles; i++) {
            Angle angle = normalAngles.get(i);
            int i1 = angle.getAtom(0).getXyzIndex() - 1;
            int i2 = angle.getAtom(1).getXyzIndex() - 1;
            int i3 = angle.getAtom(2).getXyzIndex() - 1;
            int nh = angle.nh;
            OpenMM_AmoebaAngleForce_addAngle(amoebaAngleForce, i1, i2, i3,
                    angle.angleType.angle[nh], OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant);
        }
        OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(amoebaAngleForce, AngleType.cubic);
        OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(amoebaAngleForce, AngleType.quartic);
        OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(amoebaAngleForce, AngleType.quintic);
        OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(amoebaAngleForce, AngleType.sextic);
        OpenMM_System_addForce(openMMSystem, amoebaAngleForce);
        logger.log(Level.INFO, " Added angles ({0})", nAngles);
    }

    private void addInPlaneAngles() {
        Angle angles[] = ffxForceFieldEnergy.getAngles();
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
        PointerByReference amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create();
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

    private void addUreyBradleys() {
        UreyBradley ureyBradleys[] = ffxForceFieldEnergy.getUreyBradleys();
        if (ureyBradleys == null || ureyBradleys.length < 1) {
            return;
        }
        PointerByReference amoebaBondForce = OpenMM_AmoebaBondForce_create();
        double kParameterConversion = UreyBradleyType.units * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
        int nUreys = ureyBradleys.length;
        for (int i = 0; i < nUreys; i++) {
            UreyBradley ureyBradley = ureyBradleys[i];
            int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
            int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
            UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
            OpenMM_AmoebaBondForce_addBond(amoebaBondForce, i1, i2,
                    ureyBradleyType.distance * OpenMM_NmPerAngstrom,
                    ureyBradleyType.forceConstant * kParameterConversion);
        }

        OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaBondForce,
                UreyBradleyType.cubic / OpenMM_NmPerAngstrom);
        OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaBondForce,
                UreyBradleyType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));

        OpenMM_System_addForce(openMMSystem, amoebaBondForce);
        logger.log(Level.INFO, " Added Urey-Bradleys ({0})", nUreys);
    }

    private void addOutOfPlaneBendForce(){
        PointerByReference amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create();
        OutOfPlaneBend outOfPlaneBends[] = ffxForceFieldEnergy.getOutOfPlaneBends();
        int nOutOfPlaneBends = outOfPlaneBends.length;

        for (int i = 0; i < nOutOfPlaneBends; i++){
            OutOfPlaneBend outOfPlaneBend = outOfPlaneBends[i];
            int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
            int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
            int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
            int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
            OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;

            OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(amoebaOutOfPlaneBendForce, i1, i2, i3, i4,
                    OpenMM_KJPerKcal*outOfPlaneBend.outOfPlaneBendType.forceConstant*OutOfPlaneBendType.units);
        }
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.cubic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.quartic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.quintic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.sextic);

        OpenMM_System_addForce(openMMSystem, amoebaOutOfPlaneBendForce);
        logger.log(Level.INFO, " Added Out Of Plane Bends ({0})", nOutOfPlaneBends);
    }

    private void addStretchBendForce(){
        StretchBend stretchBends[] = ffxForceFieldEnergy.getStretchBends();
        if (stretchBends == null || stretchBends.length < 1) {
            return;
        }
        int nStretchBends = stretchBends.length;
        PointerByReference amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create();
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

    private void addTorsions() {
        Torsion torsions[] = ffxForceFieldEnergy.getTorsions();
        if (torsions == null || torsions.length < 1) {
            return;
        }
        int nTorsions = torsions.length;
        PointerByReference amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create();
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

    private void addPiTorsions() {
        PiOrbitalTorsion piOrbitalTorsions[] = ffxForceFieldEnergy.getPiOrbitalTorsions();
        if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
            return;
        }
        int nPiOrbitalTorsions = piOrbitalTorsions.length;
        PointerByReference amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create();
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

    private void addTorsionTorsions() {
        TorsionTorsion torsionTorsions[] = ffxForceFieldEnergy.getTorsionTorsions();
        if (torsionTorsions == null || torsionTorsions.length < 1) {
            return;
        }
        /**
         * Load the torsion-torsions.
         */

        int nTypes = 0;
        LinkedHashMap<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();

        int nTorsionTorsions = torsionTorsions.length;
        PointerByReference amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create();
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

    private void addVDWForce() {
        VanDerWaals vdW = ffxForceFieldEnergy.getVdwNode();
        if (vdW == null) {
            return;
        }

        PointerByReference amoebaVdwForce = OpenMM_AmoebaVdwForce_create();
        OpenMM_System_addForce(openMMSystem, amoebaVdwForce);
        OpenMM_Force_setForceGroup(amoebaVdwForce, 1);

        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
        Crystal crystal = ffxForceFieldEnergy.getCrystal();

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
            OpenMM_AmoebaVdwForce_addParticle(amoebaVdwForce,
                    ired[i], OpenMM_NmPerAngstrom * vdwType.radius * radScale,
                    OpenMM_KJPerKcal * vdwType.wellDepth,
                    vdwType.reductionFactor);
        }

        // OpenMM_AmoebaVdwForce_setSigmaCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.radiusRule.name()));
        // OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.epsilonRule.name()));
        OpenMM_AmoebaVdwForce_setCutoffDistance(amoebaVdwForce, nonbondedCutoff.off * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVdwForce, OpenMM_Boolean.OpenMM_False);

        if (crystal.aperiodic()) {
            OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVdwForce,
                    OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff);
        } else {
            OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVdwForce,
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
            OpenMM_AmoebaVdwForce_setParticleExclusions(amoebaVdwForce, i, exclusions);
            OpenMM_IntArray_resize(exclusions, 0);
        }
        OpenMM_IntArray_destroy(exclusions);
        logger.log(Level.INFO, " Added van der Waals force.");
    }

    private void addMultipoleForce() {
        ParticleMeshEwald pme = ffxForceFieldEnergy.getPmeNode();
        if (pme == null) {
            return;
        }
        int axisAtom[][] = pme.getAxisAtoms();
        double dipoleConversion = OpenMM_NmPerAngstrom;
        double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
        double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                * OpenMM_NmPerAngstrom;
        double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

        PointerByReference amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();
        OpenMM_System_addForce(openMMSystem, amoebaMultipoleForce);
        OpenMM_Force_setForceGroup(amoebaMultipoleForce, 1);

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
                    polarType.polarizability * polarityConversion);
        }
        OpenMM_DoubleArray_destroy(dipoles);
        OpenMM_DoubleArray_destroy(quadrupoles);

        Crystal crystal = ffxForceFieldEnergy.getCrystal();
        if (!crystal.aperiodic()) {
            double ewaldTolerance = 1.0e-04;
            OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_PME);
            OpenMM_AmoebaMultipoleForce_setCutoffDistance(amoebaMultipoleForce,
                    pme.getEwaldCutoff() * OpenMM_NmPerAngstrom);
            OpenMM_AmoebaMultipoleForce_setAEwald(amoebaMultipoleForce,
                    pme.getEwaldCoefficient() / OpenMM_NmPerAngstrom);

            PointerByReference gridDimensions = OpenMM_IntArray_create(3);
            ReciprocalSpace recip = pme.getReciprocalSpace();
            OpenMM_IntArray_set(gridDimensions, 0, recip.getXDim());
            OpenMM_IntArray_set(gridDimensions, 1, recip.getYDim());
            OpenMM_IntArray_set(gridDimensions, 2, recip.getZDim());
            OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(amoebaMultipoleForce,
                    gridDimensions);
            OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(amoebaMultipoleForce, ewaldTolerance);
            OpenMM_IntArray_destroy(gridDimensions);
        } else {
            OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_NoCutoff);
        }

        if (pme.getPolarizationType() == Polarization.DIRECT) {
            OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Direct);
        } else {
            OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Mutual);
        }

        OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(amoebaMultipoleForce, 500);
        OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(amoebaMultipoleForce, pme.getPolarEps());

        int ip11[][] = pme.getPolarization11();
        int ip12[][] = pme.getPolarization12();
        int ip13[][] = pme.getPolarization13();

        PointerByReference covalentMap = OpenMM_IntArray_create(0);
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            for (Bond bond : ai.getBonds()) {
                int index = bond.get1_2(ai).getIndex() - 1;
                OpenMM_IntArray_append(covalentMap, index);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Angle angle : ai.getAngles()) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    int index = ak.getIndex() - 1;
                    OpenMM_IntArray_append(covalentMap, index);
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Torsion torsion : ai.getTorsions()) {
                Atom ak = torsion.get1_4(ai);
                if (ak != null) {
                    int index = ak.getIndex() - 1;
                    OpenMM_IntArray_append(covalentMap, index);
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Atom ak : ai.get1_5s()) {
                int index = ak.getIndex() - 1;
                OpenMM_IntArray_append(covalentMap, index);
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

            for (int j = 0; j < ip12[i].length; j++) {
                OpenMM_IntArray_append(covalentMap, ip12[i][j]);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent12, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (int j = 0; j < ip13[i].length; j++) {
                OpenMM_IntArray_append(covalentMap, ip13[i][j]);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent13, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent14, covalentMap);
        }

        OpenMM_IntArray_destroy(covalentMap);
        logger.log(Level.INFO, " Added polarizable multipole force.");

        GeneralizedKirkwood gk = ffxForceFieldEnergy.getGK();
        if (gk != null) {
            addGKForce();
        }

    }

    private void addGKForce() {

        GeneralizedKirkwood gk = ffxForceFieldEnergy.getGK();

        PointerByReference amoebaGeneralizedKirkwoodForce = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(amoebaGeneralizedKirkwoodForce, 78.3);
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(amoebaGeneralizedKirkwoodForce, 1.0);

        double overlapScale[] = gk.getOverlapScale();
        double baseRadii[] = gk.getBaseRadii();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(amoebaGeneralizedKirkwoodForce,
                    multipoleType.charge,
                    OpenMM_NmPerAngstrom * baseRadii[i], overlapScale[i]);
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

        VanDerWaals vdW = ffxForceFieldEnergy.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        PointerByReference amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create();

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

    private void loadPositions() {
        initialPosInNm = OpenMM_Vec3Array_create(0);
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        OpenMM_Vec3.ByValue posInNm = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            posInNm.x = atom.getX() * OpenMM_NmPerAngstrom;
            posInNm.y = atom.getY() * OpenMM_NmPerAngstrom;
            posInNm.z = atom.getZ() * OpenMM_NmPerAngstrom;
            OpenMM_Vec3Array_append(initialPosInNm, posInNm);
        }
    }

    private void setDefaultPeriodicBoxVectors() {

        OpenMM_Vec3 a = new OpenMM_Vec3();
        OpenMM_Vec3 b = new OpenMM_Vec3();
        OpenMM_Vec3 c = new OpenMM_Vec3();

        Crystal crystal = ffxForceFieldEnergy.getCrystal();

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
}
