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

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import simtk.openmm.OpenMM_Vec3;

import static simtk.openmm.OpenMMAmoebaLibrary.*;
import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.*;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.VDWType;

import static ffx.potential.parameters.ForceField.toPropertyForm;

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
         * Call the method twice to avoid a bug in OpenMM where not all platforms are
         * list after the first call.
         */
        PointerByReference platforms = OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));
        OpenMM_StringArray_destroy(platforms);
        platforms = OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));

        platforms = OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));
        int numPlatforms = OpenMM_Platform_getNumPlatforms();
        logger.log(Level.INFO, " Number of OpenMM Plugins: {0}", numPlatforms);
        for (int i = 0; i < numPlatforms; i++) {
            Pointer platformPtr = OpenMM_StringArray_get(platforms, i);
            logger.log(Level.INFO, " Plugin Library :{0}", platformPtr.getString(0));
        }
        OpenMM_StringArray_destroy(platforms);

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

        platform = OpenMM_Platform_getPlatformByName("Reference");

        if (platform == null) {
            logger.info(" OpenMM Plaform could not be created.");
        } else {
            logger.info(" Created OpenMM Reference Plaform");
        }

        // Load atoms.
        addAtoms();

        // CCOM remover.
        addCCOMRemover();

        // Add Bond Forces.
        // addBonds();
        // Reference: https://github.com/jayponder/tinker/blob/release/openmm/ommstuff.cpp
        // Add Angle Forces: to do by Mallory - see setupAmoebaAngleForce line 1952 of ommsetuff.cpp
        // Add Urey-Bradley Forces: to do by Hernan - see setupAmoebaUreyBradleyForce line 2115 of openmm-stuff.cpp
        // Add vdW force.
        addVDWForce();

        // Add multipole forces.
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
        PointerByReference amoebaBondForce = OpenMM_AmoebaBondForce_create();
        double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
        Bond bonds[] = ffxForceFieldEnergy.getBonds();
        int nBonds = bonds.length;
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

    private void addVDWForce() {
        PointerByReference amoebaVdwForce = OpenMM_AmoebaVdwForce_create();
        OpenMM_System_addForce(openMMSystem, amoebaVdwForce);
        OpenMM_Force_setForceGroup(amoebaVdwForce, 1);

        VanDerWaals vdW = ffxForceFieldEnergy.getVdwNode();
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

        OpenMM_AmoebaVdwForce_setSigmaCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.radiusRule.name()));
        OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.epsilonRule.name()));
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
    }

    private void addMultipoleForce() {

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
