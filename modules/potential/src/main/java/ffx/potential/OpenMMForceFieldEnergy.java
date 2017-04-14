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

import java.util.logging.Level;
import java.util.logging.Logger;

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import simtk.openmm.OpenMMAmoebaLibrary;
import simtk.openmm.OpenMMLibrary;
import simtk.openmm.OpenMM_Vec3;

import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.*;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.BondType;

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

    private final OpenMMLibrary openmm = OpenMMLibrary.INSTANCE;
    private final OpenMMAmoebaLibrary amoeba = OpenMMAmoebaLibrary.INSTANCE;

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
        Pointer version = openmm.OpenMM_Platform_getOpenMMVersion();
        logger.log(Level.INFO, " OpenMM Version: {0}", version.getString(0));

        // Print out the OpenMM plugin directory.
        Pointer pluginDir = openmm.OpenMM_Platform_getDefaultPluginsDirectory();
        logger.log(Level.INFO, " OpenMM Plugin Dir: {0}", pluginDir.getString(0));

        // Load plugins and print out plugins.
        PointerByReference platforms = openmm.OpenMM_Platform_loadPluginsFromDirectory(pluginDir.getString(0));
        int numPlatforms = openmm.OpenMM_Platform_getNumPlatforms();
        logger.log(Level.INFO, " Number of OpenMM Plugins: {0}", numPlatforms);
        for (int i = 0; i < numPlatforms; i++) {
            Pointer platformPtr = openmm.OpenMM_StringArray_get(platforms, i);
            logger.log(Level.INFO, " Plugin Library :{0}", platformPtr.getString(0));
        }
        openmm.OpenMM_StringArray_destroy(platforms);

        // Create the OpenMM System
        openMMSystem = openmm.OpenMM_System_create();
        logger.info(" Created OpenMM System");

        openMMIntegrator = openmm.OpenMM_LangevinIntegrator_create(300.0, 30.0, 0.001);
        logger.info(" Created OpenMM Integrator");

        platform = openmm.OpenMM_Platform_getPlatformByName("Reference");
        logger.info(" Created OpenMM Reference Plaform");

        // Load atoms.
        addAtoms();

        // CCOM remover.
        addCCOMRemover();

        // Add Bond Forces.
        addBonds();

        // Reference: https://github.com/jayponder/tinker/blob/release/openmm/ommstuff.cpp

        // Add Angle Forces: to do by Mallory - see setupAmoebaAngleForce line 1952 of ommsetuff.cpp

        // Add Urey-Bradley Forces: to do by Hernan - see setupAmoebaUreyBradleyForce line 2115 of openmm-stuff.cpp

        // Set initial position.
        loadPositions();

        // Create a context.
        context = openmm.OpenMM_Context_create_2(openMMSystem, openMMIntegrator, platform);

        // Load positions into the context.
        openmm.OpenMM_Context_setPositions(context, initialPosInNm);

        int infoMask = OpenMM_State_Positions;
        infoMask += OpenMM_State_Forces;
        infoMask += OpenMM_State_Energy;

        state = openmm.OpenMM_Context_getState(context, infoMask, 0);

        openMMForces = openmm.OpenMM_State_getForces(state);
        double openMMPotentialEnergy = openmm.OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;

        logger.log(Level.INFO, " OpenMM Energy: {0}", openMMPotentialEnergy);

        // Free the OpenMM System.
        freeOpenMM();
        logger.info(" Destroyed the Context, Integrator, and OpenMMSystem.");
    }

    private void freeOpenMM() {
        openmm.OpenMM_Context_destroy(context);
        openmm.OpenMM_Integrator_destroy(openMMIntegrator);
        openmm.OpenMM_System_destroy(openMMSystem);
    }

    private void addAtoms() {
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            openmm.OpenMM_System_addParticle(openMMSystem, atom.getMass());
        }
        logger.log(Level.INFO, " Added particles ({0})", nAtoms);
    }

    private void addCCOMRemover() {
        int frequency = 100;
        PointerByReference cMMotionRemover = openmm.OpenMM_CMMotionRemover_create(frequency);
        openmm.OpenMM_System_addForce(openMMSystem, cMMotionRemover);
        logger.log(Level.INFO, " Added center of mass motion remover (frequency: {0})", frequency);
    }

    private void addBonds() {
        PointerByReference amoebaBondForce = amoeba.OpenMM_AmoebaBondForce_create();
        double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
        Bond bonds[] = ffxForceFieldEnergy.getBonds();
        int nBonds = bonds.length;
        for (int i = 0; i < nBonds; i++) {
            Bond bond = bonds[i];
            int i1 = bond.getAtom(0).getXyzIndex() - 1;
            int i2 = bond.getAtom(1).getXyzIndex() - 1;
            BondType bondType = bond.bondType;
            amoeba.OpenMM_AmoebaBondForce_addBond(amoebaBondForce, i1, i2,
                    bond.bondType.distance * OpenMM_NmPerAngstrom,
                    kParameterConversion * bondType.forceConstant * BondType.units);

        }
        amoeba.OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaBondForce,
                BondType.cubic / OpenMM_NmPerAngstrom);
        amoeba.OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaBondForce,
                BondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));

        openmm.OpenMM_System_addForce(openMMSystem, amoebaBondForce);
        logger.log(Level.INFO, " Added bonds ({0})", nBonds);
    }

    private void loadPositions() {
        initialPosInNm = openmm.OpenMM_Vec3Array_create(0);

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        OpenMM_Vec3.ByValue posInNm = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            posInNm.x = atom.getX() * OpenMM_NmPerAngstrom;
            posInNm.y = atom.getY() * OpenMM_NmPerAngstrom;
            posInNm.z = atom.getZ() * OpenMM_NmPerAngstrom;
            openmm.OpenMM_Vec3Array_append(initialPosInNm, posInNm);
        }
    }
}
