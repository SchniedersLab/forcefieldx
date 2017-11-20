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
package ffx.algorithms;

import java.io.File;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.crystal.Crystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.OpenMMForceFieldEnergy;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;

import static ffx.algorithms.Thermostat.kB;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Hernan V. Bernabe
 */
public class OpenMMMolecularDynamics extends MolecularDynamics {

    private static final Logger logger = Logger.getLogger(OpenMMMolecularDynamics.class.getName());

    private OpenMMForceFieldEnergy openMMForceFieldEnergy;
    private PointerByReference openMMContext;
    private PointerByReference openMMIntegrator;
    private PointerByReference openMMPositions;
    private PointerByReference openMMVelocities;
    private PointerByReference openMMForces;

    /**
     * Number of OpenMM MD steps per iteration.
     */
    private int intervalSteps;

    /**
     * Random number generation.
     */
    private Random random;

    /**
     * OpenMM Integrator definition.
     */
    private String integrator;
    private double frictionCoeff;
    private double collisionFreq;

    /**
     * Flag to indicate OpenMM MD iteractions are running.
     */
    private boolean running;
    private long time;

    /**
     * Constructs an OpenMMMolecularDynamics object, to perform molecular
     * dynamics using native OpenMM routines, avoiding the cost of communicating
     * coordinates, gradients, and energies back and forth across the PCI bus.
     *
     * @param assembly MolecularAssembly to operate on
     * @param openMMForceFieldEnergy OpenMMForceFieldEnergy Potential. Cannot be
     * any other type of Potential.
     * @param properties Associated properties
     * @param listener
     * @param thermostat May have to be slightly modified for native OpenMM
     * routines
     * @param integratorMD May have to be slightly modified for native OpenMM
     * routines
     */
    public OpenMMMolecularDynamics(MolecularAssembly assembly, OpenMMForceFieldEnergy openMMForceFieldEnergy,
            CompositeConfiguration properties, AlgorithmListener listener, Thermostats thermostat, Integrators integratorMD) {
        super(assembly, openMMForceFieldEnergy, properties, listener, thermostat, integratorMD);

        /**
         * Initialization specific to OpenMMMolecularDynamics
         */
        this.openMMForceFieldEnergy = openMMForceFieldEnergy;
        openMMForceFieldEnergy.addCOMMRemover(false);
        random = new Random();
        integrator = "VERLET";
        running = false;
    }

    /**
     * UNSUPPORTED: OpenMMMolecularDynamics is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     *
     * @param system
     * @param printFrequency
     */
    @Override
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        throw new UnsupportedOperationException(" OpenMMMolecularDynamics does not support extended system variables!");
    }

    /**
     * UNSUPPORTED: OpenMMMolecularDynamics is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void detachExtendedSystem() {
        throw new UnsupportedOperationException(" OpenMMMolecularDynamics does not support extended system variables!");
    }

    /**
     * takeSteps moves the simulation forward in time a user defined number of
     * steps and integrates the equations of motion for each step. This method
     * ensures that the algorithm reports back only when the time interval
     * (steps) specified by the user is completed.
     *
     * @param intervalSteps
     */
    private void takeSteps(int intervalSteps) {
        OpenMM_Integrator_step(openMMIntegrator, intervalSteps);
    }

    /**
     * openMM_Update obtains the state of the simulation from OpenMM, getting
     * positions and velocities back from the OpenMM data structure.
     *
     * @param i
     */
    private void openMM_Update(int i, boolean running) {
        openMMContext = openMMForceFieldEnergy.getContext();
        int infoMask = OpenMM_State_Positions + OpenMM_State_Velocities + OpenMM_State_Forces + OpenMM_State_Energy;

        PointerByReference state = OpenMM_Context_getState(openMMContext, infoMask, 0);
        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;

        openMMPositions = OpenMM_State_getPositions(state);
        openMMForceFieldEnergy.getOpenMMPositions(openMMPositions, numberOfVariables, x);

        openMMVelocities = OpenMM_State_getVelocities(state);
        openMMForceFieldEnergy.getOpenMMVelocities(openMMVelocities, numberOfVariables, v);

        double e = 0.0;
        for (int ii = 0; ii < v.length; ii++) {
            double velocity = v[ii];
            double v2 = velocity * velocity;
            e += mass[ii] * v2;
        }

        int dof = v.length;
        currentTemperature = e / (kB * dof);

        openMMForces = OpenMM_State_getForces(state);
        openMMForceFieldEnergy.getOpenMMAccelerations(openMMForces, numberOfVariables, mass, a);
        openMMForceFieldEnergy.getOpenMMAccelerations(openMMForces, numberOfVariables, mass, aPrevious);

        if (running) {
            if (i == 0) {
                logger.info(format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
                logger.info(format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
                logger.info(format("  %8s %12.4f %12.4f %12.4f %8.2f",
                        "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));
            } else if (i % printFrequency == 0) {
                double simTime = i * dt * 1.0e-3;
                logger.info(format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.2f",
                        simTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * NS2SEC));
            }

            if (saveSnapshotFrequency > 0 && i % (saveSnapshotFrequency * 1000 / intervalSteps) == 0 && i != 0) {
                for (AssemblyInfo ai : assemblies) {
                    if (ai.archiveFile != null && !saveSnapshotAsPDB) {
                        if (ai.xyzFilter.writeFile(ai.archiveFile, true)) {
                            logger.info(String.format(" Appended snap shot to %s", ai.archiveFile.getName()));
                        } else {
                            logger.warning(String.format(" Appending snap shot to %s failed", ai.archiveFile.getName()));
                        }
                    } else if (saveSnapshotAsPDB) {
                        if (ai.pdbFilter.writeFile(ai.pdbFile, false)) {
                            logger.info(String.format(" Wrote PDB file to %s", ai.pdbFile.getName()));
                        } else {
                            logger.warning(String.format(" Writing PDB file to %s failed.", ai.pdbFile.getName()));
                        }
                    }
                }
            }

            /**
             * Write out restart files every saveRestartFileFrequency steps.
             */
            if (restartFrequency > 0 && i % (restartFrequency * 1000 / intervalSteps ) == 0 && i != 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(String.format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(String.format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }
        }

        OpenMM_State_destroy(state);
    }

    @Override
    public void init(int numSteps, double timeStep, double printInterval, double saveInterval,
            String fileType, double restartFrequency, double temperature, boolean initVelocities, File dyn) {
        this.targetTemperature = temperature;
        this.dt = timeStep;
        this.printFrequency = (int) printInterval;
        this.restartFile = dyn;
        this.initVelocities = initVelocities;

        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = 100;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        /**
         * Convert save interval to a save frequency.
         */
        saveSnapshotFrequency = 1000;
        if (saveInterval >= this.dt) {
            saveSnapshotFrequency = (int) (saveInterval / this.dt);
        }

        done = false;

        assemblies.stream().parallel().forEach((ainfo) -> {
            MolecularAssembly mola = ainfo.getAssembly();
            CompositeConfiguration aprops = ainfo.props;
            File file = mola.getFile();
            String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
            File archFile = ainfo.archiveFile;
            if (archFile == null) {
                archFile = new File(filename + ".arc");
                ainfo.archiveFile = XYZFilter.version(archFile);
            }
            if (ainfo.pdbFile == null) {
                String extName = FilenameUtils.getExtension(file.getName());
                if (extName.toLowerCase().startsWith("pdb")) {
                    ainfo.pdbFile = file;
                } else {
                    ainfo.pdbFile = new File(filename + ".pdb");
                }
            }
            if (ainfo.xyzFilter == null) {
                ainfo.xyzFilter = new XYZFilter(file, mola, mola.getForceField(), aprops);
            }
            if (ainfo.pdbFilter == null) {
                ainfo.pdbFilter = new PDBFilter(ainfo.pdbFile, mola, mola.getForceField(), aprops);
            }
        });

        if (!initialized) {
            updateIntegrator();
        }

        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());

        if (dyn == null) {
            this.restartFile = new File(firstFileName + ".dyn");
            loadRestart = false;
        } else {
            this.restartFile = dyn;
            loadRestart = true;
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        if (!initialized) {
            if (loadRestart) {
                Crystal crystal = molecularAssembly.getCrystal();
                // possibly add check to see if OpenMM supports this space group.
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                } else {
                    molecularAssembly.getPotentialEnergy().setCrystal(crystal);
                    openMMForceFieldEnergy.setOpenMMPositions(x, numberOfVariables);
                    openMMForceFieldEnergy.setOpenMMVelocities(v, numberOfVariables);
                }
            } else {
                openMMForceFieldEnergy.loadFFXPositionToOpenMM();
                if (initVelocities) {
                    int randomNumber = random.nextInt();
                    OpenMM_Context_setVelocitiesToTemperature(openMMContext, temperature, randomNumber);
                }
            }
        }

        int i = 0;
        running = false;
        openMM_Update(i, running);
    }

    /**
     * Start sets up context, write out file name, restart file name, sets the
     * integrator and determines whether the simulation is starting out from a
     * previous molecular dynamics run (.dyn) or if the initial velocities are
     * determined by a Maxwell Boltzmann distribution. This method then calls
     * methods openMMUPdate and takeSteps to run the molecular dynamics
     * simulation.
     *
     * @param numSteps
     */
    @Override
    public void dynamic(int numSteps, double timeStep, double printInterval, double saveInterval, double temperature, boolean initVelocities, File dyn) {
        init(numSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency, temperature, initVelocities, dyn);

        storeState();
        if (intervalSteps == 0 || intervalSteps > numSteps) {
            intervalSteps = numSteps;
        }
        running = true;

        int i = 0;
        while (i < numSteps) {
            openMM_Update(i, running);
            time = -System.nanoTime();
            takeSteps(intervalSteps);
            time += System.nanoTime();
            i += intervalSteps;
        }
        openMM_Update(i, running);
        logger.info("");
    }

    /**
     * Method to set the Integrator string
     *
     * @param integrator string name of the integrator to be used during the
     * simulation
     */
    public void setIntegratorString(String integrator) {
        this.integrator = integrator;
    }

    /**
     * Method to set the coefficient of friction (for Langevin integrator)
     *
     * @param frictionCoeff coefficient of friction that acts on the atoms
     * during simulation
     */
    public void setFrictionCoefficient(double frictionCoeff) {
        this.frictionCoeff = frictionCoeff;
    }

    /**
     * Method to set the collision frequency (for Andersen thermostat)
     *
     * @param collisionFreq rate at which atoms collide in the Andersen
     * thermostat
     */
    public void setCollisionFrequency(double collisionFreq) {
        this.collisionFreq = collisionFreq;
    }

    public void setIntervalSteps(int intervalSteps) {
        this.intervalSteps = intervalSteps;
        logger.info(String.format(" Interval Steps set at %d", intervalSteps));
    }

    public void updateIntegrator() {
        openMMForceFieldEnergy.setIntegrator(integrator, dt, frictionCoeff, targetTemperature, collisionFreq);
        openMMIntegrator = openMMForceFieldEnergy.getIntegrator();
        openMMContext = openMMForceFieldEnergy.getContext();
    }

    /**
     * Returns the OpenMM DynamicsEngine
     * @return OPENMM
     */
    @Override
    public DynamicsEngine getEngine() {
        return DynamicsEngine.OPENMM;
    }
}
