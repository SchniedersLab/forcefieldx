/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocitiesToTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_resize;

import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.DYNFilter;
import static ffx.algorithms.thermostats.Thermostat.convert;
import static ffx.algorithms.thermostats.Thermostat.kB;
import java.util.Random;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Hernan V. Bernabe
 */
public class MolecularDynamicsOpenMM extends MolecularDynamics {

    private static final Logger logger = Logger.getLogger(MolecularDynamicsOpenMM.class.getName());

    /**
     * OpenMM ForceFieldEnergy.
     */
    private ForceFieldEnergyOpenMM forceFieldEnergyOpenMM;
    /**
     * OpenMM Context.
     */
    private PointerByReference context;
    /**
     * OpenMM Integrator.
     */
    private PointerByReference integrator;
    /**
     * OpenMM Positions.
     */
    private PointerByReference positions;
    /**
     * OpenMM Velocities.
     */
    private PointerByReference velocities;
    /**
     * OpenMM Forces.
     */
    private PointerByReference forces;
    /**
     * Number of OpenMM Particles (i.e. the number of FFX atoms).
     */
    private int numParticles;
    /**
     * Number of OpenMM Degrees of Freedom.
     */
    private int dof;
    /**
     * Integrator Type.
     */
    private final IntegratorEnum integratorType;
    /**
     * Thermostat Type.
     */
    private final ThermostatEnum thermostatType;
    /**
     * Integrator String.
     */
    private String integratorString;
    /**
     * Number of OpenMM MD steps per iteration.
     */
    private int intervalSteps;
    /**
     * Flag to indicate OpenMM MD interactions are running.
     */
    private boolean running;
    /**
     * Run time.
     */
    private long time;
    /**
     * Time elapsed for each molecular dynamics trajectory along with output of
     * the calculated energies.
     */
    private long mdTime = 0;
    /**
     * Time elapsed to run the init method (will probably be removed once
     * testing is complete).
     */
    private long initTime = 0;
    /**
     * Time elapsed for first update of OpenMM for a molecular dynamics move
     * (will probably be removed once testing is complete).
     */
    private long firstUpdateTime = 0;
    /**
     * Time elapsed to take steps to generate the molecular dynamics trajectory
     * (will probably be removed once testing is complete).
     */
    private long takeStepsTime = 0;
    /**
     * Time elapsed for the second update of OpenMM for a molecular dynamics
     * move (will probably be removed once testing is complete).
     */
    private long secondUpdateTime = 0;
    /**
     * Total energy at the end of a molecular dynamics move.
     */
    private double endTotalEnergy;
    /**
     * Total number of atoms in the system, typically used in for loops.
     */
    private int natoms;
    // A change in potential energy exceeding 1E6 kcal/mol triggers a warning and snapshot dump.
    private double defaultDeltaPEThresh = 1.0E6;
    /**
     * Boolean used to signify molecular dynamics in the NVE (micro canonical)
     * ensemble.
     */
    private boolean NVE = false;
    /**
     * Boolean used to signify molecular dynamics in the NPT
     * (isothermal-isobaric) ensemble.
     */
    private boolean NPT = false;
    /**
     * Boolean used to suppress logging for creation of integrator, typically
     * used to limit the amount of output to the screen.
     */
    private boolean quiet = true;
    /**
     * Double that holds the target pressure for the barostat under NPT
     * dynamics.
     */
    private double pressure;
    /**
     * Frequency of collisions for the barostat under NPT dynamics.
     */
    private int barostatFrequency;
    /**
     * Integer used to count the number of times the context has been set up,
     * mainly used to add a thermostat to the default integrator (Velocity
     * Verlet).
     */
    private int contextCounter = 0;
    /**
     * Composite properties for the system used in the simulation.
     */
    private CompositeConfiguration properties;
    /**
     * Array that holds the velocities generated by OpenMM under specific
     * temperature.
     */
    private PointerByReference velArray = null;
    /**
     * Random number generator used to psuedo seed the OpenMM velocity generator
     * method.
     */
    private Random random;
    /**
     * Boolean to signify that we are updating the system (post-MD move).
     */
    private boolean update = false;

    /**
     * Constructs an MolecularDynamicsOpenMM object, to perform molecular
     * dynamics using native OpenMM routines, avoiding the cost of communicating
     * coordinates, gradients, and energies back and forth across the PCI bus.
     *
     * @param assembly MolecularAssembly to operate on
     * @param forceFieldEnergyOpenMM ForceFieldEnergyOpenMM Potential. Cannot be
     * any other type of Potential.
     * @param properties Associated properties
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param thermostat May have to be slightly modified for native OpenMM
     * routines
     * @param integratorMD May have to be slightly modified for native OpenMM
     * routines
     */
    public MolecularDynamicsOpenMM(MolecularAssembly assembly, ForceFieldEnergyOpenMM forceFieldEnergyOpenMM,
            CompositeConfiguration properties, AlgorithmListener listener,
            ThermostatEnum thermostat, IntegratorEnum integratorMD) {
        super(assembly, forceFieldEnergyOpenMM, properties, listener, thermostat, integratorMD);

        /**
         * Initialization specific to MolecularDynamicsOpenMM
         */
        this.forceFieldEnergyOpenMM = forceFieldEnergyOpenMM;
        numParticles = forceFieldEnergyOpenMM.getNumParticles();
        forceFieldEnergyOpenMM.addCOMMRemover(false);
        thermostatType = thermostat;
        integratorType = integratorMD;
        running = false;
        integratorToString(integratorMD);
        this.properties = properties;

        random = new Random();
        if (properties.containsKey("velRandomSeed")) {
            random.setSeed(properties.getInt("velRandomSeed", 0));
        } else {
            random.setSeed(0);
        }

        updateContext();
    }

    /**
     * Get the ForceFieldEnergyOpenMM instance used to run MD.
     *
     * @return a {@link ffx.potential.ForceFieldEnergyOpenMM} object.
     */
    public ForceFieldEnergyOpenMM getForceFieldEnergyOpenMM() {
        return forceFieldEnergyOpenMM;
    }

    /**
     * takeOpenMMSteps moves the simulation forward in time a user defined
     * number of steps and integrates the equations of motion for each step.
     * This method ensures that the algorithm reports back only when the time
     * interval (steps) specified by the user is completed.
     *
     * @param intervalSteps
     */
    private void takeOpenMMSteps(int intervalSteps) {
        OpenMM_Integrator_step(integrator, intervalSteps);
    }

    /**
     * Get from OpenMM positions, velocities, forces and energy values.
     */
    private void getOpenMMState() {
        context = forceFieldEnergyOpenMM.getContext();

        int infoMask = OpenMM_State_Positions + OpenMM_State_Velocities + OpenMM_State_Forces + OpenMM_State_Energy;

        PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * convert / (kB * dof);

        long fetchPosVelForceTime = 0;
        fetchPosVelForceTime = -System.nanoTime();
        long fetchPos = 0;
        fetchPos = -System.nanoTime();
        positions = OpenMM_State_getPositions(state);
        forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles, x);
        fetchPos += System.nanoTime();
        logger.info(String.format(" Fetched positions in %6.3f", fetchPos * NS2SEC));

        long fetchVel = 0;
        fetchVel = -System.nanoTime();
        velocities = OpenMM_State_getVelocities(state);
        forceFieldEnergyOpenMM.getOpenMMVelocities(velocities, numParticles, v);
        fetchVel += System.nanoTime();
        logger.info(String.format(" Fetched velocities in %6.3f", fetchVel * NS2SEC));

        long fetchForce = 0;
        fetchForce = -System.nanoTime();
        forces = OpenMM_State_getForces(state);
        forceFieldEnergyOpenMM.getOpenMMAccelerations(forces, numParticles, mass, a);
        fetchForce += System.nanoTime();
        logger.info(String.format(" Fetched forces in %6.3f", fetchForce * NS2SEC));

        fetchPosVelForceTime += System.nanoTime();
        logger.info(String.format(" Time to fetch positions, velocities and forces from within update %6.3f", fetchPosVelForceTime * NS2SEC));

        if (aPrevious == null || aPrevious.length != a.length) {
            aPrevious = new double[a.length];
        }
        arraycopy(a, 0, aPrevious, 0, a.length);

        OpenMM_State_destroy(state);
    }

    private void getOpenMMEnergies() {
        context = forceFieldEnergyOpenMM.getContext();

        int infoMask = OpenMM_State_Energy;

        PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * convert / (kB * dof);

        OpenMM_State_destroy(state);
    }

    private void getOpenMMEnergiesAndPositions() {
        context = forceFieldEnergyOpenMM.getContext();

        int infoMask = OpenMM_State_Positions + OpenMM_State_Energy;

        PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * convert / (kB * dof);

        positions = OpenMM_State_getPositions(state);
        forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles, x);

        OpenMM_State_destroy(state);
    }

    /**
     * Set to OpenMM positions and velocities.
     */
    private void setOpenMMState(boolean setPositions, boolean setVelocities) {
        Atom atoms[] = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;

        if (setPositions) {
            if (x == null || x.length < nAtoms * 3) {
                logger.severe(" Position vector has not been allocated.");
            }

            for (int i = 0; i < nAtoms; i++) {
                int index = i * 3;
                Atom atom = atoms[i];
                x[index] = atom.getX();
                x[index + 1] = atom.getY();
                x[index + 2] = atom.getZ();
            }
            forceFieldEnergyOpenMM.setOpenMMPositions(x, numberOfVariables);
        }

        if (setVelocities) {
            double velocity[] = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                int index = i * 3;
                Atom atom = atoms[i];
                atom.getVelocity(velocity);
                v[index] = velocity[0];
                v[index + 1] = velocity[1];
                v[index + 2] = velocity[2];
            }
            forceFieldEnergyOpenMM.setOpenMMVelocities(v, numberOfVariables);
        }
    }

    /**
     * updateFromOpenMM obtains the state of the simulation from OpenMM,
     * completes some logging, and saves restart files.
     *
     * @param i
     */
    private void updateFromOpenMM(int i, boolean running) {

        double priorPE = currentPotentialEnergy;

        if (update) {
            getOpenMMEnergiesAndPositions();
        }
        /*else {
            getOpenMMState();
        } */
        //getOpenMMState();
        detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

        if (running) {
            if (i == 0) {
                logger.info(format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
                logger.info(format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
                logger.info(format("  %8s %12.4f %12.4f %12.4f %8.2f",
                        "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));
                startingKineticEnergy = currentKineticEnergy;
                startingTotalEnergy = currentTotalEnergy;
            } else if (i % printFrequency == 0) {
                double simTime = i * dt * 1.0e-3;
                time += System.nanoTime();
                mdTime = time;
                logger.info(format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.2f",
                        simTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * NS2SEC));

                endTotalEnergy = currentTotalEnergy;

                time = -System.nanoTime();
            }

            if (saveSnapshotFrequency > 0 && i % (saveSnapshotFrequency * 1000) == 0 && i != 0) {
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
            if (saveRestartFileFrequency > 0 && i % (saveRestartFileFrequency * 1000) == 0 && i != 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void init(int numSteps, double timeStep, double printInterval, double saveInterval,
            String fileType, double restartFrequency, double temperature, boolean initVelocities, File dyn) {
        this.targetTemperature = temperature;
        this.dt = timeStep;
        this.printFrequency = (int) printInterval;
        this.restartFile = dyn;
        this.initVelocities = initVelocities;

        //logger.info(String.format(" Target Temperature %f", targetTemperature));
        // Uncomment this to get code back to normal
        //updateContext();
        switch (thermostatType) {
            case BUSSI:
            case BERENDSEN:
                if (!integratorString.equalsIgnoreCase("LANGEVIN")) {
                    logger.info(String.format(" Replacing thermostat %s with OpenMM's Andersen thermostat", thermostatType));
                    forceFieldEnergyOpenMM.addAndersenThermostat(targetTemperature);
                    if (NPT) {
                        setMonteCarloBarostat(pressure, targetTemperature, barostatFrequency);
                        long updateContextTime = 0;
                        updateContextTime = -System.nanoTime();
                        updateContext();
                        updateContextTime += System.nanoTime();
                        logger.info(String.format(" Updated context in %6.3f seconds", updateContextTime * NS2SEC));
                    } else {
                        long updateContextTime = 0;
                        updateContextTime = -System.nanoTime();
                        updateContext();
                        updateContextTime += System.nanoTime();
                        logger.info(String.format(" Updated context in %6.3f seconds", updateContextTime * NS2SEC));
                    }
                } else {
                    logger.info(" Langevin/Stochastic dynamics already has temperature control, will not be adding thermostat!");
                }

                break;
            case ADIABATIC:
                if (integratorString.equalsIgnoreCase("LANGEVIN") && NPT) {
                    setMonteCarloBarostat(pressure, targetTemperature, barostatFrequency);
                    long updateContextTime = 0;
                    updateContextTime = -System.nanoTime();
                    updateContext();
                    updateContextTime += System.nanoTime();
                    logger.info(String.format(" Updated context in %6.3f seconds", updateContextTime * NS2SEC));
                }
            default:
                break;
            // No thermostat.
        }

        // Code used to have the updateContext method call here to ensure the context is updated if the user wished to add
        // a thermostat or a barostat. It has been moved inside the different cases of the switch statement to help the 
        // performance of the the MCOSRW algorithm which does not require a thermostat
        
        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = 100;
        printInterval *= 1000; // Time step is in fsec, so convert printInterval to fsec.
        if (printInterval >= dt) {
            printFrequency = (int) (printInterval / dt);
        }

        /**
         * Convert save interval to a save frequency.
         */
        saveSnapshotFrequency = 1000;
        // saveInterval *= 1000; // Time step is in fsec, so convert saveInterval to fsec.
        if (saveInterval >= dt) {
            saveSnapshotFrequency = (int) (saveInterval / dt);
        }

        /**
         * Convert restart interval to frequency.
         */
        saveRestartFileFrequency = 1000;
        // restartFrequency *= 1000; // Time step is in fsec, so converting restartFrequency
        if (restartFrequency >= dt) {
            saveRestartFileFrequency = (int) (restartFrequency / dt);
        }

        done = false;

        assemblyInfo();

        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());

        if (dyn == null) {
            restartFile = new File(firstFileName + ".dyn");
            loadRestart = false;
        } else {
            restartFile = dyn;
            loadRestart = true;
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        dof = forceFieldEnergyOpenMM.calculateDegreesOfFreedom();

        boolean setPositions = true;
        boolean setVelocities = true;

        if (!initialized) {
            if (loadRestart) {
                Crystal crystal = molecularAssembly.getCrystal();
                // possibly add check to see if OpenMM supports this space group.
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                } else {
                    logger.info(format(" Continuing from " + dyn.getAbsolutePath()));
                    forceFieldEnergyOpenMM.setCrystal(crystal);
                    // Load positions into the main FFX data structure, move into primary unit cell, then load to OpenMM.
                    Atom[] atoms = molecularAssembly.getAtomArray();
                    double[] xyz = new double[3];
                    double[] vel = new double[3];
                    double[] acc = new double[3];
                    double[] accPrev = new double[3];
                    for (int i = 0; i < atoms.length; i++) {
                        Atom atom = atoms[i];
                        int i3 = i * 3;
                        for (int j = 0; j < 3; j++) {
                            xyz[j] = x[i3 + j];
                            vel[j] = v[i3 + j];
                            acc[j] = a[i3 + j];
                            accPrev[j] = aPrevious[i3 + j];
                        }
                        atom.setXYZ(xyz);
                        atom.setVelocity(vel);
                        atom.setAcceleration(acc);
                        atom.setPreviousAcceleration(accPrev);
                    }
                    molecularAssembly.moveAllIntoUnitCell();
                }
            } else {
                if (initVelocities) {
                    if (properties.containsKey("openMMInitVel")) {
                        int randomSeed = random.nextInt();
                        long ommSetVelTime = 0;
                        ommSetVelTime = -System.nanoTime();
                        OpenMM_Context_setVelocitiesToTemperature(context, targetTemperature, randomSeed);
                        logger.info(String.format(" OpenMM set velocities to target temperature %f with random seed %d",
                                targetTemperature, randomSeed));
                        ommSetVelTime += System.nanoTime();
                        logger.info(String.format(" Set velocites with OpenMM in %6.3f", ommSetVelTime * NS2SEC));
                        setVelocities = false;
                    } else {
                        getThermostat().setQuiet(true);
                        getThermostat().maxwell(targetTemperature);
                        Atom[] atoms = molecularAssembly.getAtomArray();
                        double[] vel = new double[3];
                        long ffxSetVelTime = 0;
                        ffxSetVelTime = -System.nanoTime();
                        for (int i = 0; i < atoms.length; i++) {
                            Atom atom = atoms[i];
                            int i3 = i * 3;
                            for (int j = 0; j < 3; j++) {
                                vel[j] = v[i3 + j];
                            }
                            atom.setVelocity(vel);
                        }
                        ffxSetVelTime += System.nanoTime();
                        logger.info(String.format(" Set velocities with FFX in %6.3f", ffxSetVelTime * NS2SEC));
                    }
                }
            }

            long setPosVel = 0;
            setPosVel = -System.nanoTime();
            setOpenMMState(setPositions, setVelocities);
            setPosVel += System.nanoTime();
            logger.info(String.format(" Set positions and velocities in %6.3f seconds", setPosVel * NS2SEC));

            // Call to retrieve the starting kinetic energy for the system.
            long retrieveEnergyTime = 0;
            retrieveEnergyTime = -System.nanoTime();
            getOpenMMEnergies();
            retrieveEnergyTime += System.nanoTime();
            logger.info(String.format(" Retrieved energies in %6.3f seconds", retrieveEnergyTime * NS2SEC));
            startingKineticEnergy = currentKineticEnergy;
        }

        saveSnapshotAsPDB = true;
        if (fileType.equalsIgnoreCase("XYZ")) {
            saveSnapshotAsPDB = false;
        } else if (!fileType.equalsIgnoreCase("PDB")) {
            logger.warning(" Snapshot file type unrecognized; saving snaphshots as PDB.\n");
        }

        Atom[] atoms = molecularAssembly.getAtomArray();
        natoms = atoms.length;

        int i = 0;
        running = false;

        // logger.info(" Calling OpenMM Update from MD Init.");
        //updateFromOpenMM(i, running);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Start sets up context, write out file name, restart file name, sets the
     * integrator and determines whether the simulation is starting out from a
     * previous molecular dynamics run (.dyn) or if the initial velocities are
     * determined by a Maxwell Boltzmann distribution. This method then calls
     * methods openMMUpdate and takeOpenMMSteps to run the molecular dynamics
     * simulation.
     */
    @Override
    public void dynamic(int numSteps, double timeStep, double printInterval, double saveInterval, double temperature, boolean initVelocities, File dyn) {

        initTime = -System.nanoTime();
        init(numSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency, temperature, initVelocities, dyn);
        initTime += System.nanoTime();

        logger.info(String.format("\n Initialized system in %6.3f sec.", initTime * NS2SEC));

        storeState();

        if (intervalSteps == 0 || intervalSteps > numSteps) {
            intervalSteps = numSteps;
        }
        running = true;

        // Update the time step in Picoseconds.
        OpenMM_Integrator_setStepSize(integrator, dt * 1.0e-3);

        int i = 0;
        time = -System.nanoTime();
        while (i < numSteps) {
            // Get an update from OpenMM.
            /**
             * firstUpdateTime = -System.nanoTime(); updateFromOpenMM(i,
             * running); firstUpdateTime += System.nanoTime();
             * logger.info(String.format("\n First update finished in %6.3f",
             * firstUpdateTime * NS2SEC));
             */

            updateFromOpenMM(i, running);

            // Take MD steps in OpenMM.
            takeStepsTime = -System.nanoTime();
            takeOpenMMSteps(intervalSteps);
            takeStepsTime += System.nanoTime();
            logger.info(String.format("\n Took steps in %6.3f", takeStepsTime * NS2SEC));

            // Update the total step count.
            i += intervalSteps;

            update = true;
            secondUpdateTime = -System.nanoTime();
            updateFromOpenMM(i, running);
            secondUpdateTime += System.nanoTime();
            logger.info(String.format("\n Update finished in %6.3f", secondUpdateTime * NS2SEC));
            update = false;
        }
        /**
         * secondUpdate = true; secondUpdateTime = -System.nanoTime();
         * updateFromOpenMM(i, running); secondUpdateTime += System.nanoTime();
         * endLoopTime += System.nanoTime(); logger.info(String.format("\n
         * Second update finished in %6.3f", secondUpdateTime * NS2SEC));
         * logger.info(String.format("\n Time from beginning of loop to end of
         * method is %6.3f", endLoopTime * NS2SEC)); secondUpdate = false;
         */
    }

    /**
     * <p>
     * integratorToString.</p>
     *
     * @param integrator a {@link ffx.algorithms.integrators.IntegratorEnum}
     * object.
     */
    public final void integratorToString(IntegratorEnum integrator) {
        if (integrator == null) {
            integratorString = "VERLET";
            logger.info(String.format(" No specified integrator, will use Verlet"));
        } else {
            switch (integratorType) {
                case STOCHASTIC:
                    integratorString = "LANGEVIN";
                    break;
                case VELOCITYVERLET:
                    integratorString = "VERLET";
                    break;
                case RESPA:
                    integratorString = "RESPA";
                    //logger.info(String.format(" In RESPA integrator case"));
                    break;
                default:
                    integratorString = "VERLET";
                    logger.warning(String.format(" Integrator %s incompatible with "
                            + "OpenMM MD integration; defaulting to %s", integratorType, integratorString));
                    break;
            }
        }
    }

    /**
     * <p>
     * setMonteCarloBarostat.</p>
     *
     * @param pressure a double.
     * @param temperature a double.
     * @param frequency a int.
     */
    public void setMonteCarloBarostat(double pressure, double temperature, int frequency) {
        forceFieldEnergyOpenMM.addMonteCarloBarostat(pressure, temperature, frequency);
    }

    /**
     * <p>
     * Setter for the field <code>pressure</code>.</p>
     *
     * @param pressure a double.
     */
    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    /**
     * <p>
     * Setter for the field <code>barostatFrequency</code>.</p>
     *
     * @param barostatFrequency a int.
     */
    public void setBarostatFrequency(int barostatFrequency) {
        this.barostatFrequency = barostatFrequency;
    }

    /**
     * <p>
     * setNPTDynamics.</p>
     */
    public void setNPTDynamics() {
        NPT = true;
    }

    /**
     * <p>
     * Setter for the field <code>intervalSteps</code>.</p>
     *
     * @param intervalSteps a int.
     */
    public void setIntervalSteps(int intervalSteps) {
        this.intervalSteps = intervalSteps;
        logger.info(String.format(" Interval Steps set at %d", intervalSteps));
    }

    /**
     * <p>
     * updateContext.</p>
     */
    public final void updateContext() {
        String currentIntegrator = forceFieldEnergyOpenMM.getIntegratorString();
        double currentTimeStp = forceFieldEnergyOpenMM.getTimeStep();
        double currentTemperature = forceFieldEnergyOpenMM.getTemperature();

        //logger.info(String.format(" Outside if statement"));
        //logger.info(String.format(" Counter is %d", counter));
        if (currentTemperature != targetTemperature || currentTimeStp != dt || !currentIntegrator.equalsIgnoreCase(integratorString) || (currentIntegrator.equalsIgnoreCase("VERLET") && contextCounter != 0)) {
            if (!quiet) {
                logger.info(String.format(" Creating OpenMM Context with step size %8.3f and target temperature %8.3f.", dt, targetTemperature));
            }
            logger.info(" Creating new OpenMM Context");
            long contextTime = 0;
            contextTime = -System.nanoTime();
            forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature);
            contextTime += System.nanoTime();
            logger.info(String.format("Created new context in %6.3f", contextTime * NS2SEC));
            integrator = forceFieldEnergyOpenMM.getIntegrator();
            context = forceFieldEnergyOpenMM.getContext();
        } else {
            integrator = forceFieldEnergyOpenMM.getIntegrator();
            context = forceFieldEnergyOpenMM.getContext();
        }
        quiet = false;
        contextCounter++;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns the OpenMM DynamicsEngine
     */
    @Override
    public DynamicsEngine getEngine() {
        return DynamicsEngine.OPENMM;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns time spent calculating the molecular dynamics trajectory on the
     * GPU
     */
    @Override
    public long getMDTime() {
        return mdTime;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getStartingTotalEnergy() {
        return startingTotalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getEndTotalEnergy() {
        return endTotalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTimeStep() {
        return dt;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getIntervalSteps() {
        return intervalSteps;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumAtoms() {
        return natoms;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    /**
     * {@inheritDoc}
     * <p>
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }

    /**
     * {@inheritDoc}
     * <p>
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void detachExtendedSystem() {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }
}
