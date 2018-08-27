/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setStepSize;
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

    private long mdTime = 0;

    private double startingTotalEnergy;

    private double endTotalEnergy;

    private int natoms;
    // A change in potential energy exceeding 1E6 kcal/mol triggers a warning and snapshot dump.
    private double defaultDeltaPEThresh = 1.0E6;

    private boolean NVE = false;
    
    private boolean NPT = false;

    private boolean quiet = true;
    
    private double pressure;
    
    private int barostatFrequency;

    /**
     * Constructs an MolecularDynamicsOpenMM object, to perform molecular
     * dynamics using native OpenMM routines, avoiding the cost of communicating
     * coordinates, gradients, and energies back and forth across the PCI bus.
     *
     * @param assembly MolecularAssembly to operate on
     * @param forceFieldEnergyOpenMM ForceFieldEnergyOpenMM Potential. Cannot be
     * any other type of Potential.
     * @param properties Associated properties
     * @param listener
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
        updateContext();
    }

    /**
     * Get the ForceFieldEnergyOpenMM instance used to run MD.
     * @return
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

        positions = OpenMM_State_getPositions(state);
        forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles, x);

        velocities = OpenMM_State_getVelocities(state);
        forceFieldEnergyOpenMM.getOpenMMVelocities(velocities, numParticles, v);

        forces = OpenMM_State_getForces(state);
        forceFieldEnergyOpenMM.getOpenMMAccelerations(forces, numParticles, mass, a);

        if (aPrevious == null || aPrevious.length != a.length) {
            aPrevious = new double[a.length];
        }
        arraycopy(a, 0, aPrevious, 0, a.length);

        OpenMM_State_destroy(state);
    }

    /**
     * Set to OpenMM positions and velocities.
     */
    private void setOpenMMState() {
        Atom atoms[] = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        if (x == null || x.length < nAtoms * 3) {
            logger.severe(" Position vector has not been allocated.");
        }
        double velocity[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            int index = i * 3;
            Atom atom = atoms[i];
            x[index] = atom.getX();
            x[index + 1] = atom.getY();
            x[index + 2] = atom.getZ();
            atom.getVelocity(velocity);
            v[index] = velocity[0];
            v[index + 1] = velocity[1];
            v[index + 2] = velocity[2];
        }
        forceFieldEnergyOpenMM.setOpenMMPositions(x, numberOfVariables);
        forceFieldEnergyOpenMM.setOpenMMVelocities(v, numberOfVariables);
    }

    /**
     * updateFromOpenMM obtains the state of the simulation from OpenMM,
     * completes some logging, and saves restart files.
     *
     * @param i
     */
    private void updateFromOpenMM(int i, boolean running) {

        double priorPE = currentPotentialEnergy;
        getOpenMMState();
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
            if (restartFrequency > 0 && i % (restartFrequency * 1000) == 0 && i != 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }
        }
    }

    @Override
    public void init(int numSteps, double timeStep, double printInterval, double saveInterval,
            String fileType, double restartFrequency, double temperature, boolean initVelocities, File dyn) {
        this.targetTemperature = temperature;
        this.dt = timeStep;
        this.printFrequency = (int) printInterval;
        this.restartFile = dyn;
        this.initVelocities = initVelocities;

        //logger.info(String.format(" Target Temperature %f", targetTemperature));
        updateContext();

        switch (thermostatType) {
            case BUSSI:
            case BERENDSEN:
                logger.info(String.format(" Replacing thermostat %s with OpenMM's Andersen thermostat", thermostatType));
                forceFieldEnergyOpenMM.addAndersenThermostat(targetTemperature);
                if (NPT){
                    setMonteCarloBarostat(pressure, targetTemperature, barostatFrequency);
                }
                break;
            case ADIABATIC:
                if (integratorString.equalsIgnoreCase("LANGEVIN") && NPT){
                    setMonteCarloBarostat(pressure, targetTemperature, barostatFrequency);
                }
            default:
                break;
            // No thermostat.
            }

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
        if (saveInterval >= this.dt) {
            saveSnapshotFrequency = (int) (saveInterval / this.dt);
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
                    getThermostat().setQuiet(true);
                    getThermostat().maxwell(targetTemperature);
                    Atom[] atoms = molecularAssembly.getAtomArray();
                    double[] vel = new double[3];
                    for (int i = 0; i < atoms.length; i++) {
                        Atom atom = atoms[i];
                        int i3 = i * 3;
                        for (int j = 0; j < 3; j++) {
                            vel[j] = v[i3 + j];
                        }
                        atom.setVelocity(vel);
                    }
                }
            }
            setOpenMMState();

            // Get the OpenMM State and then set the starting kinetic energy.
            getOpenMMState();
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
        updateFromOpenMM(i, running);
    }

    /**
     * Start sets up context, write out file name, restart file name, sets the
     * integrator and determines whether the simulation is starting out from a
     * previous molecular dynamics run (.dyn) or if the initial velocities are
     * determined by a Maxwell Boltzmann distribution. This method then calls
     * methods openMMUpdate and takeOpenMMSteps to run the molecular dynamics
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

        // Update the time step in Picoseconds.
        OpenMM_Integrator_setStepSize(integrator, dt * 1.0e-3);

        int i = 0;
        time = -System.nanoTime();
        while (i < numSteps) {
            // Get an update from OpenMM.
            updateFromOpenMM(i, running);

            // Take MD steps in OpenMM.
            takeOpenMMSteps(intervalSteps);

            // Update the total step count.
            i += intervalSteps;
        }
        updateFromOpenMM(i, running);
    }

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
    
    public void setMonteCarloBarostat(double pressure, double temperature, int frequency){
        forceFieldEnergyOpenMM.addMonteCarloBarostat(pressure, temperature, frequency);
    }
    
    public void setPressure(double pressure){
        this.pressure = pressure;
    }
    
    public void setBarostatFrequency(int barostatFrequency){
        this.barostatFrequency = barostatFrequency;
    }
    
    public void setNPTDynamics(){
        NPT = true;
    }

    public void setIntervalSteps(int intervalSteps) {
        this.intervalSteps = intervalSteps;
        logger.info(String.format(" Interval Steps set at %d", intervalSteps));
    }

    public final void updateContext() {
        String currentIntegrator = forceFieldEnergyOpenMM.getIntegratorString();
        double currentTimeStp = forceFieldEnergyOpenMM.getTimeStep();
        double currentTemperature = forceFieldEnergyOpenMM.getTemperature();
        if (currentTemperature != targetTemperature || currentTimeStp != dt || !currentIntegrator.equalsIgnoreCase(integratorString)) {
            if (!quiet) {
                logger.info(String.format(" Creating OpenMM Context with step size %8.3f and target temperature %8.3f.", dt, targetTemperature));
            }
            forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature);
            integrator = forceFieldEnergyOpenMM.getIntegrator();
            context = forceFieldEnergyOpenMM.getContext();
        } else {
            integrator = forceFieldEnergyOpenMM.getIntegrator();
            context = forceFieldEnergyOpenMM.getContext();
        }
        quiet = false;
    }

    /**
     * Returns the OpenMM DynamicsEngine
     *
     * @return OPENMM
     */
    @Override
    public DynamicsEngine getEngine() {
        return DynamicsEngine.OPENMM;
    }

    /**
     * Returns time spent calculating the molecular dynamics trajectory on the
     * GPU
     *
     * @return
     */
    @Override
    public long getMDTime() {
        return mdTime;
    }

    @Override
    public double getStartingTotalEnergy() {
        return startingTotalEnergy;
    }

    @Override
    public double getEndTotalEnergy() {
        return endTotalEnergy;
    }

    @Override
    public double getTimeStep() {
        return dt;
    }

    @Override
    public int getIntervalSteps() {
        return intervalSteps;
    }

    @Override
    public int getNumAtoms() {
        return natoms;
    }

    @Override
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    /**
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     *
     * @param system
     * @param printFrequency
     */
    @Override
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }

    /**
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void detachExtendedSystem() {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }
}
