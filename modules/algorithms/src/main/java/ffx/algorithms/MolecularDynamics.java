/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.algorithms;

import java.io.File;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;

/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class MolecularDynamics implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    private final MolecularAssembly molecularAssembly;
    private final Potential potential;
    private final CompositeConfiguration properties;
    private AlgorithmListener algorithmListener;
    private Thermostat thermostat;
    private Integrator integrator;
    private File archiveFile = null;
    private File restartFile = null;
    private File pdbFile = null;
    private XYZFilter xyzFilter = null;
    private DYNFilter dynFilter = null;
    private PDBFilter pdbFilter = null;
    private int printFrequency = 100;
    private int saveSnapshotFrequency = 1000;
    private int removeCOMMotionFrequency = 100;
    private boolean initVelocities = true;
    private boolean loadRestart = false;
    private boolean initialized = false;
    private boolean done = true;
    private boolean terminate = false;
    private final int numberOfVariables;
    private final double[] x;
    private final double[] v;
    private final double[] a;
    private final double[] aPrevious;
    private final double[] grad;
    private final double[] mass;
    private int nSteps = 1000;
    private double targetTemperature = 298.15;
    private double dt = 1.0;
    private double currentTemperature;
    private double currentKineticEnergy;
    private double currentPotentialEnergy;
    private double currentTotalEnergy;
    private boolean saveSnapshotAsPDB = true;
    private int saveRestartFileFrequency = 1000;
    private String fileType = "PDB";
    private double restartFrequency = 0.1;

    /**
     * <p>Constructor for MolecularDynamics.</p>
     *
     * @param assembly a {@link ffx.potential.bonded.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     * {@link ffx.algorithms.Thermostat.Thermostats} object.
     * @param requestedIntegrator
     */
    public MolecularDynamics(MolecularAssembly assembly,
            Potential potentialEnergy,
            CompositeConfiguration properties,
            AlgorithmListener listener,
            Thermostats requestedThermostat,
            Integrators requestedIntegrator) {
        this.molecularAssembly = assembly;
        this.algorithmListener = listener;
        this.potential = potentialEnergy;

        this.properties = properties;
        mass = potentialEnergy.getMass();
        numberOfVariables = potentialEnergy.getNumberOfVariables();
        x = new double[numberOfVariables];
        v = new double[numberOfVariables];
        a = new double[numberOfVariables];
        aPrevious = new double[numberOfVariables];
        grad = new double[numberOfVariables];

        /**
         * If an Integrator wasn't passed to the MD constructor, check for one
         * specified as a property.
         */
        if (requestedIntegrator == null) {
            String integrate = properties.getString("integrate", "beeman").trim();
            try {
                requestedIntegrator = Integrators.valueOf(integrate);
            } catch (Exception e) {
                requestedIntegrator = Integrators.BEEMAN;
            }

        }
        switch (requestedIntegrator) {
            case RESPA:
                integrator = new Respa(numberOfVariables, x, v, a, aPrevious, mass);
                break;
            case STOCHASTIC:
                double friction = properties.getDouble("friction", 91.0);
                Stochastic stochastic = new Stochastic(friction, numberOfVariables, x, v, a, mass);
                if (properties.containsKey("randomseed")) {
                    stochastic.setRandomSeed(properties.getInt("randomseed", 0));
                }
                integrator = stochastic;
                /**
                 * The stochastic dynamics integration procedure will thermostat
                 * the system. The ADIABTIC thermostat just serves to report the
                 * temperature and initialize velocities if necessary.
                 */
                requestedThermostat = Thermostats.ADIABATIC;
                break;
            case VELOCITYVERLET:
                integrator = new VelocityVerlet(numberOfVariables, x, v, a, mass);
                break;
            case BEEMAN:
            default:
                integrator = new BetterBeeman(numberOfVariables, x, v, a, aPrevious, mass);
        }


        /**
         * If a Thermostat wasn't passed to the MD constructor, check for one
         * specified as a property.
         */
        if (requestedThermostat == null) {
            String thermo = properties.getString("thermostat", "Berendsen").trim();
            try {
                requestedThermostat = Thermostats.valueOf(thermo);
            } catch (Exception e) {
                requestedThermostat = Thermostats.BERENDSEN;
            }
        }

        switch (requestedThermostat) {
            case BERENDSEN:
                double tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Berendsen(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), 300.0, tau);
                break;
            case BUSSI:
                tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Bussi(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), 300.0, tau);
                break;
            case ADIABATIC:
            default:
                thermostat = new Adiabatic(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes());
        }

        if (properties.containsKey("randomseed")) {
            thermostat.setRandomSeed(properties.getInt("randomseed", 0));
        }

        /**
         * For StochasticDynamics, center of mass motion will not be removed.
         */
        if (integrator instanceof Stochastic) {
            thermostat.removingCenterOfMassMotion(false);
        }

        done = true;
    }

    /**
     * <p>Setter for the field
     * <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link ffx.algorithms.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    /**
     * <p>Getter for the field
     * <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
    }

    /**
     * <p>Setter for the field
     * <code>x</code>.</p>
     *
     * @param x a double array to set the current parameters to.
     */
    public void setParameters(double x[]) {
        System.arraycopy(x, 0, this.x, 0, numberOfVariables);
    }

    /**
     * <p>Getter for the field
     * <code>x</code>.</p>
     *
     * @return a double array with the current parameters
     */
    public double[] getParameters() {
        return x;
    }

    /**
     * <p>Setter for the field
     * <code>archiveFile</code>.</p>
     *
     * @param archive a {@link java.io.File} object.
     */
    public void setArchiveFile(File archive) {
        this.archiveFile = archive;
    }

    /**
     * <p>Getter for the field
     * <code>archiveFile</code>.</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getArchiveFile() {
        return archiveFile;
    }

    /**
     * <p>init</p>
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param fileType
     * @param restartFrequency
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final String fileType, final double restartFrequency,
            final double temperature, final boolean initVelocities, final File dyn) {

        /**
         * Return if already running.
         */
        if (!done) {
            logger.warning(" Programming error - attempt to modify parameters of a running MolecularDynamics instance.");
            return;
        }

        this.nSteps = nSteps;
        /**
         * Convert the time step from femtoseconds to picoseconds.
         */
        dt = timeStep * 1.0e-3;

        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = 100;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        /**
         * Convert the save interval to a save frequency.
         */
        saveSnapshotFrequency = 1000;
        if (saveInterval >= this.dt) {
            saveSnapshotFrequency = (int) (saveInterval / this.dt);
        }

        /**
         * Set snapshot file type.
         */
        saveSnapshotAsPDB = true;
        if (fileType.equals("XYZ")) {
            saveSnapshotAsPDB = false;
        } else if (!fileType.equals("PDB")) {
            logger.warning("Snapshot file type unrecognized; saving snaphshots as PDB.\n");
        }

        /**
         * Convert restart frequency to steps.
         */
        saveRestartFileFrequency = 1000;
        if (restartFrequency >= this.dt) {
            saveRestartFileFrequency = (int) (restartFrequency / this.dt);
        }

        File file = molecularAssembly.getFile();
        String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
        if (archiveFile == null) {
            archiveFile = new File(filename + ".arc");
            archiveFile = XYZFilter.version(archiveFile);
        }

        if (dyn == null) {
            this.restartFile = new File(filename + ".dyn");
            loadRestart = false;
        } else {
            this.restartFile = dyn;
            loadRestart = true;
        }

        if (xyzFilter == null) {
            xyzFilter = new XYZFilter(file, molecularAssembly,
                    molecularAssembly.getForceField(), properties);
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        if (pdbFilter == null) {
            pdbFile = new File(filename + "_dyn.pdb");
            pdbFilter = new PDBFilter(new File(filename + "_dyn.pdb"), molecularAssembly, null, null);
        }

        this.targetTemperature = temperature;
        this.initVelocities = initVelocities;
    }

    /**
     * A version of init with the original method header. Redirects to the new
     * method with default values for added parameters. Needed by (at least)
     * ReplicaExchange, which calls this directly.
     * @param nSteps
     * @param timeStep
     * @param printInterval
     * @param saveInterval
     * @param temperature
     * @param initVelocities
     * @param dyn
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            final File dyn) {
        init(nSteps, timeStep, printInterval, saveInterval, "PDB", 0.1, temperature, initVelocities, dyn);
    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is
     * done.
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param fileType a String (either XYZ or PDB).
     * @param restartFrequency a double specifying the restart frequency.
     * @param dyn a {@link java.io.File} object.
     */
    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            String fileType, double restartFrequency, final File dyn) {
        this.fileType = fileType;
        this.restartFrequency = restartFrequency;
        dynamic(nSteps, timeStep, printInterval, saveInterval, temperature,
                initVelocities, dyn);
    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is
     * done.
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            final File dyn) {

        /**
         * Return if already running; Could happen if two threads call dynamic
         * on the same MolecularDynamics instance.
         */
        if (!done) {
            logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
            return;
        }

        if (integrator instanceof Stochastic) {
            logger.info(format("\n Stochastic dynamics in the NVT ensemble\n"));
        } else if (!(thermostat instanceof Adiabatic)) {
            logger.info(format("\n Molecular dynamics in the NVT ensemble\n"));
        } else {
            logger.info(format("\n Molecular dynamics in the NVE ensemble\n"));
        }

        init(nSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency, temperature, initVelocities, dyn);

        done = false;

        if (dyn != null) {
            logger.info(format(" Continuing from " + dyn.getAbsolutePath()));
        }
        logger.info(String.format(" Number of steps: %8d", nSteps));
        logger.info(String.format(" Time step:       %8.3f (fsec)", timeStep));
        logger.info(String.format(" Print interval:  %8.3f (psec)", printInterval));
        logger.info(String.format(" Save interval:   %8.3f (psec)", saveInterval));
        logger.info(String.format(" Archive file: %s", archiveFile.getName()));
        logger.info(String.format(" Restart file: %s", restartFile.getName()));

        Thread dynamicThread = new Thread(this);
        dynamicThread.start();
        synchronized (this) {
            try {
                while (dynamicThread.isAlive()) {
                    wait(100);
                }
            } catch (InterruptedException e) {
                String message = " Molecular dynamics interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }
    }

    /**
     * Methods to set file type and restartFrequency from groovy scripts
     */
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    public void setRestartFrequency(double restartFrequency) {
        this.restartFrequency = restartFrequency;
    }

    /**
     * Set the number of time steps between removal of center of mass kinetic
     * energy.
     *
     * @param removeCOMMotionFrequency Number of time steps between center of
     * mass removal.
     */
    public void setRemoveCOMMotionFrequency(int removeCOMMotionFrequency) {
        if (removeCOMMotionFrequency < 0) {
            removeCOMMotionFrequency = 0;
        }
        this.removeCOMMotionFrequency = removeCOMMotionFrequency;
        if (removeCOMMotionFrequency != 0) {
            thermostat.removingCenterOfMassMotion(true);
        } else {
            thermostat.removingCenterOfMassMotion(false);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        done = false;
        terminate = false;

        /**
         * Set the target temperature.
         */
        thermostat.setTargetTemperature(targetTemperature);
        if (integrator instanceof Stochastic) {
            Stochastic stochastic = (Stochastic) integrator;
            stochastic.setTemperature(targetTemperature);
        }

        /**
         * Set the step size.
         */
        integrator.setTimeStep(dt);

        if (!initialized) {
            /**
             * Initialize from a restart file.
             */
            if (loadRestart) {
                if (!dynFilter.readDYN(restartFile, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    return;
                }
            } else {
                /**
                 * Initialize from using current atomic coordinates.
                 */
                potential.getCoordinates(x);
                /**
                 * Initialize atomic velocities from a Maxwell-Boltzmann
                 * distribution or set to 0.
                 */
                if (initVelocities) {
                    thermostat.maxwell(targetTemperature);
                } else {
                    Arrays.fill(v, 0.0);
                }
            }
        } else {
            /**
             * If MD has already been run (ie. Annealing or RepEx), then
             * initialize velocities if requested.
             */
            if (initVelocities) {
                thermostat.maxwell(targetTemperature);
            }
        }

        /**
         * Compute the current potential energy.
         */
        potential.setScaling(null);
        currentPotentialEnergy = potential.energyAndGradient(x, grad);

        /**
         * Compute the current kinetic energy.
         */
        thermostat.kineticEnergy();
        currentKineticEnergy = thermostat.getKineticEnergy();
        currentTemperature = thermostat.getCurrentTemperature();
        currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

        /**
         * Initialize current and previous accelerations.
         */
        if (!initialized) {
            if (!loadRestart) {
                for (int i = 0; i < numberOfVariables; i++) {
                    a[i] = -Thermostat.convert * grad[i] / mass[i];
                }
                if (aPrevious != null) {
                    System.arraycopy(a, 0, aPrevious, 0, numberOfVariables);
                }
            }
            initialized = true;
        }

        logger.info(String.format("\n      Time      Kinetic    Potential        Total     Temp      CPU"));
        logger.info(String.format("      psec     kcal/mol     kcal/mol     kcal/mol        K      sec\n"));
        logger.info(String.format("          %13.4f%13.4f%13.4f %8.2f ", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));

        /**
         * Integrate Newton's equations of motion for the requested number of
         * steps, unless early termination is requested.
         */
        long time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {

            /**
             * Do the half-step thermostat operation.
             */
            thermostat.halfStep(dt);

            /**
             * Do the half-step integration operation.
             */
            integrator.preForce(potential);

            /**
             * Compute the potential energy and gradients.
             */
            currentPotentialEnergy = potential.energyAndGradient(x, grad);

            /**
             * Add the potential energy of the slow degrees of freedom.
             */
            if (integrator instanceof Respa) {
                Respa r = (Respa) integrator;
                currentPotentialEnergy += r.getHalfStepEnergy();
            }

            /**
             * Do the full-step integration operation.
             */
            integrator.postForce(grad);

            /**
             * Compute the full-step kinetic energy.
             */
            thermostat.kineticEnergy();

            /**
             * Do the full-step thermostat operation.
             */
            thermostat.fullStep(dt);

            /**
             * Recompute the kinetic energy after the full-step thermostat
             * operation.
             */
            thermostat.kineticEnergy();

            /**
             * Remove center of mass motion ever ~100 steps.
             */
            if (thermostat.removingCOM() && step % removeCOMMotionFrequency == 0) {
                thermostat.centerOfMassMotion(true, false);
            }

            /**
             * Collect current kinetic energy, temperature, and total energy.
             */
            currentKineticEnergy = thermostat.getKineticEnergy();
            currentTemperature = thermostat.getCurrentTemperature();
            currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

            /**
             * Log the current state every printFrequency steps.
             */
            if (step % printFrequency == 0) {
                double simTime = step * dt;
                time = System.nanoTime() - time;
                logger.info(String.format(" %7.3e%13.4f%13.4f%13.4f%9.2f%9.3f", simTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * 1.0e-9));
                time = System.nanoTime();
            }

            /**
             * Write out snapshots in selected format every
             * saveSnapshotFrequency steps.
             */
            if (saveSnapshotFrequency > 0 && step % saveSnapshotFrequency == 0) {
                if (archiveFile != null && saveSnapshotAsPDB == false) {
                    if (xyzFilter.writeFile(archiveFile, true)) {
                        logger.info(String.format(" Appended snap shot to " + archiveFile.getName()));
                    } else {
                        logger.warning(String.format(" Appending snap shot to " + archiveFile.getName() + " failed"));
                    }
                } else if (saveSnapshotAsPDB == true) {
                    if (pdbFilter.writeFile(pdbFile, false)) {
                        logger.info(String.format(" Wrote PDB file to " + pdbFile.getName()));
                    }
                }
            }

            /**
             * Write out restart files every saveRestartFileFrequency steps.
             */
            if (saveRestartFileFrequency > 0 && step % saveRestartFileFrequency == 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(String.format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(String.format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }

            /**
             * Notify the algorithmListener.
             */
            if (algorithmListener != null && step % printFrequency == 0) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }

            /**
             * Check for a termination request.
             */
            if (terminate) {
                logger.info(String.format("\n Terminating after %8d time steps\n", step));
                break;
            }
        }

        /**
         * Log normal completion.
         */
        if (!terminate) {
            logger.info(String.format(" Completed %8d time steps\n", nSteps));
        }

        /**
         * Reset the done and terminate flags.
         */
        done = true;
        terminate = false;
    }

    /**
     * [SDL] Added as a way for pKa.groovy to request final delta_G
     */
    public double getTotalEnergy() {
        return currentTotalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating dynamics.\n", e);
                }
            }
        }
    }
}
