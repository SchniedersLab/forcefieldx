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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import static java.lang.String.format;
import static java.util.Arrays.fill;

import javax.swing.undo.CannotUndoException;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;


/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class MolecularDynamics implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    private final MolecularAssembly molecularAssembly;
    private final List<AssemblyInfo> assemblies;
    private final Potential potential;
    private AlgorithmListener algorithmListener;
    private MonteCarloListener monteCarloListener;
    private Thermostat thermostat;
    private Integrator integrator;
    private File restartFile = null;
    private DYNFilter dynFilter = null;
    private static final int DEFAULT_PRINT_FREQ = 250;
    private int printFrequency = DEFAULT_PRINT_FREQ;
    private static final int DEFAULT_SNAP_FREQ = 10000;
    private int saveSnapshotFrequency = DEFAULT_SNAP_FREQ;
    private int removeCOMMotionFrequency = 100;
    private boolean initVelocities = true;
    private boolean loadRestart = false;
    private boolean initialized = false;
    private boolean done = true;
    private boolean terminate = false;
    private int numberOfVariables;
    private double[] x;
    private double[] v;
    private double[] a;
    private double[] aPrevious;
    private double[] grad;
    private double[] mass;
    private int nSteps = 1000;
    private double targetTemperature = 298.15;
    private double dt = 1.0;
    private double currentTemperature;
    private double currentKineticEnergy;
    private double currentPotentialEnergy;
    private double currentTotalEnergy;
    private boolean saveSnapshotAsPDB = true;
    private int saveRestartFileFrequency = 1000;
    private String fileType = "XYZ";
    private double restartFrequency = 0.1;
    private boolean notifyMonteCarlo = true;
    private ExtendedSystem extendedSystem;
    private DynamicsState dynamicsState;
    private double totalSimTime = 0.0;

    /**
     * <p>
     * Constructor for MolecularDynamics.</p>
     *
     * @param assembly a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     * {@link ffx.algorithms.Thermostat.Thermostats} object.
     * @param requestedIntegrator a
     * {@link ffx.algorithms.Integrator.Integrators} object.
     */
    public MolecularDynamics(MolecularAssembly assembly,
            Potential potentialEnergy,
            CompositeConfiguration properties,
            AlgorithmListener listener,
            Thermostats requestedThermostat,
            Integrators requestedIntegrator) {
        this.molecularAssembly = assembly;
        assemblies = new ArrayList<>();
        assemblies.add(new AssemblyInfo(assembly));
        this.algorithmListener = listener;
        this.potential = potentialEnergy;

        assemblies.get(0).props = properties;
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
            thermostat.setRemoveCenterOfMassMotion(false);
        }

        done = true;
    }

    public MolecularDynamics(MolecularAssembly assembly,
            Potential potentialEnergy,
            CompositeConfiguration properties,
            AlgorithmListener listener,
            Thermostats requestedThermostat,
            Integrators requestedIntegrator,
            ExtendedSystem esvSystem) {
        this(assembly, potentialEnergy, properties, listener,
                requestedThermostat, requestedIntegrator);
        this.extendedSystem = esvSystem;
    }

    /**
     * Reinitialize the MD engine after a chemical change.
     */
    public void reInit() {
        mass = potential.getMass();
        numberOfVariables = potential.getNumberOfVariables();
        x = potential.getCoordinates(x);
        v = potential.getVelocity(v);
        a = potential.getAcceleration(a);
        aPrevious = potential.getPreviousAcceleration(aPrevious);
        if (potential instanceof ForceFieldEnergy) {
            grad = ((ForceFieldEnergy) potential).getGradients(grad);
        } else if (grad.length < numberOfVariables) {
            grad = new double[numberOfVariables];
        }
        thermostat.setNumberOfVariables(numberOfVariables, x, v, mass, potential.getVariableTypes());
        integrator.setNumberOfVariables(numberOfVariables, x, v, a, aPrevious, mass);
    }

    /**
     * <p>
     * Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link ffx.algorithms.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    /**
     * <p>
     * Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
    }

    /**
     * <p>
     * Setter for the field <code>x</code>.</p>
     *
     * @param x a double array to set the current parameters to.
     */
    public void setParameters(double x[]) {
        System.arraycopy(x, 0, this.x, 0, numberOfVariables);
    }

    /**
     * <p>
     * Getter for the field <code>x</code>.</p>
     *
     * @return a double array with the current parameters
     */
    public double[] getParameters() {
        return x;
    }

    /**
     * <p>
     * Setter for the first archive file.
     *
     * @param archive a {@link java.io.File} object.
     */
    public void setArchiveFile(File archive) {
        //this.archiveFile = archive;
        assemblies.get(0).archiveFile = archive;
    }

    public void setArchiveFile(File archive, int pos) {
        assemblies.get(pos).archiveFile = archive;
    }

    /**
     * <p>
     * Getter for the field <code>archiveFile</code>.</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getArchiveFile() {
        //return archiveFile;
        return assemblies.get(0).archiveFile;
    }

    public List<File> getArchiveFiles() {
        ArrayList<File> aFiles = new ArrayList<>();
        assemblies.forEach((ai) -> { aFiles.add(ai.archiveFile); });
        return aFiles;
        // Below may be more thread-safe and less prone to side effects.
        // Would definitely be safer than stream().forEach(add to external list).
        //return assemblies.stream().map((AssemblyInfo ai) -> {return ai.archiveFile;}).collect(Collectors.toList());
    }

    public void addAssembly(MolecularAssembly mola) {
        addAssembly(mola, null);
    }

    public void addAssembly(MolecularAssembly mola, CompositeConfiguration props) {
        AssemblyInfo mi = new AssemblyInfo(mola);
        mi.props = props;
        assemblies.add(mi);
    }

    /**
     * Finds and removes an assembly, searching by reference equality.
     * Removes all instances of the assembly.
     * @param mola Assembly to remove.
     * @return Number of times found and removed.
     */
    public int removeAssembly(MolecularAssembly mola) {
        if (mola == null) {
            return 0;
        }
        List<AssemblyInfo> toRemove = assemblies.stream().filter((AssemblyInfo ai) -> {
            return mola == ai.getAssembly();
        }).collect(Collectors.toList());
        assemblies.removeAll(toRemove);
        return toRemove.size();
    }

    public void setNotifyMonteCarlo(boolean set) {
        notifyMonteCarlo = set;
    }

    public void setMonteCarloListener(MonteCarloListener listener) {
        monteCarloListener = listener;
    }

    public void attachExtendedSystem(ExtendedSystem system) {
        if (system != null && extendedSystem != null) {
//            logger.warning("ExtendedSystem already attached to MD.");
        }
        extendedSystem = system;
    }

    public void detachExtendedSystem() {
        extendedSystem = null;
    }

    /**
     * <p>
     * init</p>
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param fileType a String.
     * @param restartFrequency the number of steps between writing restart
     * files.
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
        totalSimTime = 0.0;
        /**
         * Convert the time step from femtoseconds to picoseconds.
         */
        dt = timeStep * 1.0e-3;

        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = DEFAULT_PRINT_FREQ;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        /**
         * Convert the save interval to a save frequency.
         */
        saveSnapshotFrequency = DEFAULT_SNAP_FREQ;
        if (saveInterval >= this.dt) {
            saveSnapshotFrequency = (int) (saveInterval / this.dt);
        }

        /**
         * Set snapshot file type.
         */
        saveSnapshotAsPDB = true;
        if (fileType.equalsIgnoreCase("XYZ")) {
            saveSnapshotAsPDB = false;
        } else if (!fileType.equalsIgnoreCase("PDB")) {
            logger.warning("Snapshot file type unrecognized; saving snaphshots as PDB.\n");
        }

        /**
         * Convert restart frequency to steps.
         */
        saveRestartFileFrequency = 1000;
        if (restartFrequency >= this.dt) {
            saveRestartFileFrequency = (int) (restartFrequency / this.dt);
        }

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
            if (ainfo.xyzFilter == null) {
                ainfo.xyzFilter = new XYZFilter(file, mola, mola.getForceField(), aprops);
            }
            if (ainfo.pdbFilter == null) {
                ainfo.pdbFilter = new PDBFilter(ainfo.pdbFile, mola, mola.getForceField(), aprops);
            }
        });

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

        this.targetTemperature = temperature;
        this.initVelocities = initVelocities;
    }

    /**
     * A version of init with the original method header. Redirects to the new
     * method with default values for added parameters. Needed by (at least)
     * ReplicaExchange, which calls this directly.
     *
     * @param nSteps the number of MD steps.
     * @param timeStep the time step.
     * @param printInterval the number of steps between loggging updates.
     * @param saveInterval the number of steps between saving snapshots.
     * @param temperature the target temperature.
     * @param initVelocities true to reset velocities from a Maxwell
     * distribution.
     * @param dyn the Dynamic restart file.
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            final File dyn) {
        init(nSteps, timeStep, printInterval, saveInterval, "XYZ", 0.1, temperature, initVelocities, dyn);
    }

    private boolean skipIntro = false;
    public void redynamic(final int nSteps, final double temperature) {
        skipIntro = true;

        this.nSteps = nSteps;
        totalSimTime = 0.0;
//        this.dt = timeStep * 1.0e-3;
//        printFrequency = (int) (printInterval / this.dt);
//        saveSnapshotFrequency = (int) (saveInterval / this.dt);
//        saveSnapshotAsPDB = true;
//        if (fileType.equals("XYZ")) {
//            saveSnapshotAsPDB = false;
//        }
//        saveRestartFileFrequency = (int) (restartFrequency / this.dt);
//        if (pdbFilter == null) {
//            logger.warning("pdbf");
//        }
        this.targetTemperature = temperature;
        thermostat.setTargetTemperature(temperature);
        this.initVelocities = false;

        done = false;
        terminate = false;
        initialized = true;

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

        init(nSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency,
                temperature, initVelocities, dyn);

        done = false;

        if (dyn != null) {
            logger.info(format(" Continuing from " + dyn.getAbsolutePath()));
        }
        logger.info(String.format(" Number of steps:     %8d", nSteps));
        logger.info(String.format(" Time step:           %8.3f (fsec)", timeStep));
        logger.info(String.format(" Print interval:      %8.3f (psec)", printInterval));
        logger.info(String.format(" Save interval:       %8.3f (psec)", saveInterval));
        //logger.info(String.format(" Archive file: %s", archiveFile.getName()));
        for (int i = 0; i < assemblies.size(); i++) {
            AssemblyInfo ai = assemblies.get(i);
            logger.info(String.format(" Archive file %3d: %s", i, ai.archiveFile.getName()));
        }
        logger.info(String.format(" Restart file:     %s", restartFile.getName()));

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
        logger.info("Done with an MD round.");
    }

    /**
     * Method to set file type from groovy scripts.
     *
     * @param fileType the type of snapshot files to write.
     */
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    /**
     * Method to set the Restart Frequency.
     *
     * @param restartFrequency the time between writing restart files.
     */
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
            thermostat.setRemoveCenterOfMassMotion(true);
        } else {
            thermostat.setRemoveCenterOfMassMotion(false);
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
                Crystal crystal = molecularAssembly.getCrystal();
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    return;
                } else {
                    molecularAssembly.getPotentialEnergy().setCrystal(crystal);
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
                    fill(v, 0.0);
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

        if (!skipIntro) {
            logger.info(String.format("\n      Time      Kinetic    Potential        Total     Temp      CPU"));
            logger.info(String.format("      psec     kcal/mol     kcal/mol     kcal/mol        K      sec\n"));
            logger.info(String.format("          %13.4f%13.4f%13.4f %8.2f ", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));
        }

        /**
         * Integrate Newton's equations of motion for the requested number of
         * steps, unless early termination is requested.
         */
        long time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {
            if (notifyMonteCarlo && monteCarloListener != null) {
                long startTime = System.nanoTime();
                monteCarloListener.mcUpdate(molecularAssembly);
                x = potential.getCoordinates(x);
                long took = (long) ((System.nanoTime() - startTime) * 1e-6);
                // logger.info(String.format(" mcUpdate() took: %d ms", took));
            }

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
            if (thermostat.getRemoveCenterOfMassMotion() && step % removeCOMMotionFrequency == 0) {
                thermostat.centerOfMassMotion(true, false);
            }

            /**
             * Collect current kinetic energy, temperature, and total energy.
             */
            currentKineticEnergy = thermostat.getKineticEnergy();
            currentTemperature = thermostat.getCurrentTemperature();
            currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

            /**
             * Update atomic velocity, acceleration and previous acceleration.
             */
            potential.setVelocity(v);
            potential.setAcceleration(a);
            potential.setPreviousAcceleration(aPrevious);

            /**
             * Update extended system variables if present.
             */
            if (extendedSystem != null) {
                extendedSystem.propagateESVs(currentTemperature, dt, step*dt);
            }

            /**
             * Log the current state every printFrequency steps.
             */
            totalSimTime += dt;
            if (step % printFrequency == 0) {
                time = System.nanoTime() - time;
                logger.info(String.format(" %7.3e%13.4f%13.4f%13.4f%9.2f%9.3f", totalSimTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * 1.0e-9));
                time = System.nanoTime();
            }

            /**
             * Write out snapshots in selected format every
             * saveSnapshotFrequency steps.
             */
            if (saveSnapshotFrequency > 0 && step % saveSnapshotFrequency == 0) {
                /*if (archiveFile != null && saveSnapshotAsPDB == false) {
                    if (xyzFilter.writeFile(archiveFile, true)) {
                        logger.info(String.format(" Appended snap shot to " + archiveFile.getName()));
                    } else {
                        logger.warning(String.format(" Appending snap shot to " + archiveFile.getName() + " failed"));
                    }
                } else if (saveSnapshotAsPDB == true) {
                    if (pdbFilter.writeFile(pdbFile, false)) {
                        logger.info(String.format(" Wrote PDB file to " + pdbFile.getName()));
                    }
                }*/
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
                //algorithmListener.algorithmUpdate(molecularAssembly);
                for (AssemblyInfo assembly : assemblies) {
                    // Probably unwise to parallelize this, so that it doesn't
                    // hit the GUI with parallel updates.
                    algorithmListener.algorithmUpdate(assembly.getAssembly());
                }
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

        if (monteCarloListener != null) {
            long startTime = System.nanoTime();
            monteCarloListener.mcUpdate(molecularAssembly);
            x = potential.getCoordinates(x);
            long took = (long) ((System.nanoTime() - startTime) * 1e-6);
            // logger.info(String.format(" mcUpdate() took: %d ms", took));
        }
    }

    /**
     * Get the total system energy.
     *
     * @return total energy.
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

    /**
     * A simple container class to hold all the infrastructure associated with
     * a MolecularAssembly for MolecularDynamics; assembly, properties, archive
     * and PDB files, PDB and XYZ filters. Direct access to package-private
     * members breaks encapsulation a bit, but the private inner class shouldn't
     * be used externally anyways.
     */
    private class AssemblyInfo {
        private final MolecularAssembly mola;
        CompositeConfiguration props = null;
        File archiveFile = null;
        File pdbFile = null;
        PDBFilter pdbFilter = null;
        XYZFilter xyzFilter = null;

        public AssemblyInfo(MolecularAssembly assembly) {
            this.mola = assembly;
        }

        public MolecularAssembly getAssembly() {
            return mola;
        }
    }

    public void storeState() {
        dynamicsState = new DynamicsState();
    }

    public void revertState() {
        if (dynamicsState == null) {
            throw new CannotUndoException();
        }
        dynamicsState.restore();
    }

    private final boolean verboseDynamicsState = System.getProperty("md-verbose") != null;
    public class DynamicsState {
        double[] xBak, vBak, aBak;
        double[] aPreviousBak, massBak, gradBak;
        double currentKineticEnergyBak, currentPotentialEnergyBak, currentTotalEnergyBak;
        double currentTemperatureBak;
        public DynamicsState() {
            xBak = x.clone();
            vBak = v.clone();
            aBak = a.clone();
            aPreviousBak = aPrevious.clone();
            massBak = mass.clone();
            gradBak = grad.clone();
            currentKineticEnergyBak = currentKineticEnergy;
            currentPotentialEnergyBak = currentPotentialEnergy;
            currentTotalEnergyBak = currentTotalEnergy;
            currentTemperatureBak = currentTemperature;
            if (verboseDynamicsState) {
                describe(" Storing State:");
            }
        }
        public void describe(String title) {
                StringBuilder sb = new StringBuilder();
                sb.append(format(title));
                sb.append(format("\nx: "));
                Arrays.stream(x).forEach(val -> sb.append(format("%.2g, ", val)));
                sb.append(format("\nv: "));
                Arrays.stream(v).forEach(val -> sb.append(format("%.2g, ", val)));
                sb.append(format("\na: "));
                Arrays.stream(a).forEach(val -> sb.append(format("%.2g, ", val)));
                sb.append(format("\naP: "));
                Arrays.stream(aPrevious).forEach(val -> sb.append(format("%.2g, ", val)));
                sb.append(format("\nm: "));
                Arrays.stream(mass).forEach(val -> sb.append(format("%.2g, ", val)));
                sb.append(format("\ng: "));
                Arrays.stream(grad).forEach(val -> sb.append(format("%.2g, ", val)));
                sb.append(format("\nK,U,E,T: %g %g %g %g\n",
                        currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature));
                logger.info(sb.toString());
        }
        public void restore() {
            if (verboseDynamicsState) {
                describe(" Reverting State (From):");
            }
            x = xBak;
            v = vBak;
            a = aBak;
            aPrevious = aPreviousBak;
            mass = massBak;
            grad = gradBak;
            currentKineticEnergy = currentKineticEnergyBak;
            currentPotentialEnergy = currentPotentialEnergyBak;
            currentTotalEnergy = currentTotalEnergyBak;
            currentTemperature = currentTemperatureBak;
            if (verboseDynamicsState) {
                describe(" Reverting State (To):");
            }
        }
    }
}
