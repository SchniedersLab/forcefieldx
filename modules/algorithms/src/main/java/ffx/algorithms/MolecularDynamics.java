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
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import org.apache.commons.collections.queue.CircularFifoQueue;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.integrators.BetterBeeman;
import ffx.algorithms.integrators.Integrator;
import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.integrators.Respa;
import ffx.algorithms.integrators.Stochastic;
import ffx.algorithms.integrators.VelocityVerlet;
import ffx.algorithms.thermostats.Adiabatic;
import ffx.algorithms.thermostats.Berendsen;
import ffx.algorithms.thermostats.Bussi;
import ffx.algorithms.thermostats.Thermostat;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;

/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MolecularDynamics implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());

    private final Potential potential;

    private final DynamicsEngine engine = DynamicsEngine.FFX;

    private MonteCarloListener monteCarloListener;
    private Thermostat thermostat;
    private Integrator integrator;

    private int printEsvFrequency = -1;
    private int removeCOMMotionFrequency = 100;
    private boolean constantPressure = false;
    private boolean terminate = false;
    private int nSteps = 1000;

    private ExtendedSystem esvSystem;
    private DynamicsState dynamicsState;
    private double totalSimTime = 0.0;
    private long time;
    private long mdTime = 0;
    private boolean quiet = false;

    protected final MolecularAssembly molecularAssembly;
    protected final List<AssemblyInfo> assemblies;
    protected File restartFile = null;
    protected boolean loadRestart = false;
    protected boolean initialized = false;
    protected boolean done = true;
    protected boolean initVelocities = true;
    protected AlgorithmListener algorithmListener;
    protected DYNFilter dynFilter = null;
    protected int numberOfVariables;
    protected double[] x;
    protected double[] v;
    protected double[] a;
    protected double[] aPrevious;
    protected double[] grad;
    protected double[] mass;
    protected double dt = 1.0;
    protected int printFrequency = 100;
    protected double restartFrequency = 0.1;
    protected int saveRestartFileFrequency = 1000;
    protected int saveSnapshotFrequency = 1000;
    protected double targetTemperature = 298.15;
    protected double currentTemperature;
    protected double currentKineticEnergy;
    protected double currentPotentialEnergy;
    protected double currentTotalEnergy;
    protected boolean saveSnapshotAsPDB = true;
    protected String fileType = "XYZ";
    protected static final double NS2SEC = 1e-9;

    /**
     * Keep some old coordinate snapshots around.
     */
    private int numSnapshotsToKeep = 0;
    // Circular FIFO queues will simply discard old elements.
    private CircularFifoQueue<DynamicsState> lastSnapshots;
    // A change in potential energy exceeding 1E6 kcal/mol triggers a warning and snapshot dump.
    private double defaultDeltaPEThresh = 1.0E6;

    private MonteCarloNotification mcNotification = MonteCarloNotification.NEVER;

    public enum MonteCarloNotification {
        NEVER, EACH_STEP, AFTER_DYNAMICS;
    }

    /**
     * Method that determines whether a dynamics is done by the java
     * implementation native to ffx or the OpenMM implementation
     *
     * @param assembly
     * @param potentialEnergy
     * @param properties
     * @param listener
     * @param requestedThermostat
     * @param requestedIntegrator
     * @return
     */
    public static MolecularDynamics dynamicsFactory(MolecularAssembly assembly,
                                                    Potential potentialEnergy,
                                                    CompositeConfiguration properties,
                                                    AlgorithmListener listener,
                                                    ThermostatEnum requestedThermostat,
                                                    IntegratorEnum requestedIntegrator) {

        return dynamicsFactory(assembly, potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator, defaultEngine(potentialEnergy));
        /*
        if (potentialEnergy instanceof ForceFieldEnergyOpenMM) {
            MolecularDynamicsOpenMM ommDynamics = new MolecularDynamicsOpenMM(assembly,
                    (ForceFieldEnergyOpenMM) potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);
            return ommDynamics;
        } else {
            MolecularDynamics mDynamics = new MolecularDynamics(assembly,
                    potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);
            return mDynamics;
        }
        */
    }

    public static MolecularDynamics dynamicsFactory(MolecularAssembly assembly,
                                                    Potential potentialEnergy,
                                                    CompositeConfiguration properties,
                                                    AlgorithmListener listener,
                                                    ThermostatEnum requestedThermostat,
                                                    IntegratorEnum requestedIntegrator,
                                                    DynamicsEngine engine) {
        switch (engine) {
            case OPENMM:
                return new MolecularDynamicsOpenMM(assembly,
                        (ForceFieldEnergyOpenMM) potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);
            case FFX:
            default:
                return new MolecularDynamics(assembly,
                        potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);
        }
    }

    private static DynamicsEngine defaultEngine(Potential potentialEnergy) {
        if (System.getProperty("MD-engine") != null) {
            if (System.getProperty("MD-engine").equalsIgnoreCase("OMM")) {
                logger.info(String.format(" Creating OpenMM Dynamics Object"));
                return DynamicsEngine.OPENMM;
            } else {
                logger.info(String.format(" Creating FFX Dynamics Object"));
                return DynamicsEngine.FFX;
            }
        } else {
            if (potentialEnergy instanceof ffx.potential.ForceFieldEnergyOpenMM) {
                return DynamicsEngine.OPENMM;
            } else {
                return DynamicsEngine.FFX;
            }
        }
    }

    /**
     * <p>
     * Constructor for MolecularDynamics.</p>
     *
     * @param assembly            a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy     a {@link ffx.numerics.Potential} object.
     * @param properties          a
     *                            {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     *                            {@link ffx.algorithms.thermostats.ThermostatEnum} object.
     * @param requestedIntegrator a
     *                            {@link ffx.algorithms.integrators.IntegratorEnum} object.
     */
    public MolecularDynamics(MolecularAssembly assembly,
                             Potential potentialEnergy,
                             CompositeConfiguration properties,
                             AlgorithmListener listener,
                             ThermostatEnum requestedThermostat,
                             IntegratorEnum requestedIntegrator) {
        this.molecularAssembly = assembly;
        assemblies = new ArrayList<>();
        assemblies.add(new AssemblyInfo(assembly));
        this.algorithmListener = listener;
        this.potential = potentialEnergy;
        if (potentialEnergy instanceof Barostat) {
            constantPressure = true;
        }

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
                requestedIntegrator = IntegratorEnum.valueOf(integrate);
            } catch (Exception e) {
                requestedIntegrator = IntegratorEnum.BEEMAN;
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
                requestedThermostat = ThermostatEnum.ADIABATIC;
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
                requestedThermostat = ThermostatEnum.valueOf(thermo);
            } catch (Exception e) {
                requestedThermostat = ThermostatEnum.BERENDSEN;
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
         * For Stochastic dynamics, center of mass motion will not be removed.
         */
        if (integrator instanceof Stochastic) {
            thermostat.setRemoveCenterOfMassMotion(false);
        }

        numSnapshotsToKeep = properties.getInteger("dynamicsSnapshotMemory", 0);
        // Cannot construct a CircularFifoQueue of zero length.
        lastSnapshots = new CircularFifoQueue<>(Math.max(numSnapshotsToKeep, 1));

        done = true;
    }

    public MolecularDynamics(MolecularAssembly assembly,
                             Potential potentialEnergy,
                             CompositeConfiguration properties,
                             AlgorithmListener listener,
                             ThermostatEnum requestedThermostat,
                             IntegratorEnum requestedIntegrator,
                             ExtendedSystem esvSystem) {
        this(assembly, potentialEnergy, properties, listener,
                requestedThermostat, requestedIntegrator);
        this.esvSystem = esvSystem;
    }

    public void setQuiet(boolean quiet) {
        this.quiet = quiet;
    }

    /**
     * Reinitialize the MD engine after a chemical change.
     */
    public void reInit() {
        numberOfVariables = potential.getNumberOfVariables();
        mass = potential.getMass();
        x = new double[numberOfVariables];
        v = new double[numberOfVariables];
        a = new double[numberOfVariables];
        aPrevious = new double[numberOfVariables];
        grad = new double[numberOfVariables];
        potential.getCoordinates(x);
        potential.getVelocity(v);
        potential.getAcceleration(a);
        potential.getPreviousAcceleration(aPrevious);
        if (potential instanceof ForceFieldEnergy) {
            grad = ((ForceFieldEnergy) potential).getGradients(grad);
        }
        thermostat.setNumberOfVariables(numberOfVariables, x, v, mass,
                potential.getVariableTypes(), true);
        integrator.setNumberOfVariables(numberOfVariables, x, v, a, aPrevious, mass);
    }

    /**
     * <p>
     * Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    /**
     * <p>
     * Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
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

    /**
     * Setter for an archive file of arbitrary position.
     *
     * @param archive A File to set as archive
     * @param pos     Index of MolecularAssembly to set this for
     */
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

    /**
     * Gets a list of all archive files.
     *
     * @return A List of archive files
     */
    public List<File> getArchiveFiles() {
        // This implementation seems to work, but I'm pretty sure it's poor practice.
        ArrayList<File> aFiles = new ArrayList<>();
        assemblies.forEach((ai) -> {
            aFiles.add(ai.archiveFile);
        });
        return aFiles;
        // Below may be more thread-safe and less prone to side effects.
        // Would definitely be safer than stream().forEach(add to external list).
        // return assemblies.stream().map((AssemblyInfo ai) -> {return ai.archiveFile;}).collect(Collectors.toList());
    }

    /**
     * Adds a MolecularAssembly to be tracked by this MolecularDynamics. Note:
     * does not affect the underlying Potential.
     *
     * @param mola A MolecularAssembly to be tracked
     */
    public void addAssembly(MolecularAssembly mola) {
        addAssembly(mola, mola.getProperties());
    }

    /**
     * Adds a MolecularAssembly to be tracked by this MolecularDynamics. Note:
     * does not affect the underlying Potential.
     *
     * @param mola  A MolecularAssembly to be tracked
     * @param props Associated CompositeConfiguration
     */
    public void addAssembly(MolecularAssembly mola, CompositeConfiguration props) {
        AssemblyInfo mi = new AssemblyInfo(mola);
        mi.props = props;
        assemblies.add(mi);
    }

    /**
     * Finds and removes an assembly, searching by reference equality. Removes
     * all instances of the assembly. Note: does not affect the underlying
     * Potential.
     *
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

    public void setMonteCarloListener(MonteCarloListener listener, MonteCarloNotification when) {
        monteCarloListener = listener;
        mcNotification = when;
    }

    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        if (esvSystem != null) {
            logger.warning("An ExtendedSystem is already attached to this MD!");
        }
        esvSystem = system;
        printEsvFrequency = printFrequency;
        logger.info(format(" Attached extended system (%s) to molecular dynamics.", esvSystem.toString()));
        reInit();
    }

    public void detachExtendedSystem() {
        logger.info(format(" Detached extended system (%s) from molecular dynamics.", esvSystem.toString()));
        esvSystem = null;
        reInit();
    }

    /**
     * <p>
     * init</p>
     *
     * @param nSteps           a int.
     * @param timeStep         a double.
     * @param printInterval    a double.
     * @param saveInterval     a double.
     * @param fileType         a String.
     * @param restartFrequency the number of steps between writing restart
     *                         files.
     * @param temperature      a double.
     * @param initVelocities   a boolean.
     * @param dyn              a {@link java.io.File} object.
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

        assemblyInfo();

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

    protected void assemblyInfo() {
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
    }

    /**
     * A version of init with the original method header. Redirects to the new
     * method with default values for added parameters. Needed by (at least)
     * ReplicaExchange, which calls this directly.
     *
     * @param nSteps         the number of MD steps.
     * @param timeStep       the time step.
     * @param printInterval  the number of steps between loggging updates.
     * @param saveInterval   the number of steps between saving snapshots.
     * @param temperature    the target temperature.
     * @param initVelocities true to reset velocities from a Maxwell
     *                       distribution.
     * @param dyn            the Dynamic restart file.
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
                     final double saveInterval, final double temperature, final boolean initVelocities,
                     final File dyn) {
        init(nSteps, timeStep, printInterval, saveInterval, "XYZ", 0.1, temperature, initVelocities, dyn);
    }

    /**
     * Causes this MolecularDynamics to take an additional set of timesteps.
     *
     * @param nSteps      Number of steps to take
     * @param temperature Temperature of simulation
     */
    public void redynamic(final int nSteps, final double temperature) {
        quiet = true;

        this.nSteps = nSteps;
        totalSimTime = 0.0;
        targetTemperature = temperature;
        thermostat.setTargetTemperature(temperature);
        initVelocities = false;

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
     * @param nSteps           a int.
     * @param timeStep         a double.
     * @param printInterval    a double.
     * @param saveInterval     a double.
     * @param temperature      a double.
     * @param initVelocities   a boolean.
     * @param fileType         a String (either XYZ or PDB).
     * @param restartFrequency a double specifying the restart frequency.
     * @param dyn              a {@link java.io.File} object.
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
     * @param nSteps         a int.
     * @param timeStep       a double.
     * @param printInterval  a double.
     * @param saveInterval   a double.
     * @param temperature    a double.
     * @param initVelocities a boolean.
     * @param dyn            a {@link java.io.File} object.
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
            if (constantPressure) {
                logger.info(format("\n Stochastic dynamics in the NPT ensemble"));
            } else {
                logger.info(format("\n Stochastic dynamics in the NVT ensemble"));
            }
        } else if (!(thermostat instanceof Adiabatic)) {
            if (constantPressure) {
                logger.info(format("\n Molecular dynamics in the NPT ensemble"));
            } else {
                logger.info(format("\n Molecular dynamics in the NVT ensemble"));
            }
        } else {
            if (constantPressure) {
                logger.severe(format("\n NPT Molecular dynamics requires a thermostat"));
            } else {
                logger.info(format("\n Molecular dynamics in the NVE ensemble"));
            }
        }

        init(nSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency,
                temperature, initVelocities, dyn);

        done = false;

        if (dyn != null) {
            logger.info(format(" Continuing from " + dyn.getAbsolutePath()));
        }

        if (!quiet) {
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
        }

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
        if (!quiet) {
            logger.info(" Done with an MD round.");
        }
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
     * @throws IllegalArgumentException If restart frequency is not a positive
     *                                  number
     */
    public void setRestartFrequency(double restartFrequency) throws IllegalArgumentException {
        if (Double.isFinite(restartFrequency) && restartFrequency > 0) {
            this.restartFrequency = restartFrequency;
        } else {
            throw new IllegalArgumentException(String.format(" Restart frequency must be positive finite, was %10.4g", restartFrequency));
        }
    }

    /**
     * Set the number of time steps between removal of center of mass kinetic
     * energy.
     *
     * @param removeCOMMotionFrequency Number of time steps between center of
     *                                 mass removal.
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
        thermostat.setQuiet(quiet);
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
        try {
            currentPotentialEnergy = potential.energyAndGradient(x, grad);
        } catch (EnergyException ex) {
            writeStoredSnapshots();
            throw ex;
        }

        /**
         * Initialize current and previous accelerations, unless they were just
         * loaded from a restart file.
         */
        if (!loadRestart || initialized) {
            for (int i = 0; i < numberOfVariables; i++) {
                a[i] = -Thermostat.convert * grad[i] / mass[i];
            }
            if (aPrevious != null) {
                arraycopy(a, 0, aPrevious, 0, numberOfVariables);
            }
        }

        /**
         * Compute the current kinetic energy.
         */
        thermostat.kineticEnergy();
        currentKineticEnergy = thermostat.getKineticEnergy();
        currentTemperature = thermostat.getCurrentTemperature();
        currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

        initialized = true;
        logger.info(format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
        logger.info(format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
        logger.info(format("  %8s %12.4f %12.4f %12.4f %8.2f",
                "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));

        /**
         * Store the initialized state.
         */
        storeState();

        /**
         * Integrate Newton's equations of motion for the requested number of
         * steps, unless early termination is requested.
         */
        time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {
            /* Notify MonteCarlo handlers such as PhMD or rotamer drivers. */
            if (monteCarloListener != null && mcNotification == MonteCarloNotification.EACH_STEP) {
                monteCarloListener.mcUpdate(thermostat.getCurrentTemperature());
            }

            /**
             * Do the half-step thermostat operation.
             */
            thermostat.halfStep(dt);

            /**
             * Do the half-step integration operation.
             */
            integrator.preForce(potential);

            double priorPE = currentPotentialEnergy;
            /**
             * Compute the potential energy and gradients.
             */
            try {
                currentPotentialEnergy = potential.energyAndGradient(x, grad);
            } catch (EnergyException ex) {
                writeStoredSnapshots();
                throw ex;
            }

            detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

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
            if (esvSystem != null) {
                esvSystem.propagateESVs(currentTemperature, dt, step * dt);
            }

            /**
             * Log the current state every printFrequency steps.
             */
            totalSimTime += dt;
            if (step % printFrequency == 0) {

                /**
                 * Original print statement
                 */
                time = System.nanoTime() - time;
                mdTime = time;
                logger.info(format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.3f",
                        totalSimTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * NS2SEC));
                time = System.nanoTime();

                // Shirts et al. logging info
                Crystal crystal = molecularAssembly.getCrystal();
                double volume = crystal.getUnitCell().volume;

                // Shirts et al. print statement
//                time = System.nanoTime() - time;
//                logger.info(format("Shirts %7.3e %12.4f %12.4f %12.4f %12.4f %8.2f %8.3f",
//                totalSimTime, currentKineticEnergy, currentPotentialEnergy,
//                currentTotalEnergy, volume, currentTemperature, time * NS2SEC));
//                time = System.nanoTime();

            }
            if (step % printEsvFrequency == 0 && esvSystem != null) {
                logger.info(format(" %7.3e %s", totalSimTime, esvSystem.getLambdaList()));
            }

            /**
             * Write out snapshots in selected format every
             * saveSnapshotFrequency steps.
             */
            if (saveSnapshotFrequency > 0 && step % saveSnapshotFrequency == 0) {
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

        if (monteCarloListener != null && mcNotification == MonteCarloNotification.AFTER_DYNAMICS) {
            monteCarloListener.mcUpdate(thermostat.getCurrentTemperature());
        }
    }

    /**
     * Detects grossly atypical potential energy values that are likely incorrect,
     * and writes snapshots to disc. "Grossly atypical" is defined as: greater
     * than 1.0E100 kcal/mol, less than -1.0E100 kcal/mol, non-finite (NaN/infinite),
     * or exceeding specified delta from the prior potential energy.
     *
     * After prior snapshots have been written to disc, the queue they are stored on
     * is now empty, preventing printing of duplicate snapshots.
     *
     * @param priorPE Potential energy prior to this step.
     * @param delPEThresh If potential energy changes by this much, trigger a snapshot write.
     * @return True if atypical energy detected.
     */
    protected boolean detectAtypicalEnergy(double priorPE, double delPEThresh) {

        // If not keeping snapshots, disable functionality.
        if (numSnapshotsToKeep < 1) {
            return false;
        }

        double deltaPE = currentPotentialEnergy - priorPE;

        DynamicsState currState = new DynamicsState();
        currState.storeState();
        lastSnapshots.add(currState);

        double maxPEThresh = 1.0E100; // 1.0E100 kcal/mol is well into the territory of the absurd.
        double absPE = Math.abs(currentPotentialEnergy);

        if (absPE > maxPEThresh || !Double.isFinite(currentPotentialEnergy) || Math.abs(deltaPE) > delPEThresh) {
            logger.info(String.format(" Unusual potential energy %12.5g detected, writing snapshots.", currentPotentialEnergy));
            writeStoredSnapshots();
            currState.revertState(); // May be unnecessary, thanks to the current state always being last on the queue.
            if (absPE > 1.0E100 || !Double.isFinite(currentPotentialEnergy)) {
                logger.severe(String.format(" Dynamics exiting with atypical potential energy of %12.5g", currentPotentialEnergy));
            }
            return true;
        }
        return false;
    }

    /**
     * Performs the inner loop of writing snapshots to disk; used by both detectAtypicalEnergy and a try-catch in dynamics.
     */
    private void writeStoredSnapshots() {
        int numSnaps = lastSnapshots.size();

        File origFile = molecularAssembly.getFile();
        String timeString = LocalDateTime.now().format(DateTimeFormatter.
                ofPattern("HH_mm_ss"));
        PotentialsFunctions ef = new PotentialsUtils();

        String filename = String.format("%s-%s-SNAP.pdb",
                FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                timeString);

        for (int is = 0; is < numSnaps; is++) {
            DynamicsState oldState = lastSnapshots.poll();
            oldState.revertState();

            ef.saveAsPDB(molecularAssembly, new File(ef.versionFile(filename)));
        }
        molecularAssembly.setFile(origFile);
    }

    /**
     * Get the total system energy (kinetic plus potential).
     *
     * @return total energy.
     */
    public double getTotalEnergy() {
        return currentTotalEnergy;
    }

    /**
     * Get the system kinetic energy.
     *
     * @return kinetic energy.
     */
    public double getKineticEnergy() {
        return currentKineticEnergy;
    }

    /**
     * Get the system potential energy.
     *
     * @return potential energy.
     */
    public double getPotentialEnergy() {
        return currentPotentialEnergy;
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
     * Returns the DynamicsEngine in use; will typically be over-ridden by
     * subclasses.
     *
     * @return FFX engine.
     */
    public DynamicsEngine getEngine() {
        return engine;
    }

    /**
     * A simple container class to hold all the infrastructure associated with a
     * MolecularAssembly for MolecularDynamics; assembly, properties, archive
     * and PDB files, PDB and XYZ filters. Direct access to package-private
     * members breaks encapsulation a bit, but the private inner class shouldn't
     * be used externally anyways.
     */
    protected class AssemblyInfo {

        private final MolecularAssembly mola;
        CompositeConfiguration props = null;
        File archiveFile = null;
        File pdbFile = null;
        PDBFilter pdbFilter = null;
        XYZFilter xyzFilter = null;

        public AssemblyInfo(MolecularAssembly assembly) {
            mola = assembly;
            pdbFile = mola.getFile();
            props = mola.getProperties();
            pdbFilter = new PDBFilter(mola.getFile(), mola, mola.getForceField(), mola.getProperties());
        }

        public MolecularAssembly getAssembly() {
            return mola;
        }
    }

    public void storeState() {
        if (dynamicsState == null) {
            dynamicsState = new DynamicsState();
        }
        dynamicsState.storeState();
    }

    public void revertState() throws Exception {
        if (dynamicsState == null) {
            throw new Exception();
        }
        dynamicsState.revertState();
    }

    private final boolean verboseDynamicsState = System.getProperty("md-verbose") != null;

    protected class DynamicsState {

        double[] xBak, vBak, aBak;
        double[] aPreviousBak, massBak, gradBak;
        double currentKineticEnergyBak, currentPotentialEnergyBak, currentTotalEnergyBak;
        double currentTemperatureBak;

        public DynamicsState() {
            xBak = new double[numberOfVariables];
            vBak = new double[numberOfVariables];
            aBak = new double[numberOfVariables];
            aPreviousBak = new double[numberOfVariables];
            massBak = new double[numberOfVariables];
            gradBak = new double[numberOfVariables];
        }

        public void storeState() {
            currentKineticEnergyBak = currentKineticEnergy;
            currentPotentialEnergyBak = currentPotentialEnergy;
            currentTotalEnergyBak = currentTotalEnergy;
            currentTemperatureBak = currentTemperature;
            arraycopy(x, 0, xBak, 0, numberOfVariables);
            arraycopy(v, 0, vBak, 0, numberOfVariables);
            arraycopy(a, 0, aBak, 0, numberOfVariables);
            arraycopy(aPrevious, 0, aPreviousBak, 0, numberOfVariables);
            arraycopy(mass, 0, massBak, 0, numberOfVariables);
            arraycopy(grad, 0, gradBak, 0, numberOfVariables);
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

        public void revertState() {
            if (verboseDynamicsState) {
                describe(" Reverting State (From):");
            }
            currentKineticEnergy = currentKineticEnergyBak;
            currentPotentialEnergy = currentPotentialEnergyBak;
            currentTotalEnergy = currentTotalEnergyBak;
            currentTemperature = currentTemperatureBak;
            arraycopy(xBak, 0, x, 0, numberOfVariables);
            arraycopy(vBak, 0, v, 0, numberOfVariables);
            arraycopy(aBak, 0, a, 0, numberOfVariables);
            arraycopy(aPreviousBak, 0, aPrevious, 0, numberOfVariables);
            arraycopy(massBak, 0, mass, 0, numberOfVariables);
            arraycopy(gradBak, 0, grad, 0, numberOfVariables);

            Atom atoms[] = molecularAssembly.getActiveAtomArray();
            if (atoms.length * 3 == numberOfVariables) {
                int nAtoms = atoms.length;
                int index = 0;
                double vec[] = new double[3];
                for (int i = 0; i < nAtoms; i++) {
                    Atom atom = atoms[i];
                    atom.moveTo(x[index], x[index + 1], x[index + 2]);
                    atom.setXYZGradient(grad[index], grad[index + 1], grad[index + 2]);
                    vec[0] = v[index];
                    vec[1] = v[index + 1];
                    vec[2] = v[index + 2];
                    atom.setVelocity(vec);
                    vec[0] = a[index];
                    vec[1] = a[index + 1];
                    vec[2] = a[index + 2];
                    atom.setAcceleration(vec);
                    vec[0] = aPrevious[index];
                    vec[1] = aPrevious[index + 1];
                    vec[2] = aPrevious[index + 2];
                    atom.setPreviousAcceleration(vec);
                    index += 3;
                    if (i == 0) {
                        //logger.info(String.format(" In revert for loop"));
                    }
                }
                //logger.info(String.format(" Coordinates were reverted"));
            } else {
                //logger.info(String.format( "WARNING: Coordinates were not reverted"));
            }
            if (verboseDynamicsState) {
                describe(" Reverting State (To):");
            }
        }
    }

    /**
     * Enumerates available molecular dynamics engines; presently limited to the
     * FFX reference engine and the OpenMM engine. Distinct from the force field
     * energy Platform, as the FFX engine can use OpenMM energies, but not
     * vice-versa.
     */
    public enum DynamicsEngine {
        FFX(true, true), OPENMM(false, true);

        // Set of supported Platforms. The EnumSet paradigm is very efficient, as it
        // is internally stored as a bit field.
        private final EnumSet<ForceFieldEnergy.Platform> platforms = EnumSet.noneOf(ForceFieldEnergy.Platform.class);

        /**
         * Constructs a DynamicsEngine using the two presently known types of
         * Platform.
         *
         * @param ffx    Add support for the FFX reference energy platform.
         * @param openMM Add support for the OpenMM energy platforms.
         */
        DynamicsEngine(boolean ffx, boolean openMM) {
            if (ffx) {
                platforms.add(ForceFieldEnergy.Platform.FFX);
            }
            if (openMM) {
                platforms.add(ForceFieldEnergy.Platform.OMM);
                platforms.add(ForceFieldEnergy.Platform.OMM_REF);
                platforms.add(ForceFieldEnergy.Platform.OMM_CUDA);
                platforms.add(ForceFieldEnergy.Platform.OMM_OPENCL);
                platforms.add(ForceFieldEnergy.Platform.OMM_OPTCPU);
            }
        }

        /**
         * Checks if this energy Platform is supported by this DynamicsEngine
         *
         * @param platform
         * @return If supported
         */
        public boolean supportsPlatform(ForceFieldEnergy.Platform platform) {
            return platforms.contains(platform);
        }

        /**
         * Gets the set of Platforms supported by this DynamicsEngine
         *
         * @return An EnumSet
         */
        public EnumSet<ForceFieldEnergy.Platform> getSupportedPlatforms() {
            return EnumSet.copyOf(platforms);
        }
    }
    
    public long getMDTime(){
        return mdTime;
    }
    
    public double getStartingTotalEnergy(){
        return 0.0;
    }
    
    public double getEndTotalEnergy(){
        return 0.0;
    }
    
    public double getTimeStep(){
        return dt;
    }
    
    public int getIntervalSteps(){
        return 1;
    }
    
    public int getNumAtoms(){
        return 1;
    }
    
}
