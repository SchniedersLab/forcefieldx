//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.dynamics;

import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import ffx.potential.ForceFieldEnergyOpenMM;
import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.algorithms.dynamics.integrators.BetterBeeman;
import ffx.algorithms.dynamics.integrators.Integrator;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.integrators.Respa;
import ffx.algorithms.dynamics.integrators.Stochastic;
import ffx.algorithms.dynamics.integrators.VelocityVerlet;
import ffx.algorithms.dynamics.thermostats.Adiabatic;
import ffx.algorithms.dynamics.thermostats.Berendsen;
import ffx.algorithms.dynamics.thermostats.Bussi;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.mc.MonteCarloListener;
import ffx.crystal.Crystal;
import ffx.numerics.Constraint;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.NS2SEC;

/**
 * Run NVE, NVT, or NPT molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MolecularDynamics implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    /**
     * MolecularAssembly to run dynamics on.
     */
    protected final MolecularAssembly molecularAssembly;
    /**
     * List of MolecularAssembly instances.
     */
    protected final List<AssemblyInfo> assemblies;
    /**
     * Propogate dynamics on this potential surface.
     */
    private final Potential potential;
    /**
     * Monte Carlo listener.
     */
    private MonteCarloListener monteCarloListener;
    /**
     * An Algorithm Listener to send updates to the GUI.
     */
    protected AlgorithmListener algorithmListener;
    /**
     * Thermostat instance.
     */
    private Thermostat thermostat;
    /**
     * Integrator instance.
     */
    private Integrator integrator;
    /**
     * Flag to indicate use of constant pressure.
     */
    boolean constantPressure = false;
    /**
     * Any Barostat that may be in use.
     */
    Barostat barostat = null;
    /**
     * Flag to indicate MD should be terminated.
     */
    private boolean terminate = false;
    /**
     * Number of MD steps to take.
     */
    private int nSteps = 1000;
    /**
     * State of the dynamics.
     */
    private DynamicsState dynamicsState;
    /**
     * Total simulation time.
     */
    private double totalSimTime = 0.0;
    /**
     * Indicates how verbose MD should be.
     */
    private VerbosityLevel verbosityLevel = VerbosityLevel.VERBOSE;
    /**
     * Log some of the more frequent messages at this level. Always at or below basicLogging.
     */
    Level intermediateLogging = Level.INFO;
    /**
     * Log basic information at this level. Always at or above intermediateLogging.
     */
    Level basicLogging = Level.INFO;
    /**
     * Dynamics restart file.
     */
    File restartFile = null;
    /**
     * Flag to indicate loading of restart file.
     */
    boolean loadRestart = false;
    /**
     * Filter to parse the dynamics restart file.
     */
    DYNFilter dynFilter = null;
    /**
     * Flag to indicate dynamics has been initialized.
     */
    boolean initialized = false;
    /**
     * Flag to indicate a run has finished.
     */
    protected boolean done;
    /**
     * Flag to indicate velocities should be initialized.
     */
    protected boolean initVelocities = true;
    /**
     * Number of dynamics variables.
     */
    int numberOfVariables;
    /**
     * Coordinates.
     */
    protected double[] x;
    /**
     * Velocities.
     */
    protected double[] v;
    /**
     * Accelerations.
     */
    protected double[] a;
    /**
     * Previous accelerations.
     */
    double[] aPrevious;
    /**
     * The gradient.
     */
    protected double[] gradient;
    /**
     * Mass for each degree of freedom.
     */
    protected double[] mass;
    /**
     * Time step (picoseconds).
     */
    protected double dt = 1.0;
    /**
     * Number of steps between logging info to the screen.
     */
    int printFrequency = 100;
    /**
     * Frequency to write out restart info.
     */
    protected double restartFrequency = 0.1;
    /**
     * Frequency to save restart files.
     */
    int saveRestartFileFrequency = 1000;
    /**
     * Frequency to save snap shot files.
     */
    int saveSnapshotFrequency = 1000;
    /**
     * Target temperature. ToDo: use the Thermostat instance.
     */
    double targetTemperature = 298.15;
    /**
     * Current temperature.
     */
    double currentTemperature;
    /**
     * Current kinetic energy.
     */
    double currentKineticEnergy;
    /**
     * Current potential energy.
     */
    double currentPotentialEnergy;
    /**
     * Current total energy.
     */
    double currentTotalEnergy;
    /**
     * Kinetic energy before taking a time step.
     */
    double startingKineticEnergy;
    /**
     * Potential energy before taking a time step.
     */
    double startingPotentialEnergy;
    /**
     * Total energy before taking a time step.
     */
    double startingTotalEnergy;
    /**
     * Save snapshots in PDB format.
     */
    boolean saveSnapshotAsPDB = true;
    /**
     * File type to use when saving files.
     */
    protected String fileType = "XYZ";
    /**
     * Keep some old coordinate snapshots around.
     */
    private int numSnapshotsToKeep;
    /**
     * Circular FIFO queues will simply discard old elements.
     */
    private CircularFifoQueue<CoordinateSnapshot> lastSnapshots;
    /**
     * MC notification flag.
     */
    private MonteCarloNotification mcNotification = MonteCarloNotification.NEVER;
    /**
     * Monte Carlo notification enumeration.
     */
    public enum MonteCarloNotification {
        NEVER, EACH_STEP, AFTER_DYNAMICS
    }
    /**
     * ESV System.
     */
    private ExtendedSystem esvSystem;
    /**
     * Frequency to print ESV info.
     */
    private int printEsvFrequency = -1;
    /**
     * Verbose logging of dynamics state.
     */
    private final boolean verboseDynamicsState;

    /**
     * By default, wait 25000 ns (25 us) in between polling the dynamics thread.
     */
    private static final int DEFAULT_DYNAMICS_SLEEP_TIME = 25000;

    /**
     * Wait this many nanoseconds in between polling the dynamics thread.
     */
    private final int dynSleepTime;

    /**
     * Method that determines whether a dynamics is done by the java
     * implementation native to ffx or the OpenMM implementation
     *
     * @param assembly            a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy     a {@link ffx.numerics.Potential} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a {@link ThermostatEnum} object.
     * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum} object.
     * @return a {@link MolecularDynamics} object.
     */
    public static MolecularDynamics dynamicsFactory(MolecularAssembly assembly,
                                                    Potential potentialEnergy,
                                                    CompositeConfiguration properties,
                                                    AlgorithmListener listener,
                                                    ThermostatEnum requestedThermostat,
                                                    IntegratorEnum requestedIntegrator) {

        return dynamicsFactory(assembly, potentialEnergy, properties, listener,
                requestedThermostat, requestedIntegrator, defaultEngine(assembly, potentialEnergy));
    }

    /**
     * <p>
     * dynamicsFactory.</p>
     *
     * @param assembly            a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy     a {@link ffx.numerics.Potential} object.
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a {@link ThermostatEnum} object.
     * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum} object.
     * @param engine              a {@link MolecularDynamics.DynamicsEngine} object.
     * @return a {@link MolecularDynamics} object.
     */
    public static MolecularDynamics dynamicsFactory(MolecularAssembly assembly,
                                                    Potential potentialEnergy, CompositeConfiguration properties,
                                                    AlgorithmListener listener, ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator,
                                                    DynamicsEngine engine) {
        switch (engine) {
            case OPENMM:
                // TODO: Replace this with calls to the leaves of a proper tree structure.
                // Unfortunately, neither Java, nor Apache Commons, nor Guava has an arbitrary tree implementing Collection.
                // Nor does javax.swing have a quick "get me the leaves" method that I was able to find.
                boolean ommLeaves = (potentialEnergy instanceof ForceFieldEnergyOpenMM ||
                        potentialEnergy.getUnderlyingPotentials().stream().anyMatch((Potential p) -> p instanceof ForceFieldEnergyOpenMM));
                if (ommLeaves) {
                    return new MolecularDynamicsOpenMM(assembly,
                            potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);
                } else {
                    throw new IllegalArgumentException(format(" Requested OpenMM engine %s, but at least one leaf of the potential %s is not an OpenMM force field!", engine, potentialEnergy));
                }
            case FFX:
            default:
                return new MolecularDynamics(assembly,
                        potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);
        }
    }

    private static DynamicsEngine defaultEngine(MolecularAssembly molecularAssembly, Potential potentialEnergy) {
        CompositeConfiguration properties = molecularAssembly.getProperties();
        String mdEngine = properties.getString("MD-engine");
        if (mdEngine != null) {
            if (mdEngine.equalsIgnoreCase("OMM")) {
                logger.info(" Creating OpenMM Dynamics Object");
                return DynamicsEngine.OPENMM;
            } else {
                logger.info(" Creating FFX Dynamics Object");
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
     * @param properties          a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a {@link ThermostatEnum} object.
     * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum} object.
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
            barostat = (Barostat) potentialEnergy;
        }

        dynSleepTime = properties.getInt("dynamics-sleep-nanos", DEFAULT_DYNAMICS_SLEEP_TIME);

        assemblies.get(0).compositeConfiguration = properties;
        mass = potentialEnergy.getMass();
        numberOfVariables = potentialEnergy.getNumberOfVariables();
        x = new double[numberOfVariables];
        v = new double[numberOfVariables];
        a = new double[numberOfVariables];
        aPrevious = new double[numberOfVariables];
        gradient = new double[numberOfVariables];

        // If an Integrator wasn't passed to the MD constructor, check for one specified as a property.
        if (requestedIntegrator == null) {
            String integrate = properties.getString("integrate", "verlet").trim();
            try {
                requestedIntegrator = IntegratorEnum.valueOf(integrate);
            } catch (Exception e) {
                requestedIntegrator = IntegratorEnum.VERLET;
            }
        }
        boolean oMMLogging = false;
        if (potential instanceof ffx.potential.ForceFieldEnergyOpenMM) {
            oMMLogging = true;
        }

        List<Constraint> constraints = potentialEnergy.getConstraints();

        switch (requestedIntegrator) {
            case RESPA:
                Respa respa = new Respa(numberOfVariables, x, v, a, aPrevious, mass);
                int in = molecularAssembly.getProperties().getInt("respa-dt", 4);
                if (in < 2) {
                    in = 2;
                }
                if (!oMMLogging) {
                    respa.setInnerTimeSteps(in);
                }
                logger.log(Level.FINE, format(" Created a RESPA integrator with %d inner time steps.", in));
                integrator = respa;
                break;
            case STOCHASTIC:
                double friction = properties.getDouble("friction", 91.0);
                logger.log(Level.INFO, format(" Friction set at %.3f collisions/picosecond", friction));

                Stochastic stochastic = new Stochastic(friction, numberOfVariables, x, v, a, mass);
                if (properties.containsKey("randomseed")) {
                    stochastic.setRandomSeed(properties.getInt("randomseed", 0));
                }
                integrator = stochastic;
                // The stochastic dynamics integration procedure will thermostat
                // the system. The ADIABTIC thermostat just serves to report the
                // temperature and initialize velocities if necessary.
                requestedThermostat = ThermostatEnum.ADIABATIC;
                break;
            case BEEMAN:
                integrator = new BetterBeeman(numberOfVariables, x, v, a, aPrevious, mass);
                break;
            case VERLET:
            case VELOCITYVERLET:
            default:
                integrator = new VelocityVerlet(numberOfVariables, x, v, a, mass);
        }

        integrator.addConstraints(constraints);

        // If a Thermostat wasn't passed to the MD constructor, check for one specified as a property.
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
                thermostat = new Berendsen(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), targetTemperature, tau, constraints);
                break;
            case BUSSI:
                tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Bussi(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), targetTemperature, tau, constraints);
                if (properties.containsKey("randomseed")) {
                    thermostat.setRandomSeed(properties.getInt("randomseed", 0));
                }
                break;
            case ADIABATIC:
            default:
                thermostat = new Adiabatic(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), constraints);
        }

        if (properties.containsKey("randomseed")) {
            thermostat.setRandomSeed(properties.getInt("randomseed", 0));
        }

        // For Stochastic dynamics, center of mass motion will not be removed.
        if (integrator instanceof Stochastic) {
            thermostat.setRemoveCenterOfMassMotion(false);
        }

        numSnapshotsToKeep = properties.getInteger("dynamicsSnapshotMemory", 0);
        // Cannot construct a CircularFifoQueue of zero length.
        lastSnapshots = new CircularFifoQueue<>(Math.max(numSnapshotsToKeep, 1));

        verboseDynamicsState = properties.getBoolean("md-verbose", false);
        done = true;
    }

    /**
     * <p>
     * Constructor for MolecularDynamics.</p>
     *
     * @param assembly            a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy     a {@link ffx.numerics.Potential} object.
     * @param properties          a
     *                            {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     * @param listener            a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a {@link ThermostatEnum} object.
     * @param requestedIntegrator a
     *                            {@link ffx.algorithms.dynamics.integrators.IntegratorEnum} object.
     * @param esvSystem           a {@link ffx.potential.extended.ExtendedSystem} object.
     */
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

    public VerbosityLevel getVerbosityLevel() {
        return verbosityLevel;
    }

    public void setVerbosityLevel(VerbosityLevel level) {
        verbosityLevel = level;
        switch (level) {
            case SILENT: {
                intermediateLogging = Level.FINE;
                basicLogging = Level.FINE;
            }
            break;
            case QUIET: {
                intermediateLogging = Level.FINE;
                basicLogging = Level.INFO;
            }
            break;
            case VERBOSE:
            default: {
                intermediateLogging = Level.INFO;
                basicLogging = Level.INFO;
            }
            break;
        }
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
        gradient = new double[numberOfVariables];
        potential.getCoordinates(x);
        potential.getVelocity(v);
        potential.getAcceleration(a);
        potential.getPreviousAcceleration(aPrevious);
        if (potential instanceof ForceFieldEnergy) {
            gradient = ((ForceFieldEnergy) potential).getGradient(gradient);
        }
        thermostat.setNumberOfVariables(numberOfVariables, x, v, mass,
                potential.getVariableTypes(), true);
        integrator.setNumberOfVariables(numberOfVariables, x, v, a, aPrevious, mass);
    }

    /**
     * Not meaningful for FFX dynamics (no need to obtain velocities/accelerations from
     * a different program, especially one running on a GPU). Is a no-op.
     *
     * @param obtainVA Not meaningful for this implementation.
     */
    public void setObtainVelAcc(boolean obtainVA) {
        // Not meaningful for FFX dynamics.
    }

    /**
     * <p>
     * Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a
     *                   {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    /**
     * <p>
     * Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
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
        mi.compositeConfiguration = props;
        assemblies.add(mi);
    }

    /**
     * <p>
     * Setter for the field <code>monteCarloListener</code>.</p>
     *
     * @param listener a {@link MonteCarloListener} object.
     * @param when     a {@link MolecularDynamics.MonteCarloNotification} object.
     */
    public void setMonteCarloListener(MonteCarloListener listener, MonteCarloNotification when) {
        monteCarloListener = listener;
        mcNotification = when;
    }

    /**
     * <p>
     * attachExtendedSystem.</p>
     *
     * @param system         a {@link ffx.potential.extended.ExtendedSystem} object.
     * @param printFrequency a int.
     */
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        if (esvSystem != null) {
            logger.warning("An ExtendedSystem is already attached to this MD!");
        }
        esvSystem = system;
        printEsvFrequency = printFrequency;
        logger.info(format(" Attached extended system (%s) to molecular dynamics.", esvSystem.toString()));
        reInit();
    }

    /**
     * <p>
     * detachExtendedSystem.</p>
     */
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

        // Return if already running.
        if (!done) {
            logger.warning(" Programming error - attempt to modify parameters of a running MolecularDynamics instance.");
            return;
        }

        if (integrator instanceof Stochastic) {
            if (constantPressure) {
                logger.log(basicLogging,"\n Stochastic dynamics in the NPT ensemble");
            } else {
                logger.log(basicLogging,"\n Stochastic dynamics in the NVT ensemble");
            }
        } else if (!(thermostat instanceof Adiabatic)) {
            if (constantPressure) {
                logger.log(basicLogging,"\n Molecular dynamics in the NPT ensemble");
            } else {
                logger.log(basicLogging,"\n Molecular dynamics in the NVT ensemble");
            }
        } else {
            if (constantPressure) {
                logger.severe("\n NPT Molecular dynamics requires a thermostat");
            } else {
                logger.log(basicLogging,"\n Molecular dynamics in the NVE ensemble");
            }
        }

        this.nSteps = nSteps;
        totalSimTime = 0.0;

        // Convert the time step from femtoseconds to picoseconds.
        dt = timeStep * 1.0e-3;

        // Convert the print interval to a print frequency.
        printFrequency = 100;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        // Convert the save interval to a save frequency.
        saveSnapshotFrequency = 1000;
        if (saveInterval >= this.dt) {
            saveSnapshotFrequency = (int) (saveInterval / this.dt);
        }

        // Set snapshot file type.
        saveSnapshotAsPDB = true;
        if (fileType.equalsIgnoreCase("XYZ")) {
            saveSnapshotAsPDB = false;
        } else if (!fileType.equalsIgnoreCase("PDB")) {
            logger.warning("Snapshot file type unrecognized; saving snaphshots as PDB.\n");
        }

        // Convert restart frequency to steps.
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
            loadRestart = dyn.exists();
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        this.targetTemperature = temperature;
        this.initVelocities = initVelocities;
        done = false;

        if (dyn != null && dyn.exists()) {
            logger.info(" Continuing from " + dyn.getAbsolutePath());
        }

        if (!verbosityLevel.isQuiet) {
            logger.info(format(" Number of steps:     %8d", nSteps));
            logger.info(format(" Time step:           %8.3f (fsec)", timeStep));
            logger.info(format(" Print interval:      %8.3f (psec)", printInterval));
            logger.info(format(" Save interval:       %8.3f (psec)", saveInterval));
            //logger.info(format(" Archive file: %s", archiveFile.getName()));
            for (int i = 0; i < assemblies.size(); i++) {
                AssemblyInfo ai = assemblies.get(i);
                logger.info(format(" Archive file %3d: %s", i, ai.archiveFile.getName()));
            }
            logger.info(format(" Restart file:     %s", restartFile.getName()));
        }
    }

    /**
     * <p>
     * assemblyInfo.</p>
     */
    void assemblyInfo() {
        assemblies.stream().parallel().forEach((ainfo) -> {
            MolecularAssembly mola = ainfo.getAssembly();
            CompositeConfiguration aprops = ainfo.compositeConfiguration;
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
    @Deprecated
    public void redynamic(final int nSteps, final double temperature) {
        if (!verbosityLevel.isQuiet()) {
            setVerbosityLevel(VerbosityLevel.QUIET);
        }

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
        // Return if already running;
        // Could happen if two threads call dynamic on the same MolecularDynamics instance.
        if (!done) {
            logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
            return;
        }

        init(nSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency,
                temperature, initVelocities, dyn);

        Thread dynamicThread = new Thread(this);
        dynamicThread.start();
        synchronized (this) {
            try {
                while (dynamicThread.isAlive()) {
                    wait(0, dynSleepTime);
                }
            } catch (InterruptedException e) {
                String message = " Molecular dynamics interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }
        if (!verbosityLevel.isQuiet) {
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
     * @throws java.lang.IllegalArgumentException If restart frequency is not a
     *                                            positive number
     */
    public void setRestartFrequency(double restartFrequency) throws IllegalArgumentException {
        if (Double.isFinite(restartFrequency) && restartFrequency > 0) {
            this.restartFrequency = restartFrequency;
        } else {
            throw new IllegalArgumentException(format(" Restart frequency must be positive finite, was %10.4g", restartFrequency));
        }
    }

    /**
     * Performs basic pre-MD operations such as loading the restart file.
     */
    void preRunOps() throws IllegalStateException {
        done = false;
        terminate = false;

        // Set the target temperature.
        thermostat.setTargetTemperature(targetTemperature);
        boolean quiet = verbosityLevel.isQuiet();
        thermostat.setQuiet(quiet);
        if (integrator instanceof Stochastic) {
            Stochastic stochastic = (Stochastic) integrator;
            stochastic.setTemperature(targetTemperature);
        }

        // Set the step size.
        integrator.setTimeStep(dt);

        if (!initialized) {
            // Initialize from a restart file.
            if (loadRestart) {
                Crystal crystal = molecularAssembly.getCrystal();
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    throw new IllegalStateException(message);
                } else {
                    molecularAssembly.getPotentialEnergy().setCrystal(crystal);
                }
            } else {
                // Initialize from using current atomic coordinates.
                potential.getCoordinates(x);
                // Initialize atomic velocities from a Maxwell-Boltzmann distribution or set to 0.
                if (initVelocities) {
                    thermostat.maxwell(targetTemperature);
                } else {
                    fill(v, 0.0);
                }
            }
        } else {
            // If MD has already been run (ie. Annealing or RepEx), then initialize velocities if requested.
            if (initVelocities) {
                thermostat.maxwell(targetTemperature);
            }
        }
    }

    /**
     * Initializes energy fields, esp. potential energy.
     */
    private void initializeEnergies() {
        // Compute the current potential energy.
        try {
            currentPotentialEnergy = potential.energyAndGradient(x, gradient);
        } catch (EnergyException ex) {
            writeStoredSnapshots();
            throw ex;
        }

        // Initialize current and previous accelerations.
        if (!loadRestart || initialized || integrator instanceof Respa) {
            // For the Respa integrator, initial accelerations are from the slowly varying forces.
            if (integrator instanceof Respa) {
                potential.setEnergyTermState(Potential.STATE.SLOW);
                potential.energyAndGradient(x, gradient);
            }

            for (int i = 0; i < numberOfVariables; i++) {
                a[i] = -KCAL_TO_GRAM_ANG2_PER_PS2 * gradient[i] / mass[i];
            }

            if (aPrevious != null) {
                arraycopy(a, 0, aPrevious, 0, numberOfVariables);
            }
        }

        // Compute the current kinetic energy.
        thermostat.kineticEnergy();
        currentKineticEnergy = thermostat.getKineticEnergy();
        startingKineticEnergy = currentKineticEnergy;
        currentTemperature = thermostat.getCurrentTemperature();
        currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

        startingPotentialEnergy = currentPotentialEnergy;
        startingTotalEnergy = currentTotalEnergy;
    }

    /**
     * Pre-run operations (mostly logging) that require knowledge of system energy.
     */
    void postInitEnergies() {
        initialized = true;
        logger.log(basicLogging, format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
        logger.log(basicLogging, format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
        logger.log(basicLogging, format("  %8s %12.4f %12.4f %12.4f %8.2f",
                "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));

        // Store the initialized state.
        storeState();
    }

    /**
     * Post-run cleanup operations.
     */
    void postRun() {
        // Add the potential energy of the slow degrees of freedom.
        if (integrator instanceof Respa) {
            potential.setEnergyTermState(Potential.STATE.BOTH);
        }

        // Log normal completion.
        if (!terminate) {
            logger.log(basicLogging, format(" Completed %8d time steps\n", nSteps));
        }

        // Reset the done and terminate flags.
        done = true;
        terminate = false;

        if (monteCarloListener != null && mcNotification == MonteCarloNotification.AFTER_DYNAMICS) {
            monteCarloListener.mcUpdate(thermostat.getCurrentTemperature());
        }
    }

    /**
     * Main loop of the run method.
     */
    private void mainLoop() {
        // Integrate Newton's equations of motion for the requested number of steps,
        // unless early termination is requested.
        long time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {
            if (step > 1) {
                List<Constraint> constraints = potential.getConstraints();
                // TODO: Replace magic numbers with named constants.
                long constraintFails = constraints.stream().
                        filter((Constraint c) -> !c.constraintSatisfied(x, v, 1E-7, 1E-7)).
                        count();
                if (constraintFails > 0) {
                    logger.info(format(" %d constraint failures in step %d", constraintFails, step));
                }
            }
            /* Notify MonteCarlo handlers such as PhMD or rotamer drivers. */
            if (monteCarloListener != null && mcNotification == MonteCarloNotification.EACH_STEP) {
                monteCarloListener.mcUpdate(thermostat.getCurrentTemperature());
            }

            // Do the half-step thermostat operation.
            thermostat.halfStep(dt);

            // Do the half-step integration operation.
            integrator.preForce(potential);

            // Compute the potential energy and gradients.
            double priorPE = currentPotentialEnergy;
            try {
                currentPotentialEnergy = potential.energyAndGradient(x, gradient);
            } catch (EnergyException ex) {
                writeStoredSnapshots();
                throw ex;
            }

            // Add the potential energy of the slow degrees of freedom.
            if (integrator instanceof Respa) {
                Respa r = (Respa) integrator;
                currentPotentialEnergy += r.getHalfStepEnergy();
            }

            double defaultDeltaPEThresh = 1.0E6;
            detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

            // Do the full-step integration operation.
            integrator.postForce(gradient);

            // Compute the full-step kinetic energy.
            thermostat.kineticEnergy();

            // Do the full-step thermostat operation.
            thermostat.fullStep(dt);

            // Recompute the kinetic energy after the full-step thermostat operation.
            thermostat.kineticEnergy();

            // Remove center of mass motion ever ~100 steps.
            int removeCOMMotionFrequency = 100;
            if (thermostat.getRemoveCenterOfMassMotion() && step % removeCOMMotionFrequency == 0) {
                thermostat.centerOfMassMotion(true, false);
            }

            // Collect current kinetic energy, temperature, and total energy.
            currentKineticEnergy = thermostat.getKineticEnergy();
            currentTemperature = thermostat.getCurrentTemperature();
            currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

            // Update atomic velocity, acceleration and previous acceleration.
            potential.setVelocity(v);
            potential.setAcceleration(a);
            potential.setPreviousAcceleration(aPrevious);

            // Update extended system variables if present.
            if (esvSystem != null) {
                esvSystem.propagateESVs(currentTemperature, dt, step * dt);
            }

            // Log the current state every printFrequency steps.
            totalSimTime += dt;
            if (step % printFrequency == 0) {
                // Original print statement
                time = System.nanoTime() - time;
                logger.log(basicLogging, format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.3f",
                        totalSimTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * NS2SEC));
                time = System.nanoTime();
            }
            if (step % printEsvFrequency == 0 && esvSystem != null) {
                logger.log(basicLogging, format(" %7.3e %s", totalSimTime, esvSystem.getLambdaList()));
            }

            // Write out snapshots in selected format every saveSnapshotFrequency steps.
            if (saveSnapshotFrequency > 0 && step % saveSnapshotFrequency == 0) {
                for (AssemblyInfo ai : assemblies) {
                    if (ai.archiveFile != null && !saveSnapshotAsPDB) {
                        if (ai.xyzFilter.writeFile(ai.archiveFile, true)) {
                            logger.log(basicLogging, format(" Appended snap shot to %s", ai.archiveFile.getName()));
                        } else {
                            logger.warning(format(" Appending snap shot to %s failed", ai.archiveFile.getName()));
                        }
                    } else if (saveSnapshotAsPDB) {
                        if (ai.pdbFilter.writeFile(ai.pdbFile, false)) {
                            logger.log(basicLogging, format(" Wrote PDB file to %s", ai.pdbFile.getName()));
                        } else {
                            logger.warning(format(" Writing PDB file to %s failed.", ai.pdbFile.getName()));
                        }
                    }
                }
            }

            // Write out restart files every saveRestartFileFrequency steps.
            if (saveRestartFileFrequency > 0 && step % saveRestartFileFrequency == 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.log(basicLogging, " Wrote dynamics restart file to " + restartFile.getName());
                } else {
                    logger.log(basicLogging, " Writing dynamics restart file to " + restartFile.getName() + " failed");
                }
            }

            // Notify the algorithmListener.
            if (algorithmListener != null && step % printFrequency == 0) {
                //algorithmListener.algorithmUpdate(molecularAssembly);
                for (AssemblyInfo assembly : assemblies) {
                    // Probably unwise to parallelize this, so that it doesn't
                    // hit the GUI with parallel updates.
                    algorithmListener.algorithmUpdate(assembly.getAssembly());
                }
            }

            // Check for a termination request.
            if (terminate) {
                logger.info(format("\n Terminating after %8d time steps\n", step));
                break;
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        try {
            preRunOps();
        } catch (IllegalStateException ise) {
            return;
        }
        initializeEnergies();
        postInitEnergies();
        mainLoop();
        postRun();
    }

    /**
     * Detects grossly atypical potential energy values that are likely
     * incorrect, and writes snapshots to disc. "Grossly atypical" is defined
     * as: greater than 1.0E100 kcal/mol, less than -1.0E100 kcal/mol,
     * non-finite (NaN/infinite), or exceeding specified delta from the prior
     * potential energy.
     * <p>
     * After prior snapshots have been written to disc, the queue they are
     * stored on is now empty, preventing printing of duplicate snapshots.
     *
     * @param priorPE     Potential energy prior to this step.
     * @param delPEThresh If potential energy changes by this much, trigger a
     */
    void detectAtypicalEnergy(double priorPE, double delPEThresh) {

        // If not keeping snapshots, disable functionality.
        if (numSnapshotsToKeep < 1) {
            return;
        }

        double deltaPE = currentPotentialEnergy - priorPE;

        CoordinateSnapshot currState = new CoordinateSnapshot();
        currState.storeState();
        lastSnapshots.add(currState);

        double maxPEThresh = 1.0E100; // 1.0E100 kcal/mol is well into the territory of the absurd.
        double absPE = Math.abs(currentPotentialEnergy);

        if (absPE > maxPEThresh || !Double.isFinite(currentPotentialEnergy) || Math.abs(deltaPE) > delPEThresh) {
            logger.info(format(" Unusual potential energy %12.5g detected, writing snapshots.", currentPotentialEnergy));
            writeStoredSnapshots();
            currState.revertState(); // May be unnecessary, thanks to the current state always being last on the queue.
            if (absPE > 1.0E100 || !Double.isFinite(currentPotentialEnergy)) {
                logger.severe(format(" Dynamics exiting with atypical potential energy of %12.5g", currentPotentialEnergy));
            }
        }
    }

    /**
     * Performs the inner loop of writing snapshots to disk; used by both
     * detectAtypicalEnergy and a try-catch in dynamics.
     */
    private void writeStoredSnapshots() {
        int numSnaps = lastSnapshots.size();

        File origFile = molecularAssembly.getFile();
        String timeString = LocalDateTime.now().format(DateTimeFormatter.
                ofPattern("HH_mm_ss"));
        PotentialsFunctions potentialsFunctions = new PotentialsUtils();

        String filename = format("%s-%s-SNAP.pdb",
                FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                timeString);

        for (int is = 0; is < numSnaps; is++) {
            CoordinateSnapshot oldState = lastSnapshots.poll();
            if (oldState != null) {
                oldState.revertState();
            }
            potentialsFunctions.saveAsPDB(molecularAssembly, new File(potentialsFunctions.versionFile(filename)));
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
     * Get the system kinetic energy.
     *
     * @return kinetic energy.
     */
    public double getStartingKineticEnergy() {
        return startingKineticEnergy;
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
     * Get the starting system potential energy.
     *
     * @return potential energy.
     */
    public double getStartingPotentialEnergy() {
        return startingPotentialEnergy;
    }

    /**
     * Get the current temperature of the system
     *
     * @return currentTemperature
     */
    public double getTemperature() {
        return currentTemperature;
    }

    /**
     * Returns the associated dynamics file.
     * @return
     */
    public File getDynFile() {
        return restartFile;
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
     * A simple container class to hold all the infrastructure associated with a
     * MolecularAssembly for MolecularDynamics; assembly, properties, archive
     * and PDB files, PDB and XYZ filters. Direct access to package-private
     * members breaks encapsulation a bit, but the private inner class shouldn't
     * be used externally anyways.
     */
    protected class AssemblyInfo {

        private final MolecularAssembly assembly;
        CompositeConfiguration compositeConfiguration;
        File archiveFile = null;
        File pdbFile;
        PDBFilter pdbFilter;
        XYZFilter xyzFilter = null;

        AssemblyInfo(MolecularAssembly assembly) {
            this.assembly = assembly;
            pdbFile = this.assembly.getFile();
            compositeConfiguration = this.assembly.getProperties();
            pdbFilter = new PDBFilter(this.assembly.getFile(), this.assembly,
                    this.assembly.getForceField(), this.assembly.getProperties());
        }

        public MolecularAssembly getAssembly() {
            return assembly;
        }
    }

    /**
     * <p>
     * storeState.</p>
     */
    public void storeState() {
        if (dynamicsState == null) {
            dynamicsState = new DynamicsState();
        }
        dynamicsState.storeState();
    }

    /**
     * <p>
     * revertState.</p>
     *
     * @throws java.lang.Exception if any.
     */
    public void revertState() throws Exception {
        if (dynamicsState == null) {
            throw new Exception();
        }
        dynamicsState.revertState();
    }

    /**
     * More limited version of a DynamicsState, storing only coordinates.
     * TODO: Make DynamicsState more flexible and let it store any combination of variables.
     */
    protected class CoordinateSnapshot {
        final double[] xBak;

        CoordinateSnapshot() {
            xBak = new double[numberOfVariables];
        }

        void storeState() {
            arraycopy(x, 0, xBak, 0, numberOfVariables);
        }

        void revertState() {
            arraycopy(xBak, 0, x, 0, numberOfVariables);
            Atom[] atoms = molecularAssembly.getActiveAtomArray();
            for (int i = 0; i < atoms.length; i++) {
                int i3 = 3*i;
                double[] newXYZ = new double[3];
                arraycopy(xBak, i3, newXYZ, 0, 3);
                atoms[i].setXYZ(newXYZ);
            }
        }
    }

    protected class DynamicsState {

        double[] xBak, vBak, aBak;
        double[] aPreviousBak, massBak, gradBak;
        double currentKineticEnergyBak, currentPotentialEnergyBak, currentTotalEnergyBak;
        double currentTemperatureBak;

        DynamicsState() {
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
            arraycopy(gradient, 0, gradBak, 0, numberOfVariables);
            if (verboseDynamicsState) {
                describe(" Storing State:");
            }
        }

        public void describe(String title) {
            StringBuilder sb = new StringBuilder();
            sb.append(title);
            sb.append("\nx: ");
            Arrays.stream(x).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\nv: ");
            Arrays.stream(v).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\na: ");
            Arrays.stream(a).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\naP: ");
            Arrays.stream(aPrevious).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\nm: ");
            Arrays.stream(mass).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\ng: ");
            Arrays.stream(gradient).forEach(val -> sb.append(format("%.2g, ", val)));
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
            arraycopy(gradBak, 0, gradient, 0, numberOfVariables);

            Atom[] atoms = molecularAssembly.getActiveAtomArray();
            if (atoms.length * 3 == numberOfVariables) {
                double[] vec = new double[3];
                int index = 0;
                for (Atom atom : atoms) {
                    atom.moveTo(x[index], x[index + 1], x[index + 2]);
                    atom.setXYZGradient(gradient[index], gradient[index + 1], gradient[index + 2]);
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
                }
            }
            if (verboseDynamicsState) {
                describe(" Reverting State (To):");
            }
        }
    }

    /**
     * Enumerates available molecular dynamics engines; presently limited to the
     * FFX reference engine and the OpenMM engine.
     * <p>
     * Distinct from the force field energy Platform, as the FFX engine can use
     * OpenMM energies, but not vice-versa.
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
         * @param platform The requested platform.
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

    /**
     * <p>
     * getStartingTotalEnergy.</p>
     *
     * @return a double.
     */
    public double getStartingTotalEnergy() {
        return 0.0;
    }

    /**
     * <p>
     * getEndTotalEnergy.</p>
     *
     * @return a double.
     */
    public double getEndTotalEnergy() {
        return 0.0;
    }

    /**
     * @param intervalSteps
     */
    public void setIntervalSteps(int intervalSteps) {
        // Not meaningful for FFX MD.
    }

    /**
     * <p>
     * getTimeStep.</p>
     *
     * @return a double.
     */
    public double getTimeStep() {
        return dt;
    }

    /**
     * <p>
     * getIntervalSteps.</p>
     *
     * @return a int.
     */
    public int getIntervalSteps() {
        return 1;
    }

    /**
     * <p>
     * getNumAtoms.</p>
     *
     * @return a int.
     */
    public int getNumAtoms() {
        return 1;
    }

    /**
     * Write out restart files.
     */
    public void writeRestart() {
        writeSnapShot();
        writeDynamicsRestart();
    }

    /**
     * Write out snapshots only if lambda is greater than the lambda writeout
     * value.
     *
     * @param lambda         Current value of lambda.
     * @param lambdaWriteOut The lambda write out cut-off.
     */
    public void writeLambdaThresholdRestart(double lambda, double lambdaWriteOut) {
        if (lambda >= lambdaWriteOut) {
            writeSnapShot();
        }
        writeDynamicsRestart();
    }

    /**
     * Write out coordinate snapshots.
     */
    private void writeSnapShot() {
        for (AssemblyInfo ai : assemblies) {
            if (ai.archiveFile != null && !saveSnapshotAsPDB) {
                if (ai.xyzFilter.writeFile(ai.archiveFile, true)) {
                    logger.info(format(" Appended snap shot to %s", ai.archiveFile.getName()));
                } else {
                    logger.warning(format(" Appending snap shot to %s failed", ai.archiveFile.getName()));
                }
            } else if (saveSnapshotAsPDB) {
                if (ai.pdbFilter.writeFile(ai.pdbFile, false)) {
                    logger.info(format(" Wrote PDB file to %s", ai.pdbFile.getName()));
                } else {
                    logger.warning(format(" Writing PDB file to %s failed.", ai.pdbFile.getName()));
                }
            }
        }
    }

    /**
     * Write out Dynamics restart.
     */
    private void writeDynamicsRestart() {
        if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
            logger.info(" Wrote dynamics restart file to " + restartFile.getName());
            if (constantPressure) {
                Crystal crystal = molecularAssembly.getCrystal();
                double currentDensity = crystal.getDensity(molecularAssembly.getTotalMass());
                logger.info(format(" Density %6.3f (g/cc) with unit cell %s.",
                        currentDensity, crystal.toShortString()));
            }
        } else {
            logger.info(" Writing dynamics restart file to " + restartFile.getName() + " failed.");
        }
    }

    public enum VerbosityLevel {
        VERBOSE(false), QUIET(true), SILENT(true);

        private boolean isQuiet;

        VerbosityLevel(boolean isQuiet) {
            this.isQuiet = isQuiet;
        }

        public boolean isQuiet() {
            return isQuiet;
        }
    }
}
