// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
// ******************************************************************************
package ffx.algorithms.dynamics;

import static ffx.utilities.Constants.FSEC_TO_PSEC;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import edu.rit.pj.Comm;
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
import ffx.crystal.Crystal;
import ffx.numerics.Constraint;
import ffx.numerics.Potential;
import ffx.numerics.math.RunningStatistics;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XPHFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FileUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * Run NVE, NVT, or NPT molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MolecularDynamics implements Runnable, Terminatable {

  private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());

  private static final double DEFAULT_LOG_INTERVAL = 0.25;
  private static final double DEFAULT_RESTART_INTERVAL = 1.0;
  private static final double DEFAULT_TRAJECTORY_INTERVAL = 10.0;
  /** By default, wait 25000 nsec (25 usec) in between polling the dynamics thread. */
  private static final int DEFAULT_DYNAMICS_SLEEP_TIME = 25000;

  /**
   * MolecularAssembly to run dynamics on.
   * <p>
   * For dual or quad topology simulations, this is an array of length 2 or 4.
   */
  protected MolecularAssembly[] molecularAssembly;
  /** Propagate dynamics on this potential surface. */
  private final Potential potential;
  /** An Algorithm Listener to send updates to the GUI. */
  protected final AlgorithmListener algorithmListener;

  /**
   * Number of dynamics variables. The length of x, v, a, aPrevious, gradient, and mass.
   */
  int numberOfVariables;
  /** Coordinates. */
  protected double[] x;
  /** Velocities. */
  protected double[] v;
  /** Accelerations. */
  protected double[] a;
  /** Previous accelerations. */
  protected double[] aPrevious;
  /** The gradient. */
  protected double[] gradient;
  /** Mass for each degree of freedom. */
  protected double[] mass;

  /** Time step (picoseconds). Always use the <code>setTimestep</code> to set! */
  protected double dt = 0.001;
  /** Number of MD steps to take. */
  private long nSteps = 1000;
  /** Target temperature. ToDo: use the Thermostat instance. */
  double targetTemperature = 298.15;
  /** Flag to indicate use of constant pressure. */
  boolean constantPressure = false;
  /** Flag to indicate velocities should be initialized. */
  protected boolean initVelocities = true;
  /**
   * Whether MD handles writing restart/trajectory files itself (true), or will be commanded by
   * another class (false) to do it. The latter is true for MC-OST, for example.
   */
  protected boolean automaticWriteouts = true;
  /** Total simulation time. */
  protected double totalSimTime = 0.0;

  /** Time between writing out restart/checkpoint files in picoseconds. */
  protected double restartInterval = DEFAULT_RESTART_INTERVAL;
  /** Timesteps between writing out restart/checkpoint files. Set by the init method. */
  protected int restartFrequency;
  /** Time between appending to the trajectory file in picoseconds. */
  private double trajectoryInterval = DEFAULT_TRAJECTORY_INTERVAL;
  /** Timesteps between adding a frame to the trajectory file. Set by the init method. */
  protected int trajectoryFrequency;
  /** Time between logging information to the screen in picoseconds. */
  private double logInterval = DEFAULT_LOG_INTERVAL;
  /** TImesteps between logging information to the screen. Set by the init method. */
  protected int logFrequency;

  /** Temperature of the system as of the last init call (and start of the last dynamics run). */
  protected double initialTemp;
  /** Kinetic energy of the system as of the last init call (and start of the last dynamics run). */
  protected double initialKinetic;
  /**
   * Potential energy of the system as of the last init call (and start of the last dynamics run).
   */
  protected double initialPotential;
  /** Total energy of the system as of the last init call (and start of the last dynamics run). */
  protected double initialTotal;
  /** Current temperature. */
  double currentTemperature;
  /** Current kinetic energy. */
  double currentKineticEnergy;
  /** Current potential energy. */
  double currentPotentialEnergy;
  /** Current total energy. */
  double currentTotalEnergy;
  private final RunningStatistics temperatureStats = new RunningStatistics();
  private final RunningStatistics potentialEnergyStats = new RunningStatistics();
  private final RunningStatistics kineticEnergyStats = new RunningStatistics();
  private final RunningStatistics totalEnergyStats = new RunningStatistics();

  /** File type to use when saving files. */
  protected String fileType = "XYZ";
  /** Save snapshots in PDB format. */
  boolean saveSnapshotAsPDB = true;

  /** Dynamics restart file. */
  File restartFile = null;
  /** Flag to indicate loading of restart file. */
  boolean loadRestart = false;
  /** Filter to parse the dynamics restart file. */
  DYNFilter dynFilter = null;
  /** If asked to perform dynamics with a null dynamics file, write here. */
  private File fallbackDynFile;

  /**
   * Log basic information at this level.
   */
  Level basicLogging = Level.INFO;
  /** Indicates how verbose MD should be. */
  private MDVerbosity verbosityLevel = MDVerbosity.VERBOSE;

  /** Flag to indicate dynamics has been initialized. */
  boolean initialized = false;
  /** Flag to indicate MD should be terminated. */
  private boolean terminate = false;
  /** Flag to indicate a run has finished. */
  protected boolean done;
  /** Wait this many nanoseconds in between polling the dynamics thread. */
  private final int dynSleepTime;

  /** Integrator instance. */
  private final Integrator integrator;
  /** Thermostat instance. */
  private Thermostat thermostat;
  /** Any Barostat that may be in use. */
  Barostat barostat = null;
  /** State of the dynamics. */
  private MDState mdState;

  /** ESV System. */
  private ExtendedSystem esvSystem;
  /** ESV Stochastic integrator. */
  private Stochastic esvIntegrator;
  /** ESV Adiabatic thermostat. */
  private Adiabatic esvThermostat;
  /** Frequency to print ESV info. */
  private int printEsvFrequency = -1;

  /**
   * Constructor for MolecularDynamics.
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param potentialEnergy a {@link ffx.numerics.Potential} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param requestedThermostat a {@link ThermostatEnum} object.
   * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
   *     object.
   */
  public MolecularDynamics(MolecularAssembly assembly, Potential potentialEnergy,
      CompositeConfiguration properties, AlgorithmListener listener,
      ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator) {
    this(assembly, potentialEnergy, properties, listener, requestedThermostat,
        requestedIntegrator, defaultFallbackDyn(assembly));
  }

  /**
   * Constructor for MolecularDynamics.
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param potentialEnergy a {@link ffx.numerics.Potential} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param requestedThermostat a {@link ThermostatEnum} object.
   * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
   *     object.
   * @param fallbackDyn File to write restarts to if none is provided by the dynamics method.
   */
  public MolecularDynamics(MolecularAssembly assembly, Potential potentialEnergy,
      CompositeConfiguration properties, AlgorithmListener listener,
      ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator, File fallbackDyn) {
    this.molecularAssembly = new MolecularAssembly[1];
    this.molecularAssembly[0] = assembly;
    this.algorithmListener = listener;
    this.potential = potentialEnergy;
    if (potentialEnergy instanceof Barostat) {
      constantPressure = true;
      barostat = (Barostat) potentialEnergy;
    }

    dynSleepTime = properties.getInt("dynamics-sleep-nanos", DEFAULT_DYNAMICS_SLEEP_TIME);

    mass = potentialEnergy.getMass();
    numberOfVariables = potentialEnergy.getNumberOfVariables();
    x = new double[numberOfVariables];
    v = new double[numberOfVariables];
    a = new double[numberOfVariables];
    aPrevious = new double[numberOfVariables];
    gradient = new double[numberOfVariables];

    // If an Integrator wasn't passed to the MD constructor, check for one specified as a property.
    if (requestedIntegrator == null) {
      String integrate = properties.getString("integrate", "VERLET").trim();
      try {
        integrate = integrate.toUpperCase().replace("-", "_");
        requestedIntegrator = IntegratorEnum.valueOf(integrate.toUpperCase());
      } catch (Exception e) {
        requestedIntegrator = IntegratorEnum.VERLET;
      }
    }

    boolean oMMLogging = potential instanceof ForceFieldEnergyOpenMM;
    List<Constraint> constraints = potentialEnergy.getConstraints();

    // Set the integrator.
    integrator = switch (requestedIntegrator) {
      case RESPA, MTS -> {
        Respa respa = new Respa(numberOfVariables, x, v, a, aPrevious, mass);
        int in = molecularAssembly[0].getProperties().getInt("respa-dt", 4);
        if (in < 2) {
          in = 2;
        }
        if (!oMMLogging) {
          respa.setInnerTimeSteps(in);
        }
        logger.log(Level.FINE, format(" Created a RESPA integrator with %d inner time steps.", in));
        yield respa;
      }
      case STOCHASTIC, LANGEVIN -> {
        double friction = properties.getDouble("friction", 91.0);
        logger.log(Level.FINE, format(" Friction set at %.3f collisions/picosecond", friction));
        Stochastic stochastic = new Stochastic(friction, numberOfVariables, x, v, a, mass);
        if (properties.containsKey("randomseed")) {
          stochastic.setRandomSeed(properties.getInt("randomseed", 0));
        }
        // The stochastic dynamics integration procedure will thermostat
        // the system. The ADIABTIC thermostat just serves to report the
        // temperature and initialize velocities if necessary.
        requestedThermostat = ThermostatEnum.ADIABATIC;
        yield stochastic;
      }
      case BEEMAN -> new BetterBeeman(numberOfVariables, x, v, a, aPrevious, mass);
      default -> new VelocityVerlet(numberOfVariables, x, v, a, mass);
    };

    integrator.addConstraints(constraints);

    // If a Thermostat wasn't passed to the MD constructor, check for one specified as a property.
    if (requestedThermostat == null) {
      String thermo = properties.getString("thermostat", "Berendsen").trim();
      try {
        thermo = thermo.toUpperCase().replace("-", "_");
        requestedThermostat = ThermostatEnum.valueOf(thermo.toUpperCase());
      } catch (Exception e) {
        requestedThermostat = ThermostatEnum.BERENDSEN;
      }
    }

    // Set the thermostat.
    thermostat = switch (requestedThermostat) {
      case BERENDSEN -> {
        double tau = properties.getDouble("tau-temperature", 0.2);
        yield new Berendsen(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(),
            targetTemperature, tau, constraints);
      }
      case BUSSI -> {
        double tau = properties.getDouble("tau-temperature", 0.2);
        Bussi bussi = new Bussi(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(),
            targetTemperature, tau, constraints);
        if (properties.containsKey("randomseed")) {
          bussi.setRandomSeed(properties.getInt("randomseed", 0));
        }
        yield bussi;
      }
      default -> new Adiabatic(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(),
          constraints);
    };

    if (properties.containsKey("randomseed")) {
      thermostat.setRandomSeed(properties.getInt("randomseed", 0));
    }

    // For Stochastic dynamics, center of mass motion will not be removed.
    if (integrator instanceof Stochastic) {
      thermostat.setRemoveCenterOfMassMotion(false);
    }

    done = true;
    fallbackDynFile = fallbackDyn;
  }

  /**
   * Constructor for MolecularDynamics.
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param potentialEnergy a {@link ffx.numerics.Potential} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param requestedThermostat a {@link ThermostatEnum} object.
   * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
   *     object.
   * @param esvSystem a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public MolecularDynamics(MolecularAssembly assembly, Potential potentialEnergy,
      CompositeConfiguration properties, AlgorithmListener listener,
      ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator,
      ExtendedSystem esvSystem) {
    this(assembly, potentialEnergy, properties, listener, requestedThermostat,
        requestedIntegrator);
    this.esvSystem = esvSystem;
  }

  /**
   * Method that determines whether a dynamics is done by the java implementation native to ffx or
   * the OpenMM implementation
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param potentialEnergy a {@link ffx.numerics.Potential} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param requestedThermostat a {@link ThermostatEnum} object.
   * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
   *     object.
   * @return a {@link MolecularDynamics} object.
   */
  public static MolecularDynamics dynamicsFactory(MolecularAssembly assembly,
      Potential potentialEnergy, CompositeConfiguration properties, AlgorithmListener listener,
      ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator) {

    return dynamicsFactory(assembly, potentialEnergy, properties, listener, requestedThermostat,
        requestedIntegrator, defaultEngine(assembly, potentialEnergy));
  }

  /**
   * dynamicsFactory.
   *
   * @param assembly a {@link ffx.potential.MolecularAssembly} object.
   * @param potentialEnergy a {@link ffx.numerics.Potential} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
   * @param requestedThermostat a {@link ThermostatEnum} object.
   * @param requestedIntegrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
   *     object.
   * @param engine a {@link MDEngine} object.
   * @return a {@link MolecularDynamics} object.
   */
  public static MolecularDynamics dynamicsFactory(MolecularAssembly assembly,
      Potential potentialEnergy, CompositeConfiguration properties, AlgorithmListener listener,
      ThermostatEnum requestedThermostat, IntegratorEnum requestedIntegrator, MDEngine engine) {

    switch (engine) {
      case OPENMM, OMM -> {
        // TODO: Replace this with calls to the leaves of a proper tree structure.
        // Unfortunately, neither Java, nor Apache Commons, nor Guava has an arbitrary tree
        // implementing Collection.
        // Nor does javax.swing have a quick "get me the leaves" method that I was able to find.
        boolean ommLeaves = potentialEnergy.getUnderlyingPotentials().stream()
            .anyMatch((Potential p) -> p instanceof ForceFieldEnergyOpenMM);
        ommLeaves = ommLeaves || potentialEnergy instanceof ForceFieldEnergyOpenMM;
        if (ommLeaves) {
          return new MolecularDynamicsOpenMM(assembly, potentialEnergy, properties, listener,
              requestedThermostat, requestedIntegrator);
        } else {
          throw new IllegalArgumentException(format(
              " Requested OpenMM engine %s, but at least one leaf of the potential %s is not an OpenMM force field!",
              engine, potentialEnergy));
        }
      }
      default -> {
        return new MolecularDynamics(assembly, potentialEnergy, properties, listener,
            requestedThermostat, requestedIntegrator);
      }
    }
  }

  private static MDEngine defaultEngine(MolecularAssembly molecularAssembly,
      Potential potentialEnergy) {
    CompositeConfiguration properties = molecularAssembly.getProperties();
    String mdEngine = properties.getString("MD-engine");
    if (mdEngine != null) {
      if (mdEngine.equalsIgnoreCase("OMM")) {
        logger.info(" Creating OpenMM Dynamics Object");
        return MDEngine.OPENMM;
      } else {
        logger.info(" Creating FFX Dynamics Object");
        return MDEngine.FFX;
      }
    } else {
      // TODO: Replace this with a better check.
      boolean ommLeaves = potentialEnergy.getUnderlyingPotentials().stream()
          .anyMatch((Potential p) -> p instanceof ForceFieldEnergyOpenMM);
      ommLeaves = ommLeaves || potentialEnergy instanceof ForceFieldEnergyOpenMM;
      if (ommLeaves) {
        return MDEngine.OPENMM;
      } else {
        return MDEngine.FFX;
      }
    }
  }

  /**
   * Creates the "default" fallback dynamics file object (does not create actual file!)
   *
   * @param assembly First assembly.
   * @return Default fallback file.
   */
  private static File defaultFallbackDyn(MolecularAssembly assembly) {
    String firstFileName = FilenameUtils.removeExtension(assembly.getFile().getAbsolutePath());
    return new File(firstFileName + ".dyn");
  }

  /**
   * Adds a MolecularAssembly to be tracked by this MolecularDynamics. Note: does not affect the
   * underlying Potential.
   *
   * @param assembly A MolecularAssembly to be tracked.
   */
  public void addAssembly(MolecularAssembly assembly) {
    molecularAssembly = Arrays.copyOf(molecularAssembly, molecularAssembly.length + 1);
    molecularAssembly[molecularAssembly.length - 1] = assembly;
  }

  /**
   * attachExtendedSystem.
   *
   * @param system a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public void attachExtendedSystem(ExtendedSystem system, double reportFreq) {
    if (esvSystem != null) {
      logger.warning("An ExtendedSystem is already attached to this MD!");
    }
    esvSystem = system;
    this.esvIntegrator = new Stochastic(esvSystem.getThetaFriction(),
        esvSystem.getNumberOfVariables(), esvSystem.getThetaPosition(), esvSystem.getThetaVelocity(),
        esvSystem.getThetaAccel(), esvSystem.getThetaMassArray());
    this.esvThermostat = new Adiabatic(esvSystem.getNumberOfVariables(),
        esvSystem.getThetaPosition(), esvSystem.getThetaVelocity(), esvSystem.getThetaMassArray(),
        potential.getVariableTypes());
    printEsvFrequency = intervalToFreq(reportFreq, "Reporting (logging) interval");
    logger.info(
        format("  Attached extended system (%s) to molecular dynamics.", esvSystem.toString()));
    logger.info(format("  Extended System Theta Friction: %f", esvSystem.getThetaFriction()));
    logger.info(format("  Extended System Theta Mass: %f", esvSystem.getThetaMass()));
    logger.info(format("  Extended System Lambda Print Frequency: %d (fsec)", printEsvFrequency));
  }

  /**
   * Blocking molecular dynamics. When this method returns, the MD run is done.
   *
   * @param nSteps Number of MD steps
   * @param timeStep Time step in femtoseconds
   * @param loggingInterval Interval between printing/logging information in picoseconds.
   * @param trajectoryInterval Interval between adding a frame to the trajectory file in
   *     picoseconds.
   * @param temperature Temperature in Kelvins.
   * @param initVelocities Initialize new velocities from a Maxwell-Boltzmann distribution.
   * @param fileType XYZ or ARC to save to .arc, PDB for .pdb files
   * @param restartInterval Interval between writing new restart files in picoseconds.
   * @param dyn A {@link java.io.File} object to write the restart file to.
   */
  public void dynamic(final long nSteps, final double timeStep, final double loggingInterval,
      final double trajectoryInterval, final double temperature, final boolean initVelocities,
      String fileType, double restartInterval, final File dyn) {
    this.fileType = fileType;
    setRestartFrequency(restartInterval);
    dynamic(nSteps, timeStep, loggingInterval, trajectoryInterval, temperature, initVelocities, dyn);
  }

  /**
   * Blocking molecular dynamics. When this method returns, the MD run is done.
   *
   * @param nSteps Number of MD steps
   * @param timeStep Time step (fsec)
   * @param loggingInterval Interval between printing/logging information in picoseconds.
   * @param trajectoryInterval Interval between adding a frame to the trajectory file in
   *     picoseconds.
   * @param temperature Temperature in Kelvins.
   * @param initVelocities Initialize new velocities from a Maxwell-Boltzmann distribution.
   * @param dyn A {@link java.io.File} object to write the restart file to.
   */
  public void dynamic(final long nSteps, final double timeStep, final double loggingInterval,
      final double trajectoryInterval, final double temperature, final boolean initVelocities,
      final File dyn) {
    // Return if already running;
    // Could happen if two threads call dynamic on the same MolecularDynamics instance.
    if (!done) {
      logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
      return;
    }

    init(nSteps, timeStep, loggingInterval, trajectoryInterval, fileType, restartInterval,
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
    if (!verbosityLevel.isQuiet()) {
      logger.info(" Done with an MD round.");
    }
  }

  /**
   * Returns the MolecularAssembly array.
   *
   * @return A copy of the MolecularAssembly array.
   */
  public MolecularAssembly[] getAssemblyArray() {
    return molecularAssembly.clone();
  }

  /**
   * Returns the associated dynamics file.
   *
   * @return Dynamics restart File.
   */
  public File getDynFile() {
    return restartFile;
  }

  /**
   * Gets the kinetic energy at the start of the last dynamics run.
   *
   * @return Kinetic energy at the start of the run.
   */
  public double getInitialKineticEnergy() {
    return initialKinetic;
  }

  /**
   * Gets the potential energy at the start of the last dynamics run.
   *
   * @return potential energy at the start of the run.
   */
  public double getInitialPotentialEnergy() {
    return initialPotential;
  }

  /**
   * Gets the temperature at the start of the last dynamics run.
   *
   * @return temperature at the start of the run.
   */
  public double getInitialTemperature() {
    return initialTemp;
  }

  /**
   * Gets the total energy at the start of the last dynamics run.
   *
   * @return total energy at the start of the run.
   */
  public double getInitialTotalEnergy() {
    return initialTotal;
  }

  /**
   * getIntervalSteps.
   *
   * @return Always 1 for this implementation.
   */
  public int getIntervalSteps() {
    return 1;
  }

  /**
   * No-op; FFX does not need to occasionally return information from FFX.
   *
   * @param intervalSteps Ignored.
   */
  public void setIntervalSteps(int intervalSteps) {
    // Not meaningful for FFX MD.
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
   * Get the current temperature of the system
   *
   * @return currentTemperature
   */
  public double getTemperature() {
    return currentTemperature;
  }

  /**
   * Getter for the field <code>thermostat</code>.
   *
   * @return a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
   */
  public Thermostat getThermostat() {
    return thermostat;
  }

  /**
   * Setter for the field <code>thermostat</code>.
   *
   * @param thermostat a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
   */
  public void setThermostat(Thermostat thermostat) {
    this.thermostat = thermostat;
  }

  /**
   * getTimeStep.
   *
   * @return Timestep in picoseconds.
   */
  public double getTimeStep() {
    return dt;
  }

  /**
   * Sets the timestep and resets frequencies based on stored intervals.
   *
   * @param step Timestep in femtoseconds.
   */
  private void setTimeStep(double step) {
    dt = step * FSEC_TO_PSEC;
    // Reset frequencies to be consistent with new timestep.
    setRestartFrequency(restartInterval);
    setLoggingFrequency(logInterval);
    setTrajectoryFrequency(trajectoryInterval);
  }

  /**
   * Get the total system energy (kinetic plus potential).
   *
   * @return total energy.
   */
  public double getTotalEnergy() {
    return currentTotalEnergy;
  }

  public MDVerbosity getVerbosityLevel() {
    return verbosityLevel;
  }

  public void setVerbosityLevel(MDVerbosity level) {
    verbosityLevel = level;
    if (level == MDVerbosity.SILENT) {
      basicLogging = Level.FINE;
    } else {
      basicLogging = Level.INFO;
    }
  }

  /**
   * init
   *
   * @param nSteps Number of MD steps
   * @param timeStep Time step in femtoseconds
   * @param loggingInterval Interval between printing/logging information in picoseconds.
   * @param trajectoryInterval Interval between adding a frame to the trajectory file in
   *     picoseconds.
   * @param fileType XYZ or ARC to save to .arc, PDB for .pdb files
   * @param restartInterval Interval between writing new restart files in picoseconds.
   * @param temperature Temperature in Kelvins.
   * @param initVelocities Initialize new velocities from a Maxwell-Boltzmann distribution.
   * @param dyn A {@link java.io.File} object to write the restart file to.
   */
  public void init(final long nSteps, final double timeStep, final double loggingInterval,
      final double trajectoryInterval, final String fileType, final double restartInterval,
      final double temperature, final boolean initVelocities, final File dyn) {

    // Return if already running.
    if (!done) {
      logger.warning(
          " Programming error - attempt to modify parameters of a running MolecularDynamics instance.");
      return;
    }

    if (integrator instanceof Stochastic) {
      if (constantPressure) {
        logger.log(basicLogging, "\n Stochastic dynamics in the NPT ensemble");
      } else {
        logger.log(basicLogging, "\n Stochastic dynamics in the NVT ensemble");
      }
    } else if (!(thermostat instanceof Adiabatic)) {
      if (constantPressure) {
        logger.log(basicLogging, "\n Molecular dynamics in the NPT ensemble");
      } else {
        logger.log(basicLogging, "\n Molecular dynamics in the NVT ensemble");
      }
    } else {
      if (constantPressure) {
        logger.severe("\n NPT Molecular dynamics requires a thermostat");
      } else {
        logger.log(basicLogging, "\n Molecular dynamics in the NVE ensemble");
      }
    }

    this.nSteps = nSteps;
    totalSimTime = 0.0;

    // Convert the time step from femtoseconds to picoseconds.
    setTimeStep(timeStep);

    // Convert intervals in ps to frequencies in timesteps.
    setLoggingFrequency(loggingInterval);
    setTrajectoryFrequency(trajectoryInterval);
    setRestartFrequency(restartInterval);

    checkFrequency("Reporting           ", logFrequency);
    if (automaticWriteouts) {
      // Only sanity check these values if MD is doing this itself.
      checkFrequency("Trajectory Write-Out", trajectoryFrequency);
      checkFrequency("Restart             ", restartFrequency);
    }

    // Set snapshot file type.
    saveSnapshotAsPDB = true;
    if (fileType.equalsIgnoreCase("XYZ") || fileType.equalsIgnoreCase("ARC")) {
      saveSnapshotAsPDB = false;
    } else if (!fileType.equalsIgnoreCase("PDB")) {
      logger.warning("Snapshot file type unrecognized; saving snapshots as PDB.\n");
    }

    assemblyInfo();

    this.restartFile = (dyn == null) ? fallbackDynFile : dyn;
    loadRestart = restartFile.exists() && !initialized;

    if (dynFilter == null) {
      dynFilter = new DYNFilter(molecularAssembly[0].getName());
    }

    this.targetTemperature = temperature;
    this.initVelocities = initVelocities;
    done = false;

    if (loadRestart) {
      logger.info("  Continuing from " + restartFile.getAbsolutePath());
    }

    if (!verbosityLevel.isQuiet()) {
      logger.info(format("  Integrator:          %15s", integrator.toString()));
      logger.info(format("  Thermostat:          %15s", thermostat.toString()));
      logger.info(format("  Number of steps:     %8d", nSteps));
      logger.info(format("  Time step:           %8.3f (fsec)", timeStep));
      logger.info(format("  Print interval:      %8.3f (psec)", loggingInterval));
      logger.info(format("  Save interval:       %8.3f (psec)", trajectoryInterval));
      if (molecularAssembly.length > 1) {
        for (int i = 0; i < molecularAssembly.length; i++) {
          File archiveFile = molecularAssembly[i].getArchiveFile();
          logger.info(format("  Archive file %3d: %s", (i + 1), archiveFile.getAbsolutePath()));
        }
      } else {
        logger.info(format("  Archive file:     %s",
            molecularAssembly[0].getArchiveFile().getAbsolutePath()));
      }
      logger.info(format("  Restart file:     %s", restartFile.getAbsolutePath()));
    }
  }

  /**
   * A version of init with the original method header. Redirects to the new method with default
   * values for added parameters. Needed by (at least) ReplicaExchange, which calls this directly.
   *
   * @param nSteps Number of MD steps
   * @param timeStep Time step in femtoseconds
   * @param loggingInterval Interval between printing/logging information in picoseconds.
   * @param trajectoryInterval Interval between adding a frame to the trajectory file in
   *     picoseconds.
   * @param temperature Temperature in Kelvins.
   * @param initVelocities Initialize new velocities from a Maxwell-Boltzmann distribution.
   * @param dyn A {@link java.io.File} object to write the restart file to.
   */
  public void init(final long nSteps, final double timeStep, final double loggingInterval,
      final double trajectoryInterval, final double temperature, final boolean initVelocities,
      final File dyn) {
    init(nSteps, timeStep, loggingInterval, trajectoryInterval, "XYZ", restartInterval, temperature,
        initVelocities, dyn);
  }

  /** Reinitialize the MD engine after a chemical change. */
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
    thermostat.setNumberOfVariables(numberOfVariables, x, v, mass, potential.getVariableTypes(),
        true);
    integrator.setNumberOfVariables(numberOfVariables, x, v, a, aPrevious, mass);
  }

  /**
   * Causes this MolecularDynamics to take an additional set of timesteps.
   *
   * @param nSteps Number of steps to take
   * @param temperature Temperature of simulation
   */
  @Deprecated
  public void redynamic(final int nSteps, final double temperature) {
    if (!verbosityLevel.isQuiet()) {
      setVerbosityLevel(MDVerbosity.QUIET);
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

  /** {@inheritDoc} */
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
   * Enable or disable automatic writeout of trajectory snapshots and restart files.
   *
   * <p>Primarily intended for use with MC-OST, which handles the timing itself.
   *
   * @param autoWriteout Whether to automatically write out trajectory and restart files.
   */
  public void setAutomaticWriteouts(boolean autoWriteout) {
    this.automaticWriteouts = autoWriteout;
  }

  /**
   * Sets the "fallback" .dyn file to write to if none is passed to the dynamic method.
   *
   * @param fallback Fallback dynamics restart file.
   */
  public void setFallbackDynFile(File fallback) {
    this.fallbackDynFile = fallback;
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
   * Not meaningful for FFX dynamics (no need to obtain velocities/accelerations from a different
   * program, especially one running on a GPU). Is a no-op.
   *
   * @param obtainVA Not meaningful for this implementation.
   */
  public void setObtainVelAcc(boolean obtainVA) {
    // Not meaningful for FFX dynamics.
  }

  /**
   * Method to set the Restart Frequency.
   *
   * @param restartInterval Time in ps between writing restart files.
   * @throws java.lang.IllegalArgumentException If restart frequency is not a positive number
   */
  public void setRestartFrequency(double restartInterval) throws IllegalArgumentException {
    restartFrequency = intervalToFreq(restartInterval, "Restart interval");
    this.restartInterval = restartInterval;
  }

  /**
   * Set the archive file for each MolecularAssembly.
   *
   * @param archiveFiles
   */
  public void setArchiveFiles(File[] archiveFiles) {
    logger.info(" Setting archive files:\n " + Arrays.toString(archiveFiles));
    int n = molecularAssembly.length;
    assert archiveFiles.length == n;
    for (int i = 0; i < n; i++) {
      molecularAssembly[i].setArchiveFile(archiveFiles[i]);
    }
  }

  /** {@inheritDoc} */
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
   * Write restart and trajectory files if the provided step matches the frequency and that file type
   * is requested.
   *
   * @param step Step to write files (if any) for.
   * @param trySnapshot If false, do not write snapshot even if the timestep is correct.
   * @param tryRestart If false, do not write a restart file even if the timestep is correct.
   * @return EnumSet of actions taken by this method.
   */
  public EnumSet<MDWriteAction> writeFilesForStep(long step, boolean trySnapshot,
      boolean tryRestart) {
    return writeFilesForStep(step, trySnapshot, tryRestart, null);
  }

  /**
   * Write restart and trajectory files if the provided step matches the frequency and that file type
   * is requested.
   *
   * @param step Step to write files (if any) for.
   * @param trySnapshot If false, do not write snapshot even if the timestep is correct.
   * @param tryRestart If false, do not write a restart file even if the timestep is correct.
   * @param extraLines Additional lines to append into the comments section of the snapshot (or
   *     null).
   * @return EnumSet of actions taken by this method.
   */
  public EnumSet<MDWriteAction> writeFilesForStep(long step, boolean trySnapshot, boolean tryRestart,
      String[] extraLines) {
    List<String> linesList =
        (extraLines == null) ? new ArrayList<>() : new ArrayList<>(Arrays.asList(extraLines));

    if (potential instanceof LambdaInterface) {
      String lamString = format("Lambda: %.8f", ((LambdaInterface) potential).getLambda());
      linesList.add(lamString);
    }

    String tempString = format("Temp: %.2f", thermostat.getTargetTemperature());
    linesList.add(tempString);

    if (esvSystem != null) {
      String pHString = format("pH: %.2f", esvSystem.getConstantPh());
      linesList.add(pHString);
    }

    Comm world = Comm.world();
    if (world != null && world.size() > 1) {
      String rankString = format("Rank: %d", world.rank());
      linesList.add(rankString);
    }

    String[] allLines = new String[linesList.size()];
    allLines = linesList.toArray(allLines);

    EnumSet<MDWriteAction> written = EnumSet.noneOf(MDWriteAction.class);
    if (step != 0) {
      // Write out snapshots in selected format every saveSnapshotFrequency steps.
      if (trySnapshot && trajectoryFrequency > 0 && step % trajectoryFrequency == 0) {

        // Log stats.
        logger.info(format("\n Average Values for the Last %d Out of %d Dynamics Steps\n",
            trajectoryFrequency, step));
        logger.info(format("  Simulation Time  %16.4f Picosecond", step * dt));
        logger.info(
            format("  Total Energy     %16.4f Kcal/mole   (+/-%9.4f)", totalEnergyStats.getMean(),
                totalEnergyStats.getStandardDeviation()));
        logger.info(format("  Potential Energy %16.4f Kcal/mole   (+/-%9.4f)",
            potentialEnergyStats.getMean(), potentialEnergyStats.getStandardDeviation()));
        logger.info(
            format("  Kinetic Energy   %16.4f Kcal/mole   (+/-%9.4f)", kineticEnergyStats.getMean(),
                kineticEnergyStats.getStandardDeviation()));
        logger.info(
            format("  Temperature      %16.4f Kelvin      (+/-%9.4f)\n", temperatureStats.getMean(),
                temperatureStats.getStandardDeviation()));
        totalEnergyStats.reset();
        potentialEnergyStats.reset();
        kineticEnergyStats.reset();
        temperatureStats.reset();

        if (esvSystem != null) {
          for (Atom atom : molecularAssembly[0].getAtomList()) {
            int atomIndex = atom.getIndex() - 1;
            atom.setOccupancy(esvSystem.getTitrationLambda(atomIndex));
            atom.setTempFactor(esvSystem.getTautomerLambda(atomIndex));
          }
        }
        appendSnapshot(allLines);
        written.add(MDWriteAction.SNAPSHOT);
      }

      // Write out restart files every saveRestartFileFrequency steps.
      if (tryRestart && restartFrequency > 0 && step % restartFrequency == 0) {
        writeRestart();
        written.add(MDWriteAction.RESTART);
      }
    }
    return written;
  }

  /** Write out a restart file. */
  public void writeRestart() {
    potential.writeAdditionalRestartInfo(true);
    String dynName = FileUtils.relativePathTo(restartFile).toString();
    if (dynFilter.writeDYN(restartFile, molecularAssembly[0].getCrystal(), x, v, a, aPrevious)) {
      logger.log(basicLogging, " Wrote dynamics restart file to " + dynName);
    } else {
      logger.log(basicLogging, " Writing dynamics restart file to " + dynName + " failed");
    }
    if (esvSystem != null) {
      esvSystem.writeRestart();
      esvSystem.writeLambdaHistogram(false);
    }
  }

  /** assemblyInfo. */
  void assemblyInfo() {
    for (MolecularAssembly assembly : molecularAssembly) {
      CompositeConfiguration properties = assembly.getProperties();
      File file = assembly.getFile();
      String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
      File archiveFile = assembly.getArchiveFile();
      if (archiveFile == null) {
        archiveFile = new File(filename + ".arc");
        assembly.setArchiveFile(XYZFilter.version(archiveFile));
      }
    }

//      if (ainfo.pdbFile == null) {
//        String extName = FilenameUtils.getExtension(file.getName());
//        if (extName.toLowerCase().startsWith("pdb")) {
//          ainfo.pdbFile = file;
//        } else {
//          ainfo.pdbFile = new File(filename + ".pdb");
//        }
//      }
//      if (ainfo.xyzFilter == null && esvSystem == null) {
//        ainfo.xyzFilter = new XYZFilter(file, mola, mola.getForceField(), aprops);
//      }
//      if (ainfo.xphFilter == null) {
//        ainfo.xphFilter = new XPHFilter(file, mola, mola.getForceField(), aprops, esvSystem);
//      }
//      if (ainfo.pdbFilter == null) {
//        ainfo.pdbFilter = new PDBFilter(ainfo.pdbFile, mola, mola.getForceField(), aprops);
//      }
  }

  /**
   * Converts an interval in ps to a frequency in timesteps.
   *
   * @param interval Interval between events in ps
   * @param describe Description of the event.
   * @return Frequency of event in timesteps per event.
   * @throws IllegalArgumentException If interval is not a positive finite value.
   */
  private int intervalToFreq(double interval, String describe) throws IllegalArgumentException {
    if (!Double.isFinite(interval) || interval <= 0) {
      throw new IllegalArgumentException(
          format(" %s must be " + "positive finite value in ps, was %10.4g", describe, interval));
    }
    if (interval >= dt) {
      return (int) (interval / dt);
    } else {
      logger.warning(format(
          " Specified %s of %.6f ps < timestep %.6f ps; " + "interval is set to once per timestep!",
          describe, interval, dt));
      return 1;
    }
  }

  /**
   * Method to set the logging/printing frequency.
   *
   * @param logInterval Time in ps between logging information to the screen.
   */
  private void setLoggingFrequency(double logInterval) {
    logFrequency = intervalToFreq(logInterval, "Reporting (logging) interval");
    this.logInterval = logInterval;
  }

  /**
   * Method to set the frequency of appending snapshots to the trajectory file.
   *
   * @param snapshotInterval Time in ps between appending snapshots to the trajectory file.
   */
  private void setTrajectoryFrequency(double snapshotInterval) {
    trajectoryFrequency = intervalToFreq(snapshotInterval, "Trajectory writeout interval");
    this.trajectoryInterval = snapshotInterval;
  }

  /** Performs basic pre-MD operations such as loading the restart file. */
  void preRunOps() throws IllegalStateException {
    done = false;
    terminate = false;

    // Set the target temperature.
    thermostat.setTargetTemperature(targetTemperature);
    boolean quiet = verbosityLevel.isQuiet();
    thermostat.setQuiet(quiet);
    if (integrator instanceof Stochastic stochastic) {
      stochastic.setTemperature(targetTemperature);
    }

    // Set the step size.
    integrator.setTimeStep(dt);

    if (!initialized) {
      // Initialize from a restart file.
      if (loadRestart) {
        Crystal crystal = molecularAssembly[0].getCrystal();
        if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
          String message = " Could not load the restart file - dynamics terminated.";
          logger.log(Level.WARNING, message);
          done = true;
          throw new IllegalStateException(message);
        } else {
          molecularAssembly[0].setCrystal(crystal);
        }
      } else {
        // Initialize using current atomic coordinates.
        potential.getCoordinates(x);
        // Initialize atomic velocities from a Maxwell-Boltzmann distribution or set to 0.
        if (initVelocities) {
          thermostat.maxwell(targetTemperature);
        } else {
          fill(v, 0.0);
        }
      }
    } else {
      // If MD has already been run (i.e., Annealing or RepEx), then initialize velocities if
      // requested.
      if (initVelocities) {
        thermostat.maxwell(targetTemperature);
      }
    }
  }

  /** Initializes energy fields, esp. potential energy. */
  private void initializeEnergies() {
    // Compute the current potential energy.

    if (esvSystem != null && potential instanceof ForceFieldEnergyOpenMM) {
      currentPotentialEnergy = ((ForceFieldEnergyOpenMM) potential).energyAndGradientFFX(x,
          gradient);
    } else {
      currentPotentialEnergy = potential.energyAndGradient(x, gradient);
    }

    initialPotential = currentPotentialEnergy;

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
    thermostat.computeKineticEnergy();
    currentKineticEnergy = thermostat.getKineticEnergy();
    initialKinetic = currentKineticEnergy;
    currentTemperature = thermostat.getCurrentTemperature();
    initialTemp = currentTemperature;
    currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;
    initialTotal = currentTotalEnergy;

    temperatureStats.reset();
    potentialEnergyStats.reset();
    kineticEnergyStats.reset();
    totalEnergyStats.reset();
  }

  /** Pre-run operations (mostly logging) that require knowledge of system energy. */
  void postInitEnergies() {
    initialized = true;

    logger.log(basicLogging,
        format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp",
            "CPU"));
    logger.log(basicLogging,
        format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K",
            "sec"));
    logger.log(basicLogging,
        format("  %8s %12.4f %12.4f %12.4f %8.2f", "", currentKineticEnergy, currentPotentialEnergy,
            currentTotalEnergy, currentTemperature));

    // Store the initialized state.
    storeState();
  }

  /** Post-run cleanup operations. */
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
  }

  /**
   * Append a snapshot to the trajectory file.
   *
   * @param extraLines Strings of meta-data to include.
   */
  protected void appendSnapshot(String[] extraLines) {
    // Loop over all molecular assemblies.
    for (MolecularAssembly assembly : molecularAssembly) {
      File archiveFile = assembly.getArchiveFile();
      ForceField forceField = assembly.getForceField();
      CompositeConfiguration properties = assembly.getProperties();

      // Save as an ARC file.
      if (archiveFile != null && !saveSnapshotAsPDB) {
        String aiName = FileUtils.relativePathTo(archiveFile).toString();
        if (esvSystem == null) {
          XYZFilter xyzFilter = new XYZFilter(archiveFile, assembly, forceField, properties);
          if (xyzFilter.writeFile(archiveFile, true, extraLines)) {
            logger.log(basicLogging, format(" Appended to archive %s", aiName));
          } else {
            logger.warning(format(" Appending to archive %s failed.", aiName));
          }
        } else {
          XPHFilter xphFilter = new XPHFilter(archiveFile, assembly, forceField, properties,
              esvSystem);
          if (xphFilter.writeFile(archiveFile, true, extraLines)) {
            logger.log(basicLogging, format(" Appended to XPH archive %s", aiName));
          } else {
            logger.warning(format(" Appending to XPH archive %s failed.", aiName));
          }
        }
      } else if (saveSnapshotAsPDB) {
        File file = assembly.getFile();
        String extName = FilenameUtils.getExtension(file.getName());
        File pdbFile;
        if (extName.toLowerCase().startsWith("pdb")) {
          pdbFile = file;
        } else {
          String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
          pdbFile = new File(filename + ".pdb");
        }
        String aiName = FileUtils.relativePathTo(pdbFile).toString();
        PDBFilter pdbFilter = new PDBFilter(pdbFile, assembly, forceField, properties);
        if (pdbFilter.writeFile(pdbFile, true, extraLines)) {
          logger.log(basicLogging, format(" Appended to PDB file %s", aiName));
        } else {
          logger.warning(format(" Appending to PDB file to %s failed.", aiName));
        }
      }
    }
  }

  /**
   * Checks if thermodynamics must be logged. If logged, current time is returned, else the time
   * passed in is returned.
   *
   * @param step Timestep to possibly log thermodynamics for.
   * @param time Clock time (in nsec) thermodynamics was last logged.
   * @return Either current time (if logged), else the time variable passed in.
   */
  protected long logThermoForTime(long step, long time) {
    if (step % logFrequency == 0) {
      return logThermodynamics(time);
    } else {
      return time;
    }
  }

  /**
   * Log thermodynamics to the screen.
   *
   * @param time Clock time (in nsec) since this was last called.
   * @return Clock time at the end of this call.
   */
  private long logThermodynamics(long time) {
    time = System.nanoTime() - time;
    logger.log(basicLogging,
        format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.3f", totalSimTime, currentKineticEnergy,
            currentPotentialEnergy, currentTotalEnergy, currentTemperature, time * NS2SEC));
    return System.nanoTime();
  }

  /** Main loop of the run method. */
  private void mainLoop() {
    // Integrate Newton's equations of motion for the requested number of steps,
    // unless early termination is requested.
    long time = System.nanoTime();
    for (long step = 1; step <= nSteps; step++) {
      if (step > 1) {
        List<Constraint> constraints = potential.getConstraints();
        // TODO: Replace magic numbers with named constants.
        long constraintFails = constraints.stream()
            .filter((Constraint c) -> !c.constraintSatisfied(x, v, 1E-7, 1E-7)).count();
        if (constraintFails > 0) {
          logger.info(format(" %d constraint failures in step %d", constraintFails, step));
        }
      }

      // Do the half-step thermostat operation.
      thermostat.halfStep(dt);

      // Do the half-step integration operation.
      integrator.preForce(potential);
      if (esvSystem != null) {
        esvIntegrator.preForce(potential);
        //preForce processes theta values after
        esvSystem.preForce();
      }

      // Compute the potential energy and gradients.

      if (esvSystem != null && potential instanceof ForceFieldEnergyOpenMM) {
        currentPotentialEnergy = ((ForceFieldEnergyOpenMM) potential).energyAndGradientFFX(x,
            gradient);
      } else {
        currentPotentialEnergy = potential.energyAndGradient(x, gradient);
      }

      // Add the potential energy of the slow degrees of freedom.
      if (integrator instanceof Respa r) {
        currentPotentialEnergy += r.getHalfStepEnergy();
      }

      // Do the full-step integration operation.
      integrator.postForce(gradient);
      if (esvSystem != null) {
        double[] dEdL = esvSystem.postForce();
        esvIntegrator.postForce(dEdL);
      }

      // Compute the full-step kinetic energy.
      thermostat.computeKineticEnergy();

      // Do the full-step thermostat operation.
      thermostat.fullStep(dt);

      // Recompute the kinetic energy after the full-step thermostat operation.
      thermostat.computeKineticEnergy();
      if (esvSystem != null) {
        //Adiabatic thermostat does nothing at half step.
        esvThermostat.computeKineticEnergy();
      }

      // Remove center of mass motion every ~100 steps.
      int removeCOMMotionFrequency = 100;
      if (thermostat.getRemoveCenterOfMassMotion() && step % removeCOMMotionFrequency == 0) {
        thermostat.centerOfMassMotion(true, false);
      }

      // Collect current kinetic energy, temperature, and total energy.
      currentKineticEnergy = thermostat.getKineticEnergy();
      if (esvSystem != null) {
        currentKineticEnergy += esvThermostat.getKineticEnergy();
      }
      currentTemperature = thermostat.getCurrentTemperature();
      currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

      // Collect running statistics.
      temperatureStats.addValue(currentTemperature);
      potentialEnergyStats.addValue(currentPotentialEnergy);
      kineticEnergyStats.addValue(currentKineticEnergy);
      totalEnergyStats.addValue(currentTotalEnergy);

      // Update atomic velocity, acceleration and previous acceleration.
      potential.setVelocity(v);
      potential.setAcceleration(a);
      potential.setPreviousAcceleration(aPrevious);

      // Log the current state every printFrequency steps.
      totalSimTime += dt;
      time = logThermoForTime(step, time);
      if (step % printEsvFrequency == 0 && esvSystem != null) {
        logger.log(basicLogging, format(" %s", esvSystem.getLambdaList()));
      }

      if (automaticWriteouts) {
        writeFilesForStep(step, true, true);
      }

      // Notify the algorithmListeners.
      if (algorithmListener != null && step % logFrequency == 0) {
        for (MolecularAssembly assembly : molecularAssembly) {
          algorithmListener.algorithmUpdate(assembly);
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
   * Perform a sanity check on a frequency to ensure it's not longer than total runtime. Currently,
   * the setter parameter is ignored due to issues with our test suite.
   *
   * @param describe Description of the frequency.
   * @param frequency Frequency in timesteps.
   */
  private void checkFrequency(String describe, int frequency) {
    if (frequency > nSteps) {
      logger.fine(
          format(" Specified %s frequency of %d is greater than the number of steps %d", describe,
              frequency, nSteps));
    }
  }

  /**
   * Set the coordinates.
   *
   * @param coords The coordinates to copy into MD coordinates array.
   */
  public void setCoordinates(double[] coords) {
    if (coords.length == x.length) {
      System.arraycopy(coords, 0, x, 0, x.length);
    }
  }

  /**
   * Get the coordinates.
   *
   * @return A copy of the current coordinates are returned.
   */
  public double[] getCoordinates() {
    return Arrays.copyOf(x, x.length);
  }

  /**
   * Store the current state of the molecular dynamics simulation in a DynamicsState record.
   */
  public void storeState() {
    mdState = new MDState(x, v, a, aPrevious, mass, gradient, currentKineticEnergy,
        currentPotentialEnergy, currentTotalEnergy, currentTemperature);
  }

  /**
   * revertState.
   *
   * @throws java.lang.Exception if any.
   */
  public void revertState() throws Exception {
    if (mdState == null) {
      throw new Exception();
    }
    revertState(mdState);
  }

  /**
   * Revert the state of the MolecularDynamics instance to the provided DynamicsState.
   *
   * @param state The DynamicsState to revert to.
   */
  private void revertState(MDState state) {
    currentPotentialEnergy = state.potentialEnergy();
    currentKineticEnergy = state.kineticEnergy();
    currentTotalEnergy = state.totalEnergy();
    currentTemperature = state.temperature();
    arraycopy(state.x(), 0, x, 0, x.length);
    arraycopy(state.v(), 0, v, 0, v.length);
    arraycopy(state.a(), 0, a, 0, a.length);
    arraycopy(state.aPrevious(), 0, aPrevious, 0, aPrevious.length);
    arraycopy(state.mass(), 0, mass, 0, mass.length);
    arraycopy(state.gradient(), 0, gradient, 0, gradient.length);

    potential.setVelocity(v);
    potential.setAcceleration(a);
    potential.setPreviousAcceleration(aPrevious);

    // To Do -- move these methods into the Potential interface.
    Atom[] atoms = molecularAssembly[0].getActiveAtomArray();
    if (atoms.length * 3 == numberOfVariables) {
      int index = 0;
      for (Atom atom : atoms) {
        atom.moveTo(x[index], x[index + 1], x[index + 2]);
        atom.setXYZGradient(gradient[index], gradient[index + 1], gradient[index + 2]);
        index += 3;
      }
    }
  }

}
