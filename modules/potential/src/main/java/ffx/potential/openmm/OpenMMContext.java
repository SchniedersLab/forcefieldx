// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.potential.openmm;

import edu.rit.pj.Comm;
import edu.uiowa.jopenmm.OpenMMLibrary;
import edu.uiowa.jopenmm.OpenMMUtils;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.openmm.Context;
import ffx.openmm.Integrator;
import ffx.openmm.MinimizationReporter;
import ffx.openmm.Platform;
import ffx.openmm.State;
import ffx.openmm.StringArray;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalEnergyMinimizer_minimize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static ffx.openmm.Platform.getNumPlatforms;
import static ffx.openmm.Platform.getOpenMMVersion;
import static ffx.openmm.Platform.getPluginLoadFailures;
import static ffx.openmm.Platform.loadPluginsFromDirectory;
import static ffx.potential.ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;
import static ffx.potential.Platform.OMM;
import static ffx.potential.Platform.OMM_CUDA;
import static ffx.potential.Platform.OMM_OPENCL;
import static ffx.potential.openmm.OpenMMEnergy.getDefaultDevice;
import static ffx.potential.openmm.OpenMMIntegrator.createIntegrator;
import static java.lang.String.format;

/**
 * Creates and manage an OpenMM Context.
 *
 * <p>A Context stores the complete state of a simulation. More specifically, it includes: The
 * current time The position of each particle The velocity of each particle The values of
 * configurable parameters defined by Force objects in the System
 *
 * <p>You can retrieve a snapshot of the current state at any time by calling getState(). This
 * allows you to record the state of the simulation at various points, either for analysis or for
 * checkpointing. getState() can also be used to retrieve the current forces on each particle and
 * the current energy of the System.
 */
public class OpenMMContext extends Context {

  private static final Logger logger = Logger.getLogger(OpenMMContext.class.getName());

  /**
   * OpenMM Platform.
   */
  private final Platform openMMPlatform;
  /**
   * OpenMM System.
   */
  private final OpenMMSystem openMMSystem;
  /**
   * Instance of the OpenMM Integrator class.
   */
  private Integrator openMMIntegrator;
  /**
   * Integrator string (default = VERLET).
   */
  private String integratorName = "VERLET";
  /**
   * Time step (default = 0.001 psec).
   */
  private double timeStep = 0.001;
  /**
   * Temperature (default = 298.15).
   */
  private double temperature = 298.15;
  /**
   *
   */
  private final int enforcePBC;
  /**
   * Array of atoms.
   */
  private final Atom[] atoms;

  /**
   * Create an OpenMM Context.
   *
   * @param platform     OpenMM Platform.
   * @param openMMSystem OpenMM System.
   * @param atoms        Array of atoms.
   */
  public OpenMMContext(Platform platform, OpenMMSystem openMMSystem, Atom[] atoms) {
    this.openMMPlatform = platform;
    this.openMMSystem = openMMSystem;

    ForceField forceField = openMMSystem.getForceField();
    boolean aperiodic = openMMSystem.getCrystal().aperiodic();
    boolean pbcEnforced = forceField.getBoolean("ENFORCE_PBC", !aperiodic);
    enforcePBC = pbcEnforced ? OpenMM_True : OpenMM_False;

    this.atoms = atoms;
    update(integratorName, timeStep, temperature, true);
  }

  /**
   * Update the Context in which to run a simulation.
   *
   * @param integratorName Requested integrator.
   * @param timeStep       Time step (psec).
   * @param temperature    Temperature (K).
   * @param forceCreation  Force creation of a new context, even if the current one matches.
   */
  public void update(String integratorName, double timeStep, double temperature, boolean forceCreation) {
    // Check if the current context is consistent with the requested context.
    if (hasContextPointer() && !forceCreation) {
      if (this.temperature == temperature && this.timeStep == timeStep
          && this.integratorName.equalsIgnoreCase(integratorName)) {
        // All requested features agree.
        return;
      }
    }

    this.integratorName = integratorName;
    this.timeStep = timeStep;
    this.temperature = temperature;

    logger.info("\n Updating OpenMM Context");

    // Update the integrator.
    if (openMMIntegrator != null) {
      openMMIntegrator.destroy();
    }
    openMMIntegrator = createIntegrator(integratorName, timeStep, temperature, openMMSystem);

    // Set lambda to 1.0 when creating a context to avoid OpenMM compiling out any terms.
    // TODO: Test on a fixed charge system.
    // double currentLambda = openMMEnergy.getLambda();
    // if (openMMEnergy.getLambdaTerm()) {
    //  openMMEnergy.setLambda(1.0);
    // }

    // Update the context.
    updateContext(openMMSystem, openMMIntegrator, openMMPlatform);

    // Revert to the current lambda value.
    // if (openMMEnergy.getLambdaTerm()) {
    //   openMMEnergy.setLambda(currentLambda);
    // }

    // Get initial positions and velocities for active atoms.
    int nVar = openMMSystem.getNumberOfVariables();
    double[] x = new double[nVar];
    double[] v = new double[nVar];
    double[] vel3 = new double[3];
    int index = 0;
    for (Atom a : atoms) {
      if (a.isActive()) {
        a.getVelocity(vel3);
        // X-axis
        x[index] = a.getX();
        v[index++] = vel3[0];
        // Y-axis
        x[index] = a.getY();
        v[index++] = vel3[1];
        // Z-axis
        x[index] = a.getZ();
        v[index++] = vel3[2];
      }
    }

    // Load the current periodic box vectors.
    Crystal crystal = openMMSystem.getCrystal();
    setPeriodicBoxVectors(crystal);

    // Load current atomic positions.
    setPositions(x);

    // Load current velocities.
    setVelocities(v);

    // Apply constraints starting from current atomic positions.
    applyConstraints(DEFAULT_CONSTRAINT_TOLERANCE);

    // Application of constraints can change coordinates and velocities.
    // Retrieve them for consistency.
    OpenMMState openMMState = getOpenMMState(OpenMM_State_Positions | OpenMM_State_Velocities);
    openMMState.getPositions(x);
    openMMState.getVelocities(v);
    openMMState.destroy();
  }

  /**
   * Update the Context if necessary.
   */
  public void update() {
    if (!hasContextPointer()) {
      logger.info(" Delayed creation of OpenMM Context.");
      update(integratorName, timeStep, temperature, true);
    }
  }

  /**
   * Get an OpenMM State from the Context.
   *
   * @param mask A mask specifying which information to retrieve.
   * @return State pointer.
   */
  public OpenMMState getOpenMMState(int mask) {
    State state = getState(mask, enforcePBC);
    return new OpenMMState(state.getPointer(), atoms, openMMSystem.getNumberOfVariables());
  }

  /**
   * Use the Context / Integrator combination to take the requested number of steps.
   *
   * @param numSteps Number of steps to take.
   */
  public void integrate(int numSteps) {
    openMMIntegrator.step(numSteps);
  }

  /**
   * Use the Context to optimize the system to the requested tolerance.
   *
   * @param eps           Convergence criteria (kcal/mole/A).
   * @param maxIterations Maximum number of iterations.
   */
  public void optimize(double eps, int maxIterations) {
    // The "report" method of MinimizationReporter cannot be overridden, so the reporter does nothing.
    MinimizationReporter reporter = new MinimizationReporter();
    OpenMM_LocalEnergyMinimizer_minimize(getPointer(), eps / (OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ),
        maxIterations, reporter.getPointer());
    reporter.destroy();
  }

  /**
   * The array x contains atomic coordinates only for active atoms.
   *
   * @param x Atomic coordinate array for only active atoms.
   */
  public void setPositions(double[] x) {
    double[] allPositions = new double[3 * atoms.length];
    double[] d = new double[3];
    int index = 0;
    for (Atom a : atoms) {
      if (a.isActive()) {
        a.moveTo(x[index], x[index + 1], x[index + 2]);
      }
      a.getXYZ(d);
      allPositions[index] = d[0] * OpenMM_NmPerAngstrom;
      allPositions[index + 1] = d[1] * OpenMM_NmPerAngstrom;
      allPositions[index + 2] = d[2] * OpenMM_NmPerAngstrom;
      index += 3;
    }
    super.setPositions(allPositions);
  }

  /**
   * The array v contains velocity values for active atomic coordinates.
   *
   * @param v Velocity array for active atoms.
   */
  public void setVelocities(double[] v) {
    double[] allVelocities = new double[3 * atoms.length];
    double[] velocity = new double[3];
    int index = 0;
    for (Atom a : atoms) {
      if (a.isActive()) {
        a.setVelocity(v[index], v[index + 1], v[index + 2]);
      } else {
        // OpenMM requires velocities for even "inactive" atoms with mass of zero.
        a.setVelocity(0.0, 0.0, 0.0);
      }
      a.getVelocity(velocity);
      allVelocities[index] = velocity[0] * OpenMM_NmPerAngstrom;
      allVelocities[index + 1] = velocity[1] * OpenMM_NmPerAngstrom;
      allVelocities[index + 2] = velocity[2] * OpenMM_NmPerAngstrom;
      index += 3;
    }
    super.setVelocities(allVelocities);
  }

  /**
   * Set the periodic box vectors for a context based on the crystal instance.
   *
   * @param crystal The crystal instance.
   */
  public void setPeriodicBoxVectors(Crystal crystal) {
    if (!crystal.aperiodic()) {
      OpenMM_Vec3 a = new OpenMM_Vec3();
      OpenMM_Vec3 b = new OpenMM_Vec3();
      OpenMM_Vec3 c = new OpenMM_Vec3();
      double[][] Ai = crystal.Ai;
      a.x = Ai[0][0] * OpenMM_NmPerAngstrom;
      a.y = Ai[0][1] * OpenMM_NmPerAngstrom;
      a.z = Ai[0][2] * OpenMM_NmPerAngstrom;
      b.x = Ai[1][0] * OpenMM_NmPerAngstrom;
      b.y = Ai[1][1] * OpenMM_NmPerAngstrom;
      b.z = Ai[1][2] * OpenMM_NmPerAngstrom;
      c.x = Ai[2][0] * OpenMM_NmPerAngstrom;
      c.y = Ai[2][1] * OpenMM_NmPerAngstrom;
      c.z = Ai[2][2] * OpenMM_NmPerAngstrom;
      setPeriodicBoxVectors(a, b, c);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return format(
        " OpenMM context with integrator %s, timestep %9.3g fsec, temperature %9.3g K",
        integratorName, timeStep, temperature);
  }

  /**
   * Load an OpenMM Platform
   *
   * @param requestedPlatform the requested OpenMM platform.
   * @param forceField        the ForceField to query for platform flags.
   * @return the loaded Platform.
   */
  public static Platform loadPlatform(ffx.potential.Platform requestedPlatform, ForceField forceField) {

    OpenMMUtils.init();

    logger.log(Level.INFO, " Loaded from:\n {0}", OpenMMLibrary.JNA_NATIVE_LIB.toString());

    // Print out the OpenMM Version.
    logger.log(Level.INFO, " Version: {0}", getOpenMMVersion());

    // Print out the OpenMM lib directory.
    String libDirectory = OpenMMUtils.getLibDirectory();
    logger.log(Level.FINE, " Lib Directory:       {0}", libDirectory);
    // Load platforms and print out their names.
    StringArray libs = loadPluginsFromDirectory(libDirectory);
    int numLibs = libs.getSize();
    logger.log(Level.FINE, " Number of libraries: {0}", numLibs);
    for (int i = 0; i < numLibs; i++) {
      logger.log(Level.FINE, "  Library: {0}", libs.get(i));
    }
    libs.destroy();

    // Print out the OpenMM plugin directory.
    String pluginDirectory = OpenMMUtils.getPluginDirectory();
    logger.log(Level.INFO, "\n Plugin Directory:  {0}", pluginDirectory);
    // Load plugins and print out their names.
    StringArray plugins = loadPluginsFromDirectory(pluginDirectory);
    int numPlugins = plugins.getSize();
    logger.log(Level.INFO, " Number of Plugins: {0}", numPlugins);
    boolean cuda = false;
    boolean opencl = false;
    for (int i = 0; i < numPlugins; i++) {
      String pluginString = plugins.get(i);
      logger.log(Level.INFO, "  Plugin: {0}", pluginString);
      if (pluginString != null) {
        pluginString = pluginString.toUpperCase();
        boolean amoebaCudaAvailable = pluginString.contains("AMOEBACUDA");
        if (amoebaCudaAvailable) {
          cuda = true;
        }
        boolean amoebaOpenCLAvailable = pluginString.contains("AMOEBAOPENCL");
        if (amoebaOpenCLAvailable) {
          opencl = true;
        }
      }
    }
    plugins.destroy();

    int numPlatforms = getNumPlatforms();
    logger.log(Level.INFO, " Number of Platforms: {0}", numPlatforms);

    if (requestedPlatform == OMM_CUDA && !cuda) {
      logger.severe(" The OMM_CUDA platform was requested, but is not available.");
    }

    if (requestedPlatform == ffx.potential.Platform.OMM_OPENCL && !opencl) {
      logger.severe(" The OMM_OPENCL platform was requested, but is not available.");
    }

    // Extra logging to print out plugins that failed to load.
    if (logger.isLoggable(Level.FINE)) {
      StringArray pluginFailures = getPluginLoadFailures();
      int numFailures = pluginFailures.getSize();
      for (int i = 0; i < numFailures; i++) {
        logger.log(Level.FINE, " Plugin load failure: {0}", pluginFailures.get(i));
      }
      pluginFailures.destroy();
    }

    String defaultPrecision = "mixed";
    String precision = forceField.getString("PRECISION", defaultPrecision).toLowerCase();
    precision = precision.replace("-precision", "");
    switch (precision) {
      case "double", "mixed", "single" -> logger.info(format(" Precision level: %s", precision));
      default -> {
        logger.info(format(" Could not interpret precision level %s, defaulting to %s", precision, defaultPrecision));
        precision = defaultPrecision;
      }
    }


    Platform openMMPlatform;
    if (cuda && (requestedPlatform == OMM_CUDA || requestedPlatform == OMM)) {
      // CUDA
      int defaultDevice = getDefaultDevice(forceField.getProperties());
      openMMPlatform = new Platform("CUDA");
      // CUDA_DEVICE is deprecated; use DeviceIndex.
      int deviceID = forceField.getInteger("CUDA_DEVICE", defaultDevice);
      deviceID = forceField.getInteger("DeviceIndex", deviceID);
      String deviceIDString = Integer.toString(deviceID);
      openMMPlatform.setPropertyDefaultValue("DeviceIndex", deviceIDString);
      openMMPlatform.setPropertyDefaultValue("Precision", precision);
      String name = openMMPlatform.getName();
      logger.info(format(" Platform: %s (Device Index %d)", name, deviceID));
    } else if (opencl && (requestedPlatform == OMM_OPENCL || requestedPlatform == OMM)) {
      // OpenCL
      int defaultDevice = getDefaultDevice(forceField.getProperties());
      openMMPlatform = new Platform("OpenCL");
      int deviceID = forceField.getInteger("DeviceIndex", defaultDevice);
      String deviceIDString = Integer.toString(deviceID);
      openMMPlatform.setPropertyDefaultValue("DeviceIndex", deviceIDString);
      int openCLPlatformIndex = forceField.getInteger("OpenCLPlatformIndex", 0);
      String openCLPlatformIndexString = Integer.toString(openCLPlatformIndex);
      openMMPlatform.setPropertyDefaultValue("DeviceIndex", deviceIDString);
      openMMPlatform.setPropertyDefaultValue("OpenCLPlatformIndex", openCLPlatformIndexString);
      openMMPlatform.setPropertyDefaultValue("Precision", precision);
      String name = openMMPlatform.getName();
      logger.info(format(" Platform: %s (Platform Index %d, Device Index %d)",
          name, openCLPlatformIndex, deviceID));
    } else {
      // Reference
      openMMPlatform = new Platform("Reference");
      String name = openMMPlatform.getName();
      logger.info(format(" Platform: %s", name));
    }

    try {
      Comm world = Comm.world();
      if (world != null) {
        logger.info(format(" Running on host %s, rank %d", world.host(), world.rank()));
      }
    } catch (IllegalStateException illegalStateException) {
      logger.fine(" Could not find the world communicator!");
    }

    return openMMPlatform;
  }

  /**
   * Free OpenMM memory for the current Context and Integrator.
   */
  public void free() {
    if (openMMIntegrator != null) {
      openMMIntegrator.destroy();
    }
    destroy();
  }
}
