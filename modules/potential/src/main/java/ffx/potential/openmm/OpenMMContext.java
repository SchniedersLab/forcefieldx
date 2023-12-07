package ffx.potential.openmm;

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;
import edu.rit.pj.Comm;
import edu.uiowa.jopenmm.OpenMMLibrary;
import edu.uiowa.jopenmm.OpenMMUtils;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergy.Platform;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_applyConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_reinitialize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalEnergyMinimizer_minimize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getNumPlatforms;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getOpenMMVersion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPlatformByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPluginLoadFailures;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_loadPluginsFromDirectory;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_getSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_destroy;
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
public class OpenMMContext {

  private static final Logger logger = Logger.getLogger(OpenMMContext.class.getName());

  /**
   * ForceFieldEnergyOpenMM instance.
   */
  private final OpenMMEnergy openMMEnergy;
  /**
   * Requested Platform (i.e. Java or an OpenMM platform).
   */
  private final Platform platform;
  /**
   * Instance of the OpenMM Integrator class.
   */
  private final OpenMMIntegrator openMMIntegrator;
  /**
   * Constraint tolerance as a fraction of the constrained bond length.
   */
  private final double constraintTolerance = ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;
  /**
   * OpenMM Context pointer.
   */
  private PointerByReference contextPointer = null;
  /**
   * Integrator string (default = VERLET).
   */
  private String integratorString = "VERLET";
  /**
   * Time step (default = 0.001 psec).
   */
  private double timeStep = 0.001;
  /**
   * OpenMM Platform pointer.
   */
  private PointerByReference platformPointer = null;
  /**
   * Temperature (default = 298.15).
   */
  private double temperature = 298.15;
  /**
   *
   */
  private int enforcePBC = 0;
  /**
   * Array of atoms.
   */
  private final Atom[] atoms;

  /**
   * Create an OpenMM Context.
   *
   * @param requestedPlatform      Platform requested.
   * @param atoms                  Array of atoms.
   * @param enforcePBC             Enforce periodic boundary conditions.
   * @param openMMEnergy ForceFieldEnergyOpenMM instance.
   */
  public OpenMMContext(ForceField forceField, Platform requestedPlatform, Atom[] atoms, int enforcePBC,
                       OpenMMEnergy openMMEnergy) {
    loadPlatform(requestedPlatform);
    platform = requestedPlatform;
    this.enforcePBC = enforcePBC;
    this.atoms = atoms;
    this.openMMEnergy = openMMEnergy;
    openMMIntegrator = new OpenMMIntegrator(forceField, constraintTolerance);
  }

  /**
   * Get a Pointer to the OpenMM Context. A Context is optionally created if none has been instantiated yet.
   *
   * @param createContext If true, create a Context if none has been instantiated yet.
   * @return Context pointer.
   */
  public PointerByReference getContextPointer(boolean createContext) {
    if (contextPointer == null && createContext) {
      openMMEnergy.createContext(integratorString, timeStep, temperature, true);
    }
    return contextPointer;
  }

  /**
   * Get a Pointer to the OpenMMContext.
   * Null is returned if no OpenMMContext has been instantiated yet.
   *
   * @return Context pointer.
   */
  public PointerByReference getContextPointer() {
    return contextPointer;
  }

  /**
   * Check if the OpenMMContext has a Context pointer.
   *
   * @return True if the Context pointer is not null.
   */
  public boolean hasContextPointer() {
    return contextPointer != null;
  }

  /**
   * Get an OpenMM State from the Context.
   *
   * @param mask A mask specifying which information to retrieve.
   * @return State pointer.
   */
  public PointerByReference getState(int mask) {
    return OpenMM_Context_getState(contextPointer, mask, enforcePBC);
  }

  public ForceFieldEnergy.Platform getPlatform() {
    return platform;
  }

  /**
   * Get a Pointer to the OpenMM Integrator.
   *
   * @return Integrator pointer.
   */
  public PointerByReference getIntegrator() {
    return openMMIntegrator.getIntegratorPointer();
  }

  /**
   * Use the Context / Integrator combination to take the requested number of steps.
   *
   * @param numSteps Number of steps to take.
   */
  public void integrate(int numSteps) {
    OpenMM_Integrator_step(openMMIntegrator.getIntegratorPointer(), numSteps);
  }

  /**
   * Use the Context to optimize the system to the requested tolerance.
   *
   * @param eps           Convergence criteria (kcal/mole/A).
   * @param maxIterations Maximum number of iterations.
   */
  public void optimize(double eps, int maxIterations) {
    OpenMM_LocalEnergyMinimizer_minimize(contextPointer,
        eps / (OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ), maxIterations);
  }

  /**
   * The array x contains atomic coordinates only for active atoms.
   *
   * @param x Atomic coordinate array for only active atoms.
   */
  public void setOpenMMPositions(double[] x) {
    PointerByReference positions = OpenMM_Vec3Array_create(0);
    OpenMM_Vec3.ByValue coords = new OpenMM_Vec3.ByValue();
    double[] d = new double[3];
    int index = 0;
    for (Atom a : atoms) {
      if (a.isActive()) {
        a.moveTo(x[index++], x[index++], x[index++]);
        a.getXYZ(d);
        coords.x = d[0] * OpenMM_NmPerAngstrom;
        coords.y = d[1] * OpenMM_NmPerAngstrom;
        coords.z = d[2] * OpenMM_NmPerAngstrom;
        OpenMM_Vec3Array_append(positions, coords);
      } else {
        // OpenMM requires coordinates for even "inactive" atoms with mass of zero.
        coords.x = a.getX() * OpenMM_NmPerAngstrom;
        coords.y = a.getY() * OpenMM_NmPerAngstrom;
        coords.z = a.getZ() * OpenMM_NmPerAngstrom;
        OpenMM_Vec3Array_append(positions, coords);
      }
    }
    OpenMM_Context_setPositions(contextPointer, positions);
    logger.finer(" Free OpenMM positions.");
    OpenMM_Vec3Array_destroy(positions);
    logger.finer(" Free OpenMM positions completed.");
  }

  /**
   * The array v contains velocity values for active atomic coordinates.
   *
   * @param v Velocity array for active atoms.
   */
  public void setOpenMMVelocities(double[] v) {
    PointerByReference velocities = OpenMM_Vec3Array_create(0);
    OpenMM_Vec3.ByValue vel = new OpenMM_Vec3.ByValue();
    int index = 0;
    double[] velocity = new double[3];
    for (Atom a : atoms) {
      if (a.isActive()) {
        a.setVelocity(v[index++], v[index++], v[index++]);
        a.getVelocity(velocity);
        vel.x = velocity[0] * OpenMM_NmPerAngstrom;
        vel.y = velocity[1] * OpenMM_NmPerAngstrom;
        vel.z = velocity[2] * OpenMM_NmPerAngstrom;
        OpenMM_Vec3Array_append(velocities, vel);
      } else {
        // OpenMM requires velocities for even "inactive" atoms with mass of zero.
        a.setVelocity(0.0, 0.0, 0.0);
        vel.x = 0.0;
        vel.y = 0.0;
        vel.z = 0.0;
        OpenMM_Vec3Array_append(velocities, vel);
      }
    }
    OpenMM_Context_setVelocities(contextPointer, velocities);
    logger.finer(" Free OpenMM velocities.");
    OpenMM_Vec3Array_destroy(velocities);
    logger.finer(" Free OpenMM velocities completed.");
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
      OpenMM_Context_setPeriodicBoxVectors(contextPointer, a, b, c);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return format(
        " OpenMM context with integrator %s, timestep %9.3g fsec, temperature %9.3g K, constraintTolerance %9.3g",
        integratorString, timeStep, temperature, constraintTolerance);
  }

  /**
   * Free OpenMM memory for the current Context and Integrator.
   */
  public void free() {
    if (openMMIntegrator != null) {
      openMMIntegrator.free();
    }
    if (contextPointer != null) {
      logger.fine(" Free OpenMM Context.");
      OpenMM_Context_destroy(contextPointer);
      logger.fine(" Free OpenMM Context completed.");
      contextPointer = null;
    }
  }

  /**
   * Construct a new Context in which to run a simulation.
   *
   * @param integratorString       Requested integrator.
   * @param timeStep               Time step (psec).
   * @param temperature            Temperature (K).
   * @param forceCreation          Force creation of a new context, even if the current one matches.
   * @param openMMEnergy ForceFieldEnergyOpenMM instance.
   * @return Pointer to the created OpenMM context.
   */
  public OpenMMContext create(String integratorString, double timeStep, double temperature,
                              boolean forceCreation, OpenMMEnergy openMMEnergy) {
    // Check if the current context is consistent with the requested context.
    if (contextPointer != null && !forceCreation) {
      if (this.temperature == temperature && this.timeStep == timeStep
          && this.integratorString.equalsIgnoreCase(integratorString)) {
        // All requested features agree.
        return this;
      }
    }

    this.integratorString = integratorString;
    this.timeStep = timeStep;
    this.temperature = temperature;

    if (contextPointer != null) {
      logger.fine(" Free OpenMM Context.");
      OpenMM_Context_destroy(contextPointer);
      logger.fine(" Free OpenMM Context completed.");
      contextPointer = null;
    }

    logger.info("\n Creating OpenMM Context");

    OpenMMSystem openMMSystem = openMMEnergy.getSystem();
    PointerByReference integratorPointer = openMMIntegrator.createIntegrator(integratorString, this.timeStep, temperature, openMMSystem);

    // Set lambda to 1.0 when creating a context to avoid OpenMM compiling out any terms.
    double currentLambda = openMMEnergy.getLambda();

    if (openMMEnergy.getLambdaTerm()) {
      openMMEnergy.setLambda(1.0);
    }

    // Create a context.
    contextPointer = OpenMM_Context_create_2(openMMSystem.getSystem(), integratorPointer, platformPointer);

    // Revert to the current lambda value.
    if (openMMEnergy.getLambdaTerm()) {
      openMMEnergy.setLambda(currentLambda);
    }

    // Get initial positions and velocities for active atoms.
    int nVar = openMMEnergy.getNumberOfVariables();
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
    Crystal crystal = openMMEnergy.getCrystal();
    setPeriodicBoxVectors(crystal);

    // Load current atomic positions.
    setOpenMMPositions(x);

    // Load current velocities.
    setOpenMMVelocities(v);

    // Apply constraints starting from current atomic positions.
    OpenMM_Context_applyConstraints(contextPointer, constraintTolerance);

    // Application of constraints can change coordinates and velocities.
    // Retrieve them for consistency.
    OpenMMState openMMState = new OpenMMState(this, OpenMM_State_Positions | OpenMM_State_Velocities,
        atoms, openMMEnergy.getNumberOfVariables());
    openMMState.getPositions(x);
    openMMState.getVelocities(v);
    openMMState.free();

    return this;
  }

  /**
   * Reinitialize the context.
   *
   * <p>When a Context is created, it may cache information about the System being simulated and
   * the Force objects contained in it. This means that, if the System or Forces are then modified,
   * the Context might not see all changes. Call reinitialize() to force the Context to rebuild its
   * internal representation of the System and pick up any changes that have been made.
   *
   * <p>This is an expensive operation, so you should try to avoid calling it too frequently.
   */
  void reinitContext() {
    if (contextPointer != null) {
      int preserveState = 1;
      OpenMM_Context_reinitialize(contextPointer, preserveState);
    }
  }

  /**
   * Set the value of an adjustable parameter defined by a Force object in the System.
   *
   * @param name  the name of the parameter to set.
   * @param value the value of the parameter.
   */
  void setParameter(String name, double value) {
    if (contextPointer != null) {
      OpenMM_Context_setParameter(contextPointer, name, value);
    }
  }

  /**
   * Load an OpenMM Platform
   */
  private void loadPlatform(ForceFieldEnergy.Platform requestedPlatform) {

    OpenMMUtils.init();

    logger.log(Level.INFO, " Loaded from:\n {0}", OpenMMLibrary.JNA_NATIVE_LIB.toString());

    // Print out the OpenMM Version.
    Pointer version = OpenMM_Platform_getOpenMMVersion();
    logger.log(Level.INFO, " Version: {0}", version.getString(0));

    // Print out the OpenMM lib directory.
    logger.log(Level.FINE, " Lib Directory:       {0}", OpenMMUtils.getLibDirectory());
    // Load platforms and print out their names.
    PointerByReference libs = OpenMM_Platform_loadPluginsFromDirectory(
        OpenMMUtils.getLibDirectory());
    int numLibs = OpenMM_StringArray_getSize(libs);
    logger.log(Level.FINE, " Number of libraries: {0}", numLibs);
    for (int i = 0; i < numLibs; i++) {
      String libString = OpenMMEnergy.stringFromArray(libs, i);
      logger.log(Level.FINE, "  Library: {0}", libString);
    }
    OpenMM_StringArray_destroy(libs);

    // Print out the OpenMM plugin directory.
    logger.log(Level.INFO, "\n Plugin Directory:  {0}", OpenMMUtils.getPluginDirectory());
    // Load plugins and print out their names.
    PointerByReference plugins = OpenMM_Platform_loadPluginsFromDirectory(
        OpenMMUtils.getPluginDirectory());
    int numPlugins = OpenMM_StringArray_getSize(plugins);
    logger.log(Level.INFO, " Number of Plugins: {0}", numPlugins);
    boolean cuda = false;
    for (int i = 0; i < numPlugins; i++) {
      String pluginString = OpenMMEnergy.stringFromArray(plugins, i);
      logger.log(Level.INFO, "  Plugin: {0}", pluginString);
      if (pluginString != null) {
        pluginString = pluginString.toUpperCase();
        boolean amoebaCudaAvailable = pluginString.contains("AMOEBACUDA");
        if (amoebaCudaAvailable) {
          cuda = true;
        }
      }
    }
    OpenMM_StringArray_destroy(plugins);

    int numPlatforms = OpenMM_Platform_getNumPlatforms();
    logger.log(Level.INFO, " Number of Platforms: {0}", numPlatforms);

    if (requestedPlatform == ForceFieldEnergy.Platform.OMM_CUDA && !cuda) {
      logger.severe(" The OMM_CUDA platform was requested, but is not available.");
    }

    // Extra logging to print out plugins that failed to load.
    if (logger.isLoggable(Level.FINE)) {
      PointerByReference pluginFailers = OpenMM_Platform_getPluginLoadFailures();
      int numFailures = OpenMM_StringArray_getSize(pluginFailers);
      for (int i = 0; i < numFailures; i++) {
        String pluginString = OpenMMEnergy.stringFromArray(pluginFailers, i);
        logger.log(Level.FINE, " Plugin load failure: {0}", pluginString);
      }
      OpenMM_StringArray_destroy(pluginFailers);
    }

    String defaultPrecision = "mixed";
    MolecularAssembly molecularAssembly = openMMEnergy.getMolecularAssembly();
    String precision = molecularAssembly.getForceField().getString("PRECISION", defaultPrecision)
        .toLowerCase();
    precision = precision.replace("-precision", "");
    switch (precision) {
      case "double", "mixed", "single" -> logger.info(format(" Precision level: %s", precision));
      default -> {
        logger.info(format(" Could not interpret precision level %s, defaulting to %s", precision, defaultPrecision));
        precision = defaultPrecision;
      }
    }

    if (cuda && requestedPlatform != ForceFieldEnergy.Platform.OMM_REF) {
      int defaultDevice = OpenMMEnergy.getDefaultDevice(molecularAssembly.getProperties());
      platformPointer = OpenMM_Platform_getPlatformByName("CUDA");
      int deviceID = molecularAssembly.getForceField().getInteger("CUDA_DEVICE", defaultDevice);
      String deviceIDString = Integer.toString(deviceID);

      OpenMMLibrary.OpenMM_Platform_setPropertyDefaultValue(platformPointer, OpenMMEnergy.pointerForString("CudaDeviceIndex"),
          OpenMMEnergy.pointerForString(deviceIDString));
      OpenMMLibrary.OpenMM_Platform_setPropertyDefaultValue(platformPointer, OpenMMEnergy.pointerForString("Precision"),
          OpenMMEnergy.pointerForString(precision));
      logger.info(format(" Platform: AMOEBA CUDA (Device ID %d)", deviceID));
      try {
        Comm world = Comm.world();
        if (world != null) {
          logger.info(format(" Running on host %s, rank %d", world.host(), world.rank()));
        }
      } catch (IllegalStateException ise) {
        logger.fine(" Could not find the world communicator!");
      }
    } else {
      platformPointer = OpenMM_Platform_getPlatformByName("Reference");
      logger.info(" Platform: AMOEBA CPU Reference");
    }
  }
}
