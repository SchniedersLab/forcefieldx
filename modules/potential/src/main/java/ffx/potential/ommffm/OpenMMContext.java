// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.ommffm;

import ffx.crystal.Crystal;
import ffx.openmm.ffm.Context;
import ffx.openmm.ffm.Integrator;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.Platform;
import ffx.openmm.ffm.Vec3;
import ffx.openmm.ffm.VerletIntegrator;
import ffx.openmm.ffm.bindings.OpenMMNative;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static ffx.openmm.ffm.OpenMMUnits.ANGSTROMS_PER_NM;
import static ffx.openmm.ffm.OpenMMUnits.KCAL_PER_KJ;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static ffx.potential.Platform.OMM;
import static ffx.potential.Platform.OMM_CUDA;
import static ffx.potential.Platform.OMM_OPENCL;
import static java.lang.String.format;
import static java.lang.foreign.MemorySegment.NULL;

/**
 * Executes an OpenMM simulation for an {@link OpenMMSystem} using the FFM backend.
 */
public class OpenMMContext extends Context {

  private static final Logger logger = Logger.getLogger(OpenMMContext.class.getName());

  private static final double DEFAULT_CONSTRAINT_TOLERANCE = 1.0e-5;

  private final OpenMMSystem openMMSystem;
  private String integratorName = "VERLET";
  private double timeStep = 0.001;
  private double temperature = 298.15;
  private final boolean enforcePBC;
  private final Atom[] atoms;

  public OpenMMContext(Platform platform, OpenMMSystem openMMSystem, Atom[] atoms) {
    super(openMMSystem, createIntegrator("VERLET", 0.001, 298.15, openMMSystem), platform);
    this.openMMSystem = openMMSystem;
    this.atoms = atoms;

    ForceField forceField = openMMSystem.getForceField();
    boolean aperiodic = openMMSystem.getCrystal().aperiodic();
    this.enforcePBC = forceField.getBoolean("ENFORCE_PBC", !aperiodic);
  }

  public void update(String integratorName, double timeStep, double temperature, boolean forceCreation) {
    if (hasContextPointer() && !forceCreation) {
      if (this.temperature == temperature && this.timeStep == timeStep
          && this.integratorName.equalsIgnoreCase(integratorName)) {
        return;
      }
    }

    this.integratorName = integratorName;
    this.timeStep = timeStep;
    this.temperature = temperature;

    logger.info("\n Updating OpenMM Context");

    Integrator newIntegrator = createIntegrator(integratorName, timeStep, temperature, openMMSystem);
    Platform newPlatform = new Platform(getPlatform().getName());
    updateContext(openMMSystem, newIntegrator, newPlatform);

    int nVar = atoms.length * 3;
    double[] x = new double[nVar];
    double[] v = new double[nVar];
    double[] vel3 = new double[3];
    int index = 0;
    for (Atom a : atoms) {
      a.getVelocity(vel3);
      x[index] = a.getX();
      v[index++] = vel3[0];
      x[index] = a.getY();
      v[index++] = vel3[1];
      x[index] = a.getZ();
      v[index++] = vel3[2];
    }

    Crystal crystal = openMMSystem.getCrystal();
    setPeriodicBoxVectors(crystal);
    setPositions(x);
    setVelocities(v);
    applyConstraints(DEFAULT_CONSTRAINT_TOLERANCE);

    getPositions(x);
    getVelocities(v);
    index = 0;
    for (Atom a : atoms) {
      a.moveTo(x[index], x[index + 1], x[index + 2]);
      vel3[0] = v[index++];
      vel3[1] = v[index++];
      vel3[2] = v[index++];
      a.setVelocity(vel3);
    }
  }

  public void integrate(int numSteps) {
    Integrator integrator = getIntegrator();
    integrator.step(numSteps);
  }

  public void optimize(double eps, int maxIterations) {
    OpenMMNative.OpenMM_LocalEnergyMinimizer_minimize(getPointer(),
        eps / (NM_PER_ANGSTROM * KCAL_PER_KJ), maxIterations, NULL);
  }

  public void update() {
    update(integratorName, timeStep, temperature, false);
  }

  public void free() {
    if (hasContextPointer()) {
      logger.fine(" Free OpenMM context.");
      destroy();
      logger.fine(" Free OpenMM context completed.");
    }
  }

  public OpenMMState getOpenMMState(int mask) {
    return getOpenMMState(mask, enforcePBC);
  }

  public OpenMMState getOpenMMState(int mask, boolean enforcePBC) {
    return new OpenMMState(OpenMMNative.OpenMM_Context_getState(getPointer(), mask, enforcePBC ? 1 : 0));
  }

  public void setPositions(double[] x) {
    int n = x.length;
    double[] xnm = new double[n];
    for (int i = 0; i < n; i++) {
      xnm[i] = x[i] * NM_PER_ANGSTROM;
    }
    super.setPositions(xnm);
  }

  public void setVelocities(double[] v) {
    int n = v.length;
    double[] vnm = new double[n];
    for (int i = 0; i < n; i++) {
      vnm[i] = v[i] * NM_PER_ANGSTROM;
    }
    super.setVelocities(vnm);
  }

  public void getPositions(double[] x) {
    try (OpenMMState state = getOpenMMState(OpenMMNative.OpenMM_State_Positions())) {
      state.getCoordinates(x);
    }
  }

  public void getActivePositions(double[] x) {
    try (OpenMMState state = getOpenMMState(OpenMMNative.OpenMM_State_Positions())) {
      state.getActiveCoordinates(x, atoms);
    }
  }

  public void getVelocities(double[] v) {
    try (OpenMMState state = getOpenMMState(OpenMMNative.OpenMM_State_Velocities())) {
      state.getVelocities(v);
    }
  }

  public void getActiveVelocities(double[] v) {
    try (OpenMMState state = getOpenMMState(OpenMMNative.OpenMM_State_Velocities())) {
      state.getActiveVelocities(v, atoms);
    }
  }

  public void setStep(long step) {
    super.setStepCount(step);
  }

  public long getStep() {
    return super.getStepCount();
  }

  public void setPeriodicBoxVectors(Crystal crystal) {
    if (!crystal.aperiodic()) {
      double[][] Ai = crystal.Ai;
      Vec3 a = new Vec3(Ai[0][0] * NM_PER_ANGSTROM, Ai[0][1] * NM_PER_ANGSTROM, Ai[0][2] * NM_PER_ANGSTROM);
      Vec3 b = new Vec3(Ai[1][0] * NM_PER_ANGSTROM, Ai[1][1] * NM_PER_ANGSTROM, Ai[1][2] * NM_PER_ANGSTROM);
      Vec3 c = new Vec3(Ai[2][0] * NM_PER_ANGSTROM, Ai[2][1] * NM_PER_ANGSTROM, Ai[2][2] * NM_PER_ANGSTROM);
      setPeriodicBoxVectors(a, b, c);
    }
  }

  public void getPeriodicBoxVectors(double[][] box) {
    try (OpenMMState state = getOpenMMState(OpenMMNative.OpenMM_State_Positions())) {
      state.getPeriodicBoxVectors(box);
    }
  }

  public static Platform loadPlatform(ffx.potential.Platform requestedPlatform, ForceField forceField) {
    OpenMMRuntime.initialize();

    logger.log(Level.INFO, " OpenMM Runtime Initialized");
    logger.log(Level.INFO, " Version: {0}", Platform.getOpenMMVersion());

    int numPlatforms = Platform.getNumPlatforms();
    logger.log(Level.INFO, " Number of Platforms: {0}", numPlatforms);
    boolean cuda = false;
    boolean opencl = false;
    for (int i = 0; i < numPlatforms; i++) {
      Platform p = Platform.getPlatform(i);
      String name = p.getName().toUpperCase();
      logger.log(Level.INFO, "  Platform: {0}", name);
      if (name.contains("CUDA")) {
        cuda = true;
      }
      if (name.contains("OPENCL")) {
        opencl = true;
      }
    }

    if (requestedPlatform == OMM_CUDA && !cuda) {
      logger.severe(" The OMM_CUDA platform was requested, but is not available.");
    }

    if (requestedPlatform == OMM_OPENCL && !opencl) {
      logger.severe(" The OMM_OPENCL platform was requested, but is not available.");
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
      int defaultDevice = getDefaultDevice(forceField.getProperties());
      openMMPlatform = new Platform("CUDA");
      int deviceID = forceField.getInteger("CUDA_DEVICE", defaultDevice);
      deviceID = forceField.getInteger("DeviceIndex", deviceID);
      openMMPlatform.setPropertyDefaultValue("DeviceIndex", Integer.toString(deviceID));
      openMMPlatform.setPropertyDefaultValue("Precision", precision);
      logger.info(format(" Platform: %s (Device Index %d)", openMMPlatform.getName(), deviceID));
    } else if (opencl && (requestedPlatform == OMM_OPENCL || requestedPlatform == OMM)) {
      int defaultDevice = getDefaultDevice(forceField.getProperties());
      openMMPlatform = new Platform("OpenCL");
      int deviceID = forceField.getInteger("DeviceIndex", defaultDevice);
      int openCLPlatformIndex = forceField.getInteger("OpenCLPlatformIndex", 0);
      openMMPlatform.setPropertyDefaultValue("DeviceIndex", Integer.toString(deviceID));
      openMMPlatform.setPropertyDefaultValue("OpenCLPlatformIndex", Integer.toString(openCLPlatformIndex));
      openMMPlatform.setPropertyDefaultValue("Precision", precision);
      logger.info(format(" Platform: %s (Platform Index %d, Device Index %d)",
          openMMPlatform.getName(), openCLPlatformIndex, deviceID));
    } else {
      openMMPlatform = new Platform("Reference");
      logger.info(format(" Platform: %s", openMMPlatform.getName()));
    }

    return openMMPlatform;
  }

  public static int getDefaultDevice(CompositeConfiguration props) {
    String availDeviceProp = props.getString("availableDevices", props.getString("CUDA_DEVICES"));
    if (availDeviceProp == null) {
      int nDevs = props.getInt("numCudaDevices", 1);
      availDeviceProp = IntStream.range(0, nDevs).mapToObj(Integer::toString)
          .collect(Collectors.joining(" "));
    }
    availDeviceProp = availDeviceProp.trim();
    String[] availDevices = availDeviceProp.split("\\s+");
    return Integer.parseInt(availDevices[0]);
  }

  private static Integrator createIntegrator(String integratorName, double timeStep,
                                             double temperature, OpenMMSystem openMMSystem) {
    return OpenMMIntegrator.createIntegrator(integratorName, timeStep, temperature, openMMSystem);
  }
}
