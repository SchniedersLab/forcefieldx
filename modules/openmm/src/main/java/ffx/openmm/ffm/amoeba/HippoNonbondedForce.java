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
package ffx.openmm.ffm.amoeba;

import ffx.openmm.ffm.Context;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.OpenMMBooleans;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.Vec3Array;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.Objects;

/**
 * HIPPO combined nonbonded interactions: permanent electrostatics, induced polarization, charge
 * transfer, dispersion, and Pauli repulsion. Add one set of particle parameters for each System
 * particle; pair-specific scale exceptions can remove or reduce interactions. {@code NoCutoff}
 * (0) is the native default; {@code PME} (1) enables periodic electrostatics and dispersion PME.
 * Updating an existing Context copies particle and exception parameters only; other force settings
 * require Context reinitialization.
 */
public class HippoNonbondedForce extends Force {

  /**
   * Particle multipole, repulsion, dispersion, and charge-transfer parameters. Dipole and
   * quadrupole arrays are copied to and from native storage; their expected lengths are three and
   * nine, respectively.
   *
   * @param charge         particle monopole charge in elementary-charge units.
   * @param dipole         molecular-frame dipole components in elementary-charge nanometers, length 3.
   * @param quadrupole     molecular-frame quadrupole components in elementary-charge nanometers
   *                       squared, length 9 in the native ordering.
   * @param coreCharge     charge assigned to the atomic core, in elementary-charge units.
   * @param alpha          width parameter for the particle electron density, in inverse nanometers.
   * @param epsilon        parameter controlling charge-transfer magnitude, in OpenMM energy units
   *                       (kilojoules per mole).
   * @param damping        charge-transfer length scale, in nanometers.
   * @param c6             coefficient of the dispersion interaction, in the native energy-length^6 units.
   * @param pauliK         coefficient of Pauli repulsion.
   * @param pauliQ         charge used in the Pauli repulsion model, in elementary-charge units.
   * @param pauliAlpha     electron-density width used for Pauli repulsion, in inverse nanometers.
   * @param polarizability atomic polarizability, in nanometers cubed.
   * @param axisType       local multipole-frame enum: ZThenX=0, Bisector=1, ZBisect=2,
   *                       ThreeFold=3, ZOnly=4, or NoAxisType=5.
   * @param multipoleAtomZ particle index defining the local-frame Z direction, or -1 if unused.
   * @param multipoleAtomX particle index defining the local-frame X direction, or -1 if unused.
   * @param multipoleAtomY particle index defining the local frame's third reference, or -1 if unused.
   */
  public record ParticleParameters(
      double charge, double[] dipole, double[] quadrupole, double coreCharge, double alpha,
      double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
      double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX,
      int multipoleAtomY) {
  }

  /**
   * Pair-specific scale factors for the six HIPPO interaction classes. Setting every factor to
   * zero completely omits the pair interaction.
   *
   * @param particle1               first particle index.
   * @param particle2               second particle index.
   * @param multipoleMultipoleScale scale for charge-charge/multipole interactions.
   * @param dipoleMultipoleScale    scale for dipole-multipole interactions.
   * @param dipoleDipoleScale       scale for dipole-dipole interactions.
   * @param dispersionScale         scale for dispersion.
   * @param repulsionScale          scale for Pauli repulsion.
   * @param chargeTransferScale     scale for charge transfer.
   */
  public record ExceptionParameters(
      int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale,
      double dipoleDipoleScale, double dispersionScale, double repulsionScale,
      double chargeTransferScale) {
  }

  /**
   * Ewald or dispersion-PME grid parameters.
   *
   * @param alpha Ewald separation parameter in inverse nanometers; zero requests automatic choice
   *              based on the Ewald error tolerance.
   * @param nx    grid points along X.
   * @param ny    grid points along Y.
   * @param nz    grid points along Z.
   */
  public record PmeParameters(double alpha, int nx, int ny, int nz) {
  }

  /**
   * Create an empty HIPPO nonbonded force. Its native nonbonded method defaults to {@code NoCutoff}.
   *
   * @throws IllegalStateException if the native force cannot be created.
   */
  public HippoNonbondedForce() {
    super(create());
  }

  /**
   * Add a pair-specific exception. Scale factors independently scale multipole-multipole,
   * dipole-multipole, dipole-dipole, dispersion, repulsion, and charge-transfer contributions.
   *
   * @param particle1               first particle index.
   * @param particle2               second particle index.
   * @param multipoleMultipoleScale multipole-multipole interaction scale.
   * @param dipoleMultipoleScale    dipole-multipole interaction scale.
   * @param dipoleDipoleScale       dipole-dipole interaction scale.
   * @param dispersionScale         dispersion interaction scale.
   * @param repulsionScale          Pauli repulsion scale.
   * @param chargeTransferScale     charge-transfer interaction scale.
   * @param replace                 if true, replace an existing exception for this pair instead of adding a
   *                                duplicate; if false, add a new exception.
   * @return index of the added or replaced exception.
   */
  public int addException(
      int particle1, int particle2, double multipoleMultipoleScale,
      double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale,
      double repulsionScale, double chargeTransferScale, boolean replace) {
    return OpenMMNative.OpenMM_HippoNonbondedForce_addException(
        getPointer(), particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale,
        dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale,
        OpenMMBooleans.toNative(replace));
  }

  /**
   * Add one HIPPO particle, copying the supplied molecular multipole arrays into native storage.
   * Call once per System particle, in particle order.
   *
   * @param parameters particle parameters; dipole must have length 3 and quadrupole length 9.
   * @return index assigned to the particle.
   * @throws NullPointerException if {@code parameters} or either multipole array is null.
   */
  public int addParticle(ParticleParameters parameters) {
    Objects.requireNonNull(parameters, "Particle parameters cannot be null.");
    try (DoubleArray dipole = toNative(parameters.dipole());
         DoubleArray quadrupole = toNative(parameters.quadrupole())) {
      return OpenMMNative.OpenMM_HippoNonbondedForce_addParticle(
          getPointer(), parameters.charge(), dipole.getPointer(), quadrupole.getPointer(),
          parameters.coreCharge(), parameters.alpha(), parameters.epsilon(), parameters.damping(),
          parameters.c6(), parameters.pauliK(), parameters.pauliQ(), parameters.pauliAlpha(),
          parameters.polarizability(), parameters.axisType(), parameters.multipoleAtomZ(),
          parameters.multipoleAtomX(), parameters.multipoleAtomY());
    }
  }

  /** Release the native force owned by this façade. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_HippoNonbondedForce_destroy);
  }

  /**
   * Get the nonbonded cutoff distance.
   *
   * @return cutoff distance in nanometers; it has no effect with {@code NoCutoff}.
   */
  public double getCutoffDistance() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getCutoffDistance(getPointer());
  }

  /**
   * Get the configured dispersion-PME parameters.
   *
   * @return copied DPME alpha and grid dimensions; alpha zero means the values are chosen
   * automatically using the Ewald error tolerance.
   */
  public PmeParameters getDPMEParameters() {
    return readPmeParameters(null, true);
  }

  /**
   * Get dispersion-PME parameters actually used by a Context.
   *
   * @param context Context whose platform-selected parameters are queried.
   * @return copied DPME alpha and grid dimensions used by that Context.
   * @throws NullPointerException if {@code context} is null.
   */
  public PmeParameters getDPMEParametersInContext(Context context) {
    return readPmeParameters(Objects.requireNonNull(context, "Context cannot be null."), true);
  }

  /**
   * Get the target Ewald relative error tolerance.
   *
   * @return dimensionless relative error tolerance.
   */
  public double getEwaldErrorTolerance() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(getPointer());
  }

  /**
   * Get copied parameters for a pair-specific exception.
   *
   * @param index exception index.
   * @return particle indices and six interaction scale factors.
   */
  public ExceptionParameters getExceptionParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment particle1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment particle2 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment[] scales = allocateDoubles(arena, 6);
      OpenMMNative.OpenMM_HippoNonbondedForce_getExceptionParameters(
          getPointer(), index, particle1, particle2, scales[0], scales[1], scales[2], scales[3],
          scales[4], scales[5]);
      return new ExceptionParameters(particle1.get(ValueLayout.JAVA_INT, 0),
          particle2.get(ValueLayout.JAVA_INT, 0), get(scales[0]), get(scales[1]), get(scales[2]),
          get(scales[3]), get(scales[4]), get(scales[5]));
    }
  }

  /**
   * Return the force-owned extrapolation coefficients as a borrowed native double-array handle.
   * Do not destroy the handle; it becomes invalid when this force is destroyed.
   *
   * @return borrowed native array handle containing the induced-dipole extrapolation coefficients.
   */
  public MemorySegment getExtrapolationCoefficients() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(getPointer());
  }

  /**
   * Retrieve induced dipoles from a Context into the supplied caller-owned vector array.
   *
   * @param context Context at which induced dipoles are evaluated.
   * @param dipoles destination native vector array, sized for the force's particles; each vector
   *     is an induced dipole in elementary-charge nanometers.
   * @throws NullPointerException if either argument is null.
   */
  public void getInducedDipoles(Context context, Vec3Array dipoles) {
    OpenMMNative.OpenMM_HippoNonbondedForce_getInducedDipoles(
        getPointer(), context.getPointer(), dipoles.getPointer());
  }

  /**
   * Retrieve permanent dipoles transformed to the laboratory frame into the supplied vector array.
   *
   * @param context Context providing particle positions and local-frame transforms.
   * @param dipoles destination native vector array, sized for the force's particles; each vector
   *     is a permanent dipole in elementary-charge nanometers.
   * @throws NullPointerException if either argument is null.
   */
  public void getLabFramePermanentDipoles(Context context, Vec3Array dipoles) {
    OpenMMNative.OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(
        getPointer(), context.getPointer(), dipoles.getPointer());
  }

  /**
   * Get the native long-range nonbonded method.
   *
   * @return {@code NoCutoff} (0) for direct nonperiodic interactions or {@code PME} (1) for
   * periodic Ewald summation.
   */
  public int getNonbondedMethod() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getNonbondedMethod(getPointer());
  }

  /**
   * Get the number of pair-specific exceptions.
   *
   * @return exception count.
   */
  public int getNumExceptions() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getNumExceptions(getPointer());
  }

  /**
   * Get the number of configured particle parameter sets.
   *
   * @return particle count.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getNumParticles(getPointer());
  }

  /**
   * Get the configured PME parameters for electrostatics.
   *
   * @return copied Ewald alpha and grid dimensions; alpha zero requests automatic values based on
   * the Ewald error tolerance.
   */
  public PmeParameters getPMEParameters() {
    return readPmeParameters(null, false);
  }

  /**
   * Get electrostatic PME parameters actually used by a Context. A platform may choose grid
   * dimensions different from the configured or automatically calculated dimensions.
   *
   * @param context Context whose platform-selected parameters are queried.
   * @return copied alpha and grid dimensions used by that Context.
   * @throws NullPointerException if {@code context} is null.
   */
  public PmeParameters getPMEParametersInContext(Context context) {
    return readPmeParameters(Objects.requireNonNull(context, "Context cannot be null."), false);
  }

  /**
   * Get copied parameters for one particle. Returned dipole and quadrupole arrays are independent
   * Java copies, not views of native storage.
   *
   * @param index particle index.
   * @return copied HIPPO particle parameters.
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena arena = Arena.ofConfined();
         DoubleArray dipole = new DoubleArray(3);
         DoubleArray quadrupole = new DoubleArray(9)) {
      MemorySegment charge = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment coreCharge = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment alpha = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment damping = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment c6 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment pauliK = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment pauliQ = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment pauliAlpha = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment polarizability = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment axisType = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment atomZ = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment atomX = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment atomY = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_HippoNonbondedForce_getParticleParameters(
          getPointer(), index, charge, dipole.getPointer(), quadrupole.getPointer(), coreCharge,
          alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType,
          atomZ, atomX, atomY);
      return new ParticleParameters(get(charge), copy(dipole), copy(quadrupole), get(coreCharge),
          get(alpha), get(epsilon), get(damping), get(c6), get(pauliK), get(pauliQ),
          get(pauliAlpha), get(polarizability), getInt(axisType), getInt(atomZ), getInt(atomX),
          getInt(atomY));
    }
  }

  /**
   * Get the distance where the switching function starts reducing repulsion and charge-transfer
   * interactions.
   *
   * @return switching distance in nanometers; it must be less than the cutoff distance.
   */
  public double getSwitchingDistance() {
    return OpenMMNative.OpenMM_HippoNonbondedForce_getSwitchingDistance(getPointer());
  }

  /**
   * Set the nonbonded cutoff distance.
   *
   * @param distance cutoff distance in nanometers; it has no effect with {@code NoCutoff}.
   */
  public void setCutoffDistance(double distance) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setCutoffDistance(getPointer(), distance);
  }

  /**
   * Configure dispersion-PME summation.
   *
   * @param alpha Ewald separation parameter in inverse nanometers; zero selects automatic values
   *              from the Ewald error tolerance.
   * @param nx    dispersion grid points along X, or zero for automatic selection.
   * @param ny    dispersion grid points along Y, or zero for automatic selection.
   * @param nz    dispersion grid points along Z, or zero for automatic selection.
   */
  public void setDPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setDPMEParameters(
        getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Set the target Ewald relative error tolerance.
   *
   * @param tolerance dimensionless relative error tolerance.
   */
  public void setEwaldErrorTolerance(double tolerance) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(getPointer(), tolerance);
  }

  /**
   * Replace an existing pair exception.
   *
   * @param index      exception index.
   * @param parameters replacement particle indices and six interaction scale factors.
   * @throws NullPointerException if {@code parameters} is null.
   */
  public void setExceptionParameters(int index, ExceptionParameters parameters) {
    Objects.requireNonNull(parameters, "Exception parameters cannot be null.");
    OpenMMNative.OpenMM_HippoNonbondedForce_setExceptionParameters(
        getPointer(), index, parameters.particle1(), parameters.particle2(),
        parameters.multipoleMultipoleScale(), parameters.dipoleMultipoleScale(),
        parameters.dipoleDipoleScale(), parameters.dispersionScale(), parameters.repulsionScale(),
        parameters.chargeTransferScale());
  }

  /**
   * Set coefficients for the induced-dipole extrapolation terms mu0, mu1, and so on. The number of
   * coefficients determines the number of extrapolation iterations. Native code copies the values.
   *
   * @param coefficients native double-array handle containing the coefficients.
   * @throws NullPointerException if {@code coefficients} is null.
   */
  public void setExtrapolationCoefficients(DoubleArray coefficients) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(
        getPointer(), coefficients.getPointer());
  }

  /**
   * Set extrapolation coefficients from a native double-array handle. The native force copies the
   * coefficients; the supplied handle remains caller-owned.
   *
   * @param coefficients valid {@code OpenMM_DoubleArray} handle.
   */
  public void setExtrapolationCoefficients(MemorySegment coefficients) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(
        getPointer(), coefficients);
  }

  /**
   * Set the long-range nonbonded method.
   *
   * @param method native enum value: NoCutoff=0 or PME=1.
   */
  public void setNonbondedMethod(int method) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setNonbondedMethod(getPointer(), method);
  }

  /**
   * Configure electrostatic PME summation.
   *
   * @param alpha Ewald separation parameter in inverse nanometers; zero selects automatic values
   *              from the Ewald error tolerance.
   * @param nx    grid points along X, or zero for automatic selection.
   * @param ny    grid points along Y, or zero for automatic selection.
   * @param nz    grid points along Z, or zero for automatic selection.
   */
  public void setPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setPMEParameters(
        getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Replace all parameters for one particle. Dipole and quadrupole arrays are copied.
   *
   * @param index      particle index.
   * @param parameters replacement particle parameters.
   * @throws NullPointerException if {@code parameters} or either multipole array is null.
   */
  public void setParticleParameters(int index, ParticleParameters parameters) {
    Objects.requireNonNull(parameters, "Particle parameters cannot be null.");
    try (DoubleArray dipole = toNative(parameters.dipole());
         DoubleArray quadrupole = toNative(parameters.quadrupole())) {
      OpenMMNative.OpenMM_HippoNonbondedForce_setParticleParameters(
          getPointer(), index, parameters.charge(), dipole.getPointer(), quadrupole.getPointer(),
          parameters.coreCharge(), parameters.alpha(), parameters.epsilon(), parameters.damping(),
          parameters.c6(), parameters.pauliK(), parameters.pauliQ(), parameters.pauliAlpha(),
          parameters.polarizability(), parameters.axisType(), parameters.multipoleAtomZ(),
          parameters.multipoleAtomX(), parameters.multipoleAtomY());
    }
  }

  /**
   * Set the switching distance for repulsion and charge-transfer interactions.
   *
   * @param distance switching distance in nanometers, strictly below the cutoff distance.
   */
  public void setSwitchingDistance(double distance) {
    OpenMMNative.OpenMM_HippoNonbondedForce_setSwitchingDistance(getPointer(), distance);
  }

  /**
   * Copy changed particle parameters and pair exceptions into an existing Context. Changes to
   * nonbonded method, cutoffs, PME settings, error tolerance, switching distance, and extrapolation
   * coefficients require Context reinitialization.
   *
   * @param context Context to update.
   * @throws NullPointerException if {@code context} is null.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_HippoNonbondedForce_updateParametersInContext(
        getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Report whether this force uses periodic boundary conditions.
   *
   * @return true when the configured method is PME.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment[] allocateDoubles(Arena arena, int count) {
    MemorySegment[] result = new MemorySegment[count];
    for (int i = 0; i < count; i++) {
      result[i] = arena.allocate(ValueLayout.JAVA_DOUBLE);
    }
    return result;
  }

  private static double[] copy(DoubleArray values) {
    double[] result = new double[values.getSize()];
    for (int i = 0; i < result.length; i++) {
      result[i] = values.get(i);
    }
    return result;
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_HippoNonbondedForce_create.makeInvoker().apply();
  }

  private static double get(MemorySegment value) {
    return value.get(ValueLayout.JAVA_DOUBLE, 0);
  }

  private static int getInt(MemorySegment value) {
    return value.get(ValueLayout.JAVA_INT, 0);
  }

  private PmeParameters readPmeParameters(Context context, boolean dpme) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment alpha = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment nx = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment ny = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment nz = arena.allocate(ValueLayout.JAVA_INT);
      if (context == null) {
        if (dpme) {
          OpenMMNative.OpenMM_HippoNonbondedForce_getDPMEParameters(
              getPointer(), alpha, nx, ny, nz);
        } else {
          OpenMMNative.OpenMM_HippoNonbondedForce_getPMEParameters(
              getPointer(), alpha, nx, ny, nz);
        }
      } else if (dpme) {
        OpenMMNative.OpenMM_HippoNonbondedForce_getDPMEParametersInContext(
            getPointer(), context.getPointer(), alpha, nx, ny, nz);
      } else {
        OpenMMNative.OpenMM_HippoNonbondedForce_getPMEParametersInContext(
            getPointer(), context.getPointer(), alpha, nx, ny, nz);
      }
      return new PmeParameters(get(alpha), getInt(nx), getInt(ny), getInt(nz));
    }
  }

  private static DoubleArray toNative(double[] values) {
    Objects.requireNonNull(values, "Values cannot be null.");
    DoubleArray result = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) {
      result.set(i, values[i]);
    }
    return result;
  }
}
