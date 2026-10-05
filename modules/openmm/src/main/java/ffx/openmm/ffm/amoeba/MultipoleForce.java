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
import ffx.openmm.ffm.IntArray;
import ffx.openmm.ffm.OpenMMBooleans;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.Vec3Array;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.Objects;

/**
 * AMOEBA permanent multipole electrostatics and polarization force. Add one multipole for each
 * System particle. Native nonbonded methods are {@code NoCutoff} (0, the default) and {@code PME}
 * (1); polarization modes are {@code Mutual} (0), {@code Direct} (1), and {@code Extrapolated}
 * (2). Particle changes can be copied to an existing Context, but other force settings require
 * Context reinitialization.
 */
public class MultipoleForce extends Force {

  /**
   * Multipole parameters for one particle. Dipole and quadrupole arrays are copied Java arrays,
   * respectively containing three and nine molecular-frame components.
   *
   * @param charge              particle charge in elementary-charge units.
   * @param molecularDipole     molecular-frame dipole components in elementary-charge nanometers,
   *                            length 3.
   * @param molecularQuadrupole molecular-frame quadrupole components in elementary-charge
   *                            nanometers squared, length 9 in the native ordering.
   * @param axisType            local-frame enum: ZThenX=0, Bisector=1, ZBisect=2, ThreeFold=3, ZOnly=4,
   *                            or NoAxisType=5.
   * @param multipoleAtomZ      index of the first reference atom defining the local frame, or -1 if
   *                            unused.
   * @param multipoleAtomX      index of the second reference atom defining the local frame, or -1 if
   *                            unused.
   * @param multipoleAtomY      index of the third reference atom defining the local frame, or -1 if
   *                            unused.
   * @param thole               Thole damping parameter.
   * @param dampingFactor       polarization damping-factor parameter.
   * @param polarity            atomic polarizability parameter.
   */
  public record MultipoleParameters(
      double charge, double[] molecularDipole, double[] molecularQuadrupole, int axisType,
      int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY,
      double thole, double dampingFactor, double polarity) {
  }

  /**
   * PME reciprocal-space parameters.
   *
   * @param alpha Ewald separation parameter in inverse nanometers; zero selects automatic values
   *              based on the Ewald error tolerance.
   * @param nx    PME grid points along X.
   * @param ny    PME grid points along Y.
   * @param nz    PME grid points along Z.
   */
  public record PmeParameters(double alpha, int nx, int ny, int nz) {
  }

  /**
   * Create an empty AMOEBA multipole force. Native defaults are NoCutoff nonbonded treatment and
   * Mutual polarization.
   *
   * @throws IllegalStateException if the native force cannot be created.
   */
  public MultipoleForce() {
    super(create());
  }

  /**
   * Add multipole information for one particle. Call once per System particle, in particle order.
   * Native code copies the molecular dipole and quadrupole values.
   *
   * @param charge              particle charge in elementary-charge units.
   * @param molecularDipole     caller-owned native double array with three molecular-frame dipole
   *                            components, in elementary-charge nanometers.
   * @param molecularQuadrupole caller-owned native double array with nine molecular-frame
   *                            quadrupole components, in elementary-charge nanometers squared.
   * @param axisType            local-frame enum value: ZThenX=0, Bisector=1, ZBisect=2, ThreeFold=3,
   *                            ZOnly=4, or NoAxisType=5.
   * @param multipoleAtomZ      index of the first local-frame reference atom, or -1 if unused.
   * @param multipoleAtomX      index of the second local-frame reference atom, or -1 if unused.
   * @param multipoleAtomY      index of the third local-frame reference atom, or -1 if unused.
   * @param thole               Thole damping parameter.
   * @param dampingFactor       polarization damping-factor parameter.
   * @param polarity            atomic polarizability parameter.
   * @return index assigned to the added multipole.
   */
  public int addMultipole(
      double charge, DoubleArray molecularDipole, DoubleArray molecularQuadrupole,
      int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY,
      double thole, double dampingFactor, double polarity) {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_addMultipole(
        getPointer(), charge, molecularDipole.getPointer(), molecularQuadrupole.getPointer(),
        axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
  }

  /** Release the native force owned by this façade. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_AmoebaMultipoleForce_destroy);
  }

  /**
   * Get the legacy Ewald alpha parameter.
   *
   * @return Ewald separation parameter in inverse nanometers.
   * @deprecated Use {@link #getPMEParameters()}.
   */
  @Deprecated
  public double getAEwald() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getAEwald(getPointer());
  }

  /**
   * Copy one covalent map into a new caller-owned integer array. Covalent type values are
   * Covalent12=0, Covalent13=1, Covalent14=2, Covalent15=3, PolarizationCovalent11=4,
   * PolarizationCovalent12=5, PolarizationCovalent13=6, and PolarizationCovalent14=7.
   *
   * @param particle     particle index.
   * @param covalentType native covalent-map category value from 0 through 7.
   * @return caller-owned {@link IntArray}; destroy it when no longer needed.
   */
  public IntArray getCovalentMap(int particle, int covalentType) {
    IntArray values = new IntArray(0);
    try {
      OpenMMNative.OpenMM_AmoebaMultipoleForce_getCovalentMap(
          getPointer(), particle, covalentType, values.getPointer());
      return values;
    } catch (RuntimeException | Error exception) {
      values.destroy();
      throw exception;
    }
  }

  /**
   * Fill a caller-supplied opaque {@code OpenMM_2D_IntArray} handle with all covalent maps for a
   * particle. The installed OpenMM C wrapper exports no constructor or accessor for this type, so
   * ordinary Java callers cannot create or inspect a suitable value; this overload is usable only
   * when a valid ABI-compatible handle is obtained externally. The native routine writes into the
   * supplied handle; this façade neither owns nor destroys it.
   *
   * @param particle     particle index.
   * @param covalentMaps valid opaque {@code OpenMM_2D_IntArray} handle.
   * @throws NullPointerException if {@code covalentMaps} is null.
   */
  public void getCovalentMaps(int particle, MemorySegment covalentMaps) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_getCovalentMaps(
        getPointer(), particle, covalentMaps);
  }

  /**
   * Get all covalent maps for a particle.
   *
   * <p>The JNA counterpart passes a one-dimensional array where the native API requires an opaque
   * {@code OpenMM_2D_IntArray}. The installed C wrapper exports no constructor or accessor for
   * that type, so this overload cannot return a meaningful array and always fails. Use
   * {@link #getCovalentMaps(int, MemorySegment)} only when an ABI-compatible handle is available.</p>
   *
   * @param particle particle index.
   * @return not available through the installed C wrapper.
   * @throws UnsupportedOperationException because the wrapper lacks 2D-array accessors.
   */
  public IntArray getCovalentMaps(int particle) {
    throw new UnsupportedOperationException(
        "OpenMM_2D_IntArray has no exported constructor or accessor in the installed C wrapper.");
  }

  /**
   * Get the nonbonded cutoff distance.
   *
   * @return cutoff distance in nanometers; it has no effect with NoCutoff.
   */
  public double getCutoffDistance() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getCutoffDistance(getPointer());
  }

  /**
   * Calculate electrostatic potential at the specified points for the supplied Context. The
   * returned native double array is owned by the caller and must be destroyed.
   *
   * @param context Context at which the potential is evaluated.
   * @param points  caller-owned native array of evaluation positions in nanometers.
   * @return caller-owned native array of electrostatic potentials in kilojoules per mole per
   *     elementary charge; destroy it when no longer needed.
   * @throws NullPointerException if either argument is null.
   */
  public DoubleArray getElectrostaticPotential(Context context, Vec3Array points) {
    DoubleArray potential = new DoubleArray(0);
    try {
      OpenMMNative.OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(
          getPointer(), points.getPointer(), context.getPointer(), potential.getPointer());
      return potential;
    } catch (RuntimeException | Error exception) {
      potential.destroy();
      throw exception;
    }
  }

  /**
   * Get the requested fractional force error tolerance for Ewald summation.
   *
   * @return dimensionless Ewald error tolerance.
   */
  public double getEwaldErrorTolerance() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(getPointer());
  }

  /**
   * Return the force-owned induced-dipole extrapolation coefficients as a borrowed native array
   * handle. The current header documents default coefficients {@code [-0.154, 0.017, 0.658,
   * 0.474]}. Do not destroy the handle; it becomes invalid when this force is destroyed.
   *
   * @return borrowed native double-array handle.
   */
  public MemorySegment getExtrapolationCoefficients() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(getPointer());
  }

  /**
   * Retrieve induced dipoles from a Context into the supplied caller-owned vector array.
   *
   * @param context        Context from which induced dipoles are obtained.
   * @param inducedDipoles destination native vector array sized for the force particles; each
   *     vector contains an induced dipole in elementary-charge nanometers.
   * @throws NullPointerException if either argument is null.
   */
  public void getInducedDipoles(Context context, Vec3Array inducedDipoles) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_getInducedDipoles(
        getPointer(), context.getPointer(), inducedDipoles.getPointer());
  }

  /**
   * Retrieve permanent dipoles transformed into the laboratory frame.
   *
   * @param context Context providing positions and local-frame orientations.
   * @param dipoles destination native vector array sized for the force particles; each vector
   *     contains a permanent dipole in elementary-charge nanometers.
   * @throws NullPointerException if either argument is null.
   */
  public void getLabFramePermanentDipoles(Context context, Vec3Array dipoles) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(
        getPointer(), context.getPointer(), dipoles.getPointer());
  }

  /**
   * Get copied multipole parameters. The returned record contains Java copies of the native
   * molecular dipole and quadrupole values.
   *
   * @param index particle/multipole index.
   * @return copied particle multipole parameters.
   */
  public MultipoleParameters getMultipoleParameters(int index) {
    try (Arena arena = Arena.ofConfined();
         DoubleArray dipole = new DoubleArray(3);
         DoubleArray quadrupole = new DoubleArray(9)) {
      MemorySegment charge = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment axis = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment atomZ = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment atomX = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment atomY = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment thole = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment damping = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment polarity = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_AmoebaMultipoleForce_getMultipoleParameters(
          getPointer(), index, charge, dipole.getPointer(), quadrupole.getPointer(), axis,
          atomZ, atomX, atomY, thole, damping, polarity);
      return new MultipoleParameters(charge.get(ValueLayout.JAVA_DOUBLE, 0), copy(dipole),
          copy(quadrupole), axis.get(ValueLayout.JAVA_INT, 0), atomZ.get(ValueLayout.JAVA_INT, 0),
          atomX.get(ValueLayout.JAVA_INT, 0), atomY.get(ValueLayout.JAVA_INT, 0),
          thole.get(ValueLayout.JAVA_DOUBLE, 0), damping.get(ValueLayout.JAVA_DOUBLE, 0),
          polarity.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the maximum number of iterations for mutually induced dipoles.
   *
   * @return iteration limit.
   */
  public int getMutualInducedMaxIterations() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(getPointer());
  }

  /**
   * Get the convergence target for mutually induced dipoles.
   *
   * @return target epsilon used by the iterative solver.
   */
  public double getMutualInducedTargetEpsilon() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(getPointer());
  }

  /**
   * Get the native long-range nonbonded method.
   *
   * @return NoCutoff (0) for direct nonperiodic interactions or PME (1) for periodic Ewald
   * summation.
   */
  public int getNonbondedMethod() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getNonbondedMethod(getPointer());
  }

  /**
   * Get the number of configured multipoles.
   *
   * @return multipole count.
   */
  public int getNumMultipoles() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getNumMultipoles(getPointer());
  }

  /**
   * Get configured PME reciprocal-space parameters. When alpha is zero (the default), the
   * configured grid values are ignored and selected automatically from the Ewald error tolerance.
   *
   * @return copied Ewald alpha and grid dimensions; alpha zero requests automatic parameter
   *     selection.
   */
  public PmeParameters getPMEParameters() {
    return readPmeParameters(null);
  }

  /**
   * Get the PME parameters actually selected for a Context. Platform restrictions may cause the
   * grid dimensions to differ from the configured or automatically selected values.
   *
   * @param context Context whose PME parameters are queried.
   * @return copied Ewald alpha and grid dimensions used by the Context.
   * @throws NullPointerException if {@code context} is null.
   */
  public PmeParameters getPMEParametersInContext(Context context) {
    return readPmeParameters(Objects.requireNonNull(context, "Context cannot be null."));
  }

  /**
   * Get the configured PME grid dimensions as a caller-owned native integer array.
   *
   * @return owned {@link IntArray} containing the X, Y, and Z dimensions; destroy it when no
   * longer needed.
   */
  public IntArray getPmeGridDimensions() {
    IntArray dimensions = new IntArray(0);
    try {
      OpenMMNative.OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(
          getPointer(), dimensions.getPointer());
      return dimensions;
    } catch (RuntimeException | Error exception) {
      dimensions.destroy();
      throw exception;
    }
  }

  /**
   * Get the PME B-spline interpolation order.
   *
   * @return B-spline order.
   */
  public int getPmeBSplineOrder() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(getPointer());
  }

  /**
   * Get the native polarization approximation. Mutual iterates induced dipoles to the configured
   * convergence target; Direct uses only fixed multipoles; Extrapolated performs a limited
   * perturbation series and extrapolates using the configured coefficients.
   *
   * @return Mutual (0), Direct (1), or Extrapolated (2).
   */
  public int getPolarizationType() {
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_getPolarizationType(getPointer());
  }

  /**
   * Retrieve the system's multipole moments from a Context into a caller-owned native array.
   *
   * @param context Context at which system moments are computed.
   * @param moments destination native double array for the moment values: charge in elementary
   *     charges, dipole in elementary-charge nanometers, and quadrupole in elementary-charge
   *     nanometers squared.
   * @throws NullPointerException if either argument is null.
   */
  public void getSystemMultipoleMoments(Context context, DoubleArray moments) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(
        getPointer(), context.getPointer(), moments.getPointer());
  }

  /**
   * Retrieve total dipoles from a Context into the supplied caller-owned vector array.
   *
   * @param context Context at which total dipoles are evaluated.
   * @param dipoles destination native vector array sized for the system particles; each vector
   *     contains a total dipole in elementary-charge nanometers.
   * @throws NullPointerException if either argument is null.
   */
  public void getTotalDipoles(Context context, Vec3Array dipoles) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_getTotalDipoles(
        getPointer(), context.getPointer(), dipoles.getPointer());
  }

  @Deprecated
  /**
   * Set the legacy Ewald alpha parameter.
   *
   * @param aewald Ewald separation parameter in inverse nanometers.
   * @deprecated Use {@link #setPMEParameters(double, int, int, int)}.
   */
  public void setAEwald(double aewald) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setAEwald(getPointer(), aewald);
  }

  /**
   * Set one covalent-neighbor map. Categories are Covalent12=0, Covalent13=1, Covalent14=2,
   * Covalent15=3, PolarizationCovalent11=4, PolarizationCovalent12=5,
   * PolarizationCovalent13=6, and PolarizationCovalent14=7.
   *
   * @param particle     particle index.
   * @param covalentType native covalent-map category value from 0 through 7.
   * @param covalentMap  caller-owned native integer array of related particle indices.
   * @throws NullPointerException if {@code covalentMap} is null.
   */
  public void setCovalentMap(int particle, int covalentType, IntArray covalentMap) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setCovalentMap(
        getPointer(), particle, covalentType, covalentMap.getPointer());
  }

  /**
   * Set the nonbonded cutoff distance.
   *
   * @param cutoff cutoff distance in nanometers; it has no effect with NoCutoff.
   */
  public void setCutoffDistance(double cutoff) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setCutoffDistance(getPointer(), cutoff);
  }

  /**
   * Set the target fractional force error tolerance for Ewald summation.
   *
   * @param tolerance dimensionless relative force error target.
   */
  public void setEwaldErrorTolerance(double tolerance) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(getPointer(), tolerance);
  }

  /**
   * Set coefficients for the induced-dipole extrapolation terms mu0, mu1, and so on. The number of
   * coefficients determines the number of extrapolation iterations; native code copies the values.
   *
   * @param coefficients caller-owned native double array containing the coefficients.
   * @throws NullPointerException if {@code coefficients} is null.
   */
  public void setExtrapolationCoefficients(DoubleArray coefficients) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(
        getPointer(), coefficients.getPointer());
  }

  /**
   * Set extrapolation coefficients from a native double-array handle. The native force copies the
   * values; the supplied handle remains caller-owned.
   *
   * @param coefficients valid native {@code OpenMM_DoubleArray} handle.
   */
  public void setExtrapolationCoefficients(MemorySegment coefficients) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(
        getPointer(), coefficients);
  }

  /**
   * Replace all parameters for one multipole using a Java value record. Multipole arrays are copied
   * into temporary native storage for the call.
   *
   * @param index      particle/multipole index.
   * @param parameters replacement parameters.
   * @throws NullPointerException if {@code parameters} or either multipole array is null.
   */
  public void setMultipoleParameters(int index, MultipoleParameters parameters) {
    Objects.requireNonNull(parameters, "Multipole parameters cannot be null.");
    try (DoubleArray dipole = toNative(parameters.molecularDipole());
         DoubleArray quadrupole = toNative(parameters.molecularQuadrupole())) {
      setMultipoleParameters(index, parameters.charge(), dipole, quadrupole,
          parameters.axisType(), parameters.multipoleAtomZ(), parameters.multipoleAtomX(),
          parameters.multipoleAtomY(), parameters.thole(), parameters.dampingFactor(),
          parameters.polarity());
    }
  }

  /**
   * Replace all parameters for one multipole using caller-owned native arrays.
   *
   * @param index               particle/multipole index.
   * @param charge              particle charge in elementary-charge units.
   * @param molecularDipole     native double array of three molecular-frame dipole components.
   * @param molecularQuadrupole native double array of nine molecular-frame quadrupole components.
   * @param axisType            local-frame enum value: ZThenX=0, Bisector=1, ZBisect=2, ThreeFold=3,
   *                            ZOnly=4, or NoAxisType=5.
   * @param multipoleAtomZ      first local-frame reference atom index, or -1 if unused.
   * @param multipoleAtomX      second local-frame reference atom index, or -1 if unused.
   * @param multipoleAtomY      third local-frame reference atom index, or -1 if unused.
   * @param thole               Thole damping parameter.
   * @param dampingFactor       polarization damping-factor parameter.
   * @param polarity            atomic polarizability parameter.
   */
  public void setMultipoleParameters(
      int index, double charge, DoubleArray molecularDipole, DoubleArray molecularQuadrupole,
      int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY,
      double thole, double dampingFactor, double polarity) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setMultipoleParameters(
        getPointer(), index, charge, molecularDipole.getPointer(), molecularQuadrupole.getPointer(),
        axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
  }

  /**
   * Set the maximum number of iterations for mutually induced dipoles.
   *
   * @param iterations maximum iteration count.
   */
  public void setMutualInducedMaxIterations(int iterations) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(getPointer(), iterations);
  }

  /**
   * Set the convergence target for mutually induced dipoles.
   *
   * @param epsilon target epsilon used by the iterative solver.
   */
  public void setMutualInducedTargetEpsilon(double epsilon) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(getPointer(), epsilon);
  }

  /**
   * Set the long-range nonbonded method.
   *
   * @param method native enum value: NoCutoff=0 or PME=1.
   */
  public void setNonbondedMethod(int method) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setNonbondedMethod(getPointer(), method);
  }

  /**
   * Set PME reciprocal-space parameters.
   *
   * @param alpha Ewald separation parameter in inverse nanometers; zero selects automatic values
   *              from the Ewald error tolerance.
   * @param nx    PME grid points along X, or zero for automatic selection.
   * @param ny    PME grid points along Y, or zero for automatic selection.
   * @param nz    PME grid points along Z, or zero for automatic selection.
   */
  public void setPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setPMEParameters(getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Set explicit PME grid dimensions.
   *
   * @param gridDimensions caller-owned native integer array containing X, Y, and Z grid sizes.
   * @throws NullPointerException if {@code gridDimensions} is null.
   */
  public void setPmeGridDimensions(IntArray gridDimensions) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(
        getPointer(), gridDimensions.getPointer());
  }

  /**
   * Select the polarization approximation.
   *
   * @param method native enum value: Mutual=0, Direct=1, or Extrapolated=2.
   */
  public void setPolarizationType(int method) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_setPolarizationType(getPointer(), method);
  }

  /**
   * Copy changed multipole parameters and covalent maps into an existing Context. Other force
   * properties, including PME/polarization settings and convergence options, require Context
   * reinitialization.
   *
   * @param context Context to update.
   * @throws NullPointerException if {@code context} is null.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_AmoebaMultipoleForce_updateParametersInContext(
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
        OpenMMNative.OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(getPointer()));
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
    return OpenMMNative.OpenMM_AmoebaMultipoleForce_create.makeInvoker().apply();
  }

  private PmeParameters readPmeParameters(Context context) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment alpha = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment nx = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment ny = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment nz = arena.allocate(ValueLayout.JAVA_INT);
      if (context == null) {
        OpenMMNative.OpenMM_AmoebaMultipoleForce_getPMEParameters(
            getPointer(), alpha, nx, ny, nz);
      } else {
        OpenMMNative.OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(
            getPointer(), context.getPointer(), alpha, nx, ny, nz);
      }
      return new PmeParameters(alpha.get(ValueLayout.JAVA_DOUBLE, 0),
          nx.get(ValueLayout.JAVA_INT, 0), ny.get(ValueLayout.JAVA_INT, 0),
          nz.get(ValueLayout.JAVA_INT, 0));
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
