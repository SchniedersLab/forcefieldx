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
package ffx.openmm.ffm;

import ffx.openmm.ffm.bindings.OpenMMNative;
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;

/**
 * Gay-Berne anisotropic nonbonded interactions between ellipsoidal particles.
 *
 * <p>Particle parameters combine a Lennard-Jones-like {@code sigma} (nm) and {@code epsilon}
 * (kJ/mol) with ellipsoid diameters {@code sx,sy,sz} (nm) and dimensionless well-depth scaling
 * factors {@code ex,ey,ez}. Lorentz-Berthelot combining rules are used for particle pairs;
 * exceptions override pair parameters. Ellipsoid axes are oriented by the positions of
 * {@code xparticle} and {@code yparticle}: {@code -1} for x denotes a sphere, and {@code -1} for
 * y denotes axial symmetry. The vector to the x-axis particle defines the ellipsoid x direction;
 * the y-axis particle defines y after its x component is removed. Particle parameter entries
 * must correspond to all system particles in system order. The force applies a Lennard-Jones-like
 * interaction based on the nearest points of the oriented ellipsoids.
 *
 * <p>The native nonbonded method values are {@code 0} ({@code NoCutoff}, the default),
 * {@code 1} ({@code CutoffNonPeriodic}), and {@code 2} ({@code CutoffPeriodic}). Cutoff and
 * switching distances are in nm. Switching is ignored with {@code NoCutoff}; when enabled for a
 * cutoff method, the switching distance must be less than the cutoff distance.
 */
public class GayBerneForce extends Force {

  /**
   * Copied particle parameters.
   *
   * @param sigma Lennard-Jones-like size parameter in nm
   * @param epsilon well depth in kJ/mol
   * @param xparticle particle whose position defines the ellipsoid x axis, or {@code -1} for a
   *     sphere
   * @param yparticle particle whose position defines the ellipsoid y axis, or {@code -1} for an
   *     axially symmetric ellipsoid
   * @param sx ellipsoid diameter along x in nm
   * @param sy ellipsoid diameter along y in nm
   * @param sz ellipsoid diameter along z in nm
   * @param ex dimensionless epsilon scale along x
   * @param ey dimensionless epsilon scale along y
   * @param ez dimensionless epsilon scale along z
   */
  public record ParticleParameters(
      double sigma, double epsilon, int xparticle, int yparticle,
      double sx, double sy, double sz, double ex, double ey, double ez) {}
  /**
   * Copied parameters for a particle-pair exception.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param sigma exception size parameter in nm
   * @param epsilon exception well depth in kJ/mol; zero omits the interaction
   */
  public record ExceptionParameters(int particle1, int particle2, double sigma, double epsilon) {}

  /** Create a Gay-Berne force; its native nonbonded method defaults to {@code NoCutoff}. */
  public GayBerneForce() {
    super(create());
  }

  /**
   * Add or replace an exception for a particle pair.
   *
   * <p>With {@code replace == false}, attempting to add a duplicate pair causes an exception;
   * with {@code true}, its existing exception is replaced. An exception with zero epsilon omits
   * that pair interaction.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param sigma exception size parameter in nm
   * @param epsilon exception well depth in kJ/mol
   * @param replace whether an existing exception for the pair is replaced
   * @return index of the added or replaced exception
   */
  public int addException(
      int particle1, int particle2, double sigma, double epsilon, boolean replace) {
    return OpenMMNative.OpenMM_GayBerneForce_addException(
        getPointer(), particle1, particle2, sigma, epsilon, OpenMMBooleans.toNative(replace));
  }

  /**
   * Add or replace an exception using the legacy integer replacement flag, forwarded unchanged
   * to the native ABI. Prefer the boolean overload when possible.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param sigma exception size parameter in nm
   * @param epsilon exception well depth in kJ/mol; zero omits the interaction
   * @param replace legacy native replacement flag
   * @return index of the added or replaced exception
   */
  public int addException(
      int particle1, int particle2, double sigma, double epsilon, int replace) {
    return OpenMMNative.OpenMM_GayBerneForce_addException(
        getPointer(), particle1, particle2, sigma, epsilon, replace);
  }

  /**
   * Add the next system particle's force-field parameters.
   *
   * <p>The Java signature intentionally takes {@code ex,ey,ez} before {@code sx,sy,sz}; the
   * native wrapper ABI receives diameters first and scale factors second. The returned record
   * orders diameters before scale factors.
   *
   * @param sigma size parameter in nm
   * @param epsilon well depth in kJ/mol
   * @param xparticle particle whose position defines the ellipsoid x axis, or {@code -1} for a
   *     sphere
   * @param yparticle particle whose position defines the ellipsoid y axis, or {@code -1} for an
   *     axially symmetric ellipsoid
   * @param ex dimensionless epsilon scale along x
   * @param ey dimensionless epsilon scale along y
   * @param ez dimensionless epsilon scale along z
   * @param sx ellipsoid diameter along x in nm
   * @param sy ellipsoid diameter along y in nm
   * @param sz ellipsoid diameter along z in nm
   * @return index of the added particle
   */
  public int addParticle(
      double sigma, double epsilon, int xparticle, int yparticle,
      double ex, double ey, double ez, double sx, double sy, double sz) {
    return OpenMMNative.OpenMM_GayBerneForce_addParticle(
        getPointer(), sigma, epsilon, xparticle, yparticle,
        sx, sy, sz, ex, ey, ez);
  }

  /** Destroy the native force. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_GayBerneForce_destroy);
  }

  /**
   * Get the cutoff distance.
   *
   * @return cutoff in nm; it has no effect when the method is {@code NoCutoff}
   */
  public double getCutoffDistance() {
    return OpenMMNative.OpenMM_GayBerneForce_getCutoffDistance(getPointer());
  }

  /**
   * Get exception parameters as copied scalar values.
   *
   * @param index index of the exception
   * @return copied particle indices, size in nm, and well depth in kJ/mol
   */
  public ExceptionParameters getExceptionParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment p1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment p2 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment sigma = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_GayBerneForce_getExceptionParameters(
          getPointer(), index, p1, p2, sigma, epsilon);
      return new ExceptionParameters(
          p1.get(ValueLayout.JAVA_INT, 0), p2.get(ValueLayout.JAVA_INT, 0),
          sigma.get(ValueLayout.JAVA_DOUBLE, 0), epsilon.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the native nonbonded-method enum value: {@code 0} NoCutoff, {@code 1}
   * CutoffNonPeriodic, or {@code 2} CutoffPeriodic.
   *
   * @return native enum value
   */
  public int getNonbondedMethod() {
    return OpenMMNative.OpenMM_GayBerneForce_getNonbondedMethod(getPointer());
  }

  /** @return number of particle-pair exceptions */
  public int getNumExceptions() {
    return OpenMMNative.OpenMM_GayBerneForce_getNumExceptions(getPointer());
  }

  /** @return number of particle parameter entries */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_GayBerneForce_getNumParticles(getPointer());
  }

  /**
   * Get a particle's parameters as copied scalar values.
   *
   * @param index particle parameter index
   * @return size (nm), well depth (kJ/mol), orientation particle indices, ellipsoid diameters
   *     (nm), then dimensionless strength factors
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment sigma = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment xparticle = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment yparticle = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment sx = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment sy = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment sz = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment ex = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment ey = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment ez = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_GayBerneForce_getParticleParameters(
          getPointer(), index, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
      return new ParticleParameters(
          sigma.get(ValueLayout.JAVA_DOUBLE, 0), epsilon.get(ValueLayout.JAVA_DOUBLE, 0),
          xparticle.get(ValueLayout.JAVA_INT, 0), yparticle.get(ValueLayout.JAVA_INT, 0),
          sx.get(ValueLayout.JAVA_DOUBLE, 0), sy.get(ValueLayout.JAVA_DOUBLE, 0),
          sz.get(ValueLayout.JAVA_DOUBLE, 0), ex.get(ValueLayout.JAVA_DOUBLE, 0),
          ey.get(ValueLayout.JAVA_DOUBLE, 0), ez.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the distance at which the switching function begins reducing the interaction.
   *
   * @return switching distance in nm; it must be less than the cutoff distance when used
   */
  public double getSwitchingDistance() {
    return OpenMMNative.OpenMM_GayBerneForce_getSwitchingDistance(getPointer());
  }

  /**
   * Get the native integer switching-function flag.
   *
   * @return native integer flag (zero is false and nonzero is true)
   */
  public int getUseSwitchingFunction() {
    return OpenMMNative.OpenMM_GayBerneForce_getUseSwitchingFunction(getPointer());
  }

  /** @return whether the native switching-function flag is nonzero */
  public boolean isUseSwitchingFunction() {
    return OpenMMBooleans.fromNative(getUseSwitchingFunction());
  }

  /**
   * Set the nonbonded cutoff distance.
   *
   * @param distance cutoff in nm; ignored with {@code NoCutoff}
   */
  public void setCutoffDistance(double distance) {
    OpenMMNative.OpenMM_GayBerneForce_setCutoffDistance(getPointer(), distance);
  }

  /**
   * Replace the parameters of an existing particle-pair exception.
   *
   * @param index exception index
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param sigma size parameter in nm
   * @param epsilon well depth in kJ/mol; zero omits the interaction
   */
  public void setExceptionParameters(
      int index, int particle1, int particle2, double sigma, double epsilon) {
    OpenMMNative.OpenMM_GayBerneForce_setExceptionParameters(
        getPointer(), index, particle1, particle2, sigma, epsilon);
  }

  /**
   * Set the native nonbonded-method enum: {@code 0} NoCutoff, {@code 1} CutoffNonPeriodic, or
   * {@code 2} CutoffPeriodic. NoCutoff is the initial native setting.
   *
   * @param method native enum value
   */
  public void setNonbondedMethod(int method) {
    OpenMMNative.OpenMM_GayBerneForce_setNonbondedMethod(getPointer(), method);
  }

  /**
   * Replace an existing particle's force-field parameters.
   *
   * <p>As in {@link #addParticle}, Java argument order is strengths {@code ex,ey,ez} then
   * diameters {@code sx,sy,sz}; the native call maps these to its own diameters-then-strengths
   * parameter order.
   *
   * @param index particle parameter index
   * @param sigma size parameter in nm
   * @param epsilon well depth in kJ/mol
   * @param xparticle particle defining the ellipsoid x axis, or {@code -1} for a sphere
   * @param yparticle particle defining the ellipsoid y axis, or {@code -1} for an axially
   *     symmetric ellipsoid
   * @param ex dimensionless epsilon scale along x
   * @param ey dimensionless epsilon scale along y
   * @param ez dimensionless epsilon scale along z
   * @param sx ellipsoid diameter along x in nm
   * @param sy ellipsoid diameter along y in nm
   * @param sz ellipsoid diameter along z in nm
   */
  public void setParticleParameters(
      int index, double sigma, double epsilon, int xparticle, int yparticle,
      double ex, double ey, double ez, double sx, double sy, double sz) {
    OpenMMNative.OpenMM_GayBerneForce_setParticleParameters(
        getPointer(), index, sigma, epsilon, xparticle, yparticle,
        sx, sy, sz, ex, ey, ez);
  }

  /**
   * Set the distance at which the switching function starts reducing the interaction.
   *
   * @param distance switching distance in nm; it must be less than the cutoff distance
   */
  public void setSwitchingDistance(double distance) {
    OpenMMNative.OpenMM_GayBerneForce_setSwitchingDistance(getPointer(), distance);
  }

  /**
   * Enable or disable the switching function. It is ignored for {@code NoCutoff}; with a cutoff
   * method, set a switching distance below the cutoff.
   *
   * @param use whether to use switching
   */
  public void setUseSwitchingFunction(boolean use) {
    OpenMMNative.OpenMM_GayBerneForce_setUseSwitchingFunction(
        getPointer(), OpenMMBooleans.toNative(use));
  }

  /**
   * Set switching using the legacy integer flag, forwarded unchanged to the native ABI. Prefer
   * the boolean overload when possible.
   *
   * @param use native integer flag
   */
  public void setUseSwitchingFunction(int use) {
    OpenMMNative.OpenMM_GayBerneForce_setUseSwitchingFunction(getPointer(), use);
  }

  /**
   * Copy supported particle and exception parameter changes into an existing context.
   *
   * <p>Only particle and exception parameters are updated. The nonbonded method, cutoff, and
   * other force settings require context reinitialization. For exceptions, only sigma and epsilon
   * may change; particle pairs are fixed. The x/y orientation-defining particle indices cannot
   * change, and no particles or exceptions can be added by this operation.
   *
   * @param context context to update; this FFM wrapper performs no update if it has no native
   *     context handle
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_GayBerneForce_updateParametersInContext(
                getPointer(), pointer));
  }

  /** @return whether the selected nonbonded method is {@code CutoffPeriodic} */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_GayBerneForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_GayBerneForce_create.makeInvoker().apply();
  }
}
