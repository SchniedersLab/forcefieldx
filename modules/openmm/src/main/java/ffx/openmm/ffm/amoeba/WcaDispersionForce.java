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
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.OpenMMBooleans;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.Objects;

/**
 * Implements the AMOEBA WCA dispersion interaction, commonly used with the Generalized Kirkwood
 * implicit-solvent force.
 *
 * <p>Add one radius and well depth for each System particle. Radii use nanometers and well depths
 * use kilojoules per mole. Only per-particle radius and epsilon values can be copied into an
 * existing Context; changes to force-wide water and dispersion settings require Context
 * reinitialization. The native force does not use periodic boundary conditions.</p>
 */
public class WcaDispersionForce extends Force {

  /**
   * Create an AMOEBA WCA dispersion force.
   *
   * @throws IllegalStateException if the native force cannot be created.
   */
  public WcaDispersionForce() {
    super(create());
  }

  /**
   * Add parameters for a particle.
   *
   * @param radius particle radius in nanometers.
   * @param epsilon particle well depth in kilojoules per mole.
   * @return index of the added particle.
   */
  public int addParticle(double radius, double epsilon) {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_addParticle(getPointer(), radius, epsilon);
  }

  /**
   * Destroy this force and release its owned native object.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_AmoebaWcaDispersionForce_destroy);
  }

  /**
   * Get the water density parameter used by the WCA dispersion model.
   *
   * @return water density parameter in the model's native units.
   */
  public double getAwater() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getAwater(getPointer());
  }

  /**
   * Get the dispersion offset.
   *
   * @return dispersion offset in nanometers.
   */
  public double getDispoff() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getDispoff(getPointer());
  }

  /**
   * Get the water hydrogen epsilon parameter.
   *
   * @return water hydrogen well depth in kilojoules per mole.
   */
  public double getEpsh() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getEpsh(getPointer());
  }

  /**
   * Get the water oxygen epsilon parameter.
   *
   * @return water oxygen well depth in kilojoules per mole.
   */
  public double getEpso() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getEpso(getPointer());
  }

  /**
   * Get the number of particle parameter sets.
   *
   * @return particle count.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getNumParticles(getPointer());
  }

  /**
   * Get parameters for a particle.
   *
   * @param particleIndex particle index.
   * @return copied particle radius in nanometers and well depth in kilojoules per mole.
   */
  public ParticleParameters getParticleParameters(int particleIndex) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment radius = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getParticleParameters(getPointer(),
          particleIndex, radius, epsilon);
      return new ParticleParameters(radius.get(ValueLayout.JAVA_DOUBLE, 0),
          epsilon.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the water hydrogen radius parameter.
   *
   * @return water hydrogen radius in nanometers.
   */
  public double getRminh() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getRminh(getPointer());
  }

  /**
   * Get the water oxygen radius parameter.
   *
   * @return water oxygen radius in nanometers.
   */
  public double getRmino() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getRmino(getPointer());
  }

  /**
   * Get the overlap correction factor.
   *
   * @return overlap correction factor in the native model's units.
   */
  public double getShctd() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getShctd(getPointer());
  }

  /**
   * Get the Levy parameter.
   *
   * @return Levy dispersion parameter in the native model's units.
   */
  public double getSlevy() {
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_getSlevy(getPointer());
  }

  /**
   * Set the water density parameter.
   *
   * @param awater water density parameter in the native model's units.
   */
  public void setAwater(double awater) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setAwater(getPointer(), awater);
  }

  /**
   * Set the dispersion offset.
   *
   * @param dispoff dispersion offset in nanometers.
   */
  public void setDispoff(double dispoff) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setDispoff(getPointer(), dispoff);
  }

  /**
   * Set the water hydrogen epsilon parameter.
   *
   * @param epsh water hydrogen well depth in kilojoules per mole.
   */
  public void setEpsh(double epsh) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setEpsh(getPointer(), epsh);
  }

  /**
   * Set the water oxygen epsilon parameter.
   *
   * @param epso water oxygen well depth in kilojoules per mole.
   */
  public void setEpso(double epso) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setEpso(getPointer(), epso);
  }

  /**
   * Replace parameters for a particle.
   *
   * @param particleIndex particle index.
   * @param radius particle radius in nanometers.
   * @param epsilon particle well depth in kilojoules per mole.
   */
  public void setParticleParameters(int particleIndex, double radius, double epsilon) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setParticleParameters(getPointer(),
        particleIndex, radius, epsilon);
  }

  /**
   * Set the water hydrogen radius parameter.
   *
   * @param rminh water hydrogen radius in nanometers.
   */
  public void setRminh(double rminh) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setRminh(getPointer(), rminh);
  }

  /**
   * Set the water oxygen radius parameter.
   *
   * @param rmino water oxygen radius in nanometers.
   */
  public void setRmino(double rmino) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setRmino(getPointer(), rmino);
  }

  /**
   * Set the overlap factor.
   *
   * @param shctd overlap correction factor in the native model's units.
   */
  public void setShctd(double shctd) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setShctd(getPointer(), shctd);
  }

  /**
   * Set the Levy parameter.
   *
   * @param slevy Levy dispersion parameter in the native model's units.
   */
  public void setSlevy(double slevy) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_setSlevy(getPointer(), slevy);
  }

  /**
   * Copy modified per-particle radius and epsilon values into an existing Context. Other
   * force-wide parameters are not changed in the Context and require Context reinitialization.
   *
   * @param context context to update.
   * @throws NullPointerException if {@code context} is null.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(getPointer(),
        Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Determine whether this force uses periodic boundary conditions.
   *
   * @return always false; the native WCA dispersion force does not use periodic boundaries.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Per-particle WCA dispersion parameters.
   *
   * @param radius particle radius in nanometers.
   * @param epsilon particle well depth in kilojoules per mole.
   */
  public record ParticleParameters(double radius, double epsilon) {
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_AmoebaWcaDispersionForce_create.makeInvoker().apply();
  }
}
