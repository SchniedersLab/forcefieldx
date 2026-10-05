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
package ffx.openmm.ffm.drude;

import ffx.openmm.ffm.Context;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.OpenMMBooleans;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;

/**
 * Applies Drude oscillator forces: an anisotropic harmonic force between each Drude particle and
 * its parent, and a Thole-screened Coulomb interaction between selected pairs of dipoles.
 *
 * <p>A Drude particle is represented by its System particle index, its parent particle index, and
 * up to three additional particles used to define anisotropic polarizability. A screened pair is a
 * pair of Drude dipoles, each consisting of a Drude particle and its parent. Add all entries before
 * constructing a Context that uses this force. Changes made later with the parameter setters affect
 * existing Contexts only after {@link #updateParametersInContext(Context)}; that operation can
 * change numeric parameters, but cannot change particle identities or add entries.</p>
 *
 * <p>Periodic boundary conditions, when enabled, are used only for screened-pair displacements,
 * never for the Drude-parent spring.</p>
 */
public class DrudeForce extends Force {

  /**
   * Copied force-field parameters for one Drude particle.
   *
   * @param particle       System index of the Drude particle.
   * @param particle1      System index of its parent particle.
   * @param particle2      System index of the second particle defining anisotropy, or {@code -1} to
   *                       disable the {@code aniso12} scale factor.
   * @param particle3      System index of the third particle defining anisotropy, or {@code -1} to
   *                       disable the {@code aniso34} scale factor.
   * @param particle4      System index of the fourth particle defining anisotropy, or {@code -1} to
   *                       disable the {@code aniso34} scale factor.
   * @param charge         Charge on the Drude particle, in elementary charge units.
   * @param polarizability Isotropic polarizability, in nm^3.
   * @param aniso12        Scale factor for polarizability along the direction defined by
   *                       {@code particle1} and {@code particle2}.
   * @param aniso34        Scale factor for polarizability along the direction defined by
   *                       {@code particle3} and {@code particle4}.
   */
  public record ParticleParameters(
      int particle, int particle1, int particle2, int particle3, int particle4,
      double charge, double polarizability, double aniso12, double aniso34) {
  }

  /**
   * Copied parameters for one screened dipole pair.
   *
   * @param particle1 Index of the first particle in this force's screened-pair interaction.
   * @param particle2 Index of the second particle in this force's screened-pair interaction.
   * @param thole     Dimensionless Thole screening factor.
   */
  public record ScreenedPairParameters(int particle1, int particle2, double thole) {
  }

  /**
   * Create an empty Drude force.
   *
   * <p>The force owns its native handle; call {@link #destroy()} when it is no longer needed.</p>
   */
  public DrudeForce() {
    super(create());
  }

  /**
   * Add a Drude particle to this force.
   *
   * @param particle       System index of the Drude particle.
   * @param particle1      System index of the particle to which the Drude particle is attached.
   * @param particle2      System index of the second particle defining anisotropic polarizability, or
   *                       {@code -1} to ignore {@code aniso12}.
   * @param particle3      System index of the third particle defining anisotropic polarizability, or
   *                       {@code -1} to ignore {@code aniso34}.
   * @param particle4      System index of the fourth particle defining anisotropic polarizability, or
   *                       {@code -1} to ignore {@code aniso34}.
   * @param charge         Charge on the Drude particle, in elementary charge units.
   * @param polarizability Isotropic polarizability, in nm^3.
   * @param aniso12        Scale factor for polarizability along the direction defined by
   *                       {@code particle1} and {@code particle2}.
   * @param aniso34        Scale factor for polarizability along the direction defined by
   *                       {@code particle3} and {@code particle4}.
   * @return Index of the added Drude particle in this force's particle list.
   */
  public int addParticle(
      int particle, int particle1, int particle2, int particle3, int particle4,
      double charge, double polarizability, double aniso12, double aniso34) {
    return OpenMMNative.OpenMM_DrudeForce_addParticle(
        getPointer(), particle, particle1, particle2, particle3, particle4,
        charge, polarizability, aniso12, aniso34);
  }

  /**
   * Add a screened interaction between two dipoles.
   *
   * @param particle1 Index within this force of the first particle in the interaction.
   * @param particle2 Index within this force of the second particle in the interaction.
   * @param thole     Dimensionless Thole screening factor.
   * @return Index of the added screened pair in this force's screened-pair list.
   */
  public int addScreenedPair(int particle1, int particle2, double thole) {
    return OpenMMNative.OpenMM_DrudeForce_addScreenedPair(
        getPointer(), particle1, particle2, thole);
  }

  /**
   * Release the native force handle.
   *
   * <p>The force must not be used after it is destroyed.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_DrudeForce_destroy);
  }

  /**
   * Get the number of Drude particles whose force-field parameters are defined.
   *
   * @return Number of Drude particles in this force.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_DrudeForce_getNumParticles(getPointer());
  }

  /**
   * Get the number of screened dipole interactions defined.
   *
   * @return Number of screened pairs in this force.
   */
  public int getNumScreenedPairs() {
    return OpenMMNative.OpenMM_DrudeForce_getNumScreenedPairs(getPointer());
  }

  /**
   * Get a snapshot of one Drude particle's force-field parameters.
   *
   * <p>The native API writes these values through output references. This FFM method copies them
   * into a {@link ParticleParameters} record; the returned values do not refer to native or arena
   * memory and remain usable after the call completes.</p>
   *
   * @param index Index of the Drude particle in this force's particle list.
   * @return Copied particle indices and force-field parameters.
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment[] particle = new MemorySegment[5];
      for (int i = 0; i < particle.length; i++) {
        particle[i] = arena.allocate(ValueLayout.JAVA_INT);
      }
      MemorySegment charge = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment polarizability = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment aniso12 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment aniso34 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_DrudeForce_getParticleParameters(
          getPointer(), index, particle[0], particle[1], particle[2], particle[3], particle[4],
          charge, polarizability, aniso12, aniso34);
      return new ParticleParameters(
          particle[0].get(ValueLayout.JAVA_INT, 0),
          particle[1].get(ValueLayout.JAVA_INT, 0),
          particle[2].get(ValueLayout.JAVA_INT, 0),
          particle[3].get(ValueLayout.JAVA_INT, 0),
          particle[4].get(ValueLayout.JAVA_INT, 0),
          charge.get(ValueLayout.JAVA_DOUBLE, 0),
          polarizability.get(ValueLayout.JAVA_DOUBLE, 0),
          aniso12.get(ValueLayout.JAVA_DOUBLE, 0),
          aniso34.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get a snapshot of one screened pair's force-field parameters.
   *
   * <p>The native API writes these values through output references. This FFM method copies them
   * into a {@link ScreenedPairParameters} record, independent of native and temporary arena
   * memory.</p>
   *
   * @param index Index of the screened pair in this force's screened-pair list.
   * @return Copied particle indices and Thole screening factor.
   */
  public ScreenedPairParameters getScreenedPairParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment particle1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment particle2 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment thole = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_DrudeForce_getScreenedPairParameters(
          getPointer(), index, particle1, particle2, thole);
      return new ScreenedPairParameters(
          particle1.get(ValueLayout.JAVA_INT, 0),
          particle2.get(ValueLayout.JAVA_INT, 0),
          thole.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace the force-field parameters for one Drude particle.
   *
   * <p>As in the native API, this updates the force object only. Existing Contexts require
   * {@link #updateParametersInContext(Context)} to receive supported changes. Particle identities
   * cannot be changed in an existing Context.</p>
   *
   * @param index          Index of the Drude particle in this force's particle list.
   * @param particle       System index of the Drude particle.
   * @param particle1      System index of its parent particle.
   * @param particle2      System index of the second particle defining anisotropic polarizability, or
   *                       {@code -1} to ignore {@code aniso12}.
   * @param particle3      System index of the third particle defining anisotropic polarizability, or
   *                       {@code -1} to ignore {@code aniso34}.
   * @param particle4      System index of the fourth particle defining anisotropic polarizability, or
   *                       {@code -1} to ignore {@code aniso34}.
   * @param charge         Charge on the Drude particle, in elementary charge units.
   * @param polarizability Isotropic polarizability, in nm^3.
   * @param aniso12        Scale factor for polarizability along the direction defined by
   *                       {@code particle1} and {@code particle2}.
   * @param aniso34        Scale factor for polarizability along the direction defined by
   *                       {@code particle3} and {@code particle4}.
   */
  public void setParticleParameters(
      int index, int particle, int particle1, int particle2, int particle3, int particle4,
      double charge, double polarizability, double aniso12, double aniso34) {
    OpenMMNative.OpenMM_DrudeForce_setParticleParameters(
        getPointer(), index, particle, particle1, particle2, particle3, particle4,
        charge, polarizability, aniso12, aniso34);
  }

  /**
   * Replace the parameters for one screened pair.
   *
   * <p>This updates the force object; call {@link #updateParametersInContext(Context)} to update
   * an existing Context. The native update operation does not permit changing the identities of
   * particles in a Context.</p>
   *
   * @param index     Index of the screened pair in this force's screened-pair list.
   * @param particle1 Index within this force of the first particle in the interaction.
   * @param particle2 Index within this force of the second particle in the interaction.
   * @param thole     Dimensionless Thole screening factor.
   */
  public void setScreenedPairParameters(
      int index, int particle1, int particle2, double thole) {
    OpenMMNative.OpenMM_DrudeForce_setScreenedPairParameters(
        getPointer(), index, particle1, particle2, thole);
  }

  /**
   * Enable or disable periodic boundary conditions for screened-pair displacements.
   *
   * <p>This setting never applies periodic boundaries to the spring between a Drude particle and
   * its parent.</p>
   *
   * @param periodic {@code true} to use periodic boundary conditions for screened pairs.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Copy supported particle and screened-pair parameter changes into an existing Context.
   *
   * <p>First change numeric parameters with the setters, then call this method. Native OpenMM
   * permits numeric parameter changes only: this does not add entries or change particle
   * identities. If the supplied wrapper has no live native Context handle, this method does
   * nothing.</p>
   *
   * @param context Context whose force parameters should be updated.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMMNative.OpenMM_DrudeForce_updateParametersInContext(
          getPointer(), context.getPointer());
    }
  }

  /**
   * Report whether this force uses periodic boundary conditions for screened pairs.
   *
   * @return {@code true} if periodic boundary conditions are enabled.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_DrudeForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_DrudeForce_create.makeInvoker().apply();
  }
}
