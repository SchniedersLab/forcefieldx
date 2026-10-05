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
 * AMOEBA implicit-solvation force using the Generalized Kirkwood/Grycuk model.
 *
 * <p>Define one parameter set for each system particle. The five-parameter form separates the
 * effective-radius and descreening radii and includes a neck correction. The number of particle
 * parameter sets must match the System particle count before Context creation. Updating an existing
 * Context copies per-particle parameters only; force-wide settings require Context reinitialization.
 * This force does not use periodic boundary conditions.</p>
 */
public class GeneralizedKirkwoodForce extends Force {

  /**
   * Generalized Kirkwood parameters for one particle.
   *
   * @param charge         particle charge in elementary-charge units.
   * @param radius         base atomic radius in nanometers.
   * @param scalingFactor  unitless factor applied to the descreen radius.
   * @param descreenRadius atomic descreening radius in nanometers.
   * @param neckFactor     unitless interstitial-neck descreening factor.
   */
  public record ParticleParameters(
      double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor) {
  }

  /**
   * Coefficients of the tanh rescaling function.
   *
   * @param b0 first tanh coefficient.
   * @param b1 second tanh coefficient.
   * @param b2 third tanh coefficient.
   */
  public record TanhParameters(double b0, double b1, double b2) {
  }

  /**
   * Create an empty Generalized Kirkwood force.
   *
   * @throws IllegalStateException if the native force cannot be created.
   */
  public GeneralizedKirkwoodForce() {
    super(create());
  }

  /**
   * Add one particle's legacy three-parameter values. Call once per System particle, in particle
   * order. This compatibility form sets the descreen radius equal to {@code radius} and sets the
   * neck factor to zero (no neck descreening).
   *
   * @param charge        particle charge in elementary-charge units.
   * @param radius        base atomic radius in nanometers.
   * @param scalingFactor unitless descreening scale.
   * @return index assigned to the added particle.
   */
  public int addParticle(double charge, double radius, double scalingFactor) {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(
        getPointer(), charge, radius, scalingFactor);
  }

  /**
   * Add one particle's complete Generalized Kirkwood parameters. The separate descreen radius
   * controls water displacement during pairwise descreening; the neck factor adds interstitial-neck
   * descreening when positive. Call once per System particle, in particle order.
   *
   * @param charge         particle charge in elementary-charge units.
   * @param radius         base atomic radius in nanometers, used in the effective-radius calculation.
   * @param scalingFactor  unitless factor applied to the descreen radius.
   * @param descreenRadius atomic radius used for descreening, in nanometers.
   * @param neckFactor     unitless interstitial-neck descreening scale.
   * @return index assigned to the added particle.
   */
  public int addParticle(
      double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor) {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(
        getPointer(), charge, radius, scalingFactor, descreenRadius, neckFactor);
  }

  /**
   * Release this force's owned native object.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_AmoebaGeneralizedKirkwoodForce_destroy);
  }

  /**
   * Get the offset used in the descreening integral.
   *
   * @return descreen offset in nanometers.
   */
  public double getDescreenOffset() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(getPointer());
  }

  /**
   * Get the offset used for the cavity contribution.
   *
   * @return dielectric offset in nanometers.
   */
  public double getDielectricOffset() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(getPointer());
  }

  /**
   * Get whether the cavity term is included.
   *
   * @return native integer flag: zero disables and nonzero enables the cavity term.
   */
  public int getIncludeCavityTerm() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(getPointer());
  }

  /**
   * Get the number of particle parameter sets.
   *
   * @return number of particles configured on this force.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(getPointer());
  }

  /**
   * Get copied particle parameters. The returned record and its scalar components are Java values
   * independent of the temporary native output storage.
   *
   * @param index particle index.
   * @return copied charge, radii, scaling factor, and neck factor.
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment[] values = new MemorySegment[5];
      for (int i = 0; i < values.length; i++) {
        values[i] = arena.allocate(ValueLayout.JAVA_DOUBLE);
      }
      OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(
          getPointer(), index, values[0], values[1], values[2], values[3], values[4]);
      return new ParticleParameters(
          values[0].get(ValueLayout.JAVA_DOUBLE, 0),
          values[1].get(ValueLayout.JAVA_DOUBLE, 0),
          values[2].get(ValueLayout.JAVA_DOUBLE, 0),
          values[3].get(ValueLayout.JAVA_DOUBLE, 0),
          values[4].get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the probe radius used for the cavity contribution.
   *
   * @return probe radius in nanometers.
   */
  public double getProbeRadius() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(getPointer());
  }

  /**
   * Get the dielectric constant of the solute.
   *
   * @return dimensionless solute dielectric constant.
   */
  public double getSoluteDielectric() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(getPointer());
  }

  /**
   * Get the dielectric constant of the solvent.
   *
   * @return dimensionless solvent dielectric constant.
   */
  public double getSolventDielectric() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(getPointer());
  }

  /**
   * Get the surface-area contribution factor.
   *
   * @return surface-area factor in kilojoules per mole per square nanometer.
   */
  public double getSurfaceAreaFactor() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(getPointer());
  }

  /**
   * Get the three coefficients of the tanh rescaling function as a copied Java record.
   *
   * @return copied tanh coefficients.
   */
  public TanhParameters getTanhParameters() {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment b0 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment b1 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment b2 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(
          getPointer(), b0, b1, b2);
      return new TanhParameters(
          b0.get(ValueLayout.JAVA_DOUBLE, 0),
          b1.get(ValueLayout.JAVA_DOUBLE, 0),
          b2.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get whether tanh rescaling is enabled.
   *
   * @return native integer flag: zero disables and nonzero enables tanh rescaling.
   */
  public int getTanhRescaling() {
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(getPointer());
  }

  /**
   * Set the offset used in the descreening integral.
   *
   * @param offset descreen offset in nanometers.
   */
  public void setDescreenOffset(double offset) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(getPointer(), offset);
  }

  /**
   * Set the offset used for the cavity contribution.
   *
   * @param offset dielectric offset in nanometers.
   */
  public void setDielectricOffset(double offset) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(getPointer(), offset);
  }

  /**
   * Set whether the cavity term is included.
   *
   * @param includeCavityTerm native integer flag; zero disables and nonzero enables the term.
   */
  public void setIncludeCavityTerm(int includeCavityTerm) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(
        getPointer(), includeCavityTerm);
  }

  /**
   * Replace all five parameters for a particle.
   *
   * @param index      particle index.
   * @param parameters copied Java parameter record to apply.
   * @throws NullPointerException if {@code parameters} is null.
   */
  public void setParticleParameters(int index, ParticleParameters parameters) {
    Objects.requireNonNull(parameters, "Particle parameters cannot be null.");
    setParticleParameters(index, parameters.charge(), parameters.radius(), parameters.scalingFactor(),
        parameters.descreenRadius(), parameters.neckFactor());
  }

  /**
   * Replace all five parameters for a particle.
   *
   * @param index          particle index.
   * @param charge         particle charge in elementary-charge units.
   * @param radius         base atomic radius in nanometers.
   * @param scalingFactor  unitless factor applied to the descreen radius.
   * @param descreenRadius atomic descreening radius in nanometers.
   * @param neckFactor     unitless interstitial-neck descreening scale.
   */
  public void setParticleParameters(
      int index, double charge, double radius, double scalingFactor,
      double descreenRadius, double neckFactor) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(
        getPointer(), index, charge, radius, scalingFactor, descreenRadius, neckFactor);
  }

  /**
   * Set the probe radius used for the cavity contribution.
   *
   * @param radius probe radius in nanometers.
   */
  public void setProbeRadius(double radius) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(getPointer(), radius);
  }

  /**
   * Set the solute dielectric constant.
   *
   * @param dielectric dimensionless solute dielectric constant.
   */
  public void setSoluteDielectric(double dielectric) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(getPointer(), dielectric);
  }

  /**
   * Set the solvent dielectric constant.
   *
   * @param dielectric dimensionless solvent dielectric constant.
   */
  public void setSolventDielectric(double dielectric) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(getPointer(), dielectric);
  }

  /**
   * Set the surface-area contribution factor.
   *
   * @param surfaceAreaFactor factor in kilojoules per mole per square nanometer.
   */
  public void setSurfaceAreaFactor(double surfaceAreaFactor) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(
        getPointer(), surfaceAreaFactor);
  }

  /**
   * Set all coefficients of the tanh rescaling function.
   *
   * @param beta0 first tanh coefficient.
   * @param beta1 second tanh coefficient.
   * @param beta2 third tanh coefficient.
   */
  public void setTanhParameters(double beta0, double beta1, double beta2) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(
        getPointer(), beta0, beta1, beta2);
  }

  /**
   * Set all coefficients of the tanh rescaling function.
   *
   * @param parameters copied Java coefficient record.
   * @throws NullPointerException if {@code parameters} is null.
   */
  public void setTanhParameters(TanhParameters parameters) {
    Objects.requireNonNull(parameters, "Tanh parameters cannot be null.");
    setTanhParameters(parameters.b0(), parameters.b1(), parameters.b2());
  }

  /**
   * Enable or disable tanh rescaling.
   *
   * @param tanhRescale native integer flag; zero disables and nonzero enables rescaling.
   */
  public void setTanhRescaling(int tanhRescale) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(
        getPointer(), tanhRescale);
  }

  /**
   * Copy changed per-particle parameters into an existing Context. Force-wide properties such as
   * dielectric constants, offsets, probe radius, tanh settings, and surface-area factor are not
   * updated by this call and require Context reinitialization.
   *
   * @param context Context to update.
   * @throws NullPointerException if {@code context} is null.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(
        getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Report whether this force uses periodic boundary conditions. The native implementation always
   * reports false.
   *
   * @return always {@code false}.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(
            getPointer()));
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_AmoebaGeneralizedKirkwoodForce_create.makeInvoker().apply();
  }
}
