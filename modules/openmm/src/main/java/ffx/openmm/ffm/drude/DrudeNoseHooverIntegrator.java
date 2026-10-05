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

import ffx.openmm.ffm.NoseHooverIntegrator;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * Nose-Hoover integrator for systems containing Drude particles, with separate thermostats for
 * ordinary and Drude internal degrees of freedom.
 *
 * <p>The ordinary thermostat acts on ordinary particles and the center-of-mass motion of each
 * Drude pair; the second acts on the relative internal displacement of each pair and is typically
 * set to a lower temperature. The System must contain a DrudeForce so the integrator can identify
 * Drude particles.</p>
 *
 * <p>A hard-wall constraint can limit the distance from each Drude particle to its parent. The
 * native header's class overview says the default is 0.02 nm, while its accessor documentation
 * says zero is the default; setting the distance to zero omits the constraint. See the native
 * implementation/version in use when relying on the initial default.</p>
 */
public class DrudeNoseHooverIntegrator extends NoseHooverIntegrator {

  /**
   * Create a Drude Nose-Hoover integrator.
   *
   * <p>The Java parameters {@code chainLength}, {@code drudeChainLength}, and {@code numMTS} are
   * forwarded positionally, not translated by name. The native header instead defines those three
   * positions as {@code chainLength}, {@code numMTS}, and {@code numYoshidaSuzuki}: it has one
   * shared chain length, no distinct Drude chain length, and documents the Yoshida-Suzuki count as
   * 1, 3, or 5. Consequently, Java's {@code drudeChainLength} supplies the native MTS count, and
   * Java's {@code numMTS} supplies the native Yoshida-Suzuki count. The native header's default
   * Yoshida-Suzuki value is 7 despite its stated 1/3/5 constraint; this constructor has no defaults
   * and passes the supplied values unchanged.</p>
   *
   * @param stepSize         Integration time step in picoseconds.
   * @param temperature      Target temperature for ordinary degrees of freedom and Drude-pair center of
   *                         mass motion, in kelvin.
   * @param drudeTemperature Target temperature for Drude internal coordinates, in kelvin.
   * @param frequency        Main heat-bath coupling frequency in ps<sup>-1</sup>.
   * @param drudeFrequency   Drude internal-coordinate heat-bath coupling frequency in ps<sup>-1</sup>.
   * @param chainLength      Value forwarded as the native shared Nose-Hoover chain length.
   * @param drudeChainLength Legacy Java parameter name; forwarded positionally as the native MTS
   *                         chain-propagation count, not as a separate Drude chain length.
   * @param numMTS           Legacy Java parameter name; forwarded positionally as the native
   *                         Yoshida-Suzuki count. The header documents 1, 3, or 5 as valid, although its declared
   *                         default is 7.
   */
  public DrudeNoseHooverIntegrator(
      double stepSize, double temperature, double drudeTemperature,
      double frequency, double drudeFrequency, int chainLength,
      int drudeChainLength, int numMTS) {
    super(create(stepSize, temperature, drudeTemperature, frequency, drudeFrequency,
        chainLength, drudeChainLength, numMTS));
  }

  /**
   * Release the native integrator handle.
   *
   * <p>The integrator must not be used after it is destroyed.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_DrudeNoseHooverIntegrator_destroy);
  }

  /**
   * Compute the kinetic energy of Drude internal degrees of freedom.
   *
   * <p>Call only after this integrator has been attached to a live OpenMM Context. The native
   * method requires the associated Context to remain alive.</p>
   *
   * @return Drude kinetic energy in kJ/mol.
   */
  public double computeDrudeKineticEnergy() {
    return OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(getPointer());
  }

  /**
   * Compute the instantaneous temperature of Drude internal coordinates.
   *
   * <p>This is calculated from the kinetic energy of the relative internal motion of Drude pairs
   * and should remain close, on average, to the configured Drude temperature. Call only while the
   * integrator is attached to a live OpenMM Context; the native method requires that Context.</p>
   *
   * @return Instantaneous Drude internal-coordinate temperature in kelvin.
   */
  public double computeDrudeTemperature() {
    return OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(getPointer());
  }

  /**
   * Compute the instantaneous temperature of ordinary system degrees of freedom.
   *
   * <p>This includes ordinary particles and the center-of-mass motion of Drude pairs, but excludes
   * their relative internal motion. On average it should be approximately the main heat-bath
   * temperature. Call only while the integrator is attached to a live OpenMM Context; the native
   * method requires that Context.</p>
   *
   * @return Instantaneous system temperature in kelvin.
   */
  public double computeSystemTemperature() {
    return OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(getPointer());
  }

  /**
   * Compute the total kinetic energy of ordinary and Drude particles.
   *
   * <p>Call only while the integrator is attached to a live OpenMM Context; the native method
   * requires that Context to remain alive.</p>
   *
   * @return Total kinetic energy in kJ/mol.
   */
  public double computeTotalKineticEnergy() {
    return OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(getPointer());
  }

  /**
   * Get the hard-wall limit on the distance between each Drude particle and its parent.
   *
   * <p>Zero disables the constraint. The native header is inconsistent about the initial default:
   * its class overview says 0.02 nm, while this accessor's documentation says zero.</p>
   *
   * @return Maximum Drude-parent distance in nm.
   */
  public double getMaxDrudeDistance() {
    return OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(getPointer());
  }

  /**
   * Set the hard-wall limit on the distance between each Drude particle and its parent.
   *
   * <p>Setting the distance to zero omits the hard-wall constraint.</p>
   *
   * @param distance Maximum Drude-parent distance in nm, or zero to omit the constraint.
   */
  public void setMaxDrudeDistance(double distance) {
    OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(
        getPointer(), distance);
  }

  private static MemorySegment create(
      double stepSize, double temperature, double drudeTemperature,
      double frequency, double drudeFrequency, int chainLength,
      int drudeChainLength, int numMTS) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_DrudeNoseHooverIntegrator_create(
        temperature, frequency, drudeTemperature, drudeFrequency, stepSize,
        chainLength, drudeChainLength, numMTS);
  }
}
