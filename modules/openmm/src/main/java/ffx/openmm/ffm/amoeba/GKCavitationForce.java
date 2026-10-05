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

import java.lang.foreign.MemorySegment;

/**
 * Compatibility façade for the legacy AMOEBA Generalized Kirkwood cavitation force.
 *
 * <p>The OpenMM AMOEBA wrapper and installed native headers do not expose this force or any of its
 * native symbols. The legacy JNA constructor printed an unsupported message and terminated the JVM
 * with {@code System.exit(-1)}; this FFM compatibility façade preserves the unsupported status
 * without terminating the process. Every instance operation, including {@link #destroy()}, throws
 * {@link UnsupportedOperationException}; no native handle is created.</p>
 */
public class GKCavitationForce extends Force {

  /**
   * Construction is unavailable because the installed OpenMM native API has no cavitation force.
   *
   * @throws UnsupportedOperationException always.
   */
  public GKCavitationForce() {
    super(unavailable());
  }

  /**
   * Legacy API entry point; no native cavitation implementation is available.
   *
   * @param radius         atomic radius, conventionally in nanometers.
   * @param surfaceTension atomic surface-tension parameter; the unavailable implementation defines
   *                       no unit contract.
   * @param isHydrogen     nonzero if the atom is hydrogen, as in the legacy API.
   * @throws UnsupportedOperationException always; no native force is available.
   */
  public void addParticle(double radius, double surfaceTension, int isHydrogen) {
    throw unavailableException();
  }

  /**
   * Legacy API entry point; no native cavitation implementation is available.
   *
   * @param index          atom index.
   * @param radius         atomic radius, conventionally in nanometers.
   * @param surfaceTension atomic surface-tension parameter; the unavailable implementation defines
   *                       no unit contract.
   * @param isHydrogen     nonzero if the atom is hydrogen, as in the legacy API.
   * @throws UnsupportedOperationException always; no native force is available.
   */
  public void setParticleParameters(int index, double radius, double surfaceTension, int isHydrogen) {
    throw unavailableException();
  }

  /**
   * Legacy API entry point; no native cavitation implementation is available.
   *
   * @param method legacy nonbonded-method value; the unavailable implementation defines no enum.
   * @throws UnsupportedOperationException always; no native force is available.
   */
  public void setNonbondedMethod(int method) {
    throw unavailableException();
  }

  /**
   * No native handle exists to destroy.
   *
   * @throws UnsupportedOperationException always; construction itself is unsupported.
   */
  @Override
  public void destroy() {
    throw unavailableException();
  }

  /**
   * Legacy API entry point; no native cavitation implementation is available.
   *
   * @param context OpenMM context; ignored because no native force exists.
   * @throws UnsupportedOperationException always; no native force is available.
   */
  public void updateParametersInContext(Context context) {
    throw unavailableException();
  }

  private static MemorySegment unavailable() {
    throw unavailableException();
  }

  private static UnsupportedOperationException unavailableException() {
    return new UnsupportedOperationException(
        "GKCavitationForce is not exposed by the installed OpenMM AMOEBA native API.");
  }
}
