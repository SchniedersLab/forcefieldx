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

/**
 * A {@link LangevinMiddleIntegrator} retained under the historical name {@code LangevinIntegrator}. The OpenMM
 * header states it is identical to {@code LangevinMiddleIntegrator} and merely subclasses it, with the two classes
 * existing for historical reasons; the header's {@code LangevinMiddleIntegrator} documentation says LFMiddle tends to
 * sample more accurately than the discretization used in {@code LangevinIntegrator}.
 *
 * <p>Wrapper/header mismatch: the native C wrapper declares {@code OpenMM_LangevinIntegrator_create}, but this class
 * inherits {@link LangevinMiddleIntegrator}'s constructor, which calls {@code OpenMM_LangevinMiddleIntegrator_create},
 * exactly as the JNA subclass does. The OpenMM C++ and C constructors take (temperature, frictionCoeff, stepSize),
 * whereas this constructor takes (dt, temp, gamma).</p>
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.LangevinIntegrator} is not.</p>
 */
public class LangevinIntegrator extends LangevinMiddleIntegrator {

  /**
   * Create a Langevin integrator by delegating to {@link LangevinMiddleIntegrator#LangevinMiddleIntegrator(double,
   * double, double)}.
   *
   * @param dt    step size with which to integrate the system, in ps.
   * @param temp  temperature of the heat bath, in K.
   * @param gamma friction coefficient coupling the system to the heat bath, in 1/ps.
   */
  public LangevinIntegrator(double dt, double temp, double gamma) {
    super(dt, temp, gamma);
  }
}
