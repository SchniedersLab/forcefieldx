// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.numerics.multipole;

import jdk.incubator.vector.DoubleVector;

import static jdk.incubator.vector.DoubleVector.SPECIES_PREFERRED;
import static jdk.incubator.vector.DoubleVector.fromArray;

/**
 * The PolarizableMultipole class defines a polarizable multipole.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PolarizableMultipoleSIMD {

  private static final double oneThird = 1.0 / 3.0;
  private static final double twoThirds = 2.0 / 3.0;

  /**
   * Partial charge.
   */
  protected DoubleVector q;
  /**
   * Dipole x-component.
   */
  protected DoubleVector dx;
  /**
   * Dipole y-component.
   */
  protected DoubleVector dy;
  /**
   * Dipole z-component.
   */
  protected DoubleVector dz;
  /**
   * Quadrupole xx-component multiplied by 1/3.
   */
  protected DoubleVector qxx;
  /**
   * Quadrupole yy-component multiplied by 1/3.
   */
  protected DoubleVector qyy;
  /**
   * Quadrupole zz-component multiplied by 1/3.
   */
  protected DoubleVector qzz;
  /**
   * Quadrupole xy-component multiplied by 2/3.
   */
  protected DoubleVector qxy;
  /**
   * Quadrupole xz-component multiplied by 2/3.
   */
  protected DoubleVector qxz;
  /**
   * Quadrupole xz-component multiplied by 2/3.
   */
  protected DoubleVector qyz;

  /**
   * Induced dipole x-component.
   */
  protected DoubleVector ux;
  /**
   * Induced dipole y-component.
   */
  protected DoubleVector uy;
  /**
   * Induced dipole z-component.
   */
  protected DoubleVector uz;
  /**
   * Induced dipole chain rule x-component.
   */
  protected DoubleVector px;
  /**
   * Induced dipole chain rule y-component.
   */
  protected DoubleVector py;
  /**
   * Induced dipole chain rule z-component.
   */
  protected DoubleVector pz;
  /**
   * Averaged induced dipole + induced dipole chain-rule x-component: sx = 0.5 * (ux + px).
   */
  protected DoubleVector sx;
  /**
   * Averaged induced dipole + induced dipole chain-rule y-component: sy = 0.5 * (uy + py).
   */
  protected DoubleVector sy;
  /**
   * Averaged induced dipole + induced dipole chain-rule z-component: sz = 0.5 * (uz + pz).
   */
  protected DoubleVector sz;

  /**
   * PolarizableMultipole constructor with zero moments.
   */
  public PolarizableMultipoleSIMD() {
  }

  /**
   * PolarizableMultipole constructor.
   *
   * @param Q   Multipoles Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u   Induced dipoles u[ux, uy, uz]
   * @param uCR Induced dipole chain-rules uCR[ux, uy, uz]
   */
  public PolarizableMultipoleSIMD(double[][] Q, double[][] u, double[][] uCR) {
    setPermanentMultipole(Q);
    setInducedDipole(u, uCR);
  }

  /**
   * Set the permanent multipole.
   * <p>
   * Note that the quadrupole trace components are multiplied by 1/3 and the
   * off-diagonal components are multiplied by 2/3.
   *
   * @param Q   Multipoles Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u   Induced dipoles u[ux, uy, uz]
   * @param uCR Induced dipole chain-rules uCR[ux, uy, uz]
   */
  public void set(double[][] Q, double[][] u, double[][] uCR) {
    setPermanentMultipole(Q);
    setInducedDipole(u, uCR);
  }

  /**
   * Set the permanent multipole.
   * <p>
   * Note that the quadrupole trace components are multiplied by 1/3 and the
   * off-diagonal components are multiplied by 2/3.
   *
   * @param Q Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   */
  public final void setPermanentMultipole(double[][] Q) {
    q = fromArray(SPECIES_PREFERRED, Q[0], 0);
    dx = fromArray(SPECIES_PREFERRED, Q[1], 0);
    dy = fromArray(SPECIES_PREFERRED, Q[2], 0);
    dz = fromArray(SPECIES_PREFERRED, Q[3], 0);
    qxx = fromArray(SPECIES_PREFERRED, Q[4], 0).mul(oneThird);
    qyy = fromArray(SPECIES_PREFERRED, Q[5], 0).mul(oneThird);
    qzz = fromArray(SPECIES_PREFERRED, Q[6], 0).mul(oneThird);
    qxy = fromArray(SPECIES_PREFERRED, Q[7], 0).mul(twoThirds);
    qxz = fromArray(SPECIES_PREFERRED, Q[8], 0).mul(twoThirds);
    qyz = fromArray(SPECIES_PREFERRED, Q[9], 0).mul(twoThirds);
  }

  /**
   * Set the induced dipole.
   *
   * @param u   Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public final void setInducedDipole(double[][] u, double[][] uCR) {
    ux = fromArray(SPECIES_PREFERRED, u[0], 0);
    uy = fromArray(SPECIES_PREFERRED, u[1], 0);
    uz = fromArray(SPECIES_PREFERRED, u[2], 0);
    px = fromArray(SPECIES_PREFERRED, uCR[0], 0);
    py = fromArray(SPECIES_PREFERRED, uCR[1], 0);
    pz = fromArray(SPECIES_PREFERRED, uCR[2], 0);
    sx = ux.add(px).mul(0.5);
    sy = uy.add(py).mul(0.5);
    sz = uz.add(pz).mul(0.5);
  }

  /**
   * Compute the scaled and averaged induced dipole.
   *
   * @param scaleInduction Induction mask scale factor.
   * @param scaleEnergy    Energy mask scale factor.
   */
  public final void applyMasks(double scaleInduction, double scaleEnergy) {
    // [Ux, Uy, Uz] resulted from induction masking rules, and we now apply the energy mask.
    // [Px, Py, Pz] resulted from energy masking rules, and we now apply the induction mask.
    sx = ux.mul(scaleEnergy).add(px.mul(scaleInduction)).mul(0.5);
    sy = uy.mul(scaleEnergy).add(py.mul(scaleInduction)).mul(0.5);
    sz = uz.mul(scaleEnergy).add(pz.mul(scaleInduction)).mul(0.5);
  }

}
