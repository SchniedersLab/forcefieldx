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
package ffx.numerics.multipole;

/**
 * The PolarizableMultipole class defines a polarizable multipole.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PolarizableMultipole {

  private static final double oneThird = 1.0 / 3.0;
  private static final double twoThirds = 2.0 / 3.0;

  /**
   * Partial charge.
   */
  protected double q;
  /**
   * Dipole x-component.
   */
  protected double dx;
  /**
   * Dipole y-component.
   */
  protected double dy;
  /**
   * Dipole z-component.
   */
  protected double dz;
  /**
   * Quadrupole xx-component multiplied by 1/3.
   */
  protected double qxx;
  /**
   * Quadrupole yy-component multiplied by 1/3.
   */
  protected double qyy;
  /**
   * Quadrupole zz-component multiplied by 1/3.
   */
  protected double qzz;
  /**
   * Quadrupole xy-component multiplied by 2/3.
   */
  protected double qxy;
  /**
   * Quadrupole xz-component multiplied by 2/3.
   */
  protected double qxz;
  /**
   * Quadrupole xz-component multiplied by 2/3.
   */
  protected double qyz;

  /**
   * Induced dipole x-component.
   */
  protected double ux;
  /**
   * Induced dipole y-component.
   */
  protected double uy;
  /**
   * Induced dipole z-component.
   */
  protected double uz;
  /**
   * Induced dipole chain rule x-component.
   */
  protected double px;
  /**
   * Induced dipole chain rule y-component.
   */
  protected double py;
  /**
   * Induced dipole chain rule z-component.
   */
  protected double pz;
  /**
   * Averaged induced dipole + induced dipole chain-rule x-component: sx = 0.5 * (ux + px).
   */
  protected double sx;
  /**
   * Averaged induced dipole + induced dipole chain-rule y-component: sy = 0.5 * (uy + py).
   */
  protected double sy;
  /**
   * Averaged induced dipole + induced dipole chain-rule z-component: sz = 0.5 * (uz + pz).
   */
  protected double sz;
  protected double Z;


  /**
   * PolarizableMultipole constructor with zero moments.
   */
  public PolarizableMultipole() {
  }

  /**
   * PolarizableMultipole constructor.
   *
   * @param Q   Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u   Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public PolarizableMultipole(double[] Q, double[] u, double[] uCR) {
    setPermanentMultipole(Q);
    setInducedDipole(u, uCR);
  }

  /**
   * PolarizableMultipole constructor.
   *
   * @param Q   Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u   Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public PolarizableMultipole(double[] Q, double[] u, double[] uCR, double Z) {
    setPermanentMultipole(Q);
    setInducedDipole(u, uCR);
    this.Z = Z;
    this.q = -Z + q;
  }

  /**
   * Set the permanent multipole.
   * <p>
   * Note that the quadrupole trace components are multiplied by 1/3 and the
   * off-diagonal components are multiplied by 2/3.
   *
   * @param Q   Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u   Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public void set(double[] Q, double[] u, double[] uCR) {
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
  public final void setPermanentMultipole(double[] Q) {
    q = Q[0];
    dx =  Q[1];
    dy =  Q[2];
    dz =  Q[3];
    qxx = Q[4] * oneThird;
    qyy = Q[5] * oneThird;
    qzz = Q[6] * oneThird;
    qxy = Q[7] * twoThirds;
    qxz = Q[8] * twoThirds;
    qyz = Q[9] * twoThirds;
  }

  /**
   * Set the induced dipole.
   *
   * @param u   Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public final void setInducedDipole(double[] u, double[] uCR) {
    ux = u[0];
    uy = u[1];
    uz = u[2];
    px = uCR[0];
    py = uCR[1];
    pz = uCR[2];
    sx = 0.5 * (ux + px);
    sy = 0.5 * (uy + py);
    sz = 0.5 * (uz + pz);
  }

  /**
   * Clear the induced dipoles.
   */
  public void clearInducedDipoles() {
    ux = 0.0;
    uy = 0.0;
    uz = 0.0;
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    sx = 0.0;
    sy = 0.0;
    sz = 0.0;
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
    sx = 0.5 * (ux * scaleEnergy + px * scaleInduction);
    sy = 0.5 * (uy * scaleEnergy + py * scaleInduction);
    sz = 0.5 * (uz * scaleEnergy + pz * scaleInduction);
  }

}
