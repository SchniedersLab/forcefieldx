// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.xray;

import static ffx.numerics.math.DoubleMath.sub;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * SolventBinaryFormFactor class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public final class SolventBinaryFormFactor implements FormFactor {

  private final double[] xyz = new double[3];
  private final double[] dxyz = new double[3];
  private final double probeRad;

  /**
   * Constructor for SolventBinaryFormFactor.
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @param probeRad a double.
   */
  public SolventBinaryFormFactor(Atom atom, double probeRad) {
    this(atom, probeRad, atom.getXYZ(null));
  }

  /**
   * Constructor for SolventBinaryFormFactor.
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @param probeRad a double.
   * @param xyz an array of double.
   */
  public SolventBinaryFormFactor(Atom atom, double probeRad, double[] xyz) {
    this.probeRad = probeRad;

    update(xyz);
  }

  /** {@inheritDoc} */
  @Override
  public double rho(double f, double lambda, double[] xyz) {
    sub(this.xyz, xyz, dxyz);
    return rho(f, lambda, DoubleMath.length(dxyz));
  }

  /**
   * rho
   *
   * @param f a double.
   * @param lambda a double.
   * @param ri a double.
   * @return a double.
   */
  public double rho(double f, double lambda, double ri) {
    if (ri <= probeRad) {
      return 0.0;
    } else {
      return f;
    }
  }

  /**
   * Derivatives are zero or infinite for the binary model.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public void rhoGrad(double[] xyz, double dfc, RefinementMode refinementmode) {}

  /** {@inheritDoc} */
  @Override
  public void update(double[] xyz) {
    update(xyz, 0.0);
  }

  /** {@inheritDoc} */
  @Override
  public void update(double[] xyz, double badd) {
    this.xyz[0] = xyz[0];
    this.xyz[1] = xyz[1];
    this.xyz[2] = xyz[2];
  }
}
