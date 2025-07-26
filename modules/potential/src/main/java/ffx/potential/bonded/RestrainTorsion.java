// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.potential.bonded;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.TorsionType;

import java.io.Serial;
import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

/**
 * RestraintTorsion is a class that restrains the torsion angle defined by four atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RestrainTorsion extends BondedTerm implements LambdaInterface {

  @Serial
  private static final long serialVersionUID = 1L;

  private final Atom[] atoms;
  public final TorsionType torsionType;
  private final boolean lambdaTerm;
  private final DoubleUnaryOperator lamMapper;
  public final double units;

  private double lambda = 1.0;
  private double dEdL = 0.0;
  private double d2EdL2 = 0.0;

  /**
   * Constructor for RestrainTorsion.
   *
   * @param a1            First torsional atom.
   * @param a2            Second torsional atom.
   * @param a3            Third torsional atom.
   * @param a4            Fourth torsional atom.
   * @param torsionType   TorsionType.
   * @param lambdaEnabled True if the lambda term is enabled.
   * @param reverseLambda True if the lambda term should be reversed.
   * @param units         Units for the energy term.
   */
  public RestrainTorsion(Atom a1, Atom a2, Atom a3, Atom a4, TorsionType torsionType,
                         boolean lambdaEnabled, boolean reverseLambda, double units) {
    atoms = new Atom[]{a1, a2, a3, a4};
    this.torsionType = torsionType;
    lambdaTerm = lambdaEnabled;
    if (this.lambdaTerm) {
      if (reverseLambda) {
        lamMapper = (double l) -> 1.0 - l;
      } else {
        lamMapper = (double l) -> l;
      }
    } else {
      lamMapper = (double l) -> 1.0;
    }
    this.units = units;
  }

  @Override
  public double energy(boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
    energy = 0.0;
    value = 0.0;
    dEdL = 0.0;
    // Only compute this term if at least one atom is being used.
    if (!getUse()) {
      return energy;
    }
    var atomA = atoms[0];
    var atomB = atoms[1];
    var atomC = atoms[2];
    var atomD = atoms[3];
    var va = atomA.getXYZ();
    var vb = atomB.getXYZ();
    var vc = atomC.getXYZ();
    var vd = atomD.getXYZ();
    var vba = vb.sub(va);
    var vcb = vc.sub(vb);
    var vdc = vd.sub(vc);
    var vt = vba.X(vcb);
    var vu = vcb.X(vdc);
    var rt2 = vt.length2();
    var ru2 = vu.length2();
    var rtru2 = rt2 * ru2;
    if (rtru2 != 0.0) {
      var rr = sqrt(rtru2);
      var rcb = vcb.length();
      var cosine = vt.dot(vu) / rr;
      var sine = vcb.dot(vt.X(vu)) / (rcb * rr);
      value = toDegrees(acos(cosine));
      if (sine < 0.0) {
        value = -value;
      }
      var amp = torsionType.amplitude;
      var tsin = torsionType.sine;
      var tcos = torsionType.cosine;
      energy = amp[0] * (1.0 + cosine * tcos[0] + sine * tsin[0]);
      var dedphi = amp[0] * (cosine * tsin[0] - sine * tcos[0]);
      var cosprev = cosine;
      var sinprev = sine;
      var n = torsionType.terms;
      for (int i = 1; i < n; i++) {
        var cosn = cosine * cosprev - sine * sinprev;
        var sinn = sine * cosprev + cosine * sinprev;
        var phi = 1.0 + cosn * tcos[i] + sinn * tsin[i];
        var dphi = (1.0 + i) * (cosn * tsin[i] - sinn * tcos[i]);
        energy = energy + amp[i] * phi;
        dedphi = dedphi + amp[i] * dphi;
        cosprev = cosn;
        sinprev = sinn;
      }
      energy = units * energy * lambda;
      dEdL = units * energy;
      if (gradient || lambdaTerm) {
        dedphi = units * dedphi;
        var vca = vc.sub(va);
        var vdb = vd.sub(vb);
        var dedt = vt.X(vcb).scaleI(dedphi / (rt2 * rcb));
        var dedu = vu.X(vcb).scaleI(-dedphi / (ru2 * rcb));
        var ga = dedt.X(vcb);
        var gb = vca.X(dedt).addI(dedu.X(vdc));
        var gc = dedt.X(vba).addI(vdb.X(dedu));
        var gd = dedu.X(vcb);
        int ia = atomA.getIndex() - 1;
        int ib = atomB.getIndex() - 1;
        int ic = atomC.getIndex() - 1;
        int id = atomD.getIndex() - 1;
        if (lambdaTerm) {
          lambdaGrad.add(threadID, ia, ga);
          lambdaGrad.add(threadID, ib, gb);
          lambdaGrad.add(threadID, ic, gc);
          lambdaGrad.add(threadID, id, gd);
        }
        if (gradient) {
          grad.add(threadID, ia, ga.scaleI(lambda));
          grad.add(threadID, ib, gb.scaleI(lambda));
          grad.add(threadID, ic, gc.scaleI(lambda));
          grad.add(threadID, id, gd.scaleI(lambda));
        }
      }
    }

    return energy;
  }

  @Override
  public double getLambda() {
    return lambda;
  }

  public Atom[] getAtoms() {
    return Arrays.copyOf(atoms, atoms.length);
  }

  @Override
  public Atom getAtom(int index) {
    return atoms[index];
  }

  @Override
  public boolean applyLambda() {
    return lambdaTerm;
  }

  @Override
  public void setLambda(double lambda) {
    this.lambda = lamMapper.applyAsDouble(lambda);
  }

  @Override
  public double getd2EdL2() {
    return d2EdL2;
  }

  @Override
  public double getdEdL() {
    return dEdL;
  }

  public double mapLambda(double lambda) {
    return lamMapper.applyAsDouble(lambda);
  }

  @Override
  public void getdEdXdL(double[] gradient) {
    // The chain rule term is at least supposedly zero.
  }

  @Override
  public String toString() {
    return format(" Restrain-Torsion %s, Angle: %.3f, Energy: %.3f", torsionType, value, energy);
  }
}
