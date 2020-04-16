// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The MultipoleTensorGlobal class computes derivatives of 1/|<b>r</b>| via recursion to arbitrary
 * order for Cartesian multipoles in either a global frame.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 *     Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 *     computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 *     Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @since 1.0
 */
public class MultipoleTensorGlobal extends MultipoleTensor {

  /**
   * Constructor for MultipoleTensorGlobal.
   *
   * @param operator a OPERATOR object.
   * @param order a int.
   * @param aewald a double.
   */
  public MultipoleTensorGlobal(OPERATOR operator, int order, double aewald) {
    super(operator, COORDINATES.GLOBAL, order, aewald);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Meaningful only for QI.
   */
  @Override
  public double getd2EdZ2() {
    return 0.0;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Meaningful only for QI.
   */
  @Override
  public double getdEdZ() {
    return 0.0;
  }

  /** {@inheritDoc} */
  @Override
  public double multipoleEnergy(double[] Fi, double[] Ti, double[] Tk) {
    multipoleIField();
    double energy = dotMultipoleK();
    multipoleKField();
    energy += dotMultipoleI();
    // Torques
    multipoleKTorque(Tk);
    multipoleKField();
    multipoleITorque(Ti);
    // Forces
    multipoleIdX();
    Fi[0] = -dotMultipoleK();
    multipoleIdY();
    Fi[1] = -dotMultipoleK();
    multipoleIdZ();
    Fi[2] = -dotMultipoleK();

    return energy;
  }

  /** {@inheritDoc} */
  @Override
  public double polarizationEnergy(
      double scaleField,
      double scaleEnergy,
      double scaleMutual,
      double[] Fi,
      double[] Ti,
      double[] Tk) {

    // Find the potential, field, etc at k due to the induced dipole i.
    inducedIField();
    // Energy of multipole k in the field of induced dipole i.
    double energy = scaleEnergy * dotMultipoleK();

    // Get the induced-induced portion of the force.
    Fi[0] = -0.5 * scaleMutual * (pxk * E200 + pyk * E110 + pzk * E101);
    Fi[1] = -0.5 * scaleMutual * (pxk * E110 + pyk * E020 + pzk * E011);
    Fi[2] = -0.5 * scaleMutual * (pxk * E101 + pyk * E011 + pzk * E002);

    // Find the potential, field, etc at i due to the induced dipole k.
    inducedKField();
    // Energy of multipole i in the field of induced dipole k.
    energy += scaleEnergy * dotMultipoleI();

    // Get the induced-induced portion of the force.
    Fi[0] += 0.5 * scaleMutual * (pxi * E200 + pyi * E110 + pzi * E101);
    Fi[1] += 0.5 * scaleMutual * (pxi * E110 + pyi * E020 + pzi * E011);
    Fi[2] += 0.5 * scaleMutual * (pxi * E101 + pyi * E011 + pzi * E002);

    /*
     Apply scale factors directly to induced dipole components for
     efficiency and convenience in computing remaining force terms and
     torques.
    */
    scaleInduced(scaleField, scaleEnergy);

    // Find the potential, field, etc at k due to (ind + indCR) at i.
    inducedIFieldCR();
    // Torque on multipole k.
    multipoleKTorque(Tk);

    // Find the potential, field, etc at i due to (ind + indCR) at k.
    inducedKFieldCR();
    // Torque on multipole i.
    multipoleITorque(Ti);

    // Forces
    inducedIdX();
    Fi[0] -= dotMultipoleK();
    inducedIdY();
    Fi[1] -= dotMultipoleK();
    inducedIdZ();
    Fi[2] -= dotMultipoleK();

    inducedKdX();
    Fi[0] -= dotMultipoleI();
    inducedKdY();
    Fi[1] -= dotMultipoleI();
    inducedKdZ();
    Fi[2] -= dotMultipoleI();

    return energy;
  }

  /** {@inheritDoc} */
  @Override
  protected boolean setR(double[] r, double lambdaFunction) {
    if (r[0] == rprev[0] && r[1] == rprev[1] && r[2] == rprev[2] && lambdaFunction == rprev[3]) {
      return true;
    }
    x = r[0];
    y = r[1];
    z = r[2] + lambdaFunction;
    r2 = (x * x + y * y + z * z);
    R = sqrt(r2);
    return false;
  }

  /** {@inheritDoc} */
  @Override
  protected double Tlmnj(
      final int l, final int m, final int n, final int j, final double[] r, final double[] T000) {
    if (m == 0 && n == 0) {
      if (l > 1) {
        return r[0] * Tlmnj(l - 1, 0, 0, j + 1, r, T000)
            + (l - 1) * Tlmnj(l - 2, 0, 0, j + 1, r, T000);
      } else if (l == 1) { // l == 1, d/dx is done.
        return r[0] * Tlmnj(0, 0, 0, j + 1, r, T000);
      } else { // l = m = n = 0. Recursion is done.
        return T000[j];
      }
    } else if (n == 0) { // m >= 1
      if (m > 1) {
        return r[1] * Tlmnj(l, m - 1, 0, j + 1, r, T000)
            + (m - 1) * Tlmnj(l, m - 2, 0, j + 1, r, T000);
      }
      return r[1] * Tlmnj(l, 0, 0, j + 1, r, T000);
    } else { // n >= 1
      if (n > 1) {
        return r[2] * Tlmnj(l, m, n - 1, j + 1, r, T000)
            + (n - 1) * Tlmnj(l, m, n - 2, j + 1, r, T000);
      }
      return r[2] * Tlmnj(l, m, 0, j + 1, r, T000);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>This method is a driver to collect elements of the Cartesian multipole tensor given the
   * recursion relationships implemented by the method "Tlmnj", which can be called directly to get
   * a single tensor element. It does not store intermediate values of the recursion, causing it to
   * scale O(order^8). For order = 5, this approach is a factor of 10 slower than recursion.
   */
  @Override
  protected void noStorageRecursion(double[] r, double[] tensor) {
    setR(r);
    source(T000);
    // 1/r
    tensor[0] = T000[0];
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    for (int l = 1; l <= order; l++) {
      tensor[ti(l, 0, 0)] = Tlmnj(l, 0, 0, 0, r, T000);
    }
    // Find (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
    for (int l = 0; l <= o1; l++) {
      for (int m = 1; m <= order - l; m++) {
        tensor[ti(l, m, 0)] = Tlmnj(l, m, 0, 0, r, T000);
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
    for (int l = 0; l <= o1; l++) {
      for (int m = 0; m <= o1 - l; m++) {
        for (int n = 1; n <= order - l - m; n++) {
          tensor[ti(l, m, n)] = Tlmnj(l, m, n, 0, r, T000);
        }
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  protected void recursion(double[] r, double[] tensor) {
    setR(r);
    source(work);
    tensor[0] = work[0];
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    // Any (d/dx) term can be formed as
    // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
    // All intermediate terms are indexed as l*il + m*im + n*in + j;
    double current;
    double previous = work[1];
    // Store the l=1 tensor T100 (d/dx)
    tensor[ti(1, 0, 0)] = x * previous;
    // Starting the loop at l=2 avoids an if statement.
    for (int l = 2; l < o1; l++) {
      // Initial condition for the inner loop is formation of T100(l-1).
      // Starting the inner loop at a=2 avoid an if statement.
      // T100(l-1) = x * T000(l)
      current = x * work[l];
      int iw = il + l - 1;
      work[iw] = current;
      for (int a = 1; a < l - 1; a++) {
        // T200(l-2) = x * T100(l-1) + (2 - 1) * T000(l-1)
        // T300(l-3) = x * T200(l-2) + (3 - 1) * T100(l-2)
        // ...
        // T(l-1)001 = x * T(l-2)002 + (l - 2) * T(l-3)002
        current = x * current + a * work[iw - il];
        iw += il - 1;
        work[iw] = current;
      }
      // Store the Tl00 tensor (d/dx)^l
      // Tl00 = x * [[ T(l-1)001 ]] + (l - 1) * T(l-2)001
      tensor[ti(l, 0, 0)] = x * current + (l - 1) * previous;
      previous = current;
    }
    // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
    // Any (d/dy) term can be formed as:
    // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
    for (int l = 0; l < order; l++) {
      // Store the m=1 tensor (d/dx)^l *(d/dy)
      // Tl10 = y * Tl001
      previous = work[l * il + 1];
      tensor[ti(l, 1, 0)] = y * previous;
      for (int m = 2; m + l < o1; m++) {
        // Tl10(m-1) = y * Tl00m;
        int iw = l * il + m;
        current = y * work[iw];
        iw += im - 1;
        work[iw] = current;
        for (int a = 1; a < m - 1; a++) {
          // Tl20(m-2) = Y * Tl10(m-1) + (2 - 1) * T100(m-1)
          // Tl30(m-3) = Y * Tl20(m-2) + (3 - 1) * Tl10(m-2)
          // ...
          // Tl(m-1)01 = Y * Tl(m-2)02 + (m - 2) * T(m-3)02
          current = y * current + a * work[iw - im];
          iw += im - 1;
          work[iw] = current;
        }
        // Store the tensor (d/dx)^l * (d/dy)^m
        // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
        tensor[ti(l, m, 0)] = y * current + (m - 1) * previous;
        previous = current;
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
    // Any (d/dz) term can be formed as:
    // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m + l < order; m++) {
        // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
        // Tlmn = z * Tlm01
        final int lm = m + l;
        final int lilmim = l * il + m * im;
        previous = work[lilmim + 1];
        tensor[ti(l, m, 1)] = z * previous;
        for (int n = 2; lm + n < o1; n++) {
          // Tlm1(n-1) = z * Tlm0n;
          int iw = lilmim + n;
          current = z * work[iw];
          iw += in - 1;
          work[iw] = current;
          final int n1 = n - 1;
          for (int a = 1; a < n1; a++) {
            // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
            // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
            // ...
            // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
            current = z * current + a * work[iw - in];
            iw += in - 1;
            work[iw] = current;
          }
          // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
          // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
          tensor[ti(l, m, n)] = z * current + n1 * previous;
          previous = current;
        }
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>This function is a driver to collect elements of the Cartesian multipole tensor. Collecting
   * all tensors scales slightly better than O(order^4).
   *
   * <p>For a multipole expansion truncated at quadrupole order, for example, up to order 5 is
   * needed for energy gradients. The number of terms this requires is binomial(5 + 3, 3) or 8! /
   * (5! * 3!), which is 56.
   *
   * <p>The packing of the tensor elements for order = 1<br>
   * tensor[0] = 1/|r| <br>
   * tensor[1] = -x/|r|^3 <br>
   * tensor[2] = -y/|r|^3 <br>
   * tensor[3] = -z/|r|^3 <br>
   *
   * <p>
   *
   * @since 1.0
   */
  @Override
  protected String codeTensorRecursion(final double[] r, final double[] tensor) {
    setR(r);
    source(work);
    StringBuilder sb = new StringBuilder();
    tensor[0] = work[0];
    if (work[0] > 0) {
      sb.append(format("%s = %s;\n", rlmn(0, 0, 0), term(0, 0, 0, 0)));
    }
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    // Any (d/dx) term can be formed as
    // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
    // All intermediate terms are indexed as l*il + m*im + n*in + j;
    double current;
    double previous = work[1];
    // Store the l=1 tensor T100 (d/dx)
    tensor[ti(1, 0, 0)] = x * previous;
    sb.append(format("%s = x * %s;\n", rlmn(1, 0, 0), term(0, 0, 0, 1)));
    // Starting the loop at l=2 avoids an if statement.
    for (int l = 2; l < o1; l++) {
      // Initial condition for the inner loop is formation of T100(l-1).
      // Starting the inner loop at a=2 avoid an if statement.
      // T100(l-1) = x * T000(l)
      current = x * work[l];
      int iw = il + l - 1;
      work[iw] = current;
      sb.append(format("double %s = x * %s;\n", term(1, 0, 0, l - 1), term(0, 0, 0, l)));
      for (int a = 1; a < l - 1; a++) {
        // T200(l-2) = x * T100(l-1) + (2 - 1) * T000(l-1)
        // T300(l-3) = x * T200(l-2) + (3 - 1) * T100(l-2)
        // ...
        // T(l-1)001 = x * T(l-2)002 + (l - 2) * T(l-3)002
        current = x * current + a * work[iw - il];
        iw += il - 1;
        work[iw] = current;
        if (a > 1) {
          sb.append(
              format(
                  "double %s = x * %s + %d * %s;\n",
                  term(a + 1, 0, 0, l - a - 1), term(a, 0, 0, l - a), a, term(a - 1, 0, 0, l - a)));
        } else {
          sb.append(
              format(
                  "double %s = x * %s + %s;\n",
                  term(a + 1, 0, 0, l - a - 1), term(a, 0, 0, l - a), term(a - 1, 0, 0, l - a)));
        }
      }
      // Store the Tl00 tensor (d/dx)^l
      // Tl00 = x * [[ T(l-1)001 ]] + (l - 1) * T(l-2)001
      tensor[ti(l, 0, 0)] = x * current + (l - 1) * previous;
      previous = current;
      if (l > 2) {
        sb.append(
            format(
                "%s = x * %s + %d * %s;\n",
                rlmn(l, 0, 0), term(l - 1, 0, 0, 1), (l - 1), term(l - 2, 0, 0, 1)));
      } else {
        sb.append(
            format(
                "%s = x * %s + %s;\n", rlmn(l, 0, 0), term(l - 1, 0, 0, 1), term(l - 2, 0, 0, 1)));
      }
    }
    // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
    // Any (d/dy) term can be formed as:
    // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
    for (int l = 0; l < order; l++) {
      // Store the m=1 tensor (d/dx)^l *(d/dy)
      // Tl10 = y * Tl001
      previous = work[l * il + 1];
      tensor[ti(l, 1, 0)] = y * previous;
      sb.append(format("%s = y * %s;\n", rlmn(l, 1, 0), term(l, 0, 0, 1)));
      for (int m = 2; m + l < o1; m++) {
        // Tl10(m-1) = y * Tl00m;
        int iw = l * il + m;
        current = y * work[iw];
        iw += im - 1;
        work[iw] = current;
        sb.append(format("double %s = y * %s;\n", term(l, 1, 0, m - 1), term(l, 0, 0, m)));
        for (int a = 1; a < m - 1; a++) {
          // Tl20(m-2) = Y * Tl10(m-1) + (2 - 1) * T100(m-1)
          // Tl30(m-3) = Y * Tl20(m-2) + (3 - 1) * Tl10(m-2)
          // ...
          // Tl(m-1)01 = Y * Tl(m-2)02 + (m - 2) * T(m-3)02
          current = y * current + a * work[iw - im];
          iw += im - 1;
          work[iw] = current;
          if (a > 1) {
            sb.append(
                format(
                    "double %s = y * %s + %d * %s;\n",
                    term(l, a + 1, 0, m - a - 1),
                    term(l, a, 0, m - a),
                    a,
                    term(l, a - 1, 0, m - a)));
          } else {
            sb.append(
                format(
                    "double %s = y * %s + %s;\n",
                    term(l, a + 1, 0, m - a - 1), term(l, a, 0, m - a), term(l, a - 1, 0, m - a)));
          }
        }
        // Store the tensor (d/dx)^l * (d/dy)^m
        // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
        tensor[ti(l, m, 0)] = y * current + (m - 1) * previous;
        previous = current;
        if (m > 2) {
          sb.append(
              format(
                  "%s = y * %s + %d * %s;\n",
                  rlmn(l, m, 0), term(l, m - 1, 0, 1), (m - 1), term(l, m - 2, 0, 1)));
        } else {
          sb.append(
              format(
                  "%s = y * %s + %s;\n",
                  rlmn(l, m, 0), term(l, m - 1, 0, 1), term(l, m - 2, 0, 1)));
        }
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
    // Any (d/dz) term can be formed as:
    // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m + l < order; m++) {
        // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
        // Tlmn = z * Tlm01
        final int lm = m + l;
        final int lilmim = l * il + m * im;
        previous = work[lilmim + 1];
        tensor[ti(l, m, 1)] = z * previous;
        sb.append(format("%s = z * %s;\n", rlmn(l, m, 1), term(l, m, 0, 1)));
        for (int n = 2; lm + n < o1; n++) {
          // Tlm1(n-1) = z * Tlm0n;
          int iw = lilmim + n;
          current = z * work[iw];
          iw += in - 1;
          work[iw] = current;
          sb.append(format("double %s = z * %s;\n", term(l, m, 1, n - 1), term(l, m, 0, n)));
          final int n1 = n - 1;
          for (int a = 1; a < n1; a++) {
            // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
            // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
            // ...
            // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
            current = z * current + a * work[iw - in];
            iw += in - 1;
            work[iw] = current;
            if (a > 1) {
              sb.append(
                  format(
                      "double %s = z * %s + %d * %s;\n",
                      term(l, m, a + 1, n - a - 1),
                      term(l, m, a, n - a),
                      a,
                      term(l, m, a - 1, n - a)));
            } else {
              sb.append(
                  format(
                      "double %s = z * %s + %s;\n",
                      term(l, m, a + 1, n - a - 1),
                      term(l, m, a, n - a),
                      term(l, m, a - 1, n - a)));
            }
          }
          // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
          // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
          tensor[ti(l, m, n)] = z * current + n1 * previous;
          previous = current;
          if (n > 2) {
            sb.append(
                format(
                    "%s = z * %s + %d * %s;\n",
                    rlmn(l, m, n), term(l, m, n - 1, 1), (n - 1), term(l, m, n - 2, 1)));
          } else {
            sb.append(
                format(
                    "%s = z * %s + %s;\n",
                    rlmn(l, m, n), term(l, m, n - 1, 1), term(l, m, n - 2, 1)));
          }
        }
      }
    }
    return sb.toString();
  }

  /** {@inheritDoc} */
  @Override
  protected void qiToGlobal(double[] Fi, double[] Ti, double[] Tk) {
    /* intentional no-op */
  }

  /**
   * {@inheritDoc}
   *
   * <p>Hard coded computation of all Cartesian multipole tensors up to 4th order, in the global
   * frame, which is sufficient for quadrupole-induced dipole forces.
   */
  @Override
  protected void order4() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    double term0004 = work[4];
    R000 = term0000;
    R100 = x * term0001;
    double term1001 = x * term0002;
    R200 = x * term1001 + term0001;
    double term1002 = x * term0003;
    double term2001 = x * term1002 + term0002;
    R300 = x * term2001 + 2 * term1001;
    double term1003 = x * term0004;
    double term2002 = x * term1003 + term0003;
    double term3001 = x * term2002 + 2 * term1002;
    R400 = x * term3001 + 3 * term2001;
    R010 = y * term0001;
    double term0101 = y * term0002;
    R020 = y * term0101 + term0001;
    double term0102 = y * term0003;
    double term0201 = y * term0102 + term0002;
    R030 = y * term0201 + 2 * term0101;
    double term0103 = y * term0004;
    double term0202 = y * term0103 + term0003;
    double term0301 = y * term0202 + 2 * term0102;
    R040 = y * term0301 + 3 * term0201;
    R110 = y * term1001;
    double term1101 = y * term1002;
    R120 = y * term1101 + term1001;
    double term1102 = y * term1003;
    double term1201 = y * term1102 + term1002;
    R130 = y * term1201 + 2 * term1101;
    R210 = y * term2001;
    double term2101 = y * term2002;
    R220 = y * term2101 + term2001;
    R310 = y * term3001;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = z * term0011 + term0001;
    double term0012 = z * term0003;
    double term0021 = z * term0012 + term0002;
    R003 = z * term0021 + 2 * term0011;
    double term0013 = z * term0004;
    double term0022 = z * term0013 + term0003;
    double term0031 = z * term0022 + 2 * term0012;
    R004 = z * term0031 + 3 * term0021;
    R011 = z * term0101;
    double term0111 = z * term0102;
    R012 = z * term0111 + term0101;
    double term0112 = z * term0103;
    double term0121 = z * term0112 + term0102;
    R013 = z * term0121 + 2 * term0111;
    R021 = z * term0201;
    double term0211 = z * term0202;
    R022 = z * term0211 + term0201;
    R031 = z * term0301;
    R101 = z * term1001;
    double term1011 = z * term1002;
    R102 = z * term1011 + term1001;
    double term1012 = z * term1003;
    double term1021 = z * term1012 + term1002;
    R103 = z * term1021 + 2 * term1011;
    R111 = z * term1101;
    double term1111 = z * term1102;
    R112 = z * term1111 + term1101;
    R121 = z * term1201;
    R201 = z * term2001;
    double term2011 = z * term2002;
    R202 = z * term2011 + term2001;
    R211 = z * term2101;
    R301 = z * term3001;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Hard coded computation of all Cartesian multipole tensors up to 5th order, in the global
   * frame, which is sufficient for quadrupole-quadrupole forces.
   */
  @Override
  protected void order5() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    double term0004 = work[4];
    double term0005 = work[5];
    R000 = term0000;
    R100 = x * term0001;
    double term1001 = x * term0002;
    R200 = x * term1001 + term0001;
    double term1002 = x * term0003;
    double term2001 = x * term1002 + term0002;
    R300 = x * term2001 + 2 * term1001;
    double term1003 = x * term0004;
    double term2002 = x * term1003 + term0003;
    double term3001 = x * term2002 + 2 * term1002;
    R400 = x * term3001 + 3 * term2001;
    double term1004 = x * term0005;
    double term2003 = x * term1004 + term0004;
    double term3002 = x * term2003 + 2 * term1003;
    double term4001 = x * term3002 + 3 * term2002;
    R500 = x * term4001 + 4 * term3001;
    R010 = y * term0001;
    double term0101 = y * term0002;
    R020 = y * term0101 + term0001;
    double term0102 = y * term0003;
    double term0201 = y * term0102 + term0002;
    R030 = y * term0201 + 2 * term0101;
    double term0103 = y * term0004;
    double term0202 = y * term0103 + term0003;
    double term0301 = y * term0202 + 2 * term0102;
    R040 = y * term0301 + 3 * term0201;
    double term0104 = y * term0005;
    double term0203 = y * term0104 + term0004;
    double term0302 = y * term0203 + 2 * term0103;
    double term0401 = y * term0302 + 3 * term0202;
    R050 = y * term0401 + 4 * term0301;
    R110 = y * term1001;
    double term1101 = y * term1002;
    R120 = y * term1101 + term1001;
    double term1102 = y * term1003;
    double term1201 = y * term1102 + term1002;
    R130 = y * term1201 + 2 * term1101;
    double term1103 = y * term1004;
    double term1202 = y * term1103 + term1003;
    double term1301 = y * term1202 + 2 * term1102;
    R140 = y * term1301 + 3 * term1201;
    R210 = y * term2001;
    double term2101 = y * term2002;
    R220 = y * term2101 + term2001;
    double term2102 = y * term2003;
    double term2201 = y * term2102 + term2002;
    R230 = y * term2201 + 2 * term2101;
    R310 = y * term3001;
    double term3101 = y * term3002;
    R320 = y * term3101 + term3001;
    R410 = y * term4001;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = z * term0011 + term0001;
    double term0012 = z * term0003;
    double term0021 = z * term0012 + term0002;
    R003 = z * term0021 + 2 * term0011;
    double term0013 = z * term0004;
    double term0022 = z * term0013 + term0003;
    double term0031 = z * term0022 + 2 * term0012;
    R004 = z * term0031 + 3 * term0021;
    double term0014 = z * term0005;
    double term0023 = z * term0014 + term0004;
    double term0032 = z * term0023 + 2 * term0013;
    double term0041 = z * term0032 + 3 * term0022;
    R005 = z * term0041 + 4 * term0031;
    R011 = z * term0101;
    double term0111 = z * term0102;
    R012 = z * term0111 + term0101;
    double term0112 = z * term0103;
    double term0121 = z * term0112 + term0102;
    R013 = z * term0121 + 2 * term0111;
    double term0113 = z * term0104;
    double term0122 = z * term0113 + term0103;
    double term0131 = z * term0122 + 2 * term0112;
    R014 = z * term0131 + 3 * term0121;
    R021 = z * term0201;
    double term0211 = z * term0202;
    R022 = z * term0211 + term0201;
    double term0212 = z * term0203;
    double term0221 = z * term0212 + term0202;
    R023 = z * term0221 + 2 * term0211;
    R031 = z * term0301;
    double term0311 = z * term0302;
    R032 = z * term0311 + term0301;
    R041 = z * term0401;
    R101 = z * term1001;
    double term1011 = z * term1002;
    R102 = z * term1011 + term1001;
    double term1012 = z * term1003;
    double term1021 = z * term1012 + term1002;
    R103 = z * term1021 + 2 * term1011;
    double term1013 = z * term1004;
    double term1022 = z * term1013 + term1003;
    double term1031 = z * term1022 + 2 * term1012;
    R104 = z * term1031 + 3 * term1021;
    R111 = z * term1101;
    double term1111 = z * term1102;
    R112 = z * term1111 + term1101;
    double term1112 = z * term1103;
    double term1121 = z * term1112 + term1102;
    R113 = z * term1121 + 2 * term1111;
    R121 = z * term1201;
    double term1211 = z * term1202;
    R122 = z * term1211 + term1201;
    R131 = z * term1301;
    R201 = z * term2001;
    double term2011 = z * term2002;
    R202 = z * term2011 + term2001;
    double term2012 = z * term2003;
    double term2021 = z * term2012 + term2002;
    R203 = z * term2021 + 2 * term2011;
    R211 = z * term2101;
    double term2111 = z * term2102;
    R212 = z * term2111 + term2101;
    R221 = z * term2201;
    R301 = z * term3001;
    double term3011 = z * term3002;
    R302 = z * term3011 + term3001;
    R311 = z * term3101;
    R401 = z * term4001;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Hard coded computation of all Cartesian multipole tensors up to 5th order, in the global
   * frame, which is sufficient for quadrupole-quadrupole forces and orthogonal space sampling.
   */
  @Override
  protected void order6() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    double term0004 = work[4];
    double term0005 = work[5];
    double term0006 = work[6];
    R000 = term0000;
    R100 = x * term0001;
    double term1001 = x * term0002;
    R200 = x * term1001 + term0001;
    double term1002 = x * term0003;
    double term2001 = x * term1002 + term0002;
    R300 = x * term2001 + 2 * term1001;
    double term1003 = x * term0004;
    double term2002 = x * term1003 + term0003;
    double term3001 = x * term2002 + 2 * term1002;
    R400 = x * term3001 + 3 * term2001;
    double term1004 = x * term0005;
    double term2003 = x * term1004 + term0004;
    double term3002 = x * term2003 + 2 * term1003;
    double term4001 = x * term3002 + 3 * term2002;
    R500 = x * term4001 + 4 * term3001;
    double term1005 = x * term0006;
    double term2004 = x * term1005 + term0005;
    double term3003 = x * term2004 + 2 * term1004;
    double term4002 = x * term3003 + 3 * term2003;
    double term5001 = x * term4002 + 4 * term3002;
    R600 = x * term5001 + 5 * term4001;
    R010 = y * term0001;
    double term0101 = y * term0002;
    R020 = y * term0101 + term0001;
    double term0102 = y * term0003;
    double term0201 = y * term0102 + term0002;
    R030 = y * term0201 + 2 * term0101;
    double term0103 = y * term0004;
    double term0202 = y * term0103 + term0003;
    double term0301 = y * term0202 + 2 * term0102;
    R040 = y * term0301 + 3 * term0201;
    double term0104 = y * term0005;
    double term0203 = y * term0104 + term0004;
    double term0302 = y * term0203 + 2 * term0103;
    double term0401 = y * term0302 + 3 * term0202;
    R050 = y * term0401 + 4 * term0301;
    double term0105 = y * term0006;
    double term0204 = y * term0105 + term0005;
    double term0303 = y * term0204 + 2 * term0104;
    double term0402 = y * term0303 + 3 * term0203;
    double term0501 = y * term0402 + 4 * term0302;
    R060 = y * term0501 + 5 * term0401;
    R110 = y * term1001;
    double term1101 = y * term1002;
    R120 = y * term1101 + term1001;
    double term1102 = y * term1003;
    double term1201 = y * term1102 + term1002;
    R130 = y * term1201 + 2 * term1101;
    double term1103 = y * term1004;
    double term1202 = y * term1103 + term1003;
    double term1301 = y * term1202 + 2 * term1102;
    R140 = y * term1301 + 3 * term1201;
    double term1104 = y * term1005;
    double term1203 = y * term1104 + term1004;
    double term1302 = y * term1203 + 2 * term1103;
    double term1401 = y * term1302 + 3 * term1202;
    R150 = y * term1401 + 4 * term1301;
    R210 = y * term2001;
    double term2101 = y * term2002;
    R220 = y * term2101 + term2001;
    double term2102 = y * term2003;
    double term2201 = y * term2102 + term2002;
    R230 = y * term2201 + 2 * term2101;
    double term2103 = y * term2004;
    double term2202 = y * term2103 + term2003;
    double term2301 = y * term2202 + 2 * term2102;
    R240 = y * term2301 + 3 * term2201;
    R310 = y * term3001;
    double term3101 = y * term3002;
    R320 = y * term3101 + term3001;
    double term3102 = y * term3003;
    double term3201 = y * term3102 + term3002;
    R330 = y * term3201 + 2 * term3101;
    R410 = y * term4001;
    double term4101 = y * term4002;
    R420 = y * term4101 + term4001;
    R510 = y * term5001;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = z * term0011 + term0001;
    double term0012 = z * term0003;
    double term0021 = z * term0012 + term0002;
    R003 = z * term0021 + 2 * term0011;
    double term0013 = z * term0004;
    double term0022 = z * term0013 + term0003;
    double term0031 = z * term0022 + 2 * term0012;
    R004 = z * term0031 + 3 * term0021;
    double term0014 = z * term0005;
    double term0023 = z * term0014 + term0004;
    double term0032 = z * term0023 + 2 * term0013;
    double term0041 = z * term0032 + 3 * term0022;
    R005 = z * term0041 + 4 * term0031;
    double term0015 = z * term0006;
    double term0024 = z * term0015 + term0005;
    double term0033 = z * term0024 + 2 * term0014;
    double term0042 = z * term0033 + 3 * term0023;
    double term0051 = z * term0042 + 4 * term0032;
    R006 = z * term0051 + 5 * term0041;
    R011 = z * term0101;
    double term0111 = z * term0102;
    R012 = z * term0111 + term0101;
    double term0112 = z * term0103;
    double term0121 = z * term0112 + term0102;
    R013 = z * term0121 + 2 * term0111;
    double term0113 = z * term0104;
    double term0122 = z * term0113 + term0103;
    double term0131 = z * term0122 + 2 * term0112;
    R014 = z * term0131 + 3 * term0121;
    double term0114 = z * term0105;
    double term0123 = z * term0114 + term0104;
    double term0132 = z * term0123 + 2 * term0113;
    double term0141 = z * term0132 + 3 * term0122;
    R015 = z * term0141 + 4 * term0131;
    R021 = z * term0201;
    double term0211 = z * term0202;
    R022 = z * term0211 + term0201;
    double term0212 = z * term0203;
    double term0221 = z * term0212 + term0202;
    R023 = z * term0221 + 2 * term0211;
    double term0213 = z * term0204;
    double term0222 = z * term0213 + term0203;
    double term0231 = z * term0222 + 2 * term0212;
    R024 = z * term0231 + 3 * term0221;
    R031 = z * term0301;
    double term0311 = z * term0302;
    R032 = z * term0311 + term0301;
    double term0312 = z * term0303;
    double term0321 = z * term0312 + term0302;
    R033 = z * term0321 + 2 * term0311;
    R041 = z * term0401;
    double term0411 = z * term0402;
    R042 = z * term0411 + term0401;
    R051 = z * term0501;
    R101 = z * term1001;
    double term1011 = z * term1002;
    R102 = z * term1011 + term1001;
    double term1012 = z * term1003;
    double term1021 = z * term1012 + term1002;
    R103 = z * term1021 + 2 * term1011;
    double term1013 = z * term1004;
    double term1022 = z * term1013 + term1003;
    double term1031 = z * term1022 + 2 * term1012;
    R104 = z * term1031 + 3 * term1021;
    double term1014 = z * term1005;
    double term1023 = z * term1014 + term1004;
    double term1032 = z * term1023 + 2 * term1013;
    double term1041 = z * term1032 + 3 * term1022;
    R105 = z * term1041 + 4 * term1031;
    R111 = z * term1101;
    double term1111 = z * term1102;
    R112 = z * term1111 + term1101;
    double term1112 = z * term1103;
    double term1121 = z * term1112 + term1102;
    R113 = z * term1121 + 2 * term1111;
    double term1113 = z * term1104;
    double term1122 = z * term1113 + term1103;
    double term1131 = z * term1122 + 2 * term1112;
    R114 = z * term1131 + 3 * term1121;
    R121 = z * term1201;
    double term1211 = z * term1202;
    R122 = z * term1211 + term1201;
    double term1212 = z * term1203;
    double term1221 = z * term1212 + term1202;
    R123 = z * term1221 + 2 * term1211;
    R131 = z * term1301;
    double term1311 = z * term1302;
    R132 = z * term1311 + term1301;
    R141 = z * term1401;
    R201 = z * term2001;
    double term2011 = z * term2002;
    R202 = z * term2011 + term2001;
    double term2012 = z * term2003;
    double term2021 = z * term2012 + term2002;
    R203 = z * term2021 + 2 * term2011;
    double term2013 = z * term2004;
    double term2022 = z * term2013 + term2003;
    double term2031 = z * term2022 + 2 * term2012;
    R204 = z * term2031 + 3 * term2021;
    R211 = z * term2101;
    double term2111 = z * term2102;
    R212 = z * term2111 + term2101;
    double term2112 = z * term2103;
    double term2121 = z * term2112 + term2102;
    R213 = z * term2121 + 2 * term2111;
    R221 = z * term2201;
    double term2211 = z * term2202;
    R222 = z * term2211 + term2201;
    R231 = z * term2301;
    R301 = z * term3001;
    double term3011 = z * term3002;
    R302 = z * term3011 + term3001;
    double term3012 = z * term3003;
    double term3021 = z * term3012 + term3002;
    R303 = z * term3021 + 2 * term3011;
    R311 = z * term3101;
    double term3111 = z * term3102;
    R312 = z * term3111 + term3101;
    R321 = z * term3201;
    R401 = z * term4001;
    double term4011 = z * term4002;
    R402 = z * term4011 + term4001;
    R411 = z * term4101;
    R501 = z * term5001;
  }

  /** {@inheritDoc} */
  @Override
  protected void multipoleIField() {
    double term000 = qi * R000;
    term000 -= dxi * R100;
    term000 -= dyi * R010;
    term000 -= dzi * R001;
    term000 += qxxi * R200;
    term000 += qyyi * R020;
    term000 += qzzi * R002;
    term000 += qxyi * R110;
    term000 += qxzi * R101;
    term000 += qyzi * R011;
    E000 = term000;
    double term100 = qi * R100;
    term100 -= dxi * R200;
    term100 -= dyi * R110;
    term100 -= dzi * R101;
    term100 += qxxi * R300;
    term100 += qyyi * R120;
    term100 += qzzi * R102;
    term100 += qxyi * R210;
    term100 += qxzi * R201;
    term100 += qyzi * R111;
    E100 = term100;
    double term010 = qi * R010;
    term010 -= dxi * R110;
    term010 -= dyi * R020;
    term010 -= dzi * R011;
    term010 += qxxi * R210;
    term010 += qyyi * R030;
    term010 += qzzi * R012;
    term010 += qxyi * R120;
    term010 += qxzi * R111;
    term010 += qyzi * R021;
    E010 = term010;
    double term001 = qi * R001;
    term001 -= dxi * R101;
    term001 -= dyi * R011;
    term001 -= dzi * R002;
    term001 += qxxi * R201;
    term001 += qyyi * R021;
    term001 += qzzi * R003;
    term001 += qxyi * R111;
    term001 += qxzi * R102;
    term001 += qyzi * R012;
    E001 = term001;
    double term200 = qi * R200;
    term200 -= dxi * R300;
    term200 -= dyi * R210;
    term200 -= dzi * R201;
    term200 += qxxi * R400;
    term200 += qyyi * R220;
    term200 += qzzi * R202;
    term200 += qxyi * R310;
    term200 += qxzi * R301;
    term200 += qyzi * R211;
    E200 = term200;
    double term020 = qi * R020;
    term020 -= dxi * R120;
    term020 -= dyi * R030;
    term020 -= dzi * R021;
    term020 += qxxi * R220;
    term020 += qyyi * R040;
    term020 += qzzi * R022;
    term020 += qxyi * R130;
    term020 += qxzi * R121;
    term020 += qyzi * R031;
    E020 = term020;
    double term002 = qi * R002;
    term002 -= dxi * R102;
    term002 -= dyi * R012;
    term002 -= dzi * R003;
    term002 += qxxi * R202;
    term002 += qyyi * R022;
    term002 += qzzi * R004;
    term002 += qxyi * R112;
    term002 += qxzi * R103;
    term002 += qyzi * R013;
    E002 = term002;
    double term110 = qi * R110;
    term110 -= dxi * R210;
    term110 -= dyi * R120;
    term110 -= dzi * R111;
    term110 += qxxi * R310;
    term110 += qyyi * R130;
    term110 += qzzi * R112;
    term110 += qxyi * R220;
    term110 += qxzi * R211;
    term110 += qyzi * R121;
    E110 = term110;
    double term101 = qi * R101;
    term101 -= dxi * R201;
    term101 -= dyi * R111;
    term101 -= dzi * R102;
    term101 += qxxi * R301;
    term101 += qyyi * R121;
    term101 += qzzi * R103;
    term101 += qxyi * R211;
    term101 += qxzi * R202;
    term101 += qyzi * R112;
    E101 = term101;
    double term011 = qi * R011;
    term011 -= dxi * R111;
    term011 -= dyi * R021;
    term011 -= dzi * R012;
    term011 += qxxi * R211;
    term011 += qyyi * R031;
    term011 += qzzi * R013;
    term011 += qxyi * R121;
    term011 += qxzi * R112;
    term011 += qyzi * R022;
    E011 = term011;
  }

  /** {@inheritDoc} */
  @Override
  protected void multipoleKField() {
    double term000 = 0.0;
    term000 += qk * R000;
    term000 += dxk * R100;
    term000 += dyk * R010;
    term000 += dzk * R001;
    term000 += qxxk * R200;
    term000 += qyyk * R020;
    term000 += qzzk * R002;
    term000 += qxyk * R110;
    term000 += qxzk * R101;
    term000 += qyzk * R011;
    E000 = term000;
    double term100 = 0.0;
    term100 += qk * R100;
    term100 += dxk * R200;
    term100 += dyk * R110;
    term100 += dzk * R101;
    term100 += qxxk * R300;
    term100 += qyyk * R120;
    term100 += qzzk * R102;
    term100 += qxyk * R210;
    term100 += qxzk * R201;
    term100 += qyzk * R111;
    E100 = term100;
    double term010 = 0.0;
    term010 += qk * R010;
    term010 += dxk * R110;
    term010 += dyk * R020;
    term010 += dzk * R011;
    term010 += qxxk * R210;
    term010 += qyyk * R030;
    term010 += qzzk * R012;
    term010 += qxyk * R120;
    term010 += qxzk * R111;
    term010 += qyzk * R021;
    E010 = term010;
    double term001 = 0.0;
    term001 += qk * R001;
    term001 += dxk * R101;
    term001 += dyk * R011;
    term001 += dzk * R002;
    term001 += qxxk * R201;
    term001 += qyyk * R021;
    term001 += qzzk * R003;
    term001 += qxyk * R111;
    term001 += qxzk * R102;
    term001 += qyzk * R012;
    E001 = term001;
    double term200 = 0.0;
    term200 += qk * R200;
    term200 += dxk * R300;
    term200 += dyk * R210;
    term200 += dzk * R201;
    term200 += qxxk * R400;
    term200 += qyyk * R220;
    term200 += qzzk * R202;
    term200 += qxyk * R310;
    term200 += qxzk * R301;
    term200 += qyzk * R211;
    E200 = term200;
    double term020 = 0.0;
    term020 += qk * R020;
    term020 += dxk * R120;
    term020 += dyk * R030;
    term020 += dzk * R021;
    term020 += qxxk * R220;
    term020 += qyyk * R040;
    term020 += qzzk * R022;
    term020 += qxyk * R130;
    term020 += qxzk * R121;
    term020 += qyzk * R031;
    E020 = term020;
    double term002 = 0.0;
    term002 += qk * R002;
    term002 += dxk * R102;
    term002 += dyk * R012;
    term002 += dzk * R003;
    term002 += qxxk * R202;
    term002 += qyyk * R022;
    term002 += qzzk * R004;
    term002 += qxyk * R112;
    term002 += qxzk * R103;
    term002 += qyzk * R013;
    E002 = term002;
    double term110 = 0.0;
    term110 += qk * R110;
    term110 += dxk * R210;
    term110 += dyk * R120;
    term110 += dzk * R111;
    term110 += qxxk * R310;
    term110 += qyyk * R130;
    term110 += qzzk * R112;
    term110 += qxyk * R220;
    term110 += qxzk * R211;
    term110 += qyzk * R121;
    E110 = term110;
    double term101 = 0.0;
    term101 += qk * R101;
    term101 += dxk * R201;
    term101 += dyk * R111;
    term101 += dzk * R102;
    term101 += qxxk * R301;
    term101 += qyyk * R121;
    term101 += qzzk * R103;
    term101 += qxyk * R211;
    term101 += qxzk * R202;
    term101 += qyzk * R112;
    E101 = term101;
    double term011 = 0.0;
    term011 += qk * R011;
    term011 += dxk * R111;
    term011 += dyk * R021;
    term011 += dzk * R012;
    term011 += qxxk * R211;
    term011 += qyyk * R031;
    term011 += qzzk * R013;
    term011 += qxyk * R121;
    term011 += qxzk * R112;
    term011 += qyzk * R022;
    E011 = term011;
  }

  /** {@inheritDoc} */
  @Override
  protected void multipoleIdX() {
    double term100 = 0.0;
    term100 += qi * R100;
    term100 -= dxi * R200;
    term100 -= dyi * R110;
    term100 -= dzi * R101;
    term100 += qxxi * R300;
    term100 += qyyi * R120;
    term100 += qzzi * R102;
    term100 += qxyi * R210;
    term100 += qxzi * R201;
    term100 += qyzi * R111;
    E000 = term100;
    double term200 = 0.0;
    term200 += qi * R200;
    term200 -= dxi * R300;
    term200 -= dyi * R210;
    term200 -= dzi * R201;
    term200 += qxxi * R400;
    term200 += qyyi * R220;
    term200 += qzzi * R202;
    term200 += qxyi * R310;
    term200 += qxzi * R301;
    term200 += qyzi * R211;
    E100 = term200;
    double term110 = 0.0;
    term110 += qi * R110;
    term110 -= dxi * R210;
    term110 -= dyi * R120;
    term110 -= dzi * R111;
    term110 += qxxi * R310;
    term110 += qyyi * R130;
    term110 += qzzi * R112;
    term110 += qxyi * R220;
    term110 += qxzi * R211;
    term110 += qyzi * R121;
    E010 = term110;
    double term101 = 0.0;
    term101 += qi * R101;
    term101 -= dxi * R201;
    term101 -= dyi * R111;
    term101 -= dzi * R102;
    term101 += qxxi * R301;
    term101 += qyyi * R121;
    term101 += qzzi * R103;
    term101 += qxyi * R211;
    term101 += qxzi * R202;
    term101 += qyzi * R112;
    E001 = term101;
    double term300 = 0.0;
    term300 += qi * R300;
    term300 -= dxi * R400;
    term300 -= dyi * R310;
    term300 -= dzi * R301;
    term300 += qxxi * R500;
    term300 += qyyi * R320;
    term300 += qzzi * R302;
    term300 += qxyi * R410;
    term300 += qxzi * R401;
    term300 += qyzi * R311;
    E200 = term300;
    double term120 = 0.0;
    term120 += qi * R120;
    term120 -= dxi * R220;
    term120 -= dyi * R130;
    term120 -= dzi * R121;
    term120 += qxxi * R320;
    term120 += qyyi * R140;
    term120 += qzzi * R122;
    term120 += qxyi * R230;
    term120 += qxzi * R221;
    term120 += qyzi * R131;
    E020 = term120;
    double term102 = 0.0;
    term102 += qi * R102;
    term102 -= dxi * R202;
    term102 -= dyi * R112;
    term102 -= dzi * R103;
    term102 += qxxi * R302;
    term102 += qyyi * R122;
    term102 += qzzi * R104;
    term102 += qxyi * R212;
    term102 += qxzi * R203;
    term102 += qyzi * R113;
    E002 = term102;
    double term210 = 0.0;
    term210 += qi * R210;
    term210 -= dxi * R310;
    term210 -= dyi * R220;
    term210 -= dzi * R211;
    term210 += qxxi * R410;
    term210 += qyyi * R230;
    term210 += qzzi * R212;
    term210 += qxyi * R320;
    term210 += qxzi * R311;
    term210 += qyzi * R221;
    E110 = term210;
    double term201 = 0.0;
    term201 += qi * R201;
    term201 -= dxi * R301;
    term201 -= dyi * R211;
    term201 -= dzi * R202;
    term201 += qxxi * R401;
    term201 += qyyi * R221;
    term201 += qzzi * R203;
    term201 += qxyi * R311;
    term201 += qxzi * R302;
    term201 += qyzi * R212;
    E101 = term201;
    double term111 = 0.0;
    term111 += qi * R111;
    term111 -= dxi * R211;
    term111 -= dyi * R121;
    term111 -= dzi * R112;
    term111 += qxxi * R311;
    term111 += qyyi * R131;
    term111 += qzzi * R113;
    term111 += qxyi * R221;
    term111 += qxzi * R212;
    term111 += qyzi * R122;
    E011 = term111;
  }

  /** {@inheritDoc} */
  @Override
  protected void multipoleIdY() {
    double term010 = 0.0;
    term010 += qi * R010;
    term010 -= dxi * R110;
    term010 -= dyi * R020;
    term010 -= dzi * R011;
    term010 += qxxi * R210;
    term010 += qyyi * R030;
    term010 += qzzi * R012;
    term010 += qxyi * R120;
    term010 += qxzi * R111;
    term010 += qyzi * R021;
    E000 = term010;
    double term110 = 0.0;
    term110 += qi * R110;
    term110 -= dxi * R210;
    term110 -= dyi * R120;
    term110 -= dzi * R111;
    term110 += qxxi * R310;
    term110 += qyyi * R130;
    term110 += qzzi * R112;
    term110 += qxyi * R220;
    term110 += qxzi * R211;
    term110 += qyzi * R121;
    E100 = term110;
    double term020 = 0.0;
    term020 += qi * R020;
    term020 -= dxi * R120;
    term020 -= dyi * R030;
    term020 -= dzi * R021;
    term020 += qxxi * R220;
    term020 += qyyi * R040;
    term020 += qzzi * R022;
    term020 += qxyi * R130;
    term020 += qxzi * R121;
    term020 += qyzi * R031;
    E010 = term020;
    double term011 = 0.0;
    term011 += qi * R011;
    term011 -= dxi * R111;
    term011 -= dyi * R021;
    term011 -= dzi * R012;
    term011 += qxxi * R211;
    term011 += qyyi * R031;
    term011 += qzzi * R013;
    term011 += qxyi * R121;
    term011 += qxzi * R112;
    term011 += qyzi * R022;
    E001 = term011;
    double term210 = 0.0;
    term210 += qi * R210;
    term210 -= dxi * R310;
    term210 -= dyi * R220;
    term210 -= dzi * R211;
    term210 += qxxi * R410;
    term210 += qyyi * R230;
    term210 += qzzi * R212;
    term210 += qxyi * R320;
    term210 += qxzi * R311;
    term210 += qyzi * R221;
    E200 = term210;
    double term030 = 0.0;
    term030 += qi * R030;
    term030 -= dxi * R130;
    term030 -= dyi * R040;
    term030 -= dzi * R031;
    term030 += qxxi * R230;
    term030 += qyyi * R050;
    term030 += qzzi * R032;
    term030 += qxyi * R140;
    term030 += qxzi * R131;
    term030 += qyzi * R041;
    E020 = term030;
    double term012 = 0.0;
    term012 += qi * R012;
    term012 -= dxi * R112;
    term012 -= dyi * R022;
    term012 -= dzi * R013;
    term012 += qxxi * R212;
    term012 += qyyi * R032;
    term012 += qzzi * R014;
    term012 += qxyi * R122;
    term012 += qxzi * R113;
    term012 += qyzi * R023;
    E002 = term012;
    double term120 = 0.0;
    term120 += qi * R120;
    term120 -= dxi * R220;
    term120 -= dyi * R130;
    term120 -= dzi * R121;
    term120 += qxxi * R320;
    term120 += qyyi * R140;
    term120 += qzzi * R122;
    term120 += qxyi * R230;
    term120 += qxzi * R221;
    term120 += qyzi * R131;
    E110 = term120;
    double term111 = 0.0;
    term111 += qi * R111;
    term111 -= dxi * R211;
    term111 -= dyi * R121;
    term111 -= dzi * R112;
    term111 += qxxi * R311;
    term111 += qyyi * R131;
    term111 += qzzi * R113;
    term111 += qxyi * R221;
    term111 += qxzi * R212;
    term111 += qyzi * R122;
    E101 = term111;
    double term021 = 0.0;
    term021 += qi * R021;
    term021 -= dxi * R121;
    term021 -= dyi * R031;
    term021 -= dzi * R022;
    term021 += qxxi * R221;
    term021 += qyyi * R041;
    term021 += qzzi * R023;
    term021 += qxyi * R131;
    term021 += qxzi * R122;
    term021 += qyzi * R032;
    E011 = term021;
  }

  /** {@inheritDoc} */
  @Override
  protected void multipoleIdZ() {
    double term001 = 0.0;
    term001 += qi * R001;
    term001 -= dxi * R101;
    term001 -= dyi * R011;
    term001 -= dzi * R002;
    term001 += qxxi * R201;
    term001 += qyyi * R021;
    term001 += qzzi * R003;
    term001 += qxyi * R111;
    term001 += qxzi * R102;
    term001 += qyzi * R012;
    E000 = term001;
    double term101 = 0.0;
    term101 += qi * R101;
    term101 -= dxi * R201;
    term101 -= dyi * R111;
    term101 -= dzi * R102;
    term101 += qxxi * R301;
    term101 += qyyi * R121;
    term101 += qzzi * R103;
    term101 += qxyi * R211;
    term101 += qxzi * R202;
    term101 += qyzi * R112;
    E100 = term101;
    double term011 = 0.0;
    term011 += qi * R011;
    term011 -= dxi * R111;
    term011 -= dyi * R021;
    term011 -= dzi * R012;
    term011 += qxxi * R211;
    term011 += qyyi * R031;
    term011 += qzzi * R013;
    term011 += qxyi * R121;
    term011 += qxzi * R112;
    term011 += qyzi * R022;
    E010 = term011;
    double term002 = 0.0;
    term002 += qi * R002;
    term002 -= dxi * R102;
    term002 -= dyi * R012;
    term002 -= dzi * R003;
    term002 += qxxi * R202;
    term002 += qyyi * R022;
    term002 += qzzi * R004;
    term002 += qxyi * R112;
    term002 += qxzi * R103;
    term002 += qyzi * R013;
    E001 = term002;
    double term201 = 0.0;
    term201 += qi * R201;
    term201 -= dxi * R301;
    term201 -= dyi * R211;
    term201 -= dzi * R202;
    term201 += qxxi * R401;
    term201 += qyyi * R221;
    term201 += qzzi * R203;
    term201 += qxyi * R311;
    term201 += qxzi * R302;
    term201 += qyzi * R212;
    E200 = term201;
    double term021 = 0.0;
    term021 += qi * R021;
    term021 -= dxi * R121;
    term021 -= dyi * R031;
    term021 -= dzi * R022;
    term021 += qxxi * R221;
    term021 += qyyi * R041;
    term021 += qzzi * R023;
    term021 += qxyi * R131;
    term021 += qxzi * R122;
    term021 += qyzi * R032;
    E020 = term021;
    double term003 = 0.0;
    term003 += qi * R003;
    term003 -= dxi * R103;
    term003 -= dyi * R013;
    term003 -= dzi * R004;
    term003 += qxxi * R203;
    term003 += qyyi * R023;
    term003 += qzzi * R005;
    term003 += qxyi * R113;
    term003 += qxzi * R104;
    term003 += qyzi * R014;
    E002 = term003;
    double term111 = 0.0;
    term111 += qi * R111;
    term111 -= dxi * R211;
    term111 -= dyi * R121;
    term111 -= dzi * R112;
    term111 += qxxi * R311;
    term111 += qyyi * R131;
    term111 += qzzi * R113;
    term111 += qxyi * R221;
    term111 += qxzi * R212;
    term111 += qyzi * R122;
    E110 = term111;
    double term102 = 0.0;
    term102 += qi * R102;
    term102 -= dxi * R202;
    term102 -= dyi * R112;
    term102 -= dzi * R103;
    term102 += qxxi * R302;
    term102 += qyyi * R122;
    term102 += qzzi * R104;
    term102 += qxyi * R212;
    term102 += qxzi * R203;
    term102 += qyzi * R113;
    E101 = term102;
    double term012 = 0.0;
    term012 += qi * R012;
    term012 -= dxi * R112;
    term012 -= dyi * R022;
    term012 -= dzi * R013;
    term012 += qxxi * R212;
    term012 += qyyi * R032;
    term012 += qzzi * R014;
    term012 += qxyi * R122;
    term012 += qxzi * R113;
    term012 += qyzi * R023;
    E011 = term012;
  }

  /** {@inheritDoc} */
  @Override
  protected void multipoleIdZ2() {
    /* intentional no-op */
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedIField() {
    double term000 = -uxi * R100;
    term000 -= uyi * R010;
    term000 -= uzi * R001;
    E000 = term000;
    double term100 = -uxi * R200;
    term100 -= uyi * R110;
    term100 -= uzi * R101;
    E100 = term100;
    double term010 = -uxi * R110;
    term010 -= uyi * R020;
    term010 -= uzi * R011;
    E010 = term010;
    double term001 = -uxi * R101;
    term001 -= uyi * R011;
    term001 -= uzi * R002;
    E001 = term001;
    double term200 = -uxi * R300;
    term200 -= uyi * R210;
    term200 -= uzi * R201;
    E200 = term200;
    double term020 = -uxi * R120;
    term020 -= uyi * R030;
    term020 -= uzi * R021;
    E020 = term020;
    double term002 = -uxi * R102;
    term002 -= uyi * R012;
    term002 -= uzi * R003;
    E002 = term002;
    double term110 = -uxi * R210;
    term110 -= uyi * R120;
    term110 -= uzi * R111;
    E110 = term110;
    double term101 = -uxi * R201;
    term101 -= uyi * R111;
    term101 -= uzi * R102;
    E101 = term101;
    double term011 = -uxi * R111;
    term011 -= uyi * R021;
    term011 -= uzi * R012;
    E011 = term011;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedKField() {
    double term000 = uxk * R100;
    term000 += uyk * R010;
    term000 += uzk * R001;
    E000 = term000;
    double term100 = uxk * R200;
    term100 += uyk * R110;
    term100 += uzk * R101;
    E100 = term100;
    double term010 = uxk * R110;
    term010 += uyk * R020;
    term010 += uzk * R011;
    E010 = term010;
    double term001 = uxk * R101;
    term001 += uyk * R011;
    term001 += uzk * R002;
    E001 = term001;
    double term200 = uxk * R300;
    term200 += uyk * R210;
    term200 += uzk * R201;
    E200 = term200;
    double term020 = uxk * R120;
    term020 += uyk * R030;
    term020 += uzk * R021;
    E020 = term020;
    double term002 = uxk * R102;
    term002 += uyk * R012;
    term002 += uzk * R003;
    E002 = term002;
    double term110 = uxk * R210;
    term110 += uyk * R120;
    term110 += uzk * R111;
    E110 = term110;
    double term101 = uxk * R201;
    term101 += uyk * R111;
    term101 += uzk * R102;
    E101 = term101;
    double term011 = uxk * R111;
    term011 += uyk * R021;
    term011 += uzk * R012;
    E011 = term011;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedIFieldCR() {
    double term000 = -sxi * R100;
    term000 -= syi * R010;
    term000 -= szi * R001;
    E000 = term000;
    double term100 = -sxi * R200;
    term100 -= syi * R110;
    term100 -= szi * R101;
    E100 = term100;
    double term010 = -sxi * R110;
    term010 -= syi * R020;
    term010 -= szi * R011;
    E010 = term010;
    double term001 = -sxi * R101;
    term001 -= syi * R011;
    term001 -= szi * R002;
    E001 = term001;
    double term200 = -sxi * R300;
    term200 -= syi * R210;
    term200 -= szi * R201;
    E200 = term200;
    double term020 = -sxi * R120;
    term020 -= syi * R030;
    term020 -= szi * R021;
    E020 = term020;
    double term002 = -sxi * R102;
    term002 -= syi * R012;
    term002 -= szi * R003;
    E002 = term002;
    double term110 = -sxi * R210;
    term110 -= syi * R120;
    term110 -= szi * R111;
    E110 = term110;
    double term101 = -sxi * R201;
    term101 -= syi * R111;
    term101 -= szi * R102;
    E101 = term101;
    double term011 = -sxi * R111;
    term011 -= syi * R021;
    term011 -= szi * R012;
    E011 = term011;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedKFieldCR() {
    double term000 = sxk * R100;
    term000 += syk * R010;
    term000 += szk * R001;
    E000 = term000;
    double term100 = sxk * R200;
    term100 += syk * R110;
    term100 += szk * R101;
    E100 = term100;
    double term010 = sxk * R110;
    term010 += syk * R020;
    term010 += szk * R011;
    E010 = term010;
    double term001 = sxk * R101;
    term001 += syk * R011;
    term001 += szk * R002;
    E001 = term001;
    double term200 = sxk * R300;
    term200 += syk * R210;
    term200 += szk * R201;
    E200 = term200;
    double term020 = sxk * R120;
    term020 += syk * R030;
    term020 += szk * R021;
    E020 = term020;
    double term002 = sxk * R102;
    term002 += syk * R012;
    term002 += szk * R003;
    E002 = term002;
    double term110 = sxk * R210;
    term110 += syk * R120;
    term110 += szk * R111;
    E110 = term110;
    double term101 = sxk * R201;
    term101 += syk * R111;
    term101 += szk * R102;
    E101 = term101;
    double term011 = sxk * R111;
    term011 += syk * R021;
    term011 += szk * R012;
    E011 = term011;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedIdX() {
    double term100 = 0.0;
    term100 -= sxi * R200;
    term100 -= syi * R110;
    term100 -= szi * R101;
    E000 = term100;
    double term200 = 0.0;
    term200 -= sxi * R300;
    term200 -= syi * R210;
    term200 -= szi * R201;
    E100 = term200;
    double term110 = 0.0;
    term110 -= sxi * R210;
    term110 -= syi * R120;
    term110 -= szi * R111;
    E010 = term110;
    double term101 = 0.0;
    term101 -= sxi * R201;
    term101 -= syi * R111;
    term101 -= szi * R102;
    E001 = term101;
    double term300 = 0.0;
    term300 -= sxi * R400;
    term300 -= syi * R310;
    term300 -= szi * R301;
    E200 = term300;
    double term120 = 0.0;
    term120 -= sxi * R220;
    term120 -= syi * R130;
    term120 -= szi * R121;
    E020 = term120;
    double term102 = 0.0;
    term102 -= sxi * R202;
    term102 -= syi * R112;
    term102 -= szi * R103;
    E002 = term102;
    double term210 = 0.0;
    term210 -= sxi * R310;
    term210 -= syi * R220;
    term210 -= szi * R211;
    E110 = term210;
    double term201 = 0.0;
    term201 -= sxi * R301;
    term201 -= syi * R211;
    term201 -= szi * R202;
    E101 = term201;
    double term111 = 0.0;
    term111 -= sxi * R211;
    term111 -= syi * R121;
    term111 -= szi * R112;
    E011 = term111;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedIdY() {
    double term010 = 0.0;
    term010 -= sxi * R110;
    term010 -= syi * R020;
    term010 -= szi * R011;
    E000 = term010;
    double term110 = 0.0;
    term110 -= sxi * R210;
    term110 -= syi * R120;
    term110 -= szi * R111;
    E100 = term110;
    double term020 = 0.0;
    term020 -= sxi * R120;
    term020 -= syi * R030;
    term020 -= szi * R021;
    E010 = term020;
    double term011 = 0.0;
    term011 -= sxi * R111;
    term011 -= syi * R021;
    term011 -= szi * R012;
    E001 = term011;
    double term210 = 0.0;
    term210 -= sxi * R310;
    term210 -= syi * R220;
    term210 -= szi * R211;
    E200 = term210;
    double term030 = 0.0;
    term030 -= sxi * R130;
    term030 -= syi * R040;
    term030 -= szi * R031;
    E020 = term030;
    double term012 = 0.0;
    term012 -= sxi * R112;
    term012 -= syi * R022;
    term012 -= szi * R013;
    E002 = term012;
    double term120 = 0.0;
    term120 -= sxi * R220;
    term120 -= syi * R130;
    term120 -= szi * R121;
    E110 = term120;
    double term111 = 0.0;
    term111 -= sxi * R211;
    term111 -= syi * R121;
    term111 -= szi * R112;
    E101 = term111;
    double term021 = 0.0;
    term021 -= sxi * R121;
    term021 -= syi * R031;
    term021 -= szi * R022;
    E011 = term021;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedIdZ() {
    double term001 = 0.0;
    term001 -= sxi * R101;
    term001 -= syi * R011;
    term001 -= szi * R002;
    E000 = term001;
    double term101 = 0.0;
    term101 -= sxi * R201;
    term101 -= syi * R111;
    term101 -= szi * R102;
    E100 = term101;
    double term011 = 0.0;
    term011 -= sxi * R111;
    term011 -= syi * R021;
    term011 -= szi * R012;
    E010 = term011;
    double term002 = 0.0;
    term002 -= sxi * R102;
    term002 -= syi * R012;
    term002 -= szi * R003;
    E001 = term002;
    double term201 = 0.0;
    term201 -= sxi * R301;
    term201 -= syi * R211;
    term201 -= szi * R202;
    E200 = term201;
    double term021 = 0.0;
    term021 -= sxi * R121;
    term021 -= syi * R031;
    term021 -= szi * R022;
    E020 = term021;
    double term003 = 0.0;
    term003 -= sxi * R103;
    term003 -= syi * R013;
    term003 -= szi * R004;
    E002 = term003;
    double term111 = 0.0;
    term111 -= sxi * R211;
    term111 -= syi * R121;
    term111 -= szi * R112;
    E110 = term111;
    double term102 = 0.0;
    term102 -= sxi * R202;
    term102 -= syi * R112;
    term102 -= szi * R103;
    E101 = term102;
    double term012 = 0.0;
    term012 -= sxi * R112;
    term012 -= syi * R022;
    term012 -= szi * R013;
    E011 = term012;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedKdX() {
    double term100 = 0.0;
    term100 += sxk * R200;
    term100 += syk * R110;
    term100 += szk * R101;
    E000 = term100;
    double term200 = 0.0;
    term200 += sxk * R300;
    term200 += syk * R210;
    term200 += szk * R201;
    E100 = term200;
    double term110 = 0.0;
    term110 += sxk * R210;
    term110 += syk * R120;
    term110 += szk * R111;
    E010 = term110;
    double term101 = 0.0;
    term101 += sxk * R201;
    term101 += syk * R111;
    term101 += szk * R102;
    E001 = term101;
    double term300 = 0.0;
    term300 += sxk * R400;
    term300 += syk * R310;
    term300 += szk * R301;
    E200 = term300;
    double term120 = 0.0;
    term120 += sxk * R220;
    term120 += syk * R130;
    term120 += szk * R121;
    E020 = term120;
    double term102 = 0.0;
    term102 += sxk * R202;
    term102 += syk * R112;
    term102 += szk * R103;
    E002 = term102;
    double term210 = 0.0;
    term210 += sxk * R310;
    term210 += syk * R220;
    term210 += szk * R211;
    E110 = term210;
    double term201 = 0.0;
    term201 += sxk * R301;
    term201 += syk * R211;
    term201 += szk * R202;
    E101 = term201;
    double term111 = 0.0;
    term111 += sxk * R211;
    term111 += syk * R121;
    term111 += szk * R112;
    E011 = term111;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedKdY() {
    double term010 = 0.0;
    term010 += sxk * R110;
    term010 += syk * R020;
    term010 += szk * R011;
    E000 = term010;
    double term110 = 0.0;
    term110 += sxk * R210;
    term110 += syk * R120;
    term110 += szk * R111;
    E100 = term110;
    double term020 = 0.0;
    term020 += sxk * R120;
    term020 += syk * R030;
    term020 += szk * R021;
    E010 = term020;
    double term011 = 0.0;
    term011 += sxk * R111;
    term011 += syk * R021;
    term011 += szk * R012;
    E001 = term011;
    double term210 = 0.0;
    term210 += sxk * R310;
    term210 += syk * R220;
    term210 += szk * R211;
    E200 = term210;
    double term030 = 0.0;
    term030 += sxk * R130;
    term030 += syk * R040;
    term030 += szk * R031;
    E020 = term030;
    double term012 = 0.0;
    term012 += sxk * R112;
    term012 += syk * R022;
    term012 += szk * R013;
    E002 = term012;
    double term120 = 0.0;
    term120 += sxk * R220;
    term120 += syk * R130;
    term120 += szk * R121;
    E110 = term120;
    double term111 = 0.0;
    term111 += sxk * R211;
    term111 += syk * R121;
    term111 += szk * R112;
    E101 = term111;
    double term021 = 0.0;
    term021 += sxk * R121;
    term021 += syk * R031;
    term021 += szk * R022;
    E011 = term021;
  }

  /** {@inheritDoc} */
  @Override
  protected void inducedKdZ() {
    double term001 = 0.0;
    term001 += sxk * R101;
    term001 += syk * R011;
    term001 += szk * R002;
    E000 = term001;
    double term101 = 0.0;
    term101 += sxk * R201;
    term101 += syk * R111;
    term101 += szk * R102;
    E100 = term101;
    double term011 = 0.0;
    term011 += sxk * R111;
    term011 += syk * R021;
    term011 += szk * R012;
    E010 = term011;
    double term002 = 0.0;
    term002 += sxk * R102;
    term002 += syk * R012;
    term002 += szk * R003;
    E001 = term002;
    double term201 = 0.0;
    term201 += sxk * R301;
    term201 += syk * R211;
    term201 += szk * R202;
    E200 = term201;
    double term021 = 0.0;
    term021 += sxk * R121;
    term021 += syk * R031;
    term021 += szk * R022;
    E020 = term021;
    double term003 = 0.0;
    term003 += sxk * R103;
    term003 += syk * R013;
    term003 += szk * R004;
    E002 = term003;
    double term111 = 0.0;
    term111 += sxk * R211;
    term111 += syk * R121;
    term111 += szk * R112;
    E110 = term111;
    double term102 = 0.0;
    term102 += sxk * R202;
    term102 += syk * R112;
    term102 += szk * R103;
    E101 = term102;
    double term012 = 0.0;
    term012 += sxk * R112;
    term012 += syk * R022;
    term012 += szk * R013;
    E011 = term012;
  }

  /** {@inheritDoc} */
  @Override
  protected void setMultipoleI(double[] Qi) {
    qi = Qi[0];
    dxi = Qi[1];
    dyi = Qi[2];
    dzi = Qi[3];
    qxxi = Qi[4] * oneThird;
    qyyi = Qi[5] * oneThird;
    qzzi = Qi[6] * oneThird;
    qxyi = Qi[7] * twoThirds;
    qxzi = Qi[8] * twoThirds;
    qyzi = Qi[9] * twoThirds;
  }

  /** {@inheritDoc} */
  @Override
  protected void setMultipoleK(double[] Qk) {
    qk = Qk[0];
    dxk = Qk[1];
    dyk = Qk[2];
    dzk = Qk[3];
    qxxk = Qk[4] * oneThird;
    qyyk = Qk[5] * oneThird;
    qzzk = Qk[6] * oneThird;
    qxyk = Qk[7] * twoThirds;
    qxzk = Qk[8] * twoThirds;
    qyzk = Qk[9] * twoThirds;
  }

  /** {@inheritDoc} */
  @Override
  protected void setDipoleI(double[] ui, double[] uiCR) {
    uxi = ui[0];
    uyi = ui[1];
    uzi = ui[2];
    pxi = uiCR[0];
    pyi = uiCR[1];
    pzi = uiCR[2];
  }

  /** {@inheritDoc} */
  @Override
  protected void setDipoleK(double[] uk, double[] ukCR) {
    uxk = uk[0];
    uyk = uk[1];
    uzk = uk[2];
    pxk = ukCR[0];
    pyk = ukCR[1];
    pzk = ukCR[2];
  }
}
