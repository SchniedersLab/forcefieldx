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

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.ScalarMath.doubleFactorial;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The abstract MultipoleTensorSIMD is extended by classes that compute derivatives of 1/|<b>r</b>| via
 * recursion to arbitrary order using Cartesian multipoles in either a global frame or a
 * quasi-internal frame.
 * <br>
 * This class serves as the abstract parent to both and defines all shared logic. Non-abstract methods
 * are declared final to disallow unnecessary overrides.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class MultipoleTensorSIMD {

  /**
   * Logger for the MultipoleTensor class.
   */
  private static final Logger logger = Logger.getLogger(MultipoleTensorSIMD.class.getName());

  /**
   * Order of the tensor recursion (5th is needed for AMOEBA forces).
   */
  protected final int order;

  /**
   * Order plus 1.
   */
  protected final int o1;

  /**
   * Order plus one.
   */
  protected final int il;

  /**
   * im = (Order plus one)^2.
   */
  protected final int im;

  /**
   * in = (Order plus one)^3.
   */
  protected final int in;

  /**
   * Size = (order + 1) * (order + 2) * (order + 3) / 6;
   */
  protected final int size;

  /**
   * The OPERATOR in use.
   */
  protected Operator operator;

  /**
   * These are the "source" terms for the recursion for the Coulomb operator (1/R).
   */
  protected final double[] coulombSource;

  /**
   * The coordinate system in use (global or QI).
   */
  protected final CoordinateSystem coordinates;

  /**
   * Separation distance.
   */
  protected DoubleVector R;

  /**
   * Separation distance squared.
   */
  protected DoubleVector r2;

  /**
   * Xk - Xi.
   */
  protected DoubleVector x;

  /**
   * Yk - Yi.
   */
  protected DoubleVector y;

  /**
   * Zk - Zi.
   */
  protected DoubleVector z;

  /**
   * A work array with auxiliary terms for the recursion.
   */
  protected final DoubleVector[] work;

  // Cartesian tensor elements (for 1/R, erfc(Beta*R)/R or DoubleVectorhole damping.
  DoubleVector R000;
  // l + m + n = 1 (3)   4
  DoubleVector R100;
  DoubleVector R010;
  DoubleVector R001;
  // l + m + n = 2 (6)  10
  DoubleVector R200;
  DoubleVector R020;
  DoubleVector R002;
  DoubleVector R110;
  DoubleVector R101;
  DoubleVector R011;
  // l + m + n = 3 (10) 20
  DoubleVector R300;
  DoubleVector R030;
  DoubleVector R003;
  DoubleVector R210;
  DoubleVector R201;
  DoubleVector R120;
  DoubleVector R021;
  DoubleVector R102;
  DoubleVector R012;
  DoubleVector R111;
  // l + m + n = 4 (15) 35
  DoubleVector R400;
  DoubleVector R040;
  DoubleVector R004;
  DoubleVector R310;
  DoubleVector R301;
  DoubleVector R130;
  DoubleVector R031;
  DoubleVector R103;
  DoubleVector R013;
  DoubleVector R220;
  DoubleVector R202;
  DoubleVector R022;
  DoubleVector R211;
  DoubleVector R121;
  DoubleVector R112;
  // l + m + n = 5 (21) 56
  DoubleVector R500;
  DoubleVector R050;
  DoubleVector R005;
  DoubleVector R410;
  DoubleVector R401;
  DoubleVector R140;
  DoubleVector R041;
  DoubleVector R104;
  DoubleVector R014;
  DoubleVector R320;
  DoubleVector R302;
  DoubleVector R230;
  DoubleVector R032;
  DoubleVector R203;
  DoubleVector R023;
  DoubleVector R311;
  DoubleVector R131;
  DoubleVector R113;
  DoubleVector R221;
  DoubleVector R212;
  DoubleVector R122;
  // l + m + n = 6 (28) 84
  DoubleVector R006;
  DoubleVector R402;
  DoubleVector R042;
  DoubleVector R204;
  DoubleVector R024;
  DoubleVector R222;
  DoubleVector R600;
  DoubleVector R060;
  DoubleVector R510;
  DoubleVector R501;
  DoubleVector R150;
  DoubleVector R051;
  DoubleVector R105;
  DoubleVector R015;
  DoubleVector R420;
  DoubleVector R240;
  DoubleVector R411;
  DoubleVector R141;
  DoubleVector R114;
  DoubleVector R330;
  DoubleVector R303;
  DoubleVector R033;
  DoubleVector R321;
  DoubleVector R231;
  DoubleVector R213;
  DoubleVector R312;
  DoubleVector R132;
  DoubleVector R123;

  // Components of the potential, field and field gradient.
  DoubleVector E000; // Potential
  // l + m + n = 1 (3)   4
  DoubleVector E100; // d/dX
  DoubleVector E010; // d/dY
  DoubleVector E001; // d/dz
  // l + m + n = 2 (6)  10
  DoubleVector E200; // d^2/dXdX
  DoubleVector E020; // d^2/dYdY
  DoubleVector E002; // d^2/dZdZ
  DoubleVector E110; // d^2/dXdY
  DoubleVector E101; // d^2/dXdZ
  DoubleVector E011; // d^2/dYdZ
  // l + m + n = 3 (10) 20
  DoubleVector E300; // d^3/dXdXdX
  DoubleVector E030; // d^3/dYdYdY
  DoubleVector E003; // d^3/dZdZdZ
  DoubleVector E210; // d^3/dXdXdY
  DoubleVector E201; // d^3/dXdXdZ
  DoubleVector E120; // d^3/dXdYdY
  DoubleVector E021; // d^3/dYdYdZ
  DoubleVector E102; // d^3/dXdZdZ
  DoubleVector E012; // d^3/dYdZdZ
  DoubleVector E111; // d^3/dXdYdZ

  /**
   * Constructor for MultipoleTensor.
   *
   * @param order       The order of the tensor.
   * @param coordinates a {@link CoordinateSystem} object.
   */
  public MultipoleTensorSIMD(CoordinateSystem coordinates, int order) {
    assert (order > 0);
    o1 = order + 1;
    il = o1;
    im = il * o1;
    in = im * o1;
    size = (order + 1) * (order + 2) * (order + 3) / 6;

    this.order = order;
    this.coordinates = coordinates;
    this.operator = Operator.COULOMB;

    // Auxiliary terms for Coulomb and Thole Screening.
    coulombSource = new double[o1];
    for (short n = 0; n <= order; n++) {
      /*
       Math.pow(-1.0, j) returns positive for all j, with -1.0 as the //
       argument rather than -1. This is a bug?
       Challacombe Eq. 21, first two factors.
      */
      coulombSource[n] = pow(-1, n) * doubleFactorial(2 * n - 1);
    }

    work = new DoubleVector[in * o1];

    // l + m + n = 0 (1)
    t000 = MultipoleUtilities.ti(0, 0, 0, order);
    // l + m + n = 1 (3)   4
    t100 = MultipoleUtilities.ti(1, 0, 0, order);
    t010 = MultipoleUtilities.ti(0, 1, 0, order);
    t001 = MultipoleUtilities.ti(0, 0, 1, order);
    // l + m + n = 2 (6)  10
    t200 = MultipoleUtilities.ti(2, 0, 0, order);
    t020 = MultipoleUtilities.ti(0, 2, 0, order);
    t002 = MultipoleUtilities.ti(0, 0, 2, order);
    t110 = MultipoleUtilities.ti(1, 1, 0, order);
    t101 = MultipoleUtilities.ti(1, 0, 1, order);
    t011 = MultipoleUtilities.ti(0, 1, 1, order);
    // l + m + n = 3 (10) 20
    t300 = MultipoleUtilities.ti(3, 0, 0, order);
    t030 = MultipoleUtilities.ti(0, 3, 0, order);
    t003 = MultipoleUtilities.ti(0, 0, 3, order);
    t210 = MultipoleUtilities.ti(2, 1, 0, order);
    t201 = MultipoleUtilities.ti(2, 0, 1, order);
    t120 = MultipoleUtilities.ti(1, 2, 0, order);
    t021 = MultipoleUtilities.ti(0, 2, 1, order);
    t102 = MultipoleUtilities.ti(1, 0, 2, order);
    t012 = MultipoleUtilities.ti(0, 1, 2, order);
    t111 = MultipoleUtilities.ti(1, 1, 1, order);
    // l + m + n = 4 (15) 35
    t400 = MultipoleUtilities.ti(4, 0, 0, order);
    t040 = MultipoleUtilities.ti(0, 4, 0, order);
    t004 = MultipoleUtilities.ti(0, 0, 4, order);
    t310 = MultipoleUtilities.ti(3, 1, 0, order);
    t301 = MultipoleUtilities.ti(3, 0, 1, order);
    t130 = MultipoleUtilities.ti(1, 3, 0, order);
    t031 = MultipoleUtilities.ti(0, 3, 1, order);
    t103 = MultipoleUtilities.ti(1, 0, 3, order);
    t013 = MultipoleUtilities.ti(0, 1, 3, order);
    t220 = MultipoleUtilities.ti(2, 2, 0, order);
    t202 = MultipoleUtilities.ti(2, 0, 2, order);
    t022 = MultipoleUtilities.ti(0, 2, 2, order);
    t211 = MultipoleUtilities.ti(2, 1, 1, order);
    t121 = MultipoleUtilities.ti(1, 2, 1, order);
    t112 = MultipoleUtilities.ti(1, 1, 2, order);
    // l + m + n = 5 (21) 56
    t500 = MultipoleUtilities.ti(5, 0, 0, order);
    t050 = MultipoleUtilities.ti(0, 5, 0, order);
    t005 = MultipoleUtilities.ti(0, 0, 5, order);
    t410 = MultipoleUtilities.ti(4, 1, 0, order);
    t401 = MultipoleUtilities.ti(4, 0, 1, order);
    t140 = MultipoleUtilities.ti(1, 4, 0, order);
    t041 = MultipoleUtilities.ti(0, 4, 1, order);
    t104 = MultipoleUtilities.ti(1, 0, 4, order);
    t014 = MultipoleUtilities.ti(0, 1, 4, order);
    t320 = MultipoleUtilities.ti(3, 2, 0, order);
    t302 = MultipoleUtilities.ti(3, 0, 2, order);
    t230 = MultipoleUtilities.ti(2, 3, 0, order);
    t032 = MultipoleUtilities.ti(0, 3, 2, order);
    t203 = MultipoleUtilities.ti(2, 0, 3, order);
    t023 = MultipoleUtilities.ti(0, 2, 3, order);
    t311 = MultipoleUtilities.ti(3, 1, 1, order);
    t131 = MultipoleUtilities.ti(1, 3, 1, order);
    t113 = MultipoleUtilities.ti(1, 1, 3, order);
    t221 = MultipoleUtilities.ti(2, 2, 1, order);
    t212 = MultipoleUtilities.ti(2, 1, 2, order);
    t122 = MultipoleUtilities.ti(1, 2, 2, order);
    // l + m + n = 6 (28) 84
    t600 = MultipoleUtilities.ti(6, 0, 0, order);
    t060 = MultipoleUtilities.ti(0, 6, 0, order);
    t006 = MultipoleUtilities.ti(0, 0, 6, order);
    t510 = MultipoleUtilities.ti(5, 1, 0, order);
    t501 = MultipoleUtilities.ti(5, 0, 1, order);
    t150 = MultipoleUtilities.ti(1, 5, 0, order);
    t051 = MultipoleUtilities.ti(0, 5, 1, order);
    t105 = MultipoleUtilities.ti(1, 0, 5, order);
    t015 = MultipoleUtilities.ti(0, 1, 5, order);
    t420 = MultipoleUtilities.ti(4, 2, 0, order);
    t402 = MultipoleUtilities.ti(4, 0, 2, order);
    t240 = MultipoleUtilities.ti(2, 4, 0, order);
    t042 = MultipoleUtilities.ti(0, 4, 2, order);
    t204 = MultipoleUtilities.ti(2, 0, 4, order);
    t024 = MultipoleUtilities.ti(0, 2, 4, order);
    t411 = MultipoleUtilities.ti(4, 1, 1, order);
    t141 = MultipoleUtilities.ti(1, 4, 1, order);
    t114 = MultipoleUtilities.ti(1, 1, 4, order);
    t330 = MultipoleUtilities.ti(3, 3, 0, order);
    t303 = MultipoleUtilities.ti(3, 0, 3, order);
    t033 = MultipoleUtilities.ti(0, 3, 3, order);
    t321 = MultipoleUtilities.ti(3, 2, 1, order);
    t231 = MultipoleUtilities.ti(2, 3, 1, order);
    t213 = MultipoleUtilities.ti(2, 1, 3, order);
    t312 = MultipoleUtilities.ti(3, 1, 2, order);
    t132 = MultipoleUtilities.ti(1, 3, 2, order);
    t123 = MultipoleUtilities.ti(1, 2, 3, order);
    t222 = MultipoleUtilities.ti(2, 2, 2, order);
  }

  /**
   * Set the separation vector.
   *
   * @param r The separation vector.
   */
  public final void setR(DoubleVector[] r) {
    setR(r[0], r[1], r[2]);
  }

  /**
   * Set the separation vector.
   *
   * @param dx Separation along the X-axis.
   * @param dy Separation along the Y-axis.
   * @param dz Separation along the Z-axis.
   */
  public abstract void setR(DoubleVector dx, DoubleVector dy, DoubleVector dz);

  /**
   * Generate the tensor using hard-coded methods.
   */
  public void generateTensor() {
    switch (order) {
      case 1 -> order1();
      case 2 -> order2();
      case 3 -> order3();
      case 4 -> order4();
      case 5 -> order5();
      case 6 -> order6();
      default -> {
        if (logger.isLoggable(Level.WARNING)) {
          logger.severe("Order " + order + " not supported.");
        }
      }
    }
  }

  /**
   * Generate source terms for the Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  protected abstract void source(DoubleVector[] T000);

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 1st order.
   */
  protected abstract void order1();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 2nd order.
   */
  protected abstract void order2();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 3rd order.
   */
  protected abstract void order3();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 4th order.
   */
  protected abstract void order4();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 5th order, which is needed for
   * quadrupole-quadrupole forces.
   */
  protected abstract void order5();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 6th order, which is needed for
   * quadrupole-quadrupole forces and orthogonal space sampling.
   */
  protected abstract void order6();

  /**
   * The index is based on the idea of filling tetrahedron.
   * <p>
   * 1/r has an index of 0.
   * <br>
   * derivatives of x are first; indices from 1..o for d/dx..(d/dx)^o
   * <br>
   * derivatives of x and y are second; base triangle of size (o+1)(o+2)/2
   * <br>
   * derivatives of x, y and z are last; total size (o+1)*(o+2)*(o+3)/6
   * <br>
   * <p>
   * This function is useful to set up masking constants:
   * <br>
   * static int Tlmn = ti(l,m,n,order)
   * <br>
   * For example the (d/dy)^2 (1/R) storage location:
   * <br>
   * static int T020 = ti(0,2,0,order)
   *
   * @param dx int The number of d/dx operations.
   * @param dy int The number of d/dy operations.
   * @param dz int The number of d/dz operations.
   * @return int in the range (0..binomial(order + 3, 3) - 1)
   */
  protected final int ti(int dx, int dy, int dz) {
    return MultipoleUtilities.ti(dx, dy, dz, order);
  }

  // l + m + n = 0 (1)
  /**
   * No derivatives.
   */
  protected final int t000;
  // l + m + n = 1 (3)   4
  /**
   * First derivative with respect to x.
   */
  protected final int t100;
  /**
   * First derivative with respect to y.
   */
  protected final int t010;
  /**
   * First derivative with respect to z.
   */
  protected final int t001;
  // l + m + n = 2 (6)  10
  /**
   * Second derivative with respect to x.
   */
  protected final int t200;
  /**
   * Second derivative with respect to y.
   */
  protected final int t020;
  /**
   * Second derivative with respect to z.
   */
  protected final int t002;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t110;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t101;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t011;
  // l + m + n = 3 (10) 20
  /**
   * Third derivative with respect to x.
   */
  protected final int t300;
  /**
   * Third derivative with respect to y.
   */
  protected final int t030;
  /**
   * Third derivative with respect to z.
   */
  protected final int t003;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t210;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t201;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t120;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t021;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t102;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t012;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t111;
  // l + m + n = 4 (15) 35
  /**
   * Fourth derivative with respect to x.
   */
  protected final int t400;
  /**
   * Fourth derivative with respect to y.
   */
  protected final int t040;
  /**
   * Fourth derivative with respect to z.
   */
  protected final int t004;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t310;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t301;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t130;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t031;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t103;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t013;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t220;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t202;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t022;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t211;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t121;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t112;

  // l + m + n = 5 (21) 56
  /**
   * Fifth derivative with respect to x.
   */
  protected final int t500;
  /**
   * Fifth derivative with respect to y.
   */
  protected final int t050;
  /**
   * Fifth derivative with respect to z.
   */
  protected final int t005;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t410;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t401;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t140;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t041;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t104;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t014;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t320;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t302;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t230;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t032;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t203;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t023;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t311;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t131;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t113;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t221;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t212;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t122;

  // l + m + n = 6 (28) 84
  /**
   * Sixth derivative with respect to x.
   */
  protected final int t600;
  /**
   * Sixth derivative with respect to y.
   */
  protected final int t060;
  /**
   * Sixth derivative with respect to z.
   */
  protected final int t006;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t510;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t501;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t150;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t051;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t105;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t015;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t420;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t402;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t240;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t042;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t204;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t024;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t411;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t141;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t114;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t330;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t303;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t033;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t321;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t231;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t213;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t312;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t132;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t123;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t222;

  @SuppressWarnings("fallthrough")
  protected abstract void multipoleIPotentialAtK(PolarizableMultipoleSIMD mI, int order);

  @SuppressWarnings("fallthrough")
  protected abstract void multipoleKPotentialAtI(PolarizableMultipoleSIMD mK, int order);
}
