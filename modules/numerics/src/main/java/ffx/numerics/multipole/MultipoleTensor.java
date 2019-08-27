//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.numerics.multipole;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Math.PI;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import static ffx.numerics.math.VectorMath.binomial;
import static ffx.numerics.math.VectorMath.doubleFactorial;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.special.Erf.erfc;

/**
 * The abstract MultipoleTensor is extended by classes that compute derivatives of 1/|<b>r</b>| via recursion
 * to arbitrary order using Cartesian multipoles in either a global frame or a
 * quasi-internal frame.
 * <br>
 * This class serves as the abstract parent to both and defines all shared logic.
 * Non-abstract methods are declared final to disallow unnecessary overrides.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class MultipoleTensor {

    private static final Logger logger = Logger.getLogger(MultipoleTensor.class.getName());

    /**
     * Operators that are supported.
     */
    public enum OPERATOR {
        COULOMB, SCREENED_COULOMB, THOLE_FIELD
    }

    /**
     * Global and Quasi-Internal (QI) coordinate systems are supported.
     */
    public enum COORDINATES {
        GLOBAL, QI
    }

    protected OPERATOR operator;
    protected COORDINATES coordinates;

    protected final int order;
    /**
     * These are the "source" terms for the recursion.
     */
    private double[] T000j;
    /**
     * These are the "source" terms for the recursion for the Coulomb operator (1/R).
     */
    protected final double[] coulomb;
    /**
     * These are the "source" terms for the recursion for the screened Coulomb operator erfc(R)/R.
     */
    private final double[] screened;
    /**
     * Ewald parameter.
     */
    protected final double beta;
    /**
     * Thole damping parameter (= min(pti,ptk)).
     */
    protected double pgamma;
    /**
     * 1/(alphai^6*alphak^6) where alpha is polarizability.
     */
    private double aiak;
    /**
     * Separation distance.
     */
    protected double R;
    /**
     * Separation distance squared.
     */
    protected double r2;
    /**
     * Xk - Xi.
     */
    protected double x;
    /**
     * Yk - Yi.
     */
    protected double y;
    /**
     * Zk - Zi.
     */
    protected double z;
    protected final int o1;
    protected final int il;
    protected final int im;
    protected final int in;
    protected final int size;
    /**
     * Store the auxillary tensor memory to avoid memory consumption.
     */
    final double[] T000;
    /**
     * Store the work array to avoid memory consumption. Note that rather than
     * use an array for intermediate values, a 4D matrix was tried. It was
     * approximately 50% slower than the linear work array.
     */
    protected final double[] work;
    /**
     * Stores previous distance (and lambdaFunction), for recycling check.
     */
    final double[] rprev = new double[4];

    /**
     * Constant <code>recycleTensors=false</code>
     */
    public static final boolean recycleTensors = false;
    private int tensorsRecycled = 0;

    /**
     * <p>
     * Constructor for MultipoleTensor.</p>
     *
     * @param operator    The tensor operator.
     * @param order       The order of the tensor.
     * @param aewald      The screening parameter.
     * @param coordinates a {@link MultipoleTensor.COORDINATES} object.
     */
    public MultipoleTensor(OPERATOR operator, COORDINATES coordinates, int order, double aewald) {
        assert (order > 0);
        o1 = order + 1;
        il = o1;
        im = il * o1;
        in = im * o1;
        size = (order + 1) * (order + 2) * (order + 3) / 6;
        work = new double[in * o1];

        this.order = order;
        this.operator = operator;
        this.coordinates = coordinates;
        this.beta = aewald;
        if (operator == OPERATOR.SCREENED_COULOMB && aewald == 0.0) {
            // logger.warning("Tried beta of zero for screened coulomb tensor.");
            // Switch to the Coulomb operator.
            operator = OPERATOR.COULOMB;
        }

        // Auxillary terms for Coulomb and Thole Screening.
        coulomb = new double[o1];
        for (int n = 0; n <= order; n++) {

            /*
              Math.pow(-1.0, j) returns positive for all j, with -1.0 as the //
              argument rather than -1. This is a bug?

              Challacombe Eq. 21, first two factors.
             */
            coulomb[n] = pow(-1, n) * doubleFactorial(2 * n - 1);
        }

        // Auxillary terms for screened Coulomb (Sagui et al. Eq. 2.28)
        screened = new double[o1];
        double prefactor = 2.0 * aewald / sqrtPI;
        double twoBeta2 = -2.0 * aewald * aewald;
        for (int n = 0; n <= order; n++) {
            screened[n] = prefactor * pow(twoBeta2, n);
        }

        setOperator(operator);

        T000 = new double[order + 1];
        // l + m + n = 0 (1)
        t000 = ti(0, 0, 0, order);
        // l + m + n = 1 (3)   4
        t100 = ti(1, 0, 0, order);
        t010 = ti(0, 1, 0, order);
        t001 = ti(0, 0, 1, order);
        // l + m + n = 2 (6)  10
        t200 = ti(2, 0, 0, order);
        t020 = ti(0, 2, 0, order);
        t002 = ti(0, 0, 2, order);
        t110 = ti(1, 1, 0, order);
        t101 = ti(1, 0, 1, order);
        t011 = ti(0, 1, 1, order);
        // l + m + n = 3 (10) 20
        t300 = ti(3, 0, 0, order);
        t030 = ti(0, 3, 0, order);
        t003 = ti(0, 0, 3, order);
        t210 = ti(2, 1, 0, order);
        t201 = ti(2, 0, 1, order);
        t120 = ti(1, 2, 0, order);
        t021 = ti(0, 2, 1, order);
        t102 = ti(1, 0, 2, order);
        t012 = ti(0, 1, 2, order);
        t111 = ti(1, 1, 1, order);
        // l + m + n = 4 (15) 35
        t400 = ti(4, 0, 0, order);
        t040 = ti(0, 4, 0, order);
        t004 = ti(0, 0, 4, order);
        t310 = ti(3, 1, 0, order);
        t301 = ti(3, 0, 1, order);
        t130 = ti(1, 3, 0, order);
        t031 = ti(0, 3, 1, order);
        t103 = ti(1, 0, 3, order);
        t013 = ti(0, 1, 3, order);
        t220 = ti(2, 2, 0, order);
        t202 = ti(2, 0, 2, order);
        t022 = ti(0, 2, 2, order);
        t211 = ti(2, 1, 1, order);
        t121 = ti(1, 2, 1, order);
        t112 = ti(1, 1, 2, order);
        // l + m + n = 5 (21) 56
        t500 = ti(5, 0, 0, order);
        t050 = ti(0, 5, 0, order);
        t005 = ti(0, 0, 5, order);
        t410 = ti(4, 1, 0, order);
        t401 = ti(4, 0, 1, order);
        t140 = ti(1, 4, 0, order);
        t041 = ti(0, 4, 1, order);
        t104 = ti(1, 0, 4, order);
        t014 = ti(0, 1, 4, order);
        t320 = ti(3, 2, 0, order);
        t302 = ti(3, 0, 2, order);
        t230 = ti(2, 3, 0, order);
        t032 = ti(0, 3, 2, order);
        t203 = ti(2, 0, 3, order);
        t023 = ti(0, 2, 3, order);
        t311 = ti(3, 1, 1, order);
        t131 = ti(1, 3, 1, order);
        t113 = ti(1, 1, 3, order);
        t221 = ti(2, 2, 1, order);
        t212 = ti(2, 1, 2, order);
        t122 = ti(1, 2, 2, order);
        // l + m + n = 6 (28) 84
        t600 = ti(6, 0, 0, order);
        t060 = ti(0, 6, 0, order);
        t006 = ti(0, 0, 6, order);
        t510 = ti(5, 1, 0, order);
        t501 = ti(5, 0, 1, order);
        t150 = ti(1, 5, 0, order);
        t051 = ti(0, 5, 1, order);
        t105 = ti(1, 0, 5, order);
        t015 = ti(0, 1, 5, order);
        t420 = ti(4, 2, 0, order);
        t402 = ti(4, 0, 2, order);
        t240 = ti(2, 4, 0, order);
        t042 = ti(0, 4, 2, order);
        t204 = ti(2, 0, 4, order);
        t024 = ti(0, 2, 4, order);
        t411 = ti(4, 1, 1, order);
        t141 = ti(1, 4, 1, order);
        t114 = ti(1, 1, 4, order);
        t330 = ti(3, 3, 0, order);
        t303 = ti(3, 0, 3, order);
        t033 = ti(0, 3, 3, order);
        t321 = ti(3, 2, 1, order);
        t231 = ti(2, 3, 1, order);
        t213 = ti(2, 1, 3, order);
        t312 = ti(3, 1, 2, order);
        t132 = ti(1, 3, 2, order);
        t123 = ti(1, 2, 3, order);
        t222 = ti(2, 2, 2, order);
    }

    /**
     * <p>isGlobal.</p>
     *
     * @return a boolean.
     */
    public final boolean isGlobal() {
        return (MultipoleTensorGlobal.class.isAssignableFrom(this.getClass()));
    }

    /**
     * <p>isQI.</p>
     *
     * @return a boolean.
     */
    public final boolean isQI() {
        return (MultipoleTensorQI.class.isAssignableFrom(this.getClass()));
    }

    /**
     * <p>permScreened.</p>
     *
     * @param r              an array of {@link double} objects.
     * @param lambdaFunction a double.
     * @param Qi             an array of {@link double} objects.
     * @param Qk             an array of {@link double} objects.
     */
    public void permScreened(double[] r, double lambdaFunction, double[] Qi, double[] Qk) {
        boolean operatorChange = operator != OPERATOR.SCREENED_COULOMB;
        if (operatorChange) {
            setOperator(OPERATOR.SCREENED_COULOMB);
        }
        boolean distanceChanged = setR(r, lambdaFunction);
        setMultipoles(Qi, Qk);
        unsetDipoles();
        unsetDamping();
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        } else {
            tensorsRecycled++;
        }
    }

    /**
     * <p>permCoulomb.</p>
     *
     * @param r              an array of {@link double} objects.
     * @param lambdaFunction a double.
     * @param Qi             an array of {@link double} objects.
     * @param Qk             an array of {@link double} objects.
     */
    public void permCoulomb(double[] r, double lambdaFunction, double[] Qi, double[] Qk) {
        boolean operatorChange = operator != OPERATOR.COULOMB;
        if (operatorChange) {
            setOperator(OPERATOR.COULOMB);
        }
        boolean distanceChanged = setR(r, lambdaFunction);
        setMultipoles(Qi, Qk);
        unsetDipoles();
        unsetDamping();
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        } else {
            tensorsRecycled++;
        }
    }

    /**
     * <p>polarScreened.</p>
     *
     * @param r    an array of {@link double} objects.
     * @param Qi   an array of {@link double} objects.
     * @param Qk   an array of {@link double} objects.
     * @param ui   an array of {@link double} objects.
     * @param uiCR an array of {@link double} objects.
     * @param uk   an array of {@link double} objects.
     * @param ukCR an array of {@link double} objects.
     */
    public void polarScreened(double[] r, double[] Qi, double[] Qk,
                              double[] ui, double[] uiCR, double[] uk, double[] ukCR) {
        boolean operatorChange = operator != OPERATOR.SCREENED_COULOMB;
        if (operatorChange) {
            setOperator(OPERATOR.SCREENED_COULOMB);
        }
        boolean distanceChanged = setR(r);
        setMultipoles(Qi, Qk);
        setDipoles(ui, uiCR, uk, ukCR);
        unsetDamping();
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        } else {
            tensorsRecycled++;
        }
    }

    /**
     * <p>polarCoulomb.</p>
     *
     * @param r    an array of {@link double} objects.
     * @param Qi   an array of {@link double} objects.
     * @param Qk   an array of {@link double} objects.
     * @param ui   an array of {@link double} objects.
     * @param uiCR an array of {@link double} objects.
     * @param uk   an array of {@link double} objects.
     * @param ukCR an array of {@link double} objects.
     */
    public void polarCoulomb(double[] r, double[] Qi, double[] Qk,
                             double[] ui, double[] uiCR, double[] uk, double[] ukCR) {
        boolean operatorChange = operator != OPERATOR.COULOMB;
        if (operatorChange) {
            setOperator(OPERATOR.COULOMB);
        }
        boolean distanceChanged = setR(r);
        setMultipoles(Qi, Qk);
        setDipoles(ui, uiCR, uk, ukCR);
        unsetDamping();
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        } else {
            tensorsRecycled++;
        }
    }

    /**
     * For Thole tensors.
     *
     * @param r      an array of {@link double} objects.
     * @param Qi     an array of {@link double} objects.
     * @param Qk     an array of {@link double} objects.
     * @param ui     an array of {@link double} objects.
     * @param uiCR   an array of {@link double} objects.
     * @param uk     an array of {@link double} objects.
     * @param ukCR   an array of {@link double} objects.
     * @param pgamma a double.
     * @param aiak   a double.
     */
    public void tholeField(double[] r, double[] Qi, double[] Qk,
                           double[] ui, double[] uiCR, double[] uk, double[] ukCR, double pgamma, double aiak) {
        boolean operatorChange = operator != OPERATOR.THOLE_FIELD;
        if (operatorChange) {
            setOperator(OPERATOR.THOLE_FIELD);
        }

        boolean distanceChanged = setR(r);
        setMultipoles(Qi, Qk);
        setDipoles(ui, uiCR, uk, ukCR);
        setTholeDamping(pgamma, aiak);
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        } else {
            tensorsRecycled++;
        }
    }

    /**
     * For the MultipoleTensorTest class and testing.
     *
     * @param r              an array of {@link double} objects.
     * @param lambdaFunction a double.
     * @param Qi             an array of {@link double} objects.
     * @param Qk             an array of {@link double} objects.
     */
    final void generateTensor(double[] r, double lambdaFunction, double[] Qi, double[] Qk) {
        boolean distanceChanged = setR(r, lambdaFunction);
        setMultipoles(Qi, Qk);
        unsetDipoles();
        unsetDamping();
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        } else {
            tensorsRecycled++;
        }
    }

    /**
     * Generate the tensor for the interaction between Qi and Qk.
     *
     * @param r    an array of {@link double} objects.
     * @param Qi   an array of {@link double} objects.
     * @param Qk   an array of {@link double} objects.
     * @param ui   an array of {@link double} objects.
     * @param uiCR an array of {@link double} objects.
     * @param uk   an array of {@link double} objects.
     * @param ukCR an array of {@link double} objects.
     */
    public final void generateTensor(double[] r, double[] Qi, double[] Qk,
                                     double[] ui, double[] uiCR, double[] uk, double[] ukCR) {
        boolean distanceChanged = setR(r);
        setMultipoles(Qi, Qk);
        setDipoles(ui, uiCR, uk, ukCR);
        unsetDamping();
        if (!recycleTensors || distanceChanged) {
            generateTensor();
        }
    }

    /**
     * <p>getRecycledCount.</p>
     *
     * @return a int.
     */
    public final int getRecycledCount() {
        return tensorsRecycled;
    }

    /**
     * <p>resetRecycledCount.</p>
     */
    public final void resetRecycledCount() {
        tensorsRecycled = 0;
    }

    /**
     * Set the Operator.
     *
     * @param operator a {@link MultipoleTensor.OPERATOR} object.
     */
    public final void setOperator(OPERATOR operator) {
        this.operator = operator;

        OPERATOR op = operator;
        if (operator == OPERATOR.SCREENED_COULOMB && beta == 0.0) {
            // logger.warning("Tried beta of zero for screened coulomb tensor.");
            // Switch to the Coulomb operator.
            op = OPERATOR.COULOMB;
        }

        switch (op) {
            case SCREENED_COULOMB:
                T000j = screened;
                break;
            default:
            case THOLE_FIELD:
            case COULOMB:
                T000j = coulomb;
        }

    }

    /**
     * Package-private for use by MultipoleTensorTest class only.
     */
    void generateTensor() {
        switch (order) {
            case 6:
                order6();
                break;
            case 5:
                order5();
                break;
            case 4:
                order4();
                break;
            default:
                double[] r = {x, y, z};
                recursion(r, work);
        }
    }

    /**
     * Set the Thole damping parameters.
     *
     * @param pgamma a double.
     * @param aiak   1/(alphai^6*alphak^6) where alpha is polarizability
     */
    final void setTholeDamping(double pgamma, double aiak) {
        this.pgamma = pgamma;
        this.aiak = aiak;        // == 1/(alphai*alphak)^6 where alpha is polarizability
    }

    /**
     * <p>checkDampingCriterion.</p>
     *
     * @param dx_local an array of {@link double} objects.
     * @param pgamma   a double.
     * @param aiak     a double.
     * @return a boolean.
     */
    public static boolean checkDampingCriterion(double[] dx_local, double pgamma, double aiak) {
        double R = r(dx_local);
        return (-pgamma * (R * aiak) * (R * aiak) * (R * aiak) > -50.0);
    }

    /**
     * <p>unsetDipoles.</p>
     */
    private void unsetDipoles() {
        uxi = uyi = uzi = pxi = pyi = pzi = 0.0;
        uxk = uyk = uzk = pxk = pyk = pzk = 0.0;
    }

    /**
     * <p>unsetDamping.</p>
     */
    private void unsetDamping() {
        pgamma = aiak = 0.0;
    }

    /**
     * <p>setR.</p>
     *
     * @param r Whether distance changed as a result.
     * @return a boolean.
     */
    final boolean setR(double[] r) {
        return setR(r, 0.0);
    }

    /**
     * <p>setR.</p>
     *
     * @param r              Whether distance changed as a result.
     * @param lambdaFunction a double.
     * @return a boolean.
     */
    protected abstract boolean setR(double[] r, double lambdaFunction);

    /**
     * Generate source terms for the Challacombe et al. recursion.
     *
     * @param T000 Location to store the source terms.
     */
    protected final void source(double[] T000) {

        OPERATOR op = operator;
        if (operator == OPERATOR.SCREENED_COULOMB && beta == 0.0) {
            // Generate tensors for the Coulomb operator.
            op = OPERATOR.COULOMB;
        }

        switch (op) {
            case SCREENED_COULOMB:
                // Sagui et al. Eq. 2.22
                double betaR = beta * R;
                double betaR2 = betaR * betaR;
                double iBetaR2 = 1.0 / (2.0 * betaR2);
                double expBR2 = exp(-betaR2);
                // Fnc(x^2) = Sqrt(PI) * erfc(x) / (2*x)
                // where x = Beta*R
                double Fnc = sqrtPI * erfc(betaR) / (2.0 * betaR);
                for (int n = 0; n < o1; n++) {
                    T000[n] = T000j[n] * Fnc;
                    // Generate F(n+1)c from Fnc (Eq. 2.24 in Sagui et al.)
                    // F(n+1)c = [(2*n+1) Fnc(x) + exp(-x)] / 2x
                    // where x = (Beta*R)^2
                    Fnc = ((2.0 * n + 1.0) * Fnc + expBR2) * iBetaR2;
                }
                break;
            case THOLE_FIELD:
                assert (order <= 4);
                double ir = 1.0 / R;
                double ir2 = ir * ir;
                for (int n = 0; n < o1; n++) {
                    T000[n] = T000j[n] * ir;
                    ir *= ir2;
                }

                // Add the Thole damping terms: edamp = exp(-damp*u^3).
                double u = R * aiak;
                double u3 = pgamma * u * u * u;
                double u6 = u3 * u3;
                double u9 = u6 * u3;
                double expU3 = exp(-u3);

                // The zeroth order term is not calculated for Thole damping.
                T000[0] = 0.0;
                T000[1] *= expU3;
                T000[2] *= (1.0 + u3) * expU3;
                T000[3] *= (1.0 + u3 + threeFifths * u6) * expU3;
                T000[4] *= (1.0 + u3 + (18.0 * u6 + 9.0 * u9) * oneThirtyFifth) * expU3;
                break;
            case COULOMB:
            default:
                // Challacombe et al. Equation 21, last factor.
                // == (1/r) * (1/r^3) * (1/r^5) * (1/r^7) * ...
                ir = 1.0 / R;
                ir2 = ir * ir;
                for (int n = 0; n < o1; n++) {
                    T000[n] = T000j[n] * ir;
                    ir *= ir2;
                }
        }
    }

    /**
     * <p>log.</p>
     *
     * @param tensor an array of {@link double} objects.
     */
    public final void log(double[] tensor) {
        log(this.operator, this.order, tensor);
    }

    /**
     * Log the tensors.
     *
     * @param operator The OPERATOR to use.
     * @param order    The tensor order.
     * @param tensor   The tensor array.
     */
    private static void log(OPERATOR operator, int order, double[] tensor) {
        final int o1 = order + 1;
        StringBuilder sb = new StringBuilder();

        sb.append(String.format("\n %s Operator to order %d:", operator, order));
        sb.append(String.format("\n%5s %4s %4s %4s %12s\n", "Index", "d/dx", "d/dy", "d/dz", "Tensor"));
        sb.append(String.format("%5d %4d %4d %4d %12.8f\n", 0, 0, 0, 0, tensor[0]));
        int count = 1;
        // Print (d/dx)^l for l = 1..order (m = 0, n = 0)
        for (int l = 1; l <= order; l++) {
            double value = tensor[ti(l, 0, 0, order)];
            if (value != 0.0) {
                sb.append(String.format("%5d %4d %4d %4d %12.8f\n", ti(l, 0, 0, order), l, 0, 0, value));
                count++;
            }
        }
        // Print (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
        for (int l = 0; l <= o1; l++) {
            for (int m = 1; m <= order - l; m++) {
                double value = tensor[ti(l, m, 0, order)];
                if (value != 0.0) {
                    sb.append(String.format("%5d %4d %4d %4d %12.8f\n", ti(l, m, 0, order), l, m, 0, value));
                    count++;
                }
            }
        }
        // Print (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
        for (int l = 0; l <= o1; l++) {
            for (int m = 0; m <= o1 - l; m++) {
                for (int n = 1; n <= order - l - m; n++) {
                    double value = tensor[ti(l, m, n, order)];
                    if (value != 0.0) {
                        sb.append(String.format("%5d %4d %4d %4d %12.8f\n", ti(l, m, n, order), l, m, n, value));
                        count++;
                    }
                }
            }
        }
        sb.append(String.format("\n Total number of active tensors: %d\n", count));
        logger.log(Level.INFO, sb.toString());
    }

    /**
     * Returns the number of tensors for derivatives to the given order.
     *
     * @param order maximum number of derivatives.
     * @return the number of tensors.
     * @since 1.0
     */
    public static int tensorCount(int order) {
        long ret = binomial(order + 3, 3);
        assert (ret < Integer.MAX_VALUE);
        return (int) ret;
    }

    /**
     * The index is based on the idea of filling tetrahedron.
     * <p>
     * 1/r has an index of 0
     * <br>
     * derivatives of x are first; indeces from 1..o for d/dx..(d/dx)^o
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
     * For example the (d/dy)^2 (1/R) storage location: <br> static int T020 =
     * ti(0,2,0,order)
     * <p>
     *
     * @param dx    int The number of d/dx operations.
     * @param dy    int The number of d/dy operations.
     * @param dz    int The number of d/dz operations.
     * @param order int The maximum tensor order (0 .LE. dx + dy + dz .LE.
     *              order).
     * @return int in the range (0..binomial(order + 3, 3) - 1)
     */
    static int ti(int dx, int dy, int dz, int order) {
        if (dx < 0 || dy < 0 || dz < 0 || dx + dy + dz > order) {
            return -1;
        }

        int size = (order + 1) * (order + 2) * (order + 3) / 6;
        /*
          We only get to the top of the tetrahedron if dz = order, otherwise
          subtract off the top, including the level of the requested tensor
          index.
         */
        int top = order + 1 - dz;
        top = top * (top + 1) * (top + 2) / 6;
        int zindex = size - top;
        /*
          Given the "dz level", dy can range from 0..order - dz) To get to the
          row for a specific value of dy, dy*(order + 1) - dy*(dy-1)/2 indeces
          are skipped. This is an operation that looks like the area of
          rectangle, minus the area of an empty triangle.
         */
        int yindex = dy * (order - dz) - (dy - 1) * (dy - 2) / 2 + 1;
        /*
          Given the dz level and dy row, dx can range from (0..order - dz - dy)
          The dx index is just walking down the dy row for "dx" steps.
         */
        return dx + yindex + zindex;
    }

    /**
     * Multidimensional arrays into 1-D; doing so should be much more efficient
     * (less cache misses from indirection).
     *
     * @see {@code MultipoleTensor::ti(intdx, int dy, int dz)}
     */
    final int ti(int dx, int dy, int dz) {
        return ti(dx, dy, dz, order);
    }

    /**
     * <p>unrolled.</p>
     *
     * @param r      an array of {@link double} objects.
     * @param tensor an array of {@link double} objects.
     */
    final void unrolled(final double[] r, final double[] tensor) {
        switch (order) {
            case 4:
                order4();
                break;
            case 5:
                order5();
                break;
            case 6:
                order6();
                break;
            default:
                throw new IllegalArgumentException();
        }
    }

    /**
     * This function is a driver to collect elements of the Cartesian multipole
     * tensor. Collecting all tensors scales slightly better than O(order^4).
     * <p>
     * For a multipole expansion truncated at quadrupole order, for example, up
     * to order 5 is needed for energy gradients. The number of terms this
     * requires is binomial(5 + 3, 3) or 8! / (5! * 3!), which is 56.
     * <p>
     * The packing of the tensor elements for order = 1<br> tensor[0] = 1/|r|
     * <br>
     * tensor[1] = -x/|r|^3 <br> tensor[2] = -y/|r|^3 <br> tensor[3] = -z/|r|^3
     * <br>
     * <p>
     *
     * @param r      double[] vector between two sites.
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     * @return Java code for the tensor recursion.
     * @since 1.0
     */
    protected String codeTensorRecursion(final double[] r, final double[] tensor) {
        setR(r);
        assert (x == 0.0 && y == 0.0);
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
        // Store the l=1 tensor T100 (d/dx)
        // Starting the loop at l=2 avoids an if statement.
        double current;
        double previous = work[1];
        tensor[ti(1, 0, 0)] = 0.0;
        for (int l = 2; l < o1; l++) {
            // Initial condition for the inner loop is formation of T100(l-1).
            // Starting the inner loop at a=1 avoids an if statement.
            // T100(l-1) = 0.0 * T000(l)
            current = 0.0;
            int iw = il + (l - 1);
            work[iw] = current;
            // sb.append(format("double %s = 0.0;\n", term(1, 0, 0, l - 1)));
            for (int a = 1; a < l - 1; a++) {
                // T200(l-2) = 0.0 * T100(l-1) + (2 - 1) * T000(l-1)
                // T300(l-3) = 0.0 * T200(l-2) + (3 - 1) * T100(l-2)
                // ...
                // T(l-1)001 = 0.0 * T(l-2)002 + (l - 2) * T(l-3)002
                // iw = (a - 1) * il + (l - a)
                current = a * work[iw - il];
                iw += il - 1;
                // iw = (a + 1) * il + (l - a - 1)
                work[iw] = current;
                if (current != 0) {
                    if (a > 2) {
                        sb.append(format("double %s = %d * %s;\n",
                                term(a + 1, 0, 0, l - a - 1), a, term(a - 1, 0, 0, l - a)));
                    } else {
                        sb.append(format("double %s = %s;\n",
                                term(a + 1, 0, 0, l - a - 1), term(a - 1, 0, 0, l - a)));
                    }
                }
            }
            // Store the Tl00 tensor (d/dx)^l
            // Tl00 = 0.0 * T(l-1)001 + (l - 1) * T(l-2)001
            int index = ti(l, 0, 0);
            tensor[index] = (l - 1) * previous;
            previous = current;
            if (tensor[index] != 0) {
                if (l > 2) {
                    sb.append(format("%s = %d * %s;\n", rlmn(l, 0, 0), (l - 1), term(l - 2, 0, 0, 1)));
                } else {
                    sb.append(format("%s = %s;\n", rlmn(l, 0, 0), term(l - 2, 0, 0, 1)));
                }
            }
        }
        // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
        // Any (d/dy) term can be formed as:
        // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
        for (int l = 0; l < order; l++) {
            // Store the m=1 tensor (d/dx)^l *(d/dy)
            // Tl10 = y * Tl001
            previous = work[l * il + 1];
            tensor[ti(l, 1, 0)] = 0.0;
            for (int m = 2; m + l < o1; m++) {
                // Tl10(m-1) = y * Tl00m;
                int iw = l * il + m;
                current = 0.0;
                iw += im - 1;
                work[iw] = current;
                for (int a = 1; a < m - 1; a++) {
                    // Tl20(m-2) = 0.0 * Tl10(m-1) + (2 - 1) * Tl00(m-1)
                    // Tl30(m-3) = 0.0 * Tl20(m-2) + (3 - 1) * Tl10(m-2)
                    // ...
                    // Tl(m-1)01 = 0.0 * Tl(m-2)02 + (m - 2) * Tl(m-3)02
                    current = a * work[iw - im];
                    iw += im - 1;
                    work[iw] = current;
                    if (current != 0) {
                        if (a > 1) {
                            sb.append(format("double %s = %d * %s;\n",
                                    term(l, a + 1, 0, m - a - 1), a, term(l, a - 1, 0, m - a)));
                        } else {
                            sb.append(format("double %s = %s;\n",
                                    term(l, a + 1, 0, m - a - 1), term(l, a - 1, 0, m - a)));
                        }

                    }
                }
                // Store the tensor (d/dx)^l * (d/dy)^m
                // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
                int index = ti(l, m, 0);
                tensor[index] = (m - 1) * previous;
                previous = current;
                if (tensor[index] != 0) {
                    if (m > 2) {
                        sb.append(format("%s = %d * %s;\n", rlmn(l, m, 0), (m - 1), term(l, m - 2, 0, 1)));
                    } else {
                        sb.append(format("%s = %s;\n", rlmn(l, m, 0), term(l, m - 2, 0, 1)));
                    }
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
                int index = ti(l, m, 1);
                tensor[index] = z * previous;
                if (tensor[index] != 0) {
                    sb.append(format("%s = z * %s;\n", rlmn(l, m, 1), term(l, m, 0, 1)));
                }
                for (int n = 2; lm + n < o1; n++) {
                    // Tlm1(n-1) = z * Tlm0n;
                    int iw = lilmim + n;
                    current = z * work[iw];
                    iw += in - 1;
                    work[iw] = current;
                    if (current != 0) {
                        sb.append(format("double %s = z * %s;\n",
                                term(l, m, 1, n - 1), term(l, m, 0, n)));
                    }
                    final int n1 = n - 1;
                    for (int a = 1; a < n1; a++) {
                        // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
                        // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
                        // ...
                        // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
                        current = z * current + a * work[iw - in];
                        iw += in - 1;
                        work[iw] = current;
                        if (current != 0) {
                            if (a > 1) {
                                sb.append(format("double %s = z * %s + %d * %s;\n",
                                        term(l, m, a + 1, n - a - 1), term(l, m, a, n - a), a, term(l, m, a - 1, n - a)));
                            } else {
                                sb.append(format("double %s = z * %s + %s;\n",
                                        term(l, m, a + 1, n - a - 1), term(l, m, a, n - a), term(l, m, a - 1, n - a)));
                            }
                        }
                    }
                    // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
                    // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
                    index = ti(l, m, n);
                    tensor[index] = z * current + n1 * previous;
                    previous = current;
                    if (tensor[index] != 0) {
                        if (n > 2) {
                            sb.append(format("%s = z * %s + %d * %s;\n",
                                    rlmn(l, m, n), term(l, m, n - 1, 1), (n - 1), term(l, m, n - 2, 1)));
                        } else {
                            sb.append(format("%s = z * %s + %s;\n",
                                    rlmn(l, m, n), term(l, m, n - 1, 1), term(l, m, n - 2, 1)));
                        }
                    }
                }
            }
        }
        return sb.toString();
    }

    /**
     * Contract multipole moments with their respective electrostatic potential
     * derivatives.
     *
     * @param T array of electrostatic potential and partial derivatives
     * @param l apply (d/dx)^l to the potential
     * @param m apply (d/dy)^l to the potential
     * @param n apply (d/dz)^l to the potential
     * @return the contracted interaction.
     */
    private double contract(double[] T, int l, int m, int n) {
        double total = 0.0;
        double total2 = 0.0;
        total += qi * T[ti(l, m, n)];
        total += dxi * T[ti(l + 1, m, n)];
        total += dyi * T[ti(l, m + 1, n)];
        total += dzi * T[ti(l, m, n + 1)];
        total += qxxi * T[ti(l + 2, m, n)];
        total += qyyi * T[ti(l, m + 2, n)];
        total += qzzi * T[ti(l, m, n + 2)];
        total2 += qxyi * T[ti(l + 1, m + 1, n)];
        total2 += qxzi * T[ti(l + 1, m, n + 1)];
        total2 += qyzi * T[ti(l, m + 1, n + 1)];
        return total + 2.0 * total2;
    }

    /**
     * Contract multipole moments with their respective electrostatic potential
     * derivatives.
     *
     * @param T  array of electrostatic potential and partial derivatives
     * @param l  apply (d/dx)^l to the potential
     * @param m  apply (d/dy)^l to the potential
     * @param n  apply (d/dz)^l to the potential
     * @param sb the code will be appended to the StringBuilfer.
     * @return the contracted interaction.
     */
    private double codeContract(double[] T, int l, int m, int n, StringBuilder sb) {
        double total = 0.0;
        String name = term(l, m, n);
        sb.append(format("double %s = 0.0;\n", name));
        StringBuilder sb1 = new StringBuilder();
        double term = qi * T[ti(l, m, n)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q000i * T[%s];\n", name, tlmn(l, m, n)));
        }
        term = dxi * T[ti(l + 1, m, n)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q100i * T[%s];\n", name, tlmn(l + 1, m, n)));
        }
        term = dyi * T[ti(l, m + 1, n)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q010i * T[%s];\n", name, tlmn(l, m + 1, n)));
        }
        term = dzi * T[ti(l, m, n + 1)];
        if (term != 0) {
            total += term;
            sb1.append(format("%s += q001i * T[%s];\n", name, tlmn(l, m, n + 1)));
        }
        StringBuilder traceSB = new StringBuilder();
        double trace = 0.0;
        term = qxxi * T[ti(l + 2, m, n)];
        if (term != 0) {
            trace += term;
            // logger.info(format(" Qxx: %16.15f T: %16.15f Term: %16.15f", q200, T[ti(l + 2, m, n)], term));
            traceSB.append(format("%s += q200i * T[%s];\n", name, tlmn(l + 2, m, n)));
        }
        term = qyyi * T[ti(l, m + 2, n)];
        if (term != 0) {
            trace += term;
            // logger.info(format(" Qyy: %16.15f T: %16.15f Term: %16.15f", q020, T[ti(l, m + 2, n)], term));
            traceSB.append(format("%s += q020i * T[%s];\n", name, tlmn(l, m + 2, n)));
        }
        term = qzzi * T[ti(l, m, n + 2)];
        if (term != 0) {
            trace += term;
            // logger.info(format(" Qzz: %16.15f T: %16.15f Term: %16.15f", q002, T[ti(l, m, n + 2)], term));
            traceSB.append(format("%s += q002i * T[%s];\n", name, tlmn(l, m, n + 2)));
        }
        total += trace;
        if (total != 0) {
            sb.append(sb1.toString());
            if (trace != 0) {
                //logger.info(format(" Trace: %16.15f", trace));
                sb.append(traceSB);
            }
        }
        StringBuilder sb2 = new StringBuilder();
        double total2 = 0.0;
        term = qxyi * T[ti(l + 1, m + 1, n)];
        if (term != 0) {
            total2 += term;
            sb2.append(format("%s2 += q110i * T[%s];\n", name, tlmn(l + 1, m + 1, n)));
        }
        term = qxzi * T[ti(l + 1, m, n + 1)];
        if (term != 0) {
            total2 += term;
            sb2.append(format("%s2 += q101i * T[%s];\n", name, tlmn(l + 1, m, n + 1)));
        }
        term = qyzi * T[ti(l, m + 1, n + 1)];
        if (term != 0) {
            total2 += term;
            sb2.append(format("%s2 += q011i * T[%s];\n", name, tlmn(l, m + 1, n + 1)));
        }
        if (total2 != 0.0) {
            sb.append(format("double %s2 = 0.0;\n", name));
            sb.append(sb2);
            total += 2.0 * total2;
            sb.append(format("%s += 2.0 * %s2;\n", name, name));
        }
        return total;
    }

    /**
     * Collect the field at R due to Q multipole moments at the origin.
     *
     * @param T Electrostatic potential and partial derivatives
     * @param l apply (d/dx)^l to the potential
     * @param m apply (d/dy)^l to the potential
     * @param n apply (d/dz)^l to the potential
     */
    final void field(double[] T, int l, int m, int n) {
        E000 = contract(T, l, m, n);
        E100 = contract(T, l + 1, m, n);
        E010 = contract(T, l, m + 1, n);
        E001 = contract(T, l, m, n + 1);
        E200 = contract(T, l + 2, m, n);
        E020 = contract(T, l, m + 2, n);
        E002 = contract(T, l, m, n + 2);
        E110 = contract(T, l + 1, n + 1, m);
        E101 = contract(T, l + 1, m, n + 1);
        E011 = contract(T, l, m + 1, n + 1);
    }

    /**
     * Collect the field at R due to Q multipole moments at the origin.
     *
     * @param T  Electrostatic potential and partial derivatives.
     * @param l  apply (d/dx)^l to the potential.
     * @param m  apply (d/dy)^l to the potential.
     * @param n  apply (d/dz)^l to the potential.
     * @param sb Append the code to the StringBuilder.
     */
    private void codeField(double[] T, int l, int m, int n, StringBuilder sb) {
        E000 = codeContract(T, l, m, n, sb);
        if (E000 != 0) {
            sb.append(format("e000 = %s;\n", term(l, m, n)));
        }
        E100 = codeContract(T, l + 1, m, n, sb);
        if (E100 != 0) {
            sb.append(format("e100 = %s;\n", term(l + 1, m, n)));
        }
        E010 = codeContract(T, l, m + 1, n, sb);
        if (E100 != 0) {
            sb.append(format("e010 = %s;\n", term(l, m + 1, n)));
        }
        E001 = codeContract(T, l, m, n + 1, sb);
        if (E001 != 0) {
            sb.append(format("e001 = %s;\n", term(l, m, n + 1)));
        }
        E200 = codeContract(T, l + 2, m, n, sb);
        if (E200 != 0) {
            sb.append(format("e200 = %s;\n", term(l + 2, m, n)));
        }
        E020 = codeContract(T, l, m + 2, n, sb);
        if (E020 != 0) {
            sb.append(format("e020 = %s;\n", term(l, m + 2, n)));
        }
        E002 = codeContract(T, l, m, n + 2, sb);
        if (E002 != 0) {
            sb.append(format("e002 = %s;\n", term(l, m, n + 2)));
        }
        E110 = codeContract(T, l + 1, m + 1, n, sb);
        if (E110 != 0) {
            sb.append(format("e110 = %s;\n", term(l + 1, m + 1, n)));
        }
        E101 = codeContract(T, l + 1, m, n + 1, sb);
        if (E101 != 0) {
            sb.append(format("e101 = %s;\n", term(l + 1, m, n + 1)));
        }
        E011 = codeContract(T, l, m + 1, n + 1, sb);
        if (E011 != 0) {
            sb.append(format("e011 = %s;\n", term(l, m + 1, n + 1)));
        }
    }

    /**
     * Code a 5th order interaction.
     *
     * @param r  Separation vector.
     * @param Qi Multipole for site i.
     * @param Qk Multipole for site k.
     * @param Fi Force for site i.
     * @param Fk Force for site k.
     * @param Ti Torque for site i.
     * @param Tk Torque for site k.
     * @return Returns 0.0.
     */
    private double codeInteract5(double[] r, double[] Qi, double[] Qk,
                                 double[] Fi, double[] Fk, double[] Ti, double[] Tk) {
        double[] T = new double[tensorCount(5)];
        recursion(r, T);

        setMultipoleI(Qi);
        setMultipoleK(Qk);

        String whetherQI = (MultipoleTensorQI.class.isAssignableFrom(getClass())) ? "QI" : "";
        StringBuilder sb = new StringBuilder("\n\npublic void E" + whetherQI + "5(double T[]) {\n");
        codeField(T, 0, 0, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void Ex").append(whetherQI).append("5(double T[]) {\n");
        codeField(T, 1, 0, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void Ey").append(whetherQI).append("5(double T[]) {\n");
        codeField(T, 0, 1, 0, sb);
        sb.append("}\n");

        sb.append("\n\npublic void Ez").append(whetherQI).append("5(double T[]) {\n");
        codeField(T, 0, 0, 1, sb);
        sb.append("}\n");
        logger.log(Level.INFO, sb.toString());

        return 0.0;
    }

    /**
     * Just the energy portion of multipoleEnergyQI().
     *
     * @return a double.
     */
    public final double mpoleIFieldDotK() {
        multipoleIField();
        return dotMultipoleK();
    }

    /**
     * <p>mpoleKFieldDotI.</p>
     *
     * @return a double.
     */
    public final double mpoleKFieldDotI() {
        multipoleKField();
        return dotMultipoleI();
    }

    /**
     * inducedIField_qi + inducedIFieldCR_qi -&gt; dotMultipoleK
     *
     * @param scaleDipoleCR a double.
     * @param scaleDipole   a double.
     * @return a double.
     */
    public final double indiFieldsDotK(double scaleDipoleCR, double scaleDipole) {
        return qk * -(uzi * scaleDipole + pzi * scaleDipoleCR) * R001
                + dxk * -(uxi * scaleDipole + pxi * scaleDipoleCR) * R200
                + dyk * -(uyi * scaleDipole + pyi * scaleDipoleCR) * R020
                + dzk * -(uzi * scaleDipole + pzi * scaleDipoleCR) * R002
                + qxxk * -(uzi * scaleDipole + pzi * scaleDipoleCR) * R201
                + qyyk * -(uzi * scaleDipole + pzi * scaleDipoleCR) * R021
                + qzzk * -(uzi * scaleDipole + pzi * scaleDipoleCR) * R003
                + qxzk * -(uxi * scaleDipole + pxi * scaleDipoleCR) * R201
                + qyzk * -(uyi * scaleDipole + pyi * scaleDipoleCR) * R021;
    }

    /**
     * inducedKField_qi + inducedKFieldCR_qi -&gt; dotMultipoleI
     *
     * @param scaleDipoleCR a double.
     * @param scaleDipole   a double.
     * @return a double.
     */
    public final double indkFieldsDotI(double scaleDipoleCR, double scaleDipole) {
        return qi * (uzk * scaleDipole + pzk * scaleDipoleCR) * R001
                - dxi * (uxk * scaleDipole + pxk * scaleDipoleCR) * R200
                - dyi * (uyk * scaleDipole + pyk * scaleDipoleCR) * R020
                - dzi * (uzk * scaleDipole + pzk * scaleDipoleCR) * R002
                + qxxi * (uzk * scaleDipole + pzk * scaleDipoleCR) * R201
                + qyyi * (uzk * scaleDipole + pzk * scaleDipoleCR) * R021
                + qzzi * (uzk * scaleDipole + pzk * scaleDipoleCR) * R003
                + qxzi * (uxk * scaleDipole + pxk * scaleDipoleCR) * R201
                + qyzi * (uyk * scaleDipole + pyk * scaleDipoleCR) * R021;
    }

    /**
     * <p>multipoleITorque.</p>
     *
     * @param torque an array of {@link double} objects.
     */
    final void multipoleITorque(double[] torque) {
        // Torque on dipole moments due to the field.
        double dx = dyi * E001 - dzi * E010;
        double dy = dzi * E100 - dxi * E001;
        double dz = dxi * E010 - dyi * E100;

        // Torque on quadrupole moments due to the gradient of the field.
        double qx = qxyi * E101 + 2.0 * qyyi * E011 + qyzi * E002
                - (qxzi * E110 + qyzi * E020 + 2.0 * qzzi * E011);
        double qy = qxzi * E200 + qyzi * E110 + 2.0 * qzzi * E101
                - (2.0 * qxxi * E101 + qxyi * E011 + qxzi * E002);
        double qz = 2.0 * qxxi * E110 + qxyi * E020 + qxzi * E011
                - (qxyi * E200 + 2.0 * qyyi * E110 + qyzi * E101);

        torque[0] = dx - qx;
        torque[1] = dy - qy;
        torque[2] = dz - qz;
    }

    /**
     * <p>multipoleKTorque.</p>
     *
     * @param torque an array of {@link double} objects.
     */
    final void multipoleKTorque(double[] torque) {
        // Torque on dipole moments due to the field.
        double dx = dyk * E001 - dzk * E010;
        double dy = dzk * E100 - dxk * E001;
        double dz = dxk * E010 - dyk * E100;

        // Torque on quadrupole moments due to the gradkent of the fkeld.
        double qx = qxyk * E101 + 2.0 * qyyk * E011 + qyzk * E002
                - (qxzk * E110 + qyzk * E020 + 2.0 * qzzk * E011);
        double qy = qxzk * E200 + qyzk * E110 + 2.0 * qzzk * E101
                - (2.0 * qxxk * E101 + qxyk * E011 + qxzk * E002);
        double qz = 2.0 * qxxk * E110 + qxyk * E020 + qxzk * E011
                - (qxyk * E200 + 2.0 * qyyk * E110 + qyzk * E101);

        torque[0] = -(dx + qx);
        torque[1] = -(dy + qy);
        torque[2] = -(dz + qz);
    }

    /**
     * <p>dotMultipoleI.</p>
     *
     * @return a double.
     */
    final double dotMultipoleI() {
        double total = qi * E000;
        total -= dxi * E100;
        total -= dyi * E010;
        total -= dzi * E001;
        total += qxxi * E200;
        total += qyyi * E020;
        total += qzzi * E002;
        total += qxyi * E110;
        total += qxzi * E101;
        total += qyzi * E011;
        return total;
    }

    /**
     * <p>dotMultipoleK.</p>
     *
     * @return a double.
     */
    final double dotMultipoleK() {
        double total = qk * E000;
        total += dxk * E100;
        total += dyk * E010;
        total += dzk * E001;
        total += qxxk * E200;
        total += qyyk * E020;
        total += qzzk * E002;
        total += qxyk * E110;
        total += qxzk * E101;
        total += qyzk * E011;
        return total;
    }

    /**
     * Convenience method for writing out intermediate terms in the recursion.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @param j the jth intermediate term.
     * @return a String of the form <code>termlmnj</code>.
     */
    protected static String term(int l, int m, int n, int j) {
        return format("term%d%d%d%d", l, m, n, j);
    }

    /**
     * Convenience method for writing out intermediate terms in the recursion.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @return a String of the form <code>termlmnj</code>.
     */
    protected static String term(int l, int m, int n) {
        return format("term%d%d%d", l, m, n);
    }

    /**
     * Convenience method for writing out tensor indices.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @return a String of the form <code>tlmn</code>.
     */
    private static String tlmn(int l, int m, int n) {
        return format("t%d%d%d", l, m, n);
    }

    /**
     * Convenience method for writing out tensor indices.
     *
     * @param l number of d/dx partial derivatives.
     * @param m number of d/dx partial derivatives.
     * @param n number of d/dx partial derivatives.
     * @return a String of the form <code>Rlmn</code>.
     */
    protected static String rlmn(int l, int m, int n) {
        return format("R%d%d%d", l, m, n);
    }

    /**
     * <p>getTensor.</p>
     *
     * @param T an array of {@link double} objects.
     */
    final void getTensor(double[] T) {
        if (T == null || order < 0) {
            return;
        }
        // l + m + n = 0 (1)
        T[t000] = R000;
        if (order < 1) {
            return;
        }
        // l + m + n = 1 (3)   4
        T[t100] = R100;
        T[t010] = R010;
        T[t001] = R001;
        if (order < 2) {
            return;
        }
        // l + m + n = 2 (6)  10
        T[t200] = R200;
        T[t020] = R020;
        T[t002] = R002;
        T[t110] = R110;
        T[t101] = R101;
        T[t011] = R011;
        if (order < 3) {
            return;
        }
        // l + m + n = 3 (10) 20
        T[t300] = R300;
        T[t030] = R030;
        T[t003] = R003;
        T[t210] = R210;
        T[t201] = R201;
        T[t120] = R120;
        T[t021] = R021;
        T[t102] = R102;
        T[t012] = R012;
        T[t111] = R111;
        if (order < 4) {
            return;
        }
        // l + m + n = 4 (15) 35
        T[t400] = R400;
        T[t040] = R040;
        T[t004] = R004;
        T[t310] = R310;
        T[t301] = R301;
        T[t130] = R130;
        T[t031] = R031;
        T[t103] = R103;
        T[t013] = R013;
        T[t220] = R220;
        T[t202] = R202;
        T[t022] = R022;
        T[t211] = R211;
        T[t121] = R121;
        T[t112] = R112;
        if (order < 5) {
            return;
        }
        // l + m + n = 5 (21) 56
        T[t500] = R500;
        T[t050] = R050;
        T[t005] = R005;
        T[t410] = R410;
        T[t401] = R401;
        T[t140] = R140;
        T[t041] = R041;
        T[t104] = R104;
        T[t014] = R014;
        T[t320] = R320;
        T[t302] = R302;
        T[t230] = R230;
        T[t032] = R032;
        T[t203] = R203;
        T[t023] = R023;
        T[t311] = R311;
        T[t131] = R131;
        T[t113] = R113;
        T[t221] = R221;
        T[t212] = R212;
        T[t122] = R122;
    }

    /**
     * <p>setTensor.</p>
     *
     * @param T an array of {@link double} objects.
     */
    final void setTensor(double[] T) {
        if (T == null || order < 0) {
            return;
        }
        // l + m + n = 0 (1)
        R000 = T[t000];
        if (order < 1) {
            return;
        }
        // l + m + n = 1 (3)   4
        R100 = T[t100];
        R010 = T[t010];
        R001 = T[t001];
        if (order < 2) {
            return;
        }
        // l + m + n = 2 (6)  10
        R200 = T[t200];
        R020 = T[t020];
        R002 = T[t002];
        R110 = T[t110];
        R101 = T[t101];
        R011 = T[t011];
        if (order < 3) {
            return;
        }
        // l + m + n = 3 (10) 20
        R300 = T[t300];
        R030 = T[t030];
        R003 = T[t003];
        R210 = T[t210];
        R201 = T[t201];
        R120 = T[t120];
        R021 = T[t021];
        R102 = T[t102];
        R012 = T[t012];
        R111 = T[t111];
        if (order < 4) {
            return;
        }
        // l + m + n = 4 (15) 35
        R400 = T[t400];
        R040 = T[t040];
        R004 = T[t004];
        R310 = T[t310];
        R301 = T[t301];
        R130 = T[t130];
        R031 = T[t031];
        R103 = T[t103];
        R013 = T[t013];
        R220 = T[t220];
        R202 = T[t202];
        R022 = T[t022];
        R211 = T[t211];
        R121 = T[t121];
        R112 = T[t112];
        if (order < 5) {
            return;
        }
        // l + m + n = 5 (21) 56
        R500 = T[t500];
        R050 = T[t050];
        R005 = T[t005];
        R410 = T[t410];
        R401 = T[t401];
        R140 = T[t140];
        R041 = T[t041];
        R104 = T[t104];
        R014 = T[t014];
        R320 = T[t320];
        R302 = T[t302];
        R230 = T[t230];
        R032 = T[t032];
        R203 = T[t203];
        R023 = T[t023];
        R311 = T[t311];
        R131 = T[t131];
        R113 = T[t113];
        R221 = T[t221];
        R212 = T[t212];
        R122 = T[t122];
    }

    /**
     * <p>setMultipoles.</p>
     *
     * @param Qi an array of {@link double} objects.
     * @param Qk an array of {@link double} objects.
     */
    private void setMultipoles(double[] Qi, double[] Qk) {
        setMultipoleI(Qi);
        setMultipoleK(Qk);
    }

    /**
     * <p>setDipoles.</p>
     *
     * @param ui   an array of {@link double} objects.
     * @param uiCR an array of {@link double} objects.
     * @param uk   an array of {@link double} objects.
     * @param ukCR an array of {@link double} objects.
     */
    private void setDipoles(double[] ui, double[] uiCR, double[] uk, double[] ukCR) {
        setDipoleI(ui, uiCR);
        setDipoleK(uk, ukCR);
    }

    /**
     * <p>multipoleEnergy.</p>
     *
     * @param Fi an array of {@link double} objects.
     * @param Ti an array of {@link double} objects.
     * @param Tk an array of {@link double} objects.
     * @return a double.
     */
    public abstract double multipoleEnergy(double[] Fi, double[] Ti, double[] Tk);

    /**
     * <p>polarizationEnergy.</p>
     *
     * @param scaleField  a double.
     * @param scaleEnergy a double.
     * @param scaleMutual a double.
     * @param Fi          an array of {@link double} objects.
     * @param Ti          an array of {@link double} objects.
     * @param Tk          an array of {@link double} objects.
     * @return a double.
     */
    public abstract double polarizationEnergy(double scaleField, double scaleEnergy, double scaleMutual,
                                              double[] Fi, double[] Ti, double[] Tk);

    /**
     * Is a no-op when not using QI.
     *
     * @param Fi an array of {@link double} objects.
     * @param Ti an array of {@link double} objects.
     * @param Tk an array of {@link double} objects.
     */
    protected abstract void qiToGlobal(double[] Fi, double[] Ti, double[] Tk);

    /**
     * <p>scaleInduced.</p>
     *
     * @param scaleField  a double.
     * @param scaleEnergy a double.
     */
    final void scaleInduced(double scaleField, double scaleEnergy) {
        uxi *= scaleEnergy;
        uyi *= scaleEnergy;
        uzi *= scaleEnergy;
        pxi *= scaleField;
        pyi *= scaleField;
        pzi *= scaleField;
        sxi = 0.5 * (uxi + pxi);
        syi = 0.5 * (uyi + pyi);
        szi = 0.5 * (uzi + pzi);
        uxk *= scaleEnergy;
        uyk *= scaleEnergy;
        uzk *= scaleEnergy;
        pxk *= scaleField;
        pyk *= scaleField;
        pzk *= scaleField;
        sxk = 0.5 * (uxk + pxk);
        syk = 0.5 * (uyk + pyk);
        szk = 0.5 * (uzk + pzk);
    }

    /**
     * This method is a driver to collect elements of the Cartesian multipole
     * tensor given the recursion relationships implemented by the method
     * "Tlmnj", which can be called directly to get a single tensor element. It
     * does not store intermediate values of the recursion, causing it to scale
     * O(order^8). For order = 5, this approach is a factor of 10 slower than
     * recursion.
     *
     * @param r      double[] vector between two sites.
     * @param tensor double[] length must be at least binomial(order + 3, 3).
     */
    protected abstract void noStorageRecursion(double[] r, double[] tensor);

    /**
     * This routine implements the recurrence relations for computation of any
     * Cartesian multipole tensor in ~O(L^8) time, where L is the total order l
     * + m + n, given the auxiliary elements T0000.
     * <br>
     * It implements the recursion relationships in brute force fashion, without
     * saving intermediate values. This is useful for finding a single tensor,
     * rather than all binomial(L + 3, 3).
     * <br>
     * The specific recursion equations (41-43) and set of auxiliary tensor
     * elements from equation (40) can be found in Challacombe et al.
     *
     * @param l    int The number of (d/dx) operations.
     * @param m    int The number of (d/dy) operations.
     * @param n    int The number of (d/dz) operations.
     * @param j    int j = 0 is the Tlmn tensor, j .GT. 0 is an intermediate.
     * @param r    double[] The {x,y,z} coordinates.
     * @param T000 double[] Initial auxiliary tensor elements from Eq. (40).
     * @return double The requested Tensor element (intermediate if j .GT. 0).
     * @since 1.0
     */
    protected abstract double Tlmnj(final int l, final int m, final int n,
                                    final int j, final double[] r, final double[] T000);

    /**
     * <p>recursion.</p>
     *
     * @param r      an array of {@link double} objects.
     * @param tensor an array of {@link double} objects.
     */
    protected abstract void recursion(final double[] r, final double[] tensor);

    /**
     * <p>setMultipoleI.</p>
     *
     * @param Qi an array of {@link double} objects.
     */
    protected abstract void setMultipoleI(double[] Qi);

    /**
     * <p>setMultipoleK.</p>
     *
     * @param Qk an array of {@link double} objects.
     */
    protected abstract void setMultipoleK(double[] Qk);

    /**
     * <p>setDipoleI.</p>
     *
     * @param ui   an array of {@link double} objects.
     * @param uiCR an array of {@link double} objects.
     */
    protected abstract void setDipoleI(double[] ui, double[] uiCR);

    /**
     * <p>setDipoleK.</p>
     *
     * @param uk   an array of {@link double} objects.
     * @param ukCR an array of {@link double} objects.
     */
    protected abstract void setDipoleK(double[] uk, double[] ukCR);

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 4th
     * order, in the global frame, which is sufficient for quadrupole-induced
     * dipole forces.
     */
    protected abstract void order4();

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 5th
     * order, in the global frame, which is sufficient for quadrupole-quadrupole
     * forces.
     */
    protected abstract void order5();

    /**
     * Hard coded computation of all Cartesian multipole tensors up to 5th
     * order, in the global frame, which is sufficient for quadrupole-quadrupole
     * forces and orthogonal space sampling.
     */
    protected abstract void order6();

    /**
     * <p>multipoleIField.</p>
     */
    protected abstract void multipoleIField();

    /**
     * <p>multipoleKField.</p>
     */
    protected abstract void multipoleKField();

    /**
     * <p>multipoleIdX.</p>
     */
    protected abstract void multipoleIdX();

    /**
     * <p>multipoleIdY.</p>
     */
    protected abstract void multipoleIdY();

    /**
     * <p>multipoleIdZ.</p>
     */
    protected abstract void multipoleIdZ();

    /**
     * Never used by Global coordinates; necessary only for lambda derivatives.
     */
    protected abstract void multipoleIdZ2();

    /**
     * <p>inducedIField.</p>
     */
    protected abstract void inducedIField();

    /**
     * <p>inducedKField.</p>
     */
    protected abstract void inducedKField();

    /**
     * <p>inducedIFieldCR.</p>
     */
    protected abstract void inducedIFieldCR();

    /**
     * <p>inducedKFieldCR.</p>
     */
    protected abstract void inducedKFieldCR();

    /**
     * <p>inducedIdX.</p>
     */
    protected abstract void inducedIdX();

    /**
     * <p>inducedIdY.</p>
     */
    protected abstract void inducedIdY();

    /**
     * <p>inducedIdZ.</p>
     */
    protected abstract void inducedIdZ();

    /**
     * <p>inducedKdX.</p>
     */
    protected abstract void inducedKdX();

    /**
     * <p>inducedKdY.</p>
     */
    protected abstract void inducedKdY();

    /**
     * <p>inducedKdZ.</p>
     */
    protected abstract void inducedKdZ();

    /**
     * <p>getdEdZ.</p>
     *
     * @return a double.
     */
    public abstract double getdEdZ();

    /**
     * <p>getd2EdZ2.</p>
     *
     * @return a double.
     */
    public abstract double getd2EdZ2();

    /**
     * Constant <code>oneThird=1.0 / 3.0</code>
     */
    protected static final double oneThird = 1.0 / 3.0;
    /**
     * Constant <code>twoThirds=2.0 / 3.0</code>
     */
    protected static final double twoThirds = 2.0 / 3.0;
    /**
     * Constant <code>threeFifths=3.0 / 5.0</code>
     */
    private static final double threeFifths = 3.0 / 5.0;
    /**
     * Constant <code>oneThirtyFifth=1.0 / 35.0</code>
     */
    private static final double oneThirtyFifth = 1.0 / 35.0;
    /**
     * Constant <code>sqrtPI = sqrt(PI)</code>
     */
    private static final double sqrtPI = sqrt(PI);

    // Multipole components for atom i.
    protected double qi;
    protected double dxi;
    protected double dyi;
    protected double dzi;
    protected double qxxi;
    protected double qyyi;
    protected double qzzi;
    protected double qxyi;
    protected double qxzi;
    protected double qyzi;

    // Induced dipole components for atom i.
    protected double uxi;
    protected double uyi;
    protected double uzi;

    // Induced dipole CR components for atom i.
    protected double pxi;
    protected double pyi;
    protected double pzi;

    // Induced dipole + Induced dipole CR components for atom i.
    protected double sxi;
    protected double syi;
    protected double szi;

    // Multipole components for atom k.
    protected double qk;
    protected double dxk;
    protected double dyk;
    protected double dzk;
    protected double qxxk;
    protected double qyyk;
    protected double qzzk;
    protected double qxyk;
    protected double qxzk;
    protected double qyzk;

    // Induced dipole components for atom k.
    protected double uxk;
    protected double uyk;
    protected double uzk;

    // Induced dipole CR components for atom k.
    protected double pxk;
    protected double pyk;
    protected double pzk;

    // Induced dipole + Induced dipole CR components for atom k.
    protected double sxk;
    protected double syk;
    protected double szk;

    // Components of the potential, field and field gradient.
    double E000; // Potential
    double E100; // X Component of the Field
    double E010; // Y Component of the Field
    double E001; // Z Component of the Field
    double E200; // XX Component of the Field Gradient
    double E020; // YY Component of the Field Gradient
    double E002; // ZZ Component of the Field Gradient
    double E110; // XY Component of the Field Gradient
    double E101; // XZ Component of the Field Gradient
    double E011; // YZ Component of the Field Gradient

    // Cartesian tensor elements (for 1/R, erfc(Beta*R)/R or Thole damping.
    // l + m + n = 0 (1)
    double R000;
    // l + m + n = 1 (3)   4
    double R100;
    double R010;
    double R001;
    // l + m + n = 2 (6)  10
    double R200;
    double R020;
    double R002;
    double R110;
    double R101;
    double R011;
    // l + m + n = 3 (10) 20
    double R300;
    double R030;
    double R003;
    double R210;
    double R201;
    double R120;
    double R021;
    double R102;
    double R012;
    double R111;
    // l + m + n = 4 (15) 35
    double R400;
    double R040;
    double R004;
    double R310;
    double R301;
    double R130;
    double R031;
    double R103;
    double R013;
    double R220;
    double R202;
    double R022;
    double R211;
    double R121;
    double R112;
    // l + m + n = 5 (21) 56
    double R500;
    double R050;
    double R005;
    double R410;
    double R401;
    double R140;
    double R041;
    double R104;
    double R014;
    double R320;
    double R302;
    double R230;
    double R032;
    double R203;
    double R023;
    double R311;
    double R131;
    double R113;
    double R221;
    double R212;
    double R122;
    // l + m + n = 6 (28) 84
    double R006;
    double R402;
    double R042;
    double R204;
    double R024;
    double R222;
    double R600;
    double R060;
    double R510;
    double R501;
    double R150;
    double R051;
    double R105;
    double R015;
    double R420;
    double R240;
    double R411;
    double R141;
    double R114;
    double R330;
    double R303;
    double R033;
    double R321;
    double R231;
    double R213;
    double R312;
    double R132;
    double R123;

    // l + m + n = 0 (1)
    protected final int t000;
    // l + m + n = 1 (3)   4
    protected final int t100;
    protected final int t010;
    protected final int t001;
    // l + m + n = 2 (6)  10
    protected final int t200;
    protected final int t020;
    protected final int t002;
    protected final int t110;
    protected final int t101;
    protected final int t011;
    // l + m + n = 3 (10) 20
    protected final int t300;
    protected final int t030;
    protected final int t003;
    protected final int t210;
    protected final int t201;
    protected final int t120;
    protected final int t021;
    protected final int t102;
    protected final int t012;
    protected final int t111;
    // l + m + n = 4 (15) 35
    private final int t400;
    private final int t040;
    private final int t004;
    private final int t310;
    private final int t301;
    private final int t130;
    private final int t031;
    private final int t103;
    private final int t013;
    private final int t220;
    private final int t202;
    private final int t022;
    private final int t211;
    private final int t121;
    private final int t112;
    // l + m + n = 5 (21) 56
    private final int t500;
    private final int t050;
    private final int t005;
    private final int t410;
    private final int t401;
    private final int t140;
    private final int t041;
    private final int t104;
    private final int t014;
    private final int t320;
    private final int t302;
    private final int t230;
    private final int t032;
    private final int t203;
    private final int t023;
    private final int t311;
    private final int t131;
    private final int t113;
    private final int t221;
    private final int t212;
    private final int t122;
    // l + m + n = 6 (28) 84
    private final int t600;
    private final int t060;
    private final int t006;
    private final int t510;
    private final int t501;
    private final int t150;
    private final int t051;
    private final int t105;
    private final int t015;
    private final int t420;
    private final int t402;
    private final int t240;
    private final int t042;
    private final int t204;
    private final int t024;
    private final int t411;
    private final int t141;
    private final int t114;
    private final int t330;
    private final int t303;
    private final int t033;
    private final int t321;
    private final int t231;
    private final int t213;
    private final int t312;
    private final int t132;
    private final int t123;
    private final int t222;
}
