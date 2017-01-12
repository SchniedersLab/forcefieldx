/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.numerics;

import org.apache.commons.math3.util.FastMath;

import static org.apache.commons.math3.util.FastMath.atan2;
import static org.apache.commons.math3.util.FastMath.cosh;
import static org.apache.commons.math3.util.FastMath.hypot;
import static org.apache.commons.math3.util.FastMath.sinh;

/**
 * <p>
 * ComplexNumber class.</p>
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
 */
public class ComplexNumber {

    private double re;
    private double im;

    /**
     * <p>
     * Constructor for ComplexNumber.</p>
     */
    public ComplexNumber() {
    }

    /**
     * <p>
     * Constructor for ComplexNumber.</p>
     *
     * @param real a double.
     * @param imag a double.
     */
    public ComplexNumber(double real, double imag) {
        re = real;
        im = imag;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (im == 0) {
            return re + "";
        }
        if (re == 0) {
            return im + "i";
        }
        if (im < 0) {
            return re + " - " + (-im) + "i";
        }
        return re + " + " + im + "i";
    }

    /**
     * <p>
     * re</p>
     *
     * @return a double.
     */
    public double re() {
        return re;
    }

    /**
     * <p>
     * re</p>
     *
     * @param re a double.
     */
    public void re(double re) {
        this.re = re;
    }

    /**
     * <p>
     * im</p>
     *
     * @return a double.
     */
    public double im() {
        return im;
    }

    /**
     * <p>
     * im</p>
     *
     * @param im a double.
     */
    public void im(double im) {
        this.im = im;
    }

    /**
     * <p>
     * copy</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void copy(ComplexNumber b) {
        ComplexNumber a = this;
        a.re = b.re;
        a.im = b.im;
    }

    /**
     * <p>
     * abs</p>
     *
     * @return a double.
     */
    public double abs() {
        return hypot(re, im);
    }

    /**
     * <p>
     * phase</p>
     *
     * @return a double.
     */
    public double phase() {
        return atan2(im, re);
    }

    /**
     * <p>
     * phaseShift</p>
     *
     * @param s a double.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber phaseShift(double s) {
        ComplexNumber sc = new ComplexNumber(FastMath.cos(s), FastMath.sin(s));
        return this.times(sc);
    }

    /**
     * <p>
     * phaseShiftIP</p>
     *
     * @param s a double.
     */
    public void phaseShiftIP(double s) {
        ComplexNumber a = this;
        double sr = FastMath.cos(s);
        double si = FastMath.sin(s);
        double real = a.re * sr - a.im * si;
        double imag = a.re * si + a.im * sr;
        a.re = real;
        a.im = imag;
    }

    /**
     * Static version of phaseShift.
     *
     * @param a a {@link ffx.numerics.ComplexNumber} object.
     * @param s a double.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public static ComplexNumber phaseShift(ComplexNumber a, double s) {
        ComplexNumber sc = new ComplexNumber(FastMath.cos(s), FastMath.sin(s));
        return a.times(sc);
    }

    /**
     * Return a new Complex object whose value is (this + b).
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber plus(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re + b.re;
        double imag = a.im + b.im;
        return new ComplexNumber(real, imag);
    }

    /**
     * <p>
     * plusIP</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void plusIP(ComplexNumber b) {
        ComplexNumber a = this;
        a.re += b.re;
        a.im += b.im;
    }

    /**
     * Return a new Complex object whose value is (this - b).
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber minus(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re - b.re;
        double imag = a.im - b.im;
        return new ComplexNumber(real, imag);
    }

    /**
     * <p>
     * minusIP</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void minusIP(ComplexNumber b) {
        ComplexNumber a = this;
        a.re -= b.re;
        a.im -= b.im;
    }

    /**
     * Return a new Complex object whose value is (this * b).
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber times(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        return new ComplexNumber(real, imag);
    }

    /**
     * <p>
     * timesIP</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void timesIP(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        a.re = real;
        a.im = imag;
    }

    /**
     * Return a new object whose value is (this * alpha).
     *
     * @param alpha a double.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber times(double alpha) {
        return new ComplexNumber(alpha * re, alpha * im);
    }

    /**
     * <p>
     * timesIP</p>
     *
     * @param alpha a double.
     */
    public void timesIP(double alpha) {
        ComplexNumber a = this;
        a.re *= alpha;
        a.im *= alpha;
    }

    /**
     *
     * Return a new Complex object whose value is the conjugate of this.
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber conjugate() {
        return new ComplexNumber(re, -im);
    }

    /**
     * <p>
     * conjugateIP</p>
     */
    public void conjugateIP() {
        this.im = -this.im;
    }

    /**
     * Return a new Complex object whose value is the reciprocal of this.
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber reciprocal() {
        double scale = re * re + im * im;
        double iscale = 1.0 / scale;
        return new ComplexNumber(re * iscale, -im * iscale);
    }

    /**
     * <p>
     * reciprocalIP</p>
     */
    public void reciprocalIP() {
        double scale = re * re + im * im;
        double iscale = 1.0 / scale;
        re *= iscale;
        im *= -iscale;
    }

    /**
     * Return a / b.
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber divides(ComplexNumber b) {
        ComplexNumber a = this;
        return a.times(b.reciprocal());
    }

    /**
     * Return a new Complex object whose value is the complex exponential of
     * this.
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber exp() {
        return new ComplexNumber(FastMath.exp(re) * FastMath.cos(im), FastMath.exp(re) * FastMath.sin(im));
    }

    /**
     * Return a new Complex object whose value is the complex sine of this.
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber sin() {
        return new ComplexNumber(FastMath.sin(re) * cosh(im), FastMath.cos(re) * sinh(im));
    }

    /**
     * Return a new Complex object whose value is the complex cosine of this.
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber cos() {
        return new ComplexNumber(FastMath.cos(re) * cosh(im), -FastMath.sin(re) * sinh(im));
    }

    /**
     * Return a new Complex object whose value is the complex tangent of this.
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber tan() {
        return sin().divides(cos());
    }
}
