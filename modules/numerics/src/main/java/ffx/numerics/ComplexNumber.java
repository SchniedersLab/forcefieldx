/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.numerics;

import static java.lang.Math.atan2;
import static java.lang.Math.cosh;
import static java.lang.Math.hypot;
import static java.lang.Math.sinh;

/**
 * <p>ComplexNumber class.</p>
 *
 * @author fennt
 * @version $Id: $
 */
public class ComplexNumber {

    private double re;
    private double im;

    /**
     * <p>Constructor for ComplexNumber.</p>
     */
    public ComplexNumber() {
    }

    /**
     * <p>Constructor for ComplexNumber.</p>
     *
     * @param real a double.
     * @param imag a double.
     */
    public ComplexNumber(double real, double imag) {
        re = real;
        im = imag;
    }

    /** {@inheritDoc} */
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
     * <p>re</p>
     *
     * @return a double.
     */
    public double re() {
        return re;
    }

    /**
     * <p>re</p>
     *
     * @param re a double.
     */
    public void re(double re) {
        this.re = re;
    }

    /**
     * <p>im</p>
     *
     * @return a double.
     */
    public double im() {
        return im;
    }

    /**
     * <p>im</p>
     *
     * @param im a double.
     */
    public void im(double im) {
        this.im = im;
    }

    /**
     * <p>copy</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void copy(ComplexNumber b){
        ComplexNumber a = this;
        a.re = b.re;
        a.im = b.im;
    }

    /**
     * <p>abs</p>
     *
     * @return a double.
     */
    public double abs() {
        return hypot(re, im);
    }

    /**
     * <p>phase</p>
     *
     * @return a double.
     */
    public double phase() {
        return atan2(im, re);
    }

    /**
     * <p>phase_shift</p>
     *
     * @param s a double.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber phase_shift(double s) {
        ComplexNumber sc = new ComplexNumber(Math.cos(s), Math.sin(s));
        return this.times(sc);
    }

    /**
     * <p>phase_shift_ip</p>
     *
     * @param s a double.
     */
    public void phase_shift_ip(double s){
        ComplexNumber a = this;
        double sr = Math.cos(s);
        double si = Math.sin(s);
        double real = a.re * sr - a.im * si;
        double imag = a.re * si + a.im * sr;
        a.re = real;
        a.im = imag;
    }

    // static version of phase_shift
    /**
     * <p>phase_shift</p>
     *
     * @param a a {@link ffx.numerics.ComplexNumber} object.
     * @param s a double.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public static ComplexNumber phase_shift(ComplexNumber a, double s) {
        ComplexNumber sc = new ComplexNumber(Math.cos(s), Math.sin(s));
        return a.times(sc);
    }

    // return a new Complex object whose value is (this + b)
    /**
     * <p>plus</p>
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
     * <p>plus_ip</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void plus_ip(ComplexNumber b) {
        ComplexNumber a = this;
        a.re += b.re;
        a.im += b.im;
    }

    // return a new Complex object whose value is (this - b)
    /**
     * <p>minus</p>
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
     * <p>minus_ip</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void minus_ip(ComplexNumber b) {
        ComplexNumber a = this;
        a.re -= b.re;
        a.im -= b.im;
    }

    // return a new Complex object whose value is (this * b)
    /**
     * <p>times</p>
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
     * <p>times_ip</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     */
    public void times_ip(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        a.re = real;
        a.im = imag;
    }

    // return a new object whose value is (this * alpha)
    /**
     * <p>times</p>
     *
     * @param alpha a double.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber times(double alpha) {
        return new ComplexNumber(alpha * re, alpha * im);
    }

    /**
     * <p>times_ip</p>
     *
     * @param alpha a double.
     */
    public void times_ip(double alpha) {
        ComplexNumber a = this;
        a.re *= alpha;
        a.im *= alpha;
    }

    // return a new Complex object whose value is the conjugate of this
    /**
     * <p>conjugate</p>
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber conjugate() {
        return new ComplexNumber(re, -im);
    }

    /**
     * <p>conjugate_ip</p>
     */
    public void conjugate_ip() {
        this.im = -this.im;
    }

    // return a new Complex object whose value is the reciprocal of this
    /**
     * <p>reciprocal</p>
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber reciprocal() {
        double scale = re * re + im * im;
        return new ComplexNumber(re / scale, -im / scale);
    }

    /**
     * <p>reciprocal_ip</p>
     */
    public void reciprocal_ip(){
        ComplexNumber a = this;
        double scale = re * re + im * im;
        a.re = re / scale;
        a.im = -im / scale;
    }

    // return a / b
    /**
     * <p>divides</p>
     *
     * @param b a {@link ffx.numerics.ComplexNumber} object.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber divides(ComplexNumber b) {
        ComplexNumber a = this;
        return a.times(b.reciprocal());
    }

    // return a new Complex object whose value is the complex exponential of this
    /**
     * <p>exp</p>
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber exp() {
        return new ComplexNumber(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
    }

    // return a new Complex object whose value is the complex sine of this
    /**
     * <p>sin</p>
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber sin() {
        return new ComplexNumber(Math.sin(re) * cosh(im), Math.cos(re) * sinh(im));
    }

    // return a new Complex object whose value is the complex cosine of this
    /**
     * <p>cos</p>
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber cos() {
        return new ComplexNumber(Math.cos(re) * cosh(im), -Math.sin(re) * sinh(im));
    }

    // return a new Complex object whose value is the complex tangent of this
    /**
     * <p>tan</p>
     *
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber tan() {
        return sin().divides(cos());
    }
}
