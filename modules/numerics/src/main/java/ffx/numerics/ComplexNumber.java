/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
 *
 * @author fennt
 */
public class ComplexNumber {

    private double re;
    private double im;

    public ComplexNumber() {
    }

    public ComplexNumber(double real, double imag) {
        re = real;
        im = imag;
    }

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

    public double re() {
        return re;
    }

    public void re(double re) {
        this.re = re;
    }

    public double im() {
        return im;
    }

    public void im(double im) {
        this.im = im;
    }

    public double abs() {
        return hypot(re, im);
    }

    public double phase() {
        return atan2(im, re);
    }

    public ComplexNumber phase_shift(double s) {
        ComplexNumber sc = new ComplexNumber(Math.cos(s), Math.sin(s));
        return this.times(sc);
    }

    // static version of phase_shift
    public static ComplexNumber phase_shift(ComplexNumber a, double s) {
        ComplexNumber sc = new ComplexNumber(Math.cos(s), Math.sin(s));
        return a.times(sc);
    }

    // return a new Complex object whose value is (this + b)
    public ComplexNumber plus(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re + b.re;
        double imag = a.im + b.im;
        return new ComplexNumber(real, imag);
    }

    // a static version of plus
    public static ComplexNumber plus(ComplexNumber a, ComplexNumber b) {
        double real = a.re + b.re;
        double imag = a.im + b.im;
        ComplexNumber sum = new ComplexNumber(real, imag);
        return sum;
    }

    // return a new Complex object whose value is (this - b)
    public ComplexNumber minus(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re - b.re;
        double imag = a.im - b.im;
        return new ComplexNumber(real, imag);
    }

    // a static version of plus
    public static ComplexNumber minus(ComplexNumber a, ComplexNumber b) {
        double real = a.re - b.re;
        double imag = a.im - b.im;
        ComplexNumber sum = new ComplexNumber(real, imag);
        return sum;
    }

    // return a new Complex object whose value is (this * b)
    public ComplexNumber times(ComplexNumber b) {
        ComplexNumber a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        return new ComplexNumber(real, imag);
    }

    // a static version of times
    public ComplexNumber times(ComplexNumber a, ComplexNumber b) {
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        return new ComplexNumber(real, imag);
    }

    // return a new object whose value is (this * alpha)
    public ComplexNumber times(double alpha) {
        return new ComplexNumber(alpha * re, alpha * im);
    }

    // a static version of times
    public static ComplexNumber times(double alpha, ComplexNumber a) {
        return new ComplexNumber(alpha * a.re, alpha * a.im);
    }

    // return a new Complex object whose value is the conjugate of this
    public ComplexNumber conjugate() {
        return new ComplexNumber(re, -im);
    }

    // return a new Complex object whose value is the reciprocal of this
    public ComplexNumber reciprocal() {
        double scale = re * re + im * im;
        return new ComplexNumber(re / scale, -im / scale);
    }

    // return a / b
    public ComplexNumber divides(ComplexNumber b) {
        ComplexNumber a = this;
        return a.times(b.reciprocal());
    }

    // return a new Complex object whose value is the complex exponential of this
    public ComplexNumber exp() {
        return new ComplexNumber(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
    }

    // return a new Complex object whose value is the complex sine of this
    public ComplexNumber sin() {
        return new ComplexNumber(Math.sin(re) * cosh(im), Math.cos(re) * sinh(im));
    }

    // return a new Complex object whose value is the complex cosine of this
    public ComplexNumber cos() {
        return new ComplexNumber(Math.cos(re) * cosh(im), -Math.sin(re) * sinh(im));
    }

    // return a new Complex object whose value is the complex tangent of this
    public ComplexNumber tan() {
        return sin().divides(cos());
    }
}
