/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.numerics.fft;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

/**
 * Compute the FFT of real, double precision data of arbitrary length n using a
 * Complex transform.
 *
 * @author Michal J. Schnieders<br> Derived from:
 * <br>
 * Bruce R. Miller bruce.miller@nist.gov
 * <br>
 * Contribution of the National Institute of Standards and Technology, not
 * subject to copyright.<br> Derived from:<br> GSL (Gnu Scientific Library) FFT
 * Code by Brian Gough bjg@vvv.lanl.gov
 *
 * @see
 * <ul>
 * <li>
 * Complex
 * </li>
 * <li>
 * <a href="http://dx.doi.org/10.1109/TASSP.1987.1165220" target="_blank">
 * Henrik V. Sorenson, Douglas L. Jones, Michael T. Heideman, and C. Sidney
 * Burrus. Real-valued fast fourier fft algorithms. IEEE Transactions on
 * Acoustics, Speech, and Signal Processing, ASSP-35(6):849–863, 1987.
 * </a>
 * </li>
 * <li>
 * <a href="http://www.jstor.org/stable/2003354" target="_blank">J. W. Cooley
 * and J. W. Tukey, Mathematics of Computation 19 (90), 297 (1965)
 * </a>
 * </li>
 * <li>
 * <a href="http://en.wikipedia.org/wiki/Fast_Fourier_transform"
 * target="_blank">FFT at Wikipedia
 * </a>
 * </li>
 * </ul>
 *
 * @since 1.0
 */
public class Real {

    private final Complex complexFFT;
    private final int n;
    private final int halfN;
    private final double theta;
    private final double cosTheta;
    private final double sinTheta;
    private final double work[];

    /**
     * Constructs a Complex FFT of length (n / 2) for real data of length n.
     *
     * @param n a int.
     */
    public Real(int n) {
        //assert (n % 2 == 0);
        this.n = n;
        halfN = n / 2;
        theta = PI / halfN;
        cosTheta = cos(theta);
        sinTheta = sin(theta);
        complexFFT = new Complex(halfN);
        work = new double[n + 2];
    }

    /**
     * <p>
     * fft</p>
     *
     * @param data an array of double.
     * @param offset a int.
     */
    public void fft(double data[], int offset) {
        complexFFT.fft(data, offset, 2);
        unpack(data, offset);
    }

    /**
     * <p>
     * ifft</p>
     *
     * @param data an array of double.
     * @param offset a int.
     */
    public void ifft(double data[], int offset) {
        pack(data, offset);
        complexFFT.ifft(data, offset, 2);
        // Renormalize the 1/2 length fft.
        int ii = offset;
        for (int i = 0; i < n; i++) {
            data[ii++] *= 2.0;
        }
    }

    /**
     * Return the normalization factor. Multiply the elements of the
     * back-transformed data to get the normalized inverse.
     *
     * @return a double.
     */
    public double normalization() {
        return 1.0 / n;
    }

    /**
     * <p>
     * inverse</p>
     *
     * @param data an array of double.
     * @param offset a int.
     */
    public void inverse(double data[], int offset) {
        ifft(data, offset);
        /**
         * normalize inverse FFT with 1/n.
         */
        double norm = normalization();
        for (int i = 0; i < n; i++) {
            final int index = offset + i;
            data[index] *= norm;
        }
    }

    private void unpack(double data[], int offset) {
        for (int i = 0; i < n; i++) {
            work[i] = data[i + offset];
        }
        double wrs = cosTheta;
        double wis = -sinTheta;
        final double d0 = work[0];
        final double d1 = work[1];
        final int on = offset + n;
        final int o1 = offset + 1;
        data[offset] = d0 + d1;
        data[o1] = 0.0;
        data[on] = d0 - d1;
        data[on + 1] = 0.0;
        double wr = wrs;
        double wi = wis;
        for (int i = 1; i < halfN; i++) {
            int i1 = 2 * i;
            double rk = work[i1];
            double ik = work[i1 + 1];
            int i2 = n - 2 * i;
            double rn = work[i2];
            double in = work[i2 + 1];
            double s1r = 0.5 * (rk + rn);
            double s1i = 0.5 * (ik - in);
            double s2r = 0.5 * (ik + in);
            double s2i = 0.5 * (rn - rk);
            double fr = s1r + wr * s2r - wi * s2i;
            double fi = s1i + wi * s2r + wr * s2i;
            data[i1 + offset] = fr;
            data[i1 + 1 + offset] = fi;
            double wrp = wr;
            wr = wr * wrs - wi * wis;
            wi = wrp * wis + wi * wrs;
        }
    }

    private void pack(double data[], int offset) {
        for (int i = 0; i < n + 2; i++) {
            work[i] = data[i + offset];
        }
        double wrs = cosTheta;
        double wis = sinTheta;
        double wr = 1.0;
        double wi = 0.0;
        for (int i = 0; i < halfN; i++) {
            int i1 = 2 * i;
            double fkr = work[i1];
            double fki = work[i1 + 1];
            int i2 = n - 2 * i;
            double fnr = work[i2];
            double fni = -work[i2 + 1];
            double s1r = fkr + fnr;
            double s1i = fki + fni;
            double s2r = fkr - fnr;
            double s2i = fki - fni;
            double fr = s2r * wr - s2i * wi;
            double fi = s2r * wi + s2i * wr;
            data[i1 + offset] = 0.5 * (s1r - fi);
            data[i1 + 1 + offset] = 0.5 * (s1i + fr);
            double wrp = wr;
            wr = wr * wrs - wi * wis;
            wi = wrp * wis + wi * wrs;
        }
    }
}
