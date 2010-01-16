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
package ffx.xray;

import ffx.numerics.Complex;

/**
 *
 * @author fennt
 */
public class RefinementData {

    public final int n;
    public final double fsigf[][];
    public final int freer[];
    public final Complex fc[];
    public final Complex fs[];
    public final Complex fctot[];
    public final double sigmaa[][];
    public final Complex fofc2[];
    public final Complex fofc1[];
    public final Complex fd[];
    // scaling coefficients
    public double solvent_k, solvent_b;
    public double anisok[] = new double[6];

    public RefinementData(int n) {
        this.n = n;
        fsigf = new double[n][2];
        freer = new int[n];
        fc = new Complex[n];
        fs = new Complex[n];
        fctot = new Complex[n];
        sigmaa = new double[n][2];
        fofc2 = new Complex[n];
        fofc1 = new Complex[n];
        fd = new Complex[n];
    }

    public void f(int i, double f) {
        fsigf[i][0] = f;
    }

    public double f(int i) {
        return fsigf[i][0];
    }

    public void sigf(int i, double sigf) {
        fsigf[i][1] = sigf;
    }

    public double sigf(int i) {
        return fsigf[i][1];
    }

    public double[] fsigf(int i) {
        return fsigf[i];
    }

    public void freer(int i, int f) {
        freer[i] = f;
    }

    public int freer(int i) {
        return freer[i];
    }

    public boolean isfreer(int i, int f) {
        return (freer[i] == f);
    }

    public boolean isfreer(int i) {
        return (freer[i] == 1);
    }

    public void fc(int i, Complex c) {
        fc[i] = c;
    }

    public Complex fc(int i) {
        return fc[i];
    }

    public void fs(int i, Complex c) {
        fs[i] = c;
    }

    public Complex fs(int i) {
        return fs[i];
    }

    public void fctot(int i, Complex c) {
        fctot[i] = c;
    }

    public Complex fctot(int i) {
        return fctot[i];
    }

    public double[] sigmaa(int i) {
        return sigmaa[i];
    }

    public void sigmaa(int i, double d[]) {
        sigmaa[i] = d;
    }

    public void sigmasigmaa(int i, double d) {
        sigmaa[i][0] = d;
    }

    public double sigmasigmaa(int i) {
        return sigmaa[i][0];
    }

    public void sigsigmaa(int i, double d) {
        sigmaa[i][1] = d;
    }

    public double sigsigmaa(int i) {
        return sigmaa[i][1];
    }

    public void fofc2(int i, Complex c) {
        fofc2[i] = c;
    }

    public Complex fofc2(int i) {
        return fofc2[i];
    }

    public void fofc1(int i, Complex c) {
        fofc1[i] = c;
    }

    public Complex fofc1(int i) {
        return fofc1[i];
    }

    public void fd(int i, Complex c) {
        fd[i] = c;
    }

    public Complex fd(int i) {
        return fd[i];
    }
}
