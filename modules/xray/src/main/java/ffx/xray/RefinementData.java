/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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

import ffx.numerics.ComplexNumber;

/**
 *
 * @author fennt
 */
public class RefinementData {

    public final int n;
    public final double fsigf[][];
    public final int freer[];
    public final double fc[][];
    public final double fs[][];
    public final double fctot[][];
    public final double fofc2[][];
    public final double fofc1[][];
    public final double dfc[][];
    public final double dfs[][];
    // spline scaling coefficients
    public final int nparams;
    public double spline[];
    public double sigmaa[];
    public double sigmaw[];
    public double fcesq[];
    public double foesq[];
    // scaling coefficients
    public double solvent_k, solvent_ueq;
    public double model_k;
    public double aniso_b[] = new double[6];

    public RefinementData(int n) {
        this(n, 10);
    }

    public RefinementData(int n, int nparams) {
        this.n = n;
        this.nparams = nparams;
        fsigf = new double[n][2];
        freer = new int[n];
        fc = new double[n][2];
        fs = new double[n][2];
        fctot = new double[n][2];
        fofc2 = new double[n][2];
        fofc1 = new double[n][2];
        dfc = new double[n][2];
        dfs = new double[n][2];

        for (int i = 0; i < n; i++) {
            fsigf[i][0] = fsigf[i][1] = Double.NaN;
        }

        spline = new double[nparams];
        sigmaa = new double[nparams];
        sigmaw = new double[nparams];
        fcesq = new double[nparams];
        foesq = new double[nparams];
        for (int i = 0; i < nparams; i++) {
            spline[i] = sigmaa[i] = fcesq[i] = foesq[i] = 1.0;
        }

        // initial guess for scaling/bulk solvent correction
        solvent_k = 0.33;
        solvent_ueq = 50.0 / (8.0 * Math.PI * Math.PI);
        model_k = 1.0;
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

    public void fsigf(int i, double f, double sigf) {
        fsigf[i][0] = f;
        fsigf[i][1] = sigf;
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

    public void fc(int i, ComplexNumber c) {
        fc[i][0] = c.re();
        fc[i][1] = c.im();
    }

    public ComplexNumber fc(int i) {
        return new ComplexNumber(fc[i][0], fc[i][1]);
    }

    public double fc_f(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.abs();
    }

    public double fc_phi(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.phase();
    }

    public void fs(int i, ComplexNumber c) {
        fs[i][0] = c.re();
        fs[i][1] = c.im();
    }

    public ComplexNumber fs(int i) {
        return new ComplexNumber(fs[i][0], fs[i][1]);
    }

    public void fctot(int i, ComplexNumber c) {
        fctot[i][0] = c.re();
        fctot[i][1] = c.im();
    }

    public ComplexNumber fctot(int i) {
        return new ComplexNumber(fctot[i][0], fctot[i][1]);
    }

    public void fofc2(int i, ComplexNumber c) {
        fofc2[i][0] = c.re();
        fofc2[i][1] = c.im();
    }

    public ComplexNumber fofc2(int i) {
        return new ComplexNumber(fofc2[i][0], fofc2[i][1]);
    }

    public void fofc1(int i, ComplexNumber c) {
        fofc1[i][0] = c.re();
        fofc1[i][1] = c.im();
    }

    public ComplexNumber fofc1(int i) {
        return new ComplexNumber(fofc1[i][0], fofc1[i][1]);
    }
}
