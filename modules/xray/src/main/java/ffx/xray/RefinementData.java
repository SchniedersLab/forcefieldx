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

import org.apache.commons.configuration.CompositeConfiguration;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.ReflectionList;
import ffx.numerics.ComplexNumber;

/**
 *
 * @author fennt
 */
public class RefinementData {

    private static final Logger logger = Logger.getLogger(CrystalStats.class.getName());
    public final int n;
    public final int scale_n;
    public final int solvent_n;
    // data
    public final double fsigf[][];
    public final int freer[];
    // calculated atomic structure factors
    public final double fc[][];
    // calculted bulk solvent structure factors
    public final double fs[][];
    // scaled sum of Fc and Fs
    public final double fctot[][];
    // figure of merit and phase
    public final double fomphi[][];
    // 2mFo - DFc coefficients
    public final double fofc2[][];
    // mFo - DFc coefficients
    public final double fofc1[][];
    // derivatives wrt Fc
    public final double dfc[][];
    // derivatives wrt Fs
    public final double dfs[][];
    // log likelihoods
    public double llkr, llkf;
    // spline scaling coefficients
    public final int nparams;
    public double spline[];
    public double sigmaa[];
    public double sigmaw[];
    public double fcesq[];
    public double foesq[];
    // bulk solvent parameters
    public double solvent_a, solvent_sd;
    // scaling coefficients
    public double solvent_k, solvent_ueq;
    public double model_k;
    public double model_b[] = new double[6];
    // settings
    public final int rfreeflag;

    public RefinementData(CompositeConfiguration properties,
            ReflectionList reflectionlist) {

        int rflag = properties.getInt("rfreeflag", 1);
        int npar = properties.getInt("nbins", 10);
        boolean bulksolvent = properties.getBoolean("bulksolvent", true);

        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            sb.append("Refinement data settings:");
            sb.append("  R Free flag: " + rflag);
            sb.append("  n bins: " + npar);
            sb.append("  bulk solvent: " + bulksolvent);
            logger.info(sb.toString());
        }

        this.n = reflectionlist.hkllist.size();
        this.scale_n = reflectionlist.crystal.scale_n;
        this.solvent_a = 11.5;
        this.solvent_sd = 0.75;
        this.solvent_n = bulksolvent ? 3 : 1;
        this.nparams = npar;
        this.rfreeflag = rflag;
        fsigf = new double[n][2];
        freer = new int[n];
        fc = new double[n][2];
        fs = new double[n][2];
        fctot = new double[n][2];
        fomphi = new double[n][2];
        fofc2 = new double[n][2];
        fofc1 = new double[n][2];
        dfc = new double[n][2];
        dfs = new double[n][2];

        for (int i = 0; i < n; i++) {
            fsigf[i][0] = fsigf[i][1] = Double.NaN;
        }

        spline = new double[nparams * 2];
        sigmaa = new double[nparams];
        sigmaw = new double[nparams];
        fcesq = new double[nparams];
        foesq = new double[nparams];
        for (int i = 0; i < nparams; i++) {
            spline[i] = spline[i + nparams] = sigmaa[i] = fcesq[i] = foesq[i] = 1.0;
            sigmaw[i] = 0.05;
        }

        // initial guess for scaling/bulk solvent correction
        solvent_k = 0.33;
        solvent_ueq = 50.0 / (8.0 * Math.PI * Math.PI);
        model_k = 0.0;
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
        return (freer[i] == rfreeflag);
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

    public double fs_f(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.abs();
    }

    public double fs_phi(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.phase();
    }

    public void fctot(int i, ComplexNumber c) {
        fctot[i][0] = c.re();
        fctot[i][1] = c.im();
    }

    public ComplexNumber fctot(int i) {
        return new ComplexNumber(fctot[i][0], fctot[i][1]);
    }

    public double fctot_f(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.abs();
    }

    public double fctot_phi(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.phase();
    }

    public void fofc2(int i, ComplexNumber c) {
        fofc2[i][0] = c.re();
        fofc2[i][1] = c.im();
    }

    public ComplexNumber fofc2(int i) {
        return new ComplexNumber(fofc2[i][0], fofc2[i][1]);
    }

    public double fofc2_f(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.abs();
    }

    public double fofc2_phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.phase();
    }

    public void fofc1(int i, ComplexNumber c) {
        fofc1[i][0] = c.re();
        fofc1[i][1] = c.im();
    }

    public ComplexNumber fofc1(int i) {
        return new ComplexNumber(fofc1[i][0], fofc1[i][1]);
    }

    public double fofc1_f(int i) {
        ComplexNumber c = new ComplexNumber(fofc1[i][0], fofc1[i][1]);

        return c.abs();
    }

    public double fofc1_phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc1[i][0], fofc1[i][1]);

        return c.phase();
    }

    public double likelihood_work() {
        return llkr;
    }

    public double likelihood_free() {
        return llkf;
    }
}
