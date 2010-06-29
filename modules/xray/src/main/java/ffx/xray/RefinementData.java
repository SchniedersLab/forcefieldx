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
import java.util.Random;

/**
 *
 * @author fennt
 */
public class RefinementData {

    private static final Logger logger = Logger.getLogger(RefinementData.class.getName());
    public final int n;
    public final int scale_n;
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
    // reciprocal space reference
    // for structure factor calculations and computing derivatives
    protected CrystalReciprocalSpace crs_fc;
    protected CrystalReciprocalSpace crs_fs;
    // spline scaling coefficients
    public final int nparams;
    public double spline[];
    public double sigmaa[];
    public double sigmaw[];
    public double fcesq[];
    public double foesq[];
    // bulk solvent parameters
    public boolean gridsearch;
    public double solvent_a, solvent_b;
    // scaling coefficients
    public double solvent_k, solvent_ueq;
    public double model_k;
    public double model_b[] = new double[6];
    // settings
    public final int rfreeflag;
    public final boolean use_3g;
    public final double xrayscaletol;
    public final double sigmaatol;
    public final double bresweight;
    public final double bmass;
    public final boolean residuebfactor;
    public final boolean addanisou;
    public final double xweight;

    public RefinementData(CompositeConfiguration properties,
            ReflectionList reflectionlist) {

        int rflag = properties.getInt("rfreeflag", 1);
        int npar = properties.getInt("nbins", 10);
        gridsearch = properties.getBoolean("gridsearch", false);
        use_3g = properties.getBoolean("use_3g", true);
        xrayscaletol = properties.getDouble("xrayscaletol", 1e-4);
        sigmaatol = properties.getDouble("sigmaatol", 1.0);
        bresweight = properties.getDouble("bresweight", 1.0);
        bmass = properties.getDouble("bmass", 5.0);
        residuebfactor = properties.getBoolean("residuebfactor", false);
        addanisou = properties.getBoolean("addanisou", false);
        xweight = properties.getDouble("xweight", 1.0);

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\nRefinement data settings:\n");
            sb.append("  using cctbx 3 Gaussians: " + use_3g + "\n");
            sb.append("  R Free flag: " + rflag + "\n");
            sb.append("  n bins: " + npar + "\n");
            sb.append("  solvent grid search: " + gridsearch + "\n");
            sb.append("  X-ray scale fit tolerance: " + xrayscaletol + "\n");
            sb.append("  sigma A fit tolerance: " + sigmaatol + "\n");
            sb.append("  B restraint weight: " + bresweight + "\n");
            sb.append("  B Lagrangian mass: " + bmass + "\n");
            sb.append("  B factors refined by residue: " + residuebfactor + "\n");
            sb.append("  add ANISOU for refinement: " + addanisou + "\n");
            sb.append("  X-ray refinement weight: " + xweight + "\n");
            logger.info(sb.toString());
        }

        this.n = reflectionlist.hkllist.size();
        this.scale_n = reflectionlist.crystal.scale_n;
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
            fctot[i][0] = fctot[i][1] = Double.NaN;
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

    public void setCrystalReciprocalSpace_fc(CrystalReciprocalSpace crs) {
        this.crs_fc = crs;
    }

    public void setCrystalReciprocalSpace_fs(CrystalReciprocalSpace crs) {
        this.crs_fs = crs;
    }

    public void generateRFree() {
        Random generator = new Random();
        int nfree = 0;
        for (int i = 0; i < n; i++) {
            if (Double.isNaN(fsigf[i][0])) {
                freer[i] = 0;
                continue;
            }

            int randomi = generator.nextInt(100);
            if (randomi < 5) {
                freer[i] = 1;
                nfree++;
            } else {
                freer[i] = 0;
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\ninternally flagging Rfree reflections\n");
            sb.append("  flagging 5% of observed data reflections\n");
            sb.append(String.format("  selected %d of %d reflections\n",
                    nfree, n));

            logger.info(sb.toString());
        }
    }

    public void set_f(int i, double f) {
        fsigf[i][0] = f;
    }

    public double get_f(int i) {
        return fsigf[i][0];
    }

    public void set_sigf(int i, double sigf) {
        fsigf[i][1] = sigf;
    }

    public double get_sigf(int i) {
        return fsigf[i][1];
    }

    public void set_fsigf(int i, double f, double sigf) {
        fsigf[i][0] = f;
        fsigf[i][1] = sigf;
    }

    public double[] get_fsigf(int i) {
        return fsigf[i];
    }

    public void set_freer(int i, int f) {
        freer[i] = f;
    }

    public int get_freer(int i) {
        return freer[i];
    }

    public boolean isfreer(int i, int f) {
        return (freer[i] == f);
    }

    public boolean isfreer(int i) {
        return (freer[i] == rfreeflag);
    }

    public void set_fc(int i, ComplexNumber c) {
        fc[i][0] = c.re();
        fc[i][1] = c.im();
    }

    public ComplexNumber get_fc(int i) {
        return new ComplexNumber(fc[i][0], fc[i][1]);
    }

    public void get_fc_ip(int i, ComplexNumber c) {
        c.re(fc[i][0]);
        c.im(fc[i][1]);
    }

    public double fc_f(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.abs();
    }

    public double fc_phi(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.phase();
    }

    public void set_fs(int i, ComplexNumber c) {
        fs[i][0] = c.re();
        fs[i][1] = c.im();
    }

    public ComplexNumber get_fs(int i) {
        return new ComplexNumber(fs[i][0], fs[i][1]);
    }

    public void get_fs_ip(int i, ComplexNumber c) {
        c.re(fs[i][0]);
        c.im(fs[i][1]);
    }

    public double fs_f(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.abs();
    }

    public double fs_phi(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.phase();
    }

    public void set_fctot(int i, ComplexNumber c) {
        fctot[i][0] = c.re();
        fctot[i][1] = c.im();
    }

    public ComplexNumber get_fctot(int i) {
        return new ComplexNumber(fctot[i][0], fctot[i][1]);
    }

    public void get_fctot_ip(int i, ComplexNumber c) {
        c.re(fctot[i][0]);
        c.im(fctot[i][1]);
    }

    public double fctot_f(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.abs();
    }

    public double fctot_phi(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.phase();
    }

    public void set_fofc2(int i, ComplexNumber c) {
        fofc2[i][0] = c.re();
        fofc2[i][1] = c.im();
    }

    public ComplexNumber get_fofc2(int i) {
        return new ComplexNumber(fofc2[i][0], fofc2[i][1]);
    }

    public void get_fofc2_ip(int i, ComplexNumber c) {
        c.re(fofc2[i][0]);
        c.im(fofc2[i][1]);
    }

    public double fofc2_f(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.abs();
    }

    public double fofc2_phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.phase();
    }

    public void set_fofc1(int i, ComplexNumber c) {
        fofc1[i][0] = c.re();
        fofc1[i][1] = c.im();
    }

    public ComplexNumber get_fofc1(int i) {
        return new ComplexNumber(fofc1[i][0], fofc1[i][1]);
    }

    public void get_fofc1_ip(int i, ComplexNumber c) {
        c.re(fofc1[i][0]);
        c.im(fofc1[i][1]);
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
