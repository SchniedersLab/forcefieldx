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
    public final double anofsigf[][];
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
    public final int nbins;
    public double spline[];
    public double sigmaa[];
    public double sigmaw[];
    public double fcesq[];
    public double foesq[];
    // bulk solvent parameters
    public boolean gridsearch;
    public double solvent_a, solvent_b;
    // scaling coefficients
    public boolean splinefit;
    public double solvent_k, solvent_ueq;
    public double model_k;
    public double model_b[] = new double[6];
    // settings
    public final double fsigfcutoff;
    public int rfreeflag;
    public final boolean use_3g;
    public final double xrayscaletol;
    public final double sigmaatol;
    public final double xweight;
    public final double bsimweight;
    public final double bnonzeroweight;
    public final double bmass;
    public final boolean residuebfactor;
    public final int nresiduebfactor;
    public final boolean addanisou;
    public final boolean refinemolocc;
    public final double occmass;

    public RefinementData(CompositeConfiguration properties,
            ReflectionList reflectionlist) {

        int rflag = properties.getInt("rfreeflag", -1);
        fsigfcutoff = properties.getDouble("fsigfcutoff", -1.0);
        gridsearch = properties.getBoolean("gridsearch", false);
        splinefit = properties.getBoolean("splinefit", true);
        use_3g = properties.getBoolean("use_3g", true);
        xrayscaletol = properties.getDouble("xrayscaletol", 1e-4);
        sigmaatol = properties.getDouble("sigmaatol", 1.0);
        xweight = properties.getDouble("xweight", 1.0);
        bsimweight = properties.getDouble("bsimweight", 1.0);
        bnonzeroweight = properties.getDouble("bnonzeroweight", 1.0);
        bmass = properties.getDouble("bmass", 5.0);
        residuebfactor = properties.getBoolean("residuebfactor", false);
        nresiduebfactor = properties.getInt("nresiduebfactor", 1);
        addanisou = properties.getBoolean("addanisou", false);
        refinemolocc = properties.getBoolean("refinemolocc", false);
        occmass = properties.getDouble("occmass", 10.0);

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\nRefinement data settings:\n");
            sb.append("  using cctbx 3 Gaussians (use_3g): " + use_3g + "\n");
            sb.append("  resolution dependent spline scale (splinefit): " + splinefit + "\n");
            sb.append("  F/sigF cutoff (fsigfcutoff): " + fsigfcutoff + "\n");
            sb.append("  R Free flag (rfreeflag) (if -1, value will be updated when data is read in): " + rflag + "\n");
            sb.append("  n bins (nbins): " + reflectionlist.nbins + "\n");
            sb.append("  solvent grid search (gridsearch): " + gridsearch + "\n");
            sb.append("  X-ray scale fit tolerance (xrayscaletol): " + xrayscaletol + "\n");
            sb.append("  sigma A fit tolerance (sigmaatol): " + sigmaatol + "\n");
            sb.append("  X-ray refinement weight (xweight): " + xweight + "\n");
            sb.append("  B similarity weight (bsimweight): " + bsimweight + "\n");
            sb.append("  B non-zero weight (bnonzeroweight): " + bnonzeroweight + "\n");
            sb.append("  B Lagrangian mass (bmass): " + bmass + "\n");
            sb.append("  B factors refined by residue (residuebfactor): " + residuebfactor + "\n");
            sb.append("    (if true, num. residues per B (nresiduebfactor): " + nresiduebfactor + ")\n");
            sb.append("  add ANISOU for refinement (addanisou): " + addanisou + "\n");
            sb.append("  refine occupancies on molecules (HETATMs - refinemolocc): " + refinemolocc + "\n");
            sb.append("  Occupancy Lagrangian mass (occmass): " + occmass + "\n");
            logger.info(sb.toString());
        }

        this.n = reflectionlist.hkllist.size();
        this.scale_n = reflectionlist.crystal.scale_n;
        this.rfreeflag = rflag;
        this.nbins = reflectionlist.nbins;
        fsigf = new double[n][2];
        anofsigf = new double[n][4];
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
            anofsigf[i][0] = anofsigf[i][1] = anofsigf[i][2] = anofsigf[i][3] = Double.NaN;
            fctot[i][0] = fctot[i][1] = Double.NaN;
        }

        spline = new double[nbins * 2];
        sigmaa = new double[nbins];
        sigmaw = new double[nbins];
        fcesq = new double[nbins];
        foesq = new double[nbins];
        for (int i = 0; i < nbins; i++) {
            spline[i] = spline[i + nbins] = sigmaa[i] = fcesq[i] = foesq[i] = 1.0;
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
        int free;
        int nonfree;
        int nfree = 0;

        if (rfreeflag == 0) {
            free = 0;
            nonfree = 1;
        } else {
            free = 1;
            nonfree = 0;
        }

        for (int i = 0; i < n; i++) {
            if (Double.isNaN(fsigf[i][0])) {
                freer[i] = nonfree;
                continue;
            }

            int randomi = generator.nextInt(100);
            if (randomi < 5) {
                freer[i] = free;
                nfree++;
            } else {
                freer[i] = nonfree;
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

    public void generate_fsigf_from_anofsigf() {
        for (int i = 0; i < n; i++) {
            if (Double.isNaN(anofsigf[i][0])
                    && Double.isNaN(anofsigf[i][2])) {
            } else if (Double.isNaN(anofsigf[i][0])) {
                fsigf[i][0] = anofsigf[i][2];
                fsigf[i][1] = anofsigf[i][3];
            } else if (Double.isNaN(anofsigf[i][2])) {
                fsigf[i][0] = anofsigf[i][0];
                fsigf[i][1] = anofsigf[i][1];
            } else {
                fsigf[i][0] = (anofsigf[i][0] + anofsigf[i][2]) / 2.0;
                fsigf[i][1] = Math.sqrt(anofsigf[i][1] * anofsigf[i][1]
                        + anofsigf[i][3] * anofsigf[i][3]);
            }
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

    public void set_ano_fplus(int i, double f) {
        anofsigf[i][0] = f;
    }

    public void set_ano_fminus(int i, double f) {
        anofsigf[i][2] = f;
    }

    public double get_ano_fplus(int i) {
        return anofsigf[i][0];
    }

    public double get_ano_fminus(int i) {
        return anofsigf[i][2];
    }

    public void set_ano_sigfplus(int i, double f) {
        anofsigf[i][1] = f;
    }

    public void set_ano_sigfminus(int i, double f) {
        anofsigf[i][3] = f;
    }

    public double get_ano_sigfplus(int i) {
        return anofsigf[i][1];
    }

    public double get_ano_sigfminus(int i) {
        return anofsigf[i][3];
    }

    public void set_ano_fsigfplus(int i, double f, double sigf) {
        anofsigf[i][0] = f;
        anofsigf[i][1] = sigf;
    }

    public void set_ano_fsigfminus(int i, double f, double sigf) {
        anofsigf[i][2] = f;
        anofsigf[i][3] = sigf;
    }

    public void set_freerflag(int i) {
        rfreeflag = i;
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
