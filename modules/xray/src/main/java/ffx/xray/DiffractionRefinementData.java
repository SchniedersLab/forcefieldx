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
package ffx.xray;

import org.apache.commons.configuration.CompositeConfiguration;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.ReflectionList;
import ffx.numerics.ComplexNumber;
import java.util.Random;

/**
 * <p>DiffractionRefinementData class.</p>
 *
 * @author fennt
 * @version $Id: $
 */
public class DiffractionRefinementData {

    private static final Logger logger = Logger.getLogger(DiffractionRefinementData.class.getName());
    /**
     * number of reflections in the data set
     */
    public final int n;
    /**
     * number of scale parameters
     */
    public final int scale_n;
    /**
     * 2D array of F/sigF data
     */
    public final double fsigf[][];
    /**
     * array of Rfree
     */
    public final int freer[];
    /**
     * calculated atomic structure factors
     */
    public final double fc[][];
    /**
     * calculated bulk solvent structure factors
     */
    public final double fs[][];
    /**
     * scaled sum of Fc and Fs
     */
    public final double fctot[][];
    /**
     * figure of merit and phase
     */
    public final double fomphi[][];
    /**
     * 2mFo - DFc coefficients
     */
    public final double fofc2[][];
    /**
     * mFo - DFc coefficients
     */
    public final double fofc1[][];
    /**
     * derivatives wrt Fc
     */
    public final double dfc[][];
    /**
     * derivatives wrt Fs
     */
    public final double dfs[][];
    /**
     * log likelihoods
     */
    public double llkr, llkf;
    /**
     * reciprocal space reference for structure factor calculations and
     * computing derivatives
     */
    protected CrystalReciprocalSpace crs_fc;
    /**
     * reciprocal space reference for bulk solvent structure factor calculations
     * and computing derivatives
     */
    protected CrystalReciprocalSpace crs_fs;
    /**
     * number of resolution bins
     */
    public final int nbins;
    /**
     * spine scaling coefficients
     */
    public double spline[];
    /**
     * sigmaA coefficient - s
     */
    public double sigmaa[];
    /**
     * sigmaA coefficient - w
     */
    public double sigmaw[];
    /**
     * scaled, E-like Fc
     */
    public double fcesq[];
    /**
     * scaled, E-like Fo
     */
    public double foesq[];
    /**
     * bulk solvent parameters
     */
    public double solvent_a, solvent_b;
    /**
     * bulk solvent scale and Ueq
     */
    public double solvent_k, solvent_ueq;
    /**
     * overall model scale
     */
    public double model_k;
    /**
     * model anisotropic B
     */
    public double model_b[] = new double[6];

    /**
     * duplicated settings - these are also in DiffractionData,
     * but duplicated here until settings are put in their own class
     */
    public int rfreeflag;
    public final double fsigfcutoff;

    /**
     * allocate data given a {@link ReflectionList}
     *
     * @param properties configuration properties
     * @param reflectionlist {@link ReflectionList} to use to allocate data
     */
    public DiffractionRefinementData(CompositeConfiguration properties,
            ReflectionList reflectionlist) {

        int rflag = properties.getInt("rfreeflag", -1);
        fsigfcutoff = properties.getDouble("fsigfcutoff", -1.0);

        this.n = reflectionlist.hkllist.size();
        this.scale_n = reflectionlist.crystal.scale_n;
        this.rfreeflag = rflag;
        this.nbins = reflectionlist.nbins;
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

    /**
     * <p>setCrystalReciprocalSpace_fc</p>
     *
     * @param crs a {@link ffx.xray.CrystalReciprocalSpace} object.
     */
    public void setCrystalReciprocalSpace_fc(CrystalReciprocalSpace crs) {
        this.crs_fc = crs;
    }

    /**
     * <p>setCrystalReciprocalSpace_fs</p>
     *
     * @param crs a {@link ffx.xray.CrystalReciprocalSpace} object.
     */
    public void setCrystalReciprocalSpace_fs(CrystalReciprocalSpace crs) {
        this.crs_fs = crs;
    }

    /**
     * generate 5% of reflections to mark for cross validation/Rfree
     */
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

    /**
     * generate average F/sigF from anomalous F/sigF
     *
     * @param anofsigf an array of double.
     */
    public void generate_fsigf_from_anofsigf(double anofsigf[][]) {
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

    /**
     * set amplitude, F
     *
     * @param i reflection to set
     * @param f value of F desired
     */
    public void set_f(int i, double f) {
        fsigf[i][0] = f;
    }

    /**
     * <p>get_f</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double get_f(int i) {
        return fsigf[i][0];
    }

    /**
     * set amplitude sigma, sigF
     *
     * @param i reflection to set
     * @param sigf value of sigF desired
     */
    public void set_sigf(int i, double sigf) {
        fsigf[i][1] = sigf;
    }

    /**
     * <p>get_sigf</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double get_sigf(int i) {
        return fsigf[i][1];
    }

    /**
     * set amplitude and sigF
     *
     * @param i reflection to set
     * @param f value of F desired
     * @param sigf value of sigF desired
     */
    public void set_fsigf(int i, double f, double sigf) {
        fsigf[i][0] = f;
        fsigf[i][1] = sigf;
    }

    /**
     * <p>get_fsigf</p>
     *
     * @param i a int.
     * @return an array of double.
     */
    public double[] get_fsigf(int i) {
        return fsigf[i];
    }

    /**
     * set FreeR value flag
     *
     * @param i if FreeR value is i, it is marked for cross validation
     */
    public void set_freerflag(int i) {
        rfreeflag = i;
    }

    /**
     * set FreeR value flag of a reflection
     *
     * @param i reflection to set
     * @param f FreeR value to set reflection to
     */
    public void set_freer(int i, int f) {
        freer[i] = f;
    }

    /**
     * <p>get_freer</p>
     *
     * @param i a int.
     * @return a int.
     */
    public int get_freer(int i) {
        return freer[i];
    }

    /**
     * <p>isfreer</p>
     *
     * @param i a int.
     * @param f a int.
     * @return a boolean.
     */
    public boolean isfreer(int i, int f) {
        return (freer[i] == f);
    }

    /**
     * <p>isfreer</p>
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isfreer(int i) {
        return (freer[i] == rfreeflag);
    }

    /**
     * set complex Fc
     *
     * @param i reflection to set
     * @param c {@link ComplexNumber} to set reflection to
     */
    public void set_fc(int i, ComplexNumber c) {
        fc[i][0] = c.re();
        fc[i][1] = c.im();
    }

    /**
     * get the complex number for a Fc reflection
     *
     * @param i reflection to get
     * @return newly allocated {@link ComplexNumber}
     */
    public ComplexNumber get_fc(int i) {
        return new ComplexNumber(fc[i][0], fc[i][1]);
    }

    /**
     * get the complex number for a Fc reflection
     *
     * @param i reflection to get
     * @param c {@link ComplexNumber} to fill
     */
    public void get_fc_ip(int i, ComplexNumber c) {
        c.re(fc[i][0]);
        c.im(fc[i][1]);
    }

    /**
     * get the amplitude of a complex Fc
     *
     * @param i reflection to get
     * @return amplitude of Fc
     */
    public double fc_f(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.abs();
    }

    /**
     * get the phase of a complex Fc
     *
     * @param i reflection to get
     * @return phase of Fc
     */
    public double fc_phi(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.phase();
    }

    /**
     * <p>set_fs</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void set_fs(int i, ComplexNumber c) {
        fs[i][0] = c.re();
        fs[i][1] = c.im();
    }

    /**
     * <p>get_fs</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber get_fs(int i) {
        return new ComplexNumber(fs[i][0], fs[i][1]);
    }

    /**
     * <p>get_fs_ip</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void get_fs_ip(int i, ComplexNumber c) {
        c.re(fs[i][0]);
        c.im(fs[i][1]);
    }

    /**
     * <p>fs_f</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fs_f(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.abs();
    }

    /**
     * <p>fs_phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fs_phi(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.phase();
    }

    /**
     * <p>set_fctot</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void set_fctot(int i, ComplexNumber c) {
        fctot[i][0] = c.re();
        fctot[i][1] = c.im();
    }

    /**
     * <p>get_fctot</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber get_fctot(int i) {
        return new ComplexNumber(fctot[i][0], fctot[i][1]);
    }

    /**
     * <p>get_fctot_ip</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void get_fctot_ip(int i, ComplexNumber c) {
        c.re(fctot[i][0]);
        c.im(fctot[i][1]);
    }

    /**
     * <p>fctot_f</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fctot_f(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.abs();
    }

    /**
     * <p>fctot_phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fctot_phi(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.phase();
    }

    /**
     * <p>set_fofc2</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void set_fofc2(int i, ComplexNumber c) {
        fofc2[i][0] = c.re();
        fofc2[i][1] = c.im();
    }

    /**
     * <p>get_fofc2</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber get_fofc2(int i) {
        return new ComplexNumber(fofc2[i][0], fofc2[i][1]);
    }

    /**
     * <p>get_fofc2_ip</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void get_fofc2_ip(int i, ComplexNumber c) {
        c.re(fofc2[i][0]);
        c.im(fofc2[i][1]);
    }

    /**
     * <p>fofc2_f</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fofc2_f(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.abs();
    }

    /**
     * <p>fofc2_phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fofc2_phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.phase();
    }

    /**
     * <p>set_fofc1</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void set_fofc1(int i, ComplexNumber c) {
        fofc1[i][0] = c.re();
        fofc1[i][1] = c.im();
    }

    /**
     * <p>get_fofc1</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber get_fofc1(int i) {
        return new ComplexNumber(fofc1[i][0], fofc1[i][1]);
    }

    /**
     * <p>get_fofc1_ip</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void get_fofc1_ip(int i, ComplexNumber c) {
        c.re(fofc1[i][0]);
        c.im(fofc1[i][1]);
    }

    /**
     * <p>fofc1_f</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fofc1_f(int i) {
        ComplexNumber c = new ComplexNumber(fofc1[i][0], fofc1[i][1]);

        return c.abs();
    }

    /**
     * <p>fofc1_phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fofc1_phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc1[i][0], fofc1[i][1]);

        return c.phase();
    }

    /**
     * return the current likelihood
     *
     * @return the work likelihood (non-Rfree based)
     */
    public double likelihood_work() {
        return llkr;
    }

    /**
     * return the current likelihood
     *
     * @return the free likelihood (Rfree based)
     */
    public double likelihood_free() {
        return llkf;
    }
}
