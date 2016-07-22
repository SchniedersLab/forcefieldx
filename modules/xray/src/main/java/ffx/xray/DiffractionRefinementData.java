/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.xray;

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Double.isNaN;

import org.apache.commons.configuration.CompositeConfiguration;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.ReflectionList;
import ffx.numerics.ComplexNumber;

/**
 * <p>
 * DiffractionRefinementData class.</p>
 *
 * @author Timothy D. Fenn
 *
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
     * duplicated settings - these are also in DiffractionData, but duplicated
     * here until settings are put in their own class
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

        int rflag = -1;
        if (properties != null) {
            rflag = properties.getInt("rfreeflag", -1);
            fsigfcutoff = properties.getDouble("fsigfcutoff", -1.0);
        } else {
            fsigfcutoff = -1.0;
        }

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
        solvent_ueq = 50.0 / (8.0 * PI * PI);
        model_k = 0.0;
    }

    /**
     * <p>
     * setCrystalReciprocalSpace_fc</p>
     *
     * @param crs a {@link ffx.xray.CrystalReciprocalSpace} object.
     */
    public void setCrystalReciprocalSpace_fc(CrystalReciprocalSpace crs) {
        this.crs_fc = crs;
    }

    /**
     * <p>
     * setCrystalReciprocalSpace_fs</p>
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
            if (isNaN(fsigf[i][0])) {
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
            sb.append("\n Internally flagging Rfree reflections:\n");
            sb.append("  Flagging 5% of observed data reflections\n");
            sb.append(String.format("  Selected %d of %d reflections\n", nfree, n));
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
            if (isNaN(anofsigf[i][0]) && isNaN(anofsigf[i][2])) {
            } else if (isNaN(anofsigf[i][0])) {
                fsigf[i][0] = anofsigf[i][2];
                fsigf[i][1] = anofsigf[i][3];
            } else if (isNaN(anofsigf[i][2])) {
                fsigf[i][0] = anofsigf[i][0];
                fsigf[i][1] = anofsigf[i][1];
            } else {
                fsigf[i][0] = (anofsigf[i][0] + anofsigf[i][2]) / 2.0;
                fsigf[i][1] = sqrt(anofsigf[i][1] * anofsigf[i][1]
                        + anofsigf[i][3] * anofsigf[i][3]);
            }
        }
    }

    /**
     * generate amplitudes from intensities. Does NOT use French and Wilson
     * scaling, just simple square root.
     */
    public void intensities_to_amplitudes() {
        double tmp;

        for (int i = 0; i < n; i++) {
            if (fsigf[i][0] > 0.0) {
                tmp = fsigf[i][0];
                fsigf[i][0] = sqrt(tmp);
                if (fsigf[i][1] < tmp) {
                    fsigf[i][1] = fsigf[i][0]
                            - sqrt(tmp - fsigf[i][1]);
                } else {
                    fsigf[i][1] = fsigf[i][0];
                }
            } else if (!isNaN(fsigf[i][0])) {
                fsigf[i][0] = 0.0;
                fsigf[i][1] = 0.0;
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\n Internally converting intensities to amplitudes;\n");
            sb.append(" note this does not use French & Wilson scaling.\n");
            logger.info(sb.toString());
        }
    }

    /**
     * set amplitude, F
     *
     * @param i reflection to set
     * @param f value of F desired
     */
    public void setF(int i, double f) {
        fsigf[i][0] = f;
    }

    /**
     * <p>
     * getF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getF(int i) {
        return fsigf[i][0];
    }

    /**
     * set amplitude sigma, sigF
     *
     * @param i reflection to set
     * @param sigf value of sigF desired
     */
    public void setSigF(int i, double sigf) {
        fsigf[i][1] = sigf;
    }

    /**
     * <p>
     * getSigF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double getSigF(int i) {
        return fsigf[i][1];
    }

    /**
     * set amplitude and sigF
     *
     * @param i reflection to set
     * @param f value of F desired
     * @param sigf value of sigF desired
     */
    public void setFSigF(int i, double f, double sigf) {
        fsigf[i][0] = f;
        fsigf[i][1] = sigf;
    }

    /**
     * <p>
     * getFSigF</p>
     *
     * @param i a int.
     * @return an array of double.
     */
    public double[] getFSigF(int i) {
        return fsigf[i];
    }

    /**
     * set FreeR value flag
     *
     * @param i if FreeR value is i, it is marked for cross validation
     */
    public void setFreeRFlag(int i) {
        rfreeflag = i;
    }

    /**
     * set FreeR value flag of a reflection
     *
     * @param i reflection to set
     * @param f FreeR value to set reflection to
     */
    public void setFreeR(int i, int f) {
        freer[i] = f;
    }

    /**
     * <p>
     * getFreeR</p>
     *
     * @param i a int.
     * @return a int.
     */
    public int getFreeR(int i) {
        return freer[i];
    }

    /**
     * <p>
     * isFreeR</p>
     *
     * @param i a int.
     * @param f a int.
     * @return a boolean.
     */
    public boolean isFreeR(int i, int f) {
        return (freer[i] == f);
    }

    /**
     * <p>
     * isFreeR</p>
     *
     * @param i a int.
     * @return a boolean.
     */
    public boolean isFreeR(int i) {
        return (freer[i] == rfreeflag);
    }

    /**
     * set complex Fc
     *
     * @param i reflection to set
     * @param c {@link ComplexNumber} to set reflection to
     */
    public void setFc(int i, ComplexNumber c) {
        fc[i][0] = c.re();
        fc[i][1] = c.im();
    }

    /**
     * get the complex number for a Fc reflection
     *
     * @param i reflection to get
     * @return newly allocated {@link ComplexNumber}
     */
    public ComplexNumber getFc(int i) {
        return new ComplexNumber(fc[i][0], fc[i][1]);
    }

    /**
     * get the complex number for a Fc reflection
     *
     * @param i reflection to get
     * @param c {@link ComplexNumber} to fill
     */
    public void getFcIP(int i, ComplexNumber c) {
        c.re(fc[i][0]);
        c.im(fc[i][1]);
    }

    /**
     * get the amplitude of a complex Fc
     *
     * @param i reflection to get
     * @return amplitude of Fc
     */
    public double fcF(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.abs();
    }

    /**
     * get the phase of a complex Fc
     *
     * @param i reflection to get
     * @return phase of Fc
     */
    public double fcPhi(int i) {
        ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);

        return c.phase();
    }

    /**
     * <p>
     * setFs</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void setFs(int i, ComplexNumber c) {
        fs[i][0] = c.re();
        fs[i][1] = c.im();
    }

    /**
     * <p>
     * getFs</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber getFs(int i) {
        return new ComplexNumber(fs[i][0], fs[i][1]);
    }

    /**
     * <p>
     * getFsIP</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void getFsIP(int i, ComplexNumber c) {
        c.re(fs[i][0]);
        c.im(fs[i][1]);
    }

    /**
     * <p>
     * fsF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fsF(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.abs();
    }

    /**
     * <p>
     * fsPhi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fsPhi(int i) {
        ComplexNumber c = new ComplexNumber(fs[i][0], fs[i][1]);

        return c.phase();
    }

    /**
     * <p>
     * setFcTot</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void setFcTot(int i, ComplexNumber c) {
        fctot[i][0] = c.re();
        fctot[i][1] = c.im();
    }

    /**
     * <p>
     * getFcTot</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber getFcTot(int i) {
        return new ComplexNumber(fctot[i][0], fctot[i][1]);
    }

    /**
     * <p>
     * getFcTotIP</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void getFcTotIP(int i, ComplexNumber c) {
        c.re(fctot[i][0]);
        c.im(fctot[i][1]);
    }

    /**
     * <p>
     * fcTotF</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fcTotF(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.abs();
    }

    /**
     * <p>
     * fcTotPhi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double fcTotPhi(int i) {
        ComplexNumber c = new ComplexNumber(fctot[i][0], fctot[i][1]);

        return c.phase();
    }

    /**
     * <p>
     * setFoFc2</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void setFoFc2(int i, ComplexNumber c) {
        fofc2[i][0] = c.re();
        fofc2[i][1] = c.im();
    }

    /**
     * <p>
     * getFoFc2</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber getFoFc2(int i) {
        return new ComplexNumber(fofc2[i][0], fofc2[i][1]);
    }

    /**
     * <p>
     * getFoFc2IP</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void getFoFc2IP(int i, ComplexNumber c) {
        c.re(fofc2[i][0]);
        c.im(fofc2[i][1]);
    }

    /**
     * <p>
     * FoFc2F</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double FoFc2F(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.abs();
    }

    /**
     * <p>
     * FoFc2Phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double FoFc2Phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc2[i][0], fofc2[i][1]);

        return c.phase();
    }

    /**
     * <p>
     * setFoFc1</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void setFoFc1(int i, ComplexNumber c) {
        fofc1[i][0] = c.re();
        fofc1[i][1] = c.im();
    }

    /**
     * <p>
     * getFoFc1</p>
     *
     * @param i a int.
     * @return a {@link ffx.numerics.ComplexNumber} object.
     */
    public ComplexNumber getFoFc1(int i) {
        return new ComplexNumber(fofc1[i][0], fofc1[i][1]);
    }

    /**
     * <p>
     * getFoFc1IP</p>
     *
     * @param i a int.
     * @param c a {@link ffx.numerics.ComplexNumber} object.
     */
    public void getFoFc1IP(int i, ComplexNumber c) {
        c.re(fofc1[i][0]);
        c.im(fofc1[i][1]);
    }

    /**
     * <p>
     * foFc1F</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double foFc1F(int i) {
        ComplexNumber c = new ComplexNumber(fofc1[i][0], fofc1[i][1]);

        return c.abs();
    }

    /**
     * <p>
     * foFc1Phi</p>
     *
     * @param i a int.
     * @return a double.
     */
    public double foFc1Phi(int i) {
        ComplexNumber c = new ComplexNumber(fofc1[i][0], fofc1[i][1]);

        return c.phase();
    }

    /**
     * return the current likelihood
     *
     * @return the work likelihood (non-Rfree based)
     */
    public double likelihoodWork() {
        return llkr;
    }

    /**
     * return the current likelihood
     *
     * @return the free likelihood (Rfree based)
     */
    public double likelihoodFree() {
        return llkf;
    }
}
