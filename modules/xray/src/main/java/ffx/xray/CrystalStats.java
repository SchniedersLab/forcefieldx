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
package ffx.xray;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.ComplexNumber;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

/**
 * Crystal statistics output/logger
 *
 * @author Timothy D. Fenn
 *
 */
public class CrystalStats {

    private static final Logger logger = Logger.getLogger(CrystalStats.class.getName());
    private final ReflectionList reflectionlist;
    private final DiffractionRefinementData refinementdata;
    private final Crystal crystal;
    private final ReflectionSpline spline;
    private final int n;
    private final double fo[][];
    private final int freer[];
    private final double fc[][];
    private final double fomphi[][];
    protected int nobshkl, highnobshkl, nobsrfree, highnobsrfree;
    protected double reshigh, reslow, highreshigh, highreslow;
    protected double completeness, highcompleteness;
    protected double rall, r, rfree, highr, highrfree;
    protected double blowdpi, blowdpih, cruickdpi, cruickdpih;
    private boolean print;

    /**
     * constructor
     *
     * @param reflectionlist {@link ReflectionList} to use for logging
     * @param refinementdata {@link DiffractionRefinementData} to use for
     * logging
     */
    public CrystalStats(ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.n = refinementdata.nbins;
        this.fo = refinementdata.fsigf;
        this.freer = refinementdata.freer;
        this.fc = refinementdata.fctot;
        this.fomphi = refinementdata.fomphi;

        this.spline = new ReflectionSpline(reflectionlist,
                refinementdata.spline.length);

        blowdpi = -1.0;
        this.print = true;
    }

    /**
     * <p>
     * getPDBHeaderString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getPDBHeaderString() {
        print = false;
        printHKLStats();
        printRStats();
        print = true;

        StringBuilder sb = new StringBuilder();

        sb.append("REMARK   3  DATA USED IN REFINEMENT\n");
        sb.append(String.format("REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : %6.2f\n", reshigh));
        sb.append(String.format("REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : %6.2f\n", reslow));
        if (refinementdata.fsigfcutoff > 0.0) {
            sb.append(String.format("REMARK   3   DATA CUTOFF            (SIGMA(F)) : %6.2f\n", refinementdata.fsigfcutoff));
        } else {
            sb.append("REMARK   3   DATA CUTOFF            (SIGMA(F)) : NONE\n");
        }
        sb.append(String.format("REMARK   3   COMPLETENESS FOR RANGE        (%%) : %6.2f\n", completeness));
        sb.append("REMARK   3   NUMBER OF REFLECTIONS             : " + nobshkl + "\n");
        sb.append("REMARK   3\n");
        sb.append("REMARK   3  FIT TO DATA USED IN REFINEMENT\n");
        sb.append("REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT\n");
        sb.append("REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM\n");
        sb.append(String.format("REMARK   3   R VALUE               (OBSERVED) : %8.6f\n", rall / 100.0));
        sb.append(String.format("REMARK   3   R VALUE            (WORKING SET) : %8.6f\n", r / 100.0));
        sb.append(String.format("REMARK   3   FREE R VALUE                     : %8.6f\n", rfree / 100.0));
        sb.append(String.format("REMARK   3   FREE R VALUE TEST SET SIZE   (%%) : %6.1f\n", (((double) nobsrfree) / nobshkl) * 100.0));
        sb.append("REMARK   3   FREE R VALUE TEST SET COUNT      : " + nobsrfree + "\n");
        sb.append("REMARK   3\n");
        sb.append("REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN\n");
        sb.append("REMARK   3   TOTAL NUMBER OF BINS USED           : " + n + "\n");
        sb.append(String.format("REMARK   3   BIN RESOLUTION RANGE HIGH           : %6.2f\n", highreshigh));
        sb.append(String.format("REMARK   3   BIN RESOLUTION RANGE LOW            : %6.2f\n", highreslow));
        sb.append("REMARK   3   REFLECTION IN BIN     (WORKING SET) : " + highnobshkl + "\n");
        sb.append(String.format("REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%%) : %6.2f\n", highcompleteness));
        sb.append(String.format("REMARK   3   BIN R VALUE           (WORKING SET) : %8.6f\n", highr / 100.0));
        sb.append("REMARK   3   BIN FREE R VALUE SET COUNT          : " + highnobsrfree + "\n");
        sb.append(String.format("REMARK   3   BIN FREE R VALUE                    : %8.6f\n", highrfree / 100.0));
        sb.append("REMARK   3\n");

        sb.append("REMARK   3  OVERALL SCALE FACTORS\n");
        sb.append(String.format("REMARK   3   SCALE: %4.2f\n",
                exp(0.25 * refinementdata.model_k)));
        sb.append(String.format("REMARK   3   ANISOTROPIC SCALE TENSOR:\n"));
        sb.append(String.format("REMARK   3    %g %g %g\n",
                refinementdata.model_b[0],
                refinementdata.model_b[3],
                refinementdata.model_b[4]));
        sb.append(String.format("REMARK   3    %g %g %g\n",
                refinementdata.model_b[3],
                refinementdata.model_b[1],
                refinementdata.model_b[5]));
        sb.append(String.format("REMARK   3    %g %g %g\n",
                refinementdata.model_b[4],
                refinementdata.model_b[5],
                refinementdata.model_b[2]));
        sb.append("REMARK   3\n");

        if (refinementdata.crs_fs.solventModel != SolventModel.NONE) {
            sb.append("REMARK   3  BULK SOLVENT MODELLING\n");
            switch (refinementdata.crs_fs.solventModel) {
                case BINARY:
                    sb.append("REMARK   3   METHOD USED: BINARY MASK\n");
                    sb.append(String.format("REMARK   3    PROBE RADIUS  : %g\n",
                            refinementdata.solvent_a));
                    sb.append(String.format("REMARK   3    SHRINK RADIUS : %g\n",
                            refinementdata.solvent_b));
                    break;
                case POLYNOMIAL:
                    sb.append("REMARK   3   METHOD USED: POLYNOMIAL SWITCH\n");
                    sb.append(String.format("REMARK   3    ATOMIC RADIUS BUFFER : %g\n",
                            refinementdata.solvent_a));
                    sb.append(String.format("REMARK   3    SWITCH RADIUS        : %g\n",
                            refinementdata.solvent_b));
                    break;
                case GAUSSIAN:
                    sb.append("REMARK   3   METHOD USED: GAUSSIAN\n");
                    sb.append(String.format("REMARK   3    ATOMIC RADIUS BUFFER : %g\n",
                            refinementdata.solvent_a));
                    sb.append(String.format("REMARK   3    STD DEV SCALE        : %g\n",
                            refinementdata.solvent_b));
                    break;
            }
            sb.append(String.format("REMARK   3    K_SOL: %g\n",
                    refinementdata.solvent_k));
            sb.append(String.format("REMARK   3    B_SOL: %g\n",
                    refinementdata.solvent_ueq * 8.0 * Math.PI * Math.PI));
            sb.append("REMARK   3\n");
        }

        if (blowdpi > 0.0) {
            sb.append("REMARK   3  ERROR ESTIMATES\n");
            sb.append("REMARK   3   ACTA CRYST (1999) D55, 583-601\n");
            sb.append("REMARK   3   ACTA CRYST (2002) D58, 792-797\n");
            sb.append(String.format("REMARK   3   BLOW DPI ALL ATOMS (EQN 7)          : %7.4f\n", blowdpih));
            sb.append(String.format("REMARK   3   BLOW DPI NONH ATOMS (EQN 7)         : %7.4f\n", blowdpi));
            sb.append(String.format("REMARK   3   CRUICKSHANK DPI ALL ATOMS (EQN 27)  : %7.4f\n", cruickdpih));
            sb.append(String.format("REMARK   3   CRUICKSHANK DPI NONH ATOMS (EQN 27) : %7.4f\n", cruickdpi));
            sb.append("REMARK   3\n");
        }

        sb.append("REMARK   3  DATA TARGET\n");
        sb.append("REMARK   3   METHOD USED: MAXIMUM LIKELIHOOD\n");
        sb.append(String.format("REMARK   3    -LOG LIKELIHOOD            : %g\n",
                refinementdata.llkr));
        sb.append(String.format("REMARK   3    -LOG LIKELIHOOD (FREE SET) : %g\n",
                refinementdata.llkf));
        sb.append("REMARK   3\n");

        return sb.toString();
    }

    /**
     * simply return the current R value
     *
     * @return r value as a percent
     */
    public double getR() {
        double numer = 0.0;
        double denom = 0.0;
        double sum = 0.0;
        double sumfo = 0.0;
        double sumall = 0.0;
        double sumfoall = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            //int b = ih.bin();

            // ignored cases
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // spline setup
            double ss = Crystal.invressq(crystal, ih);
            double fh = spline.f(ss, refinementdata.spline);

            ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
            numer = abs(abs(fo[i][0]) - fh * abs(c.abs()));
            denom = abs(fo[i][0]);
            sumall += numer;
            sumfoall += denom;
            if (!refinementdata.isFreeR(i)) {
                sum += numer;
                sumfo += denom;
            }
        }

        rall = (sumall / sumfoall) * 100.0;
        r = (sum / sumfo) * 100.0;
        return r;
    }

    /**
     * simply return the current Rfree value
     *
     * @return rfree value as a percent
     */
    public double getRFree() {
        double sum = 0.0;
        double sumfo = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();

            // ignored cases
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // spline setup
            double ss = Crystal.invressq(crystal, ih);
            double fh = spline.f(ss, refinementdata.spline);

            ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
            if (refinementdata.isFreeR(i)) {
                sum += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo += abs(fo[i][0]);
            }
        }

        rfree = (sum / sumfo) * 100.0;
        return rfree;
    }

    /**
     * simply return the current sigmaA value
     *
     * @return sigmaA
     */
    public double getSigmaA() {
        double sum = 0.0;
        int nhkl = 0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                refinementdata.sigmaa.length);

        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();

            // ignored cases
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // spline setup
            double ss = Crystal.invressq(crystal, ih);
            double fh = spline.f(ss, refinementdata.spline);
            double sa = sigmaaspline.f(ss, refinementdata.sigmaa);

            nhkl++;
            sum += (sa - sum) / nhkl;
        }

        return sum;
    }

    /**
     * simply return the current sigmaW value
     *
     * @return sigmaW
     */
    public double getSigmaW() {
        double sum = 0.0;
        int nhkl = 0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                refinementdata.sigmaa.length);

        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();

            // ignored cases
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // spline setup
            double ss = Crystal.invressq(crystal, ih);
            double fh = spline.f(ss, refinementdata.spline);
            double wa = sigmaaspline.f(ss, refinementdata.sigmaw);

            nhkl++;
            sum += (wa - sum) / nhkl;
        }

        return sum;
    }

    /**
     * output Cruickshank and Blow DPI indices
     *
     * @param natoms number of atoms in the structure
     * @param nnonhatoms number of non-H atoms in the structure
     * @see <a href="http://dx.doi.org/10.1107/S0907444998012645"
     * target="_blank"> D. W. J. Cruickshank, Acta Cryst. (1999). D55,
     * 583-601</a>
     * @see <a href="http://dx.doi.org/10.1107/S0907444902003931"
     * target="_blank"> D. M. Blow, Acta Cryst. (2002). D58, 792-797</a>
     * @see <a href="http://dx.doi.org/10.1107/S0907444998012645"
     * target="_blank"> D. W. J. Cruickshank, Acta Cryst. (1999). D55,
     * 583-601</a>
     * @see <a href="http://dx.doi.org/10.1107/S0907444902003931"
     * target="_blank"> D. M. Blow, Acta Cryst. (2002). D58, 792-797</a>
     */
    public void printDPIStats(int natoms, int nnonhatoms) {
        int nhkli = 0;
        int nhklo = refinementdata.n;
        double rfreefrac = getRFree() * 0.01;
        double res = reflectionlist.resolution.resolutionLimit();
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();

            // ignored cases
            if (Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }
            nhkli++;
        }

        double va = Math.pow(crystal.volume / crystal.spaceGroup.getNumberOfSymOps(), 0.3333);
        blowdpih = 1.28 * Math.sqrt(natoms) * va
                * Math.pow(nhkli, -0.8333) * rfreefrac;

        blowdpi = 1.28 * Math.sqrt(nnonhatoms) * va
                * Math.pow(nhkli, -0.8333) * rfreefrac;

        double natni = Math.sqrt((double) natoms / nhkli);
        double noni = Math.pow((double) nhkli / nhklo, -0.3333);
        cruickdpih = natni * noni * rfreefrac * res;

        natni = Math.sqrt((double) nnonhatoms / nhkli);
        cruickdpi = natni * noni * rfreefrac * res;

        StringBuilder sb = new StringBuilder("\n");
        sb.append(String.format(" Blow DPI for all / non-H atoms:        %7.4f / %7.4f\n", blowdpih, blowdpi));
        sb.append(String.format(" Cruickshank DPI for all / non-H atoms: %7.4f / %7.4f", cruickdpih, cruickdpi));

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * print HKL statistics/completeness info
     */
    public void printHKLStats() {
        double res[][] = new double[n][2];
        int nhkl[][] = new int[n][3];

        for (int i = 0; i < n; i++) {
            res[i][0] = Double.NEGATIVE_INFINITY;
            res[i][1] = Double.POSITIVE_INFINITY;
        }

        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // ignored cases
            if (Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                nhkl[b][2]++;
                continue;
            }

            // determine res limits of each bin
            double rh = Crystal.res(crystal, ih);
            if (rh > res[b][0]) {
                res[b][0] = rh;
            }
            if (rh < res[b][1]) {
                res[b][1] = rh;
            }

            // count the reflection
            if (freer[i] == refinementdata.rfreeflag) {
                nhkl[b][1]++;
            } else {
                nhkl[b][0]++;
            }
        }

        StringBuilder sb = new StringBuilder(String.format("\n %15s | %8s|%9s| %7s | %7s | %s\n",
                "Res. Range", " HKL (R)", " HKL (cv)", " Bin", " Miss",
                "Complete (%)"));
        for (int i = 0; i < n; i++) {
            sb.append(String.format(" %7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(String.format("%7d | %7d | %7d | %7d | ",
                    nhkl[i][0], nhkl[i][1], nhkl[i][0] + nhkl[i][1],
                    nhkl[i][2]));
            sb.append(String.format("%6.2f\n", (((double) nhkl[i][0] + nhkl[i][1])
                    / (nhkl[i][0] + nhkl[i][1] + nhkl[i][2])) * 100.0));
        }
        sb.append(String.format(" %7.3f %7.3f | ", res[0][0], res[n - 1][1]));
        int sum1 = 0;
        int sum2 = 0;
        int sum3 = 0;
        for (int i = 0; i < n; i++) {
            sum1 += nhkl[i][0];
            sum2 += nhkl[i][1];
            sum3 += nhkl[i][2];
        }
        sb.append(String.format("%7d | %7d | %7d | %7d | ",
                sum1, sum2, sum1 + sum2, sum3));
        sb.append(String.format("%6.2f\n",
                (((double) sum1 + sum2) / (sum1 + sum2 + sum3)) * 100.0));
        sb.append(String.format(" Number of reflections if complete: %10d", refinementdata.n));

        nobshkl = sum1 + sum2;
        highnobshkl = nhkl[n - 1][0] + nhkl[n - 1][1];
        nobsrfree = sum2;
        highnobsrfree = nhkl[n - 1][1];
        completeness = (((double) sum1 + sum2) / (sum1 + sum2 + sum3)) * 100.0;
        highcompleteness = (((double) nhkl[n - 1][0] + nhkl[n - 1][1])
                / (nhkl[n - 1][0] + nhkl[n - 1][1] + nhkl[n - 1][2])) * 100.0;

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * print R factors and associated statistics in a binned fashion
     */
    public void printRStats() {
        double res[][] = new double[n][2];
        double nhkl[] = new double[n + 1];
        double rb[][] = new double[n + 1][2];
        double sumfo[][] = new double[n + 1][2];
        double s[][] = new double[n + 1][4];
        double numer = 0.0;
        double denom = 0.0;
        double sumall = 0.0;
        double sumfoall = 0.0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                refinementdata.sigmaa.length);

        for (int i = 0; i < n; i++) {
            res[i][0] = Double.NEGATIVE_INFINITY;
            res[i][1] = Double.POSITIVE_INFINITY;
        }

        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // ignored cases
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // spline setup
            double ss = Crystal.invressq(crystal, ih);
            double fh = spline.f(ss, refinementdata.spline);
            double sa = sigmaaspline.f(ss, refinementdata.sigmaa);
            double wa = sigmaaspline.f(ss, refinementdata.sigmaw);
            double eoscale = sigmaaspline.f(ss, refinementdata.foesq);

            // determine res limits of each bin
            double rs = Crystal.res(crystal, ih);
            if (rs > res[b][0]) {
                res[b][0] = rs;
            }
            if (rs < res[b][1]) {
                res[b][1] = rs;
            }

            ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
            numer = abs(abs(fo[i][0]) - fh * abs(c.abs()));
            denom = abs(fo[i][0]);
            if (refinementdata.isFreeR(i)) {
                rb[b][1] += numer;
                sumfo[b][1] += denom;
                rb[n][1] += numer;
                sumfo[n][1] += denom;
                sumall += numer;
                sumfoall += denom;
            } else {
                rb[b][0] += numer;
                sumfo[b][0] += denom;
                rb[n][0] += numer;
                sumfo[n][0] += denom;
                sumall += numer;
                sumfoall += denom;
            }

            nhkl[b]++;
            nhkl[n]++;
            s[b][0] += (sa - s[b][0]) / nhkl[b];
            s[b][1] += (wa - s[b][1]) / nhkl[b];
            s[b][2] += ((wa / Math.sqrt(eoscale)) - s[b][2]) / nhkl[b];
            s[b][3] += (fomphi[i][0] - s[b][3]) / nhkl[b];
            s[n][0] += (sa - s[n][0]) / nhkl[n];
            s[n][1] += (wa - s[n][1]) / nhkl[n];
            s[n][2] += ((wa / Math.sqrt(eoscale)) - s[n][2]) / nhkl[n];
            s[n][3] += (fomphi[i][0] - s[n][3]) / nhkl[n];
        }

        StringBuilder sb = new StringBuilder(
                String.format("\n %15s | %7s | %7s | %7s | %7s | %7s | %7s\n",
                        "Res. Range", "  R", "Rfree", "s", "w(E)", "w(F)", "FOM"));
        for (int i = 0; i < n; i++) {
            sb.append(String.format(" %7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(String.format("%7.2f | %7.2f | %7.4f | %7.4f | %7.2f | %7.4f\n",
                    (rb[i][0] / sumfo[i][0]) * 100.0,
                    (rb[i][1] / sumfo[i][1]) * 100.0,
                    s[i][0], s[i][1], s[i][2], s[i][3]));
        }

        sb.append(String.format(" %7.3f %7.3f | ", res[0][0], res[n - 1][1]));
        sb.append(String.format("%7.2f | %7.2f | %7.4f | %7.4f | %7.2f | %7.4f\n",
                (rb[n][0] / sumfo[n][0]) * 100.0,
                (rb[n][1] / sumfo[n][1]) * 100.0,
                s[n][0], s[n][1], s[n][2], s[n][3]));
        sb.append(" s and w are analagous to D and sum_wc");

        reslow = res[0][0];
        reshigh = res[n - 1][1];
        highreslow = res[n - 1][0];
        highreshigh = res[n - 1][1];
        r = (rb[n][0] / sumfo[n][0]) * 100.0;
        rfree = (rb[n][1] / sumfo[n][1]) * 100.0;
        rall = (sumall / sumfoall) * 100.0;
        highr = (rb[n - 1][0] / sumfo[n - 1][0]) * 100.0;
        highrfree = (rb[n - 1][1] / sumfo[n - 1][1]) * 100.0;

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * print scaling and bulk solvent statistics
     */
    public void printScaleStats() {
        int nhkl[] = new int[n];
        double scale[] = new double[n];

        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // ignored cases
            if (Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // spline setup
            double ss = Crystal.invressq(crystal, ih);
            double fh = spline.f(ss, refinementdata.spline);

            nhkl[b]++;
            scale[b] += (fh - scale[b]) / nhkl[b];
        }

        StringBuilder sb = new StringBuilder(
                String.format(" Fc to Fo scale: %4.2f\n",
                        exp(0.25 * refinementdata.model_k)));
        sb.append(" Fc to Fo spline scale: ");
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%4.2f ", scale[i]));
        }
        sb.append(String.format("\n Aniso B tensor:\n"));
        sb.append(String.format("  %10.4f %10.4f %10.4f\n",
                refinementdata.model_b[0],
                refinementdata.model_b[3],
                refinementdata.model_b[4]));
        sb.append(String.format("  %10.4f %10.4f %10.4f\n",
                refinementdata.model_b[3],
                refinementdata.model_b[1],
                refinementdata.model_b[5]));
        sb.append(String.format("  %10.4f %10.4f %10.4f\n",
                refinementdata.model_b[4],
                refinementdata.model_b[5],
                refinementdata.model_b[2]));
        if (refinementdata.crs_fs.solventModel != SolventModel.NONE) {
            switch (refinementdata.crs_fs.solventModel) {
                case BINARY:
                    sb.append(" Bulk solvent model: Binary mask\n");
                    sb.append(String.format("  Probe radius: %8.3f\n  Shrink radius: %8.3f\n",
                            refinementdata.solvent_a,
                            refinementdata.solvent_b));
                    break;
                case POLYNOMIAL:
                    sb.append(" Bulk solvent model: Polynomial switch\n");
                    sb.append(String.format("  a:     %8.3f\n  w:     %8.3f\n",
                            refinementdata.solvent_a,
                            refinementdata.solvent_b));
                    break;
                case GAUSSIAN:
                    sb.append(" Bulk solvent model: Gaussian\n");
                    sb.append(String.format("  A: %8.3f\n  sd scale: %8.3f\n",
                            refinementdata.solvent_a,
                            refinementdata.solvent_b));
                    break;
            }
            sb.append(String.format("  Scale: %8.3f\n  B:     %8.3f\n",
                    refinementdata.solvent_k,
                    refinementdata.solvent_ueq * 8.0 * Math.PI * Math.PI));
        }
        sb.append(String.format("\n -log Likelihood: %14.3f (free set: %14.3f)",
                refinementdata.llkr, refinementdata.llkf));

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * print signal to noise ratio statistics
     */
    public void printSNStats() {
        double res[][] = new double[n][2];
        double nhkl[] = new double[n + 1];
        double sn[][] = new double[n + 1][3];

        for (int i = 0; i < n; i++) {
            res[i][0] = Double.NEGATIVE_INFINITY;
            res[i][1] = Double.POSITIVE_INFINITY;
        }

        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // ignored cases
            if (Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            // determine res limits of each bin
            double rs = Crystal.res(crystal, ih);
            if (rs > res[b][0]) {
                res[b][0] = rs;
            }
            if (rs < res[b][1]) {
                res[b][1] = rs;
            }

            // running mean
            nhkl[b]++;
            nhkl[n]++;
            sn[b][0] += (fo[i][0] - sn[b][0]) / nhkl[b];
            sn[b][1] += (fo[i][1] - sn[b][1]) / nhkl[b];
            sn[b][2] += ((fo[i][0] / fo[i][1]) - sn[b][2]) / nhkl[b];

            sn[n][0] += (fo[i][0] - sn[n][0]) / nhkl[b];
            sn[n][1] += (fo[i][1] - sn[n][1]) / nhkl[b];
            sn[n][2] += ((fo[i][0] / fo[i][1]) - sn[n][2]) / nhkl[n];
        }

        StringBuilder sb = new StringBuilder(String.format("\n %15s | %7s | %7s | %7s \n",
                "Res. Range", "Signal", "Sigma", "S/N"));
        for (int i = 0; i < n; i++) {
            sb.append(String.format(" %7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(String.format("%7.2f | %7.2f | %7.2f\n",
                    sn[i][0], sn[i][1], sn[i][2]));
        }

        sb.append(String.format(" %7.3f %7.3f | ", res[0][0], res[n - 1][1]));
        sb.append(String.format("%7.2f | %7.2f | %7.2f",
                sn[n][0], sn[n][1], sn[n][2]));

        if (print) {
            logger.info(sb.toString());
        }
    }
}
