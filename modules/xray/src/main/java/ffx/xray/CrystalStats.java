//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.xray;

import java.util.logging.Logger;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.math.ComplexNumber;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import static ffx.crystal.Crystal.invressq;
import static ffx.crystal.Crystal.res;

/**
 * Crystal statistics output/logger
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class CrystalStats {

    private static final Logger logger = Logger.getLogger(CrystalStats.class.getName());
    private final ReflectionList reflectionList;
    private final DiffractionRefinementData refinementData;
    private final Crystal crystal;
    private final int nBins;
    private final double[][] fSigF;
    private final int[] freeR;
    private final double[][] fcTot;
    private final double[][] fomPhi;
    private final ReflectionSpline spline;

    private double blowDPI = -1.0;
    private boolean print = true;

    private int nObsHKL, highnObsHKL, nObsRFree, highnObsRFree;
    private double resHigh, resLow, highResHigh, highResLow;
    private double completeness, highCompleteness;
    private double rall, r, rfree, highR, highRfree;
    private double blowDPIH, cruickDPI, cruickDPIH;


    /**
     * constructor
     *
     * @param reflectionList {@link ffx.crystal.ReflectionList} to use for logging
     * @param refinementData {@link ffx.xray.DiffractionRefinementData} to use for logging
     */
    CrystalStats(ReflectionList reflectionList,
                 DiffractionRefinementData refinementData) {
        this.reflectionList = reflectionList;
        this.refinementData = refinementData;
        this.crystal = reflectionList.crystal;
        this.nBins = refinementData.nBins;
        this.fSigF = refinementData.fSigF;
        this.freeR = refinementData.freeR;
        this.fcTot = refinementData.fcTot;
        this.fomPhi = refinementData.fomPhi;
        this.spline = new ReflectionSpline(reflectionList, refinementData.spline.length);
    }

    /**
     * <p>
     * getPDBHeaderString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    String getPDBHeaderString() {
        print = false;
        printHKLStats();
        printRStats();
        print = true;

        StringBuilder sb = new StringBuilder();

        sb.append("REMARK   3  DATA USED IN REFINEMENT\n");
        sb.append(format("REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : %6.2f\n", resHigh));
        sb.append(format("REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : %6.2f\n", resLow));
        if (refinementData.fSigFCutoff > 0.0) {
            sb.append(format("REMARK   3   DATA CUTOFF            (SIGMA(F)) : %6.2f\n", refinementData.fSigFCutoff));
        } else {
            sb.append("REMARK   3   DATA CUTOFF            (SIGMA(F)) : NONE\n");
        }
        sb.append(format("REMARK   3   COMPLETENESS FOR RANGE        (%%) : %6.2f\n", completeness));
        sb.append("REMARK   3   NUMBER OF REFLECTIONS             : " + nObsHKL + "\n");
        sb.append("REMARK   3\n");
        sb.append("REMARK   3  FIT TO DATA USED IN REFINEMENT\n");
        sb.append("REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT\n");
        sb.append("REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM\n");
        sb.append(format("REMARK   3   R VALUE               (OBSERVED) : %8.6f\n", rall / 100.0));
        sb.append(format("REMARK   3   R VALUE            (WORKING SET) : %8.6f\n", r / 100.0));
        sb.append(format("REMARK   3   FREE R VALUE                     : %8.6f\n", rfree / 100.0));
        sb.append(format("REMARK   3   FREE R VALUE TEST SET SIZE   (%%) : %6.1f\n", (((double) nObsRFree) / nObsHKL) * 100.0));
        sb.append("REMARK   3   FREE R VALUE TEST SET COUNT      : " + nObsRFree + "\n");
        sb.append("REMARK   3\n");
        sb.append("REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN\n");
        sb.append("REMARK   3   TOTAL NUMBER OF BINS USED           : " + nBins + "\n");
        sb.append(format("REMARK   3   BIN RESOLUTION RANGE HIGH           : %6.2f\n", highResHigh));
        sb.append(format("REMARK   3   BIN RESOLUTION RANGE LOW            : %6.2f\n", highResLow));
        sb.append("REMARK   3   REFLECTION IN BIN     (WORKING SET) : " + highnObsHKL + "\n");
        sb.append(format("REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%%) : %6.2f\n", highCompleteness));
        sb.append(format("REMARK   3   BIN R VALUE           (WORKING SET) : %8.6f\n", highR / 100.0));
        sb.append("REMARK   3   BIN FREE R VALUE SET COUNT          : " + highnObsRFree + "\n");
        sb.append(format("REMARK   3   BIN FREE R VALUE                    : %8.6f\n", highRfree / 100.0));
        sb.append("REMARK   3\n");

        sb.append("REMARK   3  OVERALL SCALE FACTORS\n");
        sb.append(format("REMARK   3   SCALE: %4.2f\n",
                exp(0.25 * refinementData.modelScaleK)));
        sb.append("REMARK   3   ANISOTROPIC SCALE TENSOR:\n");
        sb.append(format("REMARK   3    %g %g %g\n",
                refinementData.modelAnisoB[0],
                refinementData.modelAnisoB[3],
                refinementData.modelAnisoB[4]));
        sb.append(format("REMARK   3    %g %g %g\n",
                refinementData.modelAnisoB[3],
                refinementData.modelAnisoB[1],
                refinementData.modelAnisoB[5]));
        sb.append(format("REMARK   3    %g %g %g\n",
                refinementData.modelAnisoB[4],
                refinementData.modelAnisoB[5],
                refinementData.modelAnisoB[2]));
        sb.append("REMARK   3\n");

        if (refinementData.crystalReciprocalSpaceFs.solventModel != SolventModel.NONE) {
            sb.append("REMARK   3  BULK SOLVENT MODELLING\n");
            switch (refinementData.crystalReciprocalSpaceFs.solventModel) {
                case BINARY:
                    sb.append("REMARK   3   METHOD USED: BINARY MASK\n");
                    sb.append(format("REMARK   3    PROBE RADIUS  : %g\n",
                            refinementData.solventA));
                    sb.append(format("REMARK   3    SHRINK RADIUS : %g\n",
                            refinementData.solventB));
                    break;
                case POLYNOMIAL:
                    sb.append("REMARK   3   METHOD USED: POLYNOMIAL SWITCH\n");
                    sb.append(format("REMARK   3    ATOMIC RADIUS BUFFER : %g\n",
                            refinementData.solventA));
                    sb.append(format("REMARK   3    SWITCH RADIUS        : %g\n",
                            refinementData.solventB));
                    break;
                case GAUSSIAN:
                    sb.append("REMARK   3   METHOD USED: GAUSSIAN\n");
                    sb.append(format("REMARK   3    ATOMIC RADIUS BUFFER : %g\n",
                            refinementData.solventA));
                    sb.append(format("REMARK   3    STD DEV SCALE        : %g\n",
                            refinementData.solventB));
                    break;
            }
            sb.append(format("REMARK   3    K_SOL: %g\n",
                    refinementData.bulkSolventK));
            sb.append(format("REMARK   3    B_SOL: %g\n",
                    refinementData.bulkSolventUeq * 8.0 * Math.PI * Math.PI));
            sb.append("REMARK   3\n");
        }

        if (blowDPI > 0.0) {
            sb.append("REMARK   3  ERROR ESTIMATES\n");
            sb.append("REMARK   3   ACTA CRYST (1999) D55, 583-601\n");
            sb.append("REMARK   3   ACTA CRYST (2002) D58, 792-797\n");
            sb.append(format("REMARK   3   BLOW DPI ALL ATOMS (EQN 7)          : %7.4f\n", blowDPIH));
            sb.append(format("REMARK   3   BLOW DPI NONH ATOMS (EQN 7)         : %7.4f\n", blowDPI));
            sb.append(format("REMARK   3   CRUICKSHANK DPI ALL ATOMS (EQN 27)  : %7.4f\n", cruickDPIH));
            sb.append(format("REMARK   3   CRUICKSHANK DPI NONH ATOMS (EQN 27) : %7.4f\n", cruickDPI));
            sb.append("REMARK   3\n");
        }

        sb.append("REMARK   3  DATA TARGET\n");
        sb.append("REMARK   3   METHOD USED: MAXIMUM LIKELIHOOD\n");
        sb.append(format("REMARK   3    -LOG LIKELIHOOD            : %g\n",
                refinementData.llkR));
        sb.append(format("REMARK   3    -LOG LIKELIHOOD (FREE SET) : %g\n",
                refinementData.llkF));
        sb.append("REMARK   3\n");

        return sb.toString();
    }

    /**
     * Simply return the current R value.
     *
     * @return R value as a percent.
     */
    public double getR() {
        double numer;
        double denom;
        double sum = 0.0;
        double sumfo = 0.0;
        double sumall = 0.0;
        double sumfoall = 0.0;
        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();

            // Ignored cases
            if (isNaN(fcTot[i][0])
                    || isNaN(fSigF[i][0])
                    || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Spline setup
            double ss = invressq(crystal, ih);
            double fh = spline.f(ss, refinementData.spline);

            ComplexNumber c = new ComplexNumber(fcTot[i][0], fcTot[i][1]);
            numer = abs(abs(fSigF[i][0]) - fh * abs(c.abs()));
            denom = abs(fSigF[i][0]);
            sumall += numer;
            sumfoall += denom;
            if (!refinementData.isFreeR(i)) {
                sum += numer;
                sumfo += denom;
            }
        }

        rall = (sumall / sumfoall) * 100.0;
        r = (sum / sumfo) * 100.0;
        return r;
    }

    /**
     * Simply return the current R free value.
     *
     * @return R free value as a percent.
     */
    double getRFree() {
        double sum = 0.0;
        double sumfo = 0.0;
        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();

            // Ignored cases
            if (isNaN(fcTot[i][0])
                    || isNaN(fSigF[i][0])
                    || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Spline setup
            double ss = invressq(crystal, ih);
            double fh = spline.f(ss, refinementData.spline);

            ComplexNumber c = new ComplexNumber(fcTot[i][0], fcTot[i][1]);
            if (refinementData.isFreeR(i)) {
                sum += abs(abs(fSigF[i][0]) - fh * abs(c.abs()));
                sumfo += abs(fSigF[i][0]);
            }
        }

        rfree = (sum / sumfo) * 100.0;
        return rfree;
    }

    /**
     * Simply return the current sigmaA value.
     *
     * @return sigmaA
     */
    public double getSigmaA() {
        double sum = 0.0;
        int nhkl = 0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionList,
                refinementData.sigmaA.length);

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();

            // Ignored cases
            if (isNaN(fcTot[i][0]) || isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Spline setup
            double ss = invressq(crystal, ih);
            spline.f(ss, refinementData.spline);
            double sa = sigmaaspline.f(ss, refinementData.sigmaA);

            nhkl++;
            sum += (sa - sum) / nhkl;
        }

        return sum;
    }

    /**
     * Simply return the current sigmaW value.
     *
     * @return sigmaW
     */
    double getSigmaW() {
        double sum = 0.0;
        int nhkl = 0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionList,
                refinementData.sigmaA.length);

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();

            // Ignored cases
            if (isNaN(fcTot[i][0]) || isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Spline setup
            double ss = invressq(crystal, ih);
            spline.f(ss, refinementData.spline);
            double wa = sigmaaspline.f(ss, refinementData.sigmaW);

            nhkl++;
            sum += (wa - sum) / nhkl;
        }

        return sum;
    }

    /**
     * Output Cruickshank and Blow DPI indices.
     *
     * @param nAtoms            number of atoms in the structure
     * @param nNonHydrogenAtoms number of non-H atoms in the structure
     * @see <a href="http://dx.doi.org/10.1107/S0907444998012645"
     * target="_blank"> D. W. J. Cruickshank, Acta Cryst. (1999). D55,
     * 583-601</a>
     * @see <a href="http://dx.doi.org/10.1107/S0907444902003931"
     * target="_blank"> D. M. Blow, Acta Cryst. (2002). D58, 792-797</a>
     */
    void printDPIStats(int nAtoms, int nNonHydrogenAtoms) {
        int nhkli = 0;
        int nhklo = refinementData.n;
        double rfreefrac = getRFree() * 0.01;
        double res = reflectionList.resolution.resolutionLimit();
        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();

            // ignored cases
            if (isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                continue;
            }
            nhkli++;
        }

        double va = pow(crystal.volume / crystal.spaceGroup.getNumberOfSymOps(), 0.3333);
        blowDPIH = 1.28 * sqrt(nAtoms) * va * pow(nhkli, -0.8333) * rfreefrac;
        blowDPI = 1.28 * sqrt(nNonHydrogenAtoms) * va * pow(nhkli, -0.8333) * rfreefrac;
        double natni = sqrt((double) nAtoms / nhkli);
        double noni = pow((double) nhkli / nhklo, -0.3333);
        cruickDPIH = natni * noni * rfreefrac * res;
        natni = sqrt((double) nNonHydrogenAtoms / nhkli);
        cruickDPI = natni * noni * rfreefrac * res;

        StringBuilder sb = new StringBuilder("\n");
        sb.append(format(" Blow DPI for all / non-H atoms:        %7.4f / %7.4f\n", blowDPIH, blowDPI));
        sb.append(format(" Cruickshank DPI for all / non-H atoms: %7.4f / %7.4f", cruickDPIH, cruickDPI));

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * Print HKL statistics/completeness info.
     */
    void printHKLStats() {
        double[][] res = new double[nBins][2];
        int[][] nhkl = new int[nBins][3];

        for (int i = 0; i < nBins; i++) {
            res[i][0] = Double.NEGATIVE_INFINITY;
            res[i][1] = Double.POSITIVE_INFINITY;
        }

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // Ignored cases
            if (isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                nhkl[b][2]++;
                continue;
            }

            // Determine res limits of each bin
            double rh = res(crystal, ih);
            if (rh > res[b][0]) {
                res[b][0] = rh;
            }
            if (rh < res[b][1]) {
                res[b][1] = rh;
            }

            // Count the reflection
            if (freeR[i] == refinementData.rFreeFlag) {
                nhkl[b][1]++;
            } else {
                nhkl[b][0]++;
            }
        }

        StringBuilder sb = new StringBuilder(format("\n %15s | %8s|%9s| %7s | %7s | %s\n",
                "Res. Range", " HKL (R)", " HKL (cv)", " Bin", " Miss",
                "Complete (%)"));
        for (int i = 0; i < nBins; i++) {
            sb.append(format(" %7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(format("%7d | %7d | %7d | %7d | ", nhkl[i][0], nhkl[i][1], nhkl[i][0] + nhkl[i][1], nhkl[i][2]));
            sb.append(format("%6.2f\n", (((double) nhkl[i][0] + nhkl[i][1]) / (nhkl[i][0] + nhkl[i][1] + nhkl[i][2])) * 100.0));
        }
        sb.append(format(" %7.3f %7.3f | ", res[0][0], res[nBins - 1][1]));
        int sum1 = 0;
        int sum2 = 0;
        int sum3 = 0;
        for (int i = 0; i < nBins; i++) {
            sum1 += nhkl[i][0];
            sum2 += nhkl[i][1];
            sum3 += nhkl[i][2];
        }
        sb.append(format("%7d | %7d | %7d | %7d | ", sum1, sum2, sum1 + sum2, sum3));
        sb.append(format("%6.2f\n", (((double) sum1 + sum2) / (sum1 + sum2 + sum3)) * 100.0));
        sb.append(format(" Number of reflections if complete: %10d", refinementData.n));

        nObsHKL = sum1 + sum2;
        highnObsHKL = nhkl[nBins - 1][0] + nhkl[nBins - 1][1];
        nObsRFree = sum2;
        highnObsRFree = nhkl[nBins - 1][1];
        completeness = (((double) sum1 + sum2) / (sum1 + sum2 + sum3)) * 100.0;
        highCompleteness = (((double) nhkl[nBins - 1][0] + nhkl[nBins - 1][1])
                / (nhkl[nBins - 1][0] + nhkl[nBins - 1][1] + nhkl[nBins - 1][2])) * 100.0;

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * Print R factors and associated statistics in a binned fashion.
     */
    void printRStats() {
        double[][] res = new double[nBins][2];
        double[] nhkl = new double[nBins + 1];
        double[][] rb = new double[nBins + 1][2];
        double[][] sumfo = new double[nBins + 1][2];
        double[][] s = new double[nBins + 1][4];
        double numer;
        double denom;
        double sumall = 0.0;
        double sumfoall = 0.0;
        ReflectionSpline sigmaASpline = new ReflectionSpline(reflectionList, refinementData.sigmaA.length);

        for (int i = 0; i < nBins; i++) {
            res[i][0] = Double.NEGATIVE_INFINITY;
            res[i][1] = Double.POSITIVE_INFINITY;
        }

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // Ignored cases
            if (isNaN(fcTot[i][0]) || isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Spline setup
            double ss = invressq(crystal, ih);
            double fh = spline.f(ss, refinementData.spline);
            double sa = sigmaASpline.f(ss, refinementData.sigmaA);
            double wa = sigmaASpline.f(ss, refinementData.sigmaW);
            double eoscale = sigmaASpline.f(ss, refinementData.esqFo);

            // Determine res limits of each bin
            double rs = res(crystal, ih);
            if (rs > res[b][0]) {
                res[b][0] = rs;
            }
            if (rs < res[b][1]) {
                res[b][1] = rs;
            }

            ComplexNumber c = new ComplexNumber(fcTot[i][0], fcTot[i][1]);
            numer = abs(abs(fSigF[i][0]) - fh * abs(c.abs()));
            denom = abs(fSigF[i][0]);
            if (refinementData.isFreeR(i)) {
                rb[b][1] += numer;
                sumfo[b][1] += denom;
                rb[nBins][1] += numer;
                sumfo[nBins][1] += denom;
                sumall += numer;
                sumfoall += denom;
            } else {
                rb[b][0] += numer;
                sumfo[b][0] += denom;
                rb[nBins][0] += numer;
                sumfo[nBins][0] += denom;
                sumall += numer;
                sumfoall += denom;
            }

            nhkl[b]++;
            nhkl[nBins]++;
            s[b][0] += (sa - s[b][0]) / nhkl[b];
            s[b][1] += (wa - s[b][1]) / nhkl[b];
            s[b][2] += ((wa / sqrt(eoscale)) - s[b][2]) / nhkl[b];
            s[b][3] += (fomPhi[i][0] - s[b][3]) / nhkl[b];
            s[nBins][0] += (sa - s[nBins][0]) / nhkl[nBins];
            s[nBins][1] += (wa - s[nBins][1]) / nhkl[nBins];
            s[nBins][2] += ((wa / Math.sqrt(eoscale)) - s[nBins][2]) / nhkl[nBins];
            s[nBins][3] += (fomPhi[i][0] - s[nBins][3]) / nhkl[nBins];
        }

        StringBuilder sb = new StringBuilder(format("\n %15s | %7s | %7s | %7s | %7s | %7s | %7s\n",
                "Res. Range", "  R", "Rfree", "s", "w(E)", "w(F)", "FOM"));
        for (int i = 0; i < nBins; i++) {
            sb.append(format(" %7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(format("%7.2f | %7.2f | %7.4f | %7.4f | %7.2f | %7.4f\n",
                    (rb[i][0] / sumfo[i][0]) * 100.0,
                    (rb[i][1] / sumfo[i][1]) * 100.0,
                    s[i][0], s[i][1], s[i][2], s[i][3]));
        }

        sb.append(format(" %7.3f %7.3f | ", res[0][0], res[nBins - 1][1]));
        sb.append(format("%7.2f | %7.2f | %7.4f | %7.4f | %7.2f | %7.4f\n",
                (rb[nBins][0] / sumfo[nBins][0]) * 100.0,
                (rb[nBins][1] / sumfo[nBins][1]) * 100.0,
                s[nBins][0], s[nBins][1], s[nBins][2], s[nBins][3]));
        sb.append(" s and w are analagous to D and sum_wc");

        resLow = res[0][0];
        resHigh = res[nBins - 1][1];
        highResLow = res[nBins - 1][0];
        highResHigh = res[nBins - 1][1];
        r = (rb[nBins][0] / sumfo[nBins][0]) * 100.0;
        rfree = (rb[nBins][1] / sumfo[nBins][1]) * 100.0;
        rall = (sumall / sumfoall) * 100.0;
        highR = (rb[nBins - 1][0] / sumfo[nBins - 1][0]) * 100.0;
        highRfree = (rb[nBins - 1][1] / sumfo[nBins - 1][1]) * 100.0;

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * Print scaling and bulk solvent statistics.
     */
    void printScaleStats() {
        int[] nhkl = new int[nBins];
        double[] scale = new double[nBins];

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // Ignored cases.
            if (isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Spline setup.
            double ss = invressq(crystal, ih);
            double fh = spline.f(ss, refinementData.spline);

            nhkl[b]++;
            scale[b] += (fh - scale[b]) / nhkl[b];
        }

        StringBuilder sb = new StringBuilder(format(" Fc to Fo scale: %4.2f\n",
                exp(0.25 * refinementData.modelScaleK)));
        sb.append(" Fc to Fo spline scale: ");
        for (int i = 0; i < nBins; i++) {
            sb.append(format("%4.2f ", scale[i]));
        }
        sb.append("\n Aniso B tensor:\n");
        sb.append(format("  %10.4f %10.4f %10.4f\n",
                refinementData.modelAnisoB[0],
                refinementData.modelAnisoB[3],
                refinementData.modelAnisoB[4]));
        sb.append(format("  %10.4f %10.4f %10.4f\n",
                refinementData.modelAnisoB[3],
                refinementData.modelAnisoB[1],
                refinementData.modelAnisoB[5]));
        sb.append(format("  %10.4f %10.4f %10.4f\n",
                refinementData.modelAnisoB[4],
                refinementData.modelAnisoB[5],
                refinementData.modelAnisoB[2]));
        if (refinementData.crystalReciprocalSpaceFs.solventModel != SolventModel.NONE) {
            switch (refinementData.crystalReciprocalSpaceFs.solventModel) {
                case BINARY:
                    sb.append(" Bulk solvent model: Binary mask\n");
                    sb.append(format("  Probe radius: %8.3f\n  Shrink radius: %8.3f\n",
                            refinementData.solventA,
                            refinementData.solventB));
                    break;
                case POLYNOMIAL:
                    sb.append(" Bulk solvent model: Polynomial switch\n");
                    sb.append(format("  a:     %8.3f\n  w:     %8.3f\n",
                            refinementData.solventA,
                            refinementData.solventB));
                    break;
                case GAUSSIAN:
                    sb.append(" Bulk solvent model: Gaussian\n");
                    sb.append(format("  A: %8.3f\n  sd scale: %8.3f\n",
                            refinementData.solventA,
                            refinementData.solventB));
                    break;
            }
            sb.append(format("  Scale: %8.3f\n  B:     %8.3f\n",
                    refinementData.bulkSolventK,
                    refinementData.bulkSolventUeq * 8.0 * Math.PI * Math.PI));
        }
        sb.append(format("\n -log Likelihood: %14.3f (free set: %14.3f)",
                refinementData.llkR, refinementData.llkF));

        if (print) {
            logger.info(sb.toString());
        }
    }

    /**
     * Print signal to noise ratio statistics.
     */
    void printSNStats() {
        double[][] res = new double[nBins][2];
        double[] nhkl = new double[nBins + 1];
        double[][] sn = new double[nBins + 1][3];

        for (int i = 0; i < nBins; i++) {
            res[i][0] = Double.NEGATIVE_INFINITY;
            res[i][1] = Double.POSITIVE_INFINITY;
        }

        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();
            int b = ih.bin();

            // Ignored cases
            if (isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
                continue;
            }

            // Determine res limits of each bin
            double rs = res(crystal, ih);
            if (rs > res[b][0]) {
                res[b][0] = rs;
            }
            if (rs < res[b][1]) {
                res[b][1] = rs;
            }

            // Running mean
            nhkl[b]++;
            nhkl[nBins]++;
            sn[b][0] += (fSigF[i][0] - sn[b][0]) / nhkl[b];
            sn[b][1] += (fSigF[i][1] - sn[b][1]) / nhkl[b];
            sn[b][2] += ((fSigF[i][0] / fSigF[i][1]) - sn[b][2]) / nhkl[b];

            sn[nBins][0] += (fSigF[i][0] - sn[nBins][0]) / nhkl[b];
            sn[nBins][1] += (fSigF[i][1] - sn[nBins][1]) / nhkl[b];
            sn[nBins][2] += ((fSigF[i][0] / fSigF[i][1]) - sn[nBins][2]) / nhkl[nBins];
        }

        StringBuilder sb = new StringBuilder(format("\n %15s | %7s | %7s | %7s \n",
                "Res. Range", "Signal", "Sigma", "S/N"));
        for (int i = 0; i < nBins; i++) {
            sb.append(format(" %7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(format("%7.2f | %7.2f | %7.2f\n", sn[i][0], sn[i][1], sn[i][2]));
        }

        sb.append(format(" %7.3f %7.3f | ", res[0][0], res[nBins - 1][1]));
        sb.append(format("%7.2f | %7.2f | %7.2f", sn[nBins][0], sn[nBins][1], sn[nBins][2]));

        if (print) {
            logger.info(sb.toString());
        }
    }
}
