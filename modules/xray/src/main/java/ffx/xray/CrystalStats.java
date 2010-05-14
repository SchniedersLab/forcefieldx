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

import static java.lang.Math.abs;
import static java.lang.Math.exp;

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.ComplexNumber;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

/**
 *
 * @author fennt
 */
public class CrystalStats {

    private static final Logger logger = Logger.getLogger(CrystalStats.class.getName());
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final Crystal crystal;
    private final ReflectionSpline spline;
    private final int n;
    private final double fo[][];
    private final int freer[];
    private final double fc[][];
    private final double fomphi[][];

    public CrystalStats(ReflectionList reflectionlist,
            RefinementData refinementdata) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.n = refinementdata.nparams;
        this.fo = refinementdata.fsigf;
        this.freer = refinementdata.freer;
        this.fc = refinementdata.fctot;
        this.fomphi = refinementdata.fomphi;

        this.spline = new ReflectionSpline(reflectionlist,
                refinementdata.spline.length);
    }

    public void print_hklstats() {
        double res[][] = new double[n][2];
        int nhkl[][] = new int[n][3];

        for (int i = 0; i < n; i++) {
            res[i][0] = Double.MIN_VALUE;
            res[i][1] = Double.MAX_VALUE;
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
            double r = Crystal.res(crystal, ih);
            if (r > res[b][0]) {
                res[b][0] = r;
            }
            if (r < res[b][1]) {
                res[b][1] = r;
            }

            // count the reflection
            if (freer[i] == refinementdata.rfreeflag) {
                nhkl[b][1]++;
            } else {
                nhkl[b][0]++;
            }
        }

        StringBuffer sb = new StringBuffer("\n");
        sb.append("# reflections (for 100% complete): " + refinementdata.n + "\n");
        sb.append(String.format("%15s | %8s|%9s| %7s | %7s |%s\n",
                "res. range", "#HKL (R)", "#HKL (cv)", "#bin", "#miss",
                "%complete"));
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(String.format("%7d | %7d | %7d | %7d | ",
                    nhkl[i][0], nhkl[i][1], nhkl[i][0] + nhkl[i][1],
                    nhkl[i][2]));
            sb.append(String.format("%5.2f\n", (((double) nhkl[i][0] + nhkl[i][1])
                    / (nhkl[i][0] + nhkl[i][1] + nhkl[i][2])) * 100.0));
        }

        sb.append(String.format("%7.3f %7.3f | ", res[0][0], res[n - 1][1]));
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
        sb.append(String.format("%5.2f\n\n",
                (((double) sum1 + sum2) / (sum1 + sum2 + sum3)) * 100.0));

        logger.info(sb.toString());
    }

    public double get_r() {
        double sum = 0.0;
        double sumfo = 0.0;
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

            ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
            if (!refinementdata.isfreer(i)) {
                sum += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo += abs(fo[i][0]);
            }
        }

        return (sum / sumfo) * 100.0;
    }

    public double get_rfree() {
        double sum = 0.0;
        double sumfo = 0.0;
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

            ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
            if (refinementdata.isfreer(i)) {
                sum += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo += abs(fo[i][0]);
            }
        }

        return (sum / sumfo) * 100.0;
    }

    public double get_sigmaa() {
        double sum = 0.0;
        int nhkl = 0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                refinementdata.sigmaa.length);

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

            nhkl++;
            sum += (sa - sum) / nhkl;
        }

        return sum;
    }

    public double get_sigmaw() {
        double sum = 0.0;
        int nhkl = 0;
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                refinementdata.sigmaa.length);

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
            double wa = sigmaaspline.f(ss, refinementdata.sigmaw);

            nhkl++;
            sum += (wa - sum) / nhkl;
        }

        return sum;
    }

    public void print_rstats() {
        double res[][] = new double[n][2];
        double nhkl[] = new double[n + 1];
        double r[][] = new double[n + 1][2];
        double sumfo[][] = new double[n + 1][2];
        double s[][] = new double[n + 1][3];
        ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                refinementdata.sigmaa.length);

        for (int i = 0; i < n; i++) {
            res[i][0] = Double.MIN_VALUE;
            res[i][1] = Double.MAX_VALUE;
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

            // determine res limits of each bin
            double rs = Crystal.res(crystal, ih);
            if (rs > res[b][0]) {
                res[b][0] = rs;
            }
            if (rs < res[b][1]) {
                res[b][1] = rs;
            }

            ComplexNumber c = new ComplexNumber(fc[i][0], fc[i][1]);
            if (refinementdata.isfreer(i)) {
                r[b][1] += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo[b][1] += abs(fo[i][0]);
                r[n][1] += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo[n][1] += abs(fo[i][0]);
            } else {
                r[b][0] += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo[b][0] += abs(fo[i][0]);
                r[n][0] += abs(abs(fo[i][0]) - fh * abs(c.abs()));
                sumfo[n][0] += abs(fo[i][0]);
            }

            nhkl[b]++;
            nhkl[n]++;
            s[b][0] += (sa - s[b][0]) / nhkl[b];
            s[b][1] += (wa - s[b][1]) / nhkl[b];
            s[b][2] += (fomphi[i][0] - s[b][2]) / nhkl[b];
            s[n][0] += (sa - s[n][0]) / nhkl[n];
            s[n][1] += (wa - s[n][1]) / nhkl[n];
            s[n][2] += (fomphi[i][0] - s[n][2]) / nhkl[n];
        }

        StringBuffer sb = new StringBuffer("\n");
        sb.append(String.format("%15s | %7s | %7s | %7s | %7s | %7s\n",
                "res. range", "  R", "Rfree", "s", "w", "FOM"));
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(String.format("%7.2f | %7.2f | %7.4f | %7.4f | %7.4f\n",
                    (r[i][0] / sumfo[i][0]) * 100.0,
                    (r[i][1] / sumfo[i][1]) * 100.0,
                    s[i][0], s[i][1], s[i][2]));
        }

        sb.append(String.format("%7.3f %7.3f | ", res[0][0], res[n - 1][1]));
        sb.append(String.format("%7.2f | %7.2f | %7.4f | %7.4f | %7.4f\n\n",
                (r[n][0] / sumfo[n][0]) * 100.0,
                (r[n][1] / sumfo[n][1]) * 100.0,
                s[n][0], s[n][1], s[n][2]));

        logger.info(sb.toString());
    }

    public void print_scalestats() {
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

        StringBuffer sb = new StringBuffer("\n");
        sb.append(String.format("  Fc to Fo scale: %4.2f\n",
                exp(0.25 * refinementdata.model_k)));
        sb.append("  Fc to Fo spline scale: ");
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%4.2f ", scale[i]));
        }
        sb.append("\n");
        sb.append(String.format("  aniso B tensor:\n"));
        sb.append(String.format("    %g %g %g\n",
                refinementdata.model_b[0],
                refinementdata.model_b[3],
                refinementdata.model_b[4]));
        sb.append(String.format("    %g %g %g\n",
                refinementdata.model_b[3],
                refinementdata.model_b[1],
                refinementdata.model_b[5]));
        sb.append(String.format("    %g %g %g\n",
                refinementdata.model_b[4],
                refinementdata.model_b[5],
                refinementdata.model_b[2]));
        if (refinementdata.bulksolvent) {
            if (refinementdata.crs_fs != null) {
                switch (refinementdata.crs_fs.solventmodel) {
                    case (SolventModel.BINARY):
                        sb.append("  bulk solvent model: binary mask\n");
                        sb.append(String.format("  bulk solvent probe radius: %g shrink radius: %g\n",
                                refinementdata.solvent_a,
                                refinementdata.solvent_sd));
                        break;
                    case (SolventModel.POLYNOMIAL):
                        sb.append("  bulk solvent model: polynomial switch\n");
                        sb.append(String.format("  bulk solvent atom radius: %g window size: %g\n",
                                refinementdata.solvent_a,
                                refinementdata.solvent_sd));
                        break;
                    case (SolventModel.GAUSSIAN):
                        sb.append("  bulk solvent model: Gaussian\n");
                        sb.append(String.format("  bulk solvent A: %g sd: %g\n",
                                refinementdata.solvent_a,
                                refinementdata.solvent_sd));
                        break;
                }
            }
            sb.append(String.format("  bulk solvent scale: %g  B: %g\n",
                    refinementdata.solvent_k,
                    refinementdata.solvent_ueq * 8.0 * Math.PI * Math.PI));
        }
        sb.append(String.format("  likelihood: %g (free set: %g)\n\n",
                refinementdata.llkr, refinementdata.llkf));
        logger.info(sb.toString());
    }

    public void print_snstats() {
        double res[][] = new double[n][2];
        double nhkl[] = new double[n + 1];
        double sn[][] = new double[n + 1][3];

        for (int i = 0; i < n; i++) {
            res[i][0] = Double.MIN_VALUE;
            res[i][1] = Double.MAX_VALUE;
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
            double r = Crystal.res(crystal, ih);
            if (r > res[b][0]) {
                res[b][0] = r;
            }
            if (r < res[b][1]) {
                res[b][1] = r;
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

        StringBuffer sb = new StringBuffer("\n");
        sb.append(String.format("%15s | %8s | %8s | %8s \n",
                "res. range", "signal", "sigma", "s/n"));
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%7.3f %7.3f | ", res[i][0], res[i][1]));
            sb.append(String.format("%8.2f | %8.2f | %8.2f\n",
                    sn[i][0], sn[i][1], sn[i][2]));
        }

        sb.append(String.format("%7.3f %7.3f | ", res[0][0], res[n - 1][1]));
        sb.append(String.format("%8.2f | %8.2f | %8.2f\n\n",
                sn[n][0], sn[n][1], sn[n][2]));

        logger.info(sb.toString());
    }
}
