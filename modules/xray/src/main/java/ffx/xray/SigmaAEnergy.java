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

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tanh;

import static org.apache.commons.math.util.FastMath.exp;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.ComplexNumber;
import ffx.numerics.Potential;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.mat3mat3;
import static ffx.numerics.VectorMath.mat3symvec6;
import static ffx.numerics.VectorMath.transpose3;
import static ffx.numerics.VectorMath.vec3mat3;
import java.util.Arrays;

/**
 *
 * Optimize sigmaA coefficients (using spline coefficients) and structure factor
 * derivatives using a likelihood target function
 *
 * this target can also be used for structure refinement
 *
 * @author Timothy D. Fenn<br>
 * @see <a href="http://dx.doi.org/10.1107/S0021889804031474" target="_blank">
 * K. Cowtan, J. Appl. Cryst. (2005). 38, 193-198</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444992007352" target="_blank">
 * A. T. Brunger, Acta Cryst. (1993). D49, 24-36</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444996012255" target="_blank">
 * G. N. Murshudov, A. A. Vagin and E. J. Dodson, Acta Cryst. (1997). D53,
 * 240-255</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767386099622" target="_blank">
 * R. J. Read, Acta Cryst. (1986). A42, 140-149.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767396004370" target="_blank">
 * N. S. Pannu and R. J. Read, Acta Cryst. (1996). A52, 659-668.</a>
 *
 */
public class SigmaAEnergy implements Potential {

    private static final Logger logger = Logger.getLogger(SigmaAEnergy.class.getName());
    private static final double twoPI2 = 2.0 * PI * PI;
    private static final double sim_a = 1.639294;
    private static final double sim_b = 3.553967;
    private static final double sim_c = 2.228716;
    private static final double sim_d = 3.524142;
    private static final double sim_e = 7.107935;
    private static final double sim_A = -1.28173889903;
    private static final double sim_B = 0.69231689903;
    private static final double sim_g = 2.13643992379;
    private static final double sim_p = 0.04613803811;
    private static final double sim_q = 1.82167089029;
    private static final double sim_r = -0.74817947490;

    private final ReflectionList reflectionList;
    private final Crystal crystal;
    private final DiffractionRefinementData refinementData;
    private final ParallelTeam parallelTeam;
    private final SigmaARegion sigmaARegion;
    private final double dfscale;
    private final double fo[][];
    private final double fctot[][];
    private final double fofc1[][];
    private final double fofc2[][];
    private final double fomphi[][];
    private final double dfc[][];
    private final double dfs[][];
    private final double recipt[][];
    private final double sa[];
    private final double wa[];
    private final int n;
    protected double[] optimizationScaling = null;
    private double totalEnergy;

    /**
     * <p>
     * Constructor for SigmaAEnergy.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData a {@link ffx.xray.DiffractionRefinementData}
     * object.
     */
    public SigmaAEnergy(ReflectionList reflectionList,
            DiffractionRefinementData refinementData, ParallelTeam parallelTeam) {
        this.reflectionList = reflectionList;
        this.crystal = reflectionList.crystal;
        this.refinementData = refinementData;
        this.fo = refinementData.fsigf;
        this.fctot = refinementData.fctot;
        this.fomphi = refinementData.fomphi;
        this.fofc1 = refinementData.fofc1;
        this.fofc2 = refinementData.fofc2;
        this.dfc = refinementData.dfc;
        this.dfs = refinementData.dfs;
        this.n = refinementData.nbins;

        // Initialize parameters.
        assert (refinementData.crs_fc != null);
        double fftgrid = 2.0 * refinementData.crs_fc.getXDim()
                * refinementData.crs_fc.getYDim()
                * refinementData.crs_fc.getZDim();
        dfscale = (crystal.volume * crystal.volume) / fftgrid;
        recipt = transpose3(crystal.A);
        sa = new double[n];
        wa = new double[n];

        this.parallelTeam = parallelTeam;
        sigmaARegion = new SigmaARegion(parallelTeam.getThreadCount());
    }

    /*
     * From sim and sim_integ functions in clipper utils:
     * http://www.ysbl.york.ac.uk/~cowtan/clipper/clipper.html and from ln_of_i0
     * and i1_over_i0 functions in bessel.h in scitbx module of cctbx:
     * http://cci.lbl.gov/cctbx_sources/scitbx/math/bessel.h
     */
    /**
     * <p>
     * sim</p>
     *
     * @param x a double.
     * @return a double.
     */
    public static double sim(double x) {
        if (x >= 0.0) {
            return (((x + sim_a) * x + sim_b) * x)
                    / (((x + sim_c) * x + sim_d) * x + sim_e);
        } else {
            return -(-(-(-x + sim_a) * x + sim_b) * x)
                    / (-(-(-x + sim_c) * x + sim_d) * x + sim_e);
        }
    }

    /**
     * <p>
     * sim_integ</p>
     *
     * @param x0 a double.
     * @return a double.
     */
    public static double sim_integ(double x0) {
        double x = abs(x0);
        double z = (x + sim_p) / sim_q;

        return sim_A * log(x + sim_g) + 0.5 * sim_B * log(z * z + 1.0)
                + sim_r * atan(z) + x + 1.0;
    }

    private class SigmaARegion extends ParallelRegion {

        boolean gradient = true;
        double modelK;
        double solventK;
        double solventUEq;
        double x[];
        double g[];
        private final double resm[][] = new double[3][3];
        private final double model_b[] = new double[6];
        private final double ustar[][] = new double[3][3];
        SharedInteger nsum;
        SharedInteger nsumr;
        SharedDouble sum;
        SharedDouble sumr;
        SharedDoubleArray grad;
        SigmaALoop sigmaALoop[];

        public SigmaARegion(int nThreads) {
            sigmaALoop = new SigmaALoop[nThreads];
            nsum = new SharedInteger();
            nsumr = new SharedInteger();
            sum = new SharedDouble();
            sumr = new SharedDouble();
        }

        public void init(double x[], double g[], boolean gradient) {
            this.x = x;
            this.g = g;
            this.gradient = gradient;
        }

        @Override
        public void finish() {
            if (gradient) {
                for (int i = 0; i < g.length; i++) {
                    g[i] = grad.get(i);
                }
            }
        }

        @Override
        public void start() {
            // Zero out the gradient
            if (gradient) {
                if (grad == null) {
                    grad = new SharedDoubleArray(g.length);
                }
                for (int i = 0; i < g.length; i++) {
                    grad.set(i, 0.0);
                }
            }
            sum.set(0.0);
            nsum.set(0);
            sumr.set(0.0);
            nsumr.set(0);

            modelK = refinementData.model_k;
            solventK = refinementData.solvent_k;
            solventUEq = refinementData.solvent_ueq;
            System.arraycopy(refinementData.model_b, 0, model_b, 0, 6);

            // Generate Ustar
            mat3symvec6(crystal.A, model_b, resm);
            mat3mat3(resm, recipt, ustar);

            for (int i = 0; i < n; i++) {
                sa[i] = 1.0 + x[i];
                wa[i] = x[n + i];
            }

            // Cheap method of preventing negative w values.
            for (int i = 0; i < n; i++) {
                if (wa[i] <= 0.0) {
                    wa[i] = 1.0e-6;
                }
            }
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (sigmaALoop[ti] == null) {
                sigmaALoop[ti] = new SigmaALoop();
            }

            try {
                execute(0, reflectionList.hkllist.size() - 1, sigmaALoop[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class SigmaALoop extends IntegerForLoop {

            /**
             * Thread local work variables.
             */
            private double lsum;
            private double lsumr;
            private int lnsum;
            private int lnsumr;
            private final double lgrad[];
            private final double resv[] = new double[3];
            private final double ihc[] = new double[3];
            private final ComplexNumber resc = new ComplexNumber();
            private final ComplexNumber fcc = new ComplexNumber();
            private final ComplexNumber fsc = new ComplexNumber();
            private final ComplexNumber fct = new ComplexNumber();
            private final ComplexNumber kfct = new ComplexNumber();
            private final ComplexNumber ecc = new ComplexNumber();
            private final ComplexNumber esc = new ComplexNumber();
            private final ComplexNumber ect = new ComplexNumber();
            private final ComplexNumber kect = new ComplexNumber();
            private final ComplexNumber mfo = new ComplexNumber();
            private final ComplexNumber mfo2 = new ComplexNumber();
            private final ComplexNumber dfcc = new ComplexNumber();
            private final ReflectionSpline spline = new ReflectionSpline(reflectionList, n);
            
            public SigmaALoop() {
                lgrad = new double[2*n];                 
            }
            
            @Override
            public void start() {
                lsum = 0.0;
                lsumr = 0.0;
                lnsum = 0;
                lnsumr = 0;
                Arrays.fill(lgrad, 0.0);
            }
            
            @Override
            public void finish() {
                sum.addAndGet(lsum);
                sumr.addAndGet(lsumr);
                nsum.addAndGet(lnsum);
                nsumr.addAndGet(lnsumr);
                for (int i=0; i< lgrad.length; i++) {
                    grad.getAndAdd(i, lgrad[i]);
                }
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int j = lb; j <= ub; j++) {
                    HKL ih = reflectionList.hkllist.get(j);
                    int i = ih.index();
                    // Constants
                    ihc[0] = ih.h();
                    ihc[1] = ih.k();
                    ihc[2] = ih.l();
                    vec3mat3(ihc, ustar, resv);
                    double u = modelK - dot(resv, ihc);
                    double s = Crystal.invressq(crystal, ih);
                    double ebs = exp(-twoPI2 * solventUEq * s);
                    double ksebs = solventK * ebs;
                    double kmems = exp(0.25 * u);
                    double km2 = exp(0.5 * u);
                    double epsc = ih.epsilonc();

                    // Spline setup
                    double ecscale = spline.f(s, refinementData.fcesq);
                    double eoscale = spline.f(s, refinementData.foesq);
                    double sai = spline.f(s, sa);
                    double wai = spline.f(s, wa);
                    double sa2 = pow(sai, 2.0);

                    // Structure factors
                    refinementData.get_fc_ip(i, fcc);
                    refinementData.get_fs_ip(i, fsc);
                    fct.copy(fcc);
                    if (refinementData.crs_fs.solventModel != SolventModel.NONE) {
                        resc.copy(fsc);
                        resc.times_ip(ksebs);
                        fct.plus_ip(resc);
                    }
                    kfct.copy(fct);
                    kfct.times_ip(kmems);

                    ecc.copy(fcc);
                    ecc.times_ip(sqrt(ecscale));
                    esc.copy(fsc);
                    esc.times_ip(sqrt(ecscale));
                    ect.copy(fct);
                    ect.times_ip(sqrt(ecscale));
                    kect.copy(kfct);
                    kect.times_ip(sqrt(ecscale));
                    double eo = fo[i][0] * sqrt(eoscale);
                    double sigeo = fo[i][1] * sqrt(eoscale);
                    double eo2 = pow(eo, 2.0);
                    double kect2 = pow(kect.abs(), 2.0);

                    // FOM
                    double d = 2.0 * sigeo * sigeo + epsc * wai;
                    double d2 = d * d;
                    double fomx = 2.0 * eo * sai * kect.abs() / d;

                    double inot, dinot, cf;
                    if (ih.centric()) {
                        inot = (abs(fomx) < 10.0) ? log(cosh(fomx)) : abs(fomx) + log(0.5);
                        dinot = tanh(fomx);
                        cf = 0.5;
                    } else {
                        inot = sim_integ(fomx);
                        dinot = sim(fomx);
                        cf = 1.0;
                    }
                    double llk = cf * log(d) + (eo2 + sa2 * kect2) / d - inot;

                    // Map coefficients
                    double f = dinot * eo;
                    double phi = kect.phase();
                    fomphi[i][0] = dinot;
                    fomphi[i][1] = phi;
                    mfo.re(f * cos(phi));
                    mfo.im(f * sin(phi));
                    mfo2.re(2.0 * f * cos(phi));
                    mfo2.im(2.0 * f * sin(phi));
                    dfcc.re(sai * kect.abs() * cos(phi));
                    dfcc.im(sai * kect.abs() * sin(phi));
                    // Set up map coefficients
                    fofc1[i][0] = 0.0;
                    fofc1[i][1] = 0.0;
                    fofc2[i][0] = 0.0;
                    fofc2[i][1] = 0.0;
                    dfc[i][0] = 0.0;
                    dfc[i][1] = 0.0;
                    dfs[i][0] = 0.0;
                    dfs[i][1] = 0.0;
                    if (Double.isNaN(fctot[i][0])) {
                        if (!Double.isNaN(fo[i][0])) {
                            fofc2[i][0] = mfo.re() / sqrt(eoscale);
                            fofc2[i][1] = mfo.im() / sqrt(eoscale);
                        }
                        continue;
                    }
                    if (Double.isNaN(fo[i][0])) {
                        if (!Double.isNaN(fctot[i][0])) {
                            fofc2[i][0] = dfcc.re() / sqrt(eoscale);
                            fofc2[i][1] = dfcc.im() / sqrt(eoscale);
                        }
                        continue;
                    }
                    // Update Fctot
                    fctot[i][0] = kfct.re();
                    fctot[i][1] = kfct.im();
                    // mFo - DFc
                    resc.copy(mfo);
                    resc.minus_ip(dfcc);
                    fofc1[i][0] = resc.re() / sqrt(eoscale);
                    fofc1[i][1] = resc.im() / sqrt(eoscale);
                    // 2mFo - DFc
                    resc.copy(mfo2);
                    resc.minus_ip(dfcc);
                    fofc2[i][0] = resc.re() / sqrt(eoscale);
                    fofc2[i][1] = resc.im() / sqrt(eoscale);

                    // Derivatives
                    double dfp1 = 2.0 * sa2 * km2 * ecscale;
                    double dfp2 = 2.0 * eo * sai * kmems * sqrt(ecscale);
                    double dfcr = (dfp1 * fct.re()) / d - ((dfp2 * fct.re()) / (d * fct.abs())) * dinot;
                    double dfci = (dfp1 * fct.im()) / d - ((dfp2 * fct.im()) / (d * fct.abs())) * dinot;
                    double dfsr = ((dfp2 * ksebs * fct.re()) / (d * fct.abs())) * dinot - (dfp1 * ksebs * fct.re()) / d;
                    double dfsi = ((dfp2 * ksebs * fct.im()) / (d * fct.abs())) * dinot - (dfp1 * ksebs * fct.im()) / d;
                    double dfsa = 2.0 * sai * kect2 / d - (2.0 * eo * kect.abs() / d) * dinot;
                    double dfwa = epsc * (cf / d - (eo2 + sa2 * kect2) / d2 + (2.0 * eo * sai * kect.abs() / d2) * dinot);

                    // Partial LLK wrt Fc or Fs
                    dfc[i][0] = dfcr * dfscale;
                    dfc[i][1] = dfci * dfscale;
                    dfs[i][0] = dfsr * dfscale;
                    dfs[i][1] = dfsi * dfscale;

                    // Only use freeR flagged reflections in overall sum
                    if (refinementData.isfreer(i)) {
                        lsum += llk;
                        lnsum++;
                    } else {
                        lsumr += llk;
                        lnsumr++;
                        dfsa = dfwa = 0.0;
                    }

                    if (gradient) {
                        int i0 = spline.i0();
                        int i1 = spline.i1();
                        int i2 = spline.i2();
                        double g0 = spline.dfi0();
                        double g1 = spline.dfi1();
                        double g2 = spline.dfi2();
                        // s derivative
                        lgrad[i0] += dfsa * g0;
                        lgrad[i1] += dfsa * g1;
                        lgrad[i2] += dfsa * g2;
                        // w derivative
                        lgrad[n + i0] += dfwa * g0;
                        lgrad[n + i1] += dfwa * g1;
                        lgrad[n + i2] += dfwa * g2;
          
                    }
                }
            }
        }
    }

    /**
     * <p>
     * target</p>
     *
     * @param x an array of double.
     * @param g an array of double.
     * @param gradient a boolean.
     * @param print a boolean.
     * @return a double.
     */
    public double target(double x[], double g[],
            boolean gradient, boolean print) {

        try {
            sigmaARegion.init(x, g, gradient);
            parallelTeam.execute(sigmaARegion);
        } catch (Exception e) {
            logger.info(e.toString());
        }

        double sum = sigmaARegion.sum.get();
        double sumr = sigmaARegion.sumr.get();
        refinementData.llkr = sumr;
        refinementData.llkf = sum;

        if (print) {
            int nsum = sigmaARegion.nsum.get();
            int nsumr = sigmaARegion.nsumr.get();
            StringBuilder sb = new StringBuilder("\n");
            sb.append(" sigmaA (s and w) fit using only Rfree reflections\n");
            sb.append(String.format("      # HKL: %10d (free set) %10d (working set) %10d (total)\n",
                    nsum, nsumr, nsum + nsumr));
            sb.append(String.format("   residual: %10g (free set) %10g (working set) %10g (total)\n",
                    sum, sumr, sum + sumr));
            sb.append(String.format("    x: "));
            for (int i = 0; i < x.length; i++) {
                sb.append(String.format("%g ", x[i]));
            }
            sb.append("\n    g: ");
            for (int i = 0; i < g.length; i++) {
                sb.append(String.format("%g ", g[i]));
            }
            sb.append("\n");
            logger.info(sb.toString());
        }

        totalEnergy = sum;
        return totalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double x[]) {
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        double sum = target(x, null, false, false);

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
            }
        }

        return sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double x[], double g[]) {
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        double sum = target(x, g, true, false);

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }

        return sum;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double[] scaling) {
        if (scaling != null && scaling.length == n * 2) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return null;
    }

    @Override
    public void setEnergyTermState(STATE state) {
    }
}
