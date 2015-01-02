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

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.atan;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.cosh;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.tanh;

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
import static ffx.numerics.VectorMath.mat3Mat3;
import static ffx.numerics.VectorMath.mat3SymVec6;
import static ffx.numerics.VectorMath.transpose3;
import static ffx.numerics.VectorMath.vec3Mat3;

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
     * @param parallelTeam the ParallelTeam to execute the SigmaAEnergy.
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

        // parallelTeam = new ParallelTeam(1);
        this.parallelTeam = parallelTeam;
        sigmaARegion = new SigmaARegion(this.parallelTeam.getThreadCount());
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
    
    /**
     * Adapted from CERN's cern.jet.math.tdouble package as included with ParallelColt:
     * http://sourceforge.net/projects/parallelcolt/
     */
    private static class Bessel {
        /**
         * Analagous to Fenn's sim(x).
         */
        public static double ln_of_i0(double x) {
            return Math.log(i0(x));
        }
        /**
         * Analagous to Fenn's sim_integ(x).
         */
        public static double i1_over_i0(double x) {
            return (i1(x) / i0(x));
        }
        
        /**
         * Modified zero-order Bessel function.
         * The function is defined as i0(x) = j0( ix ).
         * The range is partitioned into the two intervals [0,8] and (8, infinity).
         * Chebyshev polynomial expansions are employed in each interval.
         * 
         * @param x the value to compute the bessel function of.
         * @return the result
         */
        public static double i0(double x) {
            double y;
            if (x < 0)
                x = -x;
            if (x <= 8.0) {
                y = (x / 2.0) - 2.0;
                return (Math.exp(x) * chbevl(y, A_i0, 30));
            }

            return (Math.exp(x) * chbevl(32.0 / x - 2.0, B_i0, 25) / Math.sqrt(x));
        }
        
        /**
         * Modified 1st-order Bessel function.
         * The function is defined as i1(x) = -i j1( ix ).
         * The range is partitioned into the two intervals [0,8] and (8, infinity).
         * Chebyshev polynomial expansions are employed in each interval.
         * 
         * @param x the value to compute the bessel function of.
         * @return the result
         */
        public static double i1(double x) {
            double y, z;

            z = Math.abs(x);
            if (z <= 8.0) {
                y = (z / 2.0) - 2.0;
                z = chbevl(y, A_i1, 29) * z * Math.exp(z);
            } else {
                z = Math.exp(z) * chbevl(32.0 / z - 2.0, B_i1, 25) / Math.sqrt(z);
            }
            if (x < 0.0)
                z = -z;
            return (z);
        }
        
        /**
         * Chebyshev coefficients for exp(-x) I0(x) in the interval [0,8].
         * lim(x->0){ exp(-x) I0(x) } = 1.
         */
        private static final double[] A_i0 = { -4.41534164647933937950E-18, 3.33079451882223809783E-17,
                -2.43127984654795469359E-16, 1.71539128555513303061E-15, -1.16853328779934516808E-14,
                7.67618549860493561688E-14, -4.85644678311192946090E-13, 2.95505266312963983461E-12,
                -1.72682629144155570723E-11, 9.67580903537323691224E-11, -5.18979560163526290666E-10,
                2.65982372468238665035E-9, -1.30002500998624804212E-8, 6.04699502254191894932E-8,
                -2.67079385394061173391E-7, 1.11738753912010371815E-6, -4.41673835845875056359E-6,
                1.64484480707288970893E-5, -5.75419501008210370398E-5, 1.88502885095841655729E-4,
                -5.76375574538582365885E-4, 1.63947561694133579842E-3, -4.32430999505057594430E-3,
                1.05464603945949983183E-2, -2.37374148058994688156E-2, 4.93052842396707084878E-2,
                -9.49010970480476444210E-2, 1.71620901522208775349E-1, -3.04682672343198398683E-1,
                6.76795274409476084995E-1 };
        
        /**
         * Chebyshev coefficients for exp(-x) sqrt(x) I0(x) in the inverted interval [8,infinity].
         * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
         */
        private static final double[] B_i0 = { -7.23318048787475395456E-18, -4.83050448594418207126E-18,
                4.46562142029675999901E-17, 3.46122286769746109310E-17, -2.82762398051658348494E-16,
                -3.42548561967721913462E-16, 1.77256013305652638360E-15, 3.81168066935262242075E-15,
                -9.55484669882830764870E-15, -4.15056934728722208663E-14, 1.54008621752140982691E-14,
                3.85277838274214270114E-13, 7.18012445138366623367E-13, -1.79417853150680611778E-12,
                -1.32158118404477131188E-11, -3.14991652796324136454E-11, 1.18891471078464383424E-11,
                4.94060238822496958910E-10, 3.39623202570838634515E-9, 2.26666899049817806459E-8,
                2.04891858946906374183E-7, 2.89137052083475648297E-6, 6.88975834691682398426E-5, 3.36911647825569408990E-3,
                8.04490411014108831608E-1 };
        
        /**
         * Chebyshev coefficients for exp(-x) I1(x) / x in the interval [0,8].
         * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
         */
        private static final double[] A_i1 = { 2.77791411276104639959E-18, -2.11142121435816608115E-17,
                1.55363195773620046921E-16, -1.10559694773538630805E-15, 7.60068429473540693410E-15,
                -5.04218550472791168711E-14, 3.22379336594557470981E-13, -1.98397439776494371520E-12,
                1.17361862988909016308E-11, -6.66348972350202774223E-11, 3.62559028155211703701E-10,
                -1.88724975172282928790E-9, 9.38153738649577178388E-9, -4.44505912879632808065E-8,
                2.00329475355213526229E-7, -8.56872026469545474066E-7, 3.47025130813767847674E-6,
                -1.32731636560394358279E-5, 4.78156510755005422638E-5, -1.61760815825896745588E-4,
                5.12285956168575772895E-4, -1.51357245063125314899E-3, 4.15642294431288815669E-3,
                -1.05640848946261981558E-2, 2.47264490306265168283E-2, -5.29459812080949914269E-2,
                1.02643658689847095384E-1, -1.76416518357834055153E-1, 2.52587186443633654823E-1 };
        
        /*
         * Chebyshev coefficients for exp(-x) sqrt(x) I1(x) in the inverted interval [8,infinity].
         * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
         */
        private static final double[] B_i1 = { 7.51729631084210481353E-18, 4.41434832307170791151E-18,
                -4.65030536848935832153E-17, -3.20952592199342395980E-17, 2.96262899764595013876E-16,
                3.30820231092092828324E-16, -1.88035477551078244854E-15, -3.81440307243700780478E-15,
                1.04202769841288027642E-14, 4.27244001671195135429E-14, -2.10154184277266431302E-14,
                -4.08355111109219731823E-13, -7.19855177624590851209E-13, 2.03562854414708950722E-12,
                1.41258074366137813316E-11, 3.25260358301548823856E-11, -1.89749581235054123450E-11,
                -5.58974346219658380687E-10, -3.83538038596423702205E-9, -2.63146884688951950684E-8,
                -2.51223623787020892529E-7, -3.88256480887769039346E-6, -1.10588938762623716291E-4,
                -9.76109749136146840777E-3, 7.78576235018280120474E-1 };
        
        /**
         * Evaluates Chebyshev polynomials (first kind) at x/2.
         * NOTE: Argument x must first be transformed to the interval (-1,1).
         * NOTE: Coefficients are in reverse; zero-order term is last.
         * 
         * @param x argument to the polynomial.
         * @param coef the coefficients of the polynomial.
         * @param N the number of coefficients.
         * @return the result
         */
        public static double chbevl(double x, double coef[], int N) {
            double b0, b1, b2;

            int p = 0;
            int i;

            b0 = coef[p++];
            b1 = 0.0;
            i = N - 1;

            do {
                b2 = b1;
                b1 = b0;
                b0 = x * b1 - b2 + coef[p++];
            } while (--i > 0);

            return (0.5 * (b0 - b2));
        }
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

            // logger.info(String.format("SigmaA Sum   %16.8f", sum.get()));
            // logger.info(String.format("SigmaA Sum R %16.8f", sumr.get()));
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
            mat3SymVec6(crystal.A, model_b, resm);
            mat3Mat3(resm, recipt, ustar);

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
                lgrad = new double[2 * n];
            }

            @Override
            public void start() {
                lsum = 0.0;
                lsumr = 0.0;
                lnsum = 0;
                lnsumr = 0;
                fill(lgrad, 0.0);
            }

            @Override
            public void finish() {
                sum.addAndGet(lsum);
                sumr.addAndGet(lsumr);
                nsum.addAndGet(lnsum);
                nsumr.addAndGet(lnsumr);
                for (int i = 0; i < lgrad.length; i++) {
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
                    vec3Mat3(ihc, ustar, resv);
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
                    double sqrtECScale = sqrt(ecscale);
                    double sqrtEOScale = sqrt(eoscale);
                    double iSqrtEOScale = 1.0 / sqrtEOScale;

                    double sai = spline.f(s, sa);
                    double wai = spline.f(s, wa);
                    double sa2 = sai * sai;

                    // Structure factors
                    refinementData.getFcIP(i, fcc);
                    refinementData.getFsIP(i, fsc);
                    fct.copy(fcc);
                    if (refinementData.crs_fs.solventModel != SolventModel.NONE) {
                        resc.copy(fsc);
                        resc.timesIP(ksebs);
                        fct.plusIP(resc);
                    }
                    kfct.copy(fct);
                    kfct.timesIP(kmems);

                    ecc.copy(fcc);
                    ecc.timesIP(sqrtECScale);
                    esc.copy(fsc);
                    esc.timesIP(sqrtECScale);
                    ect.copy(fct);
                    ect.timesIP(sqrtECScale);
                    kect.copy(kfct);
                    kect.timesIP(sqrtECScale);
                    double eo = fo[i][0] * sqrtEOScale;
                    double sigeo = fo[i][1] * sqrtEOScale;
                    double eo2 = eo * eo;
                    double akect = kect.abs();
                    double kect2 = akect * akect;

                    // FOM
                    double d = 2.0 * sigeo * sigeo + epsc * wai;
                    double id = 1.0 / d;
                    double id2 = id * id;
                    double fomx = 2.0 * eo * sai * kect.abs() * id;

                    double inot, dinot, cf;
                    if (ih.centric()) {
                        inot = (abs(fomx) < 10.0) ? log(cosh(fomx)) : abs(fomx) + log(0.5);
                        dinot = tanh(fomx);
                        cf = 0.5;
                    } else {
//                        inot = sim_integ(fomx);
//                        dinot = sim(fomx);
                        inot = Bessel.i1_over_i0(fomx);
                        dinot = Bessel.ln_of_i0(fomx);
                        cf = 1.0;
                    }
                    double llk = cf * log(d) + (eo2 + sa2 * kect2) * id - inot;

                    // Map coefficients
                    double f = dinot * eo;
                    double phi = kect.phase();
                    double sinPhi = sin(phi);
                    double cosPhi = cos(phi);
                    fomphi[i][0] = dinot;
                    fomphi[i][1] = phi;
                    mfo.re(f * cosPhi);
                    mfo.im(f * sinPhi);
                    mfo2.re(2.0 * f * cosPhi);
                    mfo2.im(2.0 * f * sinPhi);
                    akect = kect.abs();
                    dfcc.re(sai * akect * cosPhi);
                    dfcc.im(sai * akect * sinPhi);
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
                            fofc2[i][0] = mfo.re() * iSqrtEOScale;
                            fofc2[i][1] = mfo.im() * iSqrtEOScale;
                        }
                        continue;
                    }
                    if (Double.isNaN(fo[i][0])) {
                        if (!Double.isNaN(fctot[i][0])) {
                            fofc2[i][0] = dfcc.re() * iSqrtEOScale;
                            fofc2[i][1] = dfcc.im() * iSqrtEOScale;
                        }
                        continue;
                    }
                    // Update Fctot
                    fctot[i][0] = kfct.re();
                    fctot[i][1] = kfct.im();
                    // mFo - DFc
                    resc.copy(mfo);
                    resc.minusIP(dfcc);
                    fofc1[i][0] = resc.re() * iSqrtEOScale;
                    fofc1[i][1] = resc.im() * iSqrtEOScale;
                    // 2mFo - DFc
                    resc.copy(mfo2);
                    resc.minusIP(dfcc);
                    fofc2[i][0] = resc.re() * iSqrtEOScale;
                    fofc2[i][1] = resc.im() * iSqrtEOScale;

                    // Derivatives
                    double dafct = d * fct.abs();
                    double idafct = 1.0 / dafct;
                    double dfp1 = 2.0 * sa2 * km2 * ecscale;
                    double dfp2 = 2.0 * eo * sai * kmems * sqrt(ecscale);
                    double dfp1id = dfp1 * id;
                    double dfp2id = dfp2 * idafct * dinot;
                    double dfp12 = dfp1id - dfp2id;
                    double dfp21 = ksebs * (dfp2id - dfp1id);
                    double dfcr = fct.re() * dfp12;
                    double dfci = fct.im() * dfp12;
                    double dfsr = fct.re() * dfp21;
                    double dfsi = fct.im() * dfp21;
                    double dfsa = 2.0 * (sai * kect2 - eo * akect * dinot) * id;
                    double dfwa = epsc * (cf * id - (eo2 + sa2 * kect2) * id2
                            + 2.0 * eo * sai * akect * id2 * dinot);

                    // Partial LLK wrt Fc or Fs
                    dfc[i][0] = dfcr * dfscale;
                    dfc[i][1] = dfci * dfscale;
                    dfs[i][0] = dfsr * dfscale;
                    dfs[i][1] = dfsi * dfscale;

                    // Only use freeR flagged reflections in overall sum
                    if (refinementData.isFreeR(i)) {
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
