// ******************************************************************************
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
// ******************************************************************************
package ffx.xray;

import static ffx.crystal.Crystal.invressq;
import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.MatrixMath.mat3Mat3;
import static ffx.numerics.math.MatrixMath.mat3SymVec6;
import static ffx.numerics.math.MatrixMath.transpose3;
import static ffx.numerics.math.MatrixMath.vec3Mat3;
import static ffx.numerics.special.ModifiedBessel.i1OverI0;
import static ffx.numerics.special.ModifiedBessel.lnI0;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
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
import ffx.numerics.Potential;
import ffx.numerics.math.ComplexNumber;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import java.util.logging.Logger;

/**
 * Optimize SigmaA coefficients (using spline coefficients) and structure factor derivatives using a
 * likelihood target function.
 *
 * <p>This target can also be used for structure refinement.
 *
 * @author Timothy D. Fenn<br>
 * @see <a href="http://dx.doi.org/10.1107/S0021889804031474" target="_blank"> K. Cowtan, J. Appl.
 *     Cryst. (2005). 38, 193-198</a>
 * @see <a href="http://dx.doi.org/10.1107/S0907444992007352" target="_blank"> A. T. Brunger, Acta
 *     Cryst. (1993). D49, 24-36</a>
 * @see <a href="http://dx.doi.org/10.1107/S0907444996012255" target="_blank"> G. N. Murshudov, A.
 *     A. Vagin and E. J. Dodson, Acta Cryst. (1997). D53, 240-255</a>
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank"> A. T. Brunger, Acta
 *     Cryst. (1989). A45, 42-50.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0108767386099622" target="_blank"> R. J. Read, Acta
 *     Cryst. (1986). A42, 140-149.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0108767396004370" target="_blank"> N. S. Pannu and R. J.
 *     Read, Acta Cryst. (1996). A52, 659-668.</a>
 * @since 1.0
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
  private final DiffractionRefinementData refinementData;
  private final ParallelTeam parallelTeam;
  private final Crystal crystal;
  private final double[][] fSigF;
  private final double[][] fcTot;
  private final double[][] fomPhi;
  private final double[][] foFc1;
  private final double[][] foFc2;
  private final double[][] dFc;
  private final double[][] dFs;
  private final int nBins;

  /** Crystal Volume^2 / (2 * number of grid points) */
  private final double dfScale;
  /** Transpose of the matrix 'A' that converts from Cartesian to fractional coordinates. */
  private final double[][] transposeA;

  private final double[] sa;
  private final double[] wa;

  private final SigmaARegion sigmaARegion;
  private final boolean useCernBessel;
  private double[] optimizationScaling = null;
  private double totalEnergy;
  private STATE state = STATE.BOTH;

  /**
   * Constructor for SigmaAEnergy.
   *
   * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
   * @param refinementData a {@link ffx.xray.DiffractionRefinementData} object.
   * @param parallelTeam the ParallelTeam to execute the SigmaAEnergy.
   */
  SigmaAEnergy(
      ReflectionList reflectionList,
      DiffractionRefinementData refinementData,
      ParallelTeam parallelTeam) {
    this.reflectionList = reflectionList;
    this.refinementData = refinementData;
    this.parallelTeam = parallelTeam;
    this.crystal = reflectionList.crystal;
    this.fSigF = refinementData.fSigF;
    this.fcTot = refinementData.fcTot;
    this.fomPhi = refinementData.fomPhi;
    this.foFc1 = refinementData.foFc1;
    this.foFc2 = refinementData.foFc2;
    this.dFc = refinementData.dFc;
    this.dFs = refinementData.dFs;
    this.nBins = refinementData.nBins;

    // Initialize parameters.
    assert (refinementData.crystalReciprocalSpaceFc != null);
    double nGrid2 =
        2.0
            * refinementData.crystalReciprocalSpaceFc.getXDim()
            * refinementData.crystalReciprocalSpaceFc.getYDim()
            * refinementData.crystalReciprocalSpaceFc.getZDim();
    dfScale = (crystal.volume * crystal.volume) / nGrid2;
    transposeA = transpose3(crystal.A);
    sa = new double[nBins];
    wa = new double[nBins];

    sigmaARegion = new SigmaARegion(this.parallelTeam.getThreadCount());
    useCernBessel = true;
  }

  /**
   * From sim and sim_integ functions in clipper utils:
   * http://www.ysbl.york.ac.uk/~cowtan/clipper/clipper.html and from lnI0 and i1OverI0 functions in
   * bessel.h in scitbx module of cctbx: http://cci.lbl.gov/cctbx_sources/scitbx/math/bessel.h
   *
   * @param x a double.
   * @return a double.
   */
  private static double sim(double x) {
    if (x >= 0.0) {
      return (((x + sim_a) * x + sim_b) * x) / (((x + sim_c) * x + sim_d) * x + sim_e);
    } else {
      return -(-(-(-x + sim_a) * x + sim_b) * x) / (-(-(-x + sim_c) * x + sim_d) * x + sim_e);
    }
  }

  /**
   * sim_integ
   *
   * @param x0 a double.
   * @return a double.
   */
  private static double sim_integ(double x0) {
    double x = abs(x0);
    double z = (x + sim_p) / sim_q;
    return sim_A * log(x + sim_g) + 0.5 * sim_B * log(z * z + 1.0) + sim_r * atan(z) + x + 1.0;
  }

  /** {@inheritDoc} */
  @Override
  public boolean destroy() {
    // Should be destroyed upstream in DiffractionData.
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public double energy(double[] x) {
    unscaleCoordinates(x);
    double sum = target(x, null, false, false);
    scaleCoordinates(x);
    return sum;
  }

  /** {@inheritDoc} */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    unscaleCoordinates(x);
    double sum = target(x, g, true, false);
    scaleCoordinatesAndGradient(x, g);
    return sum;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public double[] getCoordinates(double[] parameters) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public STATE getEnergyTermState() {
    return state;
  }

  /** {@inheritDoc} */
  @Override
  public void setEnergyTermState(STATE state) {
    this.state = state;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getMass() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public int getNumberOfVariables() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public double[] getScaling() {
    return optimizationScaling;
  }

  /** {@inheritDoc} */
  @Override
  public void setScaling(double[] scaling) {
    if (scaling != null && scaling.length == nBins * 2) {
      optimizationScaling = scaling;
    } else {
      optimizationScaling = null;
    }
  }

  /** {@inheritDoc} */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Return a reference to each variables type.
   */
  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getVelocity(double[] velocity) {

    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public void setAcceleration(double[] acceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /** {@inheritDoc} */
  @Override
  public void setVelocity(double[] velocity) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * target
   *
   * @param x an array of double.
   * @param g an array of double.
   * @param gradient a boolean.
   * @param print a boolean.
   * @return a double.
   */
  public double target(double[] x, double[] g, boolean gradient, boolean print) {

    try {
      sigmaARegion.init(x, g, gradient);
      parallelTeam.execute(sigmaARegion);
    } catch (Exception e) {
      logger.info(e.toString());
    }

    double sum = sigmaARegion.sum.get();
    double sumR = sigmaARegion.sumR.get();
    refinementData.llkR = sumR;
    refinementData.llkF = sum;

    if (print) {
      int nSum = sigmaARegion.nSum.get();
      int nSumr = sigmaARegion.nSumR.get();
      StringBuilder sb = new StringBuilder("\n");
      sb.append(" sigmaA (s and w) fit using only R free reflections\n");
      sb.append(
          format(
              "      # HKL: %10d (free set) %10d (working set) %10d (total)\n",
              nSum, nSumr, nSum + nSumr));
      sb.append(
          format(
              "   residual: %10g (free set) %10g (working set) %10g (total)\n",
              sum, sumR, sum + sumR));
      sb.append("    x: ");
      for (double x1 : x) {
        sb.append(format("%g ", x1));
      }
      sb.append("\n    g: ");
      for (double v : g) {
        sb.append(format("%g ", v));
      }
      sb.append("\n");
      logger.info(sb.toString());
    }

    totalEnergy = sum;
    return totalEnergy;
  }

  private class SigmaARegion extends ParallelRegion {

    private final double[][] resm = new double[3][3];
    private final double[] model_b = new double[6];
    private final double[][] ustar = new double[3][3];
    boolean gradient = true;
    double modelK;
    double solventK;
    double solventUEq;
    double[] x;
    double[] g;
    SharedInteger nSum;
    SharedInteger nSumR;
    SharedDouble sum;
    SharedDouble sumR;
    SharedDoubleArray grad;
    SigmaALoop[] sigmaALoop;

    SigmaARegion(int nThreads) {
      sigmaALoop = new SigmaALoop[nThreads];
      nSum = new SharedInteger();
      nSumR = new SharedInteger();
      sum = new SharedDouble();
      sumR = new SharedDouble();
    }

    @Override
    public void finish() {
      if (gradient) {
        for (int i = 0; i < g.length; i++) {
          g[i] = grad.get(i);
        }
      }
    }

    public void init(double[] x, double[] g, boolean gradient) {
      this.x = x;
      this.g = g;
      this.gradient = gradient;
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
      nSum.set(0);
      sumR.set(0.0);
      nSumR.set(0);

      modelK = refinementData.modelScaleK;
      solventK = refinementData.bulkSolventK;
      solventUEq = refinementData.bulkSolventUeq;
      arraycopy(refinementData.modelAnisoB, 0, model_b, 0, 6);

      // Generate Ustar
      mat3SymVec6(crystal.A, model_b, resm);
      mat3Mat3(resm, transposeA, ustar);

      for (int i = 0; i < nBins; i++) {
        sa[i] = 1.0 + x[i];
        wa[i] = x[nBins + i];
      }

      // Cheap method of preventing negative w values.
      for (int i = 0; i < nBins; i++) {
        if (wa[i] <= 0.0) {
          wa[i] = 1.0e-6;
        }
      }
    }

    private class SigmaALoop extends IntegerForLoop {

      private final double[] lGrad;
      private final double[] resv = new double[3];
      private final double[] ihc = new double[3];
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
      private final ReflectionSpline spline = new ReflectionSpline(reflectionList, nBins);
      // Thread local work variables.
      private double lSum;
      private double lSumR;
      private int lSumN;
      private int lSumRN;

      SigmaALoop() {
        lGrad = new double[2 * nBins];
      }

      @Override
      public void finish() {
        sum.addAndGet(lSum);
        sumR.addAndGet(lSumR);
        nSum.addAndGet(lSumN);
        nSumR.addAndGet(lSumRN);
        for (int i = 0; i < lGrad.length; i++) {
          grad.getAndAdd(i, lGrad[i]);
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
          double s = invressq(crystal, ih);
          double ebs = exp(-twoPI2 * solventUEq * s);
          double ksebs = solventK * ebs;
          double kmems = exp(0.25 * u);
          double km2 = exp(0.5 * u);
          double epsc = ih.epsilonc();

          // Spline setup
          double ecscale = spline.f(s, refinementData.esqFc);
          double eoscale = spline.f(s, refinementData.esqFo);
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
          if (refinementData.crystalReciprocalSpaceFs.solventModel != SolventModel.NONE) {
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
          double eo = fSigF[i][0] * sqrtEOScale;
          double sigeo = fSigF[i][1] * sqrtEOScale;
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
            if (useCernBessel) {
              inot = lnI0(fomx);
              dinot = i1OverI0(fomx);
            } else {
              inot = sim_integ(fomx);
              dinot = sim(fomx);
            }
            cf = 1.0;
          }
          double llk = cf * log(d) + (eo2 + sa2 * kect2) * id - inot;

          // Map coefficients
          double f = dinot * eo;
          double phi = kect.phase();
          double sinPhi = sin(phi);
          double cosPhi = cos(phi);
          fomPhi[i][0] = dinot;
          fomPhi[i][1] = phi;
          mfo.re(f * cosPhi);
          mfo.im(f * sinPhi);
          mfo2.re(2.0 * f * cosPhi);
          mfo2.im(2.0 * f * sinPhi);
          akect = kect.abs();
          dfcc.re(sai * akect * cosPhi);
          dfcc.im(sai * akect * sinPhi);
          // Set up map coefficients
          foFc1[i][0] = 0.0;
          foFc1[i][1] = 0.0;
          foFc2[i][0] = 0.0;
          foFc2[i][1] = 0.0;
          dFc[i][0] = 0.0;
          dFc[i][1] = 0.0;
          dFs[i][0] = 0.0;
          dFs[i][1] = 0.0;
          if (isNaN(fcTot[i][0])) {
            if (!isNaN(fSigF[i][0])) {
              foFc2[i][0] = mfo.re() * iSqrtEOScale;
              foFc2[i][1] = mfo.im() * iSqrtEOScale;
            }
            continue;
          }
          if (isNaN(fSigF[i][0])) {
            if (!isNaN(fcTot[i][0])) {
              foFc2[i][0] = dfcc.re() * iSqrtEOScale;
              foFc2[i][1] = dfcc.im() * iSqrtEOScale;
            }
            continue;
          }
          // Update Fctot
          fcTot[i][0] = kfct.re();
          fcTot[i][1] = kfct.im();
          // mFo - DFc
          resc.copy(mfo);
          resc.minusIP(dfcc);
          foFc1[i][0] = resc.re() * iSqrtEOScale;
          foFc1[i][1] = resc.im() * iSqrtEOScale;
          // 2mFo - DFc
          resc.copy(mfo2);
          resc.minusIP(dfcc);
          foFc2[i][0] = resc.re() * iSqrtEOScale;
          foFc2[i][1] = resc.im() * iSqrtEOScale;

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
          double dfwa =
              epsc * (cf * id - (eo2 + sa2 * kect2) * id2 + 2.0 * eo * sai * akect * id2 * dinot);

          // Partial LLK wrt Fc or Fs
          dFc[i][0] = dfcr * dfScale;
          dFc[i][1] = dfci * dfScale;
          dFs[i][0] = dfsr * dfScale;
          dFs[i][1] = dfsi * dfScale;

          // Only use free R flagged reflections in overall sum
          if (refinementData.isFreeR(i)) {
            lSum += llk;
            lSumN++;
          } else {
            lSumR += llk;
            lSumRN++;
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
            lGrad[i0] += dfsa * g0;
            lGrad[i1] += dfsa * g1;
            lGrad[i2] += dfsa * g2;
            // w derivative
            lGrad[nBins + i0] += dfwa * g0;
            lGrad[nBins + i1] += dfwa * g1;
            lGrad[nBins + i2] += dfwa * g2;
          }
        }
      }

      @Override
      public void start() {
        lSum = 0.0;
        lSumR = 0.0;
        lSumN = 0;
        lSumRN = 0;
        fill(lGrad, 0.0);
      }
    }
  }
}
