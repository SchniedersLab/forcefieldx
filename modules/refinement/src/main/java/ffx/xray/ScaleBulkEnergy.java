// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.numerics.OptimizationInterface;
import ffx.numerics.math.ComplexNumber;
import ffx.xray.solvent.SolventModel;

import javax.annotation.Nullable;
import java.util.logging.Logger;

import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.MatrixMath.mat3Mat3Multiply;
import static ffx.numerics.math.MatrixMath.mat3SymVec6;
import static ffx.numerics.math.MatrixMath.mat3Transpose;
import static ffx.numerics.math.MatrixMath.vec3Mat3;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * Fit bulk solvent and aniso B scaling terms to correct calculated structure factors against data
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444905007894" target="_blank"> P. V. Afonine, R. W.
 * Grosse-Kunstleve and P. D. Adams, Acta Cryst. (2005). D61, 850-855</a>
 * @see <a href="http://dx.doi.org/10.1107/S0021889802008580" target="_blank"> R. W.
 * Grosse-Kunstleve and P. D. Adams, J. Appl. Cryst. (2002). 35, 477-480.</a>
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A. Grant, B. T. Pickup, A.
 * Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J. S. Jiang, A. T.
 * Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class ScaleBulkEnergy implements OptimizationInterface {

  private static final Logger logger = Logger.getLogger(ScaleBulkEnergy.class.getName());
  private static final double twopi2 = 2.0 * PI * PI;
  private static final double[] v000 = {0.0, 0.0, 0.0};
  private static final double[] v100 = {1.0, 0.0, 0.0};
  private static final double[] v010 = {0.0, 1.0, 0.0};
  private static final double[] v001 = {0.0, 0.0, 1.0};
  private static final double[][] u11 = {v100, v000, v000};
  private static final double[][] u22 = {v000, v010, v000};
  private static final double[][] u33 = {v000, v000, v001};
  private static final double[][] u12 = {v010, v100, v000};
  private static final double[][] u13 = {v001, v000, v100};
  private static final double[][] u23 = {v000, v001, v010};
  private final double[][] recipt;
  private final double[][] j11;
  private final double[][] j22;
  private final double[][] j33;
  private final double[][] j12;
  private final double[][] j13;
  private final double[][] j23;

  private final ReflectionList reflectionList;
  private final Crystal crystal;
  private final DiffractionRefinementData refinementData;
  private final double[][] fc;
  private final double[][] fcTot;
  private final double[][] fSigF;
  private final int n;
  private final int solventN;
  private final ParallelTeam parallelTeam;
  private final ScaleBulkEnergyRegion scaleBulkEnergyRegion;
  private double[] optimizationScaling = null;
  private double totalEnergy;
  private double R;
  private double Rfree;

  /**
   * Constructor for ScaleBulkEnergy.
   *
   * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
   * @param refinementData a {@link ffx.xray.DiffractionRefinementData} object.
   * @param n              an int.
   * @param parallelTeam   the ParallelTeam to execute the ScaleBulkEnergy.
   */
  ScaleBulkEnergy(ReflectionList reflectionList, DiffractionRefinementData refinementData,
      int n, ParallelTeam parallelTeam) {
    this.reflectionList = reflectionList;
    this.crystal = reflectionList.crystal;
    this.refinementData = refinementData;
    this.fc = refinementData.fc;
    this.fcTot = refinementData.fcTot;
    this.fSigF = refinementData.fSigF;
    this.n = n;
    this.solventN = n - refinementData.nScale;

    recipt = mat3Transpose(crystal.A);
    j11 = mat3Mat3Multiply(mat3Mat3Multiply(crystal.A, u11), recipt);
    j22 = mat3Mat3Multiply(mat3Mat3Multiply(crystal.A, u22), recipt);
    j33 = mat3Mat3Multiply(mat3Mat3Multiply(crystal.A, u33), recipt);
    j12 = mat3Mat3Multiply(mat3Mat3Multiply(crystal.A, u12), recipt);
    j13 = mat3Mat3Multiply(mat3Mat3Multiply(crystal.A, u13), recipt);
    j23 = mat3Mat3Multiply(mat3Mat3Multiply(crystal.A, u23), recipt);

    int threadCount = parallelTeam.getThreadCount();
    this.parallelTeam = parallelTeam;
    scaleBulkEnergyRegion = new ScaleBulkEnergyRegion(threadCount);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    // The parallelTeam should have been passed in by DiffractionData, which handles destroying it.
    return true;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    unscaleCoordinates(x);
    double sum = target(x, null, false, false);
    scaleCoordinates(x);
    return sum;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    unscaleCoordinates(x);
    double sum = target(x, g, true, false);
    scaleCoordinatesAndGradient(x, g);
    return sum;
  }

  public double getR() {
    return R;
  }
  public double getRfree() {
    return Rfree;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] parameters) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCoordinates(double[] parameters) {
    throw new UnsupportedOperationException("Not supported yet.");
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
  public double[] getScaling() {
    return optimizationScaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(@Nullable double[] scaling) {
    if (scaling != null && scaling.length == n) {
      optimizationScaling = scaling;
    } else {
      optimizationScaling = null;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * ScaleBulk target energy.
   *
   * @param x        The parameter array.
   * @param g        The gradient of the parameter array.
   * @param gradient If true, compute the gradient.
   * @param print    If true, enable verbose printing.
   * @return The target energy.
   */
  public double target(double[] x, double[] g, boolean gradient, boolean print) {
    try {
      scaleBulkEnergyRegion.init(x, g, gradient);
      parallelTeam.execute(scaleBulkEnergyRegion);
    } catch (Exception e) {
      logger.severe(e.toString());
    }

    double sum = scaleBulkEnergyRegion.sum.get();
    double sumfo = scaleBulkEnergyRegion.sumFo.get();
    double r = scaleBulkEnergyRegion.r.get();
    double rf = scaleBulkEnergyRegion.rf.get();
    double rfree = scaleBulkEnergyRegion.rFree.get();
    double rfreef = scaleBulkEnergyRegion.rFreeF.get();

    R = (r / rf) * 100.0;
    Rfree = (rfree / rfreef) * 100.0;
    totalEnergy = sum / sumfo;

    if (gradient) {
      double isumfo = 1.0 / sumfo;
      for (int i = 0; i < g.length; i++) {
        g[i] *= isumfo;
      }
    }

    if (print) {
      StringBuilder sb = new StringBuilder("\n");
      sb.append(" Bulk solvent and scale fit\n");
      sb.append(format("  Residual:  %10.5f  R:  %8.3f  Rfree:  %8.3f\n", totalEnergy, R, Rfree));
      sb.append("  Params:   ");
      for (double x1 : x) {
        sb.append(format("%16.8f ", x1));
      }
      if (gradient) {
        sb.append("\n  Gradient: ");
        for (double v : g) {
          sb.append(format("%16.8f ", v));
        }
      }
      sb.append("\n");
      logger.info(sb.toString());
    }

    return totalEnergy;
  }

  private class ScaleBulkEnergyRegion extends ParallelRegion {

    private final double[] modelB = new double[6];
    private final double[][] uStar = new double[3][3];
    private final double[][] resM = new double[3][3];
    boolean gradient = true;
    double[] x;
    double[] g;
    double solventK;
    double modelK;
    double solventUEq;
    SharedDouble r;
    SharedDouble rf;
    SharedDouble rFree;
    SharedDouble rFreeF;
    SharedDouble sum;
    SharedDouble sumFo;
    SharedDoubleArray grad;
    ScaleBulkEnergyLoop[] scaleBulkEnergyLoop;

    ScaleBulkEnergyRegion(int nThreads) {
      scaleBulkEnergyLoop = new ScaleBulkEnergyLoop[nThreads];
      r = new SharedDouble();
      rf = new SharedDouble();
      rFree = new SharedDouble();
      rFreeF = new SharedDouble();
      sum = new SharedDouble();
      sumFo = new SharedDouble();
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
    public void run() {
      int ti = getThreadIndex();
      if (scaleBulkEnergyLoop[ti] == null) {
        scaleBulkEnergyLoop[ti] = new ScaleBulkEnergyLoop();
      }

      try {
        execute(0, reflectionList.hklList.size() - 1, scaleBulkEnergyLoop[ti]);
      } catch (Exception e) {
        logger.info(e.toString());
      }
    }

    @Override
    public void start() {
      r.set(0.0);
      rf.set(0.0);
      rFree.set(0.0);
      rFreeF.set(0.0);
      sum.set(0.0);
      sumFo.set(0.0);

      for (int i = 0; i < 6; i++) {
        if (crystal.scaleB[i] >= 0) {
          modelB[i] = x[solventN + crystal.scaleB[i]];
        }
      }

      modelK = x[0];
      solventK = refinementData.bulkSolventK;
      solventUEq = refinementData.bulkSolventUeq;
      if (solventN > 1) {
        solventK = x[1];
        solventUEq = x[2];
      }

      // Generate Ustar.
      mat3SymVec6(crystal.A, modelB, resM);
      mat3Mat3Multiply(resM, recipt, uStar);

      if (gradient) {
        if (grad == null) {
          grad = new SharedDoubleArray(g.length);
        }
        for (int i = 0; i < g.length; i++) {
          grad.set(i, 0.0);
        }
      }
    }

    private class ScaleBulkEnergyLoop extends IntegerForLoop {

      private final double[] resv = new double[3];
      private final double[] ihc = new double[3];
      private final ComplexNumber resc = new ComplexNumber();
      private final ComplexNumber fcc = new ComplexNumber();
      private final ComplexNumber fsc = new ComplexNumber();
      private final ComplexNumber fct = new ComplexNumber();
      private final ComplexNumber kfct = new ComplexNumber();
      private final double[] lgrad;
      private double lr;
      private double lrf;
      private double lrfree;
      private double lrfreef;
      private double lsum;
      private double lsumfo;

      ScaleBulkEnergyLoop() {
        lgrad = new double[n];
      }

      @Override
      public void finish() {
        r.addAndGet(lr);
        rf.addAndGet(lrf);
        rFree.addAndGet(lrfree);
        rFreeF.addAndGet(lrfreef);
        sum.addAndGet(lsum);
        sumFo.addAndGet(lsumfo);
        if (gradient) {
          for (int i = 0; i < lgrad.length; i++) {
            grad.getAndAdd(i, lgrad[i]);
          }
        }
      }

      @Override
      public void run(int lb, int ub) {

        for (int j = lb; j <= ub; j++) {
          HKL ih = reflectionList.hklList.get(j);
          int i = ih.getIndex();
          if (isNaN(fc[i][0]) || isNaN(fSigF[i][0]) || fSigF[i][1] <= 0.0) {
            continue;
          }

          // Constants
          double s = crystal.invressq(ih);
          ihc[0] = ih.getH();
          ihc[1] = ih.getK();
          ihc[2] = ih.getL();
          vec3Mat3(ihc, uStar, resv);
          double u = modelK - dot(resv, ihc);
          double expBS = exp(-twopi2 * solventUEq * s);
          double ksExpBS = solventK * expBS;
          double expU = exp(0.25 * u);

          // Structure Factors
          refinementData.getFcIP(i, fcc);
          refinementData.getFsIP(i, fsc);
          fct.copy(fcc);
          if (refinementData.crystalReciprocalSpaceFs.getSolventModel() != SolventModel.NONE) {
            resc.copy(fsc);
            resc.timesIP(ksExpBS);
            fct.plusIP(resc);
          }
          kfct.copy(fct);
          kfct.timesIP(expU);

          // Total structure factor (for refinement)
          fcTot[i][0] = kfct.re();
          fcTot[i][1] = kfct.im();

          // Target
          double f1 = refinementData.getF(i);
          double akfct = kfct.abs();
          double af1 = abs(f1);
          double d = f1 - akfct;
          double d2 = d * d;
          double dr = -2.0 * d;

          lsum += d2;
          lsumfo += f1 * f1;

          if (refinementData.isFreeR(i)) {
            lrfree += abs(af1 - abs(akfct));
            lrfreef += af1;
          } else {
            lr += abs(af1 - abs(akfct));
            lrf += af1;
          }

          if (gradient) {
            // modelK/modelB - common derivative element
            double dfm = 0.25 * akfct * dr;
            // bulk solvent - common derivative element
            double afsc = fsc.abs();
            double dfb = expBS * (fcc.re() * fsc.re() + fcc.im() * fsc.im() + ksExpBS * afsc * afsc);

            // modelK derivative
            lgrad[0] += dfm;
            if (solventN > 1) {
              double iafct = 1.0 / fct.abs();
              // solventK derivative
              lgrad[1] += expU * dfb * dr * iafct;
              // solventUeq derivative
              lgrad[2] += expU * -twopi2 * s * solventK * dfb * dr * iafct;
            }

            for (int jj = 0; jj < 6; jj++) {
              if (crystal.scaleB[jj] >= 0) {
                switch (jj) {
                  case (0) -> {
                    // B11
                    vec3Mat3(ihc, j11, resv);
                    lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                  }
                  case (1) -> {
                    // B22
                    vec3Mat3(ihc, j22, resv);
                    lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                  }
                  case (2) -> {
                    // B33
                    vec3Mat3(ihc, j33, resv);
                    lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                  }
                  case (3) -> {
                    // B12
                    vec3Mat3(ihc, j12, resv);
                    lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                  }
                  case (4) -> {
                    // B13
                    vec3Mat3(ihc, j13, resv);
                    lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                  }
                  case (5) -> {
                    // B23
                    vec3Mat3(ihc, j23, resv);
                    lgrad[solventN + crystal.scaleB[jj]] += -dfm * dot(resv, ihc);
                  }
                }
              }
            }
          }
        }
      }

      @Override
      public void start() {
        lr = 0.0;
        lrf = 0.0;
        lrfree = 0.0;
        lrfreef = 0.0;
        lsum = 0.0;
        lsumfo = 0.0;
        fill(lgrad, 0.0);
      }
    }
  }
}
