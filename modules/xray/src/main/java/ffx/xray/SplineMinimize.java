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

import static java.lang.System.arraycopy;

import ffx.algorithms.Terminatable;
import ffx.crystal.ReflectionList;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch.LineSearchResult;
import ffx.numerics.optimization.OptimizationListener;
import ffx.xray.SplineEnergy.Type;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * SplineMinimize class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class SplineMinimize implements OptimizationListener, Terminatable {

  private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());

  private final SplineEnergy splineEnergy;
  private final int n;
  private final double[] x;
  private final double[] grad;
  private final double[] scaling;
  private boolean done = false;
  private boolean terminate = false;
  private double grms;
  private int nSteps;

  /**
   * Constructor for SplineMinimize.
   *
   * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
   * @param refinementData a {@link ffx.xray.DiffractionRefinementData} object.
   * @param x an array of double.
   * @param type a int.
   */
  public SplineMinimize(
      ReflectionList reflectionList,
      DiffractionRefinementData refinementData,
      double[] x,
      int type) {
    this.x = x;
    n = x.length;
    splineEnergy = new SplineEnergy(reflectionList, refinementData, n, type);
    grad = new double[n];
    scaling = new double[n];
    for (int i = 0; i < n; i++) {
      if (type == Type.FOTOESQ || type == Type.FCTOESQ) {
        x[i] = 0.1;
      } else {
        x[i] = 1.0;
      }
      scaling[i] = 1.0;
    }
  }

  /**
   * getCoordinates.
   *
   * @param x an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  public double[] getCoordinates(double x[]) {
    if (x == null) {
      x = new double[this.x.length];
    }
    arraycopy(this.x, 0, x, 0, this.x.length);
    return x;
  }

  /**
   * getNumberOfVariables.
   *
   * @return a int.
   */
  public int getNumberOfVariables() {
    return x.length;
  }

  /**
   * minimize
   *
   * @return a {@link ffx.xray.SplineEnergy} object.
   */
  public SplineEnergy minimize() {
    return minimize(0.5);
  }

  /**
   * minimize
   *
   * @param eps a double.
   * @return a {@link ffx.xray.SplineEnergy} object.
   */
  public SplineEnergy minimize(double eps) {
    return minimize(5, eps);
  }

  /**
   * minimize
   *
   * @param m a int.
   * @param eps a double.
   * @return a {@link ffx.xray.SplineEnergy} object.
   */
  public SplineEnergy minimize(int m, double eps) {

    splineEnergy.setScaling(scaling);

    double e = splineEnergy.energyAndGradient(x, grad);

    done = false;
    int status = LBFGS.minimize(n, m, x, e, grad, eps, splineEnergy, this);
    done = true;
    switch (status) {
      case 0:
        logger.fine(String.format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
        break;
      case 1:
        logger.fine(String.format("\n Optimization terminated at step %d.\n", nSteps));
        break;
      default:
        logger.warning("\n Spline Optimization failed.\n");
    }

    splineEnergy.setScaling(null);

    return splineEnergy;
  }

  /** {@inheritDoc} */
  @Override
  public boolean optimizationUpdate(
      int iter,
      int nfun,
      double grms,
      double xrms,
      double f,
      double df,
      double angle,
      LineSearchResult info) {
    this.grms = grms;
    this.nSteps = iter;
    if (terminate) {
      logger.info(" The optimization received a termination request.");
      // Tell the L-BFGS optimizer to terminate.
      return false;
    }
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public void terminate() {
    terminate = true;
    while (!done) {
      synchronized (this) {
        try {
          wait(1);
        } catch (Exception e) {
          logger.log(Level.WARNING, "Exception terminating minimization.\n", e);
        }
      }
    }
  }

  /**
   * Getter for the field <code>splineEnergy</code>.
   *
   * @return a {@link ffx.xray.SplineEnergy} object.
   */
  SplineEnergy getSplineEnergy() {
    return splineEnergy;
  }
}
