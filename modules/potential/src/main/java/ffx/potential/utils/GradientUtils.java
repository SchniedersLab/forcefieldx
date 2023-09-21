// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.utils;

import ffx.numerics.Potential;
import ffx.potential.cli.GradientOptions;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import static ffx.utilities.StringUtils.parseAtomRanges;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * GradientUtils is a utility class for testing the gradient of a potential.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GradientUtils {

  private static final Logger logger = Logger.getLogger(GradientUtils.class.getName());

  private final Potential potential;

  /**
   * Constructor.
   *
   * @param potential The Potential to test.
   */
  public GradientUtils(Potential potential) {
    this.potential = potential;
  }

  /**
   * Test the gradient of the Potential.
   *
   * @param gradientOptions The GradientOptions to use.
   * @return The number of failures.
   */
  public int testGradient(GradientOptions gradientOptions) {
    int nFailures = 0;

    // Finite-difference step size in Angstroms.
    double step = gradientOptions.getDx();
    logger.info(" Finite-difference step size:\t" + step);

    // Print out the energy for each step.
    boolean print = gradientOptions.getVerbose();
    logger.info(" Verbose printing:\t\t" + print);

    // Collect degrees of freedom to test.
    List<Integer> degreesOfFreedomToTest;
    int nAtoms = potential.getNumberOfVariables() / 3;
    String gradientAtoms = gradientOptions.getGradientAtoms();
    if (gradientAtoms.equalsIgnoreCase("NONE")) {
      logger.info(" The gradient of no atoms was evaluated.");
      return nFailures;
    } else if (gradientAtoms.equalsIgnoreCase("ALL")) {
      logger.info(" Checking gradient for all degrees of freedom.\n");
      degreesOfFreedomToTest = new ArrayList<>();
      for (int i = 0; i < nAtoms; i++) {
        degreesOfFreedomToTest.add(i);
      }
    } else {
      degreesOfFreedomToTest = parseAtomRanges(" Gradient atoms", gradientAtoms, nAtoms);
      logger.info(" Checking gradient for degrees of freedom in the range: " + gradientAtoms + "\n");
    }

    // Collect analytic gradient.
    int n = potential.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    potential.getCoordinates(x);
    potential.energyAndGradient(x, g);

    // Upper bound for a typical atomic gradient.
    double expGrad = 1000.0;
    double gradientTolerance = gradientOptions.getTolerance();
    double width = 2.0 * step;
    double avLen = 0.0;
    double avGrad = 0.0;
    double expGrad2 = expGrad * expGrad;

    int nTested = 0;
    int index = 0;
    double[] numeric = new double[3];
    for (int k = 0; k < nAtoms; k++) {
      int i0 = index++;
      int i1 = index++;
      int i2 = index++;
      if (!degreesOfFreedomToTest.contains(k)) {
        continue;
      }
      nTested++;

      // Find numeric dX
      double orig = x[i0];
      x[i0] = x[i0] + step;
      double e = potential.energy(x);
      x[i0] = orig - step;
      e -= potential.energy(x);
      x[i0] = orig;
      numeric[0] = e / width;

      // Find numeric dY
      orig = x[i1];
      x[i1] = x[i1] + step;
      e = potential.energy(x);
      x[i1] = orig - step;
      e -= potential.energy(x);
      x[i1] = orig;
      numeric[1] = e / width;

      // Find numeric dZ
      orig = x[i2];
      x[i2] = x[i2] + step;
      e = potential.energy(x);
      x[i2] = orig - step;
      e -= potential.energy(x);
      x[i2] = orig;
      numeric[2] = e / width;

      double dx = g[i0] - numeric[0];
      double dy = g[i1] - numeric[1];
      double dz = g[i2] - numeric[2];
      double len = dx * dx + dy * dy + dz * dz;
      avLen += len;
      len = sqrt(len);

      double grad2 = g[i0] * g[i0] + g[i1] * g[i1] + g[i2] * g[i2];
      avGrad += grad2;

      if (len > gradientTolerance) {
        logger.info(format(" Degree of Freedom %d\n Failed: %10.6f\n", k + 1, len) +
            format(" Analytic: (%12.4f, %12.4f, %12.4f)\n", g[i0], g[i1], g[i2]) +
            format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
        ++nFailures;
      } else {
        logger.info(format(" Degree of Freedom %d\n Passed: %10.6f\n", k + 1, len) +
            format(" Analytic: (%12.4f, %12.4f, %12.4f)\n", g[i0], g[i1], g[i2]) +
            format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]));
      }

      if (grad2 > expGrad2) {
        logger.info(format(" Degree of Freedom %d has an unusually large gradient: %10.6f", k + 1, Math.sqrt(grad2)));
      }
      logger.info("\n");
    }

    avLen = avLen / nTested;
    avLen = sqrt(avLen);
    if (avLen > gradientTolerance) {
      logger.info(format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, gradientTolerance));
    } else {
      logger.info(format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, gradientTolerance));
    }
    logger.info(format(" Number of atoms failing analytic test: %d", nFailures));

    avGrad = avGrad / nTested;
    avGrad = sqrt(avGrad);
    if (avGrad > expGrad) {
      logger.info(format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad));
    } else {
      logger.info(format(" RMS gradient: %10.6f", avGrad));
    }

    return nFailures;
  }
}
