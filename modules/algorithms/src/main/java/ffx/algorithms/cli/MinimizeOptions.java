// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.algorithms.cli;

import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that involve local energy minimization.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MinimizeOptions {

  /**
   * The ArgGroup keeps the MinimizationOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Minimization Options%n", validate = false)
  public MinimizeOptionGroup group = new MinimizeOptionGroup();

  /**
   * Convergence criteria.
   *
   * @return a double.
   */
  public double getEps() {
    return group.eps;
  }

  public void setEps(double eps) {
    group.eps = eps;
  }

  /**
   * Number of minimization steps.
   *
   * @return a int.
   */
  public int getIterations() {
    return group.iterations;
  }

  public void setIterations(int iterations) {
    group.iterations = iterations;
  }

  /**
   * Collection of Minimize Options.
   */
  private static class MinimizeOptionGroup {

    /** -i or --iterations Number of minimization steps. */
    @Option(
        names = {"-I", "--iterations"},
        paramLabel = "Unlimited",
        // Integer.MAX_VALUE = 2^31 -1 = 2147483647.
        defaultValue = "2147483647",
        description = "Number of minimization steps.")
    private int iterations;

    /** -e or --eps Convergence criteria. */
    @Option(
        names = {"-e", "--eps"},
        paramLabel = "1.0",
        defaultValue = "1.0",
        description = "Convergence criteria.")
    private double eps;
  }
}
