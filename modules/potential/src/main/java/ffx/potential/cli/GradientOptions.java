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
package ffx.potential.cli;

import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that test gradients.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GradientOptions {

  /**
   * The ArgGroup keeps the GradientOptions together when printing help.
   */
  @ArgGroup(heading = "%n Gradient Options%n", validate = false)
  public GradientOptionGroup group = new GradientOptionGroup();

  /**
   * -d or --dx Finite-difference step size.
   *
   * @return Returns the Finite-difference step size.
   */
  public double getDx() {
    return group.dx;
  }

  /**
   * --tol or --tolerance Gradient error tolerance (kcal/mol/Å).
   *
   * @return Returns the gradient error tolerance.
   */
  public double getTolerance() {
    return group.tolerance;
  }

  /**
   * --ga or --gradientAtoms Ranges of atoms to test [ALL, NONE, Range(s): 1-3,6-N].
   *
   * @return Returns the ranges of atoms to test.
   */
  public String getGradientAtoms() {
    return group.gradientAtoms;
  }

  /**
   * -v or --verbose is a flag to print out energy at each step.
   *
   * @return Returns true to print out energy at each step.
   */
  public boolean getVerbose() {
    return group.verbose;
  }

  /**
   * Collection of Gradient Options.
   */
  private static class GradientOptionGroup {

    /** -d or --dx Finite-difference step size. */
    @Option(
        names = {"-d", "--dx"},
        defaultValue = "1.0e-5",
        paramLabel = "1.0e-5 Å",
        description = "Finite-difference step size.")
    public double dx;

    /** --tol or --tolerance Gradient error tolerance (kcal/mol/Å). */
    @Option(
        names = {"--tol", "--tolerance"},
        defaultValue = "1.0e-3",
        paramLabel = "1.0e-3 kcal/mol/Å",
        description = "Gradient error tolerance.")
    public double tolerance;

    /** --ga or --gradientAtoms Ranges of atoms to test [ALL, NONE, Range(s): 1-3,6-N]. */
    @Option(
        names = {"--ga", "--gradientAtoms"},
        paramLabel = "ALL",
        defaultValue = "ALL",
        description = "Ranges of atoms to test [ALL, NONE, Range(s): 1-3,6-N].")
    public String gradientAtoms;

    /** -v or --verbose is a flag to print out energy at each step. */
    @Option(
        names = {"-v", "--verbose"},
        paramLabel = "false",
        description = "Print out the energy for each step.")
    public boolean verbose = false;
  }

}
