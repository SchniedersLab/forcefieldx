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
package ffx.algorithms.groovy;

import ffx.algorithms.groovy.test.OSTGradient;
import ffx.algorithms.misc.AlgorithmsTest;
import org.apache.commons.math3.util.FastMath;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.Arrays;
import java.util.Collection;

import static java.lang.Math.random;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.junit.Assert.assertEquals;

/** @author Michael J. Schnieders */
@RunWith(Parameterized.class)
public class OSTGradientTest extends AlgorithmsTest {

  private String info;
  private String filename;
  private int nAtoms;
  private double tolerance = 1.0e-2;

  public OSTGradientTest(String info, String filename, int nAtoms) {
    this.info = info;
    this.filename = filename;
    this.nAtoms = nAtoms;
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "C23 OST Test", // info
                "compound23.xyz", // filename
                43 // Number of atoms.
            }
        });
  }

  @Test
  public void testOSTBiasHelp() {
    // Set-up the input arguments for the script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    OSTGradient ostGradient = new OSTGradient(binding).run();
    algorithmsScript = ostGradient;
  }

  @Test
  public void testOSTGradient() {

    double lambda = random();
    int atomID = (int) floor(FastMath.random() * nAtoms) + 1;

    // Set-up the input arguments for the script.
    String[] args = {
        "--ac", "ALL",
        "-l", Double.toString(lambda),
        "--ga", Integer.toString(atomID),
        getResourcePath(filename)
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    OSTGradient ostGradient = new OSTGradient(binding).run();
    algorithmsScript = ostGradient;

    double dUdLError = ostGradient.dUdLError;
    double nFailures = ostGradient.nFailures;

    // Assert that energy is conserved at the end of the dynamics trajectory.
    assertEquals(info + ": dUdL error: ", 0.0, dUdLError, tolerance);
    assertEquals(info + ": Number of coordinate gradient errors: ", 0, nFailures, 0);
  }

  @Test
  public void testMetaDynamicsGradient() {

    double lambda = random();
    int atomID = (int) floor(FastMath.random() * nAtoms) + 1;

    // Set-up the input arguments for the script.
    String[] args = {
        "--ac", "ALL",
        "--meta",
        "-l", Double.toString(lambda),
        "--ga", Integer.toString(atomID),
        getResourcePath(filename)
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    OSTGradient ostGradient = new OSTGradient(binding).run();
    algorithmsScript = ostGradient;

    double dUdLError = ostGradient.dUdLError;
    double nFailures = ostGradient.nFailures;

    // Assert that energy is conserved at the end of the dynamics trajectory.
    assertEquals(info + ": dUdL error: ", 0.0, dUdLError, tolerance);
    assertEquals(info + ": Number of coordinate gradient errors: ", 0, nFailures, 0);
  }
}
