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
package ffx.potential.groovy;

import static org.junit.Assert.assertEquals;

import ffx.potential.groovy.test.LambdaGradient;
import ffx.potential.utils.PotentialTest;
import org.junit.Test;

/**
 * Tests test.LambdaGradient command to determine that the end state potentials and derivatives are
 * correct.
 *
 * @author Aaron J. Nessler
 * @author Mallory R. Tollefson
 */
public class LambdaGradientTest extends PotentialTest {

  /** Tests the End States of the LambdaGradient class. */
  @Test
  public void testLambdaGradient() {
    // Prepare the Binding with input arguments.
    String[] args = {"--ac", "1-26", "src/main/java/ffx/potential/structures/phenacetin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    LambdaGradient lambdaGradient = new LambdaGradient(binding).run();
    potentialScript = lambdaGradient;

    double expectedPotentialEnergyVac = -7.70864934;
    double expectedPotentialEnergyXtal = -36.22707541638729;
    double actualPotentialEnergyVac = lambdaGradient.e0;
    double actualPotentialEnergyXtal = lambdaGradient.e1;

    assertEquals(expectedPotentialEnergyVac, actualPotentialEnergyVac, 1E-6);
    assertEquals(expectedPotentialEnergyXtal, actualPotentialEnergyXtal, 1E-6);
    assertEquals(0, lambdaGradient.ndEdLFailures);
    assertEquals(0, lambdaGradient.nd2EdL2Failures);
    assertEquals(0, lambdaGradient.ndEdXdLFailures);
    assertEquals(0, lambdaGradient.ndEdXFailures);
  }

  @Test
  public void testLambdaGradientHelp() {
    // Set-up the input arguments for the LambdaGradient script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    LambdaGradient lambdaGradient = new LambdaGradient(binding).run();
    potentialScript = lambdaGradient;
  }

  /** Tests the End States of the LambdaGradient class when softcore is active. */
  @Test
  public void testLambdaGradientIntermolecularSoftcore() {
    // Set-up the input arguments for the LambdaGradient script.
    String[] args = {"--ac", "1-44", "src/main/java/ffx/potential/structures/ethylparaben.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    LambdaGradient lambdaGradient = new LambdaGradient(binding).run();
    potentialScript = lambdaGradient;

    double expectedPotentialEnergyVac = -6.67890842;
    double expectedPotentialEnergyXtal = -57.410945848200896;
    double actualPotentialEnergyVac = lambdaGradient.e0;
    double actualPotentialEnergyXtal = lambdaGradient.e1;

    assertEquals(expectedPotentialEnergyVac, actualPotentialEnergyVac, 1E-6);
    assertEquals(expectedPotentialEnergyXtal, actualPotentialEnergyXtal, 1E-6);
    assertEquals(0, lambdaGradient.ndEdLFailures);
    assertEquals(0, lambdaGradient.nd2EdL2Failures);
    assertEquals(0, lambdaGradient.ndEdXdLFailures);
    assertEquals(0, lambdaGradient.ndEdXFailures);
  }

  /** Tests the End States of the LambdaGradient class when mapped via sym op. */
  @Test
  public void testLambdaGradientSymOp() {
    // Set-up the input arguments for the LambdaGradient script.
    String[] args = {"--sf", "TRIG", "--ls", "-l", "0.5", "--ac", "3,30","--ac2","3,30",
        "src/main/java/ffx/potential/structures/roy02_P1.xyz", "src/main/java/ffx/potential/structures/roy31.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    LambdaGradient lambdaGradient = new LambdaGradient(binding).run();
    potentialScript = lambdaGradient;

    double expectedPotentialEnergyVac = -33.519321973744496;
    double expectedPotentialEnergyXtal = -41.768837155671626;
    double actualPotentialEnergyVac = lambdaGradient.e0;
    double actualPotentialEnergyXtal = lambdaGradient.e1;

    assertEquals(expectedPotentialEnergyVac, actualPotentialEnergyVac, 1E-6);
    assertEquals(expectedPotentialEnergyXtal, actualPotentialEnergyXtal, 1E-6);
    assertEquals(0, lambdaGradient.ndEdLFailures);
    assertEquals(0, lambdaGradient.nd2EdL2Failures);
    assertEquals(0, lambdaGradient.ndEdXdLFailures);
    assertEquals(0, lambdaGradient.ndEdXFailures);
  }
}
