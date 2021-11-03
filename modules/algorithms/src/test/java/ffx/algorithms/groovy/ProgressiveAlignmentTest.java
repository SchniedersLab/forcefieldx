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
package ffx.algorithms.groovy;

import static org.junit.Assert.assertEquals;

import ffx.algorithms.misc.AlgorithmsTest;
import org.junit.Test;

/**
 * Tests test.CrystalSuperpose command to verify crystals are being aligned and compared.
 *
 * @author Aaron J. Nessler
 */

public class ProgressiveAlignmentTest extends AlgorithmsTest {

  private final double tolerance = 0.001;

  // TODO: add more tests with more parameters

  /** Tests the CrystalSuperpose script with default settings.. */
  @Test
  public void testBaseProgressiveAlignment() {

    // Set-up the input arguments for the CrystalSuperpose script.
    String[] args = {"src/main/java/ffx/algorithms/structures/C23.arc"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the CrystalSuperpose script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.0851007, superposeCrystals.runningStatistics.getMean(), tolerance);
  }

  /** Tests the CrystalSuperpose script without hydrogens. */
  @Test
  public void testBaseProgressiveAlignmentNoHydrogen() {

    // Set-up the input arguments for the CrystalSuperpose script.
    String[] args = {"--nh", "src/main/java/ffx/algorithms/structures/C23.arc"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the CrystalSuperpose script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.0810508, superposeCrystals.runningStatistics.getMean(), tolerance);
  }

  /** Tests the CrystalSuperpose script in Sohncke group. */
  @Test
  public void testBaseProgressiveAlignmentSohnckeGroup() {

    // Set-up the input arguments for the CrystalSuperpose script.
    String[] args = {"--nh", "src/main/java/ffx/algorithms/structures/dap.xyz",
        "src/main/java/ffx/algorithms/structures/dap.xyz_close"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the CrystalSuperpose script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(0.066147, superposeCrystals.runningStatistics.getMean(), tolerance);
  }

  /** Tests the CrystalSuperpose script on tricky handedness case. */
  @Test
  public void testBaseProgressiveAlignmentHandedness() {

    // Set-up the input arguments for the CrystalSuperpose script.
    String[] args = {"--nh", "src/main/java/ffx/algorithms/structures/dap2.xyz",
        "src/main/java/ffx/algorithms/structures/dap2.xyz_2"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the CrystalSuperpose script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(0.114188, superposeCrystals.runningStatistics.getMean(), tolerance);
  }

  /** Tests the CrystalSuperpose script on asymmetric unit greater than one. */
  @Test
  public void testBaseProgressiveAlignmentZPrime2() {

    // Set-up the input arguments for the CrystalSuperpose script.
    String[] args = {"--nh", "src/main/java/ffx/algorithms/structures/XAFPAY02.arc"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the CrystalSuperpose script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.15354157, superposeCrystals.runningStatistics.getMean(), tolerance);
  }

  @Test
  public void testProgressiveAlignmentHelp() {
    // Set-up the input arguments for the CrystalSuperpose script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the CrystalSuperpose script.
    SuperposeCrystals SuperposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = SuperposeCrystals;
  }
}
