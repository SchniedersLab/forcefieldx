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
package ffx.algorithms.commands;

import static org.junit.Assert.assertEquals;

import ffx.algorithms.misc.AlgorithmsTest;
import org.junit.Test;

/**
 * Tests SuperposeCrystals command to verify crystals are being aligned and compared.
 *
 * @author Aaron J. Nessler
 */
public class SuperposeCrystalsTest extends AlgorithmsTest {

  private final double TOLERANCE = 0.005;

  /**
   * Tests the SuperposeCrystals script with default settings.
   */
  @Test
  public void testBaseSingleSuperposeCrystals() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "0", getResourcePath("C23.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.087248318, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.226387, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script with RMSD size of one.
   */
  @Test
  public void testBaseSingleAUSuperposeCrystals() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"--na", "1", "-l", "0", getResourcePath("C23.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.055486, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.161867, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script with small symmetric molecules (can be tricky to determine conformations).
   */
  @Test
  public void testBaseSingleMTSuperposeCrystals() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "0",
        getResourcePath("ace.arc_s"),
        getResourcePath("ace.arc_sx")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(100, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 100 comparisons (10 by 10).
    assertEquals(0.199630, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.199630, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script with small symmetric molecules (can be tricky to determine conformations).
   */
  @Test
  public void testBaseSingleMTSuperposeCrystals2() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {getResourcePath("ace.arc"), getResourcePath("ace.arc_x")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(4, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 100 comparisons (10 by 10).
    assertEquals(0.224877, superposeCrystals.runningStatistics.getMin(), TOLERANCE);
    assertEquals(0.224877, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.224877, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script with default settings and average linkage.
   */
  @Test
  public void testBaseAverageSuperposeCrystals() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "1", getResourcePath("C23.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.0830895, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.217224949, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script with hydrogens.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsHydrogen() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"--ih", "-l", "1", getResourcePath("C23.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.111014, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.297522, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script with hydrogens.
   */
  @Test
  public void testBaseSingleSuperposeCrystalsHydrogen() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"--ih", "-l", "0", getResourcePath("C23.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Only off-diagonal values.
    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.1125556, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.30190515, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script in Sohncke group.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsSohnckeGroup() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "1", getResourcePath("dap.xyz"), getResourcePath("dap.xyz_close")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(0.066569, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script in Sohncke group.
   */
  @Test
  public void testBaseSingleSuperposeCrystalsSohnckeGroup() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "0", getResourcePath("dap.xyz"), getResourcePath("dap.xyz_close")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(0.073036, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on tricky handedness case.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsHandedness() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "1", getResourcePath("dap2.xyz"), getResourcePath("dap2.xyz_2")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(0.112274, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on tricky handedness case.
   */
  @Test
  public void testBaseSingleSuperposeCrystalsHandedness() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "0", getResourcePath("dap2.xyz"), getResourcePath("dap2.xyz_2")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(0.098483, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on asymmetric unit greater than one.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsZPrime2() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "1", getResourcePath("XAFPAY02.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.151675, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.375328, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on asymmetric unit greater than one.
   */
  @Test
  public void testBaseSingleSuperposeCrystalsZPrime2() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "0", getResourcePath("XAFPAY02.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.142474326, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.355770508, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on protein crystal.
   */
  @Test
  public void testBaseSingleSuperposeCrystalsSmallProtein() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "0", getResourcePath("2olx.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.478175, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(1.383652, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on protein crystal.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsSmallProtein() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-l", "1", getResourcePath("2olx.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(6, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.441851, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(1.269418376, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script print symmetry operator functionality.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsPrintSym() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"--na", "1", "--ps", "1.0", getResourcePath("CBZ01_P1.xyz"), getResourcePath("cbz11.xyz")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(1, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.10510735, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.10510735, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script alchemical atoms functionality.
   */
  @Test
  public void testBaseAverageSuperposeCrystalsAlchemicalAtoms() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"--na", "1", "--ps", "1.0", "--ac", "4,7,9,11,13,15,17,19,21,23,25,27",
        getResourcePath("ZEYBIO01.xyz_amoebax"),
        getResourcePath("ZEYBIO04_zero.xyz_amoebax")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(1, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.
    assertEquals(0.14180274, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.14180274, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Tests the SuperposeCrystals script on long molecule (potential issue with inflation factor).
   */
  @Test
  public void testSuperposeCrystalsLongMolecule() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {getResourcePath("dameda.arc")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    assertEquals(3, superposeCrystals.runningStatistics.getCount());
    // Mean RMSD for 6 comparisons.

    assertEquals(0.200420, superposeCrystals.runningStatistics.getMean(), TOLERANCE);
    assertEquals(0.601259, superposeCrystals.runningStatistics.getMax(), TOLERANCE);
  }

  /**
   * Test the automatic symmetry operator generation.
   */
  @Test
  public void testAutoSymGeneration() {

    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"--na", "1", "--as","--ih", "--mw", getResourcePath("2olx.xyz_xtal"), getResourcePath("2onx.xyz_xtal")};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals superposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = superposeCrystals;

    // Mean RMSD for 1 comparison.
    assertEquals(1, superposeCrystals.runningStatistics.getCount());

    // Only checks last of final groups RMSD_1 and could change if final group order is different... ultimately not helpful.
    assertEquals(0.028228, superposeCrystals.runningStatistics.getMean(), TOLERANCE);

    // Check generated symmetry operators are the correct length. This ONLY checks lengths... better than nothing?
    assertEquals(1736, superposeCrystals.symOpsA.length());
    assertEquals(1737, superposeCrystals.symOpsB.length());
  }

  @Test
  public void testSuperposeCrystalsHelp() {
    // Set up the input arguments for the SuperposeCrystals script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the SuperposeCrystals script.
    SuperposeCrystals SuperposeCrystals = new SuperposeCrystals(binding).run();
    algorithmsScript = SuperposeCrystals;
  }
}
