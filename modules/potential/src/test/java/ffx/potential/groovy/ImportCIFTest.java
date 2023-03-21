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
package ffx.potential.groovy;

import ffx.potential.utils.PotentialTest;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests ImportCIF command to determine that files are being translated correctly.
 *
 * @author Aaron J. Nessler
 */

public class ImportCIFTest extends PotentialTest {

  /**
   * Test a basic CIF to XYZ conversion.
   */
  @Test
  public void testImportCIF() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"src/main/java/ffx/potential/structures/CBZ16.cif",
        "src/main/java/ffx/potential/structures/cbz.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.

    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".XYZ"));
  }

  /**
   * Test writing out a CIF file (XYZ to CIF).
   */
  @Test
  public void testImportCIFWriteAsCIF() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"--sc", "src/main/java/ffx/potential/structures/paracetamol.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".CIF"));
  }

  /**
   * Test ImportCIF when the CIF file is missing hydrogen atoms.
   */
  @Test
  public void testImportCIFNoHydrogen() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"--fl","src/main/java/ffx/potential/structures/CBZ03.cif",
        "src/main/java/ffx/potential/structures/cbz.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".XYZ"));
  }

  /**
   * Test ImportCIF when several molecules are included in asymmetric unit.
   */
  @Test
  public void testImportCIFMultipleMolecules() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"src/main/java/ffx/potential/structures/1183240.cif",
            "src/main/java/ffx/potential/structures/asplyswat.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".XYZ"));
  }

  /**
   * Test ImportCIF when similar (but different) molecules are in asymmetric unit.
   */
  @Test
  public void testImportCIFzPrimeChallenge() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"--fl","src/main/java/ffx/potential/structures/1183241.cif",
            "src/main/java/ffx/potential/structures/glulys.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".XYZ"));
  }

  /**
   * Test ImportCIF on concatenated CIF files (produces multiple ARC files).
   */
  @Test
  public void testImportCIFarc() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"--fl","src/main/java/ffx/potential/structures/cbzs.cif",
            "src/main/java/ffx/potential/structures/cbz.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(3, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".ARC"));
  }

  /**
   * Test ImportCIF across a molecular disulfide bond.
   */
  @Test
  public void testImportCIFdisulfide() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"src/main/java/ffx/potential/structures/UFAGIS01.cif",
            "src/main/java/ffx/potential/structures/uf.xyz"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".XYZ"));
  }

  /**
   * Test ImportCIF conversion to PDB.
   */
  @Test
  public void testImportCIFpdb() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"src/main/java/ffx/potential/structures/288726.cif",
        "src/main/java/ffx/potential/structures/ALA-HIE.pdb"};
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the ImportCIF script.
    ImportCIF ImportCIF = new ImportCIF(binding).run();
    potentialScript = ImportCIF;

    assertEquals(1, ImportCIF.createdFiles.length);
    assertTrue(ImportCIF.createdFiles[0].toUpperCase().contains(".PDB"));
  }

  /**
   * Print out help message.
   */
  @Test
  public void testImportCIFHelp() {
    // Set up the input arguments for the ImportCIF script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the ImportCIF script.
    potentialScript = new ImportCIF(binding).run();
  }
}