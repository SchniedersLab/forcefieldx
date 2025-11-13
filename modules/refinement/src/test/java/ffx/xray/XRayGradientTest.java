// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import ffx.algorithms.misc.AlgorithmsTest;
import ffx.xray.commands.test.Gradient;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Tests X-Ray Gradient.
 *
 * @author Michael J. Schnieders
 */
public class XRayGradientTest extends AlgorithmsTest {

  @Test
  public void testAlametXYZGradient() {
    // Set-up the input arguments.
    String[] args = {
        "-m", "coordinates",
        "--sol", "none",
        "--aRadBuffer", "3.0",
        "-G", "0.5",
        "--tol", "1.0e-2",
        "--params", "1-2",
        getResourcePath("alamet.pdb_2"),
        getResourcePath("alamet.mtz")
    };
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // alamet.pdb was refined with CNS
    // alamet.pdb_2 was refined with FFX

    // Construct and evaluate the Gradient script.
    Gradient gradient = new Gradient(binding).run();
    algorithmsScript = gradient;
    assertEquals(" Number of failed gradient components.", 0, gradient.nFailures);
  }

  @Test
  public void test1N7SXYZGradient() {
    // Set-up the input arguments.
    String[] args = {
        "-m", "coordinates",
        "--aRadBuffer", "4.0",
        "--params", "1-2",
        getResourcePath("1N7S.pdb"),
        getResourcePath("1N7S.cif")
    };
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the Gradient script.
    Gradient gradient = new Gradient(binding).run();
    algorithmsScript = gradient;
    assertEquals(" Number of failed gradient components.", 0, gradient.nFailures);
  }

  @Test
  public void test1N7SBFactorsGradient() {
    // Set-up the input arguments.
    String[] args = {
        "-m", "bfactors",
        "--aRadBuffer", "4.0",
        "--params", "1-2",
        getResourcePath("1N7S.pdb"),
        getResourcePath("1N7S.cif")
    };
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    // Construct and evaluate the Gradient script.
    Gradient gradient = new Gradient(binding).run();
    algorithmsScript = gradient;
    assertEquals(" Number of failed gradient components.", 0, gradient.nFailures);
  }

  @Test
  public void testGradientHelp() {
    // Set-up the input arguments for the Biotype script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the Gradient script.
    Gradient gradient = new Gradient(binding).run();
    algorithmsScript = gradient;
  }
}
