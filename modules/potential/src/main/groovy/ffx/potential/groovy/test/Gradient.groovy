//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy.test

import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.GradientOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.GradientUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import java.util.stream.IntStream

import static ffx.utilities.StringUtils.parseAtomRanges
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.sqrt

/**
 * The Gradient script evaluates the consistency of the energy and gradient.
 * <br>
 * Usage:
 * <br>
 * ffxc test.Gradient [options] &lt;filename&gt;
 */
@Command(description = " Test the potential energy gradient.", name = "test.Gradient")
class Gradient extends PotentialScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  @Mixin
  GradientOptions gradientOptions

  /**
   * The final argument is a single filename in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file", description = "A PDB or XYZ coordinate file.")
  String filename = null

  private ForceFieldEnergy energy
  public int nFailures = 0

  /**
   * Gradient constructor.
   */
  Gradient() {
    this(new Binding())
  }

  /**
   * Gradient constructor.
   * @param binding The Groovy Binding to use.
   */
  Gradient(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Gradient run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Testing the atomic coordinate gradient of " + filename + "\n")

    // Apply atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    energy = activeAssembly.getPotentialEnergy()
    GradientUtils gradientUtils = new GradientUtils(energy, activeAssembly.getAtomArray())

    nFailures = gradientUtils.testGradient(gradientOptions)
    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (energy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList((Potential) energy)
    }
    return potentials
  }

}
