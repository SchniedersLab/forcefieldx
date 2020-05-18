//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy.test

import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.GradientOptions
import ffx.potential.cli.PotentialScript
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
@Command(description = " Test the potential energy gradient.", name = "ffxc test.Gradient")
class Gradient extends PotentialScript {

  @Mixin
  GradientOptions gradientOptions

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

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

    if (!init()) {
      return this
    }

    String modelFilename
    if (filenames != null && filenames.size() > 0) {
      modelFilename = filenames.get(0)
      MolecularAssembly[] assemblies = potentialFunctions.openAll(modelFilename)
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      modelFilename = activeAssembly.getFile().getAbsolutePath()
    }

    logger.info("\n Testing the atomic coordinate gradient of " + modelFilename + "\n")

    energy = activeAssembly.getPotentialEnergy()
    Atom[] atoms = activeAssembly.getAtomArray()
    int nAtoms = atoms.length

    // Apply atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    // Finite-difference step size in Angstroms.
    double step = gradientOptions.dx
    logger.info(" Finite-difference step size:\t" + step)

    // Print out the energy for each step.
    boolean print = gradientOptions.verbose
    logger.info(" Verbose printing:\t\t" + print)

    // Collect atoms to test.
    List<Integer> atomsToTest
    if (gradientOptions.gradientAtoms.equalsIgnoreCase("NONE")) {
      logger.info(" The gradient of no atoms will be evaluated.")
      return this
    } else if (gradientOptions.gradientAtoms.equalsIgnoreCase("ALL")) {
      logger.info(" Checking gradient for all active atoms.\n")
      atomsToTest = new ArrayList<>()
      IntStream.range(0, nAtoms).forEach(val -> atomsToTest.add(val))
    } else {
      atomsToTest = parseAtomRanges(" Gradient atoms", gradientOptions.gradientAtoms, nAtoms)
      logger.info(
          " Checking gradient for active atoms in the range: " + gradientOptions.gradientAtoms +
              "\n")
    }

    // Map selected atom numbers to active atom numbers
    Map<Integer, Integer> allToActive = new HashMap<>()
    int nActive = 0
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i]
      if (atom.isActive()) {
        allToActive.put(i, nActive)
        nActive++
      }
    }

    // Collect analytic gradient.
    double n = energy.getNumberOfVariables()
    double[] x = new double[n]
    double[] g = new double[n]
    energy.getCoordinates(x)
    energy.energyAndGradient(x, g)
    int index = 0
    double[][] allAnalytic = new double[nAtoms][3]
    for (Atom a : atoms) {
      a.getXYZGradient(allAnalytic[index++])
    }

    // Upper bound for a typical atomic gradient.
    double expGrad = 1000.0
    double gradientTolerance = gradientOptions.tolerance
    double width = 2.0 * step
    double avLen = 0.0
    double avGrad = 0.0
    double expGrad2 = expGrad * expGrad

    int nTested = 0
    for (int k : atomsToTest) {
      Atom a0 = atoms[k]
      if (!a0.isActive()) {
        continue
      }

      nTested++
      double[] analytic = allAnalytic[k]

      // Coordinate array only includes active atoms.
      int i = allToActive.get(k)
      int i3 = i * 3
      int i0 = i3 + 0
      int i1 = i3 + 1
      int i2 = i3 + 2
      double[] numeric = new double[3]

      // Find numeric dX
      double orig = x[i0]
      x[i0] = x[i0] + step
      double e = energy.energy(x)
      x[i0] = orig - step
      e -= energy.energy(x)
      x[i0] = orig
      numeric[0] = e / width

      // Find numeric dY
      orig = x[i1]
      x[i1] = x[i1] + step
      e = energy.energy(x)
      x[i1] = orig - step
      e -= energy.energy(x)
      x[i1] = orig
      numeric[1] = e / width

      // Find numeric dZ
      orig = x[i2]
      x[i2] = x[i2] + step
      e = energy.energy(x)
      x[i2] = orig - step
      e -= energy.energy(x)
      x[i2] = orig
      numeric[2] = e / width

      double dx = analytic[0] - numeric[0]
      double dy = analytic[1] - numeric[1]
      double dz = analytic[2] - numeric[2]
      double len = dx * dx + dy * dy + dz * dz
      avLen += len
      len = sqrt(len)

      double grad2 =
          analytic[0] * analytic[0] + analytic[1] * analytic[1] + analytic[2] * analytic[2]
      avGrad += grad2

      if (len > gradientTolerance) {
        logger.info(format(" %s\n Failed: %10.6f\n", a0.toString(), len)
            + format(" Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
            + format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
        ++nFailures
      } else {
        logger.info(format(" %s\n Passed: %10.6f\n", a0.toString(), len)
            + format(" Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
            + format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
      }

      if (grad2 > expGrad2) {
        logger.info(format(" Atom %d has an unusually large gradient: %10.6f", i + 1, Math.sqrt(
            grad2)))
      }
      logger.info("\n")
    }

    avLen = avLen / nTested
    avLen = sqrt(avLen)
    if (avLen > gradientTolerance) {
      logger.info(format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen,
          gradientTolerance))
    } else {
      logger.info(format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen,
          gradientTolerance))
    }
    logger.info(format(" Number of atoms failing analytic test: %d", nFailures))

    avGrad = avGrad / nTested
    avGrad = sqrt(avGrad)
    if (avGrad > expGrad) {
      logger.info(format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad))
    } else {
      logger.info(format(" RMS gradient: %10.6f", avGrad))
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (energy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(energy)
    }
    return potentials
  }

}
