//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.OSTOptions
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.numerics.switching.PowerSwitch
import ffx.potential.DualTopologyEnergy
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The OSTBias script tests the Transition-Tempered Orthogonal Space Random Walk Potential.
 * <br>
 * Usage:
 * <br>
 * ffxc test.OSTBias [options] &lt;filename&gt;
 */
@Command(description = " Search for minimum energy polymoprhs for a given space group.", name = "ffxc test.CrystalSearch")
class OSTBias extends AlgorithmsScript {

  /**
   * -s1 or --start1 defines the first softcored atom for the first topology.
   */
  @Option(names = ['--s1', '--start1'], defaultValue = '1',
      description = 'Starting ligand atom for 1st topology')
  int ligandStart
  /**
   * -s2 or --start2 defines the first softcored atom for the second topology.
   */
  @Option(names = ['--s2', '--start2'], defaultValue = '1',
      description = 'Starting ligand atom for 2nd topology')
  int ligandStart2
  /**
   * -f1 or --final1 defines the last softcored atom for the first topology.
   */
  @Option(names = ['--f1', '--final1'], defaultValue = '-1',
      description = 'Final ligand atom for the 1st topology')
  int ligandStop
  /**
   * -f2 or --final2 defines the last softcored atom for the second topology.
   */
  @Option(names = ['--f2', '--final2'], defaultValue = '-1',
      description = 'Final ligand atom for the 2nd topology')
  int ligandStop2
  /**
   * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
   */
  @Option(names = ['--es1', '--noElecStart1'],
      defaultValue = '1', description = 'Starting no-electrostatics atom for 1st topology')
  int noElecStart
  /**
   * -es2 or --noElecStart2 defines the first atom of the second topology to have no electrostatics.
   */
  @Option(names = ['--es2', '--noElecStart2'],
      defaultValue = '1', description = 'Starting no-electrostatics atom for 2nd topology')
  int noElecStart2
  /**
   * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
   */
  @Option(names = ['--ef1', '--noElecFinal1'],
      defaultValue = '-1', description = 'Final no-electrostatics atom for 1st topology')
  int noElecStop
  /**
   * -ef2 or --noElecFinal2 defines the last atom of the second topology to have no electrostatics.
   */
  @Option(names = ['--ef2', '--noElecFinal2'],
      defaultValue = '-1', description = 'Final no-electrostatics atom for 2nd topology')
  int noElecStop2
  /**
   * -l or --lambda sets the lambda value to minimize at.
   */
  @Option(names = ['-l', '--lambda'], defaultValue = '0.5',
      description = 'Initial lambda value.')
  double initialLambda

  /**
   * A PDB or XYZ filename.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "A PDB or XYZ input file.")
  private List<String> filenames


  /**
   * CrystalSearch Constructor.
   */
  OSTBias() {
    this(new Binding())
  }

  /**
   * CrystalSearch Constructor.
   * @param binding The Groovy Binding to use.
   */
  OSTBias(Binding binding) {
    super(binding)
  }

  @Override
  OSTBias run() {

    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    MolecularAssembly molecularAssembly = getActiveAssembly(filenames.get(0))
    if (molecularAssembly == null) {
      logger.info(helpString())
      return this
    }

    String filename = molecularAssembly.getFile().getAbsolutePath()

    println("\n Testing Orthogonal Space Random Walk on " + filename)

    File structureFile = new File(FilenameUtils.normalize(filename))
    structureFile = new File(structureFile.getAbsolutePath())
    String baseFilename = FilenameUtils.removeExtension(structureFile.getName())
    File histogramRestart = new File(baseFilename + ".his")

    if (!histogramRestart.exists()) {
      logger.severe("\n Histogram restart file does not exist.")
    } else if (!histogramRestart.canRead()) {
      logger.severe("\n Histogram restart file can not be read.")
    }

    // For a single process job, try to get the restart files from the current directory.
    File lambdaRestart = new File(baseFilename + ".lam")
    File dyn = new File(baseFilename + ".dyn")
    if (!dyn.exists()) {
      dyn = null
    }

    // Turn on computation of lambda derivatives.
    System.setProperty("lambdaterm", "true")

    // Relative free energies via the DualTopologyEnergy class require different
    // default OST parameters than absolute free energies.
    if (filenames.size() == 2) {
      // Ligand vapor electrostatics are not calculated. This cancels when the
      // difference between protein and water environments is considered.
      System.setProperty("ligand-vapor-elec", "false")
    }

    // Get a reference to the first system's ForceFieldEnergy and atom array.
    ForceFieldEnergy energy = molecularAssembly.getPotentialEnergy()
    Atom[] atoms = molecularAssembly.getAtomArray()

    // Apply the no electrostatics atom selection
    if (ligandStart < 1) {
      ligandStart = 1
    }
    if (ligandStop > atoms.length) {
      ligandStop = atoms.length
    }
    // Apply the ligand atom selection
    for (int i = ligandStart; i <= ligandStop; i++) {
      Atom ai = atoms[i - 1]
      ai.setApplyLambda(true)
      ai.print()
    }

    // Apply the no electrostatics atom selection
    if (noElecStart < 1) {
      noElecStart = 1
    }
    if (noElecStop > atoms.length) {
      noElecStop = atoms.length
    }
    for (int i = noElecStart; i <= noElecStop; i++) {
      Atom ai = atoms[i - 1]
      ai.setElectrostatics(false)
      ai.print()
    }

    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    energy.getCrystal().setSpecialPositionCutoff((double) 0.0)
    // OST will be configured for either single or dual topology.
    OrthogonalSpaceTempering orthogonalSpaceTempering
    // Save a reference to the first topology.
    MolecularAssembly topology1 = molecularAssembly

    if (filenames.size() == 1) {
      orthogonalSpaceTempering =
          OSTOptions.constructOST(energy, lambdaRestart, histogramRestart, molecularAssembly,
              molecularAssembly.getProperties(), null)
    } else {
      // Open the 2nd topology.
      filename = filenames.get(1)
      MolecularAssembly active = algorithmFunctions.open(filename)
      energy = active.getPotentialEnergy()
      atoms = active.getAtomArray()

      // Apply the ligand atom selection for the 2nd topology.
      // Apply the no electrostatics atom selection
      if (ligandStart2 < 1) {
        ligandStart2 = 1
      }
      if (ligandStop2 > atoms.length) {
        ligandStop2 = atoms.length
      }
      // Apply the ligand atom selection
      for (int i = ligandStart2; i <= ligandStop2; i++) {
        Atom ai = atoms[i - 1]
        ai.setApplyLambda(true)
        ai.print()
      }

      // Apply the no electrostatics atom selection
      if (noElecStart2 < 1) {
        noElecStart2 = 1
      }
      if (noElecStop2 > atoms.length) {
        noElecStop2 = atoms.length
      }
      for (int i = noElecStart2; i <= noElecStop2; i++) {
        Atom ai = atoms[i - 1]
        ai.setElectrostatics(false)
        ai.print()
      }

      // Save a reference to the second topology.
      // Turn off checks for overlapping atoms, which is expected for lambda=0.
      energy.getCrystal().setSpecialPositionCutoff((double) 0.0)
      // Create the DualTopology potential energy.
      DualTopologyEnergy dualTopologyEnergy = new DualTopologyEnergy(topology1, active,
          new PowerSwitch())
      orthogonalSpaceTempering =
          OSTOptions.constructOST(dualTopologyEnergy, lambdaRestart, histogramRestart, active,
              active.getProperties(), null)
    }

    /**
     * Stop propagating lambda to prevent adding new Gaussians
     * to the biasing potential, which would introduce artifacts into the
     * finite-difference derivatives.
     */
    orthogonalSpaceTempering.setPropagateLambda(false)
    orthogonalSpaceTempering.setLambda(initialLambda)
    int n = orthogonalSpaceTempering.getNumberOfVariables()

    assert (n % 3 == 0)
    n = (int) (n / 3)

    // Finite-difference step size.
    double step = 1.0e-5

    double[] x = new double[3 * n]
    double[] analytic = new double[3 * n]
    double[] g = new double[3 * n]
    double[] numeric = new double[3]
    orthogonalSpaceTempering.getCoordinates(x)

    // Test Lambda gradients.
    for (int j = 0; j < 3; j++) {
      double lambda = initialLambda - 0.001 + 0.001 * j

      if (lambda - step < 0.0) {
        continue
      } else if (lambda + step > 1.0) {
        continue
      } else {
        orthogonalSpaceTempering.setLambda(lambda)
      }

      // Calculate the energy and analytic dE/dX
      double eL = orthogonalSpaceTempering.energyAndGradient(x, g)

      // Analytic dEdL
      double dEdLambda = orthogonalSpaceTempering.getTotaldEdLambda()

      // Calculate the finite-difference dEdL
      orthogonalSpaceTempering.setLambda(lambda + step)
      double lp = orthogonalSpaceTempering.energyAndGradient(x, g)

      orthogonalSpaceTempering.setLambda(lambda - step)
      double lm = orthogonalSpaceTempering.energyAndGradient(x, g)

      double dedl = (lp - lm) / (2.0 * step)

      logger.info(format(" Analytic dE/dL:   %15.8f", dEdLambda))
      logger.info(format(" Numeric  dE/dL:   %15.8f\n", dedl))

      // Calculate analytic dE/dX/dL
      orthogonalSpaceTempering.setLambda(lambda)
      double e
      double orig
      double gradientTolerance = 1.0e-3

      // Calculate finite-difference coordinate gradient
      for (int i = 0; i < n; i++) {
        int i3 = i * 3
        int i0 = i3 + 0
        int i1 = i3 + 1
        int i2 = i3 + 2

        // Calculate the analytic dE/dX
        orthogonalSpaceTempering.energyAndGradient(x, analytic)

        // Find numeric dX
        orig = x[i0]
        x[i0] = orig + step
        e = orthogonalSpaceTempering.energyAndGradient(x, g)
        x[i0] = orig - step
        e = e - orthogonalSpaceTempering.energyAndGradient(x, g)
        x[i0] = orig
        numeric[0] = e / (2.0 * step)

        // Find numeric dY
        orig = x[i1]
        x[i1] = orig + step
        e = orthogonalSpaceTempering.energyAndGradient(x, g)
        x[i1] = orig - step
        e = e - orthogonalSpaceTempering.energyAndGradient(x, g)
        x[i1] = orig
        numeric[1] = e / (2.0 * step)

        // Find numeric dZ
        orig = x[i2]
        x[i2] = orig + step
        e = orthogonalSpaceTempering.energyAndGradient(x, g)
        x[i2] = orig - step
        e = e - orthogonalSpaceTempering.energyAndGradient(x, g)
        x[i2] = orig
        numeric[2] = e / (2.0 * step)

        double dx = analytic[i0] - numeric[0]
        double dy = analytic[i1] - numeric[1]
        double dz = analytic[i2] - numeric[2]
        double len = Math.sqrt(dx * dx + dy * dy + dz * dz)
        if (len > gradientTolerance) {
          logger.info(" Atom " + i + format(" failed: %10.6f.", len)
              + format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1],
              analytic[i2])
              + format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
        } else {
          logger.info(" Atom " + i + format(" passed: %10.6f.", len)
              + format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1],
              analytic[i2])
              + format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
        }
      }
    }

    return this
  }
}
