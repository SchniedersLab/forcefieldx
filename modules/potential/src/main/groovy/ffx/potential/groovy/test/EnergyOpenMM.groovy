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
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level

import static java.lang.String.format

/**
 * The Energy script compares the energy of a system with OpenMM versus FFX platforms.
 * <br>
 * Usage:
 * <br>
 * ffxc test.EnergyOpenMM &lt;filename&gt;
 */
@Command(description = " Compute energies and gradients with both FFX and OpenMM.", name = "ffxc test.EnergyOpenMM")
class EnergyOpenMM extends PotentialScript {

  public double energy = 0.0
  public ForceFieldEnergy forceFieldEnergy = null

  /**
   * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
   */
  @Option(names = ['--es1', '--noElecStart1'], paramLabel = "1",
      description = 'Starting no-electrostatics atom for 1st topology')
  private int es1 = 1

  /**
   * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
   */
  @Option(names = ['--ef1', '--noElecFinal1'], paramLabel = "-1",
      description = 'Final no-electrostatics atom for 1st topology')
  private int ef1 = -1

  /**
   * The final argument is a single filename in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "A PDB or XYZ coordinate file.")
  String filename = null

  /**
   * EnergyOpenMM Constructor.
   */
  EnergyOpenMM() {
    this(new Binding())
  }

  /**
   * EnergyOpenMM Constructor.
   * @param binding Groovy Binding to use.
   */
  EnergyOpenMM(Binding binding) {
    super(binding)
  }

  @Override
  EnergyOpenMM run() {

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

    logger.info(" Running Energy on " + filename)

    forceFieldEnergy = activeAssembly.getPotentialEnergy()
    if (!forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
      logger.severe(" test.EnergyOpenMM requires use of an OpenMM platform!")
    }
    Atom[] atoms = activeAssembly.getAtomArray()

    // Apply the no electrostatics atom selection
    int noElecStart = es1
    noElecStart = (noElecStart < 1) ? 1 : noElecStart

    int noElecStop = ef1
    noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop

    if (noElecStart <= noElecStop) {
      logger.info(format(" Disabling electrostatics for atoms %d (%s) to %d (%s).",
          noElecStart, atoms[noElecStart - 1], noElecStop, atoms[noElecStop - 1]))
    }
    for (int i = noElecStart; i <= noElecStop; i++) {
      Atom ai = atoms[i - 1]
      ai.setElectrostatics(false)
      ai.print(Level.FINE)
    }

    int nVars = forceFieldEnergy.getNumberOfVariables()
    double[] x = new double[nVars]
    forceFieldEnergy.getCoordinates(x)

    double[] gFFX = new double[nVars]
    double[] gOMM = new double[nVars]

    ForceFieldEnergyOpenMM feOMM = (ForceFieldEnergyOpenMM) forceFieldEnergy

    double dE = feOMM.energyAndGradientFFX(x, gFFX, true)
    dE = dE - feOMM.energyAndGradient(x, gOMM, true)
    logger.info(format(" Difference in energy: %14.8g kcal/mol", dE))
    int nActAts = (int) (nVars / 3)

    for (int i = 0; i < nActAts; i++) {
      double[] dG = new double[3]
      int i3 = 3 * i
      double rmse = 0
      for (int j = 0; j < 3; j++) {
        dG[j] = gFFX[i3 + j] - gOMM[i3 + j]
        rmse += (dG[j] * dG[j])
      }
      rmse /= 3
      rmse = Math.sqrt(rmse)
      Level toPrint = (rmse > 1E-8) ? Level.INFO : Level.FINE
      logger.log(toPrint, format(" Active atom %d dG RMSE: %14.8g kcal/mol/A", i, rmse))
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (forceFieldEnergy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(forceFieldEnergy)
    }
    return potentials
  }

}
