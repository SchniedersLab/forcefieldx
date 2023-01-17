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
package ffx.realspace.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.WriteoutOptions
import ffx.realspace.cli.RealSpaceOptions
import ffx.xray.RefinementEnergy
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static ffx.utilities.TinkerUtils.version

/**
 * The Real Space Dynamics script.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.Dynamics [options] &lt;filename&gt;
 */
@Command(description = " Molecular dynamics on a Real Space target.", name = "realspace.Dynamics")
class Dynamics extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  @Mixin
  DynamicsOptions dynamicsOptions

  @Mixin
  RealSpaceOptions realSpaceOptions

  @Mixin
  WriteoutOptions writeoutOptions

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames
  private RefinementEnergy refinementEnergy

  /**
   * Dynamics constructor.
   */
  Dynamics() {
    this(new Binding())
  }

  /**
   * Dynamics constructor.
   * @param binding The Groovy Binding to use.
   */
  Dynamics(Binding binding) {
    super(binding)
  }

  @Override
  Dynamics run() {

    if (!init()) {
      return this
    }

    dynamicsOptions.init()

    String filename
    if (filenames != null && filenames.size() > 0) {
      activeAssembly = algorithmFunctions.open(filenames.get(0))
      filename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      filename = activeAssembly.getFile().getAbsolutePath()
    }
    MolecularAssembly[] molecularAssemblies = [activeAssembly] as MolecularAssembly[]

    logger.info("\n Running Real Space Dynamics on " + filename)

    atomSelectionOptions.setActiveAtoms(activeAssembly)

    // Create the refinement target.
    refinementEnergy = realSpaceOptions.toRealSpaceEnergy(filenames, molecularAssemblies)

    // Beginning force field energy.
    algorithmFunctions.energy(activeAssembly)

    // Restart File
    File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn")
    if (!dyn.exists()) {
      dyn = null
    }

    MolecularDynamics molDyn = dynamicsOptions.
        getDynamics(writeoutOptions, refinementEnergy, activeAssembly, algorithmListener)
    refinementEnergy.setThermostat(molDyn.getThermostat())

    // Reset velocities (ignored if a restart file is given)
    boolean initVelocities = true
    molDyn.dynamic(dynamicsOptions.steps, dynamicsOptions.dt, dynamicsOptions.report,
            dynamicsOptions.write, dynamicsOptions.temperature, initVelocities, dyn)

    // Final target function.
    algorithmFunctions.energy(activeAssembly)

    File file = version(new File(FilenameUtils.removeExtension(filename) + ".pdb"))
    algorithmFunctions.saveAsPDB(molecularAssemblies, file)

    return this
  }

  @Override
  List<Potential> getPotentials() {
    return refinementEnergy == null ? Collections.emptyList() :
            Collections.singletonList((Potential) refinementEnergy)
  }
}
