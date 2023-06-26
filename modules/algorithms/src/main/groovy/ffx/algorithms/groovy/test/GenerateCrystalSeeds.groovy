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
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.potential.MolecularAssembly
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import ffx.numerics.Potential


/**
 * GenerateCrystalSeeds is a Groovy script that generates a set of molecular orientations in vacuum and
 * calculates the energy of each conformation.
 * <br>
 * Usage:
 * <br>
 * ffxc test.GenerateCrystalSeeds [options] &lt;filename&gt;
 */
@Command(description = " Calculates interaction energies of different molecular orientations and saves low energy orientations.", name = "test.GenerateCrystalSeeds")
class GenerateCrystalSeeds extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    /**
     * One or two filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private String[] filenames

    /**
    * Minimize Constructor.
    */
  GenerateCrystalSeeds() {
    this(new Binding())
  }

  /**
    * Minimize Constructor.
    * @param binding The Groovy Binding to use.
    */
    GenerateCrystalSeeds(Binding binding) {
    super(binding)
  }

  /**
    * {@inheritDoc}
    */
  @Override
  GenerateCrystalSeeds run() {
      // Init the context and bind variables.
      if (!init()) {
          return this
      }

      // Load the MolecularAssemblies of the one or two input files.
      MolecularAssembly[] activeAssemblies = new MolecularAssembly[filenames.length]
      for(int i = 0; i < filenames.length; i++) {
          activeAssembly[i] = getActiveAssembly(filenames[i])
          if (activeAssembly[i] == null) {
              logger.info(helpString())
              return this
          }
      }

      // Set the filenames.
      filenames = new String[activeAssemblies.length]
      for(int i = 0; i < activeAssemblies.length; i++) {
          filenames[i] = activeAssemblies[i].getFile().getAbsolutePath()
      }

      // Run energy on each of the activeAssemblies.



      return this
  }

}

