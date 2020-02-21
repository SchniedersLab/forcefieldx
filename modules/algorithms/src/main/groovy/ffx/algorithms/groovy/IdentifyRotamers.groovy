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
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The IdentifyRotamers script attempts to identify which rotamer each Residue in the system is in.
 * <br>
 * Usage:
 * <br>
 * ffxc IdentifyRotamers [options] &lt;filename&gt;
 */
@Command(description = " Identify the rotamers a system is in.", name = "ffxc IdentifyRotamers")
class IdentifyRotamers extends AlgorithmsScript {

    @Mixin
    ManyBodyOptions mbOpts

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    @Override
    IdentifyRotamers run() {

        if (!init()) {
            return null
        }

        mbOpts.setOriginalCoordinates(false)

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = [algorithmFunctions.open(filenames.get(0))]
            activeAssembly = assemblies[0]
            if (Boolean.parseBoolean(System.getProperty("standardizeAtomNames", "false"))) {
                renameAtomsToPDBStandard(activeAssembly)
            }
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return null
        }

        activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)

        RotamerOptimization rotopt = new RotamerOptimization(
                activeAssembly, activeAssembly.getPotentialEnergy(), algorithmListener)
        mbOpts.initRotamerOptimization(rotopt, activeAssembly)

        List<Residue> residues = rotopt.getResidues()
        RotamerLibrary rLib = mbOpts.getRotamerLibrary()

        for (Residue residue : residues) {
            RotamerLibrary.RotamerGuess guess = rLib.guessRotamer(residue)
            logger.info(guess.toString())
        }

        return this
    }
}
