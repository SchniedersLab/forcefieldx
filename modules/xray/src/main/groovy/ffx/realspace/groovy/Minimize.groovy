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
package ffx.realspace.groovy

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.potential.MolecularAssembly
import ffx.potential.cli.AtomSelectionOptions
import ffx.realspace.RealSpaceData
import ffx.realspace.cli.RealSpaceOptions
import ffx.realspace.parsers.RealSpaceFile
import ffx.xray.RefinementMinimize

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The Real Space Minimize script.
 * <br>
 * Usage:
 * <br>
 * ffxc realspace.Minimize [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Minimization on a Real Space target.", name = "ffxc realspace.Minimize")
class Minimize extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    @Mixin
    RealSpaceOptions realSpaceOptions

    @Mixin
    AtomSelectionOptions atomSelectionOptions

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames
    private MolecularAssembly[] assemblies
    private RealSpaceData realspaceData

    @Override
    Minimize run() {

        if (!init()) {
            return this
        }

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            activeAssembly = algorithmFunctions.open(filenames.get(0))
            modelFilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }
        assemblies = [activeAssembly] as MolecularAssembly[]

        logger.info("\n Running Real Space Minimization on " + modelFilename)

        // Set atom (in)active selections.
        atomSelectionOptions.setActiveAtoms(activeAssembly)

        RealSpaceFile[] mapFiles = realSpaceOptions.processData(filenames, activeAssembly).toArray(new RealSpaceFile[0])
        realspaceData = new RealSpaceData(assemblies, activeAssembly.getProperties(), activeAssembly.getParallelTeam(), mapFiles)

        // Beginning target function.
        algorithmFunctions.energy(activeAssembly)

        RefinementMinimize refinementMinimize = new RefinementMinimize(realspaceData, realSpaceOptions.refinementMode)

        double eps = minimizeOptions.eps
        int maxiter = minimizeOptions.iterations
        if (eps < 0.0) {
            eps = 1.0
        }

        if (maxiter < Integer.MAX_VALUE) {
            logger.info(String.format("\n RMS gradient convergence criteria: %8.5f, Maximum iterations %d", eps, maxiter))
        } else {
            logger.info(String.format("\n RMS gradient convergence criteria: %8.5f", eps))
        }

        refinementMinimize.minimize(eps, maxiter)

        // Final target function.
        algorithmFunctions.energy(activeAssembly)

        // Suffix to append to output data
        String suffix = "_refine"
        algorithmFunctions.saveAsPDB(assemblies, new File(FilenameUtils.removeExtension(modelFilename) + suffix + ".pdb"))

        return this
    }

    @Override
    boolean destroyPotentials() {
        boolean destroyedSuccess = true
        if (assemblies != null) {
            if (realspaceData != null) {
                destroyedSuccess = realspaceData.destroy()
            }
            for (int i = 1; i < assemblies.length; i++) {
                destroyedSuccess = destroyedSuccess && assemblies[i].getPotentialEnergy().destroy()
            }
        }
        return destroyedSuccess
    }
}

