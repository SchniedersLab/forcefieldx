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

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.algorithms.optimize.MinimizeOpenMM
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.MolecularAssembly
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The Minimize script uses OpenMM accelerated L-BFGS algorithm to minimize the 
 * energy of a molecular system.
 * <br>
 * Usage:
 * <br>
 * ffxc MinimizeOpenMM [options] &lt;filename&gt;
 */
@Command(description = " Run OpenMM Accelerated L-BFGS minimization on a system.", name = "ffxc MinimizeOpenMM")
class MinimizerOpenMM extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate files in PDB or XYZ format.')

    private List<String> filenames
    private ForceFieldEnergy forceFieldEnergy

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    @Override
    MinimizerOpenMM run() {

        if (!init()) {
            return this
        }

        if (System.getProperty("platform") != null && !System.getProperty("platform").isEmpty()) {
            System.setProperty("platform", System.getProperty("platform"))
        } else {
            System.setProperty("platform", "OMM")
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        switch (forceFieldEnergy.getPlatform()) {
            case ForceFieldEnergy.Platform.OMM:
            case ForceFieldEnergy.Platform.OMM_CUDA:
            case ForceFieldEnergy.Platform.OMM_OPENCL:
            case ForceFieldEnergy.Platform.OMM_OPTCPU:
            case ForceFieldEnergy.Platform.OMM_REF:
                logger.fine(" Platform is appropriate for OpenMM Minimization.")
                break
            case ForceFieldEnergy.Platform.FFX:
            default:
                logger.severe(String.format(" Platform %s is inappropriate for OpenMM minimization. Please explicitly specify an OpenMM platform.",
                        forceFieldEnergy.getPlatform()))
                break
        }

        if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
            MinimizeOpenMM minimizeOpenMM = new MinimizeOpenMM(activeAssembly)
            minimizeOpenMM.minimize(minimizeOptions.eps, minimizeOptions.iterations)

            if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                saveDir = new File(FilenameUtils.getFullPath(modelFilename))
            }

            String dirName = saveDir.toString() + File.separator
            String fileName = FilenameUtils.getName(modelFilename)
            String ext = FilenameUtils.getExtension(fileName)
            fileName = FilenameUtils.removeExtension(fileName)

            File saveFile
            if (ext.toUpperCase().contains("XYZ")) {
                saveFile = new File(dirName + fileName + ".xyz")
                algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
            } else if (ext.toUpperCase().contains("ARC")) {
                saveFile = new File(dirName + fileName + ".arc")
                algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
            } else {
                saveFile = new File(dirName + fileName + ".pdb")
                algorithmFunctions.saveAsPDB(activeAssembly, saveFile)
            }

            SystemFilter systemFilter = algorithmFunctions.getFilter()
            saveFile = activeAssembly.getFile()

            if (systemFilter instanceof XYZFilter) {
                XYZFilter xyzFilter = (XYZFilter) systemFilter
                while (xyzFilter.readNext()) {
                    minimizeOpenMM.minimize(minimizeOptions.eps, minimizeOptions.iterations)
                    boolean append = true
                    xyzFilter.writeFile(saveFile, append)
                }
            }
        } else {
            logger.severe(" Could not start OpenMM minimization.")
        }

        return this
    }

    @Override
    List<Potential> getPotentials() {
        return forceFieldEnergy == null ? Collections.emptyList() : Collections.singletonList(forceFieldEnergy)
    }
}
