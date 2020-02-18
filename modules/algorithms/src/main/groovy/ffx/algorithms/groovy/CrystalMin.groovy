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
import static org.apache.commons.math3.util.FastMath.abs

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.algorithms.optimize.CrystalMinimize
import ffx.algorithms.optimize.Minimize
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.MolecularAssembly.FractionalMode
import ffx.potential.XtalEnergy
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The CrystalMin script uses a limited-memory BFGS algorithm to minimize the
 * energy of a crystal, including both coordinates and unit cell parameters.
 * <br>
 * Usage:
 * <br>
 * ffxc CrystalMin [options] &lt;filename&gt;
 */
@Command(description = " Minimize crystal unit cell parameters.", name = "ffxc CrystalMin")
class CrystalMin extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    /**
     * -c or --coords to cycle between lattice and coordinate optimization until both satisfy the convergence criteria.
     */
    @Option(names = ["-c", "--coords"], paramLabel = 'false',
            description = 'Cycle between lattice and coordinate optimization until both satisfy the convergence criteria.')
    boolean coords = false
    /**
     * -f or --fractional to set the optimization to maintain fractional coordinates [ATOM/MOLECULE/OFF].
     */
    @Option(names = ["-f", "--fractional"], paramLabel = 'molecule',
            description = 'Maintain fractional coordinates during lattice optimization [OFF/MOLECULE/ATOM].')
    String fractional = "MOLECULE"
    /**
     * -t or --tensor to print out partial derivatives of the energy with respect to unit cell parameters.
     */
    @Option(names = ["-t", "--tensor"], paramLabel = 'false',
            description = 'Compute partial derivatives of the energy with respect to unit cell parameters.')
    boolean tensor = false

    /**
     * The final argument(s) should be a filename.
     */
    @Parameters(arity = "1", paramLabel = "files", description = 'Atomic coordinate files in PDB or XYZ format.')
    List<String> filenames = null

    private XtalEnergy xtalEnergy
    private CrystalMinimize crystalMinimize

    @Override
    CrystalMin run() {

        if (!init()) {
            return null
        }

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = [algorithmFunctions.open(filenames.get(0))]
            activeAssembly = assemblies[0]
            modelFilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return null
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }

        logger.info("\n Running CrystalMinimize on " + modelFilename)
        logger.info("\n RMS gradient convergence criteria: " + minimizeOptions.eps)

        ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
        xtalEnergy = new XtalEnergy(forceFieldEnergy, activeAssembly)
        xtalEnergy.setFractionalCoordinateMode(FractionalMode.MOLECULE)

        SystemFilter systemFilter = algorithmFunctions.getFilter()

        // Apply fractional coordinate mode.
        if (fractional) {
            try {
                FractionalMode mode = FractionalMode.valueOf(fractional.toUpperCase())
                xtalEnergy.setFractionalCoordinateMode(mode)
            } catch (Exception e) {
                logger.info(" Unrecognized fractional coordinate mode: " + fractional)
                logger.info(" Fractional coordinate mode is set to MOLECULE.")
            }
        }

        runMinimize()

        if (tensor) {
            crystalMinimize.printTensor()
        }

        // Handle Single Topology Cases.
        modelFilename = activeAssembly.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        String ext = FilenameUtils.getExtension(fileName)
        fileName = FilenameUtils.removeExtension(fileName)
        File saveFile
        SystemFilter writeFilter
        if (ext.toUpperCase().contains("XYZ")) {
            saveFile = new File(dirName + fileName + ".xyz")
            writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
            algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
        } else if (ext.toUpperCase().contains("ARC")) {
            saveFile = new File(dirName + fileName + ".arc")
            saveFile = algorithmFunctions.versionFile(saveFile)
            writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
            algorithmFunctions.saveAsXYZ(activeAssembly, saveFile)
        } else {
            saveFile = new File(dirName + fileName + ".pdb")
            saveFile = algorithmFunctions.versionFile(saveFile)
            writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
            int numModels = systemFilter.countNumModels()
            if (numModels > 1) {
                writeFilter.setModelNumbering(0)
            }
            writeFilter.writeFile(saveFile, true, false, false)
        }

        if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
            while (systemFilter.readNext()) {
                if (systemFilter instanceof PDBFilter) {
                    saveFile.append("ENDMDL\n")
                    runMinimize()
                    PDBFilter pdbFilter = (PDBFilter) writeFilter
                    pdbFilter.writeFile(saveFile, true, false, false)
                } else if (systemFilter instanceof XYZFilter) {
                    runMinimize()
                    writeFilter.writeFile(saveFile, true)
                }
            }
            if (systemFilter instanceof PDBFilter) {
                saveFile.append("END\n")
            }
        }

        return this
    }

    void runMinimize() {
        crystalMinimize = new CrystalMinimize(activeAssembly, xtalEnergy, algorithmListener)
        crystalMinimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)
        double energy = crystalMinimize.getEnergy()

        double tolerance = 1.0e-10

        // Complete rounds of coordinate and lattice optimization.
        if (coords) {
            ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
            Minimize minimize = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
            while (true) {
                // Complete a round of coordinate optimization.
                minimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)
                double newEnergy = minimize.getEnergy()
                int status = minimize.getStatus()
                if (status != 0) {
                    break
                }
                if (abs(newEnergy - energy) <= tolerance) {
                    break
                }
                energy = newEnergy

                // Complete a round of lattice optimization.
                crystalMinimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)
                newEnergy = crystalMinimize.getEnergy()
                status = crystalMinimize.getStatus()
                if (status != 0) {
                    break
                }
                if (abs(newEnergy - energy) <= tolerance) {
                    break
                }
                energy = newEnergy
            }
        }
    }

    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (xtalEnergy == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = Collections.singletonList(xtalEnergy)
        }
        return potentials
    }

}
