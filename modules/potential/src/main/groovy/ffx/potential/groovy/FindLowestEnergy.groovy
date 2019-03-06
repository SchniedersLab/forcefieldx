package ffx.potential.groovy

import ffx.potential.ForceFieldEnergy
import ffx.potential.AssemblyState
import ffx.potential.bonded.Atom
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.Superpose
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The Superpose script superposes all molecules in an arc file to the first molecule in the arc file and reports the RMSD.
 * <br>
 * Usage:
 * <br>
 * ffxc Superpose [options] &lt;filename&gt;
 */
@Command(description = " Save the system as a PDB file.", name = "ffxc SaveAsPDB")
class FindLowestEnergy extends PotentialScript {
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

    public ForceFieldEnergy forceFieldEnergy = null

    public double energy = 0.0
    /**
     * Execute the script.
     */
    @Override
    FindLowestEnergy run() {
        if (!init()) {
            return this
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)
        energy = forceFieldEnergy.energy(x, true)
        AssemblyState assemblyState
        SystemFilter systemFilter = potentialFunctions.getFilter()
        if (systemFilter instanceof XYZFilter) {
            XYZFilter xyzFilter = (XYZFilter) systemFilter

            while (xyzFilter.readNext()) {
                forceFieldEnergy.getCoordinates(x)
                double newEnergy = forceFieldEnergy.energy(x, true)
                if (energy > newEnergy) {
                    assemblyState = new AssemblyState(activeAssembly)
                    energy = newEnergy
                }
                else if(energy < newEnergy) {
                    return
                }
            }
        }

        logger.info("The lowest energy is: "+energy)



        File saveDir = baseDir
        String modelFilename = assemblyState.mola.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
        potentialFunctions.saveAsPDB(assemblyState.mola, saveFile)

        return this
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
*/
