package ffx.algorithms.groovy

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The ApplyTorsions script is used to apply a custom set of torsions to a side 
 * chain.
 * <br>
 * Usage:
 * <br>
 * ffxc ApplyTorsions [options] &lt;filename&gt;
 * @author Jacob Litman
 * @author Michael Schnieders
 */
@Command(description = " Apply a set of torsions to a system.", name = "ffxc ApplyTorsions")
class ApplyTorsions extends AlgorithmsScript {

    /**
     * -c or --chain selects the chain name to use.
     */
    @Option(names = ["-c", "--chain"], paramLabel = ' ',
            description = 'Single character chain name (default is \" \").')
    String chain = " "
    /**
     * -r or --resid selects the residue to apply torsions to.
     */
    @Option(names = ["-r", "--resid"], paramLabel = '1',
            description = 'Residue number.')
    int resID = 1
    /**
     * -t or --torsionSets is a variable-length, colon-delimited set of comma-separated torsion sets; e.g. colons separate rotamers, commas individual values. Should not be left as final argument, as then it attempts to include the filename.
     */
    @Option(names = ["-t", "--torsionSets"], paramLabel = '0,0:180,0', arity = "1..*", split = ':',
            description = 'Torsion sets to apply (torsions comma-separated, sets colon-separated).')
    String[] torSets = null
    /**
     * -n or --nChi is the number of torsions available to this side chain; unspecified torsions are filled with 0.
     */
    @Option(names = ["-n", "--nChi"], paramLabel = '1',
            description = 'Number of torsions (unspecified torsions are filled with 0).')
    int nChi = 1
    /**
     * -v or --videoFile is the name of the file to print torsion snapshots to; defaults to filename_rots.pdb
     */
    @Option(names = ["-v", "--videoFile"],
            description = 'File to print torsion snapshots to.')
    String vidFileName = null

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private List<String> filenames

    @Override
    ApplyTorsions run() {

        if (!init()) {
            return this
        }

        String modelfilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
        }

        String newName = FilenameUtils.getBaseName(modelfilename)

        String videoFile
        if (vidFileName != null) {
            videoFile = vidFileName
        } else {
            videoFile = newName + "_rots.pdb"
        }

        File outFile = new File(newName + ".rotout.tmp")
        outFile.deleteOnExit()

        logger.info("\n Saving torsions for residue number " + resID + " of chain " + chain + ".")

        RotamerLibrary.initializeDefaultAtomicCoordinates(activeAssembly.getChains())
        Polymer polymer = activeAssembly.getChain(chain)
        if (polymer == null) {
            logger.info(" Polymer + " + chain + " does not exist.")
            return
        }
        Residue residue = polymer.getResidue(resID)
        if (residue == null) {
            logger.info(" Residue + " + resID + " does not exist.")
            return
        }

        ffx.algorithms.GenerateRotamers genr = new ffx.algorithms.GenerateRotamers(activeAssembly,
                activeAssembly.getPotentialEnergy(), residue, outFile, nChi, algorithmListener)
        genr.setVideo(videoFile)
        genr.applyAndSaveTorsions(torSets)

        return this
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
