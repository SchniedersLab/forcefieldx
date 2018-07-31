package ffx.potential.groovy.test

import ffx.numerics.Potential
import org.apache.commons.io.FilenameUtils

import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The SaveAsSIFTPDB script saves a file as a PDB file or as a SIFTPDB file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsSIFTPDB [options] &lt;filename&gt;
 */
@Command(description = " Save SIFT scores into the PDB b-factor column.", name = "ffxc SaveAsSIFTPDB")
class SaveAsSIFTPDB extends PotentialScript {

    /**
     * -f or --filename to input file of sift scores to enter
     */
    @Option(names = ['--fileName', '-f'],
            description = 'File of sift scores to enter.')
    boolean siftFilename = null

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    /**
     * Execute the script.
     */
    @Override
    SaveAsSIFTPDB run() {

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

        String filename = activeAssembly.getFile().getAbsolutePath()

        String[] data
        File siftFile = new File(siftFilename)

        if (siftFile.exists()) {
            data = siftFile.text.split('\n')
        }

        logger.info("\n Writing out PDB for " + filename)

        filename = FilenameUtils.removeExtension(filename) + ".pdb";

        if (siftFilename == null) {
            potentialFunctions.saveAsPDB(activeAssembly, new File(filename));
        } else {
            potentialFunctions.saveAsSIFTPDB(activeAssembly, new File(filename), data)
        }

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
