//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.groovy

import org.apache.commons.io.FilenameUtils

import ffx.crystal.Crystal
import ffx.crystal.SymOp
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The SaveAsXYZ script saves a file as an XYZ file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsXYZ [options] &lt;filename&gt;
 */
@Command(description = " Save the system as an XYZ file.", name = "ffxc SaveAsXYZ")
class SaveAsXYZ extends PotentialScript {

    /**
     * -p or --pos-offset to set the positive atom type offset
     */
    @Option(names = ['--pos-offset', '-p'], paramLabel = "0",
            description = 'Positive offset of atom types in the new file')
    int posOffset = 0

    /**
     * -n or --neg-offset to set the negative atom type offset
     */
    @Option(names = ['--neg-offset', '-n'], paramLabel = "0",
            description = 'Negative offset of atom types in the new file.')
    int negOffset = 0

    /**
     * -r or --random to apply a random Cartesian symmetry operator with the specified translation range -X .. X (no default).
     */
    @Option(names = ['--random', '-r'],
            description = 'Apply a random Cartesian SymOp with translation range -X .. X.')
    double scalar = -1

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

    /**
     * Execute the script.
     */
    @Override
    SaveAsXYZ run() {

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

        String modelFilename = activeAssembly.getFile().getAbsolutePath()

        int offset = 0

        // Positive offset atom types.
        if (posOffset > 0) {
            offset = posOffset;
        }

        // Negative offset atom types.
        if (negOffset > 0) {
            offset = negOffset;
            offset = -offset;
        }

        logger.info("\n Writing out XYZ for " + modelFilename)

        // Offset atom type numbers.
        if (offset != 0) {
            logger.info("\n Offset atom types by " + offset)
            ForceField forceField = activeAssembly.getForceField()
            forceField.renumberForceField(0, offset, 0)
        }

        if (scalar > 0.0) {
            SymOp symOp = SymOp.randomSymOpFactory(scalar)
            logger.info(String.format("\n Applying random Cartesian SymOp\n: %s", symOp.toString()))
            Crystal crystal = activeAssembly.getCrystal()
            Atom[] atoms = activeAssembly.getAtomArray()
            double[] xyz = new double[3];
            for (int i = 0; i < atoms.length; i++) {
                atoms[i].getXYZ(xyz)
                crystal.applyCartesianSymOp(xyz, xyz, symOp)
                atoms[i].setXYZ(xyz)
            }
        }

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.getAbsolutePath();
        String fileName = FilenameUtils.getName(modelFilename)
        fileName = FilenameUtils.removeExtension(fileName) + ".xyz"
        File modelFile = new File(dirName + File.separator + fileName)

        potentialFunctions.save(activeAssembly, modelFile)

        return this
    }
}
