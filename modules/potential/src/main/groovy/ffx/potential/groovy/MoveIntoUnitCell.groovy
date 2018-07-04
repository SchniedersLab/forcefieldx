package ffx.potential.groovy

import org.apache.commons.io.FilenameUtils

import ffx.potential.bonded.Atom
import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The MoveIntoUnitCell script moves the center of mass of each molecule into the unit cell.
 * <br>
 * Usage:
 * <br>
 * ffxc MoveIntoUnitCell &lt;filename&gt;
 */
@Command(description = " Move all molecules into the unit cell.", name = "ffxc MoveIntoUnitCell")
class MoveIntoUnitCell extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    public writeFiles = true
    public double[][] origCoordinates = null;
    public double[][] unitCellCoordinates = null;

    /**
     * Execute the script.
     */
    @Override
    MoveIntoUnitCell run() {

        if (!init()) {
            return
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            assemblies = [activeAssembly]
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        logger.info("\n Moving molecular centers of mass into the unit cell for " + modelFilename + "\n");

        // Loop over each system.
        for (int i = 0; i < assemblies.length; i++) {
            MolecularAssembly molecularAssembly = assemblies[i]

            Atom[] atoms = molecularAssembly.getAtomArray()
            int n = atoms.length
            origCoordinates = new double[n][3]
            unitCellCoordinates = new double[n][3]
            int index = 0
            for (Atom atom : atoms) {
                origCoordinates[index][0] = atom.getX()
                origCoordinates[index][1] = atom.getY()
                origCoordinates[index++][2] = atom.getZ()
            }
            molecularAssembly.moveAllIntoUnitCell()
            index = 0
            for (Atom atom : atoms) {
                unitCellCoordinates[index][0] = atom.getX()
                unitCellCoordinates[index][1] = atom.getY()
                unitCellCoordinates[index++][2] = atom.getZ()
            }
        }

        if (writeFiles) {
            String ext = FilenameUtils.getExtension(modelFilename)
            modelFilename = FilenameUtils.removeExtension(modelFilename)
            if (ext.toUpperCase().contains("XYZ")) {
                potentialFunctions.saveAsXYZ(assemblies[0], new File(modelFilename + ".xyz"))
            } else {
                potentialFunctions.saveAsPDB(assemblies, new File(modelFilename + ".pdb"))
            }
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
