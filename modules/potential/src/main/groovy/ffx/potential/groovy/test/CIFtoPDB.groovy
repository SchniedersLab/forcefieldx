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
package ffx.potential.groovy.test

import java.nio.file.Path
import java.nio.file.Paths
import static java.lang.String.format

import org.apache.commons.io.FilenameUtils
import org.rcsb.cif.CifIO
import org.rcsb.cif.model.Block
import org.rcsb.cif.model.CifFile
import org.rcsb.cif.model.Column

import ffx.crystal.Crystal
import ffx.crystal.SpaceGroup
import ffx.numerics.Potential
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import static ffx.potential.parsers.PDBFilter.toPDBAtomLine

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The CIFtoPDB script converts a CIF file to PDB format.
 * <br>
 * Usage:
 * <br>
 * ffxc test.CIFtoPDB &lt;filename.cif&gt;
 */
@Command(description = " Convert a CIF file to PDB format.", name = "ffxc test.CIFtoPDB")
class CIFtoPDB extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
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
    CIFtoPDB run() {

        if (!init()) {
            return null
        }

        CifFile cifFile
        Path path
        if (filenames != null && filenames.size() > 0) {
            path = Paths.get(filenames.get(0))
            cifFile = CifIO.readFromPath(path)
        } else {
            logger.info(helpString())
            return null
        }

        String modelFilename = path.toAbsolutePath().toString()
        logger.info("\n Opening CIF file " + path)
        Block firstBlock = cifFile.firstBlock

        // Space Group Number
        // _symmetry_Int_Tables_number      14
        // Unit Cell Parameters
        // _cell_length_a                   7.529(1)
        // _cell_length_b                   11.148(2)
        // _cell_length_c                   15.470(2)
        // _cell_angle_alpha                90
        // _cell_angle_beta                 116.17(1)
        // _cell_angle_gamma                90

        int sgNum = Integer.parseInt(firstBlock.getColumn("symmetry_Int_Tables_number").getStringData(0))
        SpaceGroup sg = SpaceGroup.spaceGroupFactory(sgNum)
        double a = toDouble(firstBlock.getColumn("cell_length_a").getStringData(0))
        double b = toDouble(firstBlock.getColumn("cell_length_b").getStringData(0))
        double c = toDouble(firstBlock.getColumn("cell_length_c").getStringData(0))
        double alpha = toDouble(firstBlock.getColumn("cell_angle_alpha").getStringData(0))
        double beta = toDouble(firstBlock.getColumn("cell_angle_beta").getStringData(0))
        double gamma = toDouble(firstBlock.getColumn("cell_angle_gamma").getStringData(0))
        Crystal crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName)
        logger.info(crystal.toString())

        int nAtoms = firstBlock.getColumn("atom_site_label").getRowCount()
        logger.info(format("\n Number of atoms: %d", nAtoms))

        // CIF file atom site columns.
        //
        // loop_
        // _atom_site_label
        // _atom_site_type_symbol
        // _atom_site_fract_x
        // _atom_site_fract_y
        // _atom_site_fract_z
        // C1 C 0.17860 0.88590 0.16870
        // ...
        Atom[] atoms = new Atom[nAtoms]
        Column labels = firstBlock.getColumn("atom_site_label")
        // Column types = firstBlock.getColumn("atom_site_type_symbol")
        Column xFract = firstBlock.getColumn("atom_site_fract_x")
        Column yFract = firstBlock.getColumn("atom_site_fract_y")
        Column zFract = firstBlock.getColumn("atom_site_fract_z")

        // Define per atom information for the PDB file.
        String resName = "CIF"
        int resID = 1
        double occupancy = 1.0
        double bfactor = 1.0
        char altLoc = ' '
        char chain = 'A'
        String segID = "A"

        // Loop over atoms.
        for (int i = 0; i < nAtoms; i++) {
            String label = labels.getStringData(i)
            // String type = types.getStringData(i)
            double x = toDouble(xFract.getStringData(i))
            double y = toDouble(yFract.getStringData(i))
            double z = toDouble(zFract.getStringData(i))
            double[] xyz = [x, y, z]
            crystal.toCartesianCoordinates(xyz, xyz)
            atoms[i] = new Atom(i + 1, label, altLoc, xyz, resName, resID, chain, occupancy, bfactor, segID)
            atoms[i].setHetero(true)
        }

        // Determine new PDB file name and location.
        File saveDir = baseDir
        File cif = new File(modelFilename)
        String cifName = cif.getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(cifName))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(cifName)
        fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        File modelFile = new File(dirName + fileName)
        File saveFile = potentialFunctions.versionFile(modelFile)

        // Write out PDB file.
        FileWriter fw = new FileWriter(saveFile, false)
        BufferedWriter bw = new BufferedWriter(fw)
        bw.write(crystal.toCRYST1())
        for (int i = 0; i < nAtoms; i++) {
            bw.write(toPDBAtomLine(atoms[i]))
        }
        bw.write("END\n")
        bw.close()

        logger.info("\n Saved PDB file to " + saveFile.getAbsolutePath())

        return this
    }

    /**
     * Remove uncertainty information from "XXX.XXX(X)", which is the part in parentheses "(X)". Then
     * convert the String to a double.
     *
     * @param string The input String, possibly with uncertainty.
     * @return A double value.
     */
    private static double toDouble(String string) {
        string = string.replace('(', ' ').replace(')', ' ').trim();
        return Double.parseDouble(string.split(" +")[0])
    }

    @Override
    List<Potential> getPotentials() {
        return new ArrayList<Potential>()
    }
}
