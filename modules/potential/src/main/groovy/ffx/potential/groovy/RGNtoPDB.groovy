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
package ffx.potential.groovy

import static java.lang.String.format

import org.apache.commons.io.FilenameUtils
import org.biojava.nbio.core.sequence.ProteinSequence
import org.biojava.nbio.core.sequence.io.FastaReaderHelper

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.ResidueEnumerations
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import static ffx.potential.parsers.PDBFilter.toPDBAtomLine

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The RGNtoPDB converts RGN output (a *.tertiary file) to PDB format.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc RGNtoPDB [options] &lt;filename.tertiary&gt; &lt;filename.fasta&gt;
 */
@Command(description = " RGNtoPDB converts RGN *tertiary file to PDB with side-chains format.", name = "ffxc RGNtoPDB")
class RGNtoPDB extends PotentialScript {

    /**
     * The final argument(s) should be an RGN output file and a sequence file in FASTA format.
     */
    @Parameters(arity = "2", paramLabel = "files",
            description = 'RGN output and a FASTA file.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    RGNtoPDB run() {

        if (!init()) {
            return null
        }

        if (filenames == null || filenames.size() < 2) {
            logger.info(helpString())
            return null
        }

        String rgnName = filenames.get(0)
        String fastaName = filenames.get(1)

        logger.info("\n Opening FASTA " + fastaName)

        LinkedHashMap<String, ProteinSequence> fastaData = FastaReaderHelper.readFastaProteinSequence(new File(fastaName))
        String sequence = fastaData.values().toArray()[0]

        logger.info("\n Opening RGN " + rgnName)

        // Read lines from RGN file.
        File rgnFile = new File(rgnName)

        String[] lines = new String[5]
        String[][] tokenizedLines = new String[5][]

        BufferedReader cr = new BufferedReader(new FileReader(rgnFile))

        int lineNumber = 0
        while (lineNumber < 5) {
            lines[lineNumber] = cr.readLine().trim()
            tokenizedLines[lineNumber] = lines[lineNumber].split(" +")
            lineNumber++
        }
        cr.close()

        int nAtoms = tokenizedLines[2].length
        int nAmino = (int) (nAtoms / 3)
        int remainder = (int) (nAtoms % 3)
        logger.info(format(" RGN record has %d atoms and %d amino acids", nAtoms, nAmino))

        if (remainder != 0) {
            logger.info(format(" RGN record %d atoms are not divisible by 3 (%d atoms left).", remainder))
            return null
        }

        File saveDir = baseDir
        rgnName = rgnFile.getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(rgnName))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(rgnName)
        fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        File modelFile = new File(dirName + fileName)
        File saveFile = potentialFunctions.versionFile(modelFile)

        FileWriter fw = new FileWriter(saveFile, false)
        BufferedWriter bw = new BufferedWriter(fw)

        int atomNumber = 1
        double[] xyz = new double[3]
        double occupancy = 1.0
        double bfactor = 1.0
        char altLoc = ' '
        char chain = 'A'
        String segID = "A"

        for (int i = 0; i < nAmino; i++) {
            int resID = i + 1
            String oneLetterResidue = sequence.charAt(i)
            String resName = convertToThreeLetter(oneLetterResidue)

            // Write N
            xyz[0] = Double.parseDouble(tokenizedLines[2][atomNumber]) / 100.0
            xyz[1] = Double.parseDouble(tokenizedLines[3][atomNumber]) / 100.0
            xyz[2] = Double.parseDouble(tokenizedLines[4][atomNumber]) / 100.0
            Atom atom = new Atom(atomNumber++, "N", altLoc, xyz, resName, resID, chain, occupancy, bfactor, segID)
            bw.write(toPDBAtomLine(atom))

            // Write CA
            xyz[0] = Double.parseDouble(tokenizedLines[2][atomNumber]) / 100.0
            xyz[1] = Double.parseDouble(tokenizedLines[3][atomNumber]) / 100.0
            xyz[2] = Double.parseDouble(tokenizedLines[4][atomNumber]) / 100.0
            atom = new Atom(atomNumber++, "CA", altLoc, xyz, resName, resID, chain, occupancy, bfactor, segID)
            bw.write(toPDBAtomLine(atom))

            // Write C
            xyz[0] = Double.parseDouble(tokenizedLines[2][atomNumber]) / 100.0
            xyz[1] = Double.parseDouble(tokenizedLines[3][atomNumber]) / 100.0
            xyz[2] = Double.parseDouble(tokenizedLines[4][atomNumber]) / 100.0
            atom = new Atom(atomNumber++, "C", altLoc, xyz, resName, resID, chain, occupancy, bfactor, segID)
            bw.write(toPDBAtomLine(atom))
        }

        bw.close()

        MolecularAssembly[] assemblies = [potentialFunctions.open(saveFile.getAbsolutePath())]
        activeAssembly = assemblies[0]
        PDBFilter pdbFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties())
        pdbFilter.writeFile(saveFile, false, false)

        return this
    }

    /**
     * This method takes in a one letter amino acid and converts it to the three letter amino acid code.
     * @param res The one letter amino acid code.
     * @return The three letter amino acid code.
     */
    private static String convertToThreeLetter(String res) {
        ResidueEnumerations.AminoAcid3 aminoAcid3 = ResidueEnumerations.getAminoAcid3From1(res)
        return aminoAcid3.toString()
    }

}
