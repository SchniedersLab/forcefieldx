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

import ffx.potential.cli.PotentialScript
import org.apache.commons.io.FilenameUtils
import org.biojava.nbio.core.sequence.ProteinSequence
import org.biojava.nbio.core.sequence.io.FastaReaderHelper
import org.biojava.nbio.core.sequence.io.FastaWriterHelper
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * Fasta outputs a sub-sequence from a FASTA file.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc Fasta [options] &lt;filename.fasta&gt;
 */
@Command(description = " Fasta outputs a sub-sequence from a FASTA file.", name = "ffxc Fasta")
class Fasta extends PotentialScript {

  /**
   * -f or --firstResidue defines the first Fasta residue to keep (index of the first residue is 1).
   */
  @Option(names = ['-f', '--firstResidue'], paramLabel = "1", defaultValue = "1",
      description = 'Define the first Fasta residue to keep (index of the first residue is 1).')
  private int firstResidue = 1

  /**
   * -l or --lastResidue defines the last Fasta residue to keep (index of the last residue is n).
   */
  @Option(names = ['-l', '--lastResidue'], paramLabel = "-1", defaultValue = "-1",
      description = 'Define the last Fasta residue to keep (index of the last residue is n).')
  private int lastResidue = -1

  /**
   * The final argument should be an Fasta file.
   */
  @Parameters(arity = "1", paramLabel = "Fasta file",
      description = 'FASTA file.')
  List<String> filenames = null

  private File baseDir = null

  private ProteinSequence proteinSequence

  void setBaseDir(File baseDir) {
    this.baseDir = baseDir
  }

  /**
   * Execute the script.
   */
  @Override
  Fasta run() {

    if (!init()) {
      return null
    }

    if (filenames === null || filenames.size() == 0) {
      logger.info(helpString())
      return null
    }

    String fastaName = filenames.get(0)

    logger.info("\n Opening FASTA " + fastaName)

    File fastaFile = new File(fastaName)
    LinkedHashMap<String, ProteinSequence> fastaData =
        FastaReaderHelper.readFastaProteinSequence(fastaFile)
    ProteinSequence sequence = fastaData.values()[0]
    String seq = sequence.sequenceAsString
    int length = seq.length()
    logger.info(format("\n %s of length: %d\n %s", sequence.getOriginalHeader(), length, seq))

    if (firstResidue < 1 || firstResidue > length) {
      firstResidue = 1
    }
    if (lastResidue < firstResidue || lastResidue > length) {
      lastResidue = length
    }

    proteinSequence = new ProteinSequence(seq.substring(firstResidue - 1, lastResidue))
    proteinSequence.setOriginalHeader(sequence.getOriginalHeader())
    length = proteinSequence.length
    logger.info(format("\n New sequence from residue %d to residue %d is of length %d: \n %s",
        firstResidue, lastResidue, length, proteinSequence.toString()))

    File saveDir = baseDir
    fastaName = fastaFile.getAbsolutePath()
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(FilenameUtils.getFullPath(fastaName))
    }
    String dirName = saveDir.toString() + File.separator
    String fileName = FilenameUtils.getName(fastaName)
    fileName = FilenameUtils.removeExtension(fileName) + ".fasta"
    File modelFile = new File(dirName + fileName)
    File saveFile = potentialFunctions.versionFile(modelFile)

    Collection<ProteinSequence> proteinSequenceCollection = new ArrayList<>()
    proteinSequenceCollection.add(proteinSequence)

    logger.info(format("\n Saving new Fasta file to: %s", saveFile.getAbsolutePath()))
    FastaWriterHelper.writeProteinSequence(saveFile, proteinSequenceCollection)

    return this
  }

}
