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
package ffx.potential.parsers;

import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.io.PDBFileReader;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

/**
 * @author Jacob Litman
 */
public class SimplePDBMatcher {

    private static final Logger logger = Logger.getLogger(PDBFileMatcher.class.getName());
    private final File[] matchFiles;
    private final File[] sourceFiles;
    private FileDoublePair[] matchedSources;

    public SimplePDBMatcher(File[] matchFiles, File[] sourceFiles) {
        this.matchFiles = matchFiles;
        this.sourceFiles = sourceFiles;
        matchedSources = new FileDoublePair[matchFiles.length];
    }

    public void match() {
        PDBFileReader reader = new PDBFileReader();
        StructurePairAligner aligner = new StructurePairAligner();
        for (int i = 0; i < matchFiles.length; i++) {
            File matchFile = matchFiles[i];
            try {
                Structure matchStructure = reader.getStructure(matchFile);
                String matchName = matchFile.getName();
                matchedSources[i] = loopOverSources(reader, aligner, matchStructure, matchName);
            } catch (IOException ex) {
                logger.warning(String.format(" Matching failed for file %s", matchFile.getName()));
            }
        }
        for (int i = 0; i < matchFiles.length; i++) {
            logger.info(String.format(" Match file %s best source: %s at %11.7f A",
                    matchFiles[i].getName(), matchedSources[i].file.getName(), matchedSources[i].rmsd));
        }
    }

    public void matchParallel() {
        try {
            new ParallelTeam().execute(new ParallelRegion() {
                @Override
                public void run() throws Exception {
                    execute(0, sourceFiles.length, new IntegerForLoop() {
                        @Override
                        public void run(int lb, int ub) {
                            PDBFileReader reader = new PDBFileReader();
                            StructurePairAligner aligner = new StructurePairAligner();
                            for (int i = lb; i <= ub; i++) {
                                File matchFile = matchFiles[i];
                                try {
                                    Structure matchStructure = reader.getStructure(matchFiles[i]);
                                    String matchName = matchFile.getName();
                                    matchedSources[i] = loopOverSources(reader, aligner, matchStructure, matchName);
                                } catch (IOException ex) {
                                    logger.warning(String.format(" Matching failed for file %s", matchFile.getName()));
                                }
                            }
                        }
                    });
                }
            });
        } catch (Exception ex) {
            logger.severe(" Matching in parallel failed.");
        }
        for (int i = 0; i < matchFiles.length; i++) {
            logger.info(String.format(" Match file %s best source: %s at %11.7f A",
                    matchFiles[i].getName(), matchedSources[i].file.getName(), matchedSources[i].rmsd));
        }
    }

    private FileDoublePair loopOverSources(PDBFileReader reader, StructurePairAligner aligner, Structure matchStructure, String matchName) {
        File bestMatch = sourceFiles[0];
        double rmsd = Double.MAX_VALUE;
        for (File sourceFile : sourceFiles) {
            try {
                Structure sourceStructure = reader.getStructure(sourceFile);
                String sourceName = sourceFile.getName();
                double fileRMSD = checkRMSD(matchStructure, sourceStructure, aligner, matchName, sourceName);
                if (fileRMSD < rmsd) {
                    bestMatch = sourceFile;
                    rmsd = fileRMSD;
                }
            } catch (IOException ex) {
                logger.warning(String.format(" Source file %s could not be read", sourceFile.toString()));
            }
        }
        logger.info(String.format(" Minimum RMSD for file %s: $11.7f to file %s", matchName, rmsd, bestMatch.getName()));
        return new FileDoublePair(bestMatch, rmsd);
    }

    private double checkRMSD(Structure matchStructure, Structure sourceStructure, StructurePairAligner aligner, String matchName, String sourceName) {
        double bestRMSD = Double.MAX_VALUE;
        try {
            aligner.align(matchStructure, sourceStructure);
            AlternativeAlignment[] alignments = aligner.getAlignments();
            for (AlternativeAlignment alignment : alignments) {
                double alignRMSD = alignment.getRmsd();
                bestRMSD = Math.min(alignRMSD, bestRMSD);
            }
            logger.info(String.format(" Minimum RMSD for match %s source %s: %11.7f", matchName, sourceName, bestRMSD));
        } catch (StructureException ex) {
            logger.warning(String.format(" Structure exception during match %s source file %s", matchName, sourceName));
        }
        return bestRMSD;
    }

    private class FileDoublePair {

        public final File file;
        public final double rmsd;

        public FileDoublePair(File file, double rmsd) {
            this.file = file;
            this.rmsd = rmsd;
        }
    }
}
