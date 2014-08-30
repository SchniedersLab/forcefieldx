/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.parsers;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.io.PDBFileReader;

/**
 *
 * @author JacobLitman
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
