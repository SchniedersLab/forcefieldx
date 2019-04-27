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
package ffx.algorithms.cli;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;

import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that use a many-body expansion
 * for global optimization.
 *
 * @author Michael J. Schnieders
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class ManyBodyOptions {

    private static final Logger logger = Logger.getLogger(ManyBodyOptions.class.getName());

    // Parameters for algorithms.

    /**
     * -L or --library Choose either Ponder and Richards (1) or Richardson (2)
     * rotamer library.
     */
    @Option(names = {"-L", "--library"}, paramLabel = "2",
            description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    int library = 2;

    /**
     * -Ln or --libraryNucleic Choose a nucleic acid library: currently only
     * Richardson available.
     */
    @Option(names = {"--Ln", "--libraryNucleic"}, paramLabel = "Richardson",
            description = "Nucleic acid library to select: [Richardson]")
    String naLibraryName = "Richardson";

    /**
     * -a or --algorithm Choices are independent residues (1), all with rotamer
     * elimination (2), all brute force (3), sliding window (4), or box
     * optimization (5).
     */
    @Option(names = {"-a", "--algorithm"}, paramLabel = "0",
            description = "Algorithm: default automatic settings (0), independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)")
    int algorithm = 0;

    /**
     * --dee or --deadEnd Use dead-end elimination criteria instead of Goldstein
     * criteria.
     */
    @Option(names = {"--dee", "--deadEnd"}, paramLabel = "false",
            description = "Use dead-end elimination criteria instead of Goldstein criteria.")
    boolean dee = false;

    /**
     * --ch or --chain Single character chain ID of the residues to optimize.
     */
    @Option(names = {"--ch", "--chain"}, paramLabel = "-1",
            description = "Single character chain ID of the residues to optimize.")
    String chain = "-1";

    /**
     * -s or --start Starting residue to perform the optimization on (-1 exits).
     * For box optimization, first box to optimize.
     */
    @Option(names = {"-s", "--start"}, paramLabel = "-1",
            description = "Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.")
    int start = -1;

    /**
     * --fi or --final Final residue to perform the optimization on (-1 exits).
     * For box optimization, final box to optimize.
     */
    @Option(names = {"--fi", "--final"}, paramLabel = "-1",
            description = "Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.")
    int finish = -1;

    /**
     * --tC or --twoBodyCutoff Cutoff distance for two-body interactions.
     */
    @Option(names = {"--tC", "--twoBodyCutoff"}, paramLabel = "3.0",
            description = "Cutoff distance for two body interactions.")
    double twoBodyCutoff = 3.0;

    /**
     * -T or --threeBody Include 3-Body interactions in the elimination
     * criteria.
     */
    @Option(names = {"-T", "--threeBody"},
            description = "Include 3-Body interactions in the elimination criteria.")
    boolean threeBody = false;

    /**
     * --thC or --threeBodyCutoff Cutoff distance for three-body interactions.
     */
    @Option(names = {"--thC", "--threeBodyCutoff"}, paramLabel = "3.0",
            description = "Cutoff distance for three-body interactions.")
    double threeBodyCutoff = 3.0;

    /**
     * --pr or --prune Prune no clashes (0), only single clashes (1), or all
     * clashes (2).
     */
    @Option(names = {"--pr", "--prune"}, paramLabel = "1",
            description = "Prune no clashes (0), only single clashes (1), or all clashes (2)")
    int prune = 1;

    /**
     * -x or --all Optimize all residues beginning from the passed value
     * (overrides other options); for box optimization, optimizes all boxes
     * beginning from the passed index. Default is to optimize all residues.
     */
    @Option(names = {"-x", "--all"}, paramLabel = "-1",
            description = "Optimize all residues beginning from the passed value (overrides other options); for box optimization, optimizes all boxes beginning from the passed index.")
    int all = -1;

    /**
     * -z or --revert Revert unfavorable changes.
     */
    @Option(names = {"-z", "--revert"}, paramLabel = "-1",
            description = "Revert unfavorable changes.")
    boolean revert = true;

    /**
     * --eR or --energyRestart Load energy restart file from a previous run
     * (requires that all parameters are the same).
     */
    @Option(names = {"--eR", "--energyRestart"}, paramLabel = "none",
            description = "Load energy restart file from a previous run (requires that all parameters are the same).")
    String energyRestart = "none";

    /**
     * -o or --noOriginal Do not include starting coordinates as their own
     * rotamer.
     */
    @Option(names = {"-O", "--noOriginal"},
            description = "Do not include starting coordinates as their own rotamer.")
    boolean noOriginal = false;

    /**
     * -E or --decompose Print energy decomposition for the input structure (no
     * optimization).
     */
    @Option(names = {"-E", "--decompose"},
            description = "Print energy decomposition for the input structure (no optimization!).")
    boolean decompose = false;

    /**
     * --lR or --listResidues Choose a list of individual residues to optimize
     * (eg. A11,A24,B40).
     */
    @Option(names = {"--lR", "--listResidues"}, paramLabel = "none",
            description = "Choose a list of individual residues to optimize (eg. A11,A24,B40).")
    String listResidues = "none";

//    /**
//     * --sO or --sequence Choose a list of individual residues to sequence
//     * optimize (example: A2.A3.A5).
//     */
//    @Option(names = {"--sO", "--sequence"}, paramLabel = "none",
//            description = "Choose a list of individual residues to sequence optimize (example: A2.A3.A5)")
//    String sequence = "none";
//    /**
//     * --tO or --titrationOptimization Optimize the titration states for a list
//     * of residues (example: H2.H3.H5).
//     */
//    @Option(names = {"--tO", "--titrationOptimization"}, paramLabel = "none",
//            description = "Optimize the titration states for a list of residues (example: H2.H3.H5).")
//    String titrationOptimization = "none";

    /**
     * --mC or --monteCarlo Follow elimination criteria with 'n' Monte Carlo
     * steps, or enumerate all remaining conformations, whichever is smaller.
     */
    @Option(names = {"--mC", "--monteCarlo"}, paramLabel = "-1",
            description = "Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.")
    int monteCarlo = -1;

    /**
     * -out or --output Save eliminated singles and eliminated pairs to a text
     * file (global and box optimization).
     */
    @Option(names = {"--out", "--output"}, paramLabel = "none",
            description = "Save eliminated singles and eliminated pairs to a text file.")
    boolean saveOutput = false;

    // Parameters for the sliding window algorithm.

    /**
     * --window Size of the sliding window with respect to adjacent residues
     * (default = 7).
     */
    @Option(names = {"--window"}, paramLabel = "7",
            description = "Size of the sliding window with respect to adjacent residues.")
    int window = 7;

    /**
     * --increment Sliding window increment (default = 3).
     */
    @Option(names = {"--increment"}, paramLabel = "3",
            description = "Sliding window increment.")
    int increment = 3;

    /**
     * --radius The sliding window cutoff radius (Angstroms).
     */
    @Option(names = {"--radius"}, paramLabel = "2.0",
            description = "The sliding window cutoff radius (Angstroms).")
    double cutoff = 2.0;

    /**
     * --clashThreshold The threshold for pruning clashes. If two self-energies on the same residue have an
     * energy difference greater than 25 kcal/mol, the high energy rotamers get pruned.
     */
    @Option(names = {"--clashThreshold"}, paramLabel = "25.0",
            description = "The threshold for pruning clashes.")
    double clashThreshold = 25.0;

    /**
     * --pairClashThreshold The threshold for pruning clashes. If two pair-energies on the same residues have an
     * energy difference greater than 25 kcal/mol, the high energy rotamers get pruned.
     */
    @Option(names = {"--pairClashThreshold"}, paramLabel = "25.0",
            description = "The threshold for pruning pair clashes.")
    double pairClashThreshold = 25.0;

    /**
     * -fR or --forceResidues Force residues in this range to be considered for
     * sliding window radii, regardless of whether they lack rotamers.
     */
    @Option(names = {"--fR", "--forceResidues"}, paramLabel = "-1,-1",
            description = "Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.")
    String forceResidues = "-1,-1";

    // Parameters for the sliding box algorithm.

    /**
     * -nB or --numBoxes Specify number of boxes along X, Y, and Z (default:
     * '3,3,3').
     */
    @Option(names = {"--nB", "--numBoxes"}, paramLabel = "3,3,3",
            description = "Specify number of boxes along X, Y, and Z (default: 3,3,3)")
    String numBoxes = "3,3,3";

    /**
     * -bB or --boxBorderSize Extent of overlap between optimization boxes in
     * Angstroms (default: 0.0).
     */
    @Option(names = {"--bB", "--boxBorderSize"}, paramLabel = "0.0",
            description = "Extent of overlap between optimization boxes in Angstroms.")
    double boxBorderSize = 0.0;

    /**
     * -bL or --approxBoxLength Approximate side lengths of boxes to be
     * constructed (over-rides numXYZBoxes). Box sizes are rounded up to make a
     * whole number of boxes along each axis (default of 0 disables this
     * function).
     */
    @Option(names = {"--bL", "--approxBoxLength"}, paramLabel = "20.0",
            description = "Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).")
    double approxBoxLength = 20.0;

    /**
     * -bC or --boxInclusionCriterion Criterion to use for adding residues to
     * boxes. (1) uses C alpha only (N1/9 for nucleic acids) (2) uses any atom.
     * (3) uses any rotamer
     */
    @Option(names = {"--bC", "--boxInclusionCriterion"}, paramLabel = "1",
            description = "Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer")
    int boxInclusionCriterion = 1;

    private RotamerOptimization rotamerOptimization;
    private RotamerLibrary rLib;
    private int allStartResID;
    private int boxStart;
    private int boxEnd;
    private int[] numXYZBoxes;
    private int forceResiduesStart;
    private int forceResiduesEnd;

    /**
     * <p>initRotamerOptimization.</p>
     *
     * @param rotamerOptimization a {@link RotamerOptimization} object.
     * @param activeAssembly      a {@link ffx.potential.MolecularAssembly} object.
     */
    public void initRotamerOptimization(RotamerOptimization rotamerOptimization, MolecularAssembly activeAssembly) {
        this.rotamerOptimization = rotamerOptimization;

        boolean useOrigCoordsRotamer = !noOriginal;
        if (decompose) {
            useOrigCoordsRotamer = true;
        }

        String rotamerFileName = activeAssembly.getFile().getName();
        rotamerFileName = FilenameUtils.removeExtension(rotamerFileName);
        rotamerFileName = rotamerFileName + ".rot";
        File rotFile = new File(rotamerFileName);

        if (rotFile.exists()) {
            logger.info(" EXPERIMENTAL: Using rotamer file " + rotamerFileName);
            rLib = new RotamerLibrary(RotamerLibrary.ProteinLibrary.None, false);
            try {
                RotamerLibrary.readRotFile(rotFile, activeAssembly);
            } catch (IOException iox) {
                logger.severe(String.format(" Exception in parsing rotamer file: %s\n%s", iox, Utilities.stackTraceToString(iox)));
            }
        } else {
            rLib = new RotamerLibrary(RotamerLibrary.ProteinLibrary.intToProteinLibrary(library), useOrigCoordsRotamer);
        }

        rotamerOptimization.setRotamerLibrary(rLib);
        rotamerOptimization.setSingletonClashThreshold(clashThreshold);
        rotamerOptimization.setPairClashThreshold(pairClashThreshold);
        rotamerOptimization.setDecomposeOriginal(decompose);

        if (algorithm == 0) {
            setAlgorithm(activeAssembly);
        }
        setSelection();
        setForcedResidue();
        setResidues(activeAssembly);
        setRotOptProperties();
    }

    /**
     * This method sets the algorithm by default. If no parameters are given,
     * the default algorithm value 0. When the default algorithm is 0, a
     * specific algorithm number (1-5) needs to be assigned. This method assigns
     * the default algorithm based on variation in input parameters (start,
     * finish, all, chain, etc.) and then determines how many amino acids are to
     * be optimized. If more than 100 amino acids are to be optimized, the
     * algorithm is set to use box optimization. If fewer than 100 amino acids
     * are to be optimized, the algorithm is set to use global optimization.
     *
     * @param activeAssembly The protein to be optimized.
     * @return The value that the algorithm should be set to since no 1-5 value
     * was assigned by the user. Either set to 2-global or 5-box optimization.
     */
    private int setAlgorithm(MolecularAssembly activeAssembly) {

        setStartAndEndDefault();

        if (allStartResID > 0) {
            Polymer[] polymers = activeAssembly.getChains();
            int nResidues = 0;
            for (Polymer polymer : polymers) {
                ArrayList<Residue> residues = polymer.getResidues();
                nResidues = residues.size() + nResidues;
            }
            if (nResidues > 100) {
                algorithm = 5;
            } else {
                algorithm = 2;
            }
        } else if (!listResidues.equalsIgnoreCase("none")) {
            List<String> resList = new ArrayList<>();
            addListResidues(resList);
            if (resList.size() > 100) {
                algorithm = 5;
            } else {
                algorithm = 2;
            }
        } else if (!chain.equalsIgnoreCase("-1")) {
            algorithm = startFinishDifference();
        } else if (start > 0 && finish > 0) {
            algorithm = startFinishDifference();
        }
        return algorithm;
    }

    /**
     * This method calculates the difference between the start and finish
     * variables and it returns the algorithm that should be used by default. If
     * more than 100 side-chains are to be optimized, the default algorithm is
     * box optimization (algorithm = 5). If fewer than 100 residues are to be
     * optimized, the default algorithm is global optimization (algorithm = 2).
     *
     * @return The default algorithm number.
     */
    private int startFinishDifference() {
        if (finish - start > 100) {
            return 5;
        } else {
            return 2;
        }
    }

    /**
     * This method sets the start and finish points by default. If no parameters
     * are given as input, the default behavior is to begin the optimization at
     * the first available amino acid and end the optimization at the last amino
     * acid. The default algorithm (i.e. box versus global) depends on how many
     * amino acids will be optimized (general cutoff of 100 amino acids).
     */
    private void setStartAndEndDefault() {
        if (start < 0 && finish < 0 && all < 0) {
            if (!listResidues.equalsIgnoreCase("none")) {
                allStartResID = -1;
                boxStart = 0;
                boxEnd = -1;
            } else {
                allStartResID = 1;
                boxStart = start - 1;
                boxEnd = finish - 1;
            }
        } else {
            allStartResID = all;
            boxStart = start - 1;
            boxEnd = finish - 1;
        }
    }

    /**
     * Set allStartResID, boxStart and boxEnd
     */
    private void setSelection() {

        // Chain, Residue and/or Box selections.
        // Internal machinery indexed 0 to (n-1)
        setStartAndEndDefault();

        if (algorithm != 5) {
            // Not Box optimization.
            if (allStartResID < 1 && listResidues.equalsIgnoreCase("none")) {
                if (finish < start || start < 0 || finish < 0) {
                    logger.warning(" FFX shutting down: no residues specified for optimization.");
                    return;
                }
            }
        } else {
            // Box optimization.
            if (allStartResID > 0) {
                // Internal machinery indexed 0 to (n-1)
                boxStart = allStartResID - 1;
            } else {
                if (boxStart < 0 || (boxEnd > -1 && boxEnd < boxStart)) {
                    logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
                    return;
                }
            }
        }

        // Box optimization options.
        numXYZBoxes = new int[3];
        if (algorithm == 5) {
            String input = numBoxes;
            Scanner boxNumInput = new java.util.Scanner(input);
            boxNumInput.useDelimiter(",");
            int inputLoopCounter = 0;
            numXYZBoxes[0] = 3; // Default
            while (inputLoopCounter < numXYZBoxes.length) {
                if (boxNumInput.hasNextInt()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextInt();
                    inputLoopCounter++;
                } else if (boxNumInput.hasNextDouble()) {
                    numXYZBoxes[inputLoopCounter] = (int) Math.floor(boxNumInput.nextDouble());
                    inputLoopCounter++;
                    logger.info(" Double input to nB truncated to integer.");
                } else if (boxNumInput.hasNext()) {
                    logger.info(" Non-numeric input to nB discarded");
                    boxNumInput.next();
                } else {
                    logger.info(" Insufficient input to nB. Non-input values assumed either equal to X or default to 3");
                    break;
                }
            }
            boxNumInput.close();
            for (int i = inputLoopCounter; i < numXYZBoxes.length; i++) {
                numXYZBoxes[i] = numXYZBoxes[0];
            }
            for (int i = 0; i < numXYZBoxes.length; i++) {
                if (numXYZBoxes[i] == 0) {
                    numXYZBoxes[i] = 3;
                    logger.info(" Input of zero to nB reset to default of three.");
                } else if (numXYZBoxes[i] < 0) {
                    numXYZBoxes[i] = -1 * numXYZBoxes[i];
                    logger.info(" Input of negative number to nB reset to positive number");
                }
            }
        }
    }

    /**
     * <p>setForcedResidue.</p>
     */
    private void setForcedResidue() {
        // Force residues.
        forceResiduesStart = -1;
        forceResiduesEnd = -1;

        List<String> resList = new ArrayList<>();
        if (!listResidues.equalsIgnoreCase("none")) {
            String tok[] = listResidues.split(",");
            for (String t : tok) {
                logger.info(" Adding " + t);
                resList.add(t);
            }
        }

        // Evaluate forced residues for the sliding window algorithm
        if (algorithm == 4 && !forceResidues.equalsIgnoreCase("-1,-1")) {
            String input = forceResidues;
            Scanner frScan = new Scanner(input);
            frScan.useDelimiter(",");
            try {
                if (!frScan.hasNextInt()) {
                    frScan.next(); // Discards extra input to indicate a negative value of frStart.
                }
                forceResiduesStart = frScan.nextInt();
                forceResiduesEnd = frScan.nextInt();
            } catch (Exception ex) {
                logger.severe(String.format(" FFX shutting down: input to -fR could not be parsed as a pair of integers: %s", input));
            }
            if (forceResiduesStart > forceResiduesEnd) {
                logger.info(" Start of range higher than ending: start flipped with end.");
                int temp = forceResiduesStart;
                forceResiduesStart = forceResiduesEnd;
                forceResiduesEnd = temp;
            }
            if (forceResiduesEnd < 1) {
                logger.severe(String.format(" FFX shutting down: end range for -fR must be at least 1; input range %d to %d",
                        forceResiduesStart, forceResiduesEnd));
            }
        }

        if (algorithm != 5) {
            if (!listResidues.equalsIgnoreCase("none")) {
                String info = "\n Evaluating rotamers for residues ";
                for (String i : resList) {
                    info += String.format("%s, ", i);
                }
                logger.info(info);
            } else if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for residues " + start + " to " + finish);
            } else {
                logger.info("\n Evaluating rotamers for all residues beginning at " + allStartResID);
            }
        } else {
            if (!listResidues.equalsIgnoreCase("none")) {
                String info = "\n Evaluating rotamers for boxes with residues ";
                for (String i : resList) {
                    info += String.format("%s, ", i);
                }
                logger.info(info);
            } else if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for boxes " + (boxStart + 1) + " to " + (boxEnd + 1));
            } else {
                logger.info("\n Evaluating rotamers for all boxes beginning at " + (boxStart + 1));
            }
        }
    }

    /**
     * <p>addListResidues.</p>
     *
     * @param resList a {@link java.util.List} object.
     */
    private void addListResidues(List<String> resList) {
        if (!listResidues.equalsIgnoreCase("none")) {
            String tok[] = listResidues.split(",");
            for (String t : tok) {
                logger.info(" Adding " + t);
                resList.add(t);
            }
        }
    }

    /**
     * <p>setResidues.</p>
     *
     * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
     */
    public void setResidues(MolecularAssembly activeAssembly) {
        List<String> resList = new ArrayList<>();
        addListResidues(resList);

        int counter = 1;
        if (algorithm != 5) {
            if (allStartResID > 0) {
                ArrayList<Residue> residueList = new ArrayList<Residue>();
                Polymer[] polymers = activeAssembly.getChains();
                for (Polymer polymer : polymers) {
                    ArrayList<Residue> residues = polymer.getResidues();
                    for (Residue residue : residues) {
                        Rotamer[] rotamers = residue.getRotamers(rLib);
                        if (rotamers != null) {
                            int nrot = rotamers.length;
                            if (nrot == 1) {
                                RotamerLibrary.applyRotamer(residue, rotamers[0]);
                            }
                            if (counter >= allStartResID) {
                                residueList.add(residue);
                            }
                        } else if (!forceResidues.equalsIgnoreCase("none")) {
                            if (counter >= allStartResID && counter >= forceResiduesStart
                                    && counter <= forceResiduesEnd) {
                                residueList.add(residue);
                            }
                        }
                        counter++;
                    }
                }
                rotamerOptimization.setResidues(residueList);
            } else if (!listResidues.equalsIgnoreCase("none")) {
                ArrayList<Residue> residueList = new ArrayList<>();
                Polymer[] polymers = activeAssembly.getChains();
                int n = 0;
                for (String s : resList) {
                    Character chainID = s.charAt(0);
                    int i = Integer.parseInt(s.substring(1));
                    for (Polymer p : polymers) {
                        if (p.getChainID() == chainID) {
                            List<Residue> rs = p.getResidues();
                            for (Residue r : rs) {
                                if (r.getResidueNumber() == i) {
                                    residueList.add(r);
                                    Rotamer[] rotamers = r.getRotamers(rLib);
                                    if (rotamers != null) {
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                }
                rotamerOptimization.setResiduesIgnoreNull(residueList);
                if (n < 1) {
                    return;
                }
            } else if (!chain.equalsIgnoreCase("-1")) {
                rotamerOptimization.setResidues(chain, start, finish);
            } else {
                rotamerOptimization.setResidues(start, finish);
            }
        } else {
            ArrayList<Residue> residueList = new ArrayList<>();
            Polymer[] polymers = activeAssembly.getChains();
            boolean ignoreNA = false;
            String ignoreNAProp = System.getProperty("ignoreNA");
            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
                ignoreNA = true;
            }
            if (!listResidues.equalsIgnoreCase("none")) {
                int n = 0;
                for (String s : resList) {
                    Character chainID = s.charAt(0);
                    int i = Integer.parseInt(s.substring(1));
                    for (Polymer p : polymers) {
                        if (p.getChainID() == chainID) {
                            List<Residue> rs = p.getResidues();
                            for (Residue r : rs) {
                                if (ignoreNA && r.getResidueType() == ResidueType.NA) {
                                    continue;
                                }
                                if (r.getResidueNumber() == i) {
                                    residueList.add(r);
                                    Rotamer[] rotamers = r.getRotamers(rLib);
                                    if (rotamers != null) {
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                }
                rotamerOptimization.setResiduesIgnoreNull(residueList);
                if (n < 1) {
                    return;
                }
            } else {
                for (Polymer p : polymers) {
                    ArrayList<Residue> rs = p.getResidues();
                    for (Residue r : rs) {
                        if (ignoreNA && r.getResidueType() == ResidueType.NA) {
                            continue;
                        }
                        Rotamer[] rotamers = r.getRotamers(rLib);
                        if (rotamers != null) {
                            int nrot = rotamers.length;
                            if (nrot == 1) {
                                RotamerLibrary.applyRotamer(r, rotamers[0]);
                            } else if (nrot > 1) {
                                residueList.add(r);
                            }
                        }
                        counter++;
                    }
                }

//            boolean ignoreNA = false;
//            String ignoreNAProp = System.getProperty("ignoreNA");
//            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
//                ignoreNA = true;
//            }
//            ArrayList<Residue> residueList = new ArrayList<>();
//            Polymer[] polymers = activeAssembly.getChains();
//            int nPolymers = polymers.length;
//            for (int p = 0; p < nPolymers; p++) {
//                Polymer polymer = polymers[p];
//                ArrayList<Residue> residues = polymer.getResidues();
//
//                System.out.print("\nresidues:\n");
//
//                int nResidues = residues.size();
//                for (int i = 0; i < nResidues; i++) {
//                    Residue residue = residues.get(i);
//
//                    System.out.print(residue+"\n");
//
//                    if (ignoreNA && residue.getResidueType() == ResidueType.NA) {
//                        continue;
//                    }
//                    Rotamer[] rotamers = residue.getRotamers(rLib);
//                    if (rotamers != null) {
//                        int nrot = rotamers.length;
//                        if (nrot == 1) {
//                            RotamerLibrary.applyRotamer(residue, rotamers[0]);
//                        } else if (nrot > 1) {
//                            residueList.add(residue);
//                        }
//                    }
//                    counter++;
//                }
            }
            rotamerOptimization.setResidues(residueList);
            rotamerOptimization.setBoxStart(boxStart);
            if (finish > 0) {
                rotamerOptimization.setBoxEnd(boxEnd);
            }
        }
    }

    /**
     * Sets the standard values for properties in rotamer optimization.
     */
    private void setRotOptProperties() {
        // General
        rotamerOptimization.setTwoBodyCutoff(twoBodyCutoff);
        rotamerOptimization.setThreeBodyCutoff(threeBodyCutoff);
        rotamerOptimization.setThreeBodyEnergy(threeBody);
        rotamerOptimization.setUseGoldstein(!dee);
        rotamerOptimization.setRevert(revert);
        rotamerOptimization.setPruning(prune);
        rotamerOptimization.setDistanceCutoff(cutoff);
        boolean monteCarloBool = false;
        if (monteCarlo > 1) {
            monteCarloBool = true;
        }
        rotamerOptimization.setMonteCarlo(monteCarloBool, monteCarlo);

        File energyRestartFile = null;
        if (!energyRestart.equalsIgnoreCase("none")) {
            energyRestartFile = new File(energyRestart);
            rotamerOptimization.setEnergyRestartFile(energyRestartFile);
        }

        // Window
        if (algorithm == 4) {
            rotamerOptimization.setWindowSize(window);
            rotamerOptimization.setIncrement(increment);
            rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
        }

        // Box
        if (algorithm == 5) {
            if (approxBoxLength < 0) {
                logger.info(" Negative box length value changed to -1 * input.");
                approxBoxLength *= -1;
            }
            rotamerOptimization.setBoxBorderSize(boxBorderSize);
            rotamerOptimization.setApproxBoxLength(approxBoxLength);
            rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
            rotamerOptimization.setBoxInclusionCriterion(boxInclusionCriterion);
        }
    }

    /**
     * Gets the restart file created during rotamer optimization.
     *
     * @return The restart file.
     */
    public File getRestartFile() {
        return rotamerOptimization.getRestartFile();
    }

    /**
     * Gets the algorithm number.
     *
     * @return Integer representing the algorithm being run (i.e. global, box
     * optimization, window, etc.)
     */
    public int getAlgorithmNumber() {
        return algorithm;
    }

//    /**
//     * Saves all eliminated rotamers to an ouput file called "eliminated.csv"
//     * when the many body command is run with the following syntax and flags:
//     * ffxc ManyBody --out ... file.pdb &gt;&gt; file.log.
//     *
//     * @throws java.io.IOException Throws an exception when output is non piped to a log
//     *                             file. The --out flag relies on the presence of a log file where output is
//     *                             piped.
//     */
//    public void saveEliminatedRotamers() throws IOException {
//        if (saveOutput) {
//            rotamerOptimization.outputEliminated();
//        }
//    }

    /**
     * List<String> sequenceOptimizationList = new ArrayList<>();
     * if(!manyBody.sequence.equalsIgnoreCase("none"))
     *
     * {
     * def tok = manyBody.sequence.tokenize('.'); for (String t : tok) {
     * logger.info(" Sequence optimizing " + t);
     * sequenceOptimizationList.add(t); } // Get the ForceFieldString properly
     * instead. //if (System.getProperty("relative-solvation") == null) {
     * //System.setProperty("relative-solvation", "AUTO"); //} }
     *
     * List<String> titrationOptimizationList = new ArrayList<>();
     * if(!manyBody.titrationOptimization.equalsIgnoreCase("none"))
     *
     * {
     * def tok = manyBody.titrationOptimization.tokenize('.'); for (String t :
     * tok) { logger.info(" Optimizing protonation states of " + t);
     * titrationOptimizationList.add(t); } }
     *
     * if(!manyBody.sequence.equalsIgnoreCase("none"))
     *
     * {
     * for (String s : sequenceOptimizationList) { Character chainID =
     * s.charAt(0); int num = Integer.parseInt(s.substring(1)); for (int i = 0;
     * i < residueList.size(); i++) { Residue res = residueList.get(i); if
     * (res.getChainID() == chainID && res.getResidueNumber() == num) {
     * MultiResidue multiRes = new MultiResidue(res,
     * activeAssembly.getForceField(), activeAssembly.getPotentialEnergy()); for
     * (Polymer polymer : activeAssembly.getChains()) { if (polymer.getChainID()
     * == chainID) { logger.info(String.format(" Adding multiresidue %s to chain
     * %c.", multiRes, chainID)); polymer.addMultiResidue(multiRes); } } for
     * (CommonAminoAcid3 aa : CommonAminoAcid3.values()) { if
     * (aa.toString().equals("PRO") || aa.toString().equals("GLY")) { continue;
     * } if (!aa.toString().equalsIgnoreCase(res.getName())) {
     * logger.info(String.format(" Adding %s to residue %s.", aa.toString(),
     * multiRes.toString())); multiRes.addResidue(new Residue(aa.toString(),
     * res.getResidueNumber(), ResidueType.AA)); } }
     * multiRes.setActiveResidue(res);
     * activeAssembly.getPotentialEnergy().reInit(); residueList.remove(i);
     * residueList.add(i, multiRes); } } } }
     *
     * if(!manyBody.titrationOptimization.equalsIgnoreCase("none"))
     *
     * {
     * ArrayList<Residue> titrating = new ArrayList<>(); for (String s :
     * titrationOptimizationList) { Character chainID = s.charAt(0); int num =
     * Integer.parseInt(s.substring(1)); for (int i = 0; i < residueList.size();
     * i++) { Residue res = residueList.get(i); if (res.getChainID() == chainID
     * && res.getResidueNumber() == num) { titrating.add(res); } } }
     * rotamerOptimization.titrationSetResidues(titrating); }
     */
}
