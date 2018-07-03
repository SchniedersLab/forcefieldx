/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms.cli;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Logger;

import ffx.algorithms.RotamerOptimization;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;

import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that use a many-body expansion for global optimization.
 *
 * @author Michael J. Schnieders
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class ManyBodyOptions {

    private static final Logger logger = Logger.getLogger(ManyBodyOptions.class.getName());

    /**
     * PARAMETERS SHARED BY ALL ALGORITHMS. 
     */

    /**
     * -L or --library Choose either Ponder and Richards (1) or Richardson (2) rotamer library.
     */
    @Option(names = {"-L", "--library"}, paramLabel = "1",
            description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    int library = 1;

    /**
     * -Ln or --libraryNucleic Choose a nucleic acid library: currently only Richardson available.
     */
    @Option(names = {"--Ln", "--libraryNucleic"}, paramLabel = "Richardson",
            description = "Nucleic acid library to select: [Richardson]")
    String naLibraryName = "Richardson";
    //RotamerLibrary.NucleicAcidLibrary naLibrary = RotamerLibrary.NucleicAcidLibrary.RICHARDSON;

    /**
     * -a or --algorithm Choices are independent residues (1), all with rotamer elimination (2),
     * all brute force (3), sliding window (4), or box optimization (5).
     */
    @Option(names = {"-a", "--algorithm"}, paramLabel = "2",
            description = "Algorithm: independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)")
    int algorithm = 2;

    /**
     * --dee or --deadEnd Use dead-end elimination criteria instead of Goldstein criteria.
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
     * -s or --start Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.
     */
    @Option(names = {"-s", "--start"}, paramLabel = "-1",
            description = "Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.")
    int start = -1;

    /**
     * --fi or --final Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.
     */
    @Option(names = {"--fi", "--final"}, paramLabel = "-1",
            description = "Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.")
    int finish = -1;

    /**
     * --tC or --twoBodyCutoff Cutoff distance for two-body interactions.
     */
    @Option(names = {"--tC", "--twoBodyCutoff"}, paramLabel = "-1.0",
            description = "Cutoff distance for two body interactions.")
    double twoBodyCutoff = -1.0;

    /**
     * -T or --threeBody Include 3-Body interactions in the elimination criteria.
     */
    @Option(names = {"-T", "--threeBody"},
            description = "Include 3-Body interactions in the elimination criteria.")
    boolean threeBody = false;

    /**
     * --thC or --threeBodyCutoff Cutoff distance for three-body interactions.
     */
    @Option(names = {"--thC", "--threeBodyCutoff"}, paramLabel = "9.0",
            description = "Cutoff distance for three-body interactions.")
    double threeBodyCutoff = 9.0;

    /**
     * --pr or --prune Prune no clashes (0), only single clashes (1), or all clashes (2).
     */
    @Option(names = {"--pr", "--prune"}, paramLabel = "2",
            description = "Prune no clashes (0), only single clashes (1), or all clashes (2)")
    int prune = 2;

    /**
     * -x or --all Optimize all residues beginning from the passed value (overrides other options);
     * for box optimization, optimizes all boxes beginning from the passed index.
     */
    @Option(names = {"-x", "--all"}, paramLabel = "-1",
            description = "Optimize all residues beginning from the passed value (overrides other options); for box optimization, optimizes all boxes beginning from the passed index.")
    int all = -1;

    /**
     * --eR or --energyRestart Load energy restart file from a previous run (requires that all parameters are the same).
     */
    @Option(names = {"--eR", "--energyRestart"}, paramLabel = "none",
            description = "Load energy restart file from a previous run (requires that all parameters are the same).")
    String energyRestart = "none";

    /**
     * -v or --verbose Prints beginning and default-conformation energies.
     */
    @Option(names = {"-v", "--verbose"},
            description = "Prints beginning and default-conformation energies.")
    boolean verbose = false;

    /**
     * -o or --original Do not include starting coordinates as their own rotamer.
     */
    @Option(names = {"-O", "--original"},
            description = "Do not include starting coordinates as their own rotamer.")
    boolean original = false;

    /**
     * -E or --decompose Print energy decomposition for the input structure (no optimization).
     */
    @Option(names = {"-E", "--decompose"},
            description = "Print energy decomposition for the input structure (no optimization!).")
    boolean decompose = false;

    /**
     * --lR or --listResidues Choose a list of individual residues to optimize (eg. A11,A24,B40).
     */
    @Option(names = {"--lR", "--listResidues"}, paramLabel = "none",
            description = "Choose a list of individual residues to optimize (eg. A11,A24,B40).")
    String listResidues = "none";

    /**
     * --sO or --sequence Choose a list of individual residues to sequence optimize (example: A2.A3.A5).
     */
    // @Option(names = {"--sO", "--sequence"}, paramLabel = "none",
    //        description = "Choose a list of individual residues to sequence optimize (example: A2.A3.A5)")
    // String sequence = "none";

    /**
     * --tO or --titrationOptimization Optimize the titration states for a list of residues (example: H2.H3.H5).
     */
    // @Option(names = {"--tO", "--titrationOptimization"}, paramLabel = "none",
    //        description = "Optimize the titration states for a list of residues (example: H2.H3.H5).")
    // String titrationOptimization = "none";

    /**
     * --mC or --monteCarlo Follow elimination criteria with 'n' Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.
     */
    @Option(names = {"--mC", "--monteCarlo"}, paramLabel = "-1",
            description = "Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.")
    int monteCarlo = -1;

    /**
     * -out or --output Save eliminated singles and eliminated pairs to a text file (global and box optimization).
     */
    @Option(names = {"--out", "--output"}, paramLabel = "none",
            description = "Save eliminated singles and eliminated pairs to a text file.")
    boolean saveOutput = false;


    /**
     * PARAMETERS SPECIFIC TO SLIDING WINDOW.
     */

    /**
     * --window Size of the sliding window with respect to adjacent residues (default = 7).
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
     * -fR or --forceResidues Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.
     */
    @Option(names = {"--fR", "--forceResidues"}, paramLabel = "-1,-1",
            description = "Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.")
    String forceResidues = "-1,-1";


    /**
     * PARAMETERS SPECIFIC TO SLIDING BOX.
     */

    /**
     * -nB or --numBoxes Specify number of boxes along X, Y, and Z (default: '3,3,3').
     */
    @Option(names = {"--nB", "--numBoxes"}, paramLabel = "3,3,3",
            description = "Specify number of boxes along X, Y, and Z (default: 3,3,3)")
    String numBoxes = "3,3,3";

    /**
     * -bB or --boxBorderSize Extent of overlap between optimization boxes in Angstroms (default: 3.0).
     */
    @Option(names = {"--bB", "--boxBorderSize"}, paramLabel = "3.0",
            description = "Extent of overlap between optimization boxes in Angstroms.")
    double boxBorderSize = 3.0;

    /**
     * -bL or --approxBoxLength Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).
     * Box sizes are rounded up to make a whole number of boxes along each axis (default of 0 disables this function).
     */
    @Option(names = {"--bL", "--approxBoxLength"}, paramLabel = "0.0",
            description = "Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).")
    double approxBoxLength = 0.0;

    /**
     * -bC or --boxInclusionCriterion Criterion to use for adding residues to boxes.
     * (1) uses C alpha only (N1/9 for nucleic acids)
     * (2) uses any atom.
     * (3) uses any rotamer
     */
    @Option(names = {"--bC", "--boxInclusionCriterion"}, paramLabel = "1",
            description = "Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer")
    int boxInclusionCriterion = 1;

    RotamerOptimization rotamerOptimization;
    RotamerLibrary rLib;

    int allStartResID;
    int boxStart;
    int boxEnd;
    int[] numXYZBoxes;
    int forceResiduesStart;
    int forceResiduesEnd;

    public void initRotamerOptimization(RotamerOptimization rotamerOptimization, MolecularAssembly activeAssembly) {
        this.rotamerOptimization = rotamerOptimization;

        /**
         * Fully initialize the Rotamer Library.
         */
        rLib = RotamerLibrary.getDefaultLibrary();
        if (library == 1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
        } else {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        }

        boolean useOrigCoordsRotamer = !original;
        if (decompose) {
            useOrigCoordsRotamer = true;
        }
        rLib.setUseOrigCoordsRotamer(useOrigCoordsRotamer);
        rotamerOptimization.setDecomposeOriginal(decompose);

        setSelection();
        setForcedResidue();
        setResidues(activeAssembly);
        setRotOptProperties();
    }

    /**
     * Set allStartResID, boxStart and boxEnd
     */
    private void setSelection() {
        /**
         * Chain, Residue and/or Box selections.
         */
        allStartResID = all;
        // Internal machinery indexed 0 to (n-1)
        boxStart = start - 1;
        boxEnd = finish - 1;

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
                if (boxStart < 0) {
                    logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
                    return;
                }
            } else {
                if (boxStart < 0 || (boxEnd > -1 && boxEnd < boxStart)) {
                    logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
                    return;
                }
            }
        }


        /**
         * Box optimization options.
         */
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
                    logger.info("Double input to nB truncated to integer.");
                } else if (boxNumInput.hasNext()) {
                    logger.info("Non-numeric input to nB discarded");
                    boxNumInput.next();
                } else {
                    logger.info("Insufficient input to nB. Non-input values assumed either equal to X or default to 3");
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
                    logger.info("Input of zero to nB reset to default of three.");
                } else if (numXYZBoxes[i] < 0) {
                    numXYZBoxes[i] = -1 * numXYZBoxes[i];
                    logger.info("Input of negative number to nB reset to positive number");
                }
            }
        }
    }

    public void setForcedResidue() {
        /**
         * Force residues.
         */
        forceResiduesStart = -1;
        forceResiduesEnd = -1;

        List<String> resList = new ArrayList<>();
        if (!listResidues.equalsIgnoreCase("none")) {
            String tok[] = listResidues.split("\\.");
            for (String t : tok) {
                logger.info(" Adding " + t);
                resList.add(t);
            }
        }

        /**
         * Evaluate forced residues for the sliding window algorithm
         */
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
            if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for boxes " + (boxStart + 1) + " to " + (boxEnd + 1));
            } else {
                logger.info("\n Evaluating rotamers for all boxes beginning at " + (boxStart + 1));
            }
        }
    }

    public void setResidues(MolecularAssembly activeAssembly) {

        List<String> resList = new ArrayList<>();
        if (!listResidues.equalsIgnoreCase("none")) {
            String tok[] = listResidues.split("\\.");
            for (String t : tok) {
                logger.info(" Adding " + t);
                resList.add(t);
            }
        }

        int counter = 1;
        if (algorithm != 5) {
            if (allStartResID > 0) {
                ArrayList<Residue> residueList = new ArrayList<Residue>();
                Polymer[] polymers = activeAssembly.getChains();
                int nPolymers = polymers.length;
                for (int p = 0; p < nPolymers; p++) {
                    Polymer polymer = polymers[p];
                    ArrayList<Residue> residues = polymer.getResidues();
                    int nResidues = residues.size();
                    for (int i = 0; i < nResidues; i++) {
                        Residue residue = residues.get(i);
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
            boolean ignoreNA = false;
            String ignoreNAProp = System.getProperty("ignoreNA");
            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
                ignoreNA = true;
            }
            ArrayList<Residue> residueList = new ArrayList<>();
            Polymer[] polymers = activeAssembly.getChains();
            int nPolymers = polymers.length;
            for (int p = 0; p < nPolymers; p++) {
                Polymer polymer = polymers[p];
                ArrayList<Residue> residues = polymer.getResidues();
                int nResidues = residues.size();
                for (int i = 0; i < nResidues; i++) {
                    Residue residue = residues.get(i);
                    if (ignoreNA && residue.getResidueType() == ResidueType.NA) {
                        continue;
                    }
                    Rotamer[] rotamers = residue.getRotamers(rLib);
                    if (rotamers != null) {
                        int nrot = rotamers.length;
                        if (nrot == 1) {
                            RotamerLibrary.applyRotamer(residue, rotamers[0]);
                        } else if (nrot > 1) {
                            residueList.add(residue);
                        }
                    }
                    counter++;
                }
            }
            rotamerOptimization.setResidues(residueList);
            rotamerOptimization.setBoxStart(boxStart);
            if (allStartResID != -1) {
                rotamerOptimization.setBoxEnd(boxEnd);
            }
        }
    }

    public void setRotOptProperties() {
        // General
        rotamerOptimization.setTwoBodyCutoff(twoBodyCutoff);
        rotamerOptimization.setThreeBodyCutoff(threeBodyCutoff);
        rotamerOptimization.setThreeBodyEnergy(threeBody);
        rotamerOptimization.setUseGoldstein(!dee);
        //rotamerOptimization.setRevert(!noRevert);
        rotamerOptimization.setPruning(prune);
        rotamerOptimization.setDistanceCutoff(cutoff);
        rotamerOptimization.setVerboseEnergies(verbose);
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

    public void saveEliminatedRotamers() throws IOException {
        if (saveOutput) {
            rotamerOptimization.outputEliminated();
        }
    }


    /**
     List<String> sequenceOptimizationList = new ArrayList<>();
     if(!manyBody.sequence.equalsIgnoreCase("none"))

     {
     def tok = manyBody.sequence.tokenize('.');
     for (String t : tok) {
     logger.info(" Sequence optimizing " + t);
     sequenceOptimizationList.add(t);
     }
     // Get the ForceFieldString properly instead.
     //if (System.getProperty("relative-solvation") == null) {
     //System.setProperty("relative-solvation", "AUTO");
     //}
     }

     List<String> titrationOptimizationList = new ArrayList<>();
     if(!manyBody.titrationOptimization.equalsIgnoreCase("none"))

     {
     def tok = manyBody.titrationOptimization.tokenize('.');
     for (String t : tok) {
     logger.info(" Optimizing protonation states of " + t);
     titrationOptimizationList.add(t);
     }
     }

     if(!manyBody.sequence.equalsIgnoreCase("none"))

     {
     for (String s : sequenceOptimizationList) {
     Character chainID = s.charAt(0);
     int num = Integer.parseInt(s.substring(1));
     for (int i = 0; i < residueList.size(); i++) {
     Residue res = residueList.get(i);
     if (res.getChainID() == chainID && res.getResidueNumber() == num) {
     MultiResidue multiRes = new MultiResidue(res, activeAssembly.getForceField(),
     activeAssembly.getPotentialEnergy());
     for (Polymer polymer : activeAssembly.getChains()) {
     if (polymer.getChainID() == chainID) {
     logger.info(String.format(" Adding multiresidue %s to chain %c.", multiRes, chainID));
     polymer.addMultiResidue(multiRes);
     }
     }
     for (CommonAminoAcid3 aa : CommonAminoAcid3.values()) {
     if (aa.toString().equals("PRO") || aa.toString().equals("GLY")) {
     continue;
     }
     if (!aa.toString().equalsIgnoreCase(res.getName())) {
     logger.info(String.format(" Adding %s to residue %s.", aa.toString(), multiRes.toString()));
     multiRes.addResidue(new Residue(aa.toString(), res.getResidueNumber(), ResidueType.AA));
     }
     }
     multiRes.setActiveResidue(res);
     activeAssembly.getPotentialEnergy().reInit();
     residueList.remove(i);
     residueList.add(i, multiRes);
     }
     }
     }
     }

     if(!manyBody.titrationOptimization.equalsIgnoreCase("none"))

     {
     ArrayList<Residue> titrating = new ArrayList<>();
     for (String s : titrationOptimizationList) {
     Character chainID = s.charAt(0);
     int num = Integer.parseInt(s.substring(1));
     for (int i = 0; i < residueList.size(); i++) {
     Residue res = residueList.get(i);
     if (res.getChainID() == chainID && res.getResidueNumber() == num) {
     titrating.add(res);
     }
     }
     }
     rotamerOptimization.titrationSetResidues(titrating);
     }
     */

}
