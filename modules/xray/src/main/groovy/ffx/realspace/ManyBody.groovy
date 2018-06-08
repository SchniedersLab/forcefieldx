package realspace

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import edu.rit.pj.Comm

import ffx.algorithms.RotamerOptimization
import ffx.potential.bonded.MultiResidue
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Residue.ResidueType
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.realspace.RealSpaceData
import ffx.realspace.RealSpaceFile
import ffx.xray.CrystalReciprocalSpace
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize
import ffx.xray.parsers.DiffractionFile

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
class ManyBody extends Script {

    /**
     * Options for the ManyBody script.
     * <br>
     * Usage:
     * <br>
     * ffxc ManyBody [options] &lt;filename&gt;
     */
    class Options {

        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * -l or --library Choose either Ponder and Richards (1) or Richardson (2) rotamer library.
         */
        @Option(shortName = 'l', longName = 'library', defaultValue = '1',
                description = 'Ponder and Richards (1) or Richardson (2) rotamer library.')
        int library
        /**
         * -a or --algorithm Choices are independent residues (1), all with rotamer elimination (2),
         * all brute force (3), sliding window (4), or box optimization (5).
         */
        @Option(shortName = 'a', longName = 'algorithm', defaultValue = '2',
                description = 'Algorithm: independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)')
        int algorithm
        /**
         * -dee or --deadEnd Use dead-end elimination criteria instead of Goldstein criteria.
         */
        @Option(shortName = 'dee', longName = 'deadEnd',
                description = 'Use dead-end elimination criteria instead of Goldstein criteria.')
        boolean dee
        /**
         * -w or --window Size of the sliding window with respect to adjacent residues (default = 7).
         */
        @Option(shortName = 'w', longName = 'window', defaultValue = '7',
                description = 'Size of the sliding window with respect to adjacent residues.')
        int window
        /**
         * -i or --increment Sliding window increment (default = 3).
         */
        @Option(shortName = 'i', longName = 'increment', defaultValue = '3',
                description = 'Sliding window increment.')
        int increment
        /**
         * -r or --cutoff The sliding window cutoff radius (Angstroms).
         */
        @Option(shortName = 'r', longName = 'cutoff', defaultValue = '2.0',
                description = 'The sliding window cutoff radius (Angstroms).')
        double cutoff
        /**
         * -c or --chain Single character chain ID of the residues to optimize.
         */
        @Option(shortName = 'c', longName = 'chain', defaultValue = '-1',
                description = 'Single character chain ID of the residues to optimize.')
        String chain;
        /**
         * -s or --start Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.
         */
        @Option(shortName = 's', longName = 'start', defaultValue = '-1',
                description = 'Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.')
        int start
        /**
         * -f or --final Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.
         */
        @Option(shortName = 'f', longName = 'final', defaultValue = '-1',
                description = 'Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.')
        int last
        /**
         * -m or --minimize Minimize the final structure to the given RMS gradient (Kcal/mole/A); the default is not
         * to minimize.
         */
        @Option(shortName = 'm', longName = 'minimize', defaultValue = '0.0',
                description = 'Minimize the final structure to the given RMS gradient (Kcal/mole/A).')
        double minimize
        /**
         * -t or --threeBody Include 3-Body interactions in the elimination criteria.
         */
        @Option(shortName = 't', longName = 'threeBody',
                description = 'Include 3-Body interactions in the elimination criteria.')
        boolean threeBody
        /**
         * -p or --prune Prune no clashes (0), only single clashes (1), or all clashes (2).
         */
        @Option(shortName = 'p', longName = 'prune', defaultValue = '2',
                description = 'Prune no clashes (0), only single clashes (1), or all clashes (2)')
        int prune
        /**
         * -x or --all Optimize all residues beginning from the passed value (overrides other options);
         * for box optimization, optimizes all boxes beginning from the passed index.
         */
        @Option(shortName = 'x', longName = 'all', defaultValue = '-1',
                description = 'Optimize all residues beginning from the passed value (overrides other options); for box optimization, optimizes all boxes beginning from the passed index.')
        int all
        /**
         * -v or --verbose Prints beginning and default-conformation energies.
         */
        @Option(shortName = 'v', longName = 'verbose',
                description = 'Prints beginning and default-conformation energies.')
        boolean verbose
        /**
         * -o or --original Do not include starting coordinates as their own rotamer.
         */
        @Option(shortName = 'o', longName = 'original',
                description = 'Do not include starting coordinates as their own rotamer.')
        boolean original
        /**
         * -d or --decompose Print energy decomposition for original-coordinates rotamers.
         */
        @Option(shortName = 'd', longName = 'decompose',
                description = 'Print energy decomposition for original-coordinates rotamers.')
        boolean decompose
        /**
         * -eR or --energyRestart Load energy restart file from a previous run (requires that all parameters are the same).
         */
        @Option(shortName = 'eR', longName = 'energyRestart', defaultValue = 'none',
                description = 'Load energy restart file from a previous run (requires that all parameters are the same).')
        String energyRestart
        /**
         * -nB or --numBoxes Specify number of boxes along X, Y, and Z (default: '3,3,3').
         */
        @Option(shortName = 'nB', longName = 'numBoxes', defaultValue = '3,3,3',
                description = 'Specify number of boxes along X, Y, and Z (default: 3,3,3)')
        String numBoxes
        /**
         * -bB or --boxBorderSize Extent of overlap between optimization boxes in Angstroms (default: 3.0).
         */
        @Option(shortName = 'bB', longName = 'boxBorderSize', defaultValue = '3.0',
                description = 'Extent of overlap between optimization boxes in Angstroms.')
        double boxBorderSize
        /**
         * -bL or --approxBoxLength Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).
         *  Box sizes are rounded up to make a whole number of boxes along each axis (default of 0 disables this function).
         */
        @Option(shortName = 'bL', longName = 'approxBoxLength', defaultValue = '0.0',
                description = 'Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).')
        double approxBoxLength
        /**
         * -bC or --boxInclusionCriterion Criterion to use for adding residues to boxes.
         *      (1) uses C alpha only (N1/9 for nucleic acids)
         *      (2) uses any atom.
         *      (3) uses any rotamer
         */
        @Option(shortName = 'bC', longName = 'boxInclusionCriterion', defaultValue = '1',
                description = 'Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer')
        int boxInclusionCriterion
        /**
         * -fR or --forceResidues Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.
         */
        @Option(shortName = 'fR', longName = 'forceResidues', defaultValue = '-1,-1',
                description = 'Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.')
        String forceResidues
        /**
         * -lR or --listResidues Choose a list of individual residues to optimize (eg. A11,A24,B40).
         */
        @Option(shortName = 'lR', longName = 'listResidues', defaultValue = 'none',
                description = 'Choose a list of individual residues to optimize (eg. A11,A24,B40).')
        String listResidues
        /**
         * -sO or --sequence Choose a list of individual residues to sequence optimize (example: A2.A3.A5).
         */
        @Option(shortName = 'sO', longName = 'sequence', defaultValue = 'none',
                description = 'Choose a list of individual residues to sequence optimize (example: A2.A3.A5)')
        String sequence
        /**
         * -tO or --titrationOptimization Optimize the titration states for a list of residues (example: H2.H3.H5).
         */
        @Option(shortName = 'tO', longName = 'titrationOptimization', defaultValue = 'none',
                description = 'Optimize the titration states for a list of residues (example: H2.H3.H5).')
        String titrationOptimization
        /**
         * -nT or --nucleicCorrectionThreshold Nucleic acid rotamers adjusted by more than a threshold distance (A) are discarded (0 disables this function).
         */
        @Option(shortName = 'nT', longName = 'nucleicCorrectionThreshold', defaultValue = '0',
                description = 'Nucleic acid rotamers adjusted by more than a threshold distance (A) are discarded (0 disables this function).')
        double nucleicCorrectionThreshold
        /**
         * -mN or --minimumAcceptedNARotamers Minimum number of NA rotamers to be accepted if a threshold distance is enabled.
         */
        @Option(shortName = 'mN', longName = 'minimumAcceptedNARotamers', defaultValue = '10',
                description = 'Minimum number of NA rotamers to be accepted if a threshold distance is enabled.')
        int minimumAcceptedNARotamers
        /**
         * -mC or --monteCarlo Follow elimination criteria with 'n' Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.
         */
        @Option(shortName = 'mC', longName = 'monteCarlo', defaultValue = '-1',
                description = 'Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.')
        int monteCarlo
        /**
         * -z or --noRevert Do not revert an unfavorable change.
         */
        @Option(shortName = 'z', longName = 'revert',
                description = 'Undo an unfavorable change.')
        boolean noRevert
        /**
         * -tC or --twoBodyCutoff Cutoff distance for two-body interactions.
         */
        @Option(shortName = 'tC', longName = 'twoBodyCutoff', defaultValue = '-1.0',
                description = 'Cutoff distance for two body interactions.')
        double twoBodyCutoff
        /**
         * -thC or --threeBodyCutoff Cutoff distance for three-body interactions.
         */
        @Option(shortName = 'thC', longName = 'threeBodyCutoff', defaultValue = '9.0',
                description = 'Cutoff distance for three-body interactions.')
        double threeBodyCutoff

        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed
        List<String> filenames;
    }

    def run() {

        def cli = new CliBuilder(usage: ' ffxc ManyBody [options] <filename> [map file]', header: ' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }

        // Read in command line.
        List<String> arguments = options.filenames
        String filename = arguments.get(0);

        /**
         * Algorithm options.
         */
        int library = options.library
        int algorithm = options.algorithm
        boolean useGoldstein = !options.dee
        boolean threeBodyTerm = options.threeBody
        int pruningType = options.prune
        boolean revert = !options.noRevert
        boolean useOrigCoordsRotamer = !options.original
        boolean verboseEnergies = options.verbose
        boolean decomposeOriginal = options.decompose
        double twoBodyCutoff = options.twoBodyCutoff
        double threeBodyCutoff = options.threeBodyCutoff
        if (decomposeOriginal) {
            useOrigCoordsRotamer = true;
        }

        // By default, rotamer optimization should silence GK warnings, because occasionally we will have unreasonable configurations.
        if (System.getProperty("gk-suppressWarnings") == null) {
            System.setProperty("gk-suppressWarnings", "true");
        }

        /**
         * Chain, Residue and/or Box selections.
         */
        String chain = options.chain;
        int startResID = options.start
        int allStartResID = options.all
        int finalResID = options.last
        // Internal machinery indexed 0 to (n-1)
        int boxStart = options.start - 1
        int boxEnd = options.last - 1

        if (algorithm != 5) {
            // Not Box optimization.
            if (allStartResID < 1 && options.listResidues.equalsIgnoreCase('none')) {
                if (finalResID < startResID || startResID < 0 || finalResID < 0) {
                    logger.warning(" FFX shutting down: no residues specified for optimization.")
                    return;
                }
            }
        } else {
            // Box optimization.
            if (allStartResID > 0) {
                boxStart = allStartResID - 1 // Internal machinery indexed 0 to (n-1)
                if (boxStart < 0) {
                    logger.warning(" FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).")
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
         * Sliding window options.
         */
        int windowSize = options.window
        int increment = options.increment
        double distance = options.cutoff

        /**
         * Box optimization options.
         */
        int[] numXYZBoxes = new int[3];
        if (algorithm == 5) {
            String input = options.numBoxes
            Scanner boxNumInput = new Scanner(input);
            boxNumInput.useDelimiter(",");
            int inputLoopCounter = 0;
            numXYZBoxes[0] = 3; // Default
            while (inputLoopCounter < numXYZBoxes.length) {
                if (boxNumInput.hasNextInt()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextInt();
                    inputLoopCounter++;
                } else if (boxNumInput.hasNextDouble()) {
                    numXYZBoxes[inputLoopCounter] = boxNumInput.nextDouble();
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

        double boxBorderSize = options.boxBorderSize
        double approxBoxLength = options.approxBoxLength
        if (approxBoxLength < 0) {
            logger.info(" Negative box length value changed to -1 * input.")
            approxBoxLength *= -1;
        }
        int boxInclusionCriterion = options.boxInclusionCriterion

        /**
         * Minimize the final structure.
         */
        boolean min = false;
        double eps = options.minimize
        if (eps > 0.0) {
            min = true
        }

        /**
         * Monte Carlo search of remaining permutations.
         */
        boolean monteCarlo = false
        int nMCSteps = options.monteCarlo;
        if (nMCSteps > 1) {
            monteCarlo = true
        }

        /**
         * Nucleic acid options.
         */
        double nucleicCorrectionThreshold = options.nucleicCorrectionThreshold
        int minimumNumberAcceptedNARotamers = options.minimumAcceptedNARotamers

        /**
         * Energy restart.
         */
        boolean useEnergyRestart = false;
        File energyRestartFile = null;
        if (!options.energyRestart.equalsIgnoreCase('none')) {
            useEnergyRestart = true;
            energyRestartFile = new File(options.energyRestart)
        }

        /**
         * Force residues.
         */
        String forceResidues = options.forceResidues
        int forceResiduesStart = -1
        int forceResiduesEnd = -1;

        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();

        List<String> resList = new ArrayList<>();
        if (!options.listResidues.equalsIgnoreCase('none')) {
            def tok = (options.listResidues).tokenize('.');
            for (String t : tok) {
                logger.info(" Adding " + t);
                resList.add(t);
            }
        }

        List<String> sequenceOptimizationList = new ArrayList<>();
        if (!options.sequence.equalsIgnoreCase('none')) {
            def tok = options.sequence.tokenize('.');
            for (String t : tok) {
                logger.info(" Sequence optimizing " + t);
                sequenceOptimizationList.add(t);
            }
            if (System.getProperty("relative-solvation") == null) {
                System.setProperty("relative-solvation", "AUTO");
            }
        }

        List<String> titrationOptimizationList = new ArrayList<>();
        if (!options.titrationOptimization.equalsIgnoreCase('none')) {
            def tok = options.titrationOptimization.tokenize('.');
            for (String t : tok) {
                logger.info(" Protonation state optimizing " + t);
                titrationOptimizationList.add(t);
            }
        }

        /**
         * Evaluate forced residues for the sliding window algorithm
         */
        if (algorithm == 4 && !forceResidues.equalsIgnoreCase('-1,-1')) {
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
            if (!options.listResidues.equalsIgnoreCase('none')) {
                String info = "\n Evaluating rotamers for residues ";
                for (String i : resList) {
                    info += String.format("%s, ", i);
                }
                logger.info(info);
            } else if (allStartResID == -1) {
                logger.info("\n Evaluating rotamers for residues " + startResID + " to " + finalResID);
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

        open(filename);
        active.getPotentialEnergy().setPrintOnFailure(false, false);

        def mapFiles = [];
        int nDiffractionData = 0;
        if (arguments.size() > 1) {
            String dataFileName = arguments.get(1);
            if (FilenameUtils.isExtension(dataFileName, "map")) {
                RealSpaceFile realspacefile = new RealSpaceFile(dataFileName, 1.0);
                mapFiles.add(realspacefile);
            } else {
                DiffractionFile diffractionFile = new DiffractionFile(dataFileName, 1.0, false);
                DiffractionData diffractionData = new DiffractionData(active, active.getProperties(), CrystalReciprocalSpace.SolventModel.POLYNOMIAL, diffractionFile);
                diffractionData.scaleBulkFit();
                diffractionData.printStats();
                String mapFileName = String.format("%s_ffx_%d", FilenameUtils.removeExtension(dataFileName), ++nDiffractionData);
                diffractionData.writeMaps(mapFileName);
                mapFiles.add(new RealSpaceFile(mapFileName + "_2fofc.map", 1.0));
            }
        } else {
            RealSpaceFile realspacefile = new RealSpaceFile(active);
            mapFiles.add(realspacefile);
        }


        RealSpaceData realSpaceData = new RealSpaceData(active,
                active.getProperties(), active.getParallelTeam(),
                mapFiles.toArray(new RealSpaceFile[mapFiles.size()]));
        RefinementEnergy refinementEnergy = new RefinementEnergy(realSpaceData, RefinementMinimize.RefinementMode.COORDINATES, null);

        double[] x = new double[refinementEnergy.getNumberOfVariables()];
        x = refinementEnergy.getCoordinates(x);
        refinementEnergy.energy(x, true);

        RotamerOptimization rotamerOptimization = new RotamerOptimization(active, refinementEnergy, sh);
        rotamerOptimization.setTwoBodyCutoff(twoBodyCutoff);
        rotamerOptimization.setThreeBodyCutoff(threeBodyCutoff);
        rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
        rotamerOptimization.setUseGoldstein(useGoldstein);
        rotamerOptimization.setRevert(revert);
        rotamerOptimization.setPruning(pruningType);
        rotamerOptimization.setWindowSize(windowSize);
        rotamerOptimization.setIncrement(increment);
        rotamerOptimization.setDistanceCutoff(distance);
        rotamerOptimization.setNucleicCorrectionThreshold(nucleicCorrectionThreshold);
        rotamerOptimization.setMinimumNumberAcceptedNARotamers(minimumNumberAcceptedNARotamers);
        rotamerOptimization.setVerboseEnergies(verboseEnergies);
        rotamerOptimization.setBoxBorderSize(boxBorderSize);
        rotamerOptimization.setApproxBoxLength(approxBoxLength);
        rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
        rotamerOptimization.setBoxInclusionCriterion(boxInclusionCriterion);
        rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
        rotamerOptimization.setMonteCarlo(monteCarlo, nMCSteps);

        if (useEnergyRestart) {
            rotamerOptimization.setEnergyRestartFile(energyRestartFile);
        }

        if (library == 1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
        } else {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
        }

        if (useOrigCoordsRotamer) {
            rLib.setUseOrigCoordsRotamer(true);
        }

        int counter = 1;
        if (algorithm != 5) {
            if (allStartResID > 0) {
                ArrayList<Residue> residueList = new ArrayList<Residue>();
                Polymer[] polymers = active.getChains();
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
                        } else if (!options.forceResidues.equalsIgnoreCase('none')) {
                            if (counter >= allStartResID && counter >= forceResiduesStart
                                    && counter <= forceResiduesEnd) {
                                residueList.add(residue);
                            }
                        }
                        counter++;
                    }
                }
                rotamerOptimization.setResidues(residueList);
            } else if (!options.listResidues.equalsIgnoreCase('none')) {
                ArrayList<Residue> residueList = new ArrayList<>();
                Polymer[] polymers = active.getChains();
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
            } else if (!chain.equalsIgnoreCase('-1')) {
                rotamerOptimization.setResidues(chain, startResID, finalResID);
            } else {
                rotamerOptimization.setResidues(startResID, finalResID);
            }
        } else {
            boolean ignoreNA = false;
            String ignoreNAProp = System.getProperty("ignoreNA");
            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
                ignoreNA = true;
            }
            ArrayList<Residue> residueList = new ArrayList<Residue>();
            Polymer[] polymers = active.getChains();
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

        ArrayList<Residue> residueList = rotamerOptimization.getResidues();

        if (!options.sequence.equalsIgnoreCase('none')) {
            for (String s : sequenceOptimizationList) {
                Character chainID = s.charAt(0);
                int num = Integer.parseInt(s.substring(1));
                for (int i = 0; i < residueList.size(); i++) {
                    Residue res = residueList.get(i);
                    if (res.getChainID() == chainID && res.getResidueNumber() == num) {
                        MultiResidue multiRes = new MultiResidue(res, active.getForceField(), active.getPotentialEnergy());
                        for (Polymer polymer : active.getChains()) {
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
                        //multiRes.requestSetActiveResidue(ResidueEnumerations.AminoAcid3.valueOf(res.getName()));
                        multiRes.setActiveResidue(res);
                        active.getPotentialEnergy().reInit();
                        residueList.remove(i);
                        residueList.add(i, multiRes);
                    }
                }
            }
        }

        if (!options.titrationOptimization.equalsIgnoreCase('none')) {
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

        boolean master = true;
        if (Comm.world().size() > 1) {
            int rank = Comm.world().rank();
            if (rank != 0) {
                master = false;
            }
        }

        x = refinementEnergy.getCoordinates(x);
        refinementEnergy.energy(x, true);

        RotamerLibrary.measureRotamers(residueList, false);

        if (decomposeOriginal) {
            rLib.setUseOrigCoordsRotamer(true);
            rotamerOptimization.setDecomposeOriginal(true);
        }

        if (algorithm == 1) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.INDEPENDENT);
        } else if (algorithm == 2) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
        } else if (algorithm == 3) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BRUTE_FORCE);
        } else if (algorithm == 4) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.WINDOW);
        } else if (algorithm == 5) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BOX);
        }

        if (master) {
            if (min) {
                minimize(eps);
            }
            logger.info(" Final Minimum Energy");
            energy();
            String ext = FilenameUtils.getExtension(filename);
            filename = FilenameUtils.removeExtension(filename);
            if (ext.toUpperCase().contains("XYZ")) {
                saveAsXYZ(new File(filename + ".xyz"));
            } else {
                saveAsPDB(new File(filename + ".pdb"));
            }
        }
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

