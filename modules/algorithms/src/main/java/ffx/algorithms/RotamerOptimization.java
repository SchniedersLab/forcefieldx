/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.algorithms;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.BiFunction;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import static java.util.Arrays.fill;

import org.apache.commons.io.FilenameUtils;

import edu.rit.mp.BooleanBuf;
import edu.rit.mp.Buf;
import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import edu.rit.pj.CommStatus;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.WorkerIteration;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;

import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.MCMove;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.NACorrectionException;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.IndexIndexPair;
import ffx.utilities.ObjectPair;

import static ffx.potential.bonded.Residue.ResidueType.AA;
import static ffx.potential.bonded.Residue.ResidueType.NA;
import static ffx.potential.bonded.RotamerLibrary.applyRotamer;

/**
 * Optimize protein side-chain conformations and nucleic acid backbone
 * conformations using rotamers.
 *
 * @author Michael J. Schnieders
 *
 * @author Jacob M. Litman
 *
 * @author Stephen D. LuCore
 *
 * @since 1.0
 */
public class RotamerOptimization implements Terminatable {

    private static final Logger logger = Logger.getLogger(RotamerOptimization.class.getName());
    /**
     * Flag to control the verbosity of printing.
     */
    private boolean print = false;
    /**
     * AlgorithmListener who should receive updates as the optimization runs.
     */
    protected AlgorithmListener algorithmListener;
    /**
     * Flag to indicate a request to terminate the optimization.
     */
    private boolean terminate = false;
    /**
     * Flag to indicate if the algorithm is running (done == false) or completed
     * (done == true).
     */
    private boolean done = true;
    /**
     * MolecularAssembly to perform rotamer optimization on.
     */
    protected final MolecularAssembly molecularAssembly;
    /**
     * List of Assemblies associated with a multi-topology Potential.
     */
    protected final List<MolecularAssembly> allAssemblies;
    /**
     * The Potential to evaluate during rotamer optimization.
     */
    protected Potential potential;
    /**
     * Number of permutations whose energy is explicitly evaluated.
     */
    private int evaluatedPermutations = 0;
    /**
     * An array of polymers from the MolecularAssembly.
     */
    private Polymer[] polymers;
    /**
     * The chain containing the residues to optimize.
     */
    private String chain;
    /**
     * List of residues to optimize; they may not be contiguous or all members
     * of the same chain.
     */
    private ArrayList<Residue> residueList;
    /**
     * Size of the sliding window.
     */
    private int windowSize = 7;
    /**
     * The distance the sliding window moves.
     */
    private int increment = 3;
    /**
     * In sliding window, whether to revert an unfavorable change.
     */
    private boolean revert = false;
    /**
     * The distance that the distance matrix checks for.
     */
    private double distance = 2.0;
    /**
     * The algorithm to use for rotamer optimization.
     */
    private Algorithm algorithm = null;
    /**
     * The sliding window direction.
     */
    private Direction direction = Direction.FORWARD;
    /**
     * Flag to indicate use of the Goldstein criteria instead of the less
     * stringent Dead-End Elimination criteria.
     */
    private boolean useGoldstein = true;
    /**
     * The number of most-favorable structures to include as output.
     */
    private int ensembleNumber = 1;
    /**
     * The energy buffer applied to each elimination criteria to affect an
     * ensemble.
     */
    private double ensembleBuffer = 0.0;
    /**
     * The step value of the energy buffer for use with ensemble search.
     */
    private double ensembleBufferStep = 0.5;
    /**
     * The energy boundary for structures to be included in the final ensemble.
     */
    private double ensembleEnergy = 0.0;
    /**
     * File to contain ensemble of structures.
     */
    private File ensembleFile;
    /**
     * PDBFilter to write out ensemble snapshots.
     */
    private PDBFilter ensembleFilter;
    /**
     * The potential energy of the system with all side-chains to be optimized
     * turned off.
     */
    protected double backboneEnergy;
    /**
     * Self-energy of each residue for each rotamer. [residue][rotamer]
     */
    protected double selfEnergy[][];
    /**
     * Pair-energies for each pair of residue and pair of rotamers.
     * [residue1][rotamer1][residue2][rotamer2]
     */
    protected double twoBodyEnergy[][][][];
    /**
     * The minimum distance between atoms of a residue pair, taking into account
     * interactions with symmetry mates.
     *
     * [residue1][rotamer1][residue2][rotamer2]
     */
    private double distanceMatrix[][][][];
    /**
     * Flag to load the distance matrix as needed; if false, matrix is prefilled
     * at the beginning of rotamer optimization.
     */
    private boolean lazyMatrix = false;
    /**
     * Trimer-energies for each trimer of rotamers.
     * [residue1][rotamer1][residue2][rotamer2][residue3][rotamer3]
     */
    protected double threeBodyEnergy[][][][][][];
    /**
     * Flag to prune clashes.
     */
    protected boolean pruneClashes = true;
    /**
     * Flag to prune pair clashes.
     */
    protected boolean prunePairClashes = true;
    /**
     * Flag to prune individual pairs (and not just entire rotamers) on pair
     * clashes. Presently set constitutively false.
     */
    private boolean pruneIndivPairs = false;
    /**
     * Clash energy threshold (kcal/mole).
     */
    private double clashThreshold = 20.0;
    /**
     * Clash energy threshold (kcal/mol) for MultiResidues, which can have much
     * more variation in self and pair energies.
     */
    private double multiResClashThreshold = 80.0;
    /**
     * Clash energy threshold (kcal/mole).
     */
    private double pairClashThreshold = 50.0;
    /**
     * Pair clash energy threshold (kcal/mol) for MultiResidues.
     */
    private double multiResPairClashAddn = 80.0;
    /**
     * Flag to control use of 3-body terms.
     */
    protected boolean threeBodyTerm = true;
    /**
     * Flag to set 3-body energies to zero outside of a cutoff.
     */
    private boolean threeBodyCutoff = true;
    /**
     * Triple cutoff distance.
     */
    private double threeBodyCutoffDist = 9.0;
    /**
     * Quad cutoff flag and distance.
     */
    private boolean quadCutoff = true;
    private double quadCutoffDist = 5.0;
    /**
     * Eliminated rotamers. [residue][rotamer]
     */
    private boolean eliminatedSingles[][];
    /**
     * Eliminated rotamer pairs. [residue1][rotamer1][residue2][rotamer2]
     */
    private boolean eliminatedPairs[][][][];
    /**
     * Eliminated rotamer pairs.
     * [residue1][rotamer1][residue2][rotamer2][residue3][rotamer3]
     */
    private boolean eliminatedTriples[][][][][][];
    /**
     * An array of atomic coordinates (length 3 * the number of atoms).
     */
    private double x[] = null;
    /**
     * A flag to indicate use of the full N-Body AMOEBA potential energy during
     * the rotamer optimization.
     */
    private boolean useFullAMOEBAEnergy = false;
    /**
     * Threshold to eliminate nucleic acid Rotamers based on excessive
     * correction distances; 0 indicates the threshold is not being implemented.
     */
    private double nucleicCorrectionThreshold = 0;
    /**
     * Minimum number of nucleic acid Rotamers to use for rotamer optimization
     * should some be eliminated by the nucleic correction threshold.
     */
    private int minNumberAcceptedNARotamers = 10;
    /**
     * Factor by which to multiply the pruning constraints for nucleic acids.
     * pairHalfPruningFactor is the arithmetic mean of 1.0 and the pruning
     * factor, and is applied for AA-NA pairs.
     */
    private double pruningFactor = 1.0;
    private double pairHalfPruningFactor = ((1.0 + pruningFactor) / 2);
    /**
     * Factor by which to multiply the singleton pruning constraints for nucleic
     * acids.
     *
     * Very important, to ensure that all possible combinations of delta(i) and
     * delta(i-1) are still represented when it comes time to calculate pair
     * energies. If this pruning factor doesn't cut it, however, it probably
     * wasn't a biologically relevant rotamer anyways.
     *
     * Not presently implemented beyond getting the value in from Groovy.
     */
    private double singletonNAPruningFactor = 1.5;
    /**
     * Flag to calculate and print additional energies (mostly for debugging).
     */
    private boolean verboseEnergies = true;
    private ArrayList<Residue> allResiduesList = null;
    private Residue allResiduesArray[] = null;
    private int numResidues = 0;
    private boolean parallelEnergies = true;
    private final Comm world;
    private final int numProc;
    private final int rank;
    private final boolean master;
    private ReceiveThread receiveThread;
    private EnergyWriterThread energyWriterThread;
    private boolean selfsDone = false, pairsDone = false, trimersDone = false, quadsDone = false;
    private boolean readyForSingles = false, readyForPairs = false, readyForTrimers = false, readyForQuads = false;
    private boolean writeEnergyRestart = true;
    private boolean loadEnergyRestart = false;
    private File energyRestartFile;
    private final HashMap<Integer, Integer[]> jobMapSingles = new HashMap<>();
    private final HashMap<Integer, Integer[]> jobMapPairs = new HashMap<>();
    private final HashMap<Integer, Integer[]> jobMapTrimers = new HashMap<>();
    private final HashMap<Integer, Integer[]> jobMapQuads = new HashMap<>();
    private List<String> energiesToWrite;
    private boolean verbose = false;
    // In X, Y, Z.
    private int[] numXYZBoxes = {3, 3, 3};
    private double boxBorderSize = 0;
    private double approxBoxLength = 0;
    private int boxInclusionCriterion = 1;
    private int boxStart = 0;
    private int boxEnd = -1;
    private boolean manualSuperbox = false;
    private double[] boxDimensions;
    private double superboxBuffer = 8.0;
    private int startForcedResidues = -1;
    private int endForcedResidues = -1;
    private boolean useForcedResidues = false;
    private double superpositionThreshold = 0.1;
    private boolean usingBoxOptimization = false;
    private int boxLoadIndex = -1;
    private int[] boxLoadCellIndices;
    /**
     * Wait time, in milliseconds, for Thread.sleep() loops in rotamerEnergies()
     * and ReceiveThread
     */
    private final int POLLING_FREQUENCY = 500;
    /**
     * ParallelJava construct used for executing the various EnergyRegions.
     * Declared as a field so as to be reused during box optimization.
     */
    private final WorkerTeam energyWorkerTeam;
    private VideoWriter videoWriter;
    private boolean writeVideo = false;
    private boolean skipEnergies = false;
    private boolean computeQuads = false;
    private boolean decomposeOriginal = false;
    private boolean addOrigRot = false; // Using original-rotamers for ALA, GLY, etc.
    private int quadMaxout = Integer.MAX_VALUE;

    /**
     * Monte Carlo parameters.
     */
    private int nMCsteps = 1000000;
    private boolean monteCarlo = false;
    private double mcTemp = 298.15;
    /**
     * Check to see if proposed move has an eliminated pair or higher-order
     * term; breaks detailed balance.
     */
    private boolean mcUseAll = false;
    // Skips brute force enumeration in favor of pure Monte Carlo. Recommended only for testing.
    private boolean mcNoEnum = false;

    /**
     * Sets whether files should be printed; true for standalone applications,
     * false for some applications which use rotamer optimization as part of a
     * larger process.
     */
    private boolean printFiles = true;
    /**
     * Stores states of each ensemble if printFiles is false.
     */
    private List<ObjectPair<ResidueState[], Double>> ensembleStates;
    protected RotamerLibrary library = RotamerLibrary.getDefaultLibrary();

    /**
     * Represents the method called to obtain the directory corresponding to
     * the current energy; will be a simple return null for potential energy
     * evaluations. While current energy calls will fill the rotamer list with
     * the current rotamers of the residue, other methods may skip applying
     * the rotamer directly.
     */
    private BiFunction<List<Residue>, List<Rotamer>, File> dirSupplier;
    /**
     * Represents the method called to obtain energy for the current rotamer or
     * state; defaults to the existing potential energy code. May discard the
     * input file.
     */
    private ToDoubleFunction<File> eFunction;

    /**
     * RotamerOptimization constructor.
     *
     * @param molecularAssembly The MolecularAssembly to search rotamers for.
     * @param potential
     * @param algorithmListener AlgorithmListener to update the GUI.
     */
    public RotamerOptimization(MolecularAssembly molecularAssembly, Potential potential,
            AlgorithmListener algorithmListener) {

        this.molecularAssembly = molecularAssembly;
        this.potential = potential;
        this.algorithmListener = algorithmListener;
        eFunction = this::currentPE;
        dirSupplier = (List<Residue> resList, List<Rotamer> rotList) -> null;
        world = Comm.world();
        numProc = world.size();
        rank = world.rank();
        master = rank == 0;
        energyWorkerTeam = new WorkerTeam(world);
        if (System.getProperty("verbose") != null) {
            if (System.getProperty("verbose").equalsIgnoreCase("true")) {
                verbose = true;
            }
        }

        // Process relevant system keys.
        String skipEnergies = System.getProperty("ro-skipEnergies");
        String undo = System.getProperty("ro-undo");
        String direction = System.getProperty("ro-direction");
        String increment = System.getProperty("ro-increment");
        String goldstein = System.getProperty("ro-goldstein");
        String parallelEnergies = System.getProperty("ro-parallelEnergies");
        String superpositionThreshold = System.getProperty("ro-superpositionThreshold");
        String ensembleNumber = System.getProperty("ro-ensembleNumber");
        String ensembleEnergy = System.getProperty("ro-ensembleEnergy");
        String ensembleBuffer = System.getProperty("ro-ensembleBuffer");
        String threeBodyCutoffDist = System.getProperty("ro-threeBodyCutoffDist");
        String pruningFactor = System.getProperty("ro-pruningFactor");
        String nucleicSinglesPruningFactor = System.getProperty("ro-nucleicSinglesPruningFactor");
        String nucleicCorrectionThreshold = System.getProperty("ro-nucleicCorrectionThreshold");
        String minimumNumberAcceptedNARotamers = System.getProperty("ro-minimumNumberAcceptedNARotamers");
        String singletonClashThreshold = System.getProperty("ro-singletonClashThreshold");
        String multiResClashThreshold = System.getProperty("ro-multiResClashThreshold");
        String pairClashThreshold = System.getProperty("ro-pairClashThreshold");
        String multiResPairClashAddition = System.getProperty("ro-multiResPairClashAddition");
        String boxDimensions = System.getProperty("ro-boxDimensions");
        String computeQuads = System.getProperty("ro-computeQuads");
        String quadCutoffDist = System.getProperty("ro-quadCutoffDist");
        String quadMaxout = System.getProperty("ro-quadMaxout");
        String lazyMatrix = System.getProperty("ro-lazyMatrix");
        String mcTemp = System.getProperty("ro-mcTemp");
        String mcUseAll = System.getProperty("ro-mcUseAll");
        String mcNoEnum = System.getProperty("ro-debug-mcNoEnum");
        String addOrigRotStr = System.getProperty("ro-addOrigRot");
        String origAtEndStr = System.getProperty("ro-origAtEnd");
        if (computeQuads != null) {
            boolean value = Boolean.parseBoolean(computeQuads);
            this.computeQuads = value;
            logger.info(String.format(" (KEY) computeQuads: %b", this.computeQuads));
        }
        if (quadCutoffDist != null) {
            double value = Double.parseDouble(quadCutoffDist);
            this.quadCutoffDist = value;
            if (this.quadCutoffDist < 0) {
                quadCutoff = false;
            }
            logger.info(String.format(" (KEY) quadCutoffDist: %.2f", this.quadCutoffDist));
        }
        if (quadMaxout != null) {
            int value = Integer.parseInt(quadMaxout);
            this.quadMaxout = value;
            logger.info(String.format(" (KEY) quadMaxout: %d", this.quadMaxout));
        }
        if (skipEnergies != null) {
            boolean value = Boolean.parseBoolean(skipEnergies);
            this.skipEnergies = value;
            logger.info(String.format(" (KEY) skipEnergies: %b", this.skipEnergies));
        }
        if (undo != null) {
            boolean value = Boolean.parseBoolean(undo);
            this.revert = value;
            logger.info(String.format(" (KEY) undo: %b", this.revert));
        }
        if (direction != null) {
            Direction value = Direction.valueOf(direction);
            this.direction = value;
            logger.info(String.format(" (KEY) direction: %s", this.direction.toString()));
        }
        if (increment != null) {
            int value = Integer.parseInt(increment);
            this.increment = value;
            logger.info(String.format(" (KEY) increment: %d", this.increment));
        }
        if (goldstein != null) {
            boolean value = Boolean.parseBoolean(goldstein);
            this.useGoldstein = value;
            logger.info(String.format(" (KEY) goldstein: %b", this.useGoldstein));
        }
        if (parallelEnergies != null) {
            boolean value = Boolean.parseBoolean(parallelEnergies);
            this.parallelEnergies = value;
            logger.info(String.format(" (KEY) parallelEnergies: %b", this.parallelEnergies));
        }
        if (superpositionThreshold != null) {
            Double value = Double.parseDouble(superpositionThreshold);
            this.superpositionThreshold = value;
            logger.info(String.format(" (KEY) superpositionThreshold: %.2f", this.superpositionThreshold));
        }
        if (ensembleNumber != null) {
            int value = Integer.parseInt(ensembleNumber);
            this.ensembleNumber = value;
            this.ensembleBuffer = 5.0;
            this.ensembleBufferStep = 0.1 * this.ensembleBuffer;
            logger.info(String.format(" (KEY) ensembleNumber: %d", this.ensembleNumber));
        }
        if (ensembleBuffer != null) {
            double value = Double.parseDouble(ensembleBuffer);
            this.ensembleBuffer = value;
            this.ensembleBufferStep = 0.1 * this.ensembleBuffer;
            logger.info(String.format(" (KEY) ensembleBuffer: %.2f", this.ensembleBuffer));
        }
        if (ensembleEnergy != null) {
            double value = Double.parseDouble(ensembleEnergy);
            this.ensembleEnergy = value;
            logger.info(String.format(" (KEY) ensembleEnergy: %.2f", this.ensembleEnergy));
        }
        if (threeBodyCutoffDist != null) {
            double value = Double.parseDouble(threeBodyCutoffDist);
            this.threeBodyCutoffDist = value;
            if (this.threeBodyCutoffDist < 0) {
                threeBodyCutoff = false;
            }
            logger.info(String.format(" (KEY) threeBodyCutoffDist: %.2f", this.threeBodyCutoffDist));
        }
        if (pruningFactor != null) {
            double value = Double.parseDouble(pruningFactor);
            this.pruningFactor = (value >= 0 ? value : 1.0);
            this.pairHalfPruningFactor = (1.0 + value) / 2;
            logger.info(String.format(" (KEY) pruningFactor: %.2f", this.pruningFactor));
        }
        if (nucleicSinglesPruningFactor != null) {
            double value = Double.parseDouble(nucleicSinglesPruningFactor);
            this.singletonNAPruningFactor = (value >= 0 ? value : 1.5);
            logger.info(String.format(" (KEY) nucleicSinglesPruningFactor: %.2f", this.singletonNAPruningFactor));
        }
        if (nucleicCorrectionThreshold != null) {
            double value = Double.parseDouble(nucleicCorrectionThreshold);
            this.nucleicCorrectionThreshold = (value >= 0 ? value : 0);
            logger.info(String.format(" (KEY) nucleicCorrectionThreshold: %.2f", this.nucleicCorrectionThreshold));
        }
        if (minimumNumberAcceptedNARotamers != null) {
            int value = Integer.parseInt(minimumNumberAcceptedNARotamers);
            this.minNumberAcceptedNARotamers = (value > 0 ? value : 10);
            logger.info(String.format(" (KEY) minimumNumberAcceptedNARotamers: %d", this.minNumberAcceptedNARotamers));
        }
        if (singletonClashThreshold != null) {
            double value = Double.parseDouble(singletonClashThreshold);
            this.clashThreshold = value;
            logger.info(String.format(" (KEY) singletonClashThreshold: %.2f", this.clashThreshold));
        }
        if (multiResClashThreshold != null) {
            double value = Double.parseDouble(multiResClashThreshold);
            this.multiResClashThreshold = value;
            logger.info(String.format(" (KEY) multiResClashThreshold: %.2f", this.multiResClashThreshold));
        }
        if (pairClashThreshold != null) {
            double value = Double.parseDouble(pairClashThreshold);
            this.pairClashThreshold = value;
            logger.info(String.format(" (KEY) pairClashThreshold: %.2f", this.pairClashThreshold));
        }
        if (multiResPairClashAddition != null) {
            double value = Double.parseDouble(multiResPairClashAddition);
            this.multiResPairClashAddn = value;
            logger.info(String.format(" (KEY) multiResPairClashAddition: %.2f", this.multiResPairClashAddn));
        }
        if (boxDimensions != null) {
            // String should be in format (buffer,xmin,xmax,ymin,ymax,zmin,zmax)
            try {
                String[] bdTokens = boxDimensions.split(",+");
                this.boxDimensions = new double[6];
                if (bdTokens.length != 7) {
                    logger.warning(" Improper number of arguments to boxDimensions; default settings used.");
                } else {
                    for (int i = 1; i < 7; i += 2) {
                        this.boxDimensions[i - 1] = Double.parseDouble(bdTokens[i]);
                        this.boxDimensions[i] = Double.parseDouble(bdTokens[i + 1]);
                        if (this.boxDimensions[i] < this.boxDimensions[i - 1]) {
                            logger.info(String.format(" Improper dimension min %8.5f > max %8.5f; max/min reversed.", this.boxDimensions[i - 1], this.boxDimensions[i]));
                            double temp = this.boxDimensions[i];
                            this.boxDimensions[i] = this.boxDimensions[i - 1];
                            this.boxDimensions[i - 1] = temp;
                        }
                    }
                    superboxBuffer = Double.parseDouble(bdTokens[0]);
                    manualSuperbox = true;
                }
            } catch (Exception ex) {
                logger.warning(String.format(" Error in parsing box dimensions: input discarded and defaults used: %s.", ex.toString()));
                manualSuperbox = false;
            }
        }
        if (mcTemp != null) {
            double value = Double.parseDouble(mcTemp);
            this.mcTemp = value;
            logIfMaster(String.format(" (KEY) mcTemp: %10.6f", this.mcTemp));
        }
        if (mcUseAll != null) {
            boolean value = Boolean.parseBoolean(mcUseAll);
            this.mcUseAll = value;
            logIfMaster(String.format(" (KEY) mcUseAll: %b", this.mcUseAll));
        }
        if (mcNoEnum != null) {
            boolean value = Boolean.parseBoolean(mcNoEnum);
            this.mcNoEnum = value;
            logIfMaster(String.format(" (KEY) debug-mcNoEnum: %b", this.mcNoEnum));
        }
        if (lazyMatrix != null) {
            boolean value = Boolean.parseBoolean(lazyMatrix);
            this.lazyMatrix = value;
            logger.info(String.format(" (KEY) lazyMatrix: %b", lazyMatrix));
        }
        if (addOrigRotStr != null) {
            boolean value = Boolean.parseBoolean(addOrigRotStr);
            this.addOrigRot = value;
            logger.info(String.format(" (KEY) addOrigRot: %b", addOrigRot));
        }
        if (origAtEndStr != null) {
            boolean value = Boolean.parseBoolean(origAtEndStr);
            // Property works in the contest of Residue class.
            logger.info(String.format(" (KEY) origAtEnd: %b", value));
        }
        allAssemblies = new ArrayList<>();
        allAssemblies.add(molecularAssembly);
    }

    public RotamerOptimization(MolecularAssembly molecularAssembly, Potential potential,
            AlgorithmListener algorithmListener, int startResID, int finalResID, Algorithm algorithm) {
        this(molecularAssembly, potential, algorithmListener);
        this.algorithm = algorithm;
    }

    public void addAdditionalAssembly(MolecularAssembly assembly) {
        this.allAssemblies.add(assembly);
    }

    /**
     * A brute-force global optimization over side-chain rotamers using a
     * recursive algorithm.
     *
     * @param molecularAssembly
     * @param residues
     * @param i
     * @param lowEnergy
     * @param optimum
     *
     * @return the current energy.
     */
    public double rotamerOptimization(MolecularAssembly molecularAssembly, Residue residues[], int i,
            double lowEnergy, int optimum[]) {

        // This is the initialization condition.
        if (i == 0) {
            evaluatedPermutations = 0;
        }

        int nResidues = residues.length;
        Residue current = residues[i];
        Rotamer[] rotamers = current.getRotamers(library);
        int lenri = rotamers.length;
        double currentEnergy = Double.MAX_VALUE;
        List<Residue> resList = Arrays.asList(residues);
        if (i < nResidues - 1) {
            /**
             * As long as there are more residues, continue the recursion for
             * each rotamer of the current residue.
             */
            int minRot = -1;
            for (int ri = 0; ri < lenri; ri++) {
                applyRotamer(current, rotamers[ri]);
                double rotEnergy = rotamerOptimization(molecularAssembly, residues, i + 1, lowEnergy, optimum);
                if (rotEnergy < currentEnergy) {
                    currentEnergy = rotEnergy;
                }
                if (rotEnergy < lowEnergy) {
                    minRot = ri;
                    lowEnergy = rotEnergy;
                }
            }
            if (minRot > -1) {
                optimum[i] = minRot;
            }
        } else {
            /**
             * At the end of the recursion, compute the potential energy for
             * each rotamer of the final residue. If a lower potential energy is
             * discovered, the rotamers of each residue will be collected as the
             * recursion returns up the chain.
             */
            for (int ri = 0; ri < lenri; ri++) {
                applyRotamer(current, rotamers[ri]);

                double rotEnergy = currentEnergy(resList);
                logger.info(String.format(" %d Energy: %16.8f", ++evaluatedPermutations, rotEnergy));
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }
                if (rotEnergy < currentEnergy) {
                    currentEnergy = rotEnergy;
                }
                if (rotEnergy < lowEnergy) {
                    lowEnergy = rotEnergy;
                    optimum[nResidues - 1] = ri;
                }
            }
        }
        return currentEnergy;
    }

    /**
     * Recursive brute-force method which uses single, pair, and potentially
     * trimer energies to calculate an optimum set of rotamers.
     *
     * @param molecularAssembly
     * @param residues Optimization window
     * @param i Current residue in the recursion.
     * @param lowEnergy Minimum energy yet found by the recursion.
     * @param optimum Optimum rotamer set yet found by the recursion.
     * @param currentRotamers Rotamer permutation under investigation.
     * @return Minimum energy found under this node in the recursion.
     */
    private double decomposedRotamerOptimization(MolecularAssembly molecularAssembly, Residue residues[], int i,
            double lowEnergy, int optimum[], int[] currentRotamers) {

        // This is the initialization condition.
        if (i == 0) {
            evaluatedPermutations = 0;
        }

        int nResidues = residues.length;
        Residue current = residues[i];
        Rotamer[] rotamers = current.getRotamers(library);
        int lenri = rotamers.length;
        double currentEnergy = Double.MAX_VALUE;
        if (i < nResidues - 1) {
            /**
             * As long as there are more residues, continue the recursion for
             * each rotamer of the current residue.
             */
            for (int ri = 0; ri < lenri; ri++) {
                currentRotamers[i] = ri;
                double rotEnergy = decomposedRotamerOptimization(molecularAssembly, residues, i + 1, lowEnergy, optimum, currentRotamers);
                if (rotEnergy < lowEnergy) {
                    lowEnergy = rotEnergy;
                }
                if (rotEnergy < currentEnergy) {
                    currentEnergy = rotEnergy;
                }
            }
        } else {
            /**
             * At the end of the recursion, compute the potential energy for
             * each rotamer of the final residue and update optimum[].
             */
            for (int ri = 0; ri < lenri; ri++) {
                currentRotamers[i] = ri;
                double rotEnergy = computeEnergy(residues, currentRotamers, false);
                ++evaluatedPermutations;
                if (rotEnergy < currentEnergy) {
                    currentEnergy = rotEnergy;
                }
                if (rotEnergy < lowEnergy) {
                    /* Because we print the rotamer set immediately on finding a
                     * more optimal structure, we have to reset the entire length
                     * of optimum instead of lazily doing it on the way down.
                     */
                    System.arraycopy(currentRotamers, 0, optimum, 0, optimum.length);
                    if (evaluatedPermutations > 1) {
                        logger.info(String.format(" Minimum energy update: %f < %f, permutation %d",
                                rotEnergy, lowEnergy, evaluatedPermutations));
                        String permutation = " Rotamer permutation: " + optimum[0];
                        for (int j = 1; j < nResidues; j++) {
                            permutation = permutation.concat(", " + optimum[j]);
                        }
                        logger.info(permutation);
                    } else {
                        logger.info(String.format(" First minimum energy (permutation 1): %f", rotEnergy));
                    }
                    lowEnergy = rotEnergy;
                }
            }
        }
        return currentEnergy;
    }

    /**
     * Runs Monte Carlo side chain optimization using the rotamer energy matrix
     * and potentially some information from dead-end or Goldstein elimination.
     * The useAllElims variable should be set false if detailed balance is to be
     * maintained. At present, no support for ensembles.
     *
     * @param residues Optimization window
     * @param optimum Array to store optimum rotamers
     * @param initialRots Array with starting rotamers
     * @param maxIters Number of MC steps to run
     * @param randomizeRots Scramble initialRots
     * @param useAllElims Use pair/triple elimination information
     * @return Lowest energy found
     */
    private double rotamerOptimizationMC(Residue[] residues, int[] optimum,
            int[] initialRots, int maxIters, boolean randomizeRots,
            boolean useAllElims) {

        long initTime = -System.nanoTime();
        if (randomizeRots) {
            randomizeRotamers(initialRots, residues, true);
        }

        int nRes = residues.length;
        System.arraycopy(initialRots, 0, optimum, 0, nRes);
        assert optimum.length == nRes;
        assert initialRots.length == nRes;

        RotamerMatrixMC rmc = new RotamerMatrixMC(initialRots, residues);
        rmc.setTemperature(mcTemp);
        RotamerMatrixMove rmove = new RotamerMatrixMove(useAllElims, initialRots, residues);
        List<MCMove> rmList = new ArrayList<>(1);
        rmList.add(rmove);

        double initialEnergy = computeEnergy(residues, initialRots, false);
        double optimumEnergy = initialEnergy;
        double currentEnergy = initialEnergy;

        int nAccept = 0;
        //int nReject = 0;

        /**
         * I have the vague idea of parallelizing down to individual threads, by
         * simply calling this method from each thread, with randomizeRots true.
         * That would require replacing logIfMaster with logger.info.
         */
        logIfMaster(String.format(" Beginning %d iterations of Monte Carlo search "
                + "starting from energy %10.6f", maxIters, initialEnergy));

        for (int i = 0; i < maxIters; i++) {
            if (rmc.mcStep(rmList, currentEnergy)) {
                currentEnergy = rmc.lastEnergy();
                ++nAccept;
                if (currentEnergy < optimumEnergy) {
                    optimumEnergy = currentEnergy;
                    System.arraycopy(initialRots, 0, optimum, 0, nRes);
                }
            } else {
                currentEnergy = rmc.lastEnergy();
                //++nReject;
            }
            /*if (threeBodyTerm) {
                logIfMaster(String.format(" %d Energy through 3-Body interactions: %16.8f",
                        i, currentEnergy));
            } else {
                logIfMaster(String.format(" %d Energy through 2-Body interactions: %16.8f",
                        i, currentEnergy));
            }*/
        }

        initTime += System.nanoTime();
        double fractAccept = ((double) nAccept) / ((double) maxIters);
        logIfMaster(String.format(" %d steps of DEE-MC completed in %10.6f seconds",
                maxIters, (initTime * 1.0E-9)));
        logIfMaster(String.format(" Number of steps accepted: %d for %10.6f of total", nAccept, fractAccept));
        logIfMaster(String.format(" Lowest energy found: %10.6f kcal/mol", optimumEnergy));
        logIfMaster(String.format(" Final energy found: %10.6f kcal/mol", currentEnergy));

        return optimumEnergy;
    }

    /**
     * Scrambles an array of rotamers.
     *
     * @param rotamers
     * @param residues
     * @param useAllElims
     */
    private void randomizeRotamers(int[] rotamers, Residue[] residues, boolean useAllElims) {
        int nRes = rotamers.length;
        for (int i = 0; i < nRes; i++) {
            Rotamer[] rotsi = residues[i].getRotamers(library);
            int lenri = rotsi.length;
            ArrayList<Integer> allowedRots = new ArrayList<>(lenri);

            for (int ri = 0; ri < lenri; ri++) {
                if (!check(i, ri)) {
                    allowedRots.add(ri);
                }
            }

            int nRots = allowedRots.size();
            if (nRots > 1) {
                boolean validMove = !useAllElims;
                int indexRI;
                do {
                    int ri = ThreadLocalRandom.current().nextInt(nRots);
                    indexRI = allowedRots.get(ri);
                    if (useAllElims) {
                        validMove = checkValidMove(i, indexRI, rotamers);
                    }
                } while (!validMove);
                rotamers[i] = indexRI;
            }
        }
    }

    /**
     * Checks the pair and triple elimination arrays to see if this permutation
     * has been eliminated. I am reasonably confident this violates detailed
     * balance.
     *
     * @param i Residue number
     * @param ri Rotamer number
     * @return If valid
     */
    private boolean checkValidMove(int i, int ri, int[] currentRots) {
        int nRes = currentRots.length;
        for (int j = 0; j < nRes; j++) {
            if (j == i) {
                continue;
            }
            if (check(j, currentRots[j], i, ri)) {
                return false;
            }
            for (int k = j + 1; k < nRes; k++) {
                if (k == i) {
                    continue;
                }
                if (check(j, currentRots[j], k, currentRots[k], i, ri)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * A global optimization over side-chain rotamers using a recursive
     * algorithm and information about eliminated rotamers, rotamer pairs and
     * rotamer triples
     *
     * @param molecularAssembly
     * @param residues
     * @param i
     * @param currentRotamers
     * @param lowEnergy
     * @param optimum Optimum set of rotamers.
     * @param permutationEnergies Energies of visited permutations or null.
     *
     * @return current energy.
     */
    public double rotamerOptimizationDEE(MolecularAssembly molecularAssembly, Residue residues[], int i,
            int currentRotamers[], double lowEnergy, int optimum[], double[] permutationEnergies) {
        // This is the initialization condition.
        if (i == 0) {
            evaluatedPermutations = 0;
        }
        int nResidues = residues.length;
        Residue residuei = residues[i];
        Rotamer[] rotamersi = residuei.getRotamers(library);
        int lenri = rotamersi.length;
        double currentEnergy = Double.MAX_VALUE;
        List<Residue> resList = Arrays.asList(residues);

        /**
         * As long as there are more residues, continue the recursion for each
         * rotamer of the current residue.
         */
        if (i < nResidues - 1) {
            /**
             * Loop over rotamers of residue i.
             */
            for (int ri = 0; ri < lenri; ri++) {
                /**
                 * Check if rotamer ri has been eliminated by DEE.
                 */
                if (check(i, ri)) {
                    continue;
                }
                /**
                 * Check if rotamer ri has been eliminated by an upstream
                 * rotamer (any residue's rotamer from j = 0 .. i-1).
                 */
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                applyRotamer(residuei, rotamersi[ri]);
                currentRotamers[i] = ri;
                double rotEnergy = rotamerOptimizationDEE(molecularAssembly, residues, i + 1,
                        currentRotamers, lowEnergy, optimum, permutationEnergies);
                if (rotEnergy < currentEnergy) {
                    currentEnergy = rotEnergy;
                }
                if (rotEnergy < lowEnergy) {
                    optimum[i] = ri;
                    lowEnergy = rotEnergy;
                }
            }
        } else {
            if (ensembleStates == null) {
                ensembleStates = new ArrayList<>();
            }

            /**
             * At the end of the recursion, compute the potential energy for
             * each rotamer of the final residue. If a lower potential energy is
             * discovered, the rotamers of each residue will be collected as the
             * recursion returns up the chain.
             */
            for (int ri = 0; ri < lenri; ri++) {
                /**
                 * Check if rotamer ri has been eliminated by DEE.
                 */
                if (check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                /**
                 * Check if rotamer ri has been eliminated by an upstream
                 * rotamer (any residue's rotamer from 0 .. i-1.
                 */
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                applyRotamer(residuei, rotamersi[ri]);
                // Compute the energy based on a 3-body approximation
                double approximateEnergy = computeEnergy(residues, currentRotamers, false);
                double comparisonEnergy = approximateEnergy;
                evaluatedPermutations++;
                // Compute the AMOEBA energy
                if (useFullAMOEBAEnergy) {
                    double amoebaEnergy = currentEnergy(resList);
                    if (permutationEnergies != null) {
                        permutationEnergies[evaluatedPermutations - 1] = amoebaEnergy;
                    }
                    comparisonEnergy = amoebaEnergy;
                    // Log current results
                    logIfMaster(String.format(" %d AMOEBA: %16.8f 3-Body: %16.8f Neglected: %16.8f",
                            evaluatedPermutations, amoebaEnergy, approximateEnergy, amoebaEnergy - approximateEnergy));
                } else {
                    if (permutationEnergies != null) {
                        permutationEnergies[evaluatedPermutations - 1] = approximateEnergy;
                    }
                    if (threeBodyTerm) {
                        logIfMaster(String.format(" %d Energy through 3-Body interactions: %16.8f",
                                evaluatedPermutations, approximateEnergy));
                    } else {
                        logIfMaster(String.format(" %d Energy through 2-Body interactions: %16.8f",
                                evaluatedPermutations, approximateEnergy));
                    }
                }

                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }

                if (ensembleNumber > 1) {
                    if (master && printFiles) {
                        try {
                            FileWriter fw = new FileWriter(ensembleFile, true);
                            BufferedWriter bw = new BufferedWriter(fw);
                            bw.write(String.format("MODEL        %d", evaluatedPermutations));
                            for (int j = 0; j < 75; j++) {
                                bw.write(" ");
                            }
                            bw.newLine();
                            bw.flush();
                            ensembleFilter.writeFile(ensembleFile, true);
                            bw.write(String.format("ENDMDL"));
                            for (int j = 0; j < 64; j++) {
                                bw.write(" ");
                            }
                            bw.newLine();
                            bw.close();
                        } catch (IOException e) {
                            logger.warning(String.format("Exception writing to file: %s", ensembleFile.getName()));
                        }
                    }
                    ResidueState[] states = ResidueState.storeAllCoordinates(residues);
                    ensembleStates.add(new ObjectPair<>(states, comparisonEnergy));
                }

                if (comparisonEnergy < currentEnergy) {
                    currentEnergy = comparisonEnergy;
                }

                if (comparisonEnergy < lowEnergy) {
                    lowEnergy = comparisonEnergy;
                    optimum[i] = ri;
                }
            }

            ensembleStates.sort(null);
        }
        return currentEnergy;
    }

    /**
     * A global optimization over side-chain rotamers using a recursive
     * algorithm and information about eliminated rotamers, rotamer pairs and
     * rotamer triples.
     *
     * @param residues
     * @param i
     * @param currentRotamers
     *
     * @return 0.
     */
    public double dryRun(Residue residues[], int i, int currentRotamers[]) {
        // This is the initialization condition.
        if (i == 0) {
            evaluatedPermutations = 0;
        }
        if (evaluatedPermutations >= 1e6) {
            if (evaluatedPermutations % 1000000 == 0) {
                logIfMaster(String.format("The permutation has reached %10.4e.", (double) evaluatedPermutations));
            }
        }
        int nResidues = residues.length;
        Residue residuei = residues[i];
        Rotamer[] rotamersi = residuei.getRotamers(library);
        int lenri = rotamersi.length;
        if (i < nResidues - 1) {
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                currentRotamers[i] = ri;
                dryRun(residues, i + 1, currentRotamers);
            }
        } else {
            /**
             * At the end of the recursion, check each rotamer of the final
             * residue.
             */
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                evaluatedPermutations++;
            }
        }
        return 0.0;
    }

    /**
     * Finds all permutations within buffer energy of GMEC.
     *
     * @param residues
     * @param i Current depth in residue/rotamer tree.
     * @param currentRotamers Current set of rotamers at this node.
     * @param gmecEnergy Minimum energy for these residues.
     * @param permutationEnergies Energy of all permutations.
     * @param permutations Contains accepted permutations.
     *
     * @return 0.
     */
    public double dryRunForEnsemble(Residue residues[], int i, int currentRotamers[],
            double gmecEnergy, double[] permutationEnergies, int[][] permutations) {
        // This is the initialization condition.
        if (i == 0) {
            evaluatedPermutations = 0;
        }
        int nResidues = residues.length;
        Residue residuei = residues[i];
        Rotamer[] rotamersi = residuei.getRotamers(library);
        int lenri = rotamersi.length;
        if (i < nResidues - 1) {
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                currentRotamers[i] = ri;
                dryRunForEnsemble(residues, i + 1,
                        currentRotamers, gmecEnergy, permutationEnergies, permutations);
            }
        } else {
            /**
             * At the end of the recursion, check each rotamer of the final
             * residue.
             */
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                if (permutationEnergies[evaluatedPermutations] - gmecEnergy < ensembleEnergy) {
                    permutations[evaluatedPermutations] = new int[nResidues];
                    System.arraycopy(currentRotamers, 0, permutations[evaluatedPermutations], 0, nResidues);
                }
                evaluatedPermutations++;
            }
        }
        return 0.0;
    }

    public void setThreeBodyEnergy(boolean threeBodyTerm) {
        this.threeBodyTerm = threeBodyTerm;
    }

    /**
     * Flag to indicate use of the full AMOEBA potential energy when ranking
     * rotamer combinations after dead-end elimination based on truncation at
     * 2-body or 3-body interactions.
     *
     * @param useFullAMOEBAEnergy
     */
    public void setUseFullAMOEBAEnergy(boolean useFullAMOEBAEnergy) {
        this.useFullAMOEBAEnergy = useFullAMOEBAEnergy;
    }

    /**
     * Uses existing backbone, self, pair, and 3-body energies from
     * rotamerEnergies() to calculate an approximate energy for a rotamer
     * permutation.
     *
     * @param residues Current window of optimization.
     * @param rotamers Set of rotamers to calculate an approximate energy for.
     * @param print Verbosity flag
     * @return Approximate permutation energy (backbone + selfs + pairs +
     * trimers).
     */
    private double computeEnergy(Residue residues[], int rotamers[], boolean print) {
        int nResidues = residues.length;
        double selfSum = 0.0;
        double pairSum = 0.0;
        double threeBodySum = 0.0;

        for (int a = 0; a < nResidues; a++) {
            int ai = rotamers[a];
            selfSum += self(a, ai);
            for (int b = a + 1; b < nResidues; b++) {
                int bi = rotamers[b];
                pairSum += pair(a, ai, b, bi);
                if (threeBodyTerm) {
                    for (int c = b + 1; c < nResidues; c++) {
                        int ci = rotamers[c];
                        threeBodySum += triple(a, ai, b, bi, c, ci);
                    }
                }
            }
        }

        double approximateEnergy = backboneEnergy + selfSum + pairSum + threeBodySum;
        if (print) {
            logger.info(String.format(" Backbone:                  %16.8f", backboneEnergy));
            logger.info(String.format(" Self Energy:               %16.8f", selfSum));
            logger.info(String.format(" Pair Energy:               %16.8f", pairSum));
            if (!threeBodyTerm) {
                logger.info(String.format(" Total Energy up to 2-Body: %16.8f", approximateEnergy));
            } else {
                logger.info(String.format(" 3-Body Energy:             %16.8f", threeBodySum));
                logger.info(String.format(" Total Energy up to 3-Body: %16.8f", approximateEnergy));
            }
        }
        return approximateEnergy;
    }

    /**
     * For decomposing the energies of the original (0th) rotamers without
     * computing the entire energy matrix.
     */
    public void decomposeOriginal() {
        allResiduesList = new ArrayList<>();
        polymers = molecularAssembly.getChains();
        for (Polymer polymer : polymers) {
            ArrayList<Residue> current = polymer.getResidues();
            for (Residue residuej : current) {
                if (residuej != null) {
                    if (residuej.getRotamers(library) != null) {
                        allResiduesList.add(residuej);
                    } else if (useForcedResidues && chain != null) {
                        Polymer setChain = molecularAssembly.getChain(chain);
                        if (setChain.equals(polymer) && checkIfForced(residuej)) {
                            allResiduesList.add(residuej);
                        }
                    } else if (useForcedResidues && checkIfForced(residuej)) {
                        allResiduesList.add(residuej);
                    }
                }
            }
        }
        numResidues = allResiduesList.size();
        allResiduesArray = allResiduesList.toArray(new Residue[numResidues]);
        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        int nRes = residues.length;
        distanceMatrix();

        double totalEnergy;
        double localBackboneEnergy = 0;
        double localSelfEnergy[] = new double[nRes];
        double pairEnergy[][] = new double[nRes][nRes];
        double triEnergy[][][] = new double[nRes][nRes][nRes];
        double residueEnergy[][] = new double[3][nRes];

        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
        totalEnergy = currentEnergy(allResiduesList);
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        // List to contain current residues.
        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));

        try {
            localBackboneEnergy = currentEnergy(rList);
        } catch (ArithmeticException ex) {
            logger.severe(String.format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            localSelfEnergy[i] = currentEnergy(rList) - localBackboneEnergy;
            logger.info(String.format(" Self %s:          %16.5f", ri, localSelfEnergy[i]));
            residueEnergy[0][i] = localSelfEnergy[i];
            turnOffAtoms(ri);
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                pairEnergy[i][j] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localBackboneEnergy;
                logger.info(String.format(" Pair %s %s:       %16.5f", ri, rj, pairEnergy[i][j]));
                double halfPair = pairEnergy[i][j] * 0.5;
                residueEnergy[1][i] += halfPair;
                residueEnergy[1][j] += halfPair;
                turnOffAtoms(rj);
            }
            turnOffAtoms(ri);
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                for (int k = j + 1; k < nRes; k++) {
                    Residue rk = residues[k];
                    rList.set(2, rk);
                    double dist = trimerDistance(i, 0, j, 0, k, 0);
                    if (dist < threeBodyCutoffDist) {
                        turnOnAtoms(rk);
                        triEnergy[i][j][k] = currentEnergy(rList)
                                - localSelfEnergy[i] - localSelfEnergy[j] - localSelfEnergy[k]
                                - pairEnergy[i][j] - pairEnergy[j][k] - pairEnergy[i][k] - localBackboneEnergy;
                        logger.info(String.format(" Tri  %s %s %s:    %16.5f", ri, rj, rk, triEnergy[i][j][k]));
                        double thirdTrimer = triEnergy[i][j][k] / 3.0;
                        residueEnergy[2][i] += thirdTrimer;
                        residueEnergy[2][j] += thirdTrimer;
                        residueEnergy[2][k] += thirdTrimer;
                        turnOffAtoms(rk);
                    } else if (dist == Double.MAX_VALUE) {
                        logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)", ri, rj, rk));
                        triEnergy[i][j][k] = 0.0;
                    } else {
                        logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms", ri, rj, rk, dist));
                        triEnergy[i][j][k] = 0.0;
                    }
                }
                turnOffAtoms(rj);
            }
            turnOffAtoms(ri);
        }

        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
        decomposePrint(residues, totalEnergy, localBackboneEnergy, residueEnergy);
    }

    /**
     * For decomposing the energies of the original (0th) rotamers without
     * computing the entire energy matrix; can accept a residue list.
     *
     * @param residues Residue array to decompose energies of.
     */
    public void decomposeOriginal(Residue[] residues) {
        allResiduesList = new ArrayList<>();
        polymers = molecularAssembly.getChains();
        for (Polymer polymer : polymers) {
            ArrayList<Residue> current = polymer.getResidues();
            for (Residue residuej : current) {
                if (residuej != null) {
                    if (residuej.getRotamers(library) != null) {
                        allResiduesList.add(residuej);
                    } else if (useForcedResidues && chain != null) {
                        Polymer setChain = molecularAssembly.getChain(chain);
                        if (setChain.equals(polymer) && checkIfForced(residuej)) {
                            allResiduesList.add(residuej);
                        }
                    } else if (useForcedResidues && checkIfForced(residuej)) {
                        allResiduesList.add(residuej);
                    }
                }
            }
        }
        numResidues = allResiduesList.size();
        allResiduesArray = allResiduesList.toArray(new Residue[numResidues]);
        int nRes = residues.length;
        distanceMatrix();

        double totalEnergy;
        double localBackboneEnergy = 0;
        double localSelfEnergy[] = new double[nRes];
        double pairEnergy[][] = new double[nRes][nRes];
        double triEnergy[][][] = new double[nRes][nRes][nRes];
        double residueEnergy[][] = new double[3][nRes];

        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
        totalEnergy = currentEnergy(allResiduesList);
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));

        try {
            localBackboneEnergy = currentEnergy(rList, false);
        } catch (ArithmeticException ex) {
            logger.severe(String.format(" FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(1, ri);
            turnOnAtoms(ri);
            localSelfEnergy[i] = currentEnergy(rList) - localBackboneEnergy;
            logger.info(String.format(" Self %s:          %16.5f", ri, localSelfEnergy[i]));
            residueEnergy[0][i] = localSelfEnergy[i];
            turnOffAtoms(ri);
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                pairEnergy[i][j] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localBackboneEnergy;
                logger.info(String.format(" Pair %s %s:       %16.5f", ri, rj, pairEnergy[i][j]));
                double halfPair = pairEnergy[i][j] * 0.5;
                residueEnergy[1][i] += halfPair;
                residueEnergy[1][j] += halfPair;
                turnOffAtoms(rj);
            }
            turnOffAtoms(ri);
        }

        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            int indexOfI = allResiduesList.indexOf(ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(1, rj);
                int indexOfJ = allResiduesList.indexOf(rj);
                turnOnAtoms(rj);
                for (int k = j + 1; k < nRes; k++) {
                    Residue rk = residues[k];
                    rList.set(2, rk);
                    int indexOfK = allResiduesList.indexOf(rk);
                    double dist = trimerDistance(indexOfI, 0, indexOfJ, 0, indexOfK, 0);
                    if (dist < threeBodyCutoffDist) {
                        turnOnAtoms(rk);
                        triEnergy[i][j][k] = currentEnergy(rList)
                                - localSelfEnergy[i] - localSelfEnergy[j] - localSelfEnergy[k]
                                - pairEnergy[i][j] - pairEnergy[j][k] - pairEnergy[i][k] - localBackboneEnergy;
                        logger.info(String.format(" Tri  %s %s %s:    %16.5f", ri, rj, rk, triEnergy[i][j][k]));
                        double thirdTrimer = triEnergy[i][j][k] / 3.0;
                        residueEnergy[2][i] += thirdTrimer;
                        residueEnergy[2][j] += thirdTrimer;
                        residueEnergy[2][k] += thirdTrimer;
                        turnOffAtoms(rk);
                    } else if (dist == Double.MAX_VALUE) {
                        logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)", ri, rj, rk));
                        triEnergy[i][j][k] = 0.0;
                    } else {
                        logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms", ri, rj, rk, dist));
                        triEnergy[i][j][k] = 0.0;
                    }
                }
                turnOffAtoms(rj);
            }
            turnOffAtoms(ri);
        }
        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
        decomposePrint(residues, totalEnergy, localBackboneEnergy, residueEnergy);
    }

    /**
     * For decomposing the energies of the original (0th) rotamers without
     * computing the entire energy matrix.
     */
    public void decomposeOriginalParallel() {
        allResiduesList = new ArrayList<>();
        polymers = molecularAssembly.getChains();
        for (Polymer polymer : polymers) {
            ArrayList<Residue> current = polymer.getResidues();
            for (Residue residuej : current) {
                if (residuej != null) {
                    if (residuej.getRotamers(library) != null) {
                        allResiduesList.add(residuej);
                    } else if (useForcedResidues && chain != null) {
                        Polymer setChain = molecularAssembly.getChain(chain);
                        if (setChain.equals(polymer) && checkIfForced(residuej)) {
                            allResiduesList.add(residuej);
                        }
                    } else if (useForcedResidues && checkIfForced(residuej)) {
                        allResiduesList.add(residuej);
                    }
                }
            }
        }
        numResidues = allResiduesList.size();
        allResiduesArray = allResiduesList.toArray(new Residue[numResidues]);
        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        int nRes = residues.length;
        distanceMatrix();

        double totalEnergy;
        double localBackboneEnergy = 0;
        double localSelfEnergy[] = new double[nRes];
        double pairEnergy[][] = new double[nRes][nRes];
        double triEnergy[][][] = new double[nRes][nRes][nRes];
        double residueEnergy[][] = new double[3][nRes];

        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
        totalEnergy = currentEnergy(allResiduesList);
        logIfMaster(String.format(" AMOEBA:   %16.5f", totalEnergy));
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));
        try {
            localBackboneEnergy = currentEnergy(rList, false);
        } catch (ArithmeticException ex) {
            logger.severe(String.format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        logIfMaster(String.format(" Backbone: %16.5f", localBackboneEnergy));

        decomposeOriginal = true;
        allocateEliminationMemory(allResiduesArray);
        rotamerEnergies(allResiduesArray);

        if (master) {
            for (int i = 0; i < nRes; i++) {
                Residue ri = residues[i];
                localSelfEnergy[i] = selfEnergy[i][0];
                residueEnergy[0][i] = localSelfEnergy[i];
                logger.info(String.format(" Self %s:          %16.5f", ri, localSelfEnergy[i]));
            }
            for (int i = 0; i < nRes; i++) {
                Residue ri = residues[i];
                for (int j = i + 1; j < nRes; j++) {
                    Residue rj = residues[j];
                    pairEnergy[i][j] = twoBodyEnergy[i][0][j][0];
                    logger.info(String.format(" Pair %s %s:       %16.5f", ri, rj, pairEnergy[i][j]));
                    double halfPair = pairEnergy[i][j] * 0.5;
                    residueEnergy[1][i] += halfPair;
                    residueEnergy[1][j] += halfPair;
                }
            }
            if (threeBodyTerm) {
                for (int i = 0; i < nRes; i++) {
                    Residue ri = residues[i];
                    for (int j = i + 1; j < nRes; j++) {
                        Residue rj = residues[j];
                        for (int k = j + 1; k < nRes; k++) {
                            Residue rk = residues[k];
                            double dist = trimerDistance(i, 0, j, 0, k, 0);
                            triEnergy[i][j][k] = threeBodyEnergy[i][0][j][0][k][0];
                            double thirdTrimer = triEnergy[i][j][k] / 3.0;
                            residueEnergy[2][i] += thirdTrimer;
                            residueEnergy[2][j] += thirdTrimer;
                            residueEnergy[2][k] += thirdTrimer;
                            if (triEnergy[i][j][k] != 0.0) {
                                logger.info(String.format(" Tri  %s %s %s:    %16.5f", ri, rj, rk, triEnergy[i][j][k]));
                            } else if (dist == Double.MAX_VALUE) {
                                logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)",
                                        ri, rj, rk));
                                triEnergy[i][j][k] = 0.0;
                            } else if (dist > threeBodyCutoffDist) {
                                logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms",
                                        ri, rj, rk, dist));
                                triEnergy[i][j][k] = 0.0;
                            } else {
                                String m = String.
                                        format(" Zero trimer energy inside cutoff: %s %s %s at %1.5f Angstroms.",
                                                ri, rj, rk, dist);
                                logger.warning(m);
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < nRes; i++) {
                turnOnAtoms(residues[i]);
            }
            decomposePrint(residues, totalEnergy, localBackboneEnergy, residueEnergy);
        }
        decomposeOriginal = false;
    }

    /**
     * Method intended to decompose energies down to quad energies. Mostly for
     * showing that quads are probably negligible.
     *
     * @param quadCutoff
     * @param maxQuads
     */
    public void decomposeOriginalQuads(double quadCutoff, int maxQuads) {
        allResiduesList = new ArrayList<>();
        polymers = molecularAssembly.getChains();
        for (Polymer polymer : polymers) {
            ArrayList<Residue> current = polymer.getResidues();
            for (Residue residuej : current) {
                if (residuej != null) {
                    if (residuej.getRotamers(library) != null) {
                        allResiduesList.add(residuej);
                    } else if (useForcedResidues && chain != null) {
                        Polymer setChain = molecularAssembly.getChain(chain);
                        if (setChain.equals(polymer) && checkIfForced(residuej)) {
                            allResiduesList.add(residuej);
                        }
                    } else if (useForcedResidues && checkIfForced(residuej)) {
                        allResiduesList.add(residuej);
                    }
                }
            }
        }
        numResidues = allResiduesList.size();
        allResiduesArray = allResiduesList.toArray(new Residue[numResidues]);
        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        int nRes = residues.length;
        distanceMatrix();

        double totalEnergy;
        double localBackboneEnergy = 0;
        double sumSelf = 0;
        double sumPair = 0;
        double sumTri = 0;
        double sumQuads = 0;
        double localSelfEnergy[] = new double[nRes];
        double pairEnergy[][] = new double[nRes][nRes];
        double triEnergy[][][] = new double[nRes][nRes][nRes];
        //double quadEnergy[][][][] = new double[nRes][][][]; This array is gigantic and unnecessary.

        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
        totalEnergy = currentEnergy(allResiduesList);
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));
        try {
            localBackboneEnergy = currentEnergy(rList, false);
        } catch (ArithmeticException ex) {
            logger.severe(String.format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            localSelfEnergy[i] = currentEnergy(rList) - localBackboneEnergy;
            logger.info(String.format(" Self %s:          %16.5f", ri, localSelfEnergy[i]));
            sumSelf += localSelfEnergy[i];
            turnOffAtoms(ri);
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(0, rj);
                turnOnAtoms(rj);
                pairEnergy[i][j] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localBackboneEnergy;
                logger.info(String.format(" Pair %s %s:       %16.5f", ri, rj, pairEnergy[i][j]));
                sumPair += pairEnergy[i][j];
                turnOffAtoms(rj);
            }
            turnOffAtoms(ri);
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                for (int k = j + 1; k < nRes; k++) {
                    Residue rk = residues[k];
                    rList.set(2, rk);
                    double dist = trimerDistance(i, 0, j, 0, k, 0);
                    if (dist < threeBodyCutoffDist) {
                        turnOnAtoms(rk);
                        triEnergy[i][j][k] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localSelfEnergy[k]
                                - pairEnergy[i][j] - pairEnergy[j][k] - pairEnergy[i][k] - localBackboneEnergy;
                        logger.info(String.format(" Tri  %s %s %s:    %16.5f", ri, rj, rk, triEnergy[i][j][k]));
                        sumTri += triEnergy[i][j][k];
                        turnOffAtoms(rk);
                    } else if (dist == Double.MAX_VALUE) {
                        logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)", ri, rj, rk));
                        triEnergy[i][j][k] = 0.0;
                    } else {
                        logger.info(String.format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms", ri, rj, rk, dist));
                        triEnergy[i][j][k] = 0.0;
                    }
                }
                turnOffAtoms(rj);
            }
            turnOffAtoms(ri);
        }

        int numQuadsEvaluated = 0;
        boolean doBreak = false;
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            //quadEnergy[i] = new double[nRes][][];
            // If for some reason storing quad energies is desired, one can allocate memory on the fly, so that only enough
            // memory is allocated for the quads you actually evaluate.
            for (int j = i + 1; j < nRes; j++) {
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                //quadEnergy[i][j] = new double[nRes][];
                for (int k = j + 1; k < nRes; k++) {
                    Residue rk = residues[k];
                    rList.set(2, rk);
                    turnOnAtoms(rk);
                    //quadEnergy[i][j][k] = new double[nRes];
                    for (int l = k + 1; l < nRes; l++) {
                        double dist = quadDistance(i, 0, j, 0, k, 0, l, 0);
                        Residue rl = residues[l];
                        rList.set(3, rl);
                        if (dist < quadCutoff) {
                            turnOnAtoms(rl);
                            /*quadEnergy[i][j][k][l] = currentEnergy() - localSelfEnergy[i] - localSelfEnergy[j]
                             - localSelfEnergy[k] - localSelfEnergy[l] - pairEnergy[i][j] - pairEnergy[i][k] -
                             pairEnergy[i][l] - pairEnergy[j][k] - pairEnergy[j][l] - pairEnergy[k][l] -
                             triEnergy[i][j][k] - triEnergy[i][j][l] - triEnergy[i][k][l] - triEnergy[j][k][l] -
                             localBackboneEnergy;*/
                            double currentQuad = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j]
                                    - localSelfEnergy[k] - localSelfEnergy[l] - pairEnergy[i][j] - pairEnergy[i][k]
                                    - pairEnergy[i][l] - pairEnergy[j][k] - pairEnergy[j][l] - pairEnergy[k][l]
                                    - triEnergy[i][j][k] - triEnergy[i][j][l] - triEnergy[i][k][l] - triEnergy[j][k][l]
                                    - localBackboneEnergy;
                            logger.info(String.format(" Quad  %s %s %s %s:    %16.5f at %1.5f Angstroms",
                                    ri, rj, rk, rl, currentQuad, dist));
                            sumQuads += currentQuad;
                            turnOffAtoms(rl);
                            if (++numQuadsEvaluated >= maxQuads) {
                                doBreak = true;
                                break;
                            }
                        } else if (dist == Double.MAX_VALUE) {
                            logger.info(String.format(" Quad  %s %s %s %s:    set to 0.0 at NaN (very long distance)",
                                    ri, rj, rk, rl));
                        } else {
                            logger.info(String.format(" Quad  %s %s %s %s:    set to 0.0 at %1.5f Angstroms",
                                    ri, rj, rk, rl, dist));
                        }
                    }
                    turnOffAtoms(rk);
                    if (doBreak) {
                        break;
                    }
                }
                turnOffAtoms(rj);
                if (doBreak) {
                    break;
                }
            }
            turnOffAtoms(ri);
            if (doBreak) {
                break;
            }
        }

        logger.info(String.format("\n\n"));
        logger.info(String.format(" Backbone:     %16.5f", localBackboneEnergy));
        logger.info(String.format(" Sum Self:     %16.5f", sumSelf));
        logger.info(String.format(" Sum Pair:     %16.5f", sumPair));
        logger.info(String.format(" Sum Tri:      %16.5f", sumTri));
        logger.info(String.format(" Sum Quad:     %16.5f", sumQuads));
        logger.info(String.format(" Neglected:    %16.5f", totalEnergy - sumSelf - sumPair - sumTri - sumQuads - localBackboneEnergy));
        logger.info(String.format(" AMOEBA:       %16.5f", totalEnergy));
        logger.info(String.format("\n\n"));
        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
    }

    private void decomposePrint(Residue residues[], double totalEnergy, double backbone, double residueEnergy[][]) {
        if (residues == null || residueEnergy == null) {
            return;
        }

        int molecule[] = molecularAssembly.getMoleculeNumbers();

        // Internal Molecule ID to Reduced Set of Molecule IDs
        HashMap<Integer, Integer> nMolecules = new HashMap<>();
        // Internal Molecule ID to Molecule Names
        HashMap<Integer, String> molNames = new HashMap<>();
        // Residue ID to Reduced Set of Molecule IDs
        HashMap<Residue, Integer> moleculeMap = new HashMap<>();

        /**
         * Sum self, pair and trimer energies.
         */
        int nRes = residues.length;
        double sumSelf = 0;
        double sumPair = 0;
        double sumTri = 0;
        int molPointer = 0;
        int reducedIndex = 0;
        for (int i = 0; i < nRes; i++) {
            sumSelf += residueEnergy[0][i];
            sumPair += residueEnergy[1][i];
            sumTri += residueEnergy[2][i];

            // Count the number of molecules.
            Residue r = residues[i];
            Atom a0 = (Atom) r.getAtomNode(0);
            int atomIndex = a0.getIndex();
            Integer moleculeIndex = molecule[atomIndex];
            if (!nMolecules.containsKey(moleculeIndex)) {
                nMolecules.put(moleculeIndex, molPointer);
                molNames.put(moleculeIndex, r.getSegID());
                reducedIndex = molPointer;
                molPointer++;
            } else {
                reducedIndex = nMolecules.get(moleculeIndex);
            }
            moleculeMap.put(r, reducedIndex);
        }

        int nMol = nMolecules.size();
        if (nMol > 1) {
            logger.info(String.format("\n"));
            logger.info(String.format(" Molecule-Based Many-Body Energy Summation\n "));
            double moleculeEnergy[][] = new double[3][nMol];
            for (int i = 0; i < nRes; i++) {
                int iMolecule = moleculeMap.get(residues[i]);
                moleculeEnergy[0][iMolecule] += residueEnergy[0][i];
                moleculeEnergy[1][iMolecule] += residueEnergy[1][i];
                moleculeEnergy[2][iMolecule] += residueEnergy[2][i];
            }
            logger.info(String.format(" %9s %9s %9s %9s %9s", "Molecule", "Self", "Pair", "3-Body", "Total"));
            for (int i = 0; i < nMol; i++) {
                String molName = molNames.get(i);
                double total = moleculeEnergy[0][i] + moleculeEnergy[1][i] + moleculeEnergy[2][i];
                logger.info(String.format(" %9s %9.3f %9.3f %9.3f %9.3f",
                        molName, moleculeEnergy[0][i], moleculeEnergy[1][i], moleculeEnergy[2][i], total));
            }
            logger.info(String.format(" %9s %9.3f %9.3f %9.3f %9.3f",
                    "Sum", sumSelf, sumPair, sumTri, sumSelf + sumPair + sumTri));
        }

        logger.info(String.format("\n"));
        logger.info(String.format(" Residue-Based Many-Body Energy Summation\n "));
        logger.info(String.format(" %9s %9s %9s %9s %9s", "Residue", "Self", "Pair", "3-Body", "Total"));

        for (int i = 0; i < nRes; i++) {
            double total = residueEnergy[0][i] + residueEnergy[1][i] + residueEnergy[2][i];
            logger.info(String.format(" %9s %9.3f %9.3f %9.3f %9.3f",
                    residues[i].toString(), residueEnergy[0][i], residueEnergy[1][i], residueEnergy[2][i], total));
        }
        logger.info(String.format(" %9s %9.3f %9.3f %9.3f %9.3f",
                "Sum", sumSelf, sumPair, sumTri, sumSelf + sumPair + sumTri));
        logger.info(String.format(" Backbone:        %9.3f", backbone));
        double target = sumSelf + sumPair + sumTri + backbone;
        logger.info(String.format(" Expansion Total: %9.3f", target));
        logger.info(String.format(" True Total:      %9.3f", totalEnergy));
        logger.info(String.format(" Neglected:       %9.3f\n", totalEnergy - target));
    }

    /**
     * Prints a summary of pair and trimer energies above [cutoff] kcal/mol.
     *
     * @param pairCutoff
     * @param trimerCutoff
     * @param dOMode
     */
    public void printLargeInteractions(double pairCutoff, double trimerCutoff, boolean dOMode) {
        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        int nRes = residues.length;

        if (dOMode) {
            logger.info(String.format(" Large pair interactions (>%.2f):", pairCutoff));
            for (int i = 0; i < nRes; i++) {
                for (int j = i + 1; j < nRes; j++) {
                    if (Math.abs(twoBodyEnergy[i][0][j][0]) >= pairCutoff) {
                        logger.info(String.format(" Large Pair %s %s:       %16.5f",
                                residues[i], residues[j], twoBodyEnergy[i][0][j][0]));
                    }
                }
            }
            logger.info(String.format("\n Large trimer interactions (>%.2f):", trimerCutoff));
            for (int i = 0; i < nRes; i++) {
                for (int j = i + 1; j < nRes; j++) {
                    for (int k = j + 1; k < nRes; k++) {
                        if (Math.abs(threeBodyEnergy[i][0][j][0][k][0]) >= trimerCutoff) {
                            logger.info(String.format(" Large Trimer  %s %s %s:    %16.5f",
                                    residues[i], residues[j], residues[k], threeBodyEnergy[i][0][j][0][k][0]));
                        }
                    }
                }
            }
            return;
        }

        logger.info(String.format(" Large pair interactions (>%.2f):", pairCutoff));
        for (int i = 0; i < nRes; i++) {
            Residue resi = residues[i];
            Rotamer roti[] = resi.getRotamers(library);
            for (int ri = 0; ri < roti.length; ri++) {
                for (int j = i + 1; j < nRes; j++) {
                    Residue resj = residues[j];
                    Rotamer rotj[] = resj.getRotamers(library);
                    for (int rj = 0; rj < rotj.length; rj++) {
                        try {
                            if (Math.abs(twoBodyEnergy[i][ri][j][rj]) >= pairCutoff) {
                                logger.info(String.format(" Large Pair %7s %-2d, %7s %-2d: %16.8f",
                                        resi, ri, resj, rj, twoBodyEnergy[i][ri][j][rj]));
                            }
                        } catch (Exception ex) {
                        }
                    }
                }
            }
        }

        logger.info(String.format("\n Large trimer interactions (>%.2f):", trimerCutoff));
        for (int i = 0; i < nRes; i++) {
            Residue resi = residues[i];
            Rotamer roti[] = resi.getRotamers(library);
            for (int ri = 0; ri < roti.length; ri++) {
                for (int j = i + 1; j < nRes; j++) {
                    Residue resj = residues[j];
                    Rotamer rotj[] = resj.getRotamers(library);
                    for (int rj = 0; rj < rotj.length; rj++) {
                        for (int k = j + 1; k < nRes; k++) {
                            Residue resk = residues[k];
                            Rotamer rotk[] = resk.getRotamers(library);
                            for (int rk = 0; rk < rotk.length; rk++) {
                                try {
                                    if (Math.abs(threeBodyEnergy[i][ri][j][rj][k][rk]) >= trimerCutoff) {
                                        logger.info(String.format(" Large Trimer %7s %-2d, %7s %-2d, %7s %-2d: %16.8f",
                                                resi, ri, resj, rj, resk, rk, threeBodyEnergy[i][ri][j][rj][k][rk]));
                                    }
                                } catch (Exception ex) {
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Set a contiguous block of residues to optimize in a specific chain.
     *
     * @param chain
     * @param startResID
     * @param finalResID
     */
    public void setResidues(String chain, int startResID, int finalResID) {
        this.chain = chain;
        this.setResidues(startResID, finalResID);
    }

    /**
     * Set a contiguous block of residues to optimize.
     *
     * @param startResID
     * @param finalResID
     */
    public void setResidues(int startResID, int finalResID) {

        Polymer polymer;
        if (chain != null) {
            polymer = molecularAssembly.getChain(chain);
        } else {
            polymers = molecularAssembly.getChains();
            polymer = polymers[0];
        }
        residueList = new ArrayList<>();
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymer.getResidue(i);
            if (residue != null) {
                Rotamer[] rotamers = residue.getRotamers(library);
                if (rotamers != null) {
                    if (rotamers.length == 1) {
                        switch (residue.getResidueType()) {
                            case NA:
                                residue.initializeDefaultAtomicCoordinates();
                                break;
                            case AA:
                            default:
                                RotamerLibrary.applyRotamer(residue, rotamers[0]);
                                break;
                        }
                        if (addOrigRot) {
                            residueList.add(residue);
                        }
                    } else {
                        residueList.add(residue);
                    }
                } else if (useForcedResidues && checkIfForced(i)) {
                    residueList.add(residue);
                }
            }
        }
    }

    public ArrayList<Residue> getResidues() {
        return residueList;
    }

    public void setResidues(ArrayList<Residue> residueList) {
        this.residueList = residueList;
    }

    public void setCoordinatesToEnsemble(int ensnum) {
        if (ensembleStates != null && !ensembleStates.isEmpty()) {
            ensnum %= ensembleStates.size();
            ResidueState.revertAllCoordinates(residueList, ensembleStates.get(ensnum).getVal());
            //ResidueState.revertAllCoordinates(residueList, ensembleStates.get(ensnum));
        } else {
            throw new IllegalArgumentException(" Ensemble states not initialized!");
        }
    }

    public List<ResidueState[]> getEnsemble() {
        if (ensembleStates == null) {
            return null;
        } else {
            List<ResidueState[]> states = new ArrayList<>(ensembleStates.size());
            ensembleStates.forEach((es) -> {
                states.add(es.getVal());
            });
            return states;
        }
    }

    /**
     * Accepts a list of residues but throws out null residues. Used by the -lR
     * flag.
     *
     * @param residueList
     * @return Added residues.
     */
    public List<Residue> setResiduesIgnoreNull(ArrayList<Residue> residueList) {
        this.residueList = new ArrayList<>();
        logger.info(" Optimizing these residues: ");
        for (Residue r : residueList) {
            if (r.getRotamers(library) != null) {
                this.residueList.add(r);
                logger.info(String.format("\t%s", r.toString()));
            } else {
                logger.info(String.format(" not \t%s", r.toString()));
            }
        }
        return new ArrayList<>(residueList);
    }

    public double optimize() {
        boolean ignoreNA = false;
        String ignoreNAProp = System.getProperty("ignoreNA");
        if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
            logger.info("(Key) Ignoring nucleic acids.");
            ignoreNA = true;
        }
        /*
         * Collect all residues in the MolecularAssembly. Use all Residues with
         * Rotamers and all forced residues if using sliding window and forced
         * residues. Forced residues is meaningless for other algorithms, and
         * will be reverted to false.
         */
        allResiduesList = new ArrayList<>();
        polymers = molecularAssembly.getChains();
        for (Polymer polymer : polymers) {
            ArrayList<Residue> current = polymer.getResidues();
            for (Residue residuej : current) {
                if (useForcedResidues) {
                    switch (algorithm) {
                        case SLIDING_WINDOW:
                            if (residuej == null) {
                                // Do nothing.
                            } else if (residuej.getRotamers(library) != null) {
                                allResiduesList.add(residuej);
                            } else {
                                int indexJ = residuej.getResidueNumber();
                                if (checkIfForced(indexJ)) {
                                    allResiduesList.add(residuej);
                                }
                            }
                            break;
                        default: // Should only trigger once before resetting useForcedResidues to false.
                            logIfMaster(" Forced residues only applicable to sliding window.", Level.WARNING);
                            useForcedResidues = false;
                            if (residuej != null && (residuej.getRotamers(library) != null)) {
                                if (!(ignoreNA && residuej.getResidueType() == Residue.ResidueType.NA)) {
                                    allResiduesList.add(residuej);
                                }
                            }
                            break;
                    }
                } else if (residuej != null && (residuej.getRotamers(library) != null)) {
                    if (!(ignoreNA && residuej.getResidueType() == Residue.ResidueType.NA)) {
                        allResiduesList.add(residuej);
                    }
                }
            }
        }

        // If -DignoreNA=true, then remove nucleic acids from residue list.
        if (ignoreNA) {
            for (int i = 0; i < residueList.size(); i++) {
                Residue res = residueList.get(i);
                if (res.getResidueType() == Residue.ResidueType.NA) {
                    residueList.remove(i);
                }
            }
        }

        RotamerLibrary.initializeDefaultAtomicCoordinates(molecularAssembly.getChains());   // for NA only
        numResidues = allResiduesList.size();
        allResiduesArray = allResiduesList.toArray(new Residue[numResidues]);

        /*
         * Distance matrix is  used to add residues to the sliding window
         * based on distance cutoff, and to automatically set some 3-body terms
         * to 0 at > 10 angstroms.
         *
         * The memory and compute overhead can be a problem for some very large
         * structures.
         */
        if (distance > 0) {
            distanceMatrix();
        }

        double e = 0.0;
        if (residueList != null) {
            done = false;
            terminate = false;
            switch (algorithm) {
                case INDEPENDENT:
                    e = independent(residueList);
                    break;
                case GLOBAL:
                    e = globalBruteForce(residueList);
                    break;
                case GLOBAL_DEE:
                    e = globalUsingEliminations(residueList);
                    break;
                case SLIDING_WINDOW:
                    e = slidingWindow(residueList, windowSize, increment, revert, distance, direction);
                    break;
                case BOX_OPTIMIZATION:
                    e = boxOptimization(residueList);
                    break;
                default:
                    break;
            }
            terminate = false;
            done = true;
        }
        return e;
    }

    /**
     * Eliminates NA backbone rotamers with corrections greater than threshold.
     * The int[] parameter allows the method to know how many Rotamers for each
     * residue have previously been pruned; currently, this means any Rotamer
     * pruned by reconcileNARotamersWithPriorResidues.
     *
     * A nucleic correction threshold of 0 skips the entire method; this check
     * is presently being performed inside the method in case it is called again
     * at some point.
     *
     * @param residues Residues to eliminate bad backbone rotamers over.
     * @param numEliminatedRotamers Number of previously eliminated rotamers per
     * residue.
     */
    private void eliminateNABackboneRotamers(Residue[] residues, int[] numEliminatedRotamers) {
        /**
         * Initialize all atoms to be used.
         */
        /* Atom atoms[] = molecularAssembly.getAtomArray();
         int nAtoms = atoms.length;
         String begin[] = new String[nAtoms];
         for (int i = 0; i < nAtoms; i++) {
         begin[i] = atoms[i].toString();
         } */
        if (nucleicCorrectionThreshold != 0) {
            logIfMaster(String.format(" Eliminating nucleic acid rotamers with correction vectors larger than %5.3f A", nucleicCorrectionThreshold));
            logIfMaster(String.format(" A minimum of %d rotamers per NA residue will carry through to energy calculations.", minNumberAcceptedNARotamers));
            ArrayList<Residue> resList = new ArrayList<>();
            resList.addAll(Arrays.asList(residues));
            ResidueState[] origCoordinates = ResidueState.storeAllCoordinates(resList);
            //double[][][] originalCoordinates = storeCoordinates(resList);
            for (int j = 0; j < residues.length; j++) {
                Residue nucleicResidue = residues[j];
                Rotamer[] rotamers = nucleicResidue.getRotamers(library);
                if (nucleicResidue.getResidueType() == NA && rotamers != null) {
                    int nrotamers = rotamers.length;
                    // Default to all rotamers that have not previously been
                    // eliminated; subtract as rotamers are rejected.
                    int numAcceptedRotamers = nrotamers - numEliminatedRotamers[j];
                    if (minNumberAcceptedNARotamers >= numAcceptedRotamers) {
                        continue;
                    }
                    ArrayList<DoubleIndexPair> rejectedRotamers = new ArrayList<>();
                    for (int i = 0; i < nrotamers; i++) {
                        if (!check(j, i)) {
                            try {
                                RotamerLibrary.applyRotamer(nucleicResidue, rotamers[i], nucleicCorrectionThreshold);
                            } catch (NACorrectionException error) {
                                double rejectedCorrection = error.getCorrection();
                                numAcceptedRotamers--;
                                DoubleIndexPair rejected = new DoubleIndexPair(i, rejectedCorrection);
                                rejectedRotamers.add(rejected);
                            }
                        }
                    }
                    int numAdditionalRotamersToAccept = minNumberAcceptedNARotamers - numAcceptedRotamers;
                    if (numAdditionalRotamersToAccept > 0) {
                        DoubleIndexPair[] rejectedArray = new DoubleIndexPair[rejectedRotamers.size()];
                        for (int i = 0; i < rejectedArray.length; i++) {
                            rejectedArray[i] = rejectedRotamers.get(i);
                        }
                        Arrays.sort(rejectedArray);
                        rejectedRotamers = new ArrayList<>();
                        rejectedRotamers.addAll(Arrays.asList(rejectedArray));
                        for (int i = 0; i < numAdditionalRotamersToAccept; i++) {
                            rejectedRotamers.remove(0);
                        }
                    }
                    for (DoubleIndexPair rotToReject : rejectedRotamers) {
                        eliminateRotamer(residues, j, rotToReject.getIndex(), print);
                        logIfMaster(String.format(" Correction magnitude was %6.4f A > %5.3f A", rotToReject.getDoubleValue(), nucleicCorrectionThreshold));
                    }
                }
                nucleicResidue.revertState(origCoordinates[j]);
                //revertSingleResidueCoordinates(nucleicResidue, originalCoordinates[j]);
                /*
                 for (int i = 0; i < nAtoms; i++) {
                 if (!begin[i].equals(atoms[i].toString())) {
                 logger.info(" ERROR IN NA BACKBONE: " + begin[i]);
                 logger.info(" ERROR IN NA BACKBONE: " + atoms[i].toString());
                 logger.severe(" Good Bye");
                 }
                 } */
            }
        }
    }

    /**
     * For NA residues inside some optimization window, prune any rotamers which
     * would be incompatible with the established rotamers of upstream NA
     * residues. Could in theory be done by self energies, but every rotamer
     * which can be eliminated without calculating a self energy makes the
     * optimization significantly faster.
     *
     * @param residues Residues to check for incompatible rotamers.
     * @return Number of rotamers eliminated for each Residue.
     */
    private int[] reconcileNARotamersWithPriorResidues(Residue[] residues) {
        int[] numEliminatedRotamers = new int[residues.length];
        for (int i = 0; i < residues.length; i++) {
            Residue residuei = residues[i];
            Rotamer[] rotamers = residuei.getRotamers(library);
            if (rotamers == null || residuei.getResidueType() != NA) {
                continue;
            }
            Residue prevResidue = residuei.getPreviousResidue();
            if (prevResidue == null || prevResidue.getResidueType() != NA) {
                continue;
            }
            boolean isInList = false;
            for (Residue residue : residues) {
                if (prevResidue.equals(residue)) {
                    isInList = true;
                    break;
                }
            }
            if (isInList) {
                continue;
            }
            double prevDelta = RotamerLibrary.measureDelta(prevResidue);
            // If between 50 and 110, assume a North pucker.
            if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                for (int j = 0; j < rotamers.length; j++) {
                    if (RotamerLibrary.checkPucker(rotamers[j].chi1) != 1) {
                        if (print) {
                            logIfMaster(String.format(" Rotamer %d of residue %s eliminated "
                                    + "for incompatibility with the sugar pucker of previous "
                                    + "residue %s outside the window.", j, residuei.toString(),
                                    prevResidue.toString()));
                        }
                        eliminateRotamer(residues, i, j, print);
                        numEliminatedRotamers[i]++;
                    }
                }
            } else {
                for (int j = 0; j < rotamers.length; j++) {
                    if (RotamerLibrary.checkPucker(rotamers[j].chi1) != 2) {
                        if (print) {
                            logIfMaster(String.format(" Rotamer %d of residue %s eliminated "
                                    + "for incompatibility with the sugar pucker of previous "
                                    + "residue %s outside the window.", j, residuei.toString(),
                                    prevResidue.toString()));
                        }
                        eliminateRotamer(residues, i, j, print);
                        numEliminatedRotamers[i]++;
                    }
                }
            } // TODO: Implement support for the DNA C3'-exo pucker.
        }
        return numEliminatedRotamers;
    }

    /**
     * For NA residues inside some optimization window, prune any rotamers which
     * would be incompatible with the established rotamers of downstream NA
     * residues. Could in theory be done by self energies, but every rotamer
     * which can be eliminated without calculating a self energy makes the
     * optimization much faster. Technically, this works by pinning it to its
     * current pucker, but if the input structure is any good, it would be the
     * current pucker anyways.
     *
     * @param residues Residues to check for incompatible rotamers.
     * @return Number of rotamers eliminated for each Residue.
     */
    private void reconcileNARotamersWithSubsequentResidues(Residue[] residues, int[] eliminatedRotamers) {
        for (int i = 0; i < residues.length; i++) {
            Residue residuei = residues[i];
            switch (residuei.getResidueType()) {
                case NA:
                    Rotamer[] rotamers = residues[i].getRotamers(library);
                    Residue nextResidue = residuei.getNextResidue();
                    if (rotamers == null || nextResidue == null) {
                        break;
                    }
                    boolean isInList = false;
                    for (Residue residue : residues) {
                        if (residue.equals(nextResidue)) {
                            isInList = true;
                            break;
                        }
                    }
                    if (isInList) {
                        break;
                    }
                    double delta = RotamerLibrary.measureDelta(residuei);
                    if (RotamerLibrary.checkPucker(delta) == 1) {
                        for (int j = 0; j < rotamers.length; j++) {
                            if (RotamerLibrary.checkPucker(rotamers[j].chi7) != 1) {
                                if (print) {
                                    logIfMaster(String.format(" Rotamer %d of residue %s eliminated "
                                            + "for incompatibility with the sugar pucker of previous "
                                            + "residue %s outside the window.", j, residuei.toString(),
                                            nextResidue.toString()));
                                }
                                eliminateRotamer(residues, i, j, print);
                                eliminatedRotamers[i]++;
                            }
                        }
                    } else {
                        for (int j = 0; j < rotamers.length; j++) {
                            if (RotamerLibrary.checkPucker(rotamers[j].chi7) != 2) {
                                if (print) {
                                    logIfMaster(String.format(" Rotamer %d of residue %s eliminated "
                                            + "for incompatibility with the sugar pucker of previous "
                                            + "residue %s outside the window.", j, residuei.toString(),
                                            nextResidue.toString()));
                                }
                                eliminateRotamer(residues, i, j, print);
                                eliminatedRotamers[i]++;
                            }
                        }
                    }
                    break;
                case AA:
                default:
                    break;
            }
        }
    }

    /**
     * Sets the number of boxes in the x, y, and z axes if the box optimization
     * is to be carried out.
     *
     * @param numXYZBoxes Int[3] of number of boxes in x, y, z.
     */
    public void setNumXYZBoxes(int[] numXYZBoxes) {
        System.arraycopy(numXYZBoxes, 0, this.numXYZBoxes, 0, this.numXYZBoxes.length);
    }

    /**
     * Sets the amount of overlap between adjacent boxes for box optimization.
     *
     * @param boxBorderSize Box overlap in Angstroms.
     */
    public void setBoxBorderSize(double boxBorderSize) {
        this.boxBorderSize = boxBorderSize;
    }

    /**
     * Sets the approximate dimensions of boxes, over-riding numXYZBoxes in
     * determining box size. Rounds box size up and number of boxes down to get
     * a whole number of boxes along each axis.
     *
     * @param approxBoxLength Optional box dimensions parameter (Angstroms).
     */
    public void setApproxBoxLength(double approxBoxLength) {
        this.approxBoxLength = approxBoxLength;
    }

    /**
     * Sets behavior for how Residues are added to boxOptCells; 1 uses just
     * reference atom (C alpha for protein, N1/9 for nucleic acids), 2 uses any
     * atom, 3 uses any atom in any rotamer.
     *
     * @param boxInclusionCriterion Criterion to use
     */
    public void setBoxInclusionCriterion(int boxInclusionCriterion) {
        this.boxInclusionCriterion = boxInclusionCriterion;
    }

    public void setBoxStart(int boxStart) {
        this.boxStart = boxStart;
    }

    public void setBoxEnd(int boxEnd) {
        // Is -1 if boxes run to completion.
        this.boxEnd = boxEnd;
    }

    public double optimize(Algorithm algorithm) {
        this.algorithm = algorithm;
        return optimize();
    }

    public void setDirection(Direction direction) {
        this.direction = direction;
    }

    public void setDistanceCutoff(double distance) {
        this.distance = distance;
    }

    public void setIncrement(int increment) {
        this.increment = increment;
    }

    public void setRevert(boolean revert) {
        this.revert = revert;
    }

    /**
     * Sets the option to use a number of Monte Carlo steps for final
     * optimization.
     *
     * @param monteCarlo If Monte Carlo is to be considered
     * @param nMCsteps Number of steps to be taken
     */
    public void setMonteCarlo(boolean monteCarlo, int nMCsteps) {
        this.monteCarlo = monteCarlo;
        this.nMCsteps = nMCsteps;
    }

    public void setNucleicCorrectionThreshold(double nucleicCorrectionThreshold) {
        this.nucleicCorrectionThreshold = nucleicCorrectionThreshold;
    }

    public void setMinimumNumberAcceptedNARotamers(int minNumberAcceptedNARotamers) {
        this.minNumberAcceptedNARotamers = minNumberAcceptedNARotamers;
    }

    /**
     * Also sets derivative pruning factors.
     *
     * @param pruningFactor
     */
    public void setPruningFactor(double pruningFactor) {
        this.pruningFactor = pruningFactor;
        this.pairHalfPruningFactor = ((1.0 + pruningFactor) / 2);
    }

    public void setSingletonNAPruningFactor(double singletonNAPruningFactor) {
        this.singletonNAPruningFactor = singletonNAPruningFactor;
    }

    public void setVerboseEnergies(Boolean verboseEnergies) {
        this.verboseEnergies = verboseEnergies;
    }

    /**
     * Sets whether rotamer optimization should print out any files, or act solely
     * to optimize a structure in memory.
     * @param printFiles
     */
    public void setPrintFiles(boolean printFiles) {
        this.printFiles = printFiles;
    }

    /**
     * Sets the use of forced residues; an endForced value of -1 indicates not
     * to use forced residues.
     *
     * @param startForced First residue to force.
     * @param endForced Last residue to force.
     */
    public void setForcedResidues(int startForced, int endForced) {
        if (endForced != -1) {
            this.startForcedResidues = startForced;
            this.endForcedResidues = endForced;
            this.useForcedResidues = true;
        }
    }

    private double independent(List<Residue> residues) {
        if (x == null) {
            Atom atoms[] = molecularAssembly.getAtomArray();
            int nAtoms = atoms.length;
            x = new double[nAtoms * 3];
        }

        double e = Double.MAX_VALUE;
        List<Residue> rList = new ArrayList<>(Collections.nCopies(1, null));
        for (Residue residue : residues) {
            rList.set(0, residue);
            logger.info(String.format(" Optimizing %s side-chain.", residue));
            Rotamer[] rotamers = residue.getRotamers(library);
            potential.getCoordinates(x);
            e = Double.MAX_VALUE;
            int bestRotamer = -1;
            for (int j = 0; j < rotamers.length; j++) {
                Rotamer rotamer = rotamers[j];
                RotamerLibrary.applyRotamer(residue, rotamer);
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }
                double newE = currentEnergy(rList);
                if (newE < e) {
                    bestRotamer = j;
                }
            }
            if (bestRotamer > -1) {
                Rotamer rotamer = rotamers[bestRotamer];
                RotamerLibrary.applyRotamer(residue, rotamer);
            }
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            if (terminate) {
                logger.info(String.format("\n Terminating after residue %s.\n", residue));
                break;
            }
        }
        return e;
    }

    /**
     * Finds the first non-eliminated rotamer permutation.
     *
     * @param residues
     * @param i
     * @param currentRotamers
     *
     * @return If valid permutation found.
     */
    public boolean firstValidPerm(Residue residues[], int i, int currentRotamers[]) {
        // This is the initialization condition.
        int nResidues = residues.length;
        Residue residuei = residues[i];
        Rotamer[] rotamersi = residuei.getRotamers(library);
        int lenri = rotamersi.length;
        if (i < nResidues - 1) {
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                    /*for (int k = j + 1; k < i; k++) { // Check triples
                        int rk = currentRotamers[k];
                        deadEnd = check(j, rj, k, rk, i, ri);
                        if (deadEnd) {
                            break;
                        }
                    }
                    if (deadEnd) {
                        break;
                    }*/
                }
                if (deadEnd) {
                    continue;
                }
                currentRotamers[i] = ri;
                if (firstValidPerm(residues, i + 1, currentRotamers)) {
                    return true;
                }
                //dryRun(residues, i + 1, currentRotamers);
            }
        } else {
            /**
             * At the end of the recursion, check each rotamer of the final
             * residue.
             */
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                return true;
            }
        }
        return false;
    }

    /**
     * The main driver for optimizing a block of residues using DEE.
     *
     * @param residueList Residues to optimize.
     * @return Final energy.
     */
    private double globalUsingEliminations(List<Residue> residueList) {
        int currentEnsemble = Integer.MAX_VALUE;
        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        int nResidues = residues.length;
        int currentRotamers[] = new int[nResidues];
        int optimum[] = new int[nResidues];
        int iterations = 0;
        boolean finalTry = false;
        int bestEnsembleTargetDiffThusFar = Integer.MAX_VALUE;
        double bestBufferThusFar = ensembleBuffer;
        double startingBuffer = ensembleBuffer;

        if (ensembleEnergy > 0.0) {
            ensembleBuffer = ensembleEnergy;
            applyEliminationCriteria(residues, true, true);
            if (x == null) {
                Atom atoms[] = molecularAssembly.getAtomArray();
                int nAtoms = atoms.length;
                x = new double[nAtoms * 3];
            }
            /**
             * Compute the number of permutations without eliminating dead-ends
             * and compute the number of permutations using singleton
             * elimination.
             */
            double permutations = 1;
            double singletonPermutations = 1;
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                int nr = rotamers.length;
                if (nr > 1) {
                    int nrot = 0;
                    for (int ri = 0; ri < nr; ri++) {
                        if (!eliminatedSingles[i][ri]) {
                            nrot++;
                        }
                    }
                    permutations *= rotamers.length;
                    if (nrot > 1) {
                        singletonPermutations *= nrot;
                    }
                }
            }
            dryRun(residues, 0, currentRotamers);
            double pairTotalElimination = singletonPermutations - (double) evaluatedPermutations;
            if (evaluatedPermutations == 0) {
                logger.severe("No valid path through rotamer space found; try recomputing without pruning or using ensemble.");
            }
            if (master && printFiles && ensembleFile == null) {
                File file = molecularAssembly.getFile();
                String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
                //ensembleFile = SystemFilter.version(new File(filename + ".ens"));
                ensembleFile = new File(filename + ".ens");
                if (ensembleFile.exists()) {
                    for (int i = 2; i < 1000; i++) {
                        ensembleFile = new File(filename + ".ens_" + i);
                        if (!ensembleFile.exists()) {
                            break;
                        }
                    }
                    if (ensembleFile.exists()) {
                        logger.warning(String.format(" Versioning failed: appending to end of file %s", ensembleFile.getName()));
                    }
                }
                ensembleFilter = new PDBFilter(new File(ensembleFile.getName()), molecularAssembly, null, null);
                logger.info(String.format(" Ensemble file: %s", ensembleFile.getName()));
            }
            logIfMaster(String.format(" Number of permutations without DEE conditions: %10.4e.", permutations));
            logIfMaster(String.format(" Number of permutations after singleton eliminations: %10.4e.", singletonPermutations));
            logIfMaster(String.format(" Number of permutations removed by pairwise eliminations: %10.4e.", pairTotalElimination));
            logIfMaster(String.format(" Number of permutations remaining: %10.4e.", (double) evaluatedPermutations));

            double e;
            if (useMonteCarlo()) {
                firstValidPerm(residues, 0, currentRotamers);
                System.arraycopy(currentRotamers, 0, optimum, 0, nResidues);
                rotamerOptimizationMC(residues, optimum, currentRotamers, nMCsteps, false, mcUseAll);

                logIfMaster(" Ensembles not currently compatible with Monte Carlo search");
                /**
                 * Not currently compatible with ensembles.
                 */
            } else {
                double[] permutationEnergies = new double[evaluatedPermutations];
                ensembleStates = new ArrayList<>();

                e = rotamerOptimizationDEE(molecularAssembly, residues, 0, currentRotamers,
                        Double.MAX_VALUE, optimum, permutationEnergies);
                int[][] acceptedPermutations = new int[evaluatedPermutations][];
                for (int i = 0; i < acceptedPermutations.length; i++) {
                    acceptedPermutations[i] = null;
                }
                logIfMaster(String.format("\n Checking permutations for distance < %5.3f kcal/mol from GMEC energy %10.8f kcal/mol", ensembleEnergy, e));
                dryRunForEnsemble(residues, 0, currentRotamers, e, permutationEnergies, acceptedPermutations);
                int numAcceptedPermutations = 0;

                for (int i = 0; i < acceptedPermutations.length; i++) {
                    if (acceptedPermutations[i] != null) {
                        ++numAcceptedPermutations;
                        logIfMaster(String.format(" Accepting permutation %d at %8.6f < %8.6f", i, permutationEnergies[i] - e, ensembleEnergy));
                        for (int j = 0; j < nResidues; j++) {
                            Residue residuej = residues[j];
                            Rotamer[] rotamersj = residuej.getRotamers(library);
                            RotamerLibrary.applyRotamer(residuej, rotamersj[acceptedPermutations[i][j]]);
                        }

                        ResidueState[] states = ResidueState.storeAllCoordinates(residues);
                        ensembleStates.add(new ObjectPair<>(states, permutationEnergies[i]));

                        if (printFiles && master) {
                            try {
                                FileWriter fw = new FileWriter(ensembleFile, true);
                                BufferedWriter bw = new BufferedWriter(fw);
                                bw.write(String.format("MODEL        %d", numAcceptedPermutations));
                                for (int j = 0; j < 75; j++) {
                                    bw.write(" ");
                                }
                                bw.newLine();
                                bw.flush();
                                ensembleFilter.writeFile(ensembleFile, true);
                                bw.write(String.format("ENDMDL"));
                                for (int j = 0; j < 64; j++) {
                                    bw.write(" ");
                                }
                                bw.newLine();
                                bw.close();
                            } catch (IOException ex) {
                                logger.warning(String.format(" Exception writing to file: %s", ensembleFile.getName()));
                            }
                        }
                    }
                }
                logIfMaster(String.format(" Number of permutations within %5.3f kcal/mol of GMEC energy: %6.4e",
                        ensembleEnergy, (double) numAcceptedPermutations));
                ensembleStates.sort(null);
            }

            logIfMaster("\n Final rotamers:");
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                int ri = optimum[i];
                Rotamer rotamer = rotamers[ri];
                logIfMaster(String.format(" %c %s %s (%d)", residue.getChainID(), residue, rotamer.toString(), ri));
                RotamerLibrary.applyRotamer(residue, rotamer);
            }

            double sumSelfEnergy = 0;
            double sumPairEnergy = 0;
            double sumTrimerEnergy = 0;
            for (int i = 0; i < nResidues; i++) {
                int ri = optimum[i];
                sumSelfEnergy += selfEnergy[i][ri];
            }
            for (int i = 0; i < nResidues - 1; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues; j++) {
                    int rj = optimum[j];
                    sumPairEnergy += twoBodyEnergy[i][ri][j][rj];
                }
            }

            e = currentEnergy(residueList);
            logIfMaster(String.format(" Self Energy:   %16.8f", sumSelfEnergy));
            logIfMaster(String.format(" Pair Energy:   %16.8f", sumPairEnergy));

            double approximateEnergy = backboneEnergy + sumSelfEnergy + sumPairEnergy;

            if (threeBodyTerm) {
                for (int i = 0; i < nResidues - 2; i++) {
                    int ri = optimum[i];
                    for (int j = i + 1; j < nResidues - 1; j++) {
                        int rj = optimum[j];
                        for (int k = j + 1; k < nResidues; k++) {
                            int rk = optimum[k];
                            sumTrimerEnergy += threeBodyEnergy[i][ri][j][rj][k][rk];
                        }
                    }
                }
                approximateEnergy += sumTrimerEnergy;
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - backboneEnergy;
                logIfMaster(String.format(" Trimer Energy: %16.8f", sumTrimerEnergy));
                logIfMaster(String.format(" Neglected:     %16.8f", higherOrderEnergy));
            } else {
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - backboneEnergy;
                logIfMaster(String.format(" Neglected:     %16.8f", higherOrderEnergy));
            }
            logIfMaster(String.format(" Approximate Energy:     %16.8f", approximateEnergy));
            return e;
        }

        /**
         * Permutations used only to set maximum bound on ensembleNumber, thus
         * it is safe here to put that value in a 32-bit int.
         */
        int nPerms = 1;
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            int nr = rotamers.length;
            if (nr > 1) {
                nPerms *= rotamers.length;
            }
            if (nPerms > ensembleNumber) {
                break;
            }
        }

        if (nPerms < ensembleNumber) {
            logger.info(String.format(" Requested an ensemble of %d, but only %d permutations exist; returning full ensemble", ensembleNumber, nPerms));
            ensembleNumber = nPerms;
        }

        while (currentEnsemble != ensembleNumber) {
            if (monteCarlo) {
                logIfMaster(" Ensemble search not currently compatible with Monte Carlo");
                ensembleNumber = 1;
            }
            if (iterations == 0) {
                applyEliminationCriteria(residues, true, true);
            } else {
                applyEliminationCriteria(residues, false, false);
            }

            if (x == null) {
                Atom atoms[] = molecularAssembly.getAtomArray();
                int nAtoms = atoms.length;
                x = new double[nAtoms * 3];
            }

            /**
             * Compute the number of permutations without eliminating dead-ends
             * and compute the number of permutations using singleton
             * elimination.
             */
            double permutations = 1;
            double singletonPermutations = 1;
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                int nr = rotamers.length;
                if (nr > 1) {
                    int nrot = 0;
                    for (int ri = 0; ri < nr; ri++) {
                        if (!eliminatedSingles[i][ri]) {
                            nrot++;
                        }
                    }
                    permutations *= rotamers.length;
                    if (nrot > 1) {
                        singletonPermutations *= nrot;
                    }
                }
            }
            dryRun(residues, 0, currentRotamers);
            double pairTotalElimination = singletonPermutations - (double) evaluatedPermutations;
            currentEnsemble = (int) evaluatedPermutations;
            if (ensembleNumber == 1 && currentEnsemble == 0) {
                logger.severe("No valid path through rotamer space found; try recomputing without pruning or using ensemble.");
                /*  PROGRAMMATIC ENSEMBLE CONTROL (dangerous)
                 if (getPruning() == 0) {
                 logger.warning(" Unable to recover a rotamer path; commencing ensemble search.");
                 setEnsemble(10, 5.0);
                 startingBuffer = 5.0;
                 continue;
                 } else {
                 logger.warning(" Pruning left no valid path through rotamer space; recomputing without pruning.");
                 setPruning(0);
                 continue;
                 }
                 */
            }
            if (ensembleNumber > 1) {
                if (master && printFiles && ensembleFile == null) {
                    File file = molecularAssembly.getFile();
                    String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
                    ensembleFile = new File(filename + ".ens");
                    if (ensembleFile.exists()) {
                        for (int i = 2; i < 1000; i++) {
                            ensembleFile = new File(filename + ".ens_" + i);
                            if (!ensembleFile.exists()) {
                                break;
                            }
                        }
                        if (ensembleFile.exists()) {
                            logger.warning(String.format(" Versioning failed: appending to end of file %s", ensembleFile.getName()));
                        }
                    }
                    ensembleFilter = new PDBFilter(new File(ensembleFile.getName()), molecularAssembly, null, null);
                    logger.info(String.format(" Ensemble file: %s", ensembleFile.getName()));
                }
                logIfMaster(String.format(" Ensemble Search Stats: (buffer: %5.3f, current: %d, target: %d)", ensembleBuffer, currentEnsemble, ensembleNumber));
            }
            if (ensembleNumber == 1 || finalTry) {
                logIfMaster(String.format(" Number of permutations without DEE conditions: %10.4e.", permutations));
                logIfMaster(String.format(" Number of permutations after singleton eliminations: %10.4e.", singletonPermutations));
                logIfMaster(String.format(" Number of permutations removed by pairwise eliminations: %10.4e.", pairTotalElimination));
                logIfMaster(String.format(" Number of permutations remaining: %10.4e.", (double) evaluatedPermutations));
                break;
            }
            if (Math.abs(currentEnsemble - ensembleNumber) < bestEnsembleTargetDiffThusFar) {
                bestEnsembleTargetDiffThusFar = Math.abs(currentEnsemble - ensembleNumber);
                bestBufferThusFar = ensembleBuffer;
            }
            if (currentEnsemble > ensembleNumber) {
                ensembleBuffer -= ensembleBufferStep;
                ensembleBufferStep -= (ensembleBufferStep * 0.01);
                iterations++;
            } else if (currentEnsemble < ensembleNumber) {
                ensembleBuffer += ensembleBufferStep;
                ensembleBufferStep -= (ensembleBufferStep * 0.01);
                iterations++;
            }
            if (iterations > 100) {
                if (currentEnsemble == 0) {
                    // TODO: Decide whether we like these next four lines.  Has the potential to produce a crazy amount of permutations.
                    logIfMaster(" Ensemble still empty; increasing buffer energy.");
                    startingBuffer = 3 * startingBuffer;
                    setEnsemble(10, startingBuffer);
                    iterations = 0;
                } else {
                    ensembleBuffer = bestBufferThusFar;
                    finalTry = true;
                }
            }
        }

        if (currentEnsemble == 0) {
            logger.warning(" No valid rotamer permutations found; results will be unreliable.  Try increasing the starting ensemble buffer.");
        }
        double[] permutationEnergyStub = null;
        if (useMonteCarlo()) {
            firstValidPerm(residues, 0, currentRotamers);
            rotamerOptimizationMC(residues, optimum, currentRotamers, nMCsteps, false, mcUseAll);
        } else {
            rotamerOptimizationDEE(molecularAssembly, residues, 0, currentRotamers,
                    Double.MAX_VALUE, optimum, permutationEnergyStub);
        }

        logIfMaster("\n Final rotamers:");
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            int ri = optimum[i];
            Rotamer rotamer = rotamers[ri];
            logIfMaster(String.format(" %c %s %s (%d)", residue.getChainID(), residue, rotamer.toString(), ri));
            RotamerLibrary.applyRotamer(residue, rotamer);
        }

        double sumSelfEnergy = 0;
        double sumPairEnergy = 0;
        double sumTrimerEnergy = 0;
        for (int i = 0; i < nResidues; i++) {
            int ri = optimum[i];
            sumSelfEnergy += selfEnergy[i][ri];
        }
        for (int i = 0; i < nResidues - 1; i++) {
            int ri = optimum[i];
            for (int j = i + 1; j < nResidues; j++) {
                int rj = optimum[j];
                sumPairEnergy += twoBodyEnergy[i][ri][j][rj];
            }
        }

        double e = currentEnergy(residueList);
        logIfMaster(String.format(" Self Energy:   %16.8f", sumSelfEnergy));
        logIfMaster(String.format(" Pair Energy:   %16.8f", sumPairEnergy));

        double approximateEnergy = backboneEnergy + sumSelfEnergy + sumPairEnergy;

        if (threeBodyTerm) {
            for (int i = 0; i < nResidues - 2; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues - 1; j++) {
                    int rj = optimum[j];
                    for (int k = j + 1; k < nResidues; k++) {
                        int rk = optimum[k];
                        sumTrimerEnergy += threeBodyEnergy[i][ri][j][rj][k][rk];
                    }
                }
            }
            approximateEnergy += sumTrimerEnergy;
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - backboneEnergy;
            logIfMaster(String.format(" Trimer Energy: %16.8f", sumTrimerEnergy));
            logIfMaster(String.format(" Neglected:     %16.8f", higherOrderEnergy));
        } else {
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - backboneEnergy;
            logIfMaster(String.format(" Neglected:     %16.8f", higherOrderEnergy));
        }

        logIfMaster(String.format(" Approximate Energy:     %16.8f", approximateEnergy));
        return e;
    }

    /**
     * Use Monte Carlo if monteCarlo specified, and either skipDEE specified or
     * nMCsteps is smaller then the remaining permutation size.
     *
     * @return Finish DEE search with Monte Carlo.
     */
    private boolean useMonteCarlo() {
        return monteCarlo && (mcNoEnum || (nMCsteps < evaluatedPermutations));
    }

    /**
     * Performs a recursive brute-force rotamer optimization over a passed list
     * of residues. Should make a final decision on how to tell the method to
     * use true brute-force or a 2/3-body energy decomposition. Temporarily uses
     * the decomposition if parallelEnergies is true.
     *
     * @param residueList Residues to be optimized.
     * @return GMEC (Global Minimum Energy Conformation) energy.
     */
    private double globalBruteForce(List<Residue> residueList) {

        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        int nResidues = residues.length;
        int optimum[] = new int[nResidues];

        if (x == null) {
            Atom atoms[] = molecularAssembly.getAtomArray();
            int nAtoms = atoms.length;
            x = new double[nAtoms * 3];
        }

        /**
         * Compute the number of permutations without eliminating dead-ends and
         * compute the number of permutations using singleton elimination.
         */
        double permutations = 1;
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            permutations *= rotamers.length;
        }

        logger.info(String.format(" Number of permutations: %16.8e.", permutations));
        double e;
        // Temporary method to make it use the 2/3-body computation.
        if (parallelEnergies) {
            useFullAMOEBAEnergy = false;
        }
        if (!useFullAMOEBAEnergy) {
            setPruning(0);
            rotamerEnergies(residues);
            int[] rotamerSet = new int[nResidues];
            fill(rotamerSet, 0);
            e = decomposedRotamerOptimization(molecularAssembly, residues, 0, Double.MAX_VALUE, optimum, rotamerSet);
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                RotamerLibrary.applyRotamer(residue, rotamers[optimum[i]]);
                turnOnAtoms(residue);
            }
            double fullEnergy = currentEnergy(residueList);
            logger.info(String.format(" Final summation of energies:    %16.5f", e));
            logger.info(String.format(" Final energy of optimized structure:    %16.5f", fullEnergy));
            logger.info(String.format(" Neglected:    %16.5f", fullEnergy - e));
        } else {
            e = rotamerOptimization(molecularAssembly, residues, 0, Double.MAX_VALUE, optimum);
        }

        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            int ri = optimum[i];
            Rotamer rotamer = rotamers[ri];
            logger.info(String.format(" %s %s (%d)", residue.getResidueNumber(), rotamer.toString(), ri));
            RotamerLibrary.applyRotamer(residue, rotamer);
            if (useFullAMOEBAEnergy) {
                e = currentEnergy(residueList);
            }
        }
        return e;
    }

    /**
     * Checks the distance matrix, finding the shortest distance between two
     * residues' rotamers under any symmetry operator; will evaluate this if
     * distance matrix not already filled.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @return Shortest distance
     */
    private double checkDistanceMatrix(int i, int ri, int j, int rj) {
        if (i > j) {
            double dist = distanceMatrix[j][rj][i][ri];
            if (dist < 0) {
                dist = evaluateDistance(j, rj, i, ri);
                distanceMatrix[j][rj][i][ri] = dist;
            }
            return dist;
        }
        double dist = distanceMatrix[i][ri][j][rj];
        if (dist < 0) {
            dist = evaluateDistance(i, ri, j, rj);
            distanceMatrix[i][ri][j][rj] = dist;
        }
        return dist;
    }

    /**
     * Returns the mean of three dimer distances.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @param k Residue k
     * @param rk Rotamer for k
     * @return mean separation distance
     */
    private double trimerDistance(int i, int ri, int j, int rj, int k, int rk) {
        double ij = checkDistanceMatrix(i, ri, j, rj);
        double ik = checkDistanceMatrix(i, ri, k, rk);
        double jk = checkDistanceMatrix(j, rj, k, rk);
        return (ij + ik + jk) / 3.0;
    }

    /**
     * Returns the mean of six dimer distances.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @param k Residue k
     * @param rk Rotamer for k
     * @param l Residue l
     * @param rl Rotamer for l
     * @return mean separation distance
     */
    private double quadDistance(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        double ij = checkDistanceMatrix(i, ri, j, rj);
        double ik = checkDistanceMatrix(i, ri, k, rk);
        double il = checkDistanceMatrix(i, ri, l, rl);
        double jk = checkDistanceMatrix(j, rj, k, rk);
        double jl = checkDistanceMatrix(j, rj, l, rl);
        double kl = checkDistanceMatrix(k, rk, l, rl);
        return (ij + ik + il + jk + jl + kl) / 6.0;
    }

    /**
     * Evaluates the pairwise distance between two residues' rotamers under any
     * symmetry operator; does "lazy loading" for the distance matrix.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @return Shortest distance
     */
    private double evaluateDistance(int i, int ri, int j, int rj) {
        Residue resi = allResiduesArray[i];
        Rotamer[] rotamersI = resi.getRotamers(library);
        Rotamer roti = rotamersI[ri];
        double[][] xi;
        if (roti.equals(resi.getRotamer())) {
            xi = resi.storeCoordinateArray();
        } else {
            ResidueState origI = resi.storeState();
            RotamerLibrary.applyRotamer(resi, roti);
            xi = resi.storeCoordinateArray();
            resi.revertState(origI);
        }

        Residue resj = allResiduesArray[j];
        Rotamer[] rotamersJ = resj.getRotamers(library);
        Rotamer rotj = rotamersJ[rj];
        double[][] xj;
        if (rotj.equals(resj.getRotamer())) {
            xj = resj.storeCoordinateArray();
        } else {
            ResidueState origJ = resj.storeState();
            RotamerLibrary.applyRotamer(resj, rotj);
            xj = resj.storeCoordinateArray();
            resj.revertState(origJ);
        }

        Crystal crystal = molecularAssembly.getCrystal();
        int nSymm = crystal.spaceGroup.getNumberOfSymOps();
        double minDist = Double.MAX_VALUE;
        for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
            SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
            double dist = interResidueDistance(xi, xj, symOp);
            minDist = dist < minDist ? dist : minDist;
        }
        return minDist;
    }

    /**
     * Checks whether an index is inside the forced residue range.
     *
     * @param i Index to check
     * @return If forced
     * @throws IllegalStateException If useForcedResidues not true.
     */
    private boolean checkIfForced(int i) throws IllegalStateException {
        if (!useForcedResidues) {
            throw new IllegalStateException("checkForcedResidues being called without useForcedResidues.");
        }
        return (startForcedResidues <= i && i <= endForcedResidues);
    }

    /**
     * A more stringent form of checkIfForced which returns true only if the
     * residue is forced, false only if it has rotamers, and otherwise throws a
     * null pointer exception.
     *
     * @param residue Residue to check
     * @return True if forced, false if it has rotamers.
     * @throws NullPointerException If non-rotameric and non-forced.
     */
    private boolean checkIfForced(Residue residue) throws NullPointerException {
        if (residue.getRotamers(library) != null) {
            return false;
        } else if (useForcedResidues && checkIfForced(residue.getResidueNumber())) {
            return true;
        } else {
            throw new NullPointerException(String.format(" Non-rotameric, non-forced residue present "
                    + "in residue list: %c %s-%d", residue.getChainID(), residue.toString(), residue.getResidueNumber()));
        }
    }

    private double slidingWindow(ArrayList<Residue> residueList, int windowSize, int increment, boolean revert,
            double distance, Direction direction) {

        long beginTime = -System.nanoTime();
        boolean incrementTruncated = false;
        boolean firstWindowSaved = false;
        int counter = 1;
        int windowEnd;

        int nOptimize = residueList.size();
        if (nOptimize < windowSize) {
            windowSize = nOptimize;
            logger.warning(" Window size too small for given residue range; truncating window size.");
        }
        switch (direction) {
            case BACKWARD:
                ArrayList<Residue> temp = new ArrayList<>();
                for (int i = nOptimize - 1; i >= 0; i--) {
                    temp.add(residueList.get(i));
                }
                residueList = temp;
            // Fall through into the FORWARD case.
            case FORWARD:
                for (int windowStart = 0; windowStart + (windowSize - 1) < nOptimize; windowStart += increment) {
                    long windowTime = -System.nanoTime();
                    windowEnd = windowStart + (windowSize - 1);
                    logIfMaster(String.format("\n Iteration %d of the sliding window.\n", counter++));
                    Residue firstResidue = residueList.get(windowStart);
                    Residue lastResidue = residueList.get(windowEnd);
                    if (firstResidue != lastResidue) {
                        logIfMaster(String.format(" Residues %s ... %s", firstResidue.toString(), lastResidue.toString()));
                    } else {
                        logIfMaster(String.format(" Residue %s", firstResidue.toString()));
                    }
                    ArrayList<Residue> currentWindow = new ArrayList<>();
                    ArrayList<Residue> onlyRotameric = new ArrayList<>(); // Not filled if useForcedResidues == false.
                    for (int i = windowStart; i <= windowEnd; i++) {
                        Residue residue = residueList.get(i);
                        if (useForcedResidues && residue.getRotamers(library) != null) {
                            onlyRotameric.add(residue);
                        }
                        currentWindow.add(residueList.get(i));
                    }

                    if (distance > 0) {
                        for (int i = windowStart; i <= windowEnd; i++) {
                            Residue residuei = residueList.get(i);
                            int indexI = allResiduesList.indexOf(residuei);
                            int lengthRi;
                            if (checkIfForced(residuei)) {
                                lengthRi = 1;
                            } else {
                                lengthRi = residuei.getRotamers(library).length;
                            }
                            for (int ri = 0; ri < lengthRi; ri++) {
                                for (int j = 0; j < numResidues; j++) {
                                    Residue residuej = allResiduesArray[j];
                                    Rotamer rotamersj[] = residuej.getRotamers(library);
                                    if (currentWindow.contains(residuej) || rotamersj == null) {
                                        continue;
                                    }
                                    int lengthRj = rotamersj.length;
                                    for (int rj = 0; rj < lengthRj; rj++) {
                                        double rotamerSeparation = checkDistanceMatrix(indexI, ri, j, rj);
                                        // if (distanceMatrix[indexI][ri][j][rj] <= distance) {
                                        if (rotamerSeparation <= distance) {
                                            if (!currentWindow.contains(residuej)) {
                                                logIfMaster(String.format(" Adding residue %s at distance %16.8f Ang from %s %d.",
                                                        residuej, rotamerSeparation, residuei, ri));
                                                currentWindow.add(residuej);
                                                if (useForcedResidues) {
                                                    onlyRotameric.add(residuej);
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    /**
                     * If the window starts with a nucleic acid, and there is a
                     * 5' NA residue, ensure that 5' NA residue has been
                     * included in the window. Otherwise, that previous residue
                     * may not have had a chance to be flexible about its sugar
                     * pucker.
                     *
                     * If window size is greater than increment, however, this
                     * has already been handled. Additionally, do not perform
                     * this for the first window (counter is already incremented
                     * by the time this check is performed, so first window's
                     * counter will be 2). Furthermore, do not include Residues
                     * with null Rotamer lists (this breaks things).
                     *
                     * The issue: if window size = increment, the last NA
                     * residue in each window will not have flexibility about
                     * its sugar pucker, because its self-energy includes the
                     * O3' (i) to P (i+1) bond, so it must remain in the
                     * original sugar pucker to meet the i+1 residue. However,
                     * this problem can be solved by ensuring that final residue
                     * is included in the next window, where it will have
                     * flexibility about its sugar pucker.
                     *
                     * If you are running successive sliding window jobs on the
                     * same file, I would suggest starting the next job on the
                     * last residue of the previous job, unless you know your
                     * settings will include it.
                     */
                    if (counter > 2 && windowSize <= increment && firstResidue.getResidueType() == NA) {
                        Residue prevResidue = firstResidue.getPreviousResidue();
                        if (prevResidue != null && prevResidue.getResidueType() == NA && !currentWindow.contains(prevResidue) && prevResidue.getRotamers(library) != null) {
                            logIfMaster(String.format(" Adding nucleic acid residue 5' of window start %s to give it flexibility about its sugar pucker.",
                                    prevResidue.toString()));
                            currentWindow.add(prevResidue);
                            if (useForcedResidues) {
                                onlyRotameric.add(prevResidue);
                            }
                        }
                    }
                    if (useForcedResidues) {
                        sortResidues(onlyRotameric);
                    } else {
                        sortResidues(currentWindow);
                    }

                    if (revert) {
                        ResidueState[] coordinates = ResidueState.storeAllCoordinates(currentWindow);
                        // x is an array of coordinates for the entire molecular assembly.
                        // If x has not yet been constructed, construct it.
                        if (x == null) {
                            Atom atoms[] = molecularAssembly.getAtomArray();
                            int nAtoms = atoms.length;
                            x = new double[nAtoms * 3];
                        }
                        double startingEnergy = currentEnergy(currentWindow);
                        if (useForcedResidues) {
                            if (onlyRotameric.size() < 1) {
                                logger.info(" Window has no rotameric residues.");
                                ResidueState.revertAllCoordinates(currentWindow, coordinates);
                            } else {
                                globalUsingEliminations(onlyRotameric);
                                double finalEnergy = currentEnergy(currentWindow);
                                if (startingEnergy <= finalEnergy) {
                                    logger.warning("Optimization did not yield a better energy. Reverting to orginal coordinates.");
                                    ResidueState.revertAllCoordinates(currentWindow, coordinates);
                                }
                            }
                        } else {
                            globalUsingEliminations(currentWindow);
                            double finalEnergy = currentEnergy(currentWindow);
                            if (startingEnergy <= finalEnergy) {
                                logger.warning("Optimization did not yield a better energy. Reverting to orginal coordinates.");
                                ResidueState.revertAllCoordinates(currentWindow, coordinates);
                            }
                        }
                    } else if (useForcedResidues) {
                        if (onlyRotameric.size() < 1) {
                            logger.info(" Window has no rotameric residues.");
                        } else {
                            globalUsingEliminations(onlyRotameric);
                        }
                    } else {
                        globalUsingEliminations(currentWindow);
                    }
                    if (!incrementTruncated) {
                        if (windowStart + (windowSize - 1) + increment > nOptimize - 1) {
                            increment = nOptimize - windowStart - windowSize;
                            if (increment == 0) {
                                break;
                            }
                            logger.warning(" Increment truncated in order to optimize entire residue range.");
                            incrementTruncated = true;
                        }
                    }

                    if (master && printFiles) {
                        File file = molecularAssembly.getFile();
                        if (firstWindowSaved) {
                            file.delete();
                        }
                        // Don't write a file if its the final iteration.
                        if (windowStart + windowSize == nOptimize) {
                            continue;
                        }
                        //   String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
                        //   File tempFile = new File(filename + ".win");
                        //   PDBFilter windowFilter = new PDBFilter(new File(filename + ".win"), molecularAssembly, null, null);
                        PDBFilter windowFilter = new PDBFilter(file, molecularAssembly, null, null);
                        //   StringBuilder header = new StringBuilder(String.format("Iteration %d of the sliding window\n", counter - 1));
                        try {
                            windowFilter.writeFile(file, false);
                            if (firstResidue != lastResidue) {
                                logger.info(String.format(" File with residues %s ... %s in window written to.", firstResidue.toString(), lastResidue.toString()));
                            } else {
                                logger.info(String.format(" File with residue %s in window written to.", firstResidue.toString()));
                            }
                        } catch (Exception e) {
                            logger.warning(String.format("Exception writing to file: %s", file.getName()));
                        }
                        firstWindowSaved = true;
                    }
                    long currentTime = System.nanoTime();
                    windowTime += currentTime;
                    logIfMaster(String.format(" Time elapsed for this iteration: %11.3f sec", windowTime * 1.0E-9));
                    logIfMaster(String.format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
                    /*for (Residue residue : residueList) {
                        if (residue instanceof MultiResidue) {
                            ((MultiResidue) residue).setDefaultResidue();
                            residue.reInitOriginalAtomList();
                        }
                    }*/
                }
                break;

            default:
                // No default case.
                break;
        }
        return 0.0;
    }

    /**
     * Returns the superbox used to generate the boxes for sliding box. If
     * superbox coordinates manually set, uses them plus the defined buffer.
     * Else, if an aperiodic system, uses maximum and minimum C alpha (or N1/9)
     * coordinates plus superboxBuffer (by default, 8A, longer than a lysine
     * side chain or N1/N9 distance to any other atom). Else, it just uses the
     * ordinary crystal.
     *
     * @param residueList List of residues to incorporate.
     * @return Superbox crystal.
     */
    private Crystal generateSuperbox(List<Residue> residueList) {
        double[] maxXYZ = new double[3];
        double[] minXYZ = new double[3];
        Crystal originalCrystal = molecularAssembly.getCrystal();
        if (manualSuperbox) {
            for (int i = 0; i < maxXYZ.length; i++) {
                int ii = 2 * i;
                minXYZ[i] = boxDimensions[ii] - superboxBuffer;
                maxXYZ[i] = boxDimensions[ii + 1] + superboxBuffer;
            }
        } else if (originalCrystal.aperiodic()) {
            if (residueList == null || residueList.isEmpty()) {
                throw new IllegalArgumentException(" Null or empty residue list when generating superbox.");
            }
            Atom initializerAtom = residueList.get(0).getReferenceAtom();
            initializerAtom.getXYZ(minXYZ);
            initializerAtom.getXYZ(maxXYZ);
            for (Residue residue : residueList) {
                Atom refAtom = residue.getReferenceAtom();
                double[] refAtomCoords = new double[3];
                refAtom.getXYZ(refAtomCoords);
                for (int i = 0; i < 3; i++) {
                    maxXYZ[i] = (refAtomCoords[i] > maxXYZ[i] ? refAtomCoords[i] : maxXYZ[i]);
                    minXYZ[i] = (refAtomCoords[i] < minXYZ[i] ? refAtomCoords[i] : minXYZ[i]);
                }
            }
            for (int i = 0; i < 3; i++) {
                minXYZ[i] -= superboxBuffer;
                maxXYZ[i] += superboxBuffer;
            }
        } else {
            return originalCrystal;
        }
        double newA = maxXYZ[0] - minXYZ[0];
        double newB = maxXYZ[1] - minXYZ[1];
        double newC = maxXYZ[2] - minXYZ[2];
        if (manualSuperbox) {
            logger.info(String.format(" Manual superbox set over (minX, maxX, minY, "
                    + "maxY, minZ, maxZ): %f, %f, %f, %f, %f, %f",
                    minXYZ[0], maxXYZ[0], minXYZ[1], maxXYZ[1], minXYZ[2], maxXYZ[2]));
            logger.info(String.format(" Buffer size (included in dimensions): %f\n", superboxBuffer));
        } else { // Crystal systems will have already returned.
            logger.info(" System is aperiodic: protein box generated over these coordinates (minX, maxX, minY, maxY, minZ, maxZ):");
            String message = " Aperiodic box dimensions: ";
            for (int i = 0; i < minXYZ.length; i++) {
                message = message.concat(String.format("%f,%f,", minXYZ[i], maxXYZ[i]));
            }
            message = message.substring(0, message.length() - 1);
            logger.info(message);
            logger.info(String.format(" Buffer size (included in dimensions): %f\n", superboxBuffer));
        }
        return new Crystal(newA, newB, newC, 90.0, 90.0, 90.0, "P1");
    }

    /**
     * Breaks down a structure into a number of overlapping boxes for
     * optimization.
     *
     * @return Potential energy of final structure.
     */
    private double boxOptimization(ArrayList<Residue> residueList) {
        this.usingBoxOptimization = true;
        long beginTime = -System.nanoTime();
        Residue residues[] = residueList.toArray(new Residue[residueList.size()]);
        boolean firstCellSaved = false;

        /*
         * A new dummy Crystal will be constructed for an aperiodic system. The
         * purpose is to avoid using the overly large dummy Crystal used for
         * Ewald purposes. Atoms are not and should not be moved into the dummy
         * Crystal boundaries; to check if an Atom is inside a cell, use an
         * array of coordinates adjusted to be 0 < coordinate < 1.0.
         */
        Crystal crystal = generateSuperbox(residueList);

        // Cells indexed by x*(YZ divisions) + y*(Z divisions) + z.
        int totalCells = getTotalCellCount(crystal); // Also initializes cell count if using -bB
        if (boxStart > totalCells - 1) {
            logger.severe(String.format(" FFX shutting down: Box optimization start is out of range of total boxes: %d > %d", (boxStart + 1), totalCells));
        }
        if (boxEnd > totalCells - 1) {
            boxEnd = totalCells - 1;
            logIfMaster(" Final box out of range: reset to last possible box.");
        } else if (boxEnd < 0) {
            boxEnd = totalCells - 1;
        }
        BoxOptCell[] cells = loadCells(crystal, residues);
        int numCells = cells.length;
        logIfMaster(String.format(" Optimizing boxes  %d  to  %d", (boxStart + 1), (boxEnd + 1)));
        /*int restartCell = -1;
         if (parallelEnergies && loadEnergyRestart) {
         restartCell = loadEnergyRestartIterations(energyRestartFile);
         if (restartCell > -1) {
         logIfMaster(String.format(" Optimization restarting from file at cell #%d", restartCell));
         }
         }*/
        for (int i = 0; i < numCells; i++) {
            /*if (restartCell > -1 && i < restartCell) {
             continue;
             }*/
            BoxOptCell celli = cells[i];
            ArrayList<Residue> residuesList = celli.getResiduesAsList();
            int[] cellIndices = celli.getXYZIndex();
            logIfMaster(String.format("\n Iteration %d of the box optimization.", (i + 1)));
            logIfMaster(String.format(" Cell index (linear): %d", (celli.getLinearIndex() + 1)));
            logIfMaster(String.format(" Cell xyz indices: x = %d, y = %d, z = %d", cellIndices[0] + 1, cellIndices[1] + 1, cellIndices[2] + 1));
            int nResidues = residuesList.size();
            if (nResidues > 0) {
                // SDL additions for writing/loading box-based restart files.
                readyForSingles = false;
                selfsDone = false;
                readyForPairs = false;
                pairsDone = false;
                readyForTrimers = false;
                trimersDone = false;
                energiesToWrite = Collections.synchronizedList(new ArrayList<String>());
                receiveThread = new ReceiveThread(residuesList.toArray(new Residue[1]));
                receiveThread.start();
                if (master && writeEnergyRestart && printFiles) {
                    if (energyWriterThread != null) {
                        int waiting = 0;
                        while (energyWriterThread.getState() != java.lang.Thread.State.TERMINATED) {
                            try {
                                if (waiting++ > 20) {
                                    logger.log(Level.SEVERE, " ReceiveThread/EnergyWriterThread from previous box locked up.");
                                }
                                logIfMaster(" Waiting for previous iteration's communication threads to shut down... ");
                                Thread.sleep(10000);
                            } catch (InterruptedException ex) {
                            }
                        }
                    }
                    energyWriterThread = new EnergyWriterThread(receiveThread, i + 1, cellIndices);
                    energyWriterThread.start();
                }

                if (loadEnergyRestart) {
                    boxLoadIndex = i + 1;
                    boxLoadCellIndices = new int[3];
                    boxLoadCellIndices[0] = cellIndices[0];
                    boxLoadCellIndices[1] = cellIndices[1];
                    boxLoadCellIndices[2] = cellIndices[2];
                }

                long boxTime = -System.nanoTime();
                Residue firstResidue = residuesList.get(0);
                Residue lastResidue = residuesList.get(nResidues - 1);
                if (firstResidue != lastResidue) {
                    logIfMaster(String.format(" Residues %s ... %s", firstResidue.toString(), lastResidue.toString()));
                } else {
                    logIfMaster(String.format(" Residue %s", firstResidue.toString()));
                }
                if (revert) {
                    ResidueState[] coordinates = ResidueState.storeAllCoordinates(residuesList);
                    // x is an array of coordinates for the entire molecular assembly.
                    // If x has not yet been constructed, construct it.
                    if (x == null) {
                        Atom atoms[] = molecularAssembly.getAtomArray();
                        int nAtoms = atoms.length;
                        x = new double[nAtoms * 3];
                    }
                    double startingEnergy = currentEnergy(residuesList);
                    globalUsingEliminations(residuesList);
                    double finalEnergy = currentEnergy(residuesList);
                    if (startingEnergy <= finalEnergy) {
                        logger.warning("Optimization did not yield a better energy. Reverting to orginal coordinates.");
                        ResidueState.revertAllCoordinates(residuesList, coordinates);
                    }
                    long currentTime = System.nanoTime();
                    boxTime += currentTime;
                    logIfMaster(String.format(" Time elapsed for this iteration: %11.3f sec", boxTime * 1.0E-9));
                    logIfMaster(String.format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
                } else {
                    globalUsingEliminations(residuesList);
                    long currentTime = System.nanoTime();
                    boxTime += currentTime;
                    logIfMaster(String.format(" Time elapsed for this iteration: %11.3f sec", boxTime * 1.0E-9));
                    logIfMaster(String.format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
                }
                if (master && printFiles) {
                    String filename = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath()) + ".partial";
                    File file = new File(filename);
                    if (firstCellSaved) {
                        file.delete();
                    }
                    // Don't write a file if it's the final iteration.
                    if (i == (numCells - 1)) {
                        continue;
                    }
                    PDBFilter windowFilter = new PDBFilter(file, molecularAssembly, null, null);
                    try {
                        windowFilter.writeFile(file, false);
                        if (firstResidue != lastResidue) {
                            logIfMaster(String.format(" File with residues %s ... %s in window written.", firstResidue.toString(), lastResidue.toString()));
                        } else {
                            logIfMaster(String.format(" File with residue %s in window written.", firstResidue.toString()));
                        }
                        firstCellSaved = true;
                    } catch (Exception e) {
                        logger.warning(String.format("Exception writing to file: %s", file.getName()));
                    }
                }
                /*for (Residue residue : residueList) {
                    if (residue instanceof MultiResidue) {
                        ((MultiResidue) residue).setDefaultResidue();
                        residue.reInitOriginalAtomList();
                    }
                }*/
            } else {
                logIfMaster(String.format(" Empty box: no residues found."));
            }
        }
        return 0.0;
    }

    /**
     * Returns the number of cells (boxes) for box optimization; if the -bB flag
     * is set, sets the final number of cells.
     *
     * @param crystal Crystal or dummy crystal being used to define boxes.
     * @return Total number of cells.
     */
    private int getTotalCellCount(Crystal crystal) {
        int numCells = 1;
        if (approxBoxLength > 0) {
            double[] boxes = new double[3];
            boxes[0] = crystal.a / approxBoxLength;
            boxes[1] = crystal.b / approxBoxLength;
            boxes[2] = crystal.c / approxBoxLength;
            for (int i = 0; i < boxes.length; i++) {
                if (boxes[i] < 1) {
                    numXYZBoxes[i] = 1;
                } else {
                    numXYZBoxes[i] = (int) boxes[i];
                }
            }
        }
        for (int i = 0; i < numXYZBoxes.length; i++) {
            numCells *= numXYZBoxes[i];
        }
        return numCells;
    }

    /**
     * Creates and fills cells (boxes) for box optimization.
     *
     * @param crystal Polymer crystal or dummy crystal
     * @param residues All residues to be optimized
     * @return Filled cells.
     */
    private BoxOptCell[] loadCells(Crystal crystal, Residue[] residues) {
        double xCellBorderFracSize = (boxBorderSize / crystal.a);
        double yCellBorderFracSize = (boxBorderSize / crystal.b);
        double zCellBorderFracSize = (boxBorderSize / crystal.c);
        logIfMaster(String.format(" Number of boxes along x: %d, y: %d, z: %d", numXYZBoxes[0], numXYZBoxes[1], numXYZBoxes[2]));

        int numCells = boxEnd - boxStart + 1;
        BoxOptCell[] cells = new BoxOptCell[numCells];
        int currentIndex = 0;
        int filledCells = 0;
        int[] xyzIndices = new int[3];
        boolean doBreak = false; // Breaks the ijk loop if the last box passed.
        /*
         * Initializes coordinates for all the boxes, indexed linearly along z,
         * then y, then x (so the box with xyz indices 2,3,2 in a crystal with
         * 4, 5, and 3 boxes along xyz would be indexed 2*5*3 + 3*3 + 2 = 41.
         * The int[] indices stores seperate x, y, and z indices.
         */
        for (int i = 0; i < numXYZBoxes[0]; i++) {
            if (doBreak) {
                break;
            }
            xyzIndices[0] = i;
            for (int j = 0; j < numXYZBoxes[1]; j++) {
                if (doBreak) {
                    break;
                }
                xyzIndices[1] = j;
                for (int k = 0; k < numXYZBoxes[2]; k++) {
                    if (currentIndex < boxStart) {
                        ++currentIndex;
                        continue;
                    } else if (currentIndex > boxEnd) {
                        doBreak = true;
                        break;
                    }
                    xyzIndices[2] = k;
                    double[] fracCoords = new double[6];
                    fracCoords[0] = (((1.0 * i) / numXYZBoxes[0]) - xCellBorderFracSize);
                    fracCoords[1] = (((1.0 * j) / numXYZBoxes[1]) - yCellBorderFracSize);
                    fracCoords[2] = (((1.0 * k) / numXYZBoxes[2]) - zCellBorderFracSize);
                    fracCoords[3] = (((1.0 + i) / numXYZBoxes[0]) + xCellBorderFracSize);
                    fracCoords[4] = (((1.0 + j) / numXYZBoxes[1]) + yCellBorderFracSize);
                    fracCoords[5] = (((1.0 + k) / numXYZBoxes[2]) + zCellBorderFracSize);
                    cells[filledCells++] = new BoxOptCell(fracCoords, xyzIndices, currentIndex);
                    ++currentIndex;
                }
            }
        }
        assignResiduesToCells(crystal, residues, cells);
        for (BoxOptCell cell : cells) {
            cell.sortResidues();
        }
        switch (direction) {
            case BACKWARD:
                BoxOptCell[] tempCells = new BoxOptCell[numCells];
                for (int i = 0; i < numCells; i++) {
                    tempCells[i] = cells[numCells - (i + 1)];
                }
                cells = tempCells;
            // Fall through into forward case (for now).
            case FORWARD:
            default:
                return cells;
        }
    }

    /**
     * Constructs the cells for box optimization and assigns them residues,
     * presently based on C alpha fractional coordinates; by default, cells are
     * sorted by global index. Presently, specifying approxBoxLength over-rides
     * numXYZBoxes, and always rounds the number of boxes down (to ensure boxes
     * are always at least the specified size).
     *
     * @param crystal Crystal group.
     * @param residues List of residues to be optimized.
     * @return Array of filled Cells
     */
    private void assignResiduesToCells(Crystal crystal, Residue[] residues, BoxOptCell[] cells) {
        // Search through residues, add them to all boxes containing their
        // fractional coordinates.
        int numCells = cells.length;
        int nResidues = residues.length;
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residues[i];
            double[] atomFracCoords = new double[3];
            boolean[] contained;
            double[][] originalCoordinates;
            switch (boxInclusionCriterion) {
                // As case 1 is default, test other cases first.
                case 2:
                    // Residue coordinates defined by any original atomic coordinate.
                    originalCoordinates = residuei.storeCoordinateArray();
                    contained = new boolean[numCells];
                    fill(contained, false);
                    // Loop over atomic coordinates in originalCoordinates.
                    for (int ai = 0; ai < originalCoordinates.length; ai++) {
                        crystal.toFractionalCoordinates(originalCoordinates[ai], atomFracCoords);
                        NeighborList.moveValuesBetweenZeroAndOne(atomFracCoords);
                        for (int j = 0; j < numCells; j++) {
                            if (!contained[j] && cells[j].checkIfContained(atomFracCoords)) {
                                cells[j].addResidue(residuei);
                                contained[j] = true;
                            }
                        }
                    }
                    break;
                case 3:
                    // Residue coordinates defined by any atomic coordinate in any rotamer.
                    //originalCoordinates = storeSingleCoordinates(residuei, true);
                    ResidueState origState = residuei.storeState();
                    contained = new boolean[numCells];
                    fill(contained, false);
                    Rotamer[] rotamersi = residuei.getRotamers(library);
                    for (Rotamer rotamer : rotamersi) {
                        RotamerLibrary.applyRotamer(residuei, rotamer);
                        double[][] currentCoordinates = residuei.storeCoordinateArray();
                        for (int ai = 0; ai < currentCoordinates.length; ai++) {
                            crystal.toFractionalCoordinates(currentCoordinates[ai], atomFracCoords);
                            NeighborList.moveValuesBetweenZeroAndOne(atomFracCoords);
                            for (int j = 0; j < numCells; j++) {
                                if (!contained[j] && cells[j].checkIfContained(atomFracCoords)) {
                                    cells[j].addResidue(residuei);
                                    contained[j] = true;
                                }
                            }
                        }
                    }
                    residuei.revertState(origState);
                    //revertSingleResidueCoordinates(residuei, originalCoordinates, true);
                    break;
                case 1:
                default:
                    // Residue coordinates defined by C alpha (protein) or N1/9
                    // (nucleic acids).
                    double[] cAlphaCoords = new double[3];
                    residuei.getReferenceAtom().getXYZ(cAlphaCoords);
                    crystal.toFractionalCoordinates(cAlphaCoords, atomFracCoords);
                    NeighborList.moveValuesBetweenZeroAndOne(atomFracCoords);
                    for (int j = 0; j < numCells; j++) {
                        if (cells[j].checkIfContained(atomFracCoords)) {
                            cells[j].addResidue(residuei);
                        }
                    }
                    break;
            }

        }
    }

    /**
     * Sorts an array of Residues by global index.
     *
     * @param residues Array of Residues to be sorted.
     */
    private Residue[] sortResidues(Residue[] residues) {
        int nResidues = residues.length;
        IndexIndexPair[] residueIndices = new IndexIndexPair[nResidues];
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residues[i];
            int indexOfI = allResiduesList.indexOf(residuei);
            residueIndices[i] = new IndexIndexPair(indexOfI, i);
        }
        Arrays.sort(residueIndices);
        Residue[] sortedList = new Residue[nResidues];
        for (int i = 0; i < nResidues; i++) {
            int indexToGet = residueIndices[i].getReferenceIndex();
            sortedList[i] = residues[indexToGet];
        }
        return sortedList;
    }

    /**
     * Sorts a passed ArrayList of Residues by global index.
     *
     * @param residues ArrayList of Residues to be sorted.
     */
    private void sortResidues(ArrayList<Residue> residues) {
        int nResidues = residues.size();
        IndexIndexPair[] residueIndices = new IndexIndexPair[nResidues];
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residues.get(i);
            int indexOfI = allResiduesList.indexOf(residuei);
            residueIndices[i] = new IndexIndexPair(indexOfI, i);
        }
        Arrays.sort(residueIndices);
        ArrayList<Residue> tempWindow = new ArrayList<>(residues);
        for (int i = 0; i < nResidues; i++) {
            int indexToGet = residueIndices[i].getReferenceIndex();
            residues.set(i, tempWindow.get(indexToGet));
        }
        /*List<ObjectPair<Residue, Integer>> sortedList = residues.parallelStream().
                map(r -> new ObjectPair<Residue, Integer>(r, allResiduesList.indexOf(r))).
                sorted().collect(Collectors.toList());
        // Pause here to avoid any concurrency issues when resetting the list.
        sortedList.stream().forEach(rip -> residues.set(rip.getKey(), rip.getVal())); // ArrayList is not thread-safe, so must set sequential.
        */
    }

    protected void applyEliminationCriteria(Residue residues[], boolean getEnergies, boolean verbose) {
        Level prevLevel = logger.getLevel();
        if (verbose == false) {
            logger.setLevel(Level.WARNING);
        }
        if (getEnergies) {
            applyEliminationCriteria(residues);
        } else {
            allocateEliminationMemory(residues);
            int i = 0;
            boolean pairEliminated = true;
            while (pairEliminated) {
                if (useGoldstein) {
                    i++;
                    logIfMaster(String.format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
                    // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                    while (goldsteinRotamerDriver(residues)) {
                        i++;
                        logIfMaster(this.toString());
                        logIfMaster(String.format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
                    }
                    i++;
                    logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                    // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                    pairEliminated = false;
                    while (goldsteinRotamerPairDriver(residues)) {
                        pairEliminated = true;
                        i++;
                        logIfMaster(this.toString());
                        logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                    }
                } else {
                    i++;
                    logIfMaster(String.format("\n Iteration %d: Applying Single DEE conditions ", i));
                    // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                    while (deeRotamerElimination(residues)) {
                        i++;
                        logIfMaster(toString());
                        logIfMaster(String.format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
                    }
                    i++;
                    logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                    // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                    pairEliminated = false;
                    while (deeRotamerPairElimination(residues)) {
                        pairEliminated = true;
                        i++;
                        logIfMaster(toString());
                        logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                    }
                }
                validateDEE(residues);
                logIfMaster(toString());
            }
            logIfMaster(" Self-consistent DEE rotamer elimination achieved.\n");
        }
        if (verbose == false) {
            logger.setLevel(prevLevel);
        }
    }

    private void logIfMaster(String msg) {
        if (master) {
            logger.info(msg);
        }
    }

    private void logIfMaster(String msg, Level level) {
        if (master) {
            logger.log(level, msg);
        }
    }

    private void applyEliminationCriteria(Residue residues[]) {
        allocateEliminationMemory(residues);
        /*
         * Must pin 5' ends of nucleic acids which are attached to nucleic acids
         * outside the window, to those prior residues' sugar puckers.  Then,
         * if a correction threshold is set, eliminate rotamers with excessive
         * correction vectors (up to a maximum defined by minNumberAcceptedNARotamers).
         */
        boolean containsNA = false;
        if (pruneClashes) {
            for (Residue residue : residues) {
                if (residue.getResidueType() == NA) {
                    containsNA = true;
                    break;
                }
            }
        }
        if (containsNA && pruneClashes) {
            logIfMaster(" Eliminating nucleic acid rotamers that conflict at their 5' end with residues outside the optimization range.");
            int[] numEliminatedRotamers = reconcileNARotamersWithPriorResidues(residues);
            // reconcileNARotamersWithSubsequentResidues(residues, numEliminatedRotamers);
            // Uncertain if I want to actually do this, since it could return bad results
            // if the input structure is no good.
            if (verboseEnergies) {
                logIfMaster(String.format("\n Beginning Energy %16.8f", currentEnergy(residues)));
            }
            // eliminateNABackboneRotamers checks for threshold != 0 internally.
            eliminateNABackboneRotamers(residues, numEliminatedRotamers);
            /* double naBackboneEnergy = currentEnergy();
             logger.info(String.format(" After Eliminate NA Backbone %16.8f", naBackboneEnergy)); */
        } else if (verboseEnergies) {
            logIfMaster(String.format("\n Beginning Energy %16.8f", currentEnergy(residues)));
        }

        rotamerEnergies(residues);

        // Beginning energy
        int currentRotamers[] = new int[residues.length];
        // Print starting energies if they're meaningful, i.e. rotamer-zeros are original coordinates.
        if (library.getUsingOrigCoordsRotamer()) {
            computeEnergy(residues, currentRotamers, master);
        }

        if (pruneClashes) {
            validateDEE(residues);
        }
        int i = 0;
        boolean pairEliminated = true;
        while (pairEliminated) {
            if (useGoldstein) {
                i++;
                logIfMaster(String.format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
                // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                while (goldsteinRotamerDriver(residues)) {
                    i++;
                    logIfMaster(this.toString());
                    logIfMaster(String.format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
                }
                i++;
                logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                pairEliminated = false;
                while (goldsteinRotamerPairDriver(residues)) {
                    pairEliminated = true;
                    i++;
                    logIfMaster(this.toString());
                    logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                }
            } else {
                i++;
                logIfMaster(String.format("\n Iteration %d: Applying Single DEE conditions ", i));
                // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                while (deeRotamerElimination(residues)) {
                    i++;
                    logIfMaster(toString());
                    logIfMaster(String.format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
                }
                i++;
                logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                pairEliminated = false;
                while (deeRotamerPairElimination(residues)) {
                    pairEliminated = true;
                    i++;
                    logIfMaster(toString());
                    logIfMaster(String.format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                }
            }
            validateDEE(residues);
            logIfMaster(toString());
        }
        logIfMaster(" Self-consistent DEE rotamer elimination achieved.\n");
    }

    private static void turnOnAtoms(Residue residue) {
        switch (residue.getResidueType()) {
            case NA:
            case AA:
                List<Atom> atomList = residue.getVariableAtoms();
                for (Atom atom : atomList) {
                    atom.setUse(true);
                }
                break;
            default:
                atomList = residue.getAtomList();
                for (Atom atom : atomList) {
                    atom.setUse(true);
                }
        }
    }

    private static void turnOffAtoms(Residue residue) {
        switch (residue.getResidueType()) {
            case NA:
            case AA:
                List<Atom> atomList = residue.getVariableAtoms();
                for (Atom atom : atomList) {
                    atom.setUse(false);
                }
                break;
            default:
                atomList = residue.getAtomList();
                for (Atom atom : atomList) {
                    atom.setUse(false);
                }
        }
    }

    private static void turnOnCBeta(Residue residue) {
        switch (residue.getResidueType()) {
            case AA:
                List<Atom> atomList = residue.getVariableAtoms();
                for (Atom atom : atomList) {
                    if (atom.getName().equals("CB")) {
                        atom.setUse(true);
                    }
                }
            default:
        }
    }

    /**
     * Calculates the energy at the current state.
     *
     * @param resArray Array of residues in current energy term.
     * @return Energy of the current state.
     */
    private double currentEnergy(Residue[] resArray) {
        return currentEnergy(Arrays.asList(resArray), true);
    }

    /**
     * Calculates the energy at the current state.
     *
     * @param resList List of residues in current energy term.
     * @return Energy of the current state.
     */
    private double currentEnergy(List<Residue> resList) {
        return currentEnergy(resList, true);
    }

    /**
     * Calculates the energy at the current state, with the option to throw
     * instead of catching exceptions in the energy calculation.
     *
     * @param resArray Array of residues in current energy term.
     * @param catchError If true, catch force field exceptions
     * @return Energy of the current state.
     */
    private double currentEnergy(Residue[] resArray, boolean catchError) {
        return currentEnergy(Arrays.asList(resArray), catchError);
    }

    /**
     * Calculates the energy at the current state, with the option to throw
     * instead of catching exceptions in the energy calculation.
     *
     * @param resList List of residues in current energy term.
     * @param catchError If true, catch force field exceptions.
     * @return Energy of the current state.
     */
    private double currentEnergy(List<Residue> resList, boolean catchError) {
        /*if (x == null) {
            int nVar = potential.getNumberOfVariables();
            x = new double[nVar * 3];
        }
        if (catchError) {
            try {
                potential.getCoordinates(x);
                return potential.energy(x);
            } catch (ArithmeticException ex) {
                logger.warning(ex.getMessage());
                return 1e100;
            }
        } else {
            potential.getCoordinates(x);
            return potential.energy(x);
        }*/
        List<Rotamer> rots = resList.stream().
                map(Residue::getRotamer).
                collect(Collectors.toList());
        File energyDir = dirSupplier.apply(resList, rots);
        if (catchError) {
            try {
                return eFunction.applyAsDouble(energyDir);
            } catch (ArithmeticException ex) {
                logger.warning(ex.getMessage());
                return 1e100;
            }
        } else {
            return eFunction.applyAsDouble(energyDir);
        }
    }

    /**
     * Default method for obtaining energy: calculates energy of the Potential.
     * @param dir Ignored, should be null
     * @return Current potential energy
     */
    private double currentPE(File dir) {
        if (x == null) {
            int nVar = potential.getNumberOfVariables();
            x = new double[nVar];
        }
        x = potential.getCoordinates(x);
        return potential.energy(x);
    }

    // Wrapper intended for use with RotamerMatrixMC.
    private double currentEnergyWrapper(List<Residue> resList) {
        return currentEnergy(resList);
    }

    /**
     * Sets the function used to calculate energy.
     * @param ef File to energy
     */
    public void setEnergyFunction(ToDoubleFunction<File> ef) {
        this.eFunction = ef;
    }

    /**
     * Sets the directory provider used.
     * @param dirProvider A function of residue list to appropriate directory
     */
    public void setDirectoryProvider(BiFunction<List<Residue>,List<Rotamer>,File> dirProvider) {
        this.dirSupplier = dirProvider;
    }

    /**
     * Turn off non-bonded contributions from all residues except for one.
     * Compute the self-energy for each residue relative to the backbone
     * contribution.
     *
     * @param residues A list of residues that we undergo rotamer optimization.
     * @return template energy
     */
    private double rotamerEnergies(Residue residues[]) {
        if (residues == null) {
            logger.warning(" Attempt to compute rotamer energies for an empty array of residues.");
            return 0.0;
        }
        int nResidues = residues.length;
        Atom atoms[] = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        /**
         * Initialize all atoms to be used.
         */
        for (int i = 0; i < nAtoms; i++) {
            atoms[i].setUse(true);
        }

        if (parallelEnergies) {
            if (!usingBoxOptimization) {
                energiesToWrite = Collections.synchronizedList(new ArrayList<String>());
                receiveThread = new ReceiveThread(residues);
                receiveThread.start();
                if (master && writeEnergyRestart && printFiles) {
                    energyWriterThread = new EnergyWriterThread(receiveThread);
                    energyWriterThread.start();
                }
            }
            int loaded = 0;
            if (loadEnergyRestart) {
                if (usingBoxOptimization) {
                    loaded = loadEnergyRestart(energyRestartFile, residues, boxLoadIndex, boxLoadCellIndices);
                } else {
                    loaded = loadEnergyRestart(energyRestartFile, residues);
                }
            }
            long energyTime = -System.nanoTime();
            SinglesEnergyRegion singlesRegion = new SinglesEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            PairsEnergyRegion pairsRegion = new PairsEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            TriplesEnergyRegion triplesRegion = new TriplesEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            QuadsEnergyRegion quadsRegion = new QuadsEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            try {
                if (loaded < 1) {
                    jobMapSingles.clear();
                    // allocate selfEnergy
                    int singleJobIndex = 0;
                    selfEnergy = new double[nResidues][];
                    for (int i = 0; i < nResidues; i++) {
                        Residue resi = residues[i];
                        Rotamer roti[] = resi.getRotamers(library);
                        selfEnergy[i] = new double[roti.length];
                        for (int ri = 0; ri < roti.length; ri++) {
                            if (!check(i, ri)) {
                            //logIfMaster(String.format("(sdl %d) singleJob %d: %d %d", BOXNUM, singleJobIndex, i, ri));
                                Integer selfJob[] = {i, ri};
                                if (decomposeOriginal && ri != 0) {
                                    continue;
                                }
                                jobMapSingles.put(singleJobIndex++, selfJob);
                            }
                        }
                    }
                }
                // broadcast that this proc is done with startup and allocation; ready for singles
                boolean thisProcReady = true;
                BooleanBuf thisProcReadyBuf = BooleanBuf.buffer(thisProcReady);
                multicastBuf(thisProcReadyBuf);
                // launch parallel singles calculation
                while (!readyForSingles) {
                    Thread.sleep(POLLING_FREQUENCY);
                }
                energyWorkerTeam.execute(singlesRegion);
                long timeOut = System.nanoTime() + energyTime;
                logIfMaster(String.format(" Time for single energies: %12.4g", (timeOut * 1.0E-9)));

                if (loaded < 2) {
                    jobMapPairs.clear();
                    // allocate twoBodyEnergy and create jobs
                    int pairJobIndex = 0;
                    twoBodyEnergy = new double[nResidues][][][];
                    for (int i = 0; i < nResidues; i++) {
                        Residue resi = residues[i];
                        Rotamer roti[] = resi.getRotamers(library);
                        twoBodyEnergy[i] = new double[roti.length][][];
                        for (int ri = 0; ri < roti.length; ri++) {
                            if (check(i, ri)) {
                                continue;
                            }
                            twoBodyEnergy[i][ri] = new double[nResidues][];
                            for (int j = i + 1; j < nResidues; j++) {
                                Residue resj = residues[j];
                                Rotamer rotj[] = resj.getRotamers(library);
                                twoBodyEnergy[i][ri][j] = new double[rotj.length];
                                for (int rj = 0; rj < rotj.length; rj++) {
                                    /*if (check(j, rj) || check(i, ri, j, rj)) {
                                        continue;
                                    }*/
                                    if (checkToJ(i,ri,j,rj)) {
                                        continue;
                                    }
                                    //logIfMaster(String.format("(sdl %d) pairJob %d: %d %d %d %d", BOXNUM, pairJobIndex, i, ri, j, rj));
                                    Integer pairJob[] = {i, ri, j, rj};
                                    if (decomposeOriginal && (ri != 0 || rj != 0)) {
                                        continue;
                                    }
                                    jobMapPairs.put(pairJobIndex++, pairJob);
                                }
                            }
                        }
                    }
                }
                // broadcast that this proc is done with pruning and allocation; ready for pairs
                multicastBuf(thisProcReadyBuf);
                // launch parallel twoBody calculation
                while (!readyForPairs) {
                    Thread.sleep(POLLING_FREQUENCY);
                }
                energyWorkerTeam.execute(pairsRegion);
                timeOut = System.nanoTime() + energyTime;
                logIfMaster(String.format(" Time for pair energies:   %12.4g", (timeOut * 1.0E-9)));

                if (threeBodyTerm) {
                    if (loaded < 3) {
                        jobMapTrimers.clear();
                        // allocate threeBodyEnergy and create jobs
                        int trimerJobIndex = 0;
                        threeBodyEnergy = new double[nResidues][][][][][];
                        for (int i = 0; i < nResidues; i++) {
                            Residue resi = residues[i];
                            Rotamer roti[] = resi.getRotamers(library);
                            threeBodyEnergy[i] = new double[roti.length][][][][];
                            for (int ri = 0; ri < roti.length; ri++) {
                                if (check(i, ri)) {
                                    continue;
                                }
                                threeBodyEnergy[i][ri] = new double[nResidues][][][];
                                for (int j = i + 1; j < nResidues; j++) {
                                    Residue resj = residues[j];
                                    Rotamer rotj[] = resj.getRotamers(library);
                                    threeBodyEnergy[i][ri][j] = new double[rotj.length][][];
                                    for (int rj = 0; rj < rotj.length; rj++) {
                                        /*if (check(j, rj) || check(i, ri, j, rj)) {
                                            continue;
                                        }*/
                                        if (checkToJ(i,ri,j,rj)) {
                                            continue;
                                        }
                                        threeBodyEnergy[i][ri][j][rj] = new double[nResidues][];
                                        for (int k = j + 1; k < nResidues; k++) {
                                            Residue resk = residues[k];
                                            Rotamer rotk[] = resk.getRotamers(library);
                                            threeBodyEnergy[i][ri][j][rj][k] = new double[rotk.length];
                                            for (int rk = 0; rk < rotk.length; rk++) {
                                                /*if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk) || check(i, ri, j, rj, k, rk)) {
                                                    continue;
                                                }*/
                                                if (checkToK(i,ri,j,rj,k,rk)) {
                                                    continue;
                                                }
                                                //logIfMaster(String.format("(sdl %d) trimerJob %d: %d %d %d %d %d %d", BOXNUM, trimerJobIndex, i, ri, j, rj, k, rk));
                                                Integer trimerJob[] = {i, ri, j, rj, k, rk};
                                                if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0)) {
                                                    continue;
                                                }
                                                jobMapTrimers.put(trimerJobIndex++, trimerJob);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // broadcast that this proc is done with pruning and allocation; ready for trimers
                    multicastBuf(thisProcReadyBuf);
                    // launch parallel threeBody calculation
                    while (!readyForTrimers) {
                        Thread.sleep(POLLING_FREQUENCY);
                    }
                    energyWorkerTeam.execute(triplesRegion);
                    timeOut = System.nanoTime() + energyTime;
                    logIfMaster(String.format(" Time for triple energies: %12.4g", (timeOut * 1.0E-9)));
                }

                if (computeQuads) {
                    logger.info(" Creating quad jobs...");
                    jobMapQuads.clear();
                    boolean maxedOut = false;
                    // create quad jobs (no memory allocation)
                    int quadJobIndex = 0;
                    for (int i = 0; i < nResidues; i++) {
                        Residue resi = residues[i];
                        Rotamer roti[] = resi.getRotamers(library);
                        for (int ri = 0; ri < roti.length; ri++) {
                            if (check(i, ri)) {
                                continue;
                            }
                            for (int j = i + 1; j < nResidues; j++) {
                                Residue resj = residues[j];
                                Rotamer rotj[] = resj.getRotamers(library);
                                for (int rj = 0; rj < rotj.length; rj++) {
                                    /*if (check(j, rj) || check(i, ri, j, rj)) {
                                        continue;
                                    }*/
                                    if (checkToJ(i,ri,j,rj)) {
                                        continue;
                                    }
                                    for (int k = j + 1; k < nResidues; k++) {
                                        Residue resk = residues[k];
                                        Rotamer rotk[] = resk.getRotamers(library);
                                        for (int rk = 0; rk < rotk.length; rk++) {
                                            /*if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk) || check(i, ri, j, rj, k, rk)) {
                                                continue;
                                            }*/
                                            if (checkToK(i,ri,j,rj,k,rk)) {
                                                continue;
                                            }
                                            for (int l = k + 1; l < nResidues; l++) {
                                                Residue resl = residues[l];
                                                Rotamer rotl[] = resl.getRotamers(library);
                                                for (int rl = 0; rl < rotl.length; rl++) {
                                                    /*if (check(l, rl) || check(i, ri, l, rl) ||
                                                            check(j, rj, l, rl) || check(k, rk, l, rl) ||
                                                            check(i, ri, j, rj, l, rl) || check(i, ri, k, rk, l, rl) ||
                                                            check(j, rj, k, rk, l, rl)) {
                                                        continue;
                                                    }*/
                                                    if (checkToL(i,ri,j,rj,k,rk,l,rl)) {
                                                        continue;
                                                    }
                                                    Integer quadJob[] = {i, ri, j, rj, k, rk, l, rl};
                                                    if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0 || rl != 0)) {
                                                        continue;
                                                    }
                                                    jobMapQuads.put(quadJobIndex++, quadJob);
                                                    if (jobMapQuads.size() > quadMaxout) {
                                                        maxedOut = true;
                                                        break;
                                                    }
                                                }
                                                if (maxedOut) {
                                                    break;
                                                }
                                            }
                                            if (maxedOut) {
                                                break;
                                            }
                                        }
                                        if (maxedOut) {
                                            break;
                                        }
                                    }
                                    if (maxedOut) {
                                        break;
                                    }
                                }
                                if (maxedOut) {
                                    break;
                                }
                            }
                            if (maxedOut) {
                                break;
                            }
                        }
                        if (maxedOut) {
                            break;
                        }
                    }

                    // broadcast that this proc is done with pruning and allocation; ready for quads
//                    logger.info(String.format(" Proc %d broadcasting ready for quads.", world.rank()));
                    multicastBuf(thisProcReadyBuf);
                    // launch parallel threeBody calculation
                    int waiting = 0;
                    while (!readyForQuads) {
                        Thread.sleep(POLLING_FREQUENCY);
                    }

                    logger.info(String.format(" Running quads: %d jobs.", jobMapQuads.size()));
                    energyWorkerTeam.execute(quadsRegion);
                    timeOut = System.nanoTime() + energyTime;
                    logIfMaster(String.format(" Time for quad energies:   %12.4g", (timeOut * 1.0E-9)));
                }
            } catch (Exception ex) {
                String message = " Exception computing rotamer energies in parallel.";
                logger.log(Level.SEVERE, message, ex);
            } finally {
                // stop receive and energyWriter threads
                // receiveThread.die();
            }

            energyTime += System.nanoTime();
            logIfMaster(String.format(" Time for all energies:    %12.4g", (energyTime * 1.0E-9)));
            // Turn on all atoms.
            for (int i = 0; i < atoms.length; i++) {
                atoms[i].setUse(true);
            }
            // Print the energy with all rotamers in their default conformation.
            if (verboseEnergies && master) {
                double defaultEnergy = currentEnergy(residues);
                logger.info(String.format(" Energy of the system with rotamers in their default conformation: %16.8f", defaultEnergy));
            }
            return backboneEnergy;
        }   // endif(parallelEnergies)

        // TODO: deprecate the rest of the rotamerEnergies() method; the parallel form should work (better) in all scenarios
        // Store coordinates: will be useful for nucleic acids.
        ArrayList<Residue> currentResidueList = new ArrayList<>();
        currentResidueList.addAll(Arrays.asList(residues));
        ResidueState[] origCoordinates = ResidueState.storeAllCoordinates(currentResidueList);
        //double[][][] originalAtomicCoordinates = storeCoordinates(currentResidueList);

        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer rotamers[] = residue.getRotamers(library);
            /*
             * Place each AA residue into its zeroth rotamer.
             *
             * The default state of NA residues will be original coordinates.
             * These original coordinates are usually BUT NOT ALWAYS default PDB
             * coordinates.  The reason for this is that the zeroth rotamer for
             * nucleic acids is frequently a very poor baseline, so one gets
             * enormous backbone and self energies.
             *
             * This would not affect the rigor of the DEE algorithm, but it
             * makes it difficult to tell when the program is going berserk and
             * throwing crazy energies.
             *
             * The primary exception: if increment is smaller than window size,
             * some NA residues will be in rotamers assigned by the prior window.
             * This is still a fairly reasonable baseline.
             */
            switch (residue.getResidueType()) {
                case NA:
                    break;
                case AA:
                default:
                    RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    break;
            }
            // RotamerLibrary.applyRotamer(residue, rotamers[0]);
            // Turn off all side-chain atoms that will be optimized.
            turnOffAtoms(residue);
        }

        /**
         * Initialize all atoms to be used.
         */
        boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        // Compute the backbone energy.
        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));
        try {
            backboneEnergy = currentEnergy(rList, false);
        } catch (ArithmeticException ex) {
            logger.severe(String.format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        logger.info(String.format(" Backbone energy:  %16.8f\n", backboneEnergy));
        // Compute the self-energy for each rotamer of each residue
        selfEnergy = new double[nResidues][];
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            rList.set(0, residue);
            Rotamer rotamers[] = residue.getRotamers(library);
            int nrot = rotamers.length;
            // Create space for this residue's rotamer self-energies.
            selfEnergy[i] = new double[nrot];
            // Turn on this residue's side-chain atoms
            turnOnAtoms(residue);
            // Loop over rotamers computing self-energies.
            double minEnergy = Double.MAX_VALUE;
            for (int ri = 0; ri < nrot; ri++) {
                if (pruneClashes && check(i, ri) && !(useOrigCoordsRot && ri == 0)) {
                    continue;
                }
                Rotamer rotamer = rotamers[ri];
                RotamerLibrary.applyRotamer(residue, rotamer);
                selfEnergy[i][ri] = currentEnergy(rList) - backboneEnergy;
                logger.info(String.format(" Self-energy %s %d: %16.8f", residue, ri, self(i, ri)));
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }
                if (self(i, ri) < minEnergy) {
                    minEnergy = self(i, ri);
                }
            }
            if (pruneClashes) {
                for (int ri = 0; ri < nrot; ri++) {
                    if (residue.getResidueType() == NA) {
                        if (!check(i, ri) && (self(i, ri) > (minEnergy + (clashThreshold * pruningFactor * singletonNAPruningFactor)))) {
                            logger.info(String.format(" Clash for rotamer %s %d: %16.8f >> %16.8f",
                                    residue, ri, self(i, ri), minEnergy));
                            eliminateRotamer(residues, i, ri, print);
                        }
                    } else if (!check(i, ri) && (self(i, ri) > (minEnergy + clashThreshold))) {
                        // don't prune orig-coords rotamers
                        if (!(useOrigCoordsRot && ri == 0)) {
                            logger.info(String.format(" Clash for rotamer %s %d: %16.8f >> %16.8f",
                                    residue, ri, self(i, ri), minEnergy));
                            eliminateRotamer(residues, i, ri, print);
                        }
                    }
                }
            }
            // Reset the residue to its default conformation.
            if (residue.getResidueType() == NA) {
                residue.revertState(origCoordinates[i]);
                //revertSingleResidueCoordinates(residue, originalAtomicCoordinates[i]);
            } else {
                RotamerLibrary.applyRotamer(residue, rotamers[0]);
            }

            // Turn off the side-chain
            turnOffAtoms(residue);
        }

        // Compute the pair-energy for each pair of rotamers
        twoBodyEnergy = new double[nResidues][][][];
        for (int i = 0; i < nResidues - 1; i++) {
            Residue residuei = residues[i];
            rList.set(0, residuei);
            Rotamer rotamersi[] = residuei.getRotamers(library);
            int indexI = allResiduesList.indexOf(residuei);
            // Turn on residue i
            turnOnAtoms(residuei);
            int ni = rotamersi.length;
            twoBodyEnergy[i] = new double[ni][][];
            for (int ri = 0; ri < ni; ri++) {
                if (pruneClashes && check(i, ri) && !(useOrigCoordsRot && ri == 0)) {
                    continue;
                }
                Rotamer rotameri = rotamersi[ri];
                RotamerLibrary.applyRotamer(residuei, rotameri);
                twoBodyEnergy[i][ri] = new double[nResidues][];
                for (int j = i + 1; j < nResidues; j++) {
                    Residue residuej = residues[j];
                    rList.set(1, residuej);
                    Rotamer rotamersj[] = residuej.getRotamers(library);
                    int indexJ = allResiduesList.indexOf(residuej);
                    turnOnAtoms(residuej);
                    // Loop over residue j's rotamers and compute pairwise energies.
                    int nj = rotamersj.length;
                    twoBodyEnergy[i][ri][j] = new double[nj];
                    for (int rj = 0; rj < nj; rj++) {
                        // don't prune orig-coords rotamers
                        if (pruneClashes && (check(j, rj) || check(i, ri, j, rj))
                                && !(useOrigCoordsRot && ri == 0 && rj == 0)) {
                            continue;
                        }
                        Rotamer rotamerj = rotamersj[rj];
                        RotamerLibrary.applyRotamer(residuej, rotamerj);
                        if (distanceMatrix != null) {
                            double dist = checkDistanceMatrix(indexI, ri, indexJ, rj);
                            if (dist < superpositionThreshold) {
                                twoBodyEnergy[i][ri][j][rj] = 1.0E100;
                                logger.info(String.format(" Pair energy %s %d, %s %d:   set to 1.0E100 at distance %13.6f Ang < %5.3f Ang",
                                        residuei, ri, residuej, rj, dist, superpositionThreshold));
                            } else {
                                twoBodyEnergy[i][ri][j][rj] = currentEnergy(rList) - self(i, ri) - self(j, rj) - backboneEnergy;
                                logger.info(String.format(" Pair energy %s %d, %s %d: %16.8f at distance %10.3f Ang",
                                        residuei, ri, residuej, rj, pair(i, ri, j, rj), dist));
                            }
                        } else {
                            twoBodyEnergy[i][ri][j][rj] = currentEnergy(rList)
                                    - self(i, ri) - self(j, rj) - backboneEnergy;
                            logger.info(String.format(" Pair energy %s %d, %s %d: %16.8f.",
                                    residuei, ri, residuej, rj, pair(i, ri, j, rj)));
                        }
                        if (algorithmListener != null) {
                            algorithmListener.algorithmUpdate(molecularAssembly);
                        }
                    }
                    // Reset the residue to its default conformation.
                    if (residuej.getResidueType() == NA) {
                        residuej.revertState(origCoordinates[j]);
                        //revertSingleResidueCoordinates(residuej, originalAtomicCoordinates[j]);
                    } else {
                        RotamerLibrary.applyRotamer(residuej, rotamersj[0]);
                    }
                    // RotamerLibrary.applyRotamer(residuej, rotamersj[0]);
                    turnOffAtoms(residuej);
                }
            }
            // Reset the residue to its default conformation.
            if (residuei.getResidueType() == NA) {
                residuei.revertState(origCoordinates[i]);
                //revertSingleResidueCoordinates(residuei, originalAtomicCoordinates[i]);
            } else {
                RotamerLibrary.applyRotamer(residuei, rotamersi[0]);
            }

            turnOffAtoms(residuei);
        }

        if (prunePairClashes) {
            prunePairClashes(residues);
        }
        // Prune each rotamer that clashes with all rotamers from a 2nd residue.
        /*if (prunePairClashes) {
         for (int i = 0; i < nResidues - 1; i++) {
         Residue resi = residues[i];
         Rotamer[] roti = resi.getRotamers(library);
         int ni = roti.length;
         for (int j = i + 1; j < nResidues; j++) {
         Residue resj = residues[j];
         Rotamer[] rotj = resj.getRotamers(library);
         int nj = rotj.length;
         double minPair = Double.MAX_VALUE;
         for (int ri = 0; ri < ni; ri++) {
         if (check(i, ri)) {
         continue;
         }
         for (int rj = 0; rj < nj; rj++) {
         if (check(j, rj)) {
         continue;
         }
         if (minPair > (pair(i, ri, j, rj) + self(i, ri) + self(j, rj))) {
         minPair = (pair(i, ri, j, rj) + self(i, ri) + self(j, rj));
         }
         }
         }
         // Check for elimination of any rotamer ri.
         for (int ri = 0; ri < ni; ri++) {
         if (check(i, ri)) {
         continue;
         }
         double eliminate = Double.MAX_VALUE;
         for (int rj = 0; rj < nj; rj++) {
         if (check(j, rj)) {
         continue;
         }
         if ((pair(i, ri, j, rj) + self(i, ri) + self(j, rj)) < eliminate) {
         eliminate = (pair(i, ri, j, rj) + self(i, ri) + self(j, rj));
         }
         }
         // Prune based on clash threshold and the appropriate
         // pruning factor (1.0 for AA pairs, pF for NA pairs,
         // arithmetic mean for AA-NA pairs).
         if (resi.getResidueType() == NA && resj.getResidueType() == NA) {
         if (eliminate > (minPair + (pairClashThreshold * pruningFactor))) {
         logger.info(String.format(
         " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
         resi, ri, resj, eliminate, minPair + (pairClashThreshold * pruningFactor)));
         eliminateRotamer(residues, i, ri, print);
         }
         } else if (resi.getResidueType() == NA || resj.getResidueType() == NA) {
         if (eliminate > (minPair + (pairClashThreshold * pairHalfPruningFactor))) {
         logger.info(String.format(
         " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
         resi, ri, resj, eliminate, minPair + (pairClashThreshold * pairHalfPruningFactor)));
         eliminateRotamer(residues, i, ri, print);
         }
         } else {
         if (eliminate > (minPair + pairClashThreshold)) {
         // don't prune orig-coords rotamers
         if (!(useOrigCoordsRot && ri == 0)) {
         logger.info(String.format(
         " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
         resi, ri, resj, eliminate, minPair + pairClashThreshold));
         eliminateRotamer(residues, i, ri, print);
         }
         }
         }
         }
         // Check for elimination of any rotamer rj.
         for (int rj = 0; rj < nj; rj++) {
         if (check(j, rj)) {
         continue;
         }
         double eliminate = Double.MAX_VALUE;
         for (int ri = 0; ri < ni; ri++) {
         if (check(i, ri)) {
         continue;
         }
         if ((pair(i, ri, j, rj) + self(i, ri) + self(j, rj)) < eliminate) {
         eliminate = (pair(i, ri, j, rj) + self(i, ri) + self(j, rj));
         }
         }
         if (resi.getResidueType() == NA && resj.getResidueType() == NA) {
         if (eliminate > (minPair + (pairClashThreshold * pruningFactor))) {
         logger.info(String.format(
         " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
         resj, rj, resi, eliminate, minPair + (pairClashThreshold * pruningFactor)));
         eliminateRotamer(residues, j, rj, print);
         }
         } else if (resi.getResidueType() == NA || resj.getResidueType() == NA) {
         if (eliminate > (minPair + (pairClashThreshold * pairHalfPruningFactor))) {
         logger.info(String.format(
         " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
         resj, rj, resi, eliminate, minPair + (pairClashThreshold * pairHalfPruningFactor)));
         eliminateRotamer(residues, j, rj, print);
         }
         } else {
         if (eliminate > (minPair + pairClashThreshold)) {
         // don't prune orig-coords rotamers
         if (!(useOrigCoordsRot && rj == 0)) {
         logger.info(String.format(
         " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
         resj, rj, resi, eliminate, minPair + pairClashThreshold));
         eliminateRotamer(residues, j, rj, print);
         }
         }
         }
         }
         }
         }
         } */

        // Compute the trimer-energy for each pair of rotamers
        if (threeBodyTerm) {
            double[][][][] localDistanceMatrix = new double[nResidues][][][];
            if (distance <= 0) {
                logger.info(" Calculating local distance matrix using non-eliminated rotamers.");
                localSequentialDistanceMatrix(residues, localDistanceMatrix);
            }
            threeBodyEnergy = new double[nResidues][][][][][];
            for (int i = 0; i < nResidues - 2; i++) {
                Residue residuei = residues[i];
                rList.set(0, residuei);
                Rotamer rotamersi[] = residuei.getRotamers(library);
                turnOnAtoms(residuei);
                int ni = rotamersi.length;
                threeBodyEnergy[i] = new double[ni][][][][];
                int indexOfI = allResiduesList.indexOf(residuei);
                for (int ri = 0; ri < ni; ri++) {
                    if (pruneClashes && check(i, ri) && !(ri == 0 && useOrigCoordsRot)) {
                        continue;
                    }
                    Rotamer rotameri = rotamersi[ri];
                    RotamerLibrary.applyRotamer(residuei, rotameri);
                    //int npairs = residues.length - (i + 1);
                    // TODO: reduce memory use.
                    threeBodyEnergy[i][ri] = new double[nResidues][][][];
                    for (int j = i + 1; j < nResidues - 1; j++) {
                        Residue residuej = residues[j];
                        rList.set(1, residuej);
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        turnOnAtoms(residuej);
                        // Loop over residue j's rotamers and compute pair-wise energies.
                        int nj = rotamersj.length;
                        threeBodyEnergy[i][ri][j] = new double[nj][][];
                        int indexOfJ = allResiduesList.indexOf(residuej);
                        for (int rj = 0; rj < nj; rj++) {
                            if ((pruneClashes && (check(j, rj)) || (prunePairClashes && check(i, ri, j, rj)))
                                    && !(useOrigCoordsRot && (ri == 0 && rj == 0))) {
                                continue;
                            }
                            Rotamer rotamerj = rotamersj[rj];
                            RotamerLibrary.applyRotamer(residuej, rotamerj);
                            threeBodyEnergy[i][ri][j][rj] = new double[nResidues][];
                            for (int k = j + 1; k < nResidues; k++) {
                                Residue residuek = residues[k];
                                rList.set(2, residuek);
                                Rotamer rotamersk[] = residuek.getRotamers(library);
                                turnOnAtoms(residuek);
                                int nk = rotamersk.length;
                                threeBodyEnergy[i][ri][j][rj][k] = new double[nk];
                                int indexOfK = allResiduesList.indexOf(residuek);
                                for (int rk = 0; rk < nk; rk++) {
                                    if ((pruneClashes
                                            && (check(k, rk)) || (prunePairClashes && check(i, ri, k, rk)) || (prunePairClashes && check(j, rj, k, rk)))
                                            && !(useOrigCoordsRot && (ri == 0 && rj == 0 && rk == 0))) {
                                        continue;
                                    }
                                    double dij;
                                    double dik;
                                    double djk;
                                    if (distance > 0) {
                                        /*
                                         * Distance matrix is asymmetric (it will have 4-6, but
                                         * not 6-4 stored). The present implementation sorts
                                         * windows beforehand, but this double-checks to ensure
                                         * the proper address in distanceMatrix is being called.
                                         */
                                        dij = checkDistanceMatrix(indexOfI, ri, indexOfJ, rj);
                                        dik = checkDistanceMatrix(indexOfI, ri, indexOfK, rk);
                                        djk = checkDistanceMatrix(indexOfJ, rj, indexOfK, rk);
                                        /*
                                         if (indexOfI > indexOfJ) {
                                         dij = distanceMatrix[indexOfJ][rj][indexOfI][ri];
                                         } else {
                                         dij = distanceMatrix[indexOfI][ri][indexOfJ][rj];
                                         }
                                         if (indexOfI > indexOfK) {
                                         dik = distanceMatrix[indexOfK][rk][indexOfI][ri];
                                         } else {
                                         dik = distanceMatrix[indexOfI][ri][indexOfK][rk];
                                         }
                                         if (indexOfJ > indexOfK) {
                                         djk = distanceMatrix[indexOfK][rk][indexOfJ][rj];
                                         } else {
                                         djk = distanceMatrix[indexOfJ][rj][indexOfK][rk];
                                         }
                                         */
                                    } else {
                                        dij = localDistanceMatrix[i][ri][j][rj];
                                        dik = localDistanceMatrix[i][ri][k][rk];
                                        djk = localDistanceMatrix[j][rj][k][rk];
                                    }

                                    double dist = Math.min(dij, Math.min(dik, djk));
                                    if (!threeBodyCutoff || (dist < threeBodyCutoffDist)) {
                                        if (dist < superpositionThreshold) {
                                            threeBodyEnergy[i][ri][j][rj][k][rk] = 1.0E100;
                                            logger.info(String.format(
                                                    " Trimer energy %s %d, %s %d, %s %d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang.",
                                                    residuei, ri, residuej, rj, residuek, rk, dist, superpositionThreshold));
                                        } else {
                                            Rotamer rotamerk = rotamersk[rk];
                                            RotamerLibrary.applyRotamer(residuek, rotamerk);
                                            threeBodyEnergy[i][ri][j][rj][k][rk] = currentEnergy(rList)
                                                    - self(i, ri) - self(j, rj) - self(k, rk)
                                                    - pair(i, ri, j, rj) - pair(i, ri, k, rk)
                                                    - pair(j, rj, k, rk) - backboneEnergy;
                                            logger.info(String.format(
                                                    " Trimer energy %s %d, %s %d, %s %d: %16.8f at %10.3f Ang.",
                                                    residuei, ri, residuej, rj, residuek, rk,
                                                    threeBodyEnergy[i][ri][j][rj][k][rk], dist));
                                            if (algorithmListener != null) {
                                                algorithmListener.algorithmUpdate(molecularAssembly);
                                            }
                                        }
                                    } else {
                                        logger.info(String.format(
                                                " Trimer energy %s %d, %s %d, %s %d: %16.8f at %10.3f Ang.",
                                                residuei, ri, residuej, rj, residuek, rk, 0.0,
                                                threeBodyCutoffDist));
                                    }
                                }
                                // Reset the residue to its default conformation.
                                if (residuek.getResidueType() == NA) {
                                    residuek.revertState(origCoordinates[k]);
                                    //revertSingleResidueCoordinates(residuek, originalAtomicCoordinates[k]);
                                } else {
                                    RotamerLibrary.applyRotamer(residuek, rotamersk[0]);
                                }
                                turnOffAtoms(residuek);
                            }
                        }
                        // Reset the residue to its default conformation.
                        if (residuej.getResidueType() == NA) {
                            residuej.revertState(origCoordinates[j]);
                            //revertSingleResidueCoordinates(residuej, originalAtomicCoordinates[j]);
                        } else {
                            RotamerLibrary.applyRotamer(residuej, rotamersj[0]);
                        }
                        turnOffAtoms(residuej);
                    }
                }
                // Reset the residue to its default conformation.
                if (residuei.getResidueType() == NA) {
                    residuei.revertState(origCoordinates[i]);
                    //revertSingleResidueCoordinates(residuei, originalAtomicCoordinates[i]);
                } else {
                    RotamerLibrary.applyRotamer(residuei, rotamersi[0]);
                }
                turnOffAtoms(residuei);
            }
        }

        // Turn on all atoms.
        for (int i = 0; i < nAtoms; i++) {
            atoms[i].setUse(true);
        }

        // Print the energy with all rotamers in their default conformation.
        if (verboseEnergies) {
            double defaultEnergy = currentEnergy(residues);
            logger.info(String.format(" Energy of the system with rotamers in their default conformation: %16.8f", defaultEnergy));
        }

        return backboneEnergy;
    }

    /**
     * Calculates the distance matrix; residue-residue distance is defined as
     * the shortest atom-atom distance in any possible rotamer-rotamer pair if
     * the residues are neighbors (central atom-central atom distances are
     * within a cutoff); else, distances are set to a default of
     * Double.MAX_VALUE.
     *
     * The intent of using a neighbor list is to avoid tediously searching
     * rotamer- rotamer pairs when two residues are so far apart we will never
     * need the exact distance. We use the distance matrix for adding residues
     * to the sliding window and determining whether to set 3-body energy to
     * 0.0. If the central atoms are too distant from each other, we can safely
     * assume no atoms will ever be close enough for addition to sliding window
     * or to cause explicit calculation of 3-body energy.
     */
    private void distanceMatrix() {

        distanceMatrix = new double[numResidues - 1][][][];
        long numDistances = 0L;
        for (int i = 0; i < (numResidues - 1); i++) {
            Residue residuei = allResiduesArray[i];
            int lengthRi;
            try {
                if (checkIfForced(residuei)) {
                    lengthRi = 1;
                } else {
                    lengthRi = residuei.getRotamers(library).length;
                }
            } catch (IndexOutOfBoundsException ex) {
                if (useForcedResidues) {
                    logger.warning(ex.toString());
                } else {
                    logIfMaster(String.format(" Non-forced Residue i %s has null rotamers.",
                            residuei.toString()), Level.WARNING);
                }
                continue;
            }
            distanceMatrix[i] = new double[lengthRi][][];
            for (int ri = 0; ri < lengthRi; ri++) {
                distanceMatrix[i][ri] = new double[numResidues][];
                for (int j = (i + 1); j < numResidues; j++) {
                    Residue residuej = allResiduesArray[j];
                    int lengthRj;
                    try {
                        if (checkIfForced(residuej)) {
                            lengthRj = 1;
                        } else {
                            lengthRj = residuej.getRotamers(library).length;
                        }
                    } catch (IndexOutOfBoundsException ex) {
                        if (useForcedResidues) {
                            logger.warning(ex.toString());
                        } else {
                            logIfMaster(String.format(" Residue j %s has null rotamers.", residuej.toString()));
                        }
                        continue;
                    }
                    distanceMatrix[i][ri][j] = new double[lengthRj];
                    numDistances += lengthRj;
                    if (!lazyMatrix) {
                        fill(distanceMatrix[i][ri][j], Double.MAX_VALUE);
                    } else {
                        fill(distanceMatrix[i][ri][j], -1.0);
                    }
                }
            }
        }

        logger.info(String.format(" Number of pairwise distances: %d", numDistances));

        if (!lazyMatrix) {
            //double orig[][][] = storeCoordinates(allResiduesList);
            ResidueState[] orig = ResidueState.storeAllCoordinates(allResiduesList);
            int nMultiRes = 0;

            /**
             * Build a list that contains one atom from each Residues: CA from
             * amino acids, C1 from nucleic acids, or the first atom otherwise.
             */
            Atom atoms[] = new Atom[numResidues];
            for (int i = 0; i < numResidues; i++) {
                Residue residuei = allResiduesArray[i];
                atoms[i] = residuei.getReferenceAtom();
                if (residuei instanceof MultiResidue) {
                    ++nMultiRes;
                }
            }
            /**
             * Use of the pre-existing ParallelTeam causes a conflict when
             * MultiResidues must re-init the force field. Temporary solution
             * for sequence optimization: if > 1 residue optimized, run on only
             * one thread.
             */
            int nThreads = 1;
            if (molecularAssembly.getPotentialEnergy().getParallelTeam() != null) {
                nThreads = (nMultiRes > 1) ? 1 : molecularAssembly.getPotentialEnergy().getParallelTeam().getThreadCount();
            } else {
                // Suggested: nThreads = (nMultiRes > 1) ? 1 : ParallelTeam.getDefaultThreadCount();
                nThreads = 16;
            }
            ParallelTeam parallelTeam = new ParallelTeam(nThreads);
            Crystal crystal = molecularAssembly.getCrystal();
            int nSymm = crystal.spaceGroup.getNumberOfSymOps();
            logger.info("\n Computing Residue Distance Matrix");

            NeighborList neighborList = new NeighborList(null, crystal,
                    atoms, distance + 25.0, 0.0, parallelTeam);

            // Expand coordinates
            double xyz[][] = new double[nSymm][3 * numResidues];
            double in[] = new double[3];
            double out[] = new double[3];
            for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
                for (int i = 0; i < numResidues; i++) {
                    int i3 = i * 3;
                    int iX = i3 + 0;
                    int iY = i3 + 1;
                    int iZ = i3 + 2;
                    Atom atom = atoms[i];
                    in[0] = atom.getX();
                    in[1] = atom.getY();
                    in[2] = atom.getZ();
                    crystal.applySymOp(in, out, symOp);
                    xyz[iSymOp][iX] = out[0];
                    xyz[iSymOp][iY] = out[1];
                    xyz[iSymOp][iZ] = out[2];
                }
            }

            // Build the residue neighbor-list.
            int lists[][][] = new int[nSymm][numResidues][];
            boolean use[] = new boolean[numResidues];
            fill(use, true);
            boolean forceRebuild = true;
            boolean printLists = false;
            long neighborTime = -System.nanoTime();
            neighborList.buildList(xyz, lists, use, forceRebuild, printLists);
            //neighborList.setDisableUpdates(false); // Should work, but be unnecessary.
            neighborTime += System.nanoTime();
            logger.info(String.format(" Built residue neighbor list:           %8.3f sec", neighborTime * 1.0e-9));

            DistanceRegion distanceRegion = new DistanceRegion(parallelTeam.getThreadCount(),
                    numResidues, crystal, lists, neighborList.getPairwiseSchedule());

            long parallelTime = -System.nanoTime();
            try {
                parallelTeam.execute(distanceRegion);
            } catch (Exception e) {
                String message = " Exception compting residue distance matrix.";
                logger.log(Level.SEVERE, message, e);
            }
            parallelTime += System.nanoTime();
            logger.info(String.format(" Pairwise distance matrix:              %8.3f sec\n", parallelTime * 1.0e-9));

            ResidueState.revertAllCoordinates(allResiduesList, orig);
            try {
                parallelTeam.shutdown();
            } catch (Exception ex) {
                logger.warning(String.format(" Exception shutting down parallel team for the distance matrix: %s", ex.toString()));
            }
        }
    }

    /**
     * Fills a local distance matrix based on residues passed to the method,
     * using only their non-eliminated Rotamers.
     *
     * @param residues Residues to create a distance matrix for
     * @param localDistanceMatrix Array to be filled.
     */
    private void localSequentialDistanceMatrix(Residue[] residues, double[][][][] localDistanceMatrix) {
        Crystal crystal = molecularAssembly.getCrystal();
        ArrayList<Residue> residuesList = new ArrayList<>();
        residuesList.addAll(Arrays.asList(residues));
        //double[][][] originalCoordinates = storeCoordinates(residuesList);
        ResidueState[] originalCoordinates = ResidueState.storeAllCoordinates(residuesList);

        int nSymm = crystal.spaceGroup.getNumberOfSymOps();
        int nResidues = residues.length;
        // isDistant used to skip over the ri loop if r-j distance is above the
        // 25 A center-to-center cutoff (ri, rj could not possibly be within
        // 9 A of each other).
        boolean[][] isDistant = new boolean[nResidues][nResidues];
        // Allocate memory and initialize isDistant.
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residues[i];
            Rotamer[] rotamersi = residuei.getRotamers(library);
            if (rotamersi == null) {
                continue;
            }
            int lengthRi = rotamersi.length;
            localDistanceMatrix[i] = new double[lengthRi][][];
            for (int ri = 0; ri < lengthRi; ri++) {
                if (pruneClashes && check(i, ri)) {
                    continue;
                }
                localDistanceMatrix[i][ri] = new double[nResidues][];
                for (int j = 0; j < nResidues; j++) {
                    if (i == j) {
                        continue;
                    }
                    Residue residuej = residues[j];
                    Rotamer[] rotamersj = residuej.getRotamers(library);
                    if (rotamersj == null) {
                        continue;
                    }
                    int lengthRj = rotamersj.length;
                    localDistanceMatrix[i][ri][j] = new double[lengthRj];
                    for (int rj = 0; rj < lengthRj; rj++) {
                        localDistanceMatrix[i][ri][j][rj] = Double.MAX_VALUE;
                    }
                }
            }
        }
        // Loop over symmetry operators.
        for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
            SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
            for (int i = 0; i < nResidues; i++) {
                fill(isDistant[i], false);
            }
            // Loop over residues i.  Testing showed no i-j vs. j-i discrepancies
            // over symmetry operators.
            for (int i = 0; i < (nResidues - 1); i++) {
                Residue residuei = residues[i];
                Rotamer rotamersi[] = residuei.getRotamers(library);
                if (rotamersi == null) {
                    continue;
                }
                int lengthRi = rotamersi.length;
                Atom atomi = residuei.getReferenceAtom();
                double[] refCoordsi = new double[3];
                atomi.getXYZ(refCoordsi);
                // Loop over Residue i's rotamers
                for (int ri = 0; ri < lengthRi; ri++) {
                    if (pruneClashes && check(i, ri)) {
                        continue;
                    }
                    Rotamer rotameri = rotamersi[ri];
                    RotamerLibrary.applyRotamer(residuei, rotameri);
                    double xi[][] = residuei.storeCoordinateArray();
                    // Loop over other residues j.
                    for (int j = (i + 1); j < nResidues; j++) {
                        if (i == j || isDistant[i][j]) {
                            continue;
                        }
                        Residue residuej = residues[j];
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        if (rotamersj == null) {
                            continue;
                        }
                        int lengthRj = rotamersj.length;
                        // Use reference atoms for a quick distance check, before
                        // looping over rj.
                        Atom atomj = residuej.getReferenceAtom();
                        double[] refCoordsj = new double[3];
                        atomj.getXYZ(refCoordsj);
                        crystal.applySymOp(refCoordsj, refCoordsj, symOp);
                        // dijRef is short for dijReferenceAtoms.
                        double[] dijRef = new double[3];
                        for (int k = 0; k < dijRef.length; k++) {
                            dijRef[k] = (refCoordsj[k] - refCoordsi[k]);
                        }
                        double dijSquare = crystal.image(dijRef);
                        // Corresponds to a distance of 25 A.
                        if (dijSquare > 625) {
                            isDistant[i][j] = true;
                            continue;
                        }
                        // Loop over the neighbor's rotamers
                        for (int rj = 0; rj < lengthRj; rj++) {
                            if ((pruneClashes && check(j, rj)) || (prunePairClashes && check(i, ri, j, rj))) {
                                continue;
                            }
                            Rotamer rotamerj = rotamersj[rj];
                            RotamerLibrary.applyRotamer(residuej, rotamerj);
                            double xj[][] = residuej.storeCoordinateArray();
                            double r = interResidueDistance(xi, xj, symOp);
                            if (r < localDistanceMatrix[i][ri][j][rj]) {
                                localDistanceMatrix[i][ri][j][rj] = r;
                            }
                        }
                    }
                }
            }
        }
        ResidueState.revertAllCoordinates(residuesList, originalCoordinates);
    }

    private void sequentialDistanceMatrix(Crystal crystal, int lists[][][]) {
        // Loop over symmetry operators.
        int nSymm = crystal.spaceGroup.getNumberOfSymOps();
        for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
            SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
            // Loop over residues.
            for (int i = 0; i < numResidues; i++) {
                Residue residuei = allResiduesArray[i];
                Rotamer rotamersi[] = residuei.getRotamers(library);
                if (rotamersi == null) {
                    continue;
                }
                int lengthRi = rotamersi.length;
                // Loop over Residue i's rotamers
                for (int ri = 0; ri < lengthRi; ri++) {
                    Rotamer rotameri = rotamersi[ri];
                    RotamerLibrary.applyRotamer(residuei, rotameri);
                    double xi[][] = residuei.storeCoordinateArray();
                    // Loop over Residue i's neighbors.
                    int list[] = lists[iSymOp][i];
                    int nList = list.length;
                    for (int k = 0; k < nList; k++) {
                        int j = list[k];
                        // Slight change: original was i == j
                        if (i >= j) {
                            continue;
                        }
                        Residue residuej = allResiduesArray[j];
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        if (rotamersj == null) {
                            continue;
                        }
                        int lengthRj = rotamersj.length;
                        // Loop over the neighbor's rotamers
                        for (int rj = 0; rj < lengthRj; rj++) {
                            Rotamer rotamerj = rotamersj[rj];
                            RotamerLibrary.applyRotamer(residuej, rotamerj);
                            double xj[][] = residuej.storeCoordinateArray();
                            double r = interResidueDistance(xi, xj, symOp);
                            if (r < distanceMatrix[i][ri][j][rj]) {
                                distanceMatrix[i][ri][j][rj] = r;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Calculates the minimum distance between two sets of coordinates in a
     * given symmetry operator.
     *
     * @param resi Coordinates of i by [atom][xyz]
     * @param resj Coordinates of j by [atom][xyz]
     * @param symOp Symmetry operator to apply
     * @return Minimum distance
     */
    private double interResidueDistance(double resi[][], double resj[][], SymOp symOp) {
        double dist = Double.MAX_VALUE;
        Crystal crystal = molecularAssembly.getCrystal();
        int ni = resi.length;
        for (int i = 0; i < ni; i++) {
            double xi[] = resi[i];
            int nj = resj.length;
            for (int j = 0; j < nj; j++) {
                double xj[] = resj[j];
                if (symOp != null) {
                    crystal.applySymOp(xj, xj, symOp);
                }
                // Generally: compare on square-of-distance, and square root only at return.
                //double r = Math.sqrt(crystal.image(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]));
                double r = crystal.image(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]);
                if (r < dist) {
                    dist = r;
                }
            }
        }
        return Math.sqrt(dist);
    }

    protected void allocateEliminationMemory(Residue[] residues) {
        int nres = residues.length;
        eliminatedSingles = new boolean[nres][];
        eliminatedPairs = new boolean[nres][][][];

        // Loop over residues.
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            Rotamer rotamersi[] = residuei.getRotamers(library);
            int lenri = rotamersi.length;  // Length rotamers i
            logIfMaster(String.format(" Residue %c %s with %d rotamers.", residuei.getChainID(), residuei.toString(), lenri));
            eliminatedSingles[i] = new boolean[lenri];
            eliminatedPairs[i] = new boolean[lenri][][];
            // Loop over the set of rotamers for residue i.
            for (int ri = 0; ri < lenri; ri++) {
                // int npairs = nres - (i + 1);
                // TODO - reduce memory by half.
                eliminatedSingles[i][ri] = false;
                eliminatedPairs[i][ri] = new boolean[nres][];
                for (int j = i + 1; j < nres; j++) {
                    Residue residuej = residues[j];
                    Rotamer rotamersj[] = residuej.getRotamers(library);
                    int lenrj = rotamersj.length;
                    eliminatedPairs[i][ri][j] = new boolean[lenrj];
                    for (int rj = 0; rj < lenrj; rj++) {
                        eliminatedPairs[i][ri][j][rj] = false;
                    }
                }
            }
        }
        if (threeBodyTerm) {
            eliminatedTriples = new boolean[nres][][][][][];
            for (int i = 0; i < nres; i++) {
                Residue residuei = residues[i];
                Rotamer rotamersi[] = residuei.getRotamers(library);
                int lenri = rotamersi.length;  // Length rotamers i
                eliminatedTriples[i] = new boolean[lenri][][][][];
                // Loop over the set of rotamers for residue i.
                for (int ri = 0; ri < lenri; ri++) {
                    eliminatedTriples[i][ri] = new boolean[nres][][][];
                    for (int j = i + 1; j < nres; j++) {
                        Residue residuej = residues[j];
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        int lenrj = rotamersj.length;
                        eliminatedTriples[i][ri][j] = new boolean[lenrj][][];
                        for (int rj = 0; rj < lenrj; rj++) {
                            eliminatedTriples[i][ri][j][rj] = new boolean[nres][];
                            for (int k = j + 1; k < nres; k++) {
                                Residue residuek = residues[k];
                                Rotamer rotamersk[] = residuek.getRotamers(library);
                                int lenrk = rotamersk.length;
                                eliminatedTriples[i][ri][j][rj][k] = new boolean[lenrk];
                                for (int rk = 0; rk < lenrk; rk++) {
                                    eliminatedTriples[i][ri][j][rj][k][rk] = false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Elimination of rotamers.
     */
    private boolean deeRotamerElimination(Residue[] residues) {
        int nres = residues.length;
        // A flag to indicate if any more rotamers or rotamer pairs were eliminated.
        boolean eliminated = false;
        // Loop over residues.
        double[] minMax = new double[2];
        double[] minEnergySingles = null;
        double[] maxEnergySingles = null;
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            Rotamer rotamersi[] = residuei.getRotamers(library);
            int lenri = rotamersi.length;
            if (minEnergySingles == null || minEnergySingles.length < lenri) {
                minEnergySingles = new double[lenri];
                maxEnergySingles = new double[lenri];
            }
            // Loop over the set of rotamers for residue i.
            for (int ri = 0; ri < lenri; ri++) {
                // Check for an eliminated single.
                if (check(i, ri)) {
                    continue;
                }
                // Start the min/max summation with the self-energy.
                minEnergySingles[ri] = selfEnergy[i][ri];
                maxEnergySingles[ri] = minEnergySingles[ri];
                for (int j = 0; j < nres; j++) {
                    if (j == i) {
                        continue;
                    }
                    if (minMaxPairEnergy(residues, minMax, i, ri, j)) {
                        minEnergySingles[ri] += minMax[0];
                        maxEnergySingles[ri] += minMax[1];
                    } else {
                        Residue residuej = residues[j];
                        logger.info(String.format(" Inconsistent Pair: %s %d, %s.",
                                residuei, ri, residuej));
                        eliminateRotamer(residues, i, ri, print);
                    }
                }
            }

            /**
             * Apply the singles elimination criteria to rotamers of residue i
             * by determining the most favorable maximum energy.
             */
            double eliminationEnergy = Double.MAX_VALUE;
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                if (maxEnergySingles[ri] < eliminationEnergy) {
                    eliminationEnergy = maxEnergySingles[ri];
                }
            }

            /**
             * Eliminate rotamers whose minimum energy is greater than the worst
             * case for another rotamer.
             */
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                if (minEnergySingles[ri] > eliminationEnergy + ensembleBuffer) {
                    if (eliminateRotamer(residues, i, ri, print)) {
                        eliminated = true;
                    }
                }
            }
        }
        return eliminated;
    }

    /**
     * Elimination of rotamers.
     */
    private boolean deeRotamerPairElimination(Residue[] residues) {
        int nres = residues.length;
        double minMax[] = new double[2];
        // A flag to indicate if any more rotamers or rotamer pairs were eliminated.
        boolean eliminated = false;

        double maxEnergyDoubles[] = null;
        double minEnergyDoubles[] = null;

        // Loop over residues.
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            Rotamer rotamersi[] = residuei.getRotamers(library);
            int lenri = rotamersi.length;
            // Loop over the set of rotamers for residue i.
            for (int ri = 0; ri < lenri; ri++) {
                // Check for an eliminated single.
                if (check(i, ri)) {
                    continue;
                }
                for (int j = i + 1; j < nres; j++) {
                    Residue residuej = residues[j];
                    Rotamer rotamersj[] = residuej.getRotamers(library);
                    int lenrj = rotamersj.length;

                    if (maxEnergyDoubles == null || maxEnergyDoubles.length < lenrj) {
                        maxEnergyDoubles = new double[lenrj];
                        minEnergyDoubles = new double[lenrj];
                    }

                    // Loop over residue j's rotamers.
                    for (int rj = 0; rj < lenrj; rj++) {
                        // Check for an eliminated single or pair.
                        if (check(j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        // Start the min/max summation with the "pair" self-energy.
                        minEnergyDoubles[rj] = selfEnergy[i][ri] + selfEnergy[j][rj] + pair(i, ri, j, rj);
                        maxEnergyDoubles[rj] = minEnergyDoubles[rj];
                        // Loop over the third residue.
                        for (int k = 0; k < nres; k++) {
                            if (k == i || k == j) {
                                continue;
                            }
                            if (minMaxTripleEnergy(residues, minMax, i, ri, j, rj, k)) {
                                minEnergyDoubles[rj] += minMax[0];
                                maxEnergyDoubles[rj] += minMax[1];
                            } else {
                                Residue residuek = residues[k];
                                logger.info(String.format(" Inconsistent triple: %s %d, %s %d, %s.",
                                        residuei, ri, residuej, rj, residuek));
                                eliminateRotamerPair(residues, i, ri, j, rj, print);
                            }
                        }
                    }

                    /**
                     * Apply the double elimination criteria to the rotamer pair
                     * by determining the most favorable maximum energy.
                     */
                    double pairEliminationEnergy = Double.MAX_VALUE;
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (check(j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        if (maxEnergyDoubles[rj] < pairEliminationEnergy) {
                            pairEliminationEnergy = maxEnergyDoubles[rj];
                        }
                    }
                    /**
                     * Eliminate rotamer pairs whose minimum energy is higher
                     * than the worst case for an alternative pair.
                     */
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (check(j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        if (minEnergyDoubles[rj] > pairEliminationEnergy + ensembleBuffer) {
                            logger.info(String.format(" Eliminating rotamer pair: %s %d, %s %d (%16.8f > %16.8f + %6.6f)",
                                    residuei, ri, residuej, rj,
                                    minEnergyDoubles[rj], pairEliminationEnergy, ensembleBuffer));
                            if (eliminateRotamerPair(residues, i, ri, j, rj, print)) {
                                eliminated = true;
                            }
                            // Check if any of i's rotamers are left to interact with residue j's rotamer rj?
                            boolean singleton = true;
                            for (int rii = 0; rii < lenri; rii++) {
                                if (!check(i, rii, j, rj)) {
                                    singleton = false;
                                }
                            }
                            // If not, then this rotamer is completely eliminated.
                            if (singleton) {
                                if (eliminateRotamer(residues, j, rj, print)) {
                                    eliminated = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        return eliminated;
    }

    private boolean minMaxPairEnergy(Residue[] residues, double minMax[], int i, int ri, int j) {
        Residue residuej = residues[j];
        Rotamer rotamersj[] = residuej.getRotamers(library);
        int lenrj = rotamersj.length;
        boolean valid = false;
        minMax[0] = Double.MAX_VALUE;
        minMax[1] = Double.MIN_VALUE;

        // Loop over the third residues' rotamers.
        for (int rj = 0; rj < lenrj; rj++) {
            // Check for an eliminated single or eliminated pair.
            if (check(i, ri) || check(j, rj) || check(i, ri, j, rj)) {
                continue;
            }
            double current = pair(i, ri, j, rj);
            if (threeBodyTerm) {
                double minMaxTriple[] = new double[2];
                // Loop over residue k to find the min/max triple energy.
                boolean validPair = minMax2BodySum(residues, minMaxTriple, i, ri, j, rj);
                if (!validPair) {
                    // Eliminate Rotamer Pair
                    Residue residuei = residues[i];
                    logger.info(String.format(" Inconsistent Pair: %s %d, %s %d.",
                            residuei, ri, residuej, rj));
                    /*
                     eliminatedPairs[i][ri][j][rj] = true;
                     eliminateRotamerTriples(residues, i, ri, j, rj);
                     */
                    continue;
                }
            }
            valid = true;
            if (current < minMax[0]) {
                minMax[0] = current;
            }
            if (current > minMax[1]) {
                minMax[1] = current;
            }
        }
        return valid;
    }

    private boolean minMaxTripleEnergy(Residue[] residues, double minMax[], int i, int ri, int j, int rj, int k) {
        Residue residuek = residues[k];
        Rotamer[] romatersk = residuek.getRotamers(library);
        int lenrk = romatersk.length;
        boolean valid = false;
        minMax[0] = Double.MAX_VALUE;
        minMax[1] = Double.MIN_VALUE;
        // Loop over the third residues' rotamers.
        for (int rk = 0; rk < lenrk; rk++) {
            // Check for an eliminated single or eliminated pair.
            if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk)) {
                continue;
            }
            valid = true;
            double current = pair(i, ri, k, rk) + pair(j, rj, k, rk);
            if (threeBodyTerm) {
                current += triple(i, ri, j, rj, k, rk);
                double minMaxTriple[] = new double[2];
                /**
                 * Loop over residue l to sum the min/max triple energy returns
                 * Sum_l [ min_rl [E_3(i,ri,k,rk,l,rl) + E_3(j,rj,k,rk,l,rl)]]
                 */
                boolean valid3Body = minMax3BodySum(residues, minMaxTriple, i, ri, j, rj, k, rk);
                if (!valid3Body) {
                    // Eliminate Rotamer Pair
//                    Residue residuei = residues[i];
//                    logger.info(String.format(" Eliminating Inconsistent Pair: %s %d, %s %d.",
//                            residuei, ri, residuek, rk));
//                    eliminatedPairs[i][ri][k][rk] = true;
//                    eliminateRotamerTriples(residues, i, ri, k, rk);
                    continue;
                }
                double currentMin = current + minMaxTriple[0];
                double currentMax = current + minMaxTriple[1];

                if (currentMin < minMax[0]) {
                    minMax[0] = currentMin;
                }
                if (currentMax > minMax[1]) {
                    minMax[1] = currentMax;
                }
            } else {
                if (current < minMax[0]) {
                    minMax[0] = current;
                }
                if (current > minMax[1]) {
                    minMax[1] = current;
                }
            }
        }
        return valid;
    }

    /**
     * Find the min/max of the 2-body energy.
     *
     * @param residues The residue array.
     * @param minMax The bound on the 3-body energy (minMax[0] = min, minMax[1]
     * = max.
     * @param i Residue i
     * @param ri Rotamer ri of Residue i
     * @param j Residue j
     * @param rj Rotamer rj of Residue j
     *
     * @return true if this term is valid.
     */
    private boolean minMax2BodySum(Residue[] residues, double minMax[], int i, int ri, int j, int rj) {
        int nres = residues.length;
        double minSum = 0.0;
        double maxSum = 0.0;
        for (int k = 0; k < nres; k++) {
            if (k == i || k == j) {
                continue;
            }
            Residue residuek = residues[k];
            Rotamer[] romatersk = residuek.getRotamers(library);
            int lenrk = romatersk.length;
            double currentMin = Double.MAX_VALUE;
            double currentMax = Double.MIN_VALUE;
            // Loop over the third residues' rotamers.
            boolean valid = false;
            for (int rk = 0; rk < lenrk; rk++) {
                // Check for an eliminated triple.
                if (check(i, ri, j, rj, k, rk)) {
                    continue;
                }
                valid = true;
                double current = triple(i, ri, j, rj, k, rk);
                if (current < currentMin) {
                    currentMin = current;
                }
                if (current > currentMax) {
                    currentMax = current;
                }
            }
            // Must find at least 1 valid rotamer.
            if (!valid) {
                return false;
            }
            minSum += currentMin;
            maxSum += currentMax;
        }
        minMax[0] = minSum;
        minMax[1] = maxSum;
        return true;
    }

    /**
     * Find the min/max of the 3-body energy.
     *
     * @param residues The residue array.
     * @param minMax The bound on the 3-body energy (minMax[0] = min, minMax[1]
     * = max.
     * @param i Residue i
     * @param ri Rotamer ri of Residue i
     * @param j Residue j
     * @param rj Rotamer rj of Residue j
     * @param k Residue k
     * @return true if this term is valid.
     */
    private boolean minMax3BodySum(Residue[] residues, double minMax[], int i, int ri, int j, int rj, int k, int rk) {
        int nres = residues.length;
        double minSum = 0.0;
        double maxSum = 0.0;
        for (int l = 0; l < nres; l++) {
            if (l == i || l == j || l == k) {
                continue;
            }
            Residue residuel = residues[l];
            Rotamer[] romatersl = residuel.getRotamers(library);
            int lenrl = romatersl.length;
            double currentMin = Double.MAX_VALUE;
            double currentMax = Double.MIN_VALUE;
            // Loop over the third residues' rotamers.
            boolean valid = false;
            for (int rl = 0; rl < lenrl; rl++) {
                // Check for an eliminated triple.
                if (check(i, ri, k, rk, l, rl) || check(j, rj, k, rk, l, rl)) {
                    continue;
                }
                valid = true;
                // Collect the 3-body energies and 4-body energy (TODO - returns 0.0 now)
                double current = triple(i, ri, k, rk, l, rl) + triple(j, rj, k, rk, l, rl)
                        + quad(i, ri, j, rj, k, rk, l, rl);
                if (current < currentMin) {
                    currentMin = current;
                }
                if (current > currentMax) {
                    currentMax = current;
                }
            }
            // Must find at least 1 valid rotamer.
            if (!valid) {
                return false;
            }
            minSum += currentMin;
            maxSum += currentMax;
        }
        minMax[0] = minSum;
        minMax[1] = maxSum;
        return true;
    }

    private boolean goldsteinRotamerDriver(Residue[] residues) {
        int nres = residues.length;
        // A flag to indicate if a rotamer is eliminated.
        boolean eliminated = false;
        // Loop over residue i.
        for (int i = 0; i < nres; i++) {
            Residue resi = residues[i];
            Rotamer roti[] = resi.getRotamers(library);
            int nri = roti.length;
            // Loop over the set of rotamers for residue i.
            for (int riA = 0; riA < nri; riA++) {
                // Continue if Residue i, Rotamer Ri_a is eliminated single.
                if (check(i, riA)) {
                    continue;
                }
                for (int riB = 0; riB < nri; riB++) {
                    if (check(i, riB) || riA == riB) {
                        continue;
                    }
                    if (goldsteinRotamerElimination(residues, i, riA, riB)) {
                        eliminated = true;
                        break;
                    }
                }
            }
        }
        return eliminated;
    }

    private boolean goldsteinRotamerElimination(Residue residues[], int i, int riA, int riB) {
        int nres = residues.length;
        Residue resi = residues[i];
        // Initialize Goldstein inequality.
        double goldsteinEnergy = selfEnergy[i][riA] - selfEnergy[i][riB];
        // Loop over residue j.
        for (int j = 0; j < nres; j++) {
            if (j == i) {
                continue;
            }
            Residue resj = residues[j];
            Rotamer rotj[] = resj.getRotamers(library);
            int nrj = rotj.length;
            double minForResJ = Double.MAX_VALUE;
            int rjEvals = 0;
            for (int rj = 0; rj < nrj; rj++) {
                if (check(j, rj) || check(i, riA, j, rj)) {
                    continue;
                }
                // This following is NOT checked in Osprey.
                //if (check(i, riB, j, rj)) {
                //    continue;
                //}
                rjEvals++;
                double sumOverK = pair(i, riA, j, rj) - pair(i, riB, j, rj);
                // Include three-body interactions.
                if (threeBodyTerm) {
                    for (int k = 0; k < nres; k++) {
                        if (k == i || k == j) {
                            continue;
                        }
                        Residue resk = residues[k];
                        Rotamer rotk[] = resk.getRotamers(library);
                        int nrk = rotk.length;
                        int rkEvals = 0;
                        double minForResK = Double.MAX_VALUE;
                        for (int rk = 0; rk < nrk; rk++) {
                            if (check(k, rk) || check(j, rj, k, rk)
                                    || check(i, riA, k, rk)) {
                                continue;
                            }
                            // By anaology, the following is NOT checked in Osprey.
                            // if (check(i, riB, k, rk)) {
                            //    continue;
                            // }
                            rkEvals++;
                            double e = triple(i, riA, j, rj, k, rk) - triple(i, riB, j, rj, k, rk);
                            if (e < minForResK) {
                                minForResK = e;
                            }
                        }
                        if (rkEvals == 0) {
                            // Old: return false;
                            // New: Osprey sets the minForResK to 0.0
                            minForResK = 0.0;
                        }
                        sumOverK += minForResK;
                    }
                }
                if (sumOverK < minForResJ) {
                    minForResJ = sumOverK;
                }
            }
            if (rjEvals == 0) {
                // Old: return false;
                // New: Osprey sets the minForResJ to 0.0
                minForResJ = 0.0;
            }
            goldsteinEnergy += minForResJ;
        }
        if (goldsteinEnergy > ensembleBuffer) {
            logIfMaster(String.format(" Goldstein rotamer elimination of %s %d by %d: %16.8f > %6.2f",
                    resi, riA, riB, goldsteinEnergy, ensembleBuffer));
            if (eliminateRotamer(residues, i, riA, print)) {
                return true;
            }
        }
        return false;
    }

    /**
     * The Goldstein Rotamer Pair driver routine generates pairs of rotamers to
     * be evaluated by the 3-body form of the Goldstein criteria.
     *
     * @param residues The list of residues to be optimized.
     *
     * @return true if a residue is eliminated.
     */
    private boolean goldsteinRotamerPairDriver(Residue[] residues) {
        int nres = residues.length;
        // A flag to indicate if any more rotamers or rotamer pairs were eliminated.
        boolean eliminated = false;
        // Loop over residue i.
        for (int i = 0; i < nres; i++) {
            Residue resi = residues[i];
            Rotamer roti[] = resi.getRotamers(library);
            int nri = roti.length;
            // Loop over the set of rotamers for residue i.
            for (int riA = 0; riA < nri; riA++) {
                if (check(i, riA)) {
                    continue;
                }
                // A 2nd loop over the set of rotamers for residue i.
                for (int riB = 0; riB < nri; riB++) {
                    if (riA == riB || check(i, riB)) {
                        continue;
                    }
                    // Loop over residue j.
                    for (int j = 0; j < nres; j++) {
                        if (j == i) {
                            continue;
                        }
                        Residue resj = residues[j];
                        Rotamer rotj[] = resj.getRotamers(library);
                        int nrj = rotj.length;
                        // Loop over the set of rotamers for residue j.
                        boolean breakOut = false;
                        for (int rjC = 0; rjC < nrj; rjC++) {
                            if (breakOut) {
                                break;
                            }
                            if (check(j, rjC) || check(i, riA, j, rjC)) {
                                continue;
                            }
                            // A 2nd loop over the set of rotamers for residue j.
                            for (int rjD = 0; rjD < nrj; rjD++) {
                                if (breakOut) {
                                    break;
                                }
                                if (rjC == rjD || check(j, rjD) || check(i, riB, j, rjD)) {
                                    continue;
                                }
                                // Try to eliminate R_i(riA) & R_j(rjC) using R_i(riB) & R_j(rjD)
                                if (goldsteinRotamerPairElimination(residues, i, riA, riB, j, rjC, rjD)) {
                                    eliminated = true;
                                    breakOut = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        return eliminated;
    }

    private boolean goldsteinRotamerPairElimination(Residue residues[],
            int i, int riA, int riB, int j, int rjC, int rjD) {
        int nres = residues.length;
        // Initialize the Goldstein energy.
        double goldsteinEnergy = selfEnergy[i][riA] + selfEnergy[j][rjC] + pair(i, riA, j, rjC)
                - selfEnergy[i][riB] - selfEnergy[j][rjD] - pair(i, riB, j, rjD);
        double sumOverK = 0.0;
        // Loop over a 3rd residue k.
        for (int k = 0; k < nres; k++) {
            if (k == j || k == i) {
                continue;
            }
            double minForResK = Double.MAX_VALUE;
            Residue resk = residues[k];
            Rotamer rotk[] = resk.getRotamers(library);
            int nrk = rotk.length;
            int rkEvals = 0;
            // Loop over residue k's rotamers.
            for (int rk = 0; rk < nrk; rk++) {
                if (check(k, rk) || check(i, riA, k, rk) || check(j, rjC, k, rk)) {
                    continue;
                }
                // These are not checked in Osprey.
                // if (check(i, rjB, k, rk) || check(j, rjD, k, rk)) {
                //    continue;
                // }
                rkEvals++;
                double currentResK = pair(i, riA, k, rk) + pair(j, rjC, k, rk)
                        - pair(i, riB, k, rk) - pair(j, rjD, k, rk);
                // Include 3-body effects.
                if (threeBodyTerm) {
                    double sumOverL = (triple(i, riA, j, rjC, k, rk) - triple(i, riB, j, rjD, k, rk));
                    // Loop over a 4th residue l.
                    for (int l = 0; l < nres; l++) {
                        if (l == k || l == i || l == j) {
                            continue;
                        }
                        Residue residuel = residues[l];
                        Rotamer rotamersl[] = residuel.getRotamers(library);
                        int nrl = rotamersl.length;
                        int rlEvaluations = 0;
                        double minForResL = Double.MAX_VALUE;
                        // Loop over rotamers for residue l.
                        for (int rl = 0; rl < nrl; rl++) {
                            if (check(l, rl) || check(i, riA, l, rl)
                                    || check(k, rk, l, rl) || check(i, riB, l, rl)) {
                                continue;
                            }
                            // By analogy, the following are not checked in Osprey.
                            //if (check(j, rjC, l, rl) || check(j, rjD, l, rl)) {
                            //    continue;
                            //}
                            rlEvaluations++;
                            double e = triple(i, riA, k, rk, l, rl) + triple(j, rjC, k, rk, l, rl)
                                    - triple(i, riB, k, rk, l, rl) - triple(j, rjD, k, rk, l, rl);
                            if (e < minForResL) {
                                minForResL = e;
                            }
                        }
                        if (rlEvaluations == 0) {
                            // Old: return false;
                            // New: The min over residue L does not contribute.
                            minForResL = 0.0;
                        }
                        sumOverL += minForResL;
                    }
                    currentResK += sumOverL;
                }
                if (currentResK < minForResK) {
                    minForResK = currentResK;
                }
            }
            if (rkEvals == 0) {
                // Old: return false;
                // New: The min over residue K does not contribute:
                minForResK = 0.0;
            }
            sumOverK += minForResK;
        }
        goldsteinEnergy += sumOverK;
        if (goldsteinEnergy > ensembleBuffer) {
            logIfMaster(String.format(" Goldstein pair %s %s elimination of %d %d by %d %d: %16.8f > %6.2f",
                    residues[i], residues[j], riA, rjC, riB, rjD, goldsteinEnergy, ensembleBuffer));
            if (eliminateRotamerPair(residues, i, riA, j, rjC, print)) {
                return true;
            }
        }
        return false;
    }

    private boolean eliminateRotamer(Residue[] residues, int i, int ri, boolean verbose) {
        if (!eliminatedSingles[i][ri] && !(library.getUsingOrigCoordsRotamer() && ri == 0)) {
            Residue residue = residues[i];
            logIfMaster(String.format(" Pruning rotamer: %s %d", residue, ri));
            eliminatedSingles[i][ri] = true;
            int pruned = eliminateRotamerPairs(residues, i, ri, verbose);
            if (pruned > 0) {
                logIfMaster(String.format("  Pruned %d rotamer pairs.", pruned));
            }
            return true;
        } else {
            return false;
        }
    }

    private boolean eliminateRotamerPair(Residue[] residues, int i, int ri, int j, int rj, boolean verbose) {
        if (!check(i, ri, j, rj) && !(library.getUsingOrigCoordsRotamer() && ri == 0 && rj == 0)) {
            if (i > j) {
                int ii = i;
                int iri = ri;
                i = j;
                ri = rj;
                j = ii;
                rj = iri;
            }
            Residue residuei = residues[i];
            Residue residuej = residues[j];
            if (verbose) {
                logIfMaster(String.format("  Pruning rotamer pair: %s %d %s %d",
                        residuei, ri, residuej, rj));
            }
            eliminatedPairs[i][ri][j][rj] = true;
            if (threeBodyTerm) {
                int pruned = eliminateRotamerTriples(residues, i, ri, j, rj, verbose);
                if (pruned > 0 && verbose) {
                    logIfMaster(String.format("   Pruned %d rotamer triples.", pruned));
                }
            }
            return true;
        } else {
            return false;
        }
    }

    private boolean eliminateRotamerTriple(Residue[] residues, int i, int ri, int j, int rj, int k, int rk,
            boolean verbose) {
        if (!check(i, ri, j, rj, k, rk) && !(library.getUsingOrigCoordsRotamer() && ri == 0 && rj == 0 && rk == 0)) {
            if (j < i) {
                int ii = i;
                int iri = ri;
                i = j;
                ri = rj;
                j = ii;
                rj = iri;
            }
            if (verbose) {
                logIfMaster(String.format("   Pruning rotamer triple: %s %d, %s %d, %s %d",
                        residues[i], ri, residues[j], rj, residues[k], rk));
            }
            eliminatedTriples[i][ri][j][rj][k][rk] = true;
            return true;
        } else {
            return false;
        }

    }

    private int eliminateRotamerPairs(Residue[] residues, int i, int ri, boolean verbose) {
        int nres = residues.length;
        int pruned = 0;
        for (int j = i + 1; j < nres; j++) {
            Residue residuej = residues[j];
            Rotamer rotamersj[] = residuej.getRotamers(library);
            int lenrj = rotamersj.length;
            for (int rj = 0; rj < lenrj; rj++) {
                if (eliminateRotamerPair(residues, i, ri, j, rj, verbose)) {
                    pruned++;
                }
            }
        }
        return pruned;
    }

    private int eliminateRotamerTriples(Residue[] residues, int i, int ri, int j, int rj, boolean verbose) {
        int nres = residues.length;
        int pruned = 0;
        for (int k = j + 1; k < nres; k++) {
            Residue residuek = residues[k];
            Rotamer[] rotamers = residuek.getRotamers(library);
            int lenrk = rotamers.length;
            for (int rk = 0; rk < lenrk; rk++) {
                if (eliminateRotamerTriple(residues, i, ri, j, rj, k, rk, verbose)) {
                    pruned++;
                }
            }
        }
        return pruned;
    }

    private double self(int i, int ri) {
        return selfEnergy[i][ri];
    }

    private double pair(int i, int ri, int j, int rj) {
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        double ret = 0;
        try {
            ret = twoBodyEnergy[i][ri][j][rj];
        } catch (NullPointerException ex) {
            String message = String.format("NPE for pair %d-%d, %d-%d", i, ri, j, rj);
            logger.log(Level.SEVERE, message, ex);
        }
        return ret;
    }

    private double triple(int i, int ri, int j, int rj, int k, int rk) {

        if (!threeBodyTerm) {
            return 0.0;
        }
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        if (k < i) {
            int ii = i;
            int iri = ri;
            i = k;
            ri = rk;
            k = ii;
            rk = iri;
        }
        if (k < j) {
            int jj = j;
            int jrj = rj;
            j = k;
            rj = rk;
            k = jj;
            rk = jrj;
        }

        return threeBodyEnergy[i][ri][j][rj][k][rk];
    }

    public void setThreeBodyCutoffDist(double dist) {
        this.threeBodyCutoffDist = dist;
        if (threeBodyCutoffDist < 0) {
            threeBodyCutoff = false;
        }
    }

    public void setSuperpositionThreshold(double superpositionThreshold) {
        this.superpositionThreshold = superpositionThreshold;
    }

    public void setGoldstein(boolean set) {
        this.useGoldstein = set;
    }

    /**
     * Sets level of pruning: 0 for fully off, 1 for only singles, 2 for single
     * and pair pruning.
     *
     * @param set Pruning option.
     */
    public void setPruning(int set) {
        if (set == 0) {
            this.pruneClashes = false;
            this.prunePairClashes = false;
        } else if (set == 1) {
            this.pruneClashes = true;
            this.prunePairClashes = false;
        } else if (set == 2) {
            this.pruneClashes = true;
            this.prunePairClashes = true;
        }
    }

    public void setSingletonClashThreshold(double singletonClashThreshold) {
        this.clashThreshold = singletonClashThreshold;
    }

    public void setPairClashThreshold(double pairClashThreshold) {
        this.pairClashThreshold = pairClashThreshold;
    }

    public void setWindowSize(int windowSize) {
        this.windowSize = windowSize;
        if (this.increment > windowSize) {
            logger.info(String.format(" Decreasing increment to match window size %d", windowSize));
            this.increment = windowSize;
        }
    }

    public void setVideoWriter(boolean writeVideo, boolean ignoreInactiveAtoms, boolean skipEnergies) {
        this.writeVideo = writeVideo;
        this.videoWriter = new VideoWriter(molecularAssembly, ignoreInactiveAtoms, skipEnergies);
    }

    public int getPruning() {
        if (pruneClashes && prunePairClashes) {
            return 2;
        } else if (pruneClashes && !prunePairClashes) {
            return 1;
        } else if (!pruneClashes && !prunePairClashes) {
            return 0;
        } else {
            logger.warning(" Invalid pruning state.");
            return -1;
        }
    }

    public void setEnsemble(int ensemble) {
        setEnsemble(ensemble, 5.0);
    }

    public void setEnsemble(int ensemble, double ensembleBuffer) {
        this.ensembleNumber = ensemble;
        this.ensembleBuffer = ensembleBuffer;
        this.ensembleBufferStep = 0.1 * ensembleBuffer;
        if (ensemble > 1) {
            setPruning(0);
        }
    }

    public void setEnsembleTarget(double ensembleTarget) {
        this.ensembleEnergy = ensembleTarget;
    }

    /**
     * TODO: Implement Quad energy.
     */
    private double quad(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        return 0.0;
    }

    /**
     * Check for eliminated rotamer; true if eliminated.
     *
     * @param i Residue i.
     * @param ri Rotamer ri.
     * @return True if rotamer eliminated.
     */
    protected boolean check(int i, int ri) {
        return eliminatedSingles[i][ri];
    }

    /**
     * Check for eliminated rotamer pair; true if eliminated.
     *
     * @param i Residue i.
     * @param ri Rotamer ri.
     * @param j Residue j.
     * @param rj Rotamer rj.
     * @return True if eliminated pair.
     */
    protected boolean check(int i, int ri, int j, int rj) {
        // If j is an earlier residue than i, swap j with i, as eliminated
        // rotamers are stored with the earlier residue listed first.
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        return eliminatedPairs[i][ri][j][rj];
    }

    /**
     * Check for eliminated combination of three rotamers; true if eliminated.
     *
     * @param i Residue i.
     * @param ri Rotamer ri.
     * @param j Residue j.
     * @param rj Rotamer rj.
     * @param k Residue k.
     * @param rk Rotamer rk.
     * @return True if eliminated triple.
     */
    protected boolean check(int i, int ri, int j, int rj, int k, int rk) {
        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        if (k < i) {
            int ii = i;
            int iri = ri;
            i = k;
            ri = rk;
            k = ii;
            rk = iri;
        }
        if (k < j) {
            int jj = j;
            int jrj = rj;
            j = k;
            rj = rk;
            k = jj;
            rk = jrj;
        }
        return eliminatedTriples[i][ri][j][rj][k][rk];
    }

    /**
     * Presently constitutively returns false, as we do not eliminate quads at
     * this point.
     * @param i Residue i.
     * @param ri Rotamer ri.
     * @param j Residue j.
     * @param rj Rotamer rj.
     * @param k Residue k.
     * @param rk Rotamer rk.
     * @param l Residue l
     * @param rl Rotamer rl
     * @return False (quad eliminations unimplemented)
     */
    protected boolean check(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        return false;
    }

    protected boolean checkAll(int... iris) {
        int nVals = iris.length;
        int i = iris[0];
        int ri = iris[1];
        // Fall-through switch behavior.
        switch (nVals) {
            case 8:
                int j = iris[2];
                int rj = iris[3];
                int k = iris[4];
                int rk = iris[5];
                int l = iris[6];
                int rl = iris[7];
                if (checkToL(i,ri,j,rj,k,rk,l,rl)) {
                    return true;
                }
            case 6:
                j = iris[2];
                rj = iris[3];
                k = iris[4];
                rk = iris[5];
                if (checkToK(i,ri,j,rj,k,rk)) {
                    return true;
                }
            case 4:
                j = iris[2];
                rj = iris[3];
                if (checkToJ(i,ri,j,rj)) {
                    return true;
                }
            case 2:
                return check(i,ri);
            default:
                throw new IllegalArgumentException(" checkAll function only implemented for selfs, pairs, triples, or quads");
        }
    }

    protected boolean checkTo(int... iris) {
        int nVals = iris.length;
        switch (nVals) {
            case 2:
                int i = iris[0];
                int ri = iris[1];
                return check(i,ri);
            case 4:
                i = iris[0];
                ri = iris[1];
                int j = iris[2];
                int rj = iris[3];
                return checkToJ(i,ri,j,rj);
            case 6:
                i = iris[0];
                ri = iris[1];
                j = iris[2];
                rj = iris[3];
                int k = iris[4];
                int rk = iris[5];
                return checkToK(i,ri,j,rj,k,rk);
            case 8:
                i = iris[0];
                ri = iris[1];
                j = iris[2];
                rj = iris[3];
                k = iris[4];
                rk = iris[5];
                int l = iris[6];
                int rl = iris[7];
                return checkToL(i,ri,j,rj,k,rk,l,rl);
            default:
                throw new IllegalArgumentException(" checkTo function only implemented for selfs, pairs, triples, or quads");
        }
    }

    /**
     * Wrapper for check(i,ri); ideally just call check(i,ri). Checks if residue-
     * rotamer self i,ri has been eliminated
     * @param i Residue i
     * @param ri Rotamer ri
     * @return Eliminated
     */
    protected boolean checkToI(int i, int ri) {
        return check(i,ri);
    }

    /**
     * Checks to see if any eliminations with j,rj have occurred; assumes i,ri
     * self has already been checked. Checks j,rj self and i,ri,j,rj pair. The
     * intent is to be part of a loop over i,ri,j,rj, and check for eliminations
     * at the j,rj point.
     *
     * @param i Residue i
     * @param ri Rotamer ri
     * @param j Residue j
     * @param rj Rotamer rj
     * @return j eliminated with i
     */
    protected boolean checkToJ(int i, int ri, int j, int rj) {
        return (check(j, rj) || check(i,ri,j,rj));
    }

    /**
     * Checks to see if any eliminations with k,rk have occurred; assumes
     * i,ri,j,rj pair has already been checked. Checks the k,rk self, all pairs
     * with k,rk, and the i,ri,j,rj,k,rk triple. The intent is to be part of a
     * loop over i,ri,j,rj,k,rk, and check for eliminations at the k,rk point.
     *
     * @param i Residue i
     * @param ri Rotamer ri
     * @param j Residue j
     * @param rj Rotamer rj
     * @param k Residue k
     * @param rk Rotamer rk
     * @return k eliminated with i,j
     */
    protected boolean checkToK(int i, int ri, int j, int rj, int k, int rk) {
        return (check(k,rk) || check(i,ri,k,rk) || check(j,rj,k,rk) || check(i,ri,j,rj,k,rk));
    }

    /**
     * Checks to see if any eliminations with l,rl have occurred; assumes
     * i,ri,j,rj,k,rk triple has already been checked. Checks the l,rl self, all
     * pairs with l,rl, all triples with l,rl, and the quad. The intent is to be
     * part of a loop over i,ri,j,rj,k,rk,l,rl, and check for eliminations at
     * the l,rl point.
     *
     * @param i Residue i
     * @param ri Rotamer ri
     * @param j Residue j
     * @param rj Rotamer rj
     * @param k Residue k
     * @param rk Rotamer rk
     * @param l Residue l
     * @param rl Rotamer rl
     * @return l eliminated with i,j,k
     */
    protected boolean checkToL(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        return (check(l,rl) || check(i,ri,l,rl) || check(j,rj,l,rl) || check(k,rk,l,rl) || check(i,ri,j,rj,l,rl) || check(i,ri,k,rk,l,rl) || check(j,rj,k,rk,l,rl) || check(i,ri,j,rj,k,rk,l,rl));
    }

    private boolean validateDEE(Residue residues[]) {
        int nres = eliminatedSingles.length;
        // Validate residues
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            int ni = eliminatedSingles[i].length;
            boolean valid = false;
            for (int ri = 0; ri < ni; ri++) {
                if (!check(i, ri)) {
                    valid = true;
                }
            }
            if (!valid) {
                logger.severe(String.format(" Coding error: all %d rotamers for residue %s eliminated.", ni, residuei));
            }
        }

        // Validate pairs
        for (int i = 0; i < nres; i++) {
            Residue residuei = residues[i];
            Rotamer rotamersi[] = residuei.getRotamers(library);
            int ni = rotamersi.length;
            for (int j = i + 1; j < nres; j++) {
                Residue residuej = residues[j];
                Rotamer rotamersj[] = residuej.getRotamers(library);
                int nj = rotamersj.length;
                boolean valid = false;
                for (int ri = 0; ri < ni; ri++) {
                    for (int rj = 0; rj < nj; rj++) {
                        if (!check(i, ri, j, rj)) {
                            valid = true;
                        }
                    }
                }
                if (!valid) {
                    logger.severe(String.format(" Coding error: all pairs for %s with residue %s eliminated.",
                            residuei, residuej));
                }
            }
        }

        if (threeBodyTerm) {
            // Validate triples
            for (int i = 0; i < nres; i++) {
                Residue residuei = residues[i];
                Rotamer rotamersi[] = residuei.getRotamers(library);
                int ni = rotamersi.length;
                for (int j = i + 1; j < nres; j++) {
                    Residue residuej = residues[j];
                    Rotamer rotamersj[] = residuej.getRotamers(library);
                    int nj = rotamersj.length;
                    for (int k = j + 1; k < nres; k++) {
                        Residue residuek = residues[k];
                        Rotamer rotamersk[] = residuek.getRotamers(library);
                        int nk = rotamersk.length;
                        boolean valid = false;
                        for (int ri = 0; ri < ni; ri++) {
                            for (int rj = 0; rj < nj; rj++) {
                                for (int rk = 0; rk < nk; rk++) {
                                    if (!check(i, ri, j, rj, k, rk)) {
                                        valid = true;
                                    }
                                }
                            }
                        }
                        if (!valid) {
                            logger.severe(String.format("Coding error: all triples for residues %s, %s and %s eliminated",
                                    residuei, residuej, residuek));
                        }
                    }
                }
            }
        }
        return true;
    }

    @Override
    public String toString() {
        int rotamerCount = 0;
        int pairCount = 0;
        int singles = 0;
        int pairs = 0;
        int nres = eliminatedSingles.length;
        for (int i = 0; i < nres; i++) {
            int nroti = eliminatedSingles[i].length;
            rotamerCount += nroti;
            for (int ri = 0; ri < nroti; ri++) {
                if (eliminatedSingles[i][ri]) {
                    singles++;
                }
                for (int j = i + 1; j < nres; j++) {
                    int nrotj = eliminatedPairs[i][ri][j].length;
                    pairCount += nrotj;
                    for (int rj = 0; rj < nrotj; rj++) {
                        if (eliminatedPairs[i][ri][j][rj]) {
                            pairs++;
                        }
                    }
                }
            }
        }
        StringBuilder sb = new StringBuilder(String.format(" %d out of %d rotamers eliminated.\n", singles, rotamerCount));
        sb.append(String.format(" %d out of %d rotamer pairs eliminated.", pairs, pairCount));
        return sb.toString();
    }

    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (InterruptedException e) {
                    logger.log(Level.WARNING, " Exception terminating rotamer optimization.\n", e);

                }
            }
        }
    }

    /**
     * Contains a cell used for box optimization, its residues, the fractional
     * coordinates within a crystal it takes up, its overall (linear) index, and
     * its indices along the a, b, and c crystal axes.
     */
    public class BoxOptCell {

        private ArrayList<Residue> residues = new ArrayList<>();
        // fracCoords indexed by 1-3 min x,y,z, 4-6 are max x,y,z
        private final double[] fracCoords = new double[6];
        private final int[] indexXYZ = new int[3];
        private final int linearIndex;

        /**
         * Constructs a BoxOptCell object, which takes up a set of fractional
         * coordinates within the Crystal, the Residues contained within, and
         * the index of the cell along the crystal's a, b, and c axes.
         *
         * @param fractionalCoordinates Fractional coordinates contained,
         * indexed by 1-3 min x,y,z, 4-6 max x,y,z
         * @param indices Index of cell along a, b, and c (x, y, and z).
         * @param linearIndex Index of box in linear box array.
         */
        public BoxOptCell(double[] fractionalCoordinates, int[] indices, int linearIndex) {
            System.arraycopy(fractionalCoordinates, 0, fracCoords, 0, fracCoords.length);
            System.arraycopy(indices, 0, indexXYZ, 0, indexXYZ.length);
            this.linearIndex = linearIndex;
        }

        /**
         * Constructs a BoxOptCell object, which takes up a set of fractional
         * coordinates within the Crystal, the Residues contained within, and
         * the index of the cell along the crystal's a, b, and c axes.
         *
         * @param fractionalCoordinates Fractional coordinates contained,
         * indexed by 1-3 min x,y,z, 4-6 max x,y,z
         * @param indices Index of cell along a, b, and c (x, y, and z).
         * @param linearIndex Index of box in linear box array.
         * @param residuesIn Array of Residues to initialize the BoxOptCell
         * with.
         */
        public BoxOptCell(double[] fractionalCoordinates, int[] indices, int linearIndex, Residue[] residuesIn) {
            System.arraycopy(fractionalCoordinates, 0, fracCoords, 0, fracCoords.length);
            System.arraycopy(indices, 0, indexXYZ, 0, indexXYZ.length);
            this.linearIndex = linearIndex;
            for (Residue residue : residuesIn) {
                this.residues.add(residue);
            }
        }

        /**
         * Returns an array of the Residues contained within the cell.
         *
         * @return Array of Residues.
         */
        public Residue[] getResidues() {
            Residue[] residueList = (Residue[]) residues.toArray();
            return residueList;
        }

        /**
         * Directly returns the residue list. DOES NOT RETURN A COPY.
         *
         * @return The residue list.
         */
        public ArrayList<Residue> getResidueList() {
            return residues;
        }

        /**
         * Returns a copy of the ArrayList of residues.
         *
         * @return ArrayList of Residues in the cell.
         */
        public ArrayList<Residue> getResiduesAsList() {
            ArrayList<Residue> returnedList = new ArrayList<>();
            for (Residue residue : residues) {
                returnedList.add(residue);
            }
            return returnedList;
        }

        /**
         * Get a residue in the list by index.
         *
         * @param index Index of residue to be returned.
         * @return Residue
         */
        public Residue getResidue(int index) {
            return residues.get(index);
        }

        /**
         * Gets the index of a residue inside the cell. Throws
         * IndexOutOfBoundsException if the residue is not contained within the
         * cell.
         *
         * @param residue Residue to be searched for.
         * @return Index of residue within the cell's residue list.
         * @throws IndexOutOfBoundsException If cell does not contain the
         * residue.
         */
        public int getResidueIndex(Residue residue) throws IndexOutOfBoundsException {
            for (int i = 0; i < residues.size(); i++) {
                if (residue.equals(residueList.get(i))) {
                    return i;
                }
            }
            throw new IndexOutOfBoundsException(" Residue list does not contain "
                    + "the specified residue.");
        }

        /**
         * Add a residue to the box.
         *
         * @param residue Residue to be added.
         */
        public void addResidue(Residue residue) {
            residues.add(residue);
        }

        /**
         * Adds an array of Residues to the cell.
         *
         * @param residues Array of Residues to be added.
         */
        public void addResidues(Residue[] residues) {
            for (Residue residue : residues) {
                this.residues.add(residue);
            }
        }

        /**
         * Adds an ArrayList of Residues to the cell.
         *
         * @param residueList ArrayList of Residues to be added.
         */
        public void addResidues(ArrayList<Residue> residueList) {
            this.residues.addAll(residueList);
        }

        /**
         * Get the fractional coordinates which this box contains.
         *
         * @return Fractional coordinates contained.
         */
        public double[] getFracCoords() {
            double[] returnedCoords = new double[6];
            System.arraycopy(fracCoords, 0, returnedCoords, 0, returnedCoords.length);
            return returnedCoords;
        }

        /**
         * Returns the x, y, and z indices of this box.
         *
         * @return Box indices.
         */
        public int[] getXYZIndex() {
            int[] returnedIndices = new int[3];
            System.arraycopy(indexXYZ, 0, returnedIndices, 0, returnedIndices.length);
            return returnedIndices;
        }

        /**
         * Returns the linear index of this Box.
         *
         * @return Linear index.
         */
        public int getLinearIndex() {
            return linearIndex;
        }

        /**
         * Sorts residues in the cell by index in allResiduesList. Identical to
         * RotamerOptimization.sortResidues().
         */
        public void sortResidues() {
            int nResidues = residues.size();
            IndexIndexPair[] residueIndices = new IndexIndexPair[nResidues];
            for (int i = 0; i < nResidues; i++) {
                Residue residuei = residues.get(i);
                int indexOfI = allResiduesList.indexOf(residuei);
                residueIndices[i] = new IndexIndexPair(indexOfI, i);
            }
            Arrays.sort(residueIndices);
            ArrayList<Residue> tempWindow = new ArrayList<>(residues);
            for (int i = 0; i < nResidues; i++) {
                int indexToGet = residueIndices[i].getReferenceIndex();
                residues.set(i, tempWindow.get(indexToGet));
            }
            /*residues = residues.parallelStream().
                    map(r -> new ObjectPair<Residue, Integer>(r, allResiduesList.indexOf(r))).
                    sorted().
                    map(rip -> rip.getVal()).
                    collect(Collectors.toCollection(ArrayList<Residue>::new));*/
        }

        /**
         * Clears the list of residues in the cell.
         */
        public void clearResidues() {
            residues = new ArrayList<>();
        }

        /**
         * Removes a residue in the cell by index in the cell.
         *
         * @param index Residue to be removed.
         */
        public void removeResidue(int index) {
            residues.remove(index);
        }

        /**
         * Removes a specified residue in the cell.
         *
         * @param residue Removed from cell.
         */
        public void removeResidue(Residue residue) {
            residues.remove(residue);
        }

        /**
         * Check if an atom's fractional coordinates would be contained by the
         * box.
         *
         * @param atomFracCoords Atomic fractional coordinates
         * @return If contained.
         */
        public boolean checkIfContained(double[] atomFracCoords) {
            for (int i = 0; i < 3; i++) {
                if ((fracCoords[i] > atomFracCoords[i]) || (fracCoords[i + 3] < atomFracCoords[i])) {
                    return false;
                }
            }
            return true;
        }
    }

    public void setParallelEnergies(boolean set) {
        this.parallelEnergies = set;
    }

    public void setEnergyRestartFile(File file) {
        loadEnergyRestart = true;
        energyRestartFile = file;
    }

    public int nameToNumber(String residueString, Residue residues[]) throws NumberFormatException {
        int ret = -1;
        for (int x = 0; x < residues.length; x++) {
            if (residueString.equals(residues[x].toString())) {
                ret = x;
                break;
            }
        }
        if (ret == -1) {
            throw new NumberFormatException();
        }
        return ret;
    }

    /**
     * Loads the number of iterations a box optimization or sliding window has
     * gone through if using an energy restart file: NOT IMPLEMENTED.
     *
     * @param restartFile
     * @return Number of iterations of the sliding window or box optimization.
     */
    public int loadEnergyRestartIterations(File restartFile) {
        int numIterations = -1;
        try {
            Path path = Paths.get(restartFile.getCanonicalPath());
            List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
            for (String line : lines) {
                String tok[] = line.split("\\s");
                if (tok[0].startsWith("Iteration")) {
                    numIterations = Integer.parseInt(tok[1]);
                    return numIterations;
                }
            }
        } catch (IOException ex) {
            logIfMaster(String.format(" Error in loading number of iterations of box optimization or sliding window: %s", ex.toString()), Level.WARNING);
        }
        return numIterations;
    }

    public int loadEnergyRestart(File restartFile, Residue residues[]) {
        return loadEnergyRestart(restartFile, residues, -1, null);
    }

    public int loadEnergyRestart(File restartFile, Residue residues[], int boxIteration, int[] cellIndices) {
        try {
            int nResidues = residues.length;
            Path path = Paths.get(restartFile.getCanonicalPath());
            List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
            List<String> linesThisBox = new ArrayList<>();

            if (usingBoxOptimization && boxIteration >= 0) {
                boolean foundBox = false;
                for (int i = 0; i < lines.size(); i++) {
                    String line = lines.get(i);
                    if (line.startsWith("Box")) {
                        String tok[] = line.replaceAll("Box", "").replaceAll(":", ",").replaceAll(" ", "").split(",");
                        int readIteration = Integer.parseInt(tok[0]);
                        int readCellIndexX = Integer.parseInt(tok[1]);
                        int readCellIndexY = Integer.parseInt(tok[2]);
                        int readCellIndexZ = Integer.parseInt(tok[3]);
                        if (readIteration == boxIteration
                                && readCellIndexX == cellIndices[0]
                                && readCellIndexY == cellIndices[1]
                                && readCellIndexZ == cellIndices[2]) {
                            foundBox = true;
                            for (int j = i + 1; j < lines.size(); j++) {
                                String l = lines.get(j);
                                if (l.startsWith("Box")) {
                                    break;
                                }
                                linesThisBox.add(l);
                            }
                            break;
                        }
                    }
                }
                if (!foundBox) {
                    logIfMaster(String.format(" Didn't find restart energies for Box %d: %d,%d,%d",
                            boxIteration, cellIndices[0], cellIndices[1], cellIndices[2]));
                    return 0;
                } else if (linesThisBox.size() == 0) {
                    return 0;
                } else {
                    lines = linesThisBox;
                }
            }

            List<String> singleLines = new ArrayList<>();
            List<String> pairLines = new ArrayList<>();
            List<String> tripleLines = new ArrayList<>();
            for (String line : lines) {
                String tok[] = line.split("\\s");
                if (tok[0].startsWith("Self")) {
                    singleLines.add(line);
                } else if (tok[0].startsWith("Pair")) {
                    pairLines.add(line);
                } else if (tok[0].startsWith("Triple")) {
                    tripleLines.add(line);
                }
            }
            int loaded = 0;
            if (tripleLines.size() > 0) {
                loaded = 3;
            } else if (pairLines.size() > 0) {
                loaded = 2;
            } else if (singleLines.size() > 0) {
                loaded = 1;
            } else {
                logger.warning(String.format("Empty or unreadable energy restart file: %s.", restartFile.getCanonicalPath()));
            }
            if (loaded >= 1) {
                jobMapSingles.clear();
                // allocate selfEnergy array and create self jobs
                HashMap<String, Integer> reverseJobMapSingles = new HashMap<>();
                int singleJobIndex = 0;
                selfEnergy = new double[nResidues][];
                for (int i = 0; i < nResidues; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    selfEnergy[i] = new double[roti.length];
                    for (int ri = 0; ri < roti.length; ri++) {
                        Integer selfJob[] = {i, ri};
                        if (decomposeOriginal && ri != 0) {
                            continue;
                        }
                        jobMapSingles.put(singleJobIndex, selfJob);
                        String revKey = String.format("%d %d", i, ri);
                        reverseJobMapSingles.put(revKey, singleJobIndex);
                        singleJobIndex++;
                    }
                }
                // fill in self-energies from file while removing the corresponding jobs from jobMapSingles
                for (String line : singleLines) {
                    try {
                        String tok[] = line.replace(",", "").replace(":", "").split("\\s+");
                        int i;
                        if (tok[1].contains("-")) {
                            i = nameToNumber(tok[1], residues);
                        } else {
                            i = Integer.parseInt(tok[1]);
                        }
                        int ri = Integer.parseInt(tok[2]);
                        double energy = Double.parseDouble(tok[3]);
                        try {
                            selfEnergy[i][ri] = energy;
                        } catch (ArrayIndexOutOfBoundsException ex) {
                            logger.log(Level.SEVERE, String.format("Restart file contained an out-of-bounds index.  Offending line: %s", line), ex);
                        }
                        if (verbose) {
                            logIfMaster(String.format(" From restart file: Self energy %3d %-2d: %7s %-2d: %16.8f", i, ri, residues[i], ri, energy));
                        }
                        // remove that job from the pool
                        String revKey = String.format("%d %d", i, ri);
                        Integer ret[] = jobMapSingles.remove(reverseJobMapSingles.get(revKey));
                        if (ret == null) {
                            //logIfMaster(String.format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, String.format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded self energies from restart file.");
                // prune singles
                if (pruneClashes) {
                    pruneSingleClashes(residues);
                }
            }
            if (loaded >= 2) {
                if (jobMapSingles.size() > 0) {
                    if (master) {
                        logger.warning("Double-check that parameters match original run!  Found pairs in restart file, but singles job queue is non-empty.");
                    }
                }
                jobMapPairs.clear();
                // allocated twoBodyEnergy array and create pair jobs
                HashMap<String, Integer> reverseJobMapPairs = new HashMap<>();
                int pairJobIndex = 0;
                twoBodyEnergy = new double[nResidues][][][];
                for (int i = 0; i < nResidues; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    twoBodyEnergy[i] = new double[roti.length][][];
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (pruneClashes && check(i, ri)) {
                            continue;
                        }
                        twoBodyEnergy[i][ri] = new double[nResidues][];
                        for (int j = i + 1; j < nResidues; j++) {
                            Residue resj = residues[j];
                            Rotamer rotj[] = resj.getRotamers(library);
                            twoBodyEnergy[i][ri][j] = new double[rotj.length];
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if ((pruneClashes && check(j, rj)) || (prunePairClashes && check(i, ri, j, rj))) {
                                    continue;
                                }
                                Integer pairJob[] = {i, ri, j, rj};
                                if (decomposeOriginal && (ri != 0 || rj != 0)) {
                                    continue;
                                }
                                jobMapPairs.put(pairJobIndex, pairJob);
                                String revKey = String.format("%d %d %d %d", i, ri, j, rj);
                                reverseJobMapPairs.put(revKey, pairJobIndex);
                                pairJobIndex++;
                            }
                        }
                    }
                }
                // fill in pair-energies from file while removing the corresponding jobs from jobMapPairs
                for (String line : pairLines) {
                    try {
                        String tok[] = line.replace(",", "").replace(":", "").split("\\s+");
                        int i;
                        if (tok[1].contains("-")) {
                            i = nameToNumber(tok[1], residues);
                        } else {
                            i = Integer.parseInt(tok[1]);
                        }
                        int ri = Integer.parseInt(tok[2]);
                        int j;
                        if (tok[3].contains("-")) {
                            j = nameToNumber(tok[3], residues);
                        } else {
                            j = Integer.parseInt(tok[3]);
                        }
                        int rj = Integer.parseInt(tok[4]);
                        double energy = Double.parseDouble(tok[5]);
                        try {
                            twoBodyEnergy[i][ri][j][rj] = energy;
                        } catch (ArrayIndexOutOfBoundsException ex) {
                            logger.log(Level.SEVERE, String.format("Restart file contained an out-of-bounds index.  Offending line: %s", line), ex);
                        }
                        if (verbose) {
                            logIfMaster(String.format(" From restart file: Pair energy %3d %-2d, %3d %-2d: %16.8f", i, ri, j, rj, energy));
                        }
                        // remove that job from the pool
                        String revKey = String.format("%d %d %d %d", i, ri, j, rj);
                        Integer ret[] = jobMapPairs.remove(reverseJobMapPairs.get(revKey));
                        if (ret == null) {
                            //logIfMaster(String.format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, String.format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded pair energies from restart file.");
                // prune pairs
                if (prunePairClashes) {
                    prunePairClashes(residues);
                }
            }
            if (loaded >= 3) {
                if (jobMapPairs.size() > 0) {
                    if (master) {
                        logger.warning("Double-check that parameters match original run!  Found trimers in restart file, but pairs job queue is non-empty.");
                    }
                }
                HashMap<String, Integer> reverseJobMapTrimers = new HashMap<>();
                jobMapTrimers.clear();
                // allocate threeBodyEnergy array, fill in triple-energies from file
                int trimerJobIndex = 0;
                threeBodyEnergy = new double[nResidues][][][][][];
                for (int i = 0; i < nResidues; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    threeBodyEnergy[i] = new double[roti.length][][][][];
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (pruneClashes && check(i, ri)) {
                            continue;
                        }
                        threeBodyEnergy[i][ri] = new double[nResidues][][][];
                        for (int j = i + 1; j < nResidues; j++) {
                            Residue resj = residues[j];
                            Rotamer rotj[] = resj.getRotamers(library);
                            threeBodyEnergy[i][ri][j] = new double[rotj.length][][];
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if ((pruneClashes && check(j, rj)) || (prunePairClashes && check(i, ri, j, rj))) {
                                    continue;
                                }
                                threeBodyEnergy[i][ri][j][rj] = new double[nResidues][];
                                for (int k = j + 1; k < nResidues; k++) {
                                    Residue resk = residues[k];
                                    Rotamer rotk[] = resk.getRotamers(library);
                                    threeBodyEnergy[i][ri][j][rj][k] = new double[rotk.length];
                                    for (int rk = 0; rk < rotk.length; rk++) {
                                        if ((pruneClashes && check(k, rk)) || (prunePairClashes && (check(i, ri, k, rk) || check(j, rj, k, rk)))) {
                                            continue;
                                        }
                                        Integer trimerJob[] = {i, ri, j, rj, k, rk};
                                        if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0)) {
                                            continue;
                                        }
                                        jobMapTrimers.put(trimerJobIndex, trimerJob);
                                        String revKey = String.format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                                        reverseJobMapTrimers.put(revKey, trimerJobIndex);
                                        trimerJobIndex++;
                                    }
                                }
                            }
                        }
                    }
                }
                // fill in triple-energies from file while removing the corresponding jobs from jobMapTrimers
                for (String line : tripleLines) {
                    try {
                        String tok[] = line.replace(",", "").replace(":", "").split("\\s+");
                        int i;
                        if (tok[1].contains("-")) {
                            i = nameToNumber(tok[1], residues);
                        } else {
                            i = Integer.parseInt(tok[1]);
                        }
                        int ri = Integer.parseInt(tok[2]);
                        int j;
                        if (tok[3].contains("-")) {
                            j = nameToNumber(tok[3], residues);
                        } else {
                            j = Integer.parseInt(tok[3]);
                        }
                        int rj = Integer.parseInt(tok[4]);
                        int k;
                        if (tok[5].contains("-")) {
                            k = nameToNumber(tok[5], residues);
                        } else {
                            k = Integer.parseInt(tok[5]);
                        }
                        int rk = Integer.parseInt(tok[6]);
                        double energy = Double.parseDouble(tok[7]);
                        try {
                            threeBodyEnergy[i][ri][j][rj][k][rk] = energy;
                        } catch (ArrayIndexOutOfBoundsException ex) {
                            logger.log(Level.SEVERE, String.format("Restart file contained an out-of-bounds index.  Offending line: %s", line), ex);
                        }
                        if (verbose) {
                            logIfMaster(String.format(" From restart file: Trimer energy %3d %-2d, %3d %-2d, %3d %-2d: %16.8f", i, ri, j, rj, k, rk, energy));
                        }
                        // remove that job from the pool
                        String revKey = String.format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                        Integer ret[] = jobMapTrimers.remove(reverseJobMapTrimers.get(revKey));
                        if (ret == null) {
                            //logIfMaster(String.format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, String.format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded trimer energies from restart file.");
            }
            return loaded;
        } catch (IOException ex) {
            logger.log(Level.WARNING, "Exception while loading energy restart file.", ex);
        }
        return 0;
    }

    /**
     * Writes files that allow restarting from partially-completed call to
     * rotamerEnergies(). Spawned only in a parallel environment and only by the
     * master process.
     */
    private class EnergyWriterThread extends Thread {

        private ReceiveThread receiveThread;
        private File restartFile;
        private BufferedWriter bw;
        private final int writeFrequency = 100;
        private String boxHeader = null;

        public EnergyWriterThread(ReceiveThread receiveThread) {
            this.receiveThread = receiveThread;
            if (loadEnergyRestart) {
                restartFile = energyRestartFile;
            } else {
                File file = molecularAssembly.getFile();
                String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
                Path restartPath = Paths.get(filename + ".restart");
                // if there's already one there, back it up before starting a new one
                // this happens when wrapping globalUsingEliminations with e.g. box optimization
//                if (Files.exists(restartPath)) {
//                    int i = 0;
//                    Path backupRestartPath = Paths.get(filename + ".resBak");
//                    while (Files.exists(backupRestartPath)) {
//                        backupRestartPath = Paths.get(filename + ".resBak" + ++i);
//                    }
//                    try {
//                        FileUtils.moveFile(restartPath.toFile(), backupRestartPath.toFile());
//                    } catch (IOException ex) {
//                        logger.warning("I/O exception while backing up existing restart file.");
//                    }
//                }
                restartFile = restartPath.toFile();
            }
            try {
                bw = new BufferedWriter(new FileWriter(restartFile, true));
            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Couldn't open energy restart file.", ex);
            }
            logger.info(String.format(" Energy restart file: %s", restartFile.getName()));
        }

        public EnergyWriterThread(ReceiveThread receiveThread, int iteration, int[] cellIndices) {
            this.receiveThread = receiveThread;
            if (loadEnergyRestart) {
                restartFile = energyRestartFile;
            } else {
                File file = molecularAssembly.getFile();
                String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
                Path restartPath = Paths.get(filename + ".restart");
                restartFile = restartPath.toFile();
            }
            try {
                bw = new BufferedWriter(new FileWriter(restartFile, true));
                boxHeader = String.format("Box %d: %d,%d,%d", iteration, cellIndices[0], cellIndices[1], cellIndices[2]);
            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Couldn't open energy restart file.", ex);
            }
            logger.info(String.format(" Energy restart file: %s", restartFile.getName()));
        }

        @Override
        public void run() {
            boolean die = false;
            List<String> writing = new ArrayList<>();
            while (!die) {
                if (receiveThread.getState() == java.lang.Thread.State.TERMINATED) {
                    die = true;
                }
                if (energiesToWrite.size() >= writeFrequency || die) {
                    // copy energies to new local array; avoid blocking the receive thread for too long
                    synchronized (energiesToWrite) {
                        writing.addAll(energiesToWrite);
                        energiesToWrite.clear();
                    }
                    try {
                        if (boxHeader != null && !writing.isEmpty()) {
                            bw.append(boxHeader);
                            bw.newLine();
                            bw.flush();
                            boxHeader = null;
                        }
                        for (String line : writing) {
                            bw.append(line);
                            bw.newLine();
                        }
                        bw.flush();
                    } catch (IOException ex) {
                        logger.log(Level.SEVERE, "Exception writing energy restart file.", ex);
                    }
                    writing.clear();
                }
            }
            try {
                bw.close();
            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Exception while closing energy restart file.", ex);
            }
            if (verbose) {
                logger.info(" EnergyWriterThread shutting down.");
            }
        }
    }

    private class ReceiveThread extends Thread {

        private double incSelf[], incPair[], incTriple[];
        private DoubleBuf incSelfBuf, incPairBuf, incTripleBuf;
        private boolean alive = true;

        public ReceiveThread(Residue residues[]) {
            incSelf = new double[3];
            incSelfBuf = DoubleBuf.buffer(incSelf);
            incPair = new double[5];
            incPairBuf = DoubleBuf.buffer(incPair);
            incTriple = new double[7];
            incTripleBuf = DoubleBuf.buffer(incTriple);
        }

        public void die() {
            logger.info(" ReceiveThread registered die() command.");
            alive = false;
        }

        @Override
        public void run() {
            // Barrier; wait for everyone to be done starting up and allocating singleEnergy arrays.
            int procsDone = 0;
            boolean readyForNext = false;
            BooleanBuf readyForNextBuf = BooleanBuf.buffer(readyForNext);
            while (alive) {
                try {
                    CommStatus cs = world.receive(null, readyForNextBuf);
                    if (cs.length > 0) {
                        procsDone++;    // boolean check unnecessary; message is only sent by ready processes
                    }
                } catch (IOException ex) {
                    logger.log(Level.WARNING, ex.getMessage(), ex);
                }
                if (procsDone >= numProc) {
                    readyForSingles = true;
                    break;
                }
            }

            // Receive all singles energies.
            procsDone = 0;
            while (alive) {
                try {
                    CommStatus cs = world.receive(null, incSelfBuf);
                    if (cs.length > 0) {
                        int resi = (int) incSelf[0];
                        int roti = (int) incSelf[1];
                        double energy = incSelf[2];
                        // check for "process finished" announcements
                        if (resi < 0 && roti < 0) {
                            procsDone++;
                        } else {
                            selfEnergy[resi][roti] = energy;
                            if (writeEnergyRestart && printFiles) {
                                energiesToWrite.add(String.format("Self %d %d: %16.8f", resi, roti, energy));
                            }
                        }
                    }
                } catch (IOException ex) {
                    logger.log(Level.WARNING, ex.getMessage(), ex);
                }

                if (procsDone >= numProc) {
                    selfsDone = true;
                    break;
                }
            }

            // Barrier; wait for everyone to be done pruning and allocating twoBodyEnergy arrays before computing pairs.
            procsDone = 0;
            while (alive) {
                try {
                    CommStatus cs = world.receive(null, readyForNextBuf);
                    if (cs.length > 0) {
                        procsDone++;    // boolean check unnecessary; message is only sent by ready processes
                    }
                } catch (IOException ex) {
                    logger.log(Level.WARNING, ex.getMessage(), ex);
                }
                if (procsDone >= numProc) {
                    readyForPairs = true;
                    break;
                }
            }

            // Receive all pair energies.
            procsDone = 0;
            while (alive) {
                try {
                    CommStatus cs = world.receive(null, incPairBuf);
                    if (cs.length > 0) {
                        int resi = (int) incPair[0];
                        int roti = (int) incPair[1];
                        int resj = (int) incPair[2];
                        int rotj = (int) incPair[3];
                        double energy = incPair[4];
                        // check for "process finished" announcements
                        if (resi < 0 && roti < 0 && resj < 0 && rotj < 0) {
                            procsDone++;
                        } else {
                            twoBodyEnergy[resi][roti][resj][rotj] = energy;
                            if (writeEnergyRestart && printFiles) {
                                energiesToWrite.add(String.format("Pair %d %d, %d %d: %16.8f", resi, roti, resj, rotj, energy));
                            }
                        }
                    }
                } catch (IOException ex) {
                    logger.log(Level.WARNING, ex.getMessage(), ex);
                }

                if (procsDone >= numProc) {
                    pairsDone = true;
                    break;
                }
            }

            if (threeBodyTerm) {
                // Barrier; wait for everyone to be done pruning and allocating threeBodyEnergy arrays before computing trimers.
                procsDone = 0;
                while (alive) {
                    readyForNext = false;
                    try {
                        CommStatus cs = world.receive(null, readyForNextBuf);
                        if (cs.length > 0) {
                            procsDone++;    // boolean check unnecessary; message is only sent by ready processes
                        }
                    } catch (IOException ex) {
                        logger.log(Level.WARNING, ex.getMessage(), ex);
                    }
                    if (procsDone >= numProc) {
                        readyForTrimers = true;
                        break;
                    }
                }

                // Receive all trimer energies.
                procsDone = 0;
                while (alive) {
                    try {
                        CommStatus cs = world.receive(null, incTripleBuf);
                        if (cs.length > 0) {
                            int resi = (int) incTriple[0];
                            int roti = (int) incTriple[1];
                            int resj = (int) incTriple[2];
                            int rotj = (int) incTriple[3];
                            int resk = (int) incTriple[4];
                            int rotk = (int) incTriple[5];
                            double energy = incTriple[6];
                            // check for "process finished" announcements
                            if (resi < 0 && roti < 0 && resj < 0 && rotj < 0 && resk < 0 && rotk < 0) {
                                procsDone++;
                            } else {
                                threeBodyEnergy[resi][roti][resj][rotj][resk][rotk] = energy;
                                if (writeEnergyRestart && printFiles) {
                                    energiesToWrite.add(String.format("Triple %d %d, %d %d, %d %d: %16.8f", resi, roti, resj, rotj, resk, rotk, energy));
                                }
                            }
                        }
                    } catch (IOException ex) {
                        logger.log(Level.WARNING, ex.getMessage(), ex);
                    }

                    if (procsDone >= numProc) {
                        trimersDone = true;
                        break;
                    }
                }
            }

            if (computeQuads) {
                // Barrier; wait for everyone to be done pruning before starting quad energies.
                procsDone = 0;
                while (alive) {
                    readyForNext = false;
                    try {
                        CommStatus cs = world.receive(null, readyForNextBuf);
                        if (cs.length > 0) {
                            procsDone++;    // boolean check unnecessary; message is only sent by ready processes
                        }
                    } catch (IOException ex) {
                        logger.log(Level.WARNING, ex.getMessage(), ex);
                    }
                    if (procsDone >= numProc) {
                        readyForQuads = true;
                        break;
                    }
                }
                // No receive mechanism for quad energies since we only print them.
            }
            if (verbose) {
                logger.info(" Receive thread shutting down.");
            }
        }
    }

    /**
     * Uses calculated energies to prune rotamers based on a threshold distance
     * from that residue's minimum energy rotamer (by default 20 kcal/mol). The
     * threshold can be modulated by presence of nucleic acids or MultiResidues,
     * which require more generous pruning criteria.
     *
     * @param residues Residues to prune rotamers over.
     */
    public void pruneSingleClashes(Residue residues[]) {
        if (!pruneClashes) {
            return;
        }
        for (int i = 0; i < residues.length; i++) {
            Residue residue = residues[i];
            Rotamer rotamers[] = residue.getRotamers(library);
            int nrot = rotamers.length;
            double minEnergy = Double.MAX_VALUE;
            for (int ri = 0; ri < nrot; ri++) {
                if (!check(i, ri) && self(i, ri) < minEnergy) {
                    minEnergy = self(i, ri);
                }
            }
            /**
             * Regular: ep = minEnergy + clashThreshold Nucleic acids: ep =
             * minEnergy + (clashThreshold * factor * factor) MultiResidues: ep
             * = minEnergy + multiResClashThreshold
             *
             * Nucleic acids are bigger than amino acids, and MultiResidues can
             * have wild swings in energy on account of chemical perturbation.
             */
            double energyToPrune = (residue instanceof MultiResidue) ? multiResClashThreshold : clashThreshold;
            energyToPrune = (residue.getResidueType() == NA) ? energyToPrune * singletonNAPruningFactor * pruningFactor : energyToPrune;
            energyToPrune += minEnergy;

            for (int ri = 0; ri < nrot; ri++) {
                if (!check(i, ri) && (self(i, ri) > energyToPrune)) {
                    if (ri != 0 || !library.getUsingOrigCoordsRotamer()) {
                        eliminateRotamer(residues, i, ri, print);
                        logIfMaster(String.format(" Clash for rotamer %s %d: %16.8f >> %16.8f",
                                residue, ri, self(i, ri), minEnergy));
                    }
                }
            }
        }
    }

    /**
     * Prunes rotamer ri of residue i if all ri-j pair energies are worse than
     * the best i-j pair by some threshold value; additionally prunes ri-rj
     * pairs if they exceed the best i-j pair by a greater threshold value;
     * additionally performs this in reverse (searches over j-i).
     *
     * @param residues Residues whose rotamers are to be pruned.
     */
    public void prunePairClashes(Residue residues[]) {
        if (!prunePairClashes) {
            return;
        }
        int nResidues = residues.length;
        for (int i = 0; i < nResidues - 1; i++) {
            Residue resi = residues[i];
            Rotamer[] roti = resi.getRotamers(library);
            int ni = roti.length;
            for (int j = i + 1; j < nResidues; j++) {
                Residue resj = residues[j];
                Rotamer[] rotj = resj.getRotamers(library);
                int nj = rotj.length;
                double minPair = Double.MAX_VALUE;
                for (int ri = 0; ri < ni; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    for (int rj = 0; rj < nj; rj++) {
                        if (check(j, rj)) {
                            continue;
                        }
                        if (minPair > (pair(i, ri, j, rj) + self(i, ri) + self(j, rj))) {
                            minPair = (pair(i, ri, j, rj) + self(i, ri) + self(j, rj));
                        }
                    }
                }
                //double pruneThreshold = pairClashThreshold;
                double threshold = pairClashThreshold;
                if (resi instanceof MultiResidue) {
                    threshold += multiResPairClashAddn;
                }
                if (resj instanceof MultiResidue) {
                    threshold += multiResPairClashAddn;
                }

                int numNARes = (resi.getResidueType() == NA ? 1 : 0) + (resj.getResidueType() == NA ? 1 : 0);
                switch (numNARes) {
                    case 0:
                        break;
                    case 1:
                        threshold *= pairHalfPruningFactor;
                        //pruneThreshold *= pairHalfPruningFactor;
                        break;
                    case 2:
                        threshold *= pruningFactor;
                        //pruneThreshold *= pruningFactor;
                        break;
                    default:
                        throw new ArithmeticException(" RotamerOptimization.prunePairClashes() has somehow "
                                + "found less than zero or more than two nucleic acid residues in a pair of"
                                + " residues. This result should be impossible.");
                }
                threshold += minPair;

                // Check for elimination of any rotamer ri.
                for (int ri = 0; ri < ni; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    double eliminate = Double.MAX_VALUE;
                    for (int rj = 0; rj < nj; rj++) {
                        if (check(j, rj)) {
                            continue;
                        }
                        if ((pair(i, ri, j, rj) + self(i, ri) + self(j, rj)) < eliminate) {
                            eliminate = (pair(i, ri, j, rj) + self(i, ri) + self(j, rj));
                        }
                    }
                    // Prune based on clash threshold and the appropriate
                    // pruning factor (1.0 for AA pairs, pF for NA pairs,
                    // arithmetic mean for AA-NA pairs).
                    if (eliminate > threshold) {
                        // don't prune orig-coords rotamers
                        if (ri != 0 || !library.getUsingOrigCoordsRotamer()) {
                            eliminateRotamer(residues, i, ri, print);
                            logIfMaster(String.format(
                                    " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                                    resi, ri, resj, eliminate, threshold));
                        }
                    }
                    /*if (resi.getResidueType() == NA && resj.getResidueType() == NA) {
                     if (eliminate > (minPair + (pairClashThreshold * pruningFactor))) {
                     if (ri != 0 || !RotamerLibrary.getUsingOrigCoordsRotamer()) {
                     eliminateRotamer(residues, i, ri, print);
                     logIfMaster(String.format(
                     " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                     resi, ri, resj, eliminate, minPair + (pairClashThreshold * pruningFactor)));
                     }
                     }
                     } else if (resi.getResidueType() == NA || resj.getResidueType() == NA) {
                     if (eliminate > (minPair + (pairClashThreshold * pairHalfPruningFactor))) {
                     if (ri != 0 || !RotamerLibrary.getUsingOrigCoordsRotamer()) {
                     eliminateRotamer(residues, i, ri, print);
                     logIfMaster(String.format(
                     " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                     resi, ri, resj, eliminate, minPair + (pairClashThreshold * pairHalfPruningFactor)));
                     }
                     }
                     } else {
                     if (eliminate > (minPair + pairClashThreshold)) {
                     // don't prune orig-coords rotamers
                     if (ri != 0 || !RotamerLibrary.getUsingOrigCoordsRotamer()) {
                     eliminateRotamer(residues, i, ri, print);
                     logIfMaster(String.format(
                     " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                     resi, ri, resj, eliminate, minPair + pairClashThreshold));
                     }
                     }
                     }*/
                }
                // Check for elimination of any rotamer rj.
                for (int rj = 0; rj < nj; rj++) {
                    if (check(j, rj)) {
                        continue;
                    }
                    double eliminate = Double.MAX_VALUE;
                    for (int ri = 0; ri < ni; ri++) {
                        if (check(i, ri)) {
                            continue;
                        }
                        if ((pair(i, ri, j, rj) + self(i, ri) + self(j, rj)) < eliminate) {
                            eliminate = (pair(i, ri, j, rj) + self(i, ri) + self(j, rj));
                        }
                    }
                    if (eliminate > threshold) {
                        if (rj != 0 || !library.getUsingOrigCoordsRotamer()) {
                            eliminateRotamer(residues, j, rj, print);
                            logIfMaster(String.format(
                                    " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                                    resj, rj, resi, eliminate, threshold));
                        }
                    }
                    /*if (resi.getResidueType() == NA && resj.getResidueType() == NA) {
                     if (eliminate > (minPair + (pairClashThreshold * pruningFactor))) {
                     if (rj != 0 || !RotamerLibrary.getUsingOrigCoordsRotamer()) {
                     eliminateRotamer(residues, j, rj, print);
                     logIfMaster(String.format(
                     " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                     resj, rj, resi, eliminate, minPair + (pairClashThreshold * pruningFactor)));
                     }
                     }
                     } else if (resi.getResidueType() == NA || resj.getResidueType() == NA) {
                     if (eliminate > (minPair + (pairClashThreshold * pairHalfPruningFactor))) {
                     if (rj != 0 || !RotamerLibrary.getUsingOrigCoordsRotamer()) {
                     eliminateRotamer(residues, j, rj, print);
                     logIfMaster(String.format(
                     " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                     resj, rj, resi, eliminate, minPair + (pairClashThreshold * pairHalfPruningFactor)));
                     }
                     }
                     } else {
                     if (eliminate > (minPair + pairClashThreshold)) {
                     // don't prune orig-coords rotamers
                     if (rj != 0 || !RotamerLibrary.getUsingOrigCoordsRotamer()) {
                     eliminateRotamer(residues, j, rj, print);
                     logIfMaster(String.format(
                     " Pruning rotamer %s %d that clashes with all %s rotamers %16.8f >> %16.8f.",
                     resj, rj, resi, eliminate, minPair + pairClashThreshold));
                     }
                     }
                     }*/
                }
                /* Presently consitutively false
                 if (pruneIndivPairs) {
                 for (int ri = 0; ri < ni; ni++) {
                 if (check(i, ri)) {
                 continue;
                 }
                 for (int rj = 0; rj < nj; nj++) {
                 if (check(j, rj) || check(i, ri, j, rj)) {
                 continue;
                 }
                 double currentPairEnergy = pair(i, ri, j, rj);
                 if (currentPairEnergy > indivPruneThreshold) {
                 eliminateRotamerPair(residues, i, ri, j, ri, print);
                 logIfMaster(String.format(" Pruning rotamer pair %s %d %s %d that clashes with best "
                 + "%s %s rotamer pair %16.8f >> %16.8f.", resi, ri, resj, rj, resi, resj,
                 currentPairEnergy, indivPruneThreshold));
                 }
                 }
                 }
                 }*/
            }
        }
    }

    private class SinglesEnergyRegion extends WorkerRegion {

        private final SinglesEnergyLoop energyLoop;
        private final Residue residues[];
        private final int nResidues;
        private final boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        //private double originalAtomicCoordinates[][][];

        public SinglesEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new SinglesEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;

        }

        @Override
        public void start() {
            // setup and compute E(BB)
            /*ArrayList<Residue> currentResidueList = new ArrayList<>();
             currentResidueList.addAll(Arrays.asList(residues));
             originalAtomicCoordinates = storeCoordinates(currentResidueList);*/
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer rotamers[] = residue.getRotamers(library);
                switch (residue.getResidueType()) {
                    case NA:
                        break;
                    case AA:
                    default:
//                        logIfMaster(String.format(" Applying to residue %s, rotamer %s.", residue.toString(), rotamers[0].toString()));
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                        break;
                }
                turnOffAtoms(residue);
            }
            try {
                backboneEnergy = currentEnergy(new ArrayList<>(), false);
            } catch (ArithmeticException ex) {
                logger.severe(String.format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
            }
            logIfMaster(String.format(" Backbone energy:  %16.8f\n", backboneEnergy));
        }

        @Override
        public void run() throws Exception {
            if (!jobMapSingles.isEmpty()) {
//                Integer jobKeysArray[] = jobMapSingles.keySet().toArray(new Integer[1]);
//                int maxJobKey = Integer.MIN_VALUE;
//                for (int i = 0; i < jobKeysArray.length; i++) {
//                    maxJobKey = (jobKeysArray[i] > maxJobKey) ? jobKeysArray[i] : maxJobKey;
//                }
//                execute(0, maxJobKey, energyLoop);
                execute(jobMapSingles.keySet().iterator(), energyLoop);
            }
        }

        @Override
        public void finish() {
            // broadcast the "I'm finished" signal
            double finished[] = new double[3];
            for (int x = 0; x < finished.length; x++) {
                finished[x] = -1.0;
            }
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            // wait for everyone else to send in their self energies
            int waiting = 0;
            while (!selfsDone) {
                try {
                    Thread.sleep(POLLING_FREQUENCY);
                    if (waiting++ == 1000) {
                        logger.warning(String.format("Process %d experiencing long wait for others' trimer energies.", Comm.world().rank()));
                    }
                } catch (InterruptedException ex) {
                }
            }

            // Prune clashes for all singles (not just the ones this node did).
            if (pruneClashes) {
                pruneSingleClashes(residues);
            }

            // Print what we've got so far.
            if (master && verbose) {
                for (int i = 0; i < residues.length; i++) {
                    Residue residue = residues[i];
                    Rotamer rotamers[] = residue.getRotamers(library);
                    for (int ri = 0; ri < rotamers.length; ri++) {
                        logger.info(String.format(" (recap) Self energy %7s %-2d: %16.8f", residues[i], ri, self(i, ri)));
                    }
                }
            }
        }

        private class SinglesEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                // Compute the self-energy for each rotamer
                // Rotamer-level job indexing method
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!jobMapSingles.keySet().contains(jobKey)) {
                        //logger.warning(String.format("(sbox %d) Unexpected: singles jobKey not contained in map (key = %d).", BOXNUM, jobKey));
                        //logger.warning(String.format("Unexpected: singles jobKey not contained in map (key = %d).", jobKey));
                        continue;
                    }
                    Integer job[] = jobMapSingles.get(jobKey);
                    int i = job[0];
                    int ri = job[1];
                    if (!(useOrigCoordsRot && ri == 0) && pruneClashes && check(i, ri)) {
                        continue;
                    }
                    Residue resi = residues[i];
                    Rotamer roti = resi.getRotamers(library)[ri];
                    //double[][] resiOriginalCoordinates = (resi.getResidueType() == NA ? storeSingleCoordinates(resi) : null);
                    ResidueState resiOriginal = (resi.getResidueType() == NA ? resi.storeState() : null);
                    turnOnAtoms(resi);
                    RotamerLibrary.applyRotamer(resi, roti);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                    if (writeVideo) {
                        videoWriter.write(String.format("%03d-%03d", i, ri));
//                        if (videoWriter.skipEnergies) {
//                            videoWriter.broadcastDummy(i, ri);
//                            continue;
//                        }
                    }

                    double selfEnergy;
                    if (writeVideo || skipEnergies) {
                        selfEnergy = 0;
                    } else {
                        List<Residue> rList = Arrays.asList(resi);
                        selfEnergy = currentEnergy(rList) - backboneEnergy;
                    }

                    double anEnergy[] = new double[3];
                    anEnergy[0] = i;
                    anEnergy[1] = ri;
                    anEnergy[2] = selfEnergy;
                    DoubleBuf anEnergyBuf = DoubleBuf.buffer(anEnergy);

                    multicastBuf(anEnergyBuf);
                    logger.info(String.format(" Self %7s %-2d: %16.8f", resi, ri, selfEnergy));
                    if (resi.getResidueType() == NA) {
                        //revertSingleResidueCoordinates(resi, resiOriginalCoordinates);
                        resi.revertState(resiOriginal);
                    } else {
                        RotamerLibrary.applyRotamer(resi, resi.getRotamers(library)[0]);
                    }
                    turnOffAtoms(resi);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                }
            }
        }
    }

    public void multicastBuf(Buf message) {
        for (int p = 0; p < numProc; p++) {
            try {
                world.send(p, message);
            } catch (IOException ex) {
                logger.log(Level.WARNING, ex.getMessage(), ex);
            }
        }
    }

    private class PairsEnergyRegion extends WorkerRegion {

        private final PairsEnergyLoop energyLoop;
        private final Residue residues[];
        private final int nResidues;
        private final boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        //private double originalAtomicCoordinates[][][];

        public PairsEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new PairsEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
        }

        @Override
        public void run() throws Exception {
            if (!jobMapPairs.isEmpty()) {
                execute(jobMapPairs.keySet().iterator(), energyLoop);
            }
        }

        @Override
        public void finish() {
            // broadcast the "I'm finished" signal
            double finished[] = new double[5];
            for (int x = 0; x < finished.length; x++) {
                finished[x] = -1.0;
            }
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            // wait for everyone else to send in their pair energies
            int waiting = 0;
            while (!pairsDone) {
                try {
                    Thread.sleep(POLLING_FREQUENCY);
                    if (waiting++ == 1000) {
                        logger.warning(String.format("Process %d experiencing long wait for others' pair energies.", Comm.world().rank()));
                    }
                } catch (InterruptedException ex) {
                }
            }

            // Prune each rotamer that clashes with all rotamers from a 2nd residue.
            if (prunePairClashes) {
                prunePairClashes(residues);
            }

            // Print what we've got so far.
            if (master && verbose) {
                for (int i = 0; i < nResidues; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (pruneClashes && check(i, ri) && !(useOrigCoordsRot && ri == 0)) {
                            continue;
                        }
                        for (int j = i + 1; j < nResidues; j++) {
                            Residue resj = residues[j];
                            Rotamer rotj[] = resj.getRotamers(library);
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if (pruneClashes && (check(j, rj) || check(i, ri, j, rj))
                                        && !(useOrigCoordsRot && ri == 0 && rj == 0)) {
                                    continue;
                                }
                                logger.info(String.format(" (recap) Pair energy %7s %-2d, %7s %-2d: %16.8f", residues[i], ri, residues[j], rj, pair(i, ri, j, rj)));
                            }
                        }
                    }
                }
            }
        }

        private class PairsEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                // Compute the pair-energy for each pair of rotamers
                // Pair-level job indexing method
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!jobMapPairs.keySet().contains(jobKey)) {
                        //logger.warning(String.format("(sbox %d) Unexpected: pairs jobKey not contained in map (key = %d).", BOXNUM, jobKey));
                        //logger.warning(String.format("Unexpected: pairs jobKey not contained in map (key = %d).", jobKey));
                        continue;
                    }
                    Integer job[] = jobMapPairs.get(jobKey);
                    int i = job[0];
                    int ri = job[1];
                    int j = job[2];
                    int rj = job[3];
                    if (!(useOrigCoordsRot && ri == 0 && rj == 0)
                            && pruneClashes && (check(i, ri) || check(j, rj) || check(i, ri, j, rj))) {
                        continue;
                    }
                    Residue resi = residues[i];
                    Rotamer roti = resi.getRotamers(library)[ri];
                    Residue resj = residues[j];
                    Rotamer rotj = resj.getRotamers(library)[rj];
                    /*double[][] resiOriginalCoordinates = (resi.getResidueType() == NA ? storeSingleCoordinates(resi) : null);
                     double[][] resjOriginalCoordinates = (resj.getResidueType() == NA ? storeSingleCoordinates(resj) : null);*/
                    ResidueState resiOriginal = resi.getResidueType() == NA ? resi.storeState() : null;
                    ResidueState resjOriginal = resj.getResidueType() == NA ? resj.storeState() : null;

                    // turn on, apply rot
                    turnOnAtoms(resi);
                    turnOnAtoms(resj);
                    RotamerLibrary.applyRotamer(resi, roti);
                    RotamerLibrary.applyRotamer(resj, rotj);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                    if (writeVideo) {
                        videoWriter.write(String.format("%03d-%03d_%03d-%03d", i, ri, j, rj));
//                        if (videoWriter.skipEnergies) {
//                            videoWriter.broadcastDummy(i, ri, j, rj);
//                            continue;
//                        }
                    }
                    double twoBodyEnergy;

                    List<Residue> rList = Arrays.asList(new Residue[]{resi, resj});
                    if (writeVideo || skipEnergies) {
                        twoBodyEnergy = 0;
                    } else if (distanceMatrix != null) {
                        int indexI = allResiduesList.indexOf(resi);
                        int indexJ = allResiduesList.indexOf(resj);
                        double dist = checkDistanceMatrix(indexI, ri, indexJ, rj);
                        if (dist < superpositionThreshold) {
                            twoBodyEnergy = 1.0E100;
                            logger.info(String.format(" Pair %7s %-2d, %7s %-2d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang", resi, ri, resj, rj, dist, superpositionThreshold));
                        } else {
                            twoBodyEnergy = currentEnergy(rList) - self(i, ri) - self(j, rj) - backboneEnergy;
                            String distString = (dist < Double.MAX_VALUE) ? String.format("%10.3f", dist) : String.format("     large");
                            logger.info(String.format(" Pair %7s %-2d, %7s %-2d: %16.8f at %s Ang", resi, ri, resj, rj, twoBodyEnergy, distString));
                        }
                    } else {
                        twoBodyEnergy = currentEnergy(rList)
                                - self(i, ri) - self(j, rj) - backboneEnergy;
                        logger.info(String.format(" Pair %7s %-2d, %7s %-2d: %16.8f", resi, ri, resj, rj, twoBodyEnergy));
                    }

                    // get energy and broadcast it
                    double commEnergy[] = new double[5];
                    commEnergy[0] = i;
                    commEnergy[1] = ri;
                    commEnergy[2] = j;
                    commEnergy[3] = rj;
                    commEnergy[4] = twoBodyEnergy;
                    DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
                    multicastBuf(commEnergyBuf);

                    // move back, turn off
                    if (resi.getResidueType() == NA) {
                        //revertSingleResidueCoordinates(resi, resiOriginalCoordinates);
                        resi.revertState(resiOriginal);
                    } else {
                        RotamerLibrary.applyRotamer(resi, resi.getRotamers(library)[0]);
                    }
                    if (resj.getResidueType() == NA) {
                        //revertSingleResidueCoordinates(resj, resjOriginalCoordinates);
                        resj.revertState(resjOriginal);
                    } else {
                        RotamerLibrary.applyRotamer(resj, resj.getRotamers(library)[0]);
                    }
                    turnOffAtoms(resi);
                    turnOffAtoms(resj);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                }
            }
        }
    }

    private class TriplesEnergyRegion extends WorkerRegion {

        private final TriplesEnergyLoop energyLoop;
        private final Residue residues[];
        private final int nResidues;
        private final boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        //private double originalAtomicCoordinates[][][];
        private double localDistanceMatrix[][][][];

        public TriplesEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new TriplesEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
            localDistanceMatrix = new double[nResidues][][][];
        }

        @Override
        public void start() {
            if (distance <= 0) {
                logger.info(" Calculating local distance matrix using non-eliminated rotamers.");
                // TODO: check on the location of this call - might need to be done per trimer-job
                localSequentialDistanceMatrix(residues, localDistanceMatrix);
            }
        }

        @Override
        public void run() throws Exception {
            if (!jobMapTrimers.isEmpty()) {
                execute(jobMapTrimers.keySet().iterator(), energyLoop);
            }
        }

        @Override
        public void finish() {
            // broadcast the "I'm finished" signal
            double finished[] = new double[7];
            for (int x = 0; x < finished.length; x++) {
                finished[x] = -1.0;
            }
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            // wait for everyone else to send in their trimer energies
            int waiting = 0;
            while (!trimersDone) {
                try {
                    Thread.sleep(POLLING_FREQUENCY);
                    if (waiting++ == 1000) {
                        logger.warning(String.format("Process %d experiencing long wait for others' trimer energies.", Comm.world().rank()));
                    }
                } catch (InterruptedException ex) {
                }
            }

            // Print what we've got so far.
            if (master && verbose) {
                for (int i = 0; i < nResidues; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (pruneClashes && check(i, ri) && !(useOrigCoordsRot && ri == 0)) {
                            continue;
                        }
                        for (int j = i + 1; j < nResidues; j++) {
                            Residue resj = residues[j];
                            Rotamer rotj[] = resj.getRotamers(library);
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if (pruneClashes && (check(j, rj) || check(i, ri, j, rj))
                                        && !(useOrigCoordsRot && ri == 0 && rj == 0)) {
                                    continue;
                                }
                                for (int k = j + 1; k < nResidues; k++) {
                                    Residue resk = residues[k];
                                    Rotamer rotk[] = resk.getRotamers(library);
                                    for (int rk = 0; rk < rotk.length; rk++) {
                                        if (pruneClashes && (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk))
                                                && !(useOrigCoordsRot && ri == 0 && rj == 0 && rk == 0)) {
                                            continue;
                                        }
                                        logger.info(String.format(" (recap) Trimer energy %7s %-2d, %7s %-2d, %7s %-2d: %16.8f", resi, ri, resj, rj, resk, rk, triple(i, ri, j, rj, k, rk)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        private class TriplesEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                // Trimer-level job indexing method
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!jobMapTrimers.keySet().contains(jobKey)) {
                        //logger.warning(String.format("(sbox %d) Unexpected: trimers jobKey not contained in map (key = %d).", BOXNUM, jobKey));
                        //logger.warning(String.format("Unexpected: trimers jobKey not contained in map (key = %d).", jobKey));
                        continue;
                    }
                    Integer job[] = jobMapTrimers.get(jobKey);
                    int i = job[0];
                    int ri = job[1];
                    int j = job[2];
                    int rj = job[3];
                    int k = job[4];
                    int rk = job[5];
                    if (!(useOrigCoordsRot && ri == 0 && rj == 0 && rk == 0)
                            && pruneClashes && (check(i, ri) || check(j, rj) || check(k, rk) || check(i, ri, j, rj) || check(i, ri, k, rk) || check(j, rj, k, rk))) {
                        continue;
                    }
                    Residue resi = residues[i];
                    Rotamer roti = resi.getRotamers(library)[ri];
                    Residue resj = residues[j];
                    Rotamer rotj = resj.getRotamers(library)[rj];
                    Residue resk = residues[k];
                    Rotamer rotk = resk.getRotamers(library)[rk];
                    /*double[][] resiOriginalCoordinates = (resi.getResidueType() == NA ? storeSingleCoordinates(resi) : null);
                     double[][] resjOriginalCoordinates = (resj.getResidueType() == NA ? storeSingleCoordinates(resj) : null);
                     double[][] reskOriginalCoordinates = (resk.getResidueType() == NA ? storeSingleCoordinates(resk) : null);*/
                    ResidueState resiOriginal = resi.getResidueType() == NA ? resi.storeState() : null;
                    ResidueState resjOriginal = resj.getResidueType() == NA ? resj.storeState() : null;
                    ResidueState reskOriginal = resk.getResidueType() == NA ? resk.storeState() : null;

                    int indexOfI = allResiduesList.indexOf(resi);
                    int indexOfJ = allResiduesList.indexOf(resj);
                    int indexOfK = allResiduesList.indexOf(resk);
                    double dij, dik, djk;
                    if (distance > 0) {
                        // Distance matrix is asymmetric, but in present implementation i < j < k.
                        dij = checkDistanceMatrix(indexOfI, ri, indexOfJ, rj);
                        dik = checkDistanceMatrix(indexOfI, ri, indexOfK, rk);
                        djk = checkDistanceMatrix(indexOfJ, rj, indexOfK, rk);
                    } else {
                        dij = localDistanceMatrix[i][ri][j][rj];
                        dik = localDistanceMatrix[i][ri][k][rk];
                        djk = localDistanceMatrix[j][rj][k][rk];
                    }
                    double dist = Math.min(dij, Math.min(dik, djk));
                    double threeBodyEnergy;
                    List<Residue> rList = Arrays.asList(new Residue[]{resi, resj, resk});
                    if (!threeBodyCutoff || (dist < threeBodyCutoffDist)) {
                        if (dist < superpositionThreshold) {
                            threeBodyEnergy = 1.0E100;
                            logger.info(String.format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang.",
                                    resi, ri, resj, rj, resk, rk, dist, superpositionThreshold));
                        } else {
                            // turn on, apply rotamers
                            turnOnAtoms(resi);
                            turnOnAtoms(resj);
                            turnOnAtoms(resk);
                            RotamerLibrary.applyRotamer(resi, roti);
                            RotamerLibrary.applyRotamer(resj, rotj);
                            RotamerLibrary.applyRotamer(resk, rotk);
                            if (algorithmListener != null) {
                                algorithmListener.algorithmUpdate(molecularAssembly);
                            }
                            if (writeVideo) {
                                videoWriter.write(String.format("%03d-%03d_%03d-%03d_%03d-%03d", i, ri, j, rj, k, rk));
//                                if (videoWriter.skipEnergies) {
//                                    videoWriter.broadcastDummy(i, ri, j, rj, k, rk);
//                                    continue;
//                                }
                            }

                            // compute trimer and broadcast it
                            if (writeVideo || skipEnergies) {
                                threeBodyEnergy = 0;
//                                videoWriter.broadcastDummy(i, ri, j, rj, k, rk);
                            } else {
                                threeBodyEnergy = currentEnergy(rList)
                                        - self(i, ri) - self(j, rj) - self(k, rk)
                                        - pair(i, ri, j, rj) - pair(i, ri, k, rk)
                                        - pair(j, rj, k, rk) - backboneEnergy;
                                String distString = (dist < Double.MAX_VALUE) ? String.format("%10.3f", dist) : String.format("     large");
                                logger.info(String.format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %s Ang.",
                                        resi, ri, resj, rj, resk, rk, threeBodyEnergy, distString));
                            }

                            // revert rotamers, turn off
                            if (resi.getResidueType() == NA) {
                                //revertSingleResidueCoordinates(resi, resiOriginalCoordinates);
                                resi.revertState(resiOriginal);
                            } else {
                                RotamerLibrary.applyRotamer(resi, resi.getRotamers(library)[0]);
                            }
                            if (resj.getResidueType() == NA) {
                                //revertSingleResidueCoordinates(resj, resjOriginalCoordinates);
                                resj.revertState(resjOriginal);
                            } else {
                                RotamerLibrary.applyRotamer(resj, resj.getRotamers(library)[0]);
                            }
                            if (resk.getResidueType() == NA) {
                                //revertSingleResidueCoordinates(resk, reskOriginalCoordinates);
                                resk.revertState(reskOriginal);
                            } else {
                                RotamerLibrary.applyRotamer(resk, resk.getRotamers(library)[0]);
                            }
                            turnOffAtoms(resi);
                            turnOffAtoms(resj);
                            turnOffAtoms(resk);
                        }
                        double commEnergy[] = new double[7];
                        commEnergy[0] = i;
                        commEnergy[1] = ri;
                        commEnergy[2] = j;
                        commEnergy[3] = rj;
                        commEnergy[4] = k;
                        commEnergy[5] = rk;
                        commEnergy[6] = threeBodyEnergy;
                        DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
                        multicastBuf(commEnergyBuf);
                        if (algorithmListener != null) {
                            algorithmListener.algorithmUpdate(molecularAssembly);
                        }
                    } else {
                        logger.info(String.format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %10.3f Ang.",
                                resi, ri, resj, rj, resk, rk, 0.0, threeBodyCutoffDist));
                    }
                }
            }
        }
    }

    private class QuadsEnergyRegion extends WorkerRegion {

        private final QuadsEnergyLoop energyLoop;
        private final Residue residues[];
        private final int nResidues;
        private final boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        //private double originalAtomicCoordinates[][][];
        private double localDistanceMatrix[][][][];

        public QuadsEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new QuadsEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
            localDistanceMatrix = new double[nResidues][][][];
        }

        @Override
        public void start() {
            if (distance <= 0) {
                logger.info(" Calculating local distance matrix using non-eliminated rotamers.");
                // TODO: check on the location of this call - might need to be done per quad-job
                localSequentialDistanceMatrix(residues, localDistanceMatrix);
            }
        }

        @Override
        public void run() throws Exception {
            if (!jobMapQuads.isEmpty()) {
                execute(jobMapQuads.keySet().iterator(), energyLoop);
            }
        }

        @Override
        public void finish() {
            // no "I'm finished" signal for quads
        }

        private class QuadsEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                // Quad-level job indexing method
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!jobMapQuads.keySet().contains(jobKey)) {
                        continue;
                    }
                    Integer job[] = jobMapQuads.get(jobKey);
                    int i = job[0];
                    int ri = job[1];
                    int j = job[2];
                    int rj = job[3];
                    int k = job[4];
                    int rk = job[5];
                    int l = job[6];
                    int rl = job[7];
                    if (!(useOrigCoordsRot && ri == 0 && rj == 0 && rk == 0 && rl == 0)
                            && pruneClashes
                            && (check(i, ri) || check(j, rj) || check(k, rk) || check(l, rl)
                            || check(i, ri, j, rj) || check(i, ri, k, rk) || check(j, rj, k, rk)
                            || check(i, ri, l, rl) || check(j, rj, l, rl) || check(k, rk, l, rl))) {
                        continue;
                    }
                    Residue resi = residues[i];
                    Rotamer roti = resi.getRotamers(library)[ri];
                    Residue resj = residues[j];
                    Rotamer rotj = resj.getRotamers(library)[rj];
                    Residue resk = residues[k];
                    Rotamer rotk = resk.getRotamers(library)[rk];
                    Residue resl = residues[l];
                    Rotamer rotl = resl.getRotamers(library)[rl];
                    ResidueState resiOriginalCoordinates = (resi.getResidueType() == NA ? resi.storeState() : null);
                    ResidueState resjOriginalCoordinates = (resj.getResidueType() == NA ? resj.storeState() : null);
                    ResidueState reskOriginalCoordinates = (resk.getResidueType() == NA ? resk.storeState() : null);
                    ResidueState reslOriginalCoordinates = (resl.getResidueType() == NA ? resl.storeState() : null);
                    /*double[][] resiOriginalCoordinates = (resi.getResidueType() == NA ? storeSingleCoordinates(resi) : null);
                     double[][] resjOriginalCoordinates = (resj.getResidueType() == NA ? storeSingleCoordinates(resj) : null);
                     double[][] reskOriginalCoordinates = (resk.getResidueType() == NA ? storeSingleCoordinates(resk) : null);
                     double[][] reslOriginalCoordinates = (resl.getResidueType() == NA ? storeSingleCoordinates(resl) : null);*/

                    int indexOfI = allResiduesList.indexOf(resi);
                    int indexOfJ = allResiduesList.indexOf(resj);
                    int indexOfK = allResiduesList.indexOf(resk);
                    int indexOfL = allResiduesList.indexOf(resl);
                    double dij, dik, djk, dil, djl, dkl;
                    List<Residue> rList = Arrays.asList(new Residue[]{resi, resj, resk, resl});
                    if (distance > 0) {
                        // Distance matrix is asymmetric, but in present implementation i < j < k.
                        dij = checkDistanceMatrix(indexOfI, ri, indexOfJ, rj);
                        dik = checkDistanceMatrix(indexOfI, ri, indexOfK, rk);
                        djk = checkDistanceMatrix(indexOfJ, rj, indexOfK, rk);
                        dil = checkDistanceMatrix(indexOfI, ri, indexOfL, rl);
                        djl = checkDistanceMatrix(indexOfJ, rj, indexOfL, rl);
                        dkl = checkDistanceMatrix(indexOfK, rk, indexOfL, rl);
                    } else {
                        dij = localDistanceMatrix[i][ri][j][rj];
                        dik = localDistanceMatrix[i][ri][k][rk];
                        djk = localDistanceMatrix[j][rj][k][rk];
                        dil = localDistanceMatrix[i][ri][l][rl];
                        djl = localDistanceMatrix[j][rj][l][rl];
                        dkl = localDistanceMatrix[k][rk][l][rl];
                    }
                    double dist = Math.min(dkl, Math.min(djl, Math.min(dil, Math.min(dij, Math.min(dik, djk)))));
                    double quadEnergy;
                    if (!quadCutoff || (dist < quadCutoffDist)) {
                        if (dist < superpositionThreshold) {
                            quadEnergy = 1.0E100;
                            logger.info(String.format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang.",
                                    resi, ri, resj, rj, resk, rk, resl, rl, dist, superpositionThreshold));
                        } else {
                            // turn on, apply rotamers
                            turnOnAtoms(resi);
                            turnOnAtoms(resj);
                            turnOnAtoms(resk);
                            turnOnAtoms(resl);
                            RotamerLibrary.applyRotamer(resi, roti);
                            RotamerLibrary.applyRotamer(resj, rotj);
                            RotamerLibrary.applyRotamer(resk, rotk);
                            RotamerLibrary.applyRotamer(resl, rotl);
                            if (algorithmListener != null) {
                                algorithmListener.algorithmUpdate(molecularAssembly);
                            }

                            // compute quad and broadcast it
                            if (writeVideo || skipEnergies) {
                                quadEnergy = 0;
                            } else {
                                quadEnergy = currentEnergy(rList)
                                        - self(i, ri) - self(j, rj) - self(k, rk) - self(l, rl)
                                        - pair(i, ri, j, rj) - pair(i, ri, k, rk) - pair(i, ri, l, rl)
                                        - pair(j, rj, k, rk) - pair(j, rj, l, rl) - pair(k, rk, l, rl)
                                        - triple(i, ri, j, rj, k, rk) - triple(i, ri, j, rj, l, rl) - triple(i, ri, k, rk, l, rl) - triple(j, rj, k, rk, l, rl)
                                        - backboneEnergy;
                                String distString = (dist < Double.MAX_VALUE) ? String.format("%10.3f", dist) : String.format("     large");
                                if (Math.abs(quadEnergy) > 1.0) {
                                    StringBuilder sb = new StringBuilder();
                                    sb.append(String.format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %s Ang.\n",
                                            resi, ri, resj, rj, resk, rk, resl, rl, quadEnergy, distString));
                                    sb.append(String.format("   Explain: (ref %d) \n", jobKey));
                                    sb.append(String.format("     Self %3d %3d:                  %.3f\n", i, ri, self(i, ri)));
                                    sb.append(String.format("     Self %3d %3d:                  %.3f\n", j, rj, self(j, rj)));
                                    sb.append(String.format("     Self %3d %3d:                  %.3f\n", k, rk, self(k, rk)));
                                    sb.append(String.format("     Self %3d %3d:                  %.3f\n", l, rl, self(l, rl)));
                                    sb.append(String.format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, j, rj, pair(i, ri, j, rj)));
                                    sb.append(String.format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, k, rk, pair(i, ri, k, rk)));
                                    sb.append(String.format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, l, rl, pair(i, ri, l, rl)));
                                    sb.append(String.format("     Pair %3d %3d %3d %3d:          %.3f\n", j, rj, k, rk, pair(j, rj, k, rk)));
                                    sb.append(String.format("     Pair %3d %3d %3d %3d:          %.3f\n", j, rj, l, rl, pair(j, rj, l, rl)));
                                    sb.append(String.format("     Pair %3d %3d %3d %3d:          %.3f\n", k, rk, l, rl, pair(k, rk, l, rl)));
                                    sb.append(String.format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, j, rj, k, rk, triple(i, ri, j, rj, k, rk)));
                                    sb.append(String.format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, j, rj, l, rl, triple(i, ri, j, rj, l, rl)));
                                    sb.append(String.format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, k, rk, l, rl, triple(i, ri, k, rk, l, rl)));
                                    sb.append(String.format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", j, rj, k, rk, l, rl, triple(j, rj, k, rk, l, rl)));
                                    sb.append(String.format("     backbone:                      %.3f\n", backboneEnergy));
                                    sb.append(String.format("     quadEnergy:                 %.3f\n", quadEnergy));
                                    sb.append(String.format("     --s--\n"));
                                    sb.append(String.format("     Active residues:\n"));
                                    for (int debug = 0; debug < residues.length; debug++) {
                                        if (residues[debug].getSideChainAtoms().get(0).getUse()) {
                                            sb.append(String.format("       %s\n", residues[debug].toString()));
                                        }
                                    }
                                    sb.append(String.format("     --f--\n"));
                                    logger.info(sb.toString());
//                                    String filename = FilenameUtils.removeExtension(molecularAssembly.getFile().toString()) + "." + jobKey;
//                                    File file = new File(filename);
//                                    PDBFilter filt = new PDBFilter(file, molecularAssembly, null, null);
//                                    filt.writeFile(file, false);
                                } else {
                                    logger.info(String.format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %s Ang.",
                                            resi, ri, resj, rj, resk, rk, resl, rl, quadEnergy, distString));
                                }
                            }

                            // revert rotamers, turn off
                            if (resi.getResidueType() == NA) {
                                resi.revertState(resiOriginalCoordinates);
                                //revertSingleResidueCoordinates(resi, resiOriginalCoordinates);
                            } else {
                                RotamerLibrary.applyRotamer(resi, resi.getRotamers(library)[0]);
                            }
                            if (resj.getResidueType() == NA) {
                                resj.revertState(resjOriginalCoordinates);
                                //revertSingleResidueCoordinates(resj, resjOriginalCoordinates);
                            } else {
                                RotamerLibrary.applyRotamer(resj, resj.getRotamers(library)[0]);
                            }
                            if (resk.getResidueType() == NA) {
                                resk.revertState(reskOriginalCoordinates);
                                //revertSingleResidueCoordinates(resk, reskOriginalCoordinates);
                            } else {
                                RotamerLibrary.applyRotamer(resk, resk.getRotamers(library)[0]);
                            }
                            if (resl.getResidueType() == NA) {
                                resl.revertState(reslOriginalCoordinates);
                                //revertSingleResidueCoordinates(resl, reslOriginalCoordinates);
                            } else {
                                RotamerLibrary.applyRotamer(resl, resl.getRotamers(library)[0]);
                            }
                            turnOffAtoms(resi);
                            turnOffAtoms(resj);
                            turnOffAtoms(resk);
                            turnOffAtoms(resl);
                        }
                        /* No broadcasting of quad energies, just print!
                         double commEnergy[] = new double[9];
                         commEnergy[0] = i;
                         commEnergy[1] = ri;
                         commEnergy[2] = j;
                         commEnergy[3] = rj;
                         commEnergy[4] = k;
                         commEnergy[5] = rk;
                         commEnergy[6] = l;
                         commEnergy[7] = rl;
                         commEnergy[8] = quadEnergy;
                         DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
                         multicastBuf(commEnergyBuf);
                         if (algorithmListener != null) {
                         algorithmListener.algorithmUpdate(molecularAssembly);
                         }
                         */
                    } else {
                        logger.info(String.format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %10.3f Ang.",
                                resi, ri, resj, rj, resk, rk, resl, rl, 0.0, quadCutoffDist));
                    }
                }
            }
        }
    }

    private class DistanceRegion extends ParallelRegion {

        private final DistanceLoop distanceLoops[];
        private final int nResidues;
        private final Crystal crystal;
        private final int nSymm;
        private final int lists[][][];
        private final IntegerSchedule pairwiseSchedule;

        public DistanceRegion(int nt, int nResidues, Crystal crystal,
                int lists[][][], IntegerSchedule schedule) {
            distanceLoops = new DistanceLoop[nt];
            this.nResidues = nResidues;
            this.crystal = crystal;
            this.nSymm = crystal.spaceGroup.getNumberOfSymOps();
            this.lists = lists;
            for (int i = 0; i < nt; i++) {
                distanceLoops[i] = new DistanceLoop();
            }
            pairwiseSchedule = schedule;
        }

        @Override
        public void run() throws Exception {
            try {
                int threadID = getThreadIndex();
                execute(0, nResidues - 1, distanceLoops[threadID]);
            } catch (Exception e) {
                String message = " Exception computing residue-residue distances.";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class DistanceLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return pairwiseSchedule;
            }

            @Override
            public void run(int lb, int ub) {
                // Loop over symmetry operators.
                for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                    SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
                    // Loop over residues.
                    for (int i = lb; i <= ub; i++) {
                        Residue residuei = allResiduesArray[i];
                        Rotamer rotamersi[] = residuei.getRotamers(library);
                        int lengthRi;
                        boolean forcedResidueI = false;
                        try {
                            if (checkIfForced(residuei)) {
                                forcedResidueI = true;
                                lengthRi = 1;
                            } else {
                                lengthRi = rotamersi.length;
                            }
                        } catch (IndexOutOfBoundsException ex) {
                            logger.warning(String.format(" Exception in distance loop: %s", ex.toString()));
                            continue;
                        }
                        int list[] = lists[iSymOp][i];
                        int nList = list.length;

                        // Loop over Residue i's rotamers
                        for (int ri = 0; ri < lengthRi; ri++) {
                            double xi[][];
                            synchronized (residuei) {
                                if (!forcedResidueI) {
                                    Rotamer rotameri = rotamersi[ri];
                                    RotamerLibrary.applyRotamer(residuei, rotameri);
                                    //xi = storeSingleCoordinates(residuei);
                                    xi = residuei.storeCoordinateArray();
                                } else {
                                    //xi = storeSingleCoordinates(residuei);
                                    xi = residuei.storeCoordinateArray();
                                }
                            }
                            // Loop over Residue i's neighbors.
                            for (int k = 0; k < nList; k++) {
                                int j = list[k];
                                if (i == j) {
                                    continue;
                                }
                                Residue residuej = allResiduesArray[j];
                                Rotamer rotamersj[] = residuej.getRotamers(library);
                                int lengthRj;
                                boolean forcedResidueJ = false;
                                try {
                                    if (checkIfForced(residuej)) {
                                        forcedResidueJ = true;
                                        lengthRj = 1;
                                    } else {
                                        lengthRj = rotamersj.length;
                                    }
                                } catch (IndexOutOfBoundsException ex) {
                                    continue;
                                }
                                // Loop over the neighbor's rotamers
                                for (int rj = 0; rj < lengthRj; rj++) {
                                    double xj[][];
                                    synchronized (residuej) {
                                        if (!forcedResidueJ) {
                                            Rotamer rotamerj = rotamersj[rj];
                                            RotamerLibrary.applyRotamer(residuej, rotamerj);
                                            //xj = storeSingleCoordinates(residuej);
                                            xj = residuej.storeCoordinateArray();
                                        } else {
                                            //xj = storeSingleCoordinates(residuej);
                                            xj = residuej.storeCoordinateArray();
                                        }
                                    }

                                    if (getThreadIndex() == 0 && algorithmListener != null) {
                                        algorithmListener.algorithmUpdate(molecularAssembly);
                                    }

                                    double r = interResidueDistance(xi, xj, symOp);
                                    //logger.info(String.format(" %s %d %s %d %16.8f", residuei.toString(), i, residuej.toString(), j, r));
                                    if (i < j) {
                                        if (r < distanceMatrix[i][ri][j][rj]) {
                                            distanceMatrix[i][ri][j][rj] = r;
                                            // Would apply symmetry to the distance matrix:
                                            // distanceMatrix[j][rj][i][ri] = r;
                                        }
                                    } else if (r < distanceMatrix[j][rj][i][ri]) {
                                        distanceMatrix[j][rj][i][ri] = r;
                                        // Would apply symmetry to the distance matrix:
                                        // distanceMatrix[j][rj][i][ri] = r;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private class VideoWriter {

        private final MolecularAssembly molecularAssembly;
        private final PDBFilter filter;
        private boolean ignoreInactiveAtoms = false;
        private final String videoDir;
        private int snapshotNum = 1;
        public final boolean skipEnergies;

        public VideoWriter(MolecularAssembly molecularAssembly, boolean ignoreInactiveAtoms, boolean skipEnergies) {
            this.molecularAssembly = molecularAssembly;
            this.ignoreInactiveAtoms = ignoreInactiveAtoms;
            this.skipEnergies = skipEnergies;
            String molAssFile = molecularAssembly.getFile().getAbsolutePath();
            this.videoDir = FilenameUtils.separatorsToSystem(FilenameUtils.getFullPath(molAssFile) + "video/");
            this.filter = new PDBFilter(new File(videoDir + String.format("%08d", snapshotNum)), molecularAssembly, null, null);
            logger.info(String.format(" SDL: video snapshot directory: %s", videoDir));
            if (ignoreInactiveAtoms) {
                filter.setIgnoreInactiveAtoms(true);
            }
        }

        public void broadcastDummy(int i, int ri) {
            double commEnergy[] = new double[3];
            commEnergy[0] = i;
            commEnergy[1] = ri;
            commEnergy[2] = 0;
            DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
            multicastBuf(commEnergyBuf);
        }

        public void broadcastDummy(int i, int ri, int j, int rj) {
            double commEnergy[] = new double[5];
            commEnergy[0] = i;
            commEnergy[1] = ri;
            commEnergy[2] = j;
            commEnergy[3] = rj;
            commEnergy[4] = 0;
            DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
            multicastBuf(commEnergyBuf);
        }

        public void broadcastDummy(int i, int ri, int j, int rj, int k, int rk) {
            double commEnergy[] = new double[7];
            commEnergy[0] = i;
            commEnergy[1] = ri;
            commEnergy[2] = j;
            commEnergy[3] = rj;
            commEnergy[4] = k;
            commEnergy[5] = rk;
            commEnergy[6] = 0;
            DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
            multicastBuf(commEnergyBuf);
        }

        public void write(String name) {
            File previous = molecularAssembly.getFile();
            if (!name.endsWith(".pdb")) {
                name += ".pdb";
            }
            File file = new File(videoDir + name);
            filter.setFile(file);
            filter.writeFile(file, false);
            molecularAssembly.setFile(previous);
        }

        public void write() {
            File previous = molecularAssembly.getFile();
            String name = String.format("%08d.pdb", snapshotNum);
            File file = new File(videoDir + name);
            filter.setFile(file);
            filter.writeFile(file, false);
            molecularAssembly.setFile(previous);
            snapshotNum++;
        }
    }

    /**
     * Monte Carlo driver for DEE-MC.
     */
    private class RotamerMatrixMC extends BoltzmannMC {

        private final int[] currentRots;
        private final int[] oldRots;
        private final int nRes;
        private final Residue[] residues;

        // Strongly considering adding another internal class to represent the
        // residues-rotamers data, because I'm paranoid somebody won't remember
        // that everybody has to be using the same int[] rotamers.
        /**
         * The rotamers array MUST be the same array as passed to any MCMove
         * objects used, and NOT a copy.
         *
         * @param rotamers
         * @param residues
         */
        RotamerMatrixMC(int[] rotamers, Residue[] residues) {
            currentRots = rotamers; // This is intentional.
            nRes = rotamers.length;
            oldRots = new int[nRes];
            System.arraycopy(rotamers, 0, oldRots, 0, nRes);
            this.residues = residues;
        }

        @Override
        public void revertStep() {
            System.arraycopy(oldRots, 0, currentRots, 0, nRes);
        }

        /**
         * If useFullAMOEBAEnergy is set to true, explicitly evaluates energy,
         * else computes energy from the rotamer energy matrices.
         *
         * @return Energy at the current state
         */
        @Override
        protected double currentEnergy() {
            try {
                return useFullAMOEBAEnergy ? currentEnergyWrapper(Arrays.asList(residues)) : computeEnergy(residues, currentRots, false);
            } catch (NullPointerException ex) {
                // If using the rotamer energy matrix, and there is some missing
                // energy term, just return a default, very large energy.
                return 1e100;
            }
            //return computeEnergy(residues, currentRots, false);
        }

        @Override
        protected void storeState() {
            System.arraycopy(currentRots, 0, oldRots, 0, nRes);
        }
    }

    /**
     * This implements single-rotamer changes in the framework of the rotamer
     * energy matrices.
     */
    private class RotamerMatrixMove implements MCMove {

        private final boolean useAllElims;
        /**
         * currentRots should point to the same array as being used in the
         * overlying MetropolisMC implementation.
         */
        private final int[] currentRots;
        private final int nRes;
        private final List<Integer> allowedRes;
        private final List<List<Integer>> allowedRots;
        private final int nAllowed;
        /**
         * When we take a step, we need to remember which rotamer of which
         * residue was changed.
         */
        private int changedRes;
        private int changedRot;

        /**
         * Constructs the RotamerMatrixMove set; at present, a new object must
         * be made if rotamers or residues are changed outside the scope of this
         * class.
         *
         * @param useAllElims Use eliminated pair/triple info
         * @param rots Initial rotamer set
         */
        RotamerMatrixMove(boolean useAllElims, int[] rots, Residue[] residues) {
            this.useAllElims = useAllElims;
            nRes = rots.length;
            currentRots = rots; // Again, very intentional.

            allowedRes = new ArrayList<>(nRes);
            allowedRots = new ArrayList<>(nRes);

            for (int i = 0; i < nRes; i++) {
                ArrayList<Integer> resAllowed = new ArrayList<>();

                int lenri = residues[i].getRotamers(library).length;
                for (int ri = 0; ri < lenri; ri++) {
                    if (!check(i, ri)) {
                        resAllowed.add(ri);
                    }
                }

                if (resAllowed.size() > 1) {
                    resAllowed.trimToSize();
                    allowedRes.add(i);
                    allowedRots.add(resAllowed);
                }
            }

            ((ArrayList) allowedRes).trimToSize();
            nAllowed = allowedRes.size();
        }

        @Override
        public void move() {
            boolean validMove = !useAllElims;
            int indexI;
            int indexRI;
            do {
                /**
                 * resi and roti correspond to their positions in allowedRes and
                 * allowedRots. indexI and indexRI correspond to their numbers
                 * in the rotamer matrix.
                 */
                int resi = ThreadLocalRandom.current().nextInt(nAllowed);
                indexI = allowedRes.get(resi);
                List<Integer> rotsi = allowedRots.get(resi);
                int lenri = rotsi.size();

                int roti = ThreadLocalRandom.current().nextInt(lenri);
                indexRI = rotsi.get(roti);
                if (useAllElims) {
                    validMove = checkValidMove(indexI, indexRI, currentRots);
                }
            } while (!validMove);

            changedRes = indexI;
            changedRot = currentRots[indexI];

            currentRots[indexI] = indexRI;
        }

        @Override
        public void revertMove() {
            currentRots[changedRes] = changedRot;
        }

        /**
         * This method should not be necessary, because currentRots should be
         * the same array as used in the MetropolisMC.
         *
         * @return The current rotamers array
         */
        private int[] getCurrentRots() {
            return currentRots;
        }

        @Override
        public String toString() {
            return "Rotamer moves utlizing a rotamer energy matrix";
        }
    }

    /**
     * Describes an integer array of rotamers and its energy for DEE-MC.
     */
    private class RotamerPermutation implements Comparable {

        private final int[] rotamers;
        private final double energy;
        private final int nRots;

        private RotamerPermutation(int[] rotamers, double energy) {
            nRots = rotamers.length;
            this.rotamers = new int[nRots];
            System.arraycopy(rotamers, 0, this.rotamers, 0, nRots);
            this.energy = energy;
        }

        private int[] getMCRotamers() {
            int[] ret = new int[nRots];
            System.arraycopy(rotamers, 0, ret, 0, nRots);
            return ret;
        }

        private double getEnergy() {
            return energy;
        }

        /**
         * Equal if the rotamer arrays are identical.
         *
         * @param o
         * @return
         */
        @Override
        public boolean equals(Object o) {
            if (o == null) {
                return false;
            }
            if (!(o instanceof RotamerPermutation)) {
                return false;
            }
            RotamerPermutation other = (RotamerPermutation) o;
            int[] otherRots = other.getMCRotamers();
            if (otherRots.length != nRots) {
                return false;
            }
            for (int i = 0; i < nRots; i++) {
                if (rotamers[i] < otherRots[i]) {
                    return false;
                } else if (rotamers[i] > otherRots[i]) {
                    return false;
                }
            }
            return true;
        }

        @Override
        public int hashCode() {
            int hash = 7;
            hash = 37 * hash + Arrays.hashCode(this.rotamers);
            hash = 37 * hash + this.nRots;
            return hash;
        }

        /**
         * Compares based on the rotamer permutation. Consistent with equals.
         * Much more likely, however, to compare using the rotPermEnComparator.
         *
         * @param o To compare to
         * @return
         */
        @Override
        public int compareTo(Object o) {
            if (o == null) {
                return 0;
            }
            if (!(o instanceof RotamerPermutation)) {
                return 0;
            }
            RotamerPermutation other = (RotamerPermutation) o;
            int[] otherRots = other.getMCRotamers();
            if (otherRots.length != nRots) {
                return 0;
            }
            for (int i = 0; i < nRots; i++) {
                if (rotamers[i] < otherRots[i]) {
                    return -1;
                } else if (rotamers[i] > otherRots[i]) {
                    return 1;
                }
            }
            return 0;
            /*if (energy < other.getEnergy()) {
                return -1;
            } else if (energy > other.getEnergy()) {
                return 1;
            } else {
                return 0;
            }*/
        }
    }

    /**
     * Compares RotamerPermutation objects based on energy instead of the
     * natural ordering. Inconsistent with equals.
     */
    private class rotPermEcomparator implements Comparator<RotamerPermutation> {

        @Override
        public int compare(RotamerPermutation o1, RotamerPermutation o2) {
            if (o1.getEnergy() < o2.getEnergy()) {
                return -1;
            } else if (o1.getEnergy() > o2.getEnergy()) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    public enum Algorithm {

        INDEPENDENT, GLOBAL, GLOBAL_DEE, SLIDING_WINDOW, SLIDING_WINDOW_DEE, BOX_OPTIMIZATION
    }

    public enum Direction {

        FORWARD, BACKWARD
    }
}
