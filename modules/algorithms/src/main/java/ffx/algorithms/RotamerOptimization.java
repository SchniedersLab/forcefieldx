/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.algorithms;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.BiFunction;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import static java.lang.String.format;
import static java.util.Arrays.fill;

import edu.rit.mp.BooleanBuf;
import edu.rit.mp.Buf;
import edu.rit.mp.DoubleBuf;
import edu.rit.pj.*;
import edu.rit.pj.reduction.SharedDouble;

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.MCMove;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.IndexIndexPair;
import ffx.utilities.ObjectPair;

import static ffx.potential.bonded.Residue.ResidueType.NA;
import static ffx.potential.bonded.RotamerLibrary.applyRotamer;

/**
 * Optimize protein side-chain conformations and nucleic acid backbone
 * conformations using rotamers.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @author Stephen D. LuCore
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
    private int evaluatedPermutationsPrint = 0;
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
     * <p>
     * [residue1][rotamer1][residue2][rotamer2]
     */
    private double distanceMatrix[][][][];
    /**
     * Flag to load the distance matrix as needed; if false, matrix is prefilled
     * at the beginning of rotamer optimization.
     */
    private boolean lazyMatrix = false;

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
    private double pairClashThreshold = 20.0;
    /**
     * Pair clash energy threshold (kcal/mol) for MultiResidues.
     */
    private double multiResPairClashAddn = 80.0;
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
     * <p>
     * Very important, to ensure that all possible combinations of delta(i) and
     * delta(i-1) are still represented when it comes time to calculate pair
     * energies. If this pruning factor doesn't cut it, however, it probably
     * wasn't a biologically relevant rotamer anyways.
     * <p>
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

    /**
     * Parallel Java Variables.
     */
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
    private final HashMap<Integer, Integer[]> singlesMap = new HashMap<>();
    private final HashMap<Integer, Integer[]> pairsMap = new HashMap<>();
    private final HashMap<Integer, Integer[]> trimersMap = new HashMap<>();
    private final HashMap<Integer, Integer[]> quadsMap = new HashMap<>();
    private List<String> energiesToWrite;

    private ParallelTeam parallelTeam;
    private GoldsteinPairRegion goldsteinPairRegion;
    private EnergyRegion energyRegion;

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
            logger.info(format(" (KEY) computeQuads: %b", this.computeQuads));
        }
        if (quadCutoffDist != null) {
            double value = Double.parseDouble(quadCutoffDist);
            this.quadCutoffDist = value;
            if (this.quadCutoffDist < 0) {
                quadCutoff = false;
            }
            logger.info(format(" (KEY) quadCutoffDist: %.2f", this.quadCutoffDist));
        }
        if (quadMaxout != null) {
            int value = Integer.parseInt(quadMaxout);
            this.quadMaxout = value;
            logger.info(format(" (KEY) quadMaxout: %d", this.quadMaxout));
        }
        if (skipEnergies != null) {
            boolean value = Boolean.parseBoolean(skipEnergies);
            this.skipEnergies = value;
            logger.info(format(" (KEY) skipEnergies: %b", this.skipEnergies));
        }
        if (undo != null) {
            boolean value = Boolean.parseBoolean(undo);
            this.revert = value;
            logger.info(format(" (KEY) undo: %b", this.revert));
        }
        if (direction != null) {
            Direction value = Direction.valueOf(direction);
            this.direction = value;
            logger.info(format(" (KEY) direction: %s", this.direction.toString()));
        }
        if (increment != null) {
            int value = Integer.parseInt(increment);
            this.increment = value;
            logger.info(format(" (KEY) increment: %d", this.increment));
        }
        if (goldstein != null) {
            boolean value = Boolean.parseBoolean(goldstein);
            this.useGoldstein = value;
            logger.info(format(" (KEY) goldstein: %b", this.useGoldstein));
        }
        if (parallelEnergies != null) {
            boolean value = Boolean.parseBoolean(parallelEnergies);
            this.parallelEnergies = value;
            logger.info(format(" (KEY) parallelEnergies: %b", this.parallelEnergies));
        }
        if (superpositionThreshold != null) {
            Double value = Double.parseDouble(superpositionThreshold);
            this.superpositionThreshold = value;
            logger.info(format(" (KEY) superpositionThreshold: %.2f", this.superpositionThreshold));
        }
        if (ensembleNumber != null) {
            int value = Integer.parseInt(ensembleNumber);
            this.ensembleNumber = value;
            this.ensembleBuffer = 5.0;
            this.ensembleBufferStep = 0.1 * this.ensembleBuffer;
            logger.info(format(" (KEY) ensembleNumber: %d", this.ensembleNumber));
        }
        if (ensembleBuffer != null) {
            double value = Double.parseDouble(ensembleBuffer);
            this.ensembleBuffer = value;
            this.ensembleBufferStep = 0.1 * this.ensembleBuffer;
            logger.info(format(" (KEY) ensembleBuffer: %.2f", this.ensembleBuffer));
        }
        if (ensembleEnergy != null) {
            double value = Double.parseDouble(ensembleEnergy);
            this.ensembleEnergy = value;
            logger.info(format(" (KEY) ensembleEnergy: %.2f", this.ensembleEnergy));
        }
        if (threeBodyCutoffDist != null) {
            double value = Double.parseDouble(threeBodyCutoffDist);
            this.threeBodyCutoffDist = value;
            if (this.threeBodyCutoffDist < 0) {
                threeBodyCutoff = false;
            }
            logger.info(format(" (KEY) threeBodyCutoffDist: %.2f", this.threeBodyCutoffDist));
        }
        if (pruningFactor != null) {
            double value = Double.parseDouble(pruningFactor);
            this.pruningFactor = (value >= 0 ? value : 1.0);
            this.pairHalfPruningFactor = (1.0 + value) / 2;
            logger.info(format(" (KEY) pruningFactor: %.2f", this.pruningFactor));
        }
        if (nucleicSinglesPruningFactor != null) {
            double value = Double.parseDouble(nucleicSinglesPruningFactor);
            this.singletonNAPruningFactor = (value >= 0 ? value : 1.5);
            logger.info(format(" (KEY) nucleicSinglesPruningFactor: %.2f", this.singletonNAPruningFactor));
        }
        if (nucleicCorrectionThreshold != null) {
            double value = Double.parseDouble(nucleicCorrectionThreshold);
            this.nucleicCorrectionThreshold = (value >= 0 ? value : 0);
            logger.info(format(" (KEY) nucleicCorrectionThreshold: %.2f", this.nucleicCorrectionThreshold));
        }
        if (minimumNumberAcceptedNARotamers != null) {
            int value = Integer.parseInt(minimumNumberAcceptedNARotamers);
            this.minNumberAcceptedNARotamers = (value > 0 ? value : 10);
            logger.info(format(" (KEY) minimumNumberAcceptedNARotamers: %d", this.minNumberAcceptedNARotamers));
        }
        if (singletonClashThreshold != null) {
            double value = Double.parseDouble(singletonClashThreshold);
            this.clashThreshold = value;
            logger.info(format(" (KEY) singletonClashThreshold: %.2f", this.clashThreshold));
        }
        if (multiResClashThreshold != null) {
            double value = Double.parseDouble(multiResClashThreshold);
            this.multiResClashThreshold = value;
            logger.info(format(" (KEY) multiResClashThreshold: %.2f", this.multiResClashThreshold));
        }
        if (pairClashThreshold != null) {
            double value = Double.parseDouble(pairClashThreshold);
            this.pairClashThreshold = value;
            logger.info(format(" (KEY) pairClashThreshold: %.2f", this.pairClashThreshold));
        }
        if (multiResPairClashAddition != null) {
            double value = Double.parseDouble(multiResPairClashAddition);
            this.multiResPairClashAddn = value;
            logger.info(format(" (KEY) multiResPairClashAddition: %.2f", this.multiResPairClashAddn));
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
                            logger.info(format(" Improper dimension min %8.5f > max %8.5f; max/min reversed.", this.boxDimensions[i - 1], this.boxDimensions[i]));
                            double temp = this.boxDimensions[i];
                            this.boxDimensions[i] = this.boxDimensions[i - 1];
                            this.boxDimensions[i - 1] = temp;
                        }
                    }
                    superboxBuffer = Double.parseDouble(bdTokens[0]);
                    manualSuperbox = true;
                }
            } catch (Exception ex) {
                logger.warning(format(" Error in parsing box dimensions: input discarded and defaults used: %s.", ex.toString()));
                manualSuperbox = false;
            }
        }
        if (mcTemp != null) {
            double value = Double.parseDouble(mcTemp);
            this.mcTemp = value;
            logIfMaster(format(" (KEY) mcTemp: %10.6f", this.mcTemp));
        }
        if (mcUseAll != null) {
            boolean value = Boolean.parseBoolean(mcUseAll);
            this.mcUseAll = value;
            logIfMaster(format(" (KEY) mcUseAll: %b", this.mcUseAll));
        }
        if (mcNoEnum != null) {
            boolean value = Boolean.parseBoolean(mcNoEnum);
            this.mcNoEnum = value;
            logIfMaster(format(" (KEY) debug-mcNoEnum: %b", this.mcNoEnum));
        }
        if (lazyMatrix != null) {
            boolean value = Boolean.parseBoolean(lazyMatrix);
            this.lazyMatrix = value;
            logger.info(format(" (KEY) lazyMatrix: %b", lazyMatrix));
        }
        if (addOrigRotStr != null) {
            boolean value = Boolean.parseBoolean(addOrigRotStr);
            this.addOrigRot = value;
            logger.info(format(" (KEY) addOrigRot: %b", addOrigRot));
        }
        if (origAtEndStr != null) {
            boolean value = Boolean.parseBoolean(origAtEndStr);
            // Property works in the contest of Residue class.
            logger.info(format(" (KEY) origAtEnd: %b", value));
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

                double rotEnergy = 1E100;
                try {
                    rotEnergy = currentEnergy(resList);
                    logger.info(format(" %d Energy: %16.8f", ++evaluatedPermutations, rotEnergy));
                } catch (ArithmeticException ex) {
                    logger.info(String.format(" %d Energy set to %10.6g (unreasonable conformation)", ++evaluatedPermutations, rotEnergy));
                }
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
     * @param residues          Optimization window
     * @param i                 Current residue in the recursion.
     * @param lowEnergy         Minimum energy yet found by the recursion.
     * @param optimum           Optimum rotamer set yet found by the recursion.
     * @param currentRotamers   Rotamer permutation under investigation.
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
                        logger.info(format(" Minimum energy update: %f < %f, permutation %d",
                                rotEnergy, lowEnergy, evaluatedPermutations));
                        String permutation = " Rotamer permutation: " + optimum[0];
                        for (int j = 1; j < nResidues; j++) {
                            permutation = permutation.concat(", " + optimum[j]);
                        }
                        logger.info(permutation);
                    } else {
                        logger.info(format(" First minimum energy (permutation 1): %f", rotEnergy));
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
     * @param residues      Optimization window
     * @param optimum       Array to store optimum rotamers
     * @param initialRots   Array with starting rotamers
     * @param maxIters      Number of MC steps to run
     * @param randomizeRots Scramble initialRots
     * @param useAllElims   Use pair/triple elimination information
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

        /**
         * I have the vague idea of parallelizing down to individual threads, by
         * simply calling this method from each thread, with randomizeRots true.
         * That would require replacing logIfMaster with logger.info.
         */
        logIfMaster(format(" Beginning %d iterations of Monte Carlo search "
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
        }

        initTime += System.nanoTime();
        double fractAccept = ((double) nAccept) / ((double) maxIters);
        logIfMaster(format(" %d steps of DEE-MC completed in %10.6f seconds",
                maxIters, (initTime * 1.0E-9)));
        logIfMaster(format(" Number of steps accepted: %d for %10.6f of total", nAccept, fractAccept));
        logIfMaster(format(" Lowest energy found: %10.6f kcal/mol", optimumEnergy));
        logIfMaster(format(" Final energy found: %10.6f kcal/mol", currentEnergy));

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
     * Checks the pair elimination array to see if this permutation
     * has been eliminated.
     *
     * @param i  Residue number
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
     * @param optimum             Optimum set of rotamers.
     * @param permutationEnergies Energies of visited permutations or null.
     * @return current energy.
     */
    private double rotamerOptimizationDEE(MolecularAssembly molecularAssembly, Residue residues[], int i,
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
                    double amoebaEnergy = 1E100;
                    try {
                        amoebaEnergy = currentEnergy(resList);
                    } catch (ArithmeticException ex) {
                        logger.warning(String.format(" Exception %s in calculating full AMOEBA energy for permutation %d", ex.toString(), evaluatedPermutations));
                    }
                    comparisonEnergy = amoebaEnergy;
                }
                if (permutationEnergies != null) {
                    permutationEnergies[evaluatedPermutations - 1] = comparisonEnergy;
                }
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }
                if (ensembleNumber > 1) {
                    if (master && printFiles) {
                        try {
                            FileWriter fw = new FileWriter(ensembleFile, true);
                            BufferedWriter bw = new BufferedWriter(fw);
                            bw.write(format("MODEL        %d", evaluatedPermutations));
                            for (int j = 0; j < 75; j++) {
                                bw.write(" ");
                            }
                            bw.newLine();
                            bw.flush();
                            ensembleFilter.writeFile(ensembleFile, true);
                            bw.write(format("ENDMDL"));
                            for (int j = 0; j < 64; j++) {
                                bw.write(" ");
                            }
                            bw.newLine();
                            bw.close();
                        } catch (IOException e) {
                            logger.warning(format("Exception writing to file: %s", ensembleFile.getName()));
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

                if (useFullAMOEBAEnergy) {
                    // Log current results
                    logIfMaster(format(" %6e AMOEBA: %12.4f 3-Body: %12.4f Neglected: %12.4f (%12.4f)",
                            (double) evaluatedPermutations, comparisonEnergy, approximateEnergy,
                            comparisonEnergy - approximateEnergy, lowEnergy));
                } else {
                    if (threeBodyTerm) {
                        logIfMaster(format(" %6e 3-Body: %12.4f (%12.4f)",
                                (double) evaluatedPermutations, approximateEnergy, lowEnergy));
                    } else {
                        logIfMaster(format(" %6e 2-Body: %12.4f (%12.4f)",
                                (double) evaluatedPermutations, approximateEnergy, lowEnergy));
                    }
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
     * @return 0.
     */
    private double dryRun(Residue residues[], int i, int currentRotamers[]) {
        // This is the initialization condition.
        if (i == 0) {
            evaluatedPermutations = 0;
            evaluatedPermutationsPrint = 1000;
        }

        if (evaluatedPermutations >= evaluatedPermutationsPrint) {
            if (evaluatedPermutations % evaluatedPermutationsPrint == 0) {
                logIfMaster(format(" The permutations have reached %10.4e.", (double) evaluatedPermutationsPrint));
                evaluatedPermutationsPrint *= 10;
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
             * At the end of the recursion, check each rotamer of the final residue.
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
                if (!deadEnd) {
                    evaluatedPermutations++;
                }
            }
        }

        return 0.0;
    }

    /**
     * Finds all permutations within buffer energy of GMEC.
     *
     * @param residues
     * @param i                   Current depth in residue/rotamer tree.
     * @param currentRotamers     Current set of rotamers at this node.
     * @param gmecEnergy          Minimum energy for these residues.
     * @param permutationEnergies Energy of all permutations.
     * @param permutations        Contains accepted permutations.
     * @return 0.
     */
    private double dryRunForEnsemble(Residue residues[], int i, int currentRotamers[],
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
     * @param print    Verbosity flag
     * @return Approximate permutation energy (backbone + selfs + pairs +
     * trimers).
     */
    private double computeEnergy(Residue residues[], int rotamers[], boolean print) {
        double selfSum = 0.0;
        double pairSum = 0.0;
        double threeBodySum = 0.0;
        try {
            if (parallelTeam == null) {
                parallelTeam = new ParallelTeam();
            }
            if (energyRegion == null) {
                energyRegion = new EnergyRegion(parallelTeam.getThreadCount());
            }
            energyRegion.init(residues, rotamers);
            parallelTeam.execute(energyRegion);
            selfSum = energyRegion.getSelf();
            pairSum = energyRegion.getPair();
            threeBodySum = energyRegion.getThreeBody();
        } catch (Exception e) {
            logger.log(Level.WARNING, " Exception in EnergyRegion.", e);
        }

        double approximateEnergy = backboneEnergy + selfSum + pairSum + threeBodySum;
        if (print) {
            logger.info(format(" Backbone:                  %16.8f", backboneEnergy));
            logger.info(format(" Self Energy:               %16.8f", selfSum));
            logger.info(format(" Pair Energy:               %16.8f", pairSum));
            if (!threeBodyTerm) {
                logger.info(format(" Total Energy up to 2-Body: %16.8f", approximateEnergy));
            } else {
                logger.info(format(" 3-Body Energy:             %16.8f", threeBodySum));
                logger.info(format(" Total Energy up to 3-Body: %16.8f", approximateEnergy));
            }
        }
        return approximateEnergy;
    }

    private class EnergyRegion extends ParallelRegion {

        private SharedDouble self;
        private SharedDouble pair;
        private SharedDouble threeBody;
        private EnergyLoop energyLoops[];
        private Residue residues[];
        private int rotamers[];
        private int nResidues;

        public EnergyRegion(int nThreads) {
            self = new SharedDouble();
            pair = new SharedDouble();
            threeBody = new SharedDouble();
            energyLoops = new EnergyLoop[nThreads];
        }

        public void init(Residue residues[], int rotamers[]) {
            this.residues = residues;
            this.rotamers = rotamers;
            this.nResidues = residues.length;
        }

        public void start() {
            self.set(0.0);
            pair.set(0.0);
            threeBody.set(0.0);
        }

        public double getSelf() {
            return self.get();
        }

        public double getPair() {
            return pair.get();
        }

        public double getThreeBody() {
            return threeBody.get();
        }

        @Override
        public void run() {
            int threadID = getThreadIndex();
            if (energyLoops[threadID] == null) {
                energyLoops[threadID] = new EnergyLoop();
            }
            try {
                execute(0, nResidues - 1, energyLoops[threadID]);
            } catch (Exception e) {
                logger.log(Level.WARNING, " Exception in EnergyLoop.", e);
            }
        }

        private class EnergyLoop extends IntegerForLoop {
            private double selfSum;
            private double pairSum;
            private double threeBodySum;

            @Override
            public void start() {
                selfSum = 0.0;
                pairSum = 0.0;
                threeBodySum = 0.0;
            }

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.dynamic();
            }

            @Override
            public void finish() {
                self.addAndGet(selfSum);
                pair.addAndGet(pairSum);
                threeBody.addAndGet(threeBodySum);
            }

            @Override
            public void run(int lb, int ub) {
                for (int a = lb; a <= ub; a++) {
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
            }
        }
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
        try {
            totalEnergy = currentEnergy(allResiduesList);
        } catch (ArithmeticException ex) {
            totalEnergy = 0;
            logger.severe(String.format(" Exception %s in calculating total energy at the start of decompose-original", ex.toString()));
        }
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        // List to contain current residues.
        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));

        try {
            localBackboneEnergy = currentEnergy(rList);
        } catch (ArithmeticException ex) {
            logger.severe(format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        for (int i = 0; i < nRes; i++) {
            long time = -System.nanoTime();
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            // Unchecked use because the original structure should be reasonable.
            localSelfEnergy[i] = currentEnergy(rList) - localBackboneEnergy;
            residueEnergy[0][i] = localSelfEnergy[i];
            turnOffAtoms(ri);
            time += System.nanoTime();
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                long time = -System.nanoTime();
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                // Unchecked use because the original structure should be reasonable.
                pairEnergy[i][j] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localBackboneEnergy;
                double halfPair = pairEnergy[i][j] * 0.5;
                residueEnergy[1][i] += halfPair;
                residueEnergy[1][j] += halfPair;
                turnOffAtoms(rj);
                time += System.nanoTime();
                logger.info(format(" Pair %s %s:       %16.5f in %6.4f (sec).", ri, rj, pairEnergy[i][j], time * 1.0e-9));
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
                        long time = -System.nanoTime();
                        turnOnAtoms(rk);
                        // Unchecked use because the original structure should be reasonable.
                        triEnergy[i][j][k] = currentEnergy(rList)
                                - localSelfEnergy[i] - localSelfEnergy[j] - localSelfEnergy[k]
                                - pairEnergy[i][j] - pairEnergy[j][k] - pairEnergy[i][k] - localBackboneEnergy;
                        double thirdTrimer = triEnergy[i][j][k] / 3.0;
                        residueEnergy[2][i] += thirdTrimer;
                        residueEnergy[2][j] += thirdTrimer;
                        residueEnergy[2][k] += thirdTrimer;
                        turnOffAtoms(rk);
                        time += System.nanoTime();
                        logger.info(format(" Tri  %s %s %s:    %16.5f in %6.4f (sec).", ri, rj, rk, triEnergy[i][j][k], time * 1.0e-9));
                    } else if (dist == Double.MAX_VALUE) {
                        logger.info(format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)", ri, rj, rk));
                        triEnergy[i][j][k] = 0.0;
                    } else {
                        logger.info(format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms", ri, rj, rk, dist));
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
        try {
            totalEnergy = currentEnergy(allResiduesList);
        } catch (ArithmeticException ex) {
            totalEnergy = 0;
            logger.severe(String.format(" Exception %s in calculating total energy for decomposition; FFX shutting down", ex.toString()));
        }
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));

        try {
            localBackboneEnergy = currentEnergy(rList);
        } catch (ArithmeticException ex) {
            logger.severe(format(" FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        for (int i = 0; i < nRes; i++) {
            long time = -System.nanoTime();
            Residue ri = residues[i];
            rList.set(1, ri);
            turnOnAtoms(ri);
            // Unchecked use because the original structure should be reasonable.
            localSelfEnergy[i] = currentEnergy(rList) - localBackboneEnergy;
            residueEnergy[0][i] = localSelfEnergy[i];
            turnOffAtoms(ri);
            time += System.nanoTime();
            logger.info(format(" Self %s:          %16.5f in %6.4f (sec).", ri, localSelfEnergy[i], time * 1.0e-9));
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            for (int j = i + 1; j < nRes; j++) {
                long time = -System.nanoTime();
                Residue rj = residues[j];
                rList.set(1, rj);
                turnOnAtoms(rj);
                // Unchecked use because the original structure should be reasonable.
                pairEnergy[i][j] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localBackboneEnergy;
                double halfPair = pairEnergy[i][j] * 0.5;
                residueEnergy[1][i] += halfPair;
                residueEnergy[1][j] += halfPair;
                turnOffAtoms(rj);
                time += System.nanoTime();
                logger.info(format(" Pair %s %s:       %16.5f in %6.4f (sec).", ri, rj, pairEnergy[i][j], time * 1.0e-9));
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
                        long time = -System.nanoTime();
                        turnOnAtoms(rk);
                        // Unchecked use because the original structure should be reasonable.
                        triEnergy[i][j][k] = currentEnergy(rList)
                                - localSelfEnergy[i] - localSelfEnergy[j] - localSelfEnergy[k]
                                - pairEnergy[i][j] - pairEnergy[j][k] - pairEnergy[i][k] - localBackboneEnergy;
                        double thirdTrimer = triEnergy[i][j][k] / 3.0;
                        residueEnergy[2][i] += thirdTrimer;
                        residueEnergy[2][j] += thirdTrimer;
                        residueEnergy[2][k] += thirdTrimer;
                        turnOffAtoms(rk);
                        time += System.nanoTime();
                        logger.info(format(" Tri  %s %s %s:    %16.5f in %6.4f (sec).", ri, rj, rk, triEnergy[i][j][k], time * 1.0e-9));
                    } else if (dist == Double.MAX_VALUE) {
                        logger.info(format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)", ri, rj, rk));
                        triEnergy[i][j][k] = 0.0;
                    } else {
                        logger.info(format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms", ri, rj, rk, dist));
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
        try {
            totalEnergy = currentEnergy(allResiduesList);
        } catch (ArithmeticException ex) {
            totalEnergy = 0;
            logger.severe(String.format(" Exception %s in calculating total energy for decomposition; FFX shutting down", ex.toString()));
        }
        logIfMaster(format(" AMOEBA:   %16.5f", totalEnergy));
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));
        try {
            localBackboneEnergy = currentEnergy(rList);
        } catch (ArithmeticException ex) {
            logger.severe(format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        logIfMaster(format(" Backbone: %16.5f", localBackboneEnergy));

        decomposeOriginal = true;
        allocateEliminationMemory(allResiduesArray);
        rotamerEnergies(allResiduesArray);

        if (master) {
            for (int i = 0; i < nRes; i++) {
                Residue ri = residues[i];
                localSelfEnergy[i] = selfEnergy[i][0];
                residueEnergy[0][i] = localSelfEnergy[i];
                logger.info(format(" Self %s:          %16.5f", ri, localSelfEnergy[i]));
            }
            for (int i = 0; i < nRes; i++) {
                Residue ri = residues[i];
                for (int j = i + 1; j < nRes; j++) {
                    Residue rj = residues[j];
                    pairEnergy[i][j] = twoBodyEnergy[i][0][j][0];
                    logger.info(format(" Pair %s %s:       %16.5f", ri, rj, pairEnergy[i][j]));
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
                                logger.info(format(" Tri  %s %s %s:    %16.5f", ri, rj, rk, triEnergy[i][j][k]));
                            } else if (dist == Double.MAX_VALUE) {
                                logger.info(format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)",
                                        ri, rj, rk));
                                triEnergy[i][j][k] = 0.0;
                            } else if (dist > threeBodyCutoffDist) {
                                logger.info(format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms",
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
     * Test the self-energy elimination by setting
     * 2-body and 3-body interactions to zero.
     */
    public void testSelfEnergyElimination(Residue residues[]) {
        int nRes = residues.length;
        for (int i = 0; i < nRes; i++) {
            Residue resI = residues[i];
            Rotamer[] rotI = resI.getRotamers(library);
            int nI = rotI.length;
            for (int ri = 0; ri < nI; ri++) {
                for (int j = i + 1; j < nRes; j++) {
                    Residue resJ = residues[j];
                    Rotamer[] rotJ = resJ.getRotamers(library);
                    int nJ = rotJ.length;
                    for (int rj = 0; rj < nJ; rj++) {
                        try {
                            twoBodyEnergy[i][ri][j][rj] = 0.0;
                        } catch (Exception e) {
                            // catch NPE.
                        }
                        if (threeBodyTerm) {
                            for (int k = j + 1; k < nRes; k++) {
                                Residue resK = residues[k];
                                Rotamer[] rotK = resK.getRotamers(library);
                                int nK = rotK.length;
                                for (int rk = 0; rk < nK; rk++) {
                                    try {
                                        threeBodyEnergy[i][ri][j][rj][k][rk] = 0.0;
                                    } catch (Exception e) {
                                        // catch NPE.
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Test the self-energy elimination by setting
     * 2-body and 3-body interactions to zero.
     */
    public void testPairEnergyElimination(Residue residues[], int resID) {
        int nRes = residues.length;

        if (resID >= nRes) {
            return;
        }

        for (int i = 0; i < nRes; i++) {
            Residue resI = residues[i];
            Rotamer[] rotI = resI.getRotamers(library);
            int nI = rotI.length;
            for (int ri = 0; ri < nI; ri++) {
                try {
                    selfEnergy[i][ri] = 0.0;
                } catch (Exception e) {
                    // catch NPE.
                }
                for (int j = i + 1; j < nRes; j++) {
                    Residue resJ = residues[j];
                    Rotamer[] rotJ = resJ.getRotamers(library);
                    int nJ = rotJ.length;
                    for (int rj = 0; rj < nJ; rj++) {
                        if (i != resID && j != resID) {
                            try {
                                twoBodyEnergy[i][ri][j][rj] = 0.0;
                            } catch (Exception e) {
                                // catch NPE.
                            }
                        }
                        if (threeBodyTerm) {
                            for (int k = j + 1; k < nRes; k++) {
                                Residue resK = residues[k];
                                Rotamer[] rotK = resK.getRotamers(library);
                                int nK = rotK.length;
                                for (int rk = 0; rk < nK; rk++) {
                                    try {
                                        threeBodyEnergy[i][ri][j][rj][k][rk] = 0.0;
                                    } catch (Exception e) {
                                        // catch NPE.
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
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
        try {
            totalEnergy = currentEnergy(allResiduesList);
        } catch (ArithmeticException ex) {
            totalEnergy = 0;
            logger.severe(String.format(" Exception %s in calculating total energy for decomposition; FFX shutting down", ex.toString()));
        }
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }

        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));
        try {
            localBackboneEnergy = currentEnergy(rList);
        } catch (ArithmeticException ex) {
            logger.severe(format("FFX shutting down: error in calculation of backbone energy %s", ex.getMessage()));
        }
        for (int i = 0; i < nRes; i++) {
            Residue ri = residues[i];
            rList.set(0, ri);
            turnOnAtoms(ri);
            // Unchecked use because the original structure should be reasonable.
            localSelfEnergy[i] = currentEnergy(rList) - localBackboneEnergy;
            logger.info(format(" Self %s:          %16.5f", ri, localSelfEnergy[i]));
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
                // Unchecked use because the original structure should be reasonable.
                pairEnergy[i][j] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localBackboneEnergy;
                logger.info(format(" Pair %s %s:       %16.5f", ri, rj, pairEnergy[i][j]));
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
                        // Unchecked use because the original structure should be reasonable.
                        triEnergy[i][j][k] = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j] - localSelfEnergy[k]
                                - pairEnergy[i][j] - pairEnergy[j][k] - pairEnergy[i][k] - localBackboneEnergy;
                        logger.info(format(" Tri  %s %s %s:    %16.5f", ri, rj, rk, triEnergy[i][j][k]));
                        sumTri += triEnergy[i][j][k];
                        turnOffAtoms(rk);
                    } else if (dist == Double.MAX_VALUE) {
                        logger.info(format(" Tri  %s %s %s:    set to 0.0 at NaN (very long distance)", ri, rj, rk));
                        triEnergy[i][j][k] = 0.0;
                    } else {
                        logger.info(format(" Tri  %s %s %s:    set to 0.0 at %1.5f Angstroms", ri, rj, rk, dist));
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
                            // Unchecked use because the original structure should be reasonable.
                            double currentQuad = currentEnergy(rList) - localSelfEnergy[i] - localSelfEnergy[j]
                                    - localSelfEnergy[k] - localSelfEnergy[l] - pairEnergy[i][j] - pairEnergy[i][k]
                                    - pairEnergy[i][l] - pairEnergy[j][k] - pairEnergy[j][l] - pairEnergy[k][l]
                                    - triEnergy[i][j][k] - triEnergy[i][j][l] - triEnergy[i][k][l] - triEnergy[j][k][l]
                                    - localBackboneEnergy;
                            logger.info(format(" Quad  %s %s %s %s:    %16.5f at %1.5f Angstroms",
                                    ri, rj, rk, rl, currentQuad, dist));
                            sumQuads += currentQuad;
                            turnOffAtoms(rl);
                            if (++numQuadsEvaluated >= maxQuads) {
                                doBreak = true;
                                break;
                            }
                        } else if (dist == Double.MAX_VALUE) {
                            logger.info(format(" Quad  %s %s %s %s:    set to 0.0 at NaN (very long distance)",
                                    ri, rj, rk, rl));
                        } else {
                            logger.info(format(" Quad  %s %s %s %s:    set to 0.0 at %1.5f Angstroms",
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

        logger.info(format("\n\n"));
        logger.info(format(" Backbone:     %16.5f", localBackboneEnergy));
        logger.info(format(" Sum Self:     %16.5f", sumSelf));
        logger.info(format(" Sum Pair:     %16.5f", sumPair));
        logger.info(format(" Sum Tri:      %16.5f", sumTri));
        logger.info(format(" Sum Quad:     %16.5f", sumQuads));
        logger.info(format(" Neglected:    %16.5f", totalEnergy - sumSelf - sumPair - sumTri - sumQuads - localBackboneEnergy));
        logger.info(format(" AMOEBA:       %16.5f", totalEnergy));
        logger.info(format("\n\n"));
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
            logger.info(format("\n"));
            logger.info(format(" Molecule-Based Many-Body Energy Summation\n "));
            double moleculeEnergy[][] = new double[3][nMol];
            for (int i = 0; i < nRes; i++) {
                int iMolecule = moleculeMap.get(residues[i]);
                moleculeEnergy[0][iMolecule] += residueEnergy[0][i];
                moleculeEnergy[1][iMolecule] += residueEnergy[1][i];
                moleculeEnergy[2][iMolecule] += residueEnergy[2][i];
            }
            logger.info(format(" %9s %9s %9s %9s %9s", "Molecule", "Self", "Pair", "3-Body", "Total"));
            for (int i = 0; i < nMol; i++) {
                String molName = molNames.get(i);
                double total = moleculeEnergy[0][i] + moleculeEnergy[1][i] + moleculeEnergy[2][i];
                logger.info(format(" %9s %9.3f %9.3f %9.3f %9.3f",
                        molName, moleculeEnergy[0][i], moleculeEnergy[1][i], moleculeEnergy[2][i], total));
            }
            logger.info(format(" %9s %9.3f %9.3f %9.3f %9.3f",
                    "Sum", sumSelf, sumPair, sumTri, sumSelf + sumPair + sumTri));
        }

        logger.info(format("\n"));
        logger.info(format(" Residue-Based Many-Body Energy Summation\n "));
        logger.info(format(" %9s %9s %9s %9s %9s", "Residue", "Self", "Pair", "3-Body", "Total"));

        for (int i = 0; i < nRes; i++) {
            double total = residueEnergy[0][i] + residueEnergy[1][i] + residueEnergy[2][i];
            logger.info(format(" %9s %9.3f %9.3f %9.3f %9.3f",
                    residues[i].toString(), residueEnergy[0][i], residueEnergy[1][i], residueEnergy[2][i], total));
        }
        logger.info(format(" %9s %9.3f %9.3f %9.3f %9.3f",
                "Sum", sumSelf, sumPair, sumTri, sumSelf + sumPair + sumTri));
        logger.info(format(" Backbone:        %9.3f", backbone));
        double target = sumSelf + sumPair + sumTri + backbone;
        logger.info(format(" Expansion Total: %9.3f", target));
        logger.info(format(" True Total:      %9.3f", totalEnergy));
        logger.info(format(" Neglected:       %9.3f\n", totalEnergy - target));
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
            logger.info(format(" Large pair interactions (>%.2f):", pairCutoff));
            for (int i = 0; i < nRes; i++) {
                for (int j = i + 1; j < nRes; j++) {
                    if (Math.abs(twoBodyEnergy[i][0][j][0]) >= pairCutoff) {
                        logger.info(format(" Large Pair %s %s:       %16.5f",
                                residues[i], residues[j], twoBodyEnergy[i][0][j][0]));
                    }
                }
            }
            logger.info(format("\n Large trimer interactions (>%.2f):", trimerCutoff));
            for (int i = 0; i < nRes; i++) {
                for (int j = i + 1; j < nRes; j++) {
                    for (int k = j + 1; k < nRes; k++) {
                        if (Math.abs(threeBodyEnergy[i][0][j][0][k][0]) >= trimerCutoff) {
                            logger.info(format(" Large Trimer  %s %s %s:    %16.5f",
                                    residues[i], residues[j], residues[k], threeBodyEnergy[i][0][j][0][k][0]));
                        }
                    }
                }
            }
            return;
        }

        logger.info(format(" Large pair interactions (>%.2f):", pairCutoff));
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
                                logger.info(format(" Large Pair %7s %-2d, %7s %-2d: %16.8f",
                                        resi, ri, resj, rj, twoBodyEnergy[i][ri][j][rj]));
                            }
                        } catch (Exception ex) {
                        }
                    }
                }
            }
        }

        logger.info(format("\n Large trimer interactions (>%.2f):", trimerCutoff));
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
                                        logger.info(format(" Large Trimer %7s %-2d, %7s %-2d, %7s %-2d: %16.8f",
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
                logger.info(format("\t%s", r.toString()));
            } else {
                logger.info(format(" not \t%s", r.toString()));
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

        logger.info(format("\n Rotamer Library:     %s", library.getLibrary()));
        logger.info(format(" Algorithm:           %s", algorithm));
        logger.info(format(" Goldstein Criteria:  %b", useGoldstein));
        logger.info(format(" Three-Body Energies: %b\n", threeBodyTerm));

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
                        case WINDOW:
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
                case BRUTE_FORCE:
                    e = bruteForce(residueList);
                    break;
                case ALL:
                    e = globalOptimization(residueList);
                    break;
                case WINDOW:
                    e = slidingWindowOptimization(residueList, windowSize, increment, revert, distance, direction);
                    break;
                case BOX:
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
     * <p>
     * A nucleic correction threshold of 0 skips the entire method; this check
     * is presently being performed inside the method in case it is called again
     * at some point.
     *
     * @param residues              Residues to eliminate bad backbone rotamers over.
     * @param numEliminatedRotamers Number of previously eliminated rotamers per
     *                              residue.
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
            logIfMaster(format(" Eliminating nucleic acid rotamers with correction vectors larger than %5.3f A", nucleicCorrectionThreshold));
            logIfMaster(format(" A minimum of %d rotamers per NA residue will carry through to energy calculations.", minNumberAcceptedNARotamers));
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
                        logIfMaster(format(" Correction magnitude was %6.4f A > %5.3f A", rotToReject.getDoubleValue(), nucleicCorrectionThreshold));
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
                            logIfMaster(format(" Rotamer %d of residue %s eliminated "
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
                            logIfMaster(format(" Rotamer %d of residue %s eliminated "
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
                                    logIfMaster(format(" Rotamer %d of residue %s eliminated "
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
                                    logIfMaster(format(" Rotamer %d of residue %s eliminated "
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

    public void setUseGoldstein(boolean useGoldstein) {
        this.useGoldstein = useGoldstein;
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
     * @param nMCsteps   Number of steps to be taken
     */
    public void setMonteCarlo(boolean monteCarlo, int nMCsteps) {
        this.monteCarlo = monteCarlo;
        this.nMCsteps = nMCsteps;
    }

    public void setNucleicCorrectionThreshold(double nucleicCorrectionThreshold) {
        if (nucleicCorrectionThreshold >= 0) {
            this.nucleicCorrectionThreshold = nucleicCorrectionThreshold;
        } else {
            logger.warning("\n Correction threshold must be >= 0. Setting to default of 0 (threshold inactive).\n");
            this.nucleicCorrectionThreshold = 0;
        }
    }

    public void setMinimumNumberAcceptedNARotamers(int minNumberAcceptedNARotamers) {
        if (minNumberAcceptedNARotamers > 0) {
            this.minNumberAcceptedNARotamers = minNumberAcceptedNARotamers;
        } else {
            logger.warning("\n Minimum number of accepted NA rotamers must be a positive integer.\n Setting to default value 10.\n");
            this.minNumberAcceptedNARotamers = 10;
        }
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
     *
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
     * @param endForced   Last residue to force.
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
            logger.info(format(" Optimizing %s side-chain.", residue));
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
                double newE = 1E100;
                try {
                    newE = currentEnergy(rList);
                } catch (ArithmeticException ex) {
                    logger.fine(String.format(" Exception %s in energy calculations during independent for %s-%d", ex.toString(), residue, j));
                }
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
                logger.info(format("\n Terminating after residue %s.\n", residue));
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
     * @return If valid permutation found.
     */
    private boolean firstValidPerm(Residue residues[], int i, int currentRotamers[]) {
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
                if (firstValidPerm(residues, i + 1, currentRotamers)) {
                    return true;
                }
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
    private double globalOptimization(List<Residue> residueList) {
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
             * Compute the number of permutations without eliminating dead-ends and compute the number
             * of permutations using singleton elimination.
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
                logger.severe(" No valid path through rotamer space found; try recomputing without pruning or using ensemble.");
            }
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
                        logger.warning(format(" Versioning failed: appending to end of file %s", ensembleFile.getName()));
                    }
                }
                ensembleFilter = new PDBFilter(new File(ensembleFile.getName()), molecularAssembly, null, null);
                logger.info(format(" Ensemble file: %s", ensembleFile.getName()));
            }
            logIfMaster(format(" Number of permutations without DEE conditions: %10.4e.", permutations));
            logIfMaster(format(" Number of permutations after singleton eliminations: %10.4e.", singletonPermutations));
            logIfMaster(format(" Number of permutations removed by pairwise eliminations: %10.4e.", pairTotalElimination));
            logIfMaster(format(" Number of permutations remaining: %10.4e.", (double) evaluatedPermutations));

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
                logIfMaster(format("\n Checking permutations for distance < %5.3f kcal/mol from GMEC energy %10.8f kcal/mol", ensembleEnergy, e));
                dryRunForEnsemble(residues, 0, currentRotamers, e, permutationEnergies, acceptedPermutations);
                int numAcceptedPermutations = 0;

                for (int i = 0; i < acceptedPermutations.length; i++) {
                    if (acceptedPermutations[i] != null) {
                        ++numAcceptedPermutations;
                        logIfMaster(format(" Accepting permutation %d at %8.6f < %8.6f", i, permutationEnergies[i] - e, ensembleEnergy));
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
                                bw.write(format("MODEL        %d", numAcceptedPermutations));
                                for (int j = 0; j < 75; j++) {
                                    bw.write(" ");
                                }
                                bw.newLine();
                                bw.flush();
                                ensembleFilter.writeFile(ensembleFile, true);
                                bw.write(format("ENDMDL"));
                                for (int j = 0; j < 64; j++) {
                                    bw.write(" ");
                                }
                                bw.newLine();
                                bw.close();
                            } catch (IOException ex) {
                                logger.warning(format(" Exception writing to file: %s", ensembleFile.getName()));
                            }
                        }
                    }
                }
                logIfMaster(format(" Number of permutations within %5.3f kcal/mol of GMEC energy: %6.4e",
                        ensembleEnergy, (double) numAcceptedPermutations));
                ensembleStates.sort(null);
            }

            logIfMaster("\n Final rotamers:");
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                int ri = optimum[i];
                Rotamer rotamer = rotamers[ri];
                logIfMaster(format(" %c (%7s,2d) %s", residue.getChainID(), residue, ri, rotamer.toAngleString()));
                RotamerLibrary.applyRotamer(residue, rotamer);
            }

            double sumSelfEnergy = 0;
            double sumPairEnergy = 0;
            double sumTrimerEnergy = 0;
            for (int i = 0; i < nResidues; i++) {
                int ri = optimum[i];
                sumSelfEnergy += self(i, ri);
                logIfMaster(format(" Final self Energy (%7s,%2d): %12.4f", residues[i], ri, self(i, ri)));
            }
            for (int i = 0; i < nResidues - 1; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues; j++) {
                    int rj = optimum[j];
                    sumPairEnergy += pair(i, ri, j, rj);
                    if (pair(i, ri, j, rj) > 10.0) {
                        logIfMaster(format(" Large Final Pair Energy (%7s,%2d) (%7s,%2d): %12.4f", residues[i], ri,
                                residues[j], rj, pair(i, ri, j, rj)));
                    }
                }
            }

            try {
                e = currentEnergy(residueList);
            } catch (ArithmeticException ex) {
                e = 1E100;
                logger.severe(String.format(" Exception %s in calculating current energy at the end of triples", ex.toString()));
            }
            logIfMaster(format(" Self Energy:   %16.8f", sumSelfEnergy));
            logIfMaster(format(" Pair Energy:   %16.8f", sumPairEnergy));

            double approximateEnergy = backboneEnergy + sumSelfEnergy + sumPairEnergy;

            if (threeBodyTerm) {
                for (int i = 0; i < nResidues - 2; i++) {
                    int ri = optimum[i];
                    for (int j = i + 1; j < nResidues - 1; j++) {
                        int rj = optimum[j];
                        for (int k = j + 1; k < nResidues; k++) {
                            int rk = optimum[k];
                            sumTrimerEnergy += triple(i, ri, j, rj, k, rk);
                        }
                    }
                }
                approximateEnergy += sumTrimerEnergy;
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - backboneEnergy;
                logIfMaster(format(" Trimer Energy: %16.8f", sumTrimerEnergy));
                logIfMaster(format(" Neglected:     %16.8f", higherOrderEnergy));
            } else {
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - backboneEnergy;
                logIfMaster(format(" Neglected:     %16.8f", higherOrderEnergy));
            }
            logIfMaster(format(" Approximate Energy:     %16.8f", approximateEnergy));
            return e;
        }

        /**
         * Permutations used only to set maximum bound on ensembleNumber, thus it is safe here to put that value in a 32-bit int.
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
            logger.info(format(" Requested an ensemble of %d, but only %d permutations exist; returning full ensemble", ensembleNumber, nPerms));
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

            logIfMaster(format(" Collecting permutations."));
            dryRun(residues, 0, currentRotamers);

            double pairTotalElimination = singletonPermutations - (double) evaluatedPermutations;
            currentEnsemble = (int) evaluatedPermutations;
            if (ensembleNumber == 1 && currentEnsemble == 0) {
                logger.severe(" No valid path through rotamer space found; try recomputing without pruning or using ensemble.");
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
                            logger.warning(format(" Versioning failed: appending to end of file %s", ensembleFile.getName()));
                        }
                    }
                    ensembleFilter = new PDBFilter(new File(ensembleFile.getName()), molecularAssembly, null, null);
                    logger.info(format(" Ensemble file: %s", ensembleFile.getName()));
                }
                logIfMaster(format(" Ensemble Search Stats: (buffer: %5.3f, current: %d, target: %d)", ensembleBuffer, currentEnsemble, ensembleNumber));
            }
            if (ensembleNumber == 1 || finalTry) {
                logIfMaster(format(" Number of permutations without DEE conditions: %10.4e.", permutations));
                logIfMaster(format(" Number of permutations after singleton eliminations: %10.4e.", singletonPermutations));
                logIfMaster(format(" Number of permutations removed by pairwise eliminations: %10.4e.", pairTotalElimination));
                logIfMaster(format(" Number of permutations remaining: %10.4e.", (double) evaluatedPermutations));
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

        double residueEnergy[] = new double[nResidues];

        double sumSelfEnergy = 0;
        double sumLowSelfEnergy = 0;
        for (int i = 0; i < nResidues; i++) {
            int ri = optimum[i];
            Residue residue = residues[i];
            Rotamer rotamers[] = residue.getRotamers(library);
            turnOnAtoms(residue);
            RotamerLibrary.applyRotamer(residue, rotamers[ri]);
            double self = self(i, ri);
            residueEnergy[i] = self;
            sumSelfEnergy += self;
            double lowest = lowestSelfEnergy(residues, i);
            sumLowSelfEnergy += lowest;
            if (self - lowest > 10.0) {
                logIfMaster(format(" Self-Energy (%7s,%2d): %12.4f (Lowest: %12.4f).",
                        residues[i], ri, self, lowest));
            }
        }

        double sumPairEnergy = 0.0;
        double sumLowPairEnergy = 0.0;
        double resPairEnergy[] = new double[nResidues];
        double lowPairEnergy[] = new double[nResidues];
        for (int i = 0; i < nResidues - 1; i++) {
            StringBuilder sb = new StringBuilder();
            int ri = optimum[i];
            double sumPairEnergyI = 0;
            double sumLowPairEnergyI = 0;
            for (int j = i + 1; j < nResidues; j++) {
                int rj = optimum[j];
                double pair = pair(i, ri, j, rj);
                residueEnergy[i] += 0.5 * pair;
                residueEnergy[j] += 0.5 * pair;
                sumPairEnergy += pair;
                sumPairEnergyI += pair;
                double lowest = lowestPairEnergy(residues, i, ri, j);
                sumLowPairEnergy += lowest;
                sumLowPairEnergyI += lowest;
                resPairEnergy[i] = 0.5 * pair;
                resPairEnergy[j] = 0.5 * pair;
                lowPairEnergy[i] = 0.5 * lowest;
                lowPairEnergy[j] = 0.5 * lowest;
                if (resPairEnergy[i] - lowPairEnergy[i] > 10.0) {
                    sb.append(format("  Pair Energy (%7s,%2d) (%7s,%2d): %12.4f (Lowest: %12.4f).\n",
                            residues[i], ri, residues[j], rj, pair, lowest));
                }
            }
            if (sumPairEnergyI - sumLowPairEnergyI > 10.0) {
                logIfMaster(format(" Total Pair Energy (%7s,%2d):        %12.4f (Lowest: %12.4f).",
                        residues[i], ri, sumPairEnergyI, sumLowPairEnergyI));
                sb.trimToSize();
                if (!sb.toString().isEmpty()) {
                    logIfMaster(sb.toString());
                }
            }
        }

        double e = 1E100;
        try {
            e = currentEnergy(residueList);
        } catch (ArithmeticException ex) {
            logger.severe(String.format(" Exception %s in calculating current energy at the end of self and pairs", ex.toString()));
        }
        logIfMaster(format(" Backbone Energy:    %16.8f.", backboneEnergy));
        logIfMaster(format(" Self Energy:        %16.8f (Lowest possible: %16.8f).", sumSelfEnergy, sumLowSelfEnergy));
        logIfMaster(format(" Pair Energy:        %16.8f (Lowest possible: %16.8f).", sumPairEnergy, sumLowPairEnergy));

        double approximateEnergy = backboneEnergy + sumSelfEnergy + sumPairEnergy;

        double sumTrimerEnergy = 0;
        if (threeBodyTerm) {
            for (int i = 0; i < nResidues - 2; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues - 1; j++) {
                    int rj = optimum[j];
                    for (int k = j + 1; k < nResidues; k++) {
                        int rk = optimum[k];
                        double triple = triple(i, ri, j, rj, k, rj);
                        residueEnergy[i] += triple / 3.0;
                        residueEnergy[j] += triple / 3.0;
                        residueEnergy[k] += triple / 3.0;
                        sumTrimerEnergy += triple;
                    }
                }
            }
            approximateEnergy += sumTrimerEnergy;
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - backboneEnergy;
            logIfMaster(format(" Trimer Energy:      %16.8f", sumTrimerEnergy));
            logIfMaster(format(" Neglected:          %16.8f", higherOrderEnergy));
        } else {
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - backboneEnergy;
            logIfMaster(format(" Neglected:          %16.8f", higherOrderEnergy));
        }

        logIfMaster(format(" Approximate Energy: %16.8f", approximateEnergy));

        logIfMaster("\n Final rotamers:");
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            int ri = optimum[i];
            Rotamer rotamer = rotamers[ri];
            logIfMaster(format(" %3d %c (%7s,%2d) %s %12.4f",
                    i + 1, residue.getChainID(), residue, ri, rotamer.toAngleString(), residueEnergy[i]));
            RotamerLibrary.applyRotamer(residue, rotamer);
        }

        return e;
    }

    /**
     * Return the lowest self-energy for residue i.
     *
     * @param residues
     * @param i
     * @return
     */
    private double lowestSelfEnergy(Residue residues[], int i) {
        if (residues == null) {
            return 0.0;
        }
        int n = residues.length;
        if (i < 0 || i >= n) {
            return 0.0;
        }
        Rotamer rotamers[] = residues[i].getRotamers(this.library);
        int nr = rotamers.length;
        double energy = Double.MAX_VALUE;
        for (int ni = 0; ni < nr; ni++) {
            try {
                double e = self(i, ni);
                if (e < energy) {
                    energy = e;
                }
            } catch (Exception e) {
                continue;
            }

        }
        return energy;
    }

    /**
     * Return the lowest pair-energy for residue (i,ri) with residue j.
     *
     * @param residues
     * @param i
     * @param ri
     * @param j
     * @return
     */
    private double lowestPairEnergy(Residue residues[], int i, int ri, int j) {
        if (residues == null) {
            return 0.0;
        }
        int n = residues.length;
        if (i < 0 || i >= n) {
            return 0.0;
        }
        if (j < 0 || j >= n) {
            return 0.0;
        }

        Rotamer rotamers[] = residues[j].getRotamers(this.library);
        int nr = rotamers.length;
        double energy = Double.MAX_VALUE;
        for (int jr = 0; jr < nr; jr++) {
            try {
                double e = pair(i, ri, j, jr);
                if (e < energy) {
                    energy = e;
                }
            } catch (Exception e) {
                continue;
            }
        }
        return energy;
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
    private double bruteForce(List<Residue> residueList) {

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

        logger.info(format(" Number of permutations: %16.8e.", permutations));
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

            double fullEnergy = 0;
            try {
                fullEnergy = currentEnergy(residueList);
            } catch (Exception ex) {
                logger.severe(String.format(" Exception %s in calculating full energy; FFX shutting down", ex.toString()));
            }

            logger.info(format(" Final summation of energies:    %16.5f", e));
            logger.info(format(" Final energy of optimized structure:    %16.5f", fullEnergy));
            logger.info(format(" Neglected:    %16.5f", fullEnergy - e));
        } else {
            e = rotamerOptimization(molecularAssembly, residues, 0, Double.MAX_VALUE, optimum);
        }

        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            int ri = optimum[i];
            Rotamer rotamer = rotamers[ri];
            logger.info(format(" %s %s (%d)", residue.getResidueNumber(), rotamer.toString(), ri));
            RotamerLibrary.applyRotamer(residue, rotamer);
            if (useFullAMOEBAEnergy) {
                try {
                    e = currentEnergy(residueList);
                } catch (ArithmeticException ex) {
                    logger.fine(String.format(" Exception %s in calculating full AMOEBA energy at the end of brute force", ex.toString()));
                }
            }
        }
        return e;
    }

    /**
     * Checks the distance matrix, finding the shortest distance between two
     * residues' rotamers under any symmetry operator; will evaluate this if
     * distance matrix not already filled.
     *
     * @param i  Residue i
     * @param ri Rotamer for i
     * @param j  Residue j
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
     * @param i  Residue i
     * @param ri Rotamer for i
     * @param j  Residue j
     * @param rj Rotamer for j
     * @param k  Residue k
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
     * @param i  Residue i
     * @param ri Rotamer for i
     * @param j  Residue j
     * @param rj Rotamer for j
     * @param k  Residue k
     * @param rk Rotamer for k
     * @param l  Residue l
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
     * @param i  Residue i
     * @param ri Rotamer for i
     * @param j  Residue j
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
            throw new NullPointerException(format(" Non-rotameric, non-forced residue present "
                    + "in residue list: %c %s-%d", residue.getChainID(), residue.toString(), residue.getResidueNumber()));
        }
    }

    private double slidingWindowOptimization(ArrayList<Residue> residueList, int windowSize, int increment, boolean revert,
                                             double distance, Direction direction) {

        long beginTime = -System.nanoTime();
        boolean incrementTruncated = false;
        boolean firstWindowSaved = false;
        int counter = 1;
        int windowEnd;

        int nOptimize = residueList.size();
        if (nOptimize < windowSize) {
            windowSize = nOptimize;
            logger.info(" Window size too small for given residue range; truncating window size.");
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
                    logIfMaster(format("\n Iteration %d of the sliding window.\n", counter++));
                    Residue firstResidue = residueList.get(windowStart);
                    Residue lastResidue = residueList.get(windowEnd);
                    if (firstResidue != lastResidue) {
                        logIfMaster(format(" Residues %s ... %s", firstResidue.toString(), lastResidue.toString()));
                    } else {
                        logIfMaster(format(" Residue %s", firstResidue.toString()));
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
                                                logIfMaster(format(" Adding residue %s at distance %16.8f Ang from %s %d.",
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
                            logIfMaster(format(" Adding nucleic acid residue 5' of window start %s to give it flexibility about its sugar pucker.",
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
                        double startingEnergy = 1E100;
                        try {
                            startingEnergy = currentEnergy(currentWindow);
                        } catch (ArithmeticException ex) {
                            logger.severe(String.format(" Exception %s in calculating starting energy of a window; FFX shutting down", ex.toString()));
                        }
                        if (useForcedResidues) {
                            if (onlyRotameric.size() < 1) {
                                logger.info(" Window has no rotameric residues.");
                                ResidueState.revertAllCoordinates(currentWindow, coordinates);
                            } else {
                                globalOptimization(onlyRotameric);
                                double finalEnergy = 1E100;
                                try {
                                    finalEnergy = currentEnergy(currentWindow);
                                } catch (ArithmeticException ex) {
                                    logger.severe(String.format(" Exception %s in calculating final energy of a window; FFX shutting down", ex.toString()));
                                }
                                if (startingEnergy <= finalEnergy) {
                                    logger.warning("Optimization did not yield a better energy. Reverting to orginal coordinates.");
                                    ResidueState.revertAllCoordinates(currentWindow, coordinates);
                                }
                            }
                        } else {
                            globalOptimization(currentWindow);
                            double finalEnergy = 1E100;
                            try {
                                finalEnergy = currentEnergy(currentWindow);
                            } catch (ArithmeticException ex) {
                                logger.severe(String.format(" Exception %s in calculating final energy of a window; FFX shutting down", ex.toString()));
                            }
                            if (startingEnergy <= finalEnergy) {
                                logger.warning("Optimization did not yield a better energy. Reverting to orginal coordinates.");
                                ResidueState.revertAllCoordinates(currentWindow, coordinates);
                            }
                        }
                    } else if (useForcedResidues) {
                        if (onlyRotameric.size() < 1) {
                            logger.info(" Window has no rotameric residues.");
                        } else {
                            globalOptimization(onlyRotameric);
                        }
                    } else {
                        globalOptimization(currentWindow);
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
                        //   StringBuilder header = new StringBuilder(format("Iteration %d of the sliding window\n", counter - 1));
                        try {
                            windowFilter.writeFile(file, false);
                            if (firstResidue != lastResidue) {
                                logger.info(format(" File with residues %s ... %s in window written to.", firstResidue.toString(), lastResidue.toString()));
                            } else {
                                logger.info(format(" File with residue %s in window written to.", firstResidue.toString()));
                            }
                        } catch (Exception e) {
                            logger.warning(format("Exception writing to file: %s", file.getName()));
                        }
                        firstWindowSaved = true;
                    }
                    long currentTime = System.nanoTime();
                    windowTime += currentTime;
                    logIfMaster(format(" Time elapsed for this iteration: %11.3f sec", windowTime * 1.0E-9));
                    logIfMaster(format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
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
            logger.info(format(" Manual superbox set over (minX, maxX, minY, "
                            + "maxY, minZ, maxZ): %f, %f, %f, %f, %f, %f",
                    minXYZ[0], maxXYZ[0], minXYZ[1], maxXYZ[1], minXYZ[2], maxXYZ[2]));
            logger.info(format(" Buffer size (included in dimensions): %f\n", superboxBuffer));
        } else { // Crystal systems will have already returned.
            logger.info(" System is aperiodic: protein box generated over these coordinates (minX, maxX, minY, maxY, minZ, maxZ):");
            String message = " Aperiodic box dimensions: ";
            for (int i = 0; i < minXYZ.length; i++) {
                message = message.concat(format("%f,%f,", minXYZ[i], maxXYZ[i]));
            }
            message = message.substring(0, message.length() - 1);
            logger.info(message);
            logger.info(format(" Buffer size (included in dimensions): %f\n", superboxBuffer));
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
            logger.severe(format(" FFX shutting down: Box optimization start is out of range of total boxes: %d > %d", (boxStart + 1), totalCells));
        }
        if (boxEnd > totalCells - 1) {
            boxEnd = totalCells - 1;
            logIfMaster(" Final box out of range: reset to last possible box.");
        } else if (boxEnd < 0) {
            boxEnd = totalCells - 1;
        }
        BoxOptCell[] cells = loadCells(crystal, residues);
        int numCells = cells.length;
        logIfMaster(format(" Optimizing boxes  %d  to  %d", (boxStart + 1), (boxEnd + 1)));
        /*int restartCell = -1;
         if (parallelEnergies && loadEnergyRestart) {
         restartCell = loadEnergyRestartIterations(energyRestartFile);
         if (restartCell > -1) {
         logIfMaster(format(" Optimization restarting from file at cell #%d", restartCell));
         }
         }*/
        for (int i = 0; i < numCells; i++) {
            /*if (restartCell > -1 && i < restartCell) {
             continue;
             }*/
            BoxOptCell celli = cells[i];
            ArrayList<Residue> residuesList = celli.getResiduesAsList();
            int[] cellIndices = celli.getXYZIndex();
            logIfMaster(format("\n Iteration %d of the box optimization.", (i + 1)));
            logIfMaster(format(" Cell index (linear): %d", (celli.getLinearIndex() + 1)));
            logIfMaster(format(" Cell xyz indices: x = %d, y = %d, z = %d", cellIndices[0] + 1, cellIndices[1] + 1, cellIndices[2] + 1));
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
                    logIfMaster(format(" Residues %s ... %s", firstResidue.toString(), lastResidue.toString()));
                } else {
                    logIfMaster(format(" Residue %s", firstResidue.toString()));
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

                    double startingEnergy = 0;
                    double finalEnergy = 0;
                    try {
                        startingEnergy = currentEnergy(residuesList);
                    } catch (ArithmeticException ex) {
                        logger.severe(String.format(" Exception %s in calculating starting energy of a box; FFX shutting down", ex.toString()));
                    }
                    globalOptimization(residuesList);
                    try {
                        finalEnergy = currentEnergy(residuesList);
                    } catch (ArithmeticException ex) {
                        logger.severe(String.format(" Exception %s in calculating starting energy of a box; FFX shutting down", ex.toString()));
                    }

                    if (startingEnergy <= finalEnergy) {
                        logger.warning("Optimization did not yield a better energy. Reverting to orginal coordinates.");
                        ResidueState.revertAllCoordinates(residuesList, coordinates);
                    }
                    long currentTime = System.nanoTime();
                    boxTime += currentTime;
                    logIfMaster(format(" Time elapsed for this iteration: %11.3f sec", boxTime * 1.0E-9));
                    logIfMaster(format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
                } else {
                    globalOptimization(residuesList);
                    long currentTime = System.nanoTime();
                    boxTime += currentTime;
                    logIfMaster(format(" Time elapsed for this iteration: %11.3f sec", boxTime * 1.0E-9));
                    logIfMaster(format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
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
                            logIfMaster(format(" File with residues %s ... %s in window written.", firstResidue.toString(), lastResidue.toString()));
                        } else {
                            logIfMaster(format(" File with residue %s in window written.", firstResidue.toString()));
                        }
                        firstCellSaved = true;
                    } catch (Exception e) {
                        logger.warning(format("Exception writing to file: %s", file.getName()));
                    }
                }
                /*for (Residue residue : residueList) {
                    if (residue instanceof MultiResidue) {
                        ((MultiResidue) residue).setDefaultResidue();
                        residue.reInitOriginalAtomList();
                    }
                }*/
            } else {
                logIfMaster(format(" Empty box: no residues found."));
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
     * @param crystal  Polymer crystal or dummy crystal
     * @param residues All residues to be optimized
     * @return Filled cells.
     */
    private BoxOptCell[] loadCells(Crystal crystal, Residue[] residues) {
        double xCellBorderFracSize = (boxBorderSize / crystal.a);
        double yCellBorderFracSize = (boxBorderSize / crystal.b);
        double zCellBorderFracSize = (boxBorderSize / crystal.c);
        logIfMaster(format(" Number of boxes along x: %d, y: %d, z: %d", numXYZBoxes[0], numXYZBoxes[1], numXYZBoxes[2]));

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
     * @param crystal  Crystal group.
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
                    logIfMaster(format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
                    // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                    while (goldsteinDriver(residues)) {
                        i++;
                        logIfMaster(this.toString());
                        logIfMaster(format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
                    }
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                    // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                    pairEliminated = false;
                    while (goldsteinPairDriver(residues)) {
                        pairEliminated = true;
                        i++;
                        logIfMaster(this.toString());
                        logIfMaster(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                    }
                } else {
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Single DEE conditions ", i));
                    // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                    while (deeRotamerElimination(residues)) {
                        i++;
                        logIfMaster(toString());
                        logIfMaster(format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
                    }
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                    // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                    pairEliminated = false;
                    while (deeRotamerPairElimination(residues)) {
                        pairEliminated = true;
                        i++;
                        logIfMaster(toString());
                        logIfMaster(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
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
                try {
                    logIfMaster(format("\n Beginning Energy %16.8f", currentEnergy(residues)));
                } catch (ArithmeticException ex) {
                    logger.severe(String.format(" Exception %s in calculating beginning energy; FFX shutting down.", ex.toString()));
                }
            }
            // eliminateNABackboneRotamers checks for threshold != 0 internally.
            eliminateNABackboneRotamers(residues, numEliminatedRotamers);
            /* double naBackboneEnergy = currentEnergy();
             logger.info(format(" After Eliminate NA Backbone %16.8f", naBackboneEnergy)); */
        } else if (verboseEnergies) {
            try {
                logIfMaster(format("\n Beginning Energy %16.8f", currentEnergy(residues)));
            } catch (ArithmeticException ex) {
                logger.severe(String.format(" Exception %s in calculating beginning energy; FFX shutting down.", ex.toString()));
            }
        }

        rotamerEnergies(residues);

        // testSelfEnergyElimination(residues);

        // testPairEnergyElimination(residues, 19);

        // Beginning energy
        int currentRotamers[] = new int[residues.length];

        if (pruneClashes) {
            validateDEE(residues);
        }

        int i = 0;
        boolean pairEliminated = true;
        while (pairEliminated) {
            if (useGoldstein) {
                pairEliminated = false;
                i++;
                logIfMaster(format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
                // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                while (goldsteinDriver(residues)) {
                    i++;
                    logIfMaster(this.toString());
                    logIfMaster(format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
                }
                i++;
                logIfMaster(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                while (goldsteinPairDriver(residues)) {
                    pairEliminated = true;
                    i++;
                    logIfMaster(this.toString());
                    logIfMaster(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                }
            } else {
                i++;
                logIfMaster(format("\n Iteration %d: Applying Single DEE conditions ", i));
                // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                while (deeRotamerElimination(residues)) {
                    i++;
                    logIfMaster(toString());
                    logIfMaster(format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
                }
                i++;
                logIfMaster(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                pairEliminated = false;
                while (deeRotamerPairElimination(residues)) {
                    pairEliminated = true;
                    i++;
                    logIfMaster(toString());
                    logIfMaster(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
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

    /**
     * Calculates the energy at the current state.
     *
     * @param resArray Array of residues in current energy term.
     * @return Energy of the current state.
     */
    private double currentEnergy(Residue[] resArray) throws ArithmeticException {
        return currentEnergy(Arrays.asList(resArray));
    }

    /**
     * Calculates the energy at the current state.
     *
     * @param resList List of residues in current energy term.
     * @return Energy of the current state.
     */
    private double currentEnergy(List<Residue> resList) throws ArithmeticException {
        List<Rotamer> rots = resList.stream().
                filter(res -> res != null ).
                map(Residue::getRotamer).
                collect(Collectors.toList());
        File energyDir = dirSupplier.apply(resList, rots);
        return eFunction.applyAsDouble(energyDir);
    }

    /**
     * Default method for obtaining energy: calculates energy of the Potential.
     *
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
    private double currentEnergyWrapper(List<Residue> resList) throws ArithmeticException {
        return currentEnergy(resList);
    }

    /**
     * Sets the function used to calculate energy.
     *
     * @param ef File to energy
     */
    public void setEnergyFunction(ToDoubleFunction<File> ef) {
        this.eFunction = ef;
    }

    /**
     * Sets the directory provider used.
     *
     * @param dirProvider A function of residue list to appropriate directory
     */
    public void setDirectoryProvider(BiFunction<List<Residue>, List<Rotamer>, File> dirProvider) {
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

            long energyStartTime = System.nanoTime();
            SinglesEnergyRegion singlesRegion = new SinglesEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            PairsEnergyRegion pairsRegion = new PairsEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            TriplesEnergyRegion triplesRegion = new TriplesEnergyRegion(energyWorkerTeam.getThreadCount(), residues);
            QuadsEnergyRegion quadsRegion = new QuadsEnergyRegion(energyWorkerTeam.getThreadCount(), residues);

            try {
                if (loaded < 1) {
                    singlesMap.clear();
                    // allocate selfEnergy
                    int singleJobIndex = 0;
                    selfEnergy = new double[nResidues][];
                    for (int i = 0; i < nResidues; i++) {
                        Residue resi = residues[i];
                        Rotamer roti[] = resi.getRotamers(library);
                        selfEnergy[i] = new double[roti.length];
                        for (int ri = 0; ri < roti.length; ri++) {
                            if (!check(i, ri)) {
                                Integer selfJob[] = {i, ri};
                                if (decomposeOriginal && ri != 0) {
                                    continue;
                                }
                                singlesMap.put(singleJobIndex++, selfJob);
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
                long singlesTime = System.nanoTime() - energyStartTime;
                logIfMaster(format(" Time for single energies: %12.4g", (singlesTime * 1.0E-9)));

                if (loaded < 2) {
                    pairsMap.clear();
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
                                    if (checkToJ(i, ri, j, rj)) {
                                        continue;
                                    }
                                    Integer pairJob[] = {i, ri, j, rj};
                                    if (decomposeOriginal && (ri != 0 || rj != 0)) {
                                        continue;
                                    }
                                    pairsMap.put(pairJobIndex++, pairJob);
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
                long pairsTime = System.nanoTime() - (singlesTime + energyStartTime);
                long triplesTime = 0;
                long quadsTime = 0;
                logIfMaster(format(" Time for pair energies:   %12.4g", (pairsTime * 1.0E-9)));

                if (threeBodyTerm) {
                    if (loaded < 3) {
                        trimersMap.clear();
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
                                        if (checkToJ(i, ri, j, rj)) {
                                            continue;
                                        }
                                        threeBodyEnergy[i][ri][j][rj] = new double[nResidues][];
                                        for (int k = j + 1; k < nResidues; k++) {
                                            Residue resk = residues[k];
                                            Rotamer rotk[] = resk.getRotamers(library);
                                            threeBodyEnergy[i][ri][j][rj][k] = new double[rotk.length];
                                            for (int rk = 0; rk < rotk.length; rk++) {
                                                if (checkToK(i, ri, j, rj, k, rk)) {
                                                    continue;
                                                }
                                                Integer trimerJob[] = {i, ri, j, rj, k, rk};
                                                if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0)) {
                                                    continue;
                                                }
                                                trimersMap.put(trimerJobIndex++, trimerJob);
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
                    triplesTime = System.nanoTime() - (pairsTime + singlesTime + energyStartTime);
                    logIfMaster(format(" Time for triple energies: %12.4g", (triplesTime * 1.0E-9)));
                }

                if (computeQuads) {
                    logger.info(" Creating quad jobs...");
                    quadsMap.clear();
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
                                    if (checkToJ(i, ri, j, rj)) {
                                        continue;
                                    }
                                    for (int k = j + 1; k < nResidues; k++) {
                                        Residue resk = residues[k];
                                        Rotamer rotk[] = resk.getRotamers(library);
                                        for (int rk = 0; rk < rotk.length; rk++) {
                                            /*if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk) || check(i, ri, j, rj, k, rk)) {
                                                continue;
                                            }*/
                                            if (checkToK(i, ri, j, rj, k, rk)) {
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
                                                    if (checkToL(i, ri, j, rj, k, rk, l, rl)) {
                                                        continue;
                                                    }
                                                    Integer quadJob[] = {i, ri, j, rj, k, rk, l, rl};
                                                    if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0 || rl != 0)) {
                                                        continue;
                                                    }
                                                    quadsMap.put(quadJobIndex++, quadJob);
                                                    if (quadsMap.size() > quadMaxout) {
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
//                    logger.info(format(" Proc %d broadcasting ready for quads.", world.rank()));
                    multicastBuf(thisProcReadyBuf);
                    // launch parallel threeBody calculation
                    int waiting = 0;
                    while (!readyForQuads) {
                        Thread.sleep(POLLING_FREQUENCY);
                    }

                    logger.info(format(" Running quads: %d jobs.", quadsMap.size()));
                    energyWorkerTeam.execute(quadsRegion);
                    quadsTime = System.nanoTime() - (triplesTime + pairsTime + singlesTime + energyStartTime);
                    logIfMaster(format(" Time for quad energies:   %12.4g", quadsTime * 1.0E-9));
                }
                long allTime = singlesTime + pairsTime + triplesTime + quadsTime;
                logIfMaster(format(" Time for all energies:    %12.4g", allTime * 1.0E-9));
            } catch (Exception ex) {
                String message = " Exception computing rotamer energies in parallel.";
                logger.log(Level.SEVERE, message, ex);
            } finally {
                // stop receive and energyWriter threads
                // receiveThread.die();
            }

            // Turn on all atoms.
            for (int i = 0; i < atoms.length; i++) {
                atoms[i].setUse(true);
            }
            // Print the energy with all rotamers in their default conformation.
            if (verboseEnergies && master) {
                try {
                    double defaultEnergy = currentEnergy(residues);
                    logger.info(format(" Energy of the system with rotamers in their default conformation: %16.8f", defaultEnergy));
                } catch (ArithmeticException ex) {
                    logger.severe(String.format(" Exception %s in calculating default energy; FFX shutting down", ex.toString()));
                }
            }
            return backboneEnergy;
        }   // endif(parallelEnergies)

        // TODO: deprecate the rest of the rotamerEnergies() method; the parallel form should work (better) in all scenarios
        // Store coordinates: will be useful for nucleic acids.
        ArrayList<Residue> currentResidueList = new ArrayList<>();
        currentResidueList.addAll(Arrays.asList(residues));
        ResidueState[] origCoordinates = ResidueState.storeAllCoordinates(currentResidueList);

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

        /**
         * Compute the backbone energy.
         */
        List<Residue> rList = new ArrayList<>(Collections.nCopies(4, null));
        try {
            backboneEnergy = currentEnergy(rList);
        } catch (ArithmeticException ex) {
            logger.severe(format(" Error in calculation of backbone energy %s", ex.getMessage()));
        }
        logger.info(format(" Backbone energy:  %16.8f\n", backboneEnergy));

        /**
         * Compute the self-energy for each rotamer of each residue.
         */
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
                long time = -System.nanoTime();
                // This entire section of code is deprecated.
                selfEnergy[i][ri] = currentEnergy(rList) - backboneEnergy;
                time += System.nanoTime();
                logger.info(format(" Self-energy %s %d: %16.8f in %6.4 sec",
                        residue, ri, self(i, ri), time * 1.0e-9));
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
                            logger.info(format(" Clash for rotamer %s %d: %16.8f >> %16.8f",
                                    residue, ri, self(i, ri), minEnergy));
                            eliminateRotamer(residues, i, ri, print);
                        }
                    } else if (!check(i, ri) && (self(i, ri) > (minEnergy + clashThreshold))) {
                        // don't prune orig-coords rotamers
                        if (!(useOrigCoordsRot && ri == 0)) {
                            logger.info(format(" Clash for rotamer %s %d: %16.8f >> %16.8f",
                                    residue, ri, self(i, ri), minEnergy));
                            eliminateRotamer(residues, i, ri, print);
                        }
                    }
                }
            }

            // Reset the residue to its default conformation.
            if (residue.getResidueType() == NA) {
                residue.revertState(origCoordinates[i]);
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
                                logger.info(format(" Pair energy %s %d, %s %d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang",
                                        residuei, ri, residuej, rj, dist, superpositionThreshold));
                            } else {
                                long time = -System.nanoTime();
                                twoBodyEnergy[i][ri][j][rj] = currentEnergy(rList) - self(i, ri) - self(j, rj) - backboneEnergy;
                                time += System.nanoTime();
                                logger.info(format(" Pair energy %s %d, %s %d: %16.8f at %10.3f Ang in %5.3f sec",
                                        residuei, ri, residuej, rj, pair(i, ri, j, rj), dist, time * 1.0e-9));
                            }
                        } else {
                            long time = -System.nanoTime();
                            twoBodyEnergy[i][ri][j][rj] = currentEnergy(rList)
                                    - self(i, ri) - self(j, rj) - backboneEnergy;
                            time += System.nanoTime();
                            logger.info(format(" Pair energy %s %d, %s %d: %16.8f in %5.3f sec.",
                                    residuei, ri, residuej, rj, pair(i, ri, j, rj), time * 1.0e-9));
                        }
                        if (algorithmListener != null) {
                            algorithmListener.algorithmUpdate(molecularAssembly);
                        }
                    }
                    if (residuej.getResidueType() == NA) {
                        residuej.revertState(origCoordinates[j]);
                    } else {
                        RotamerLibrary.applyRotamer(residuej, rotamersj[0]);
                    }
                    turnOffAtoms(residuej);
                }
            }

            /**
             * Reset the residue to its default conformation.
             */
            if (residuei.getResidueType() == NA) {
                residuei.revertState(origCoordinates[i]);
            } else {
                RotamerLibrary.applyRotamer(residuei, rotamersi[0]);
            }
            turnOffAtoms(residuei);
        }

        if (prunePairClashes) {
            prunePairClashes(residues);
        }

        /**
         * Compute the trimer-energy for each pair of rotamers.
         */
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
                    threeBodyEnergy[i][ri] = new double[nResidues][][][];
                    for (int j = i + 1; j < nResidues - 1; j++) {
                        Residue residuej = residues[j];
                        rList.set(1, residuej);
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        turnOnAtoms(residuej);
                        // Loop over residue j's rotamers.
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
                                    if ((pruneClashes && (check(k, rk))
                                            || (prunePairClashes && check(i, ri, k, rk))
                                            || (prunePairClashes && check(j, rj, k, rk)))
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
                                    } else {
                                        dij = localDistanceMatrix[i][ri][j][rj];
                                        dik = localDistanceMatrix[i][ri][k][rk];
                                        djk = localDistanceMatrix[j][rj][k][rk];
                                    }

                                    double dist = Math.min(dij, Math.min(dik, djk));
                                    if (!threeBodyCutoff || (dist < threeBodyCutoffDist)) {
                                        if (dist < superpositionThreshold) {
                                            threeBodyEnergy[i][ri][j][rj][k][rk] = 1.0E100;
                                            logger.info(format(
                                                    " Trimer energy %s %d, %s %d, %s %d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang.",
                                                    residuei, ri, residuej, rj, residuek, rk, dist, superpositionThreshold));
                                        } else {
                                            Rotamer rotamerk = rotamersk[rk];
                                            RotamerLibrary.applyRotamer(residuek, rotamerk);
                                            long time = -System.nanoTime();
                                            threeBodyEnergy[i][ri][j][rj][k][rk] = currentEnergy(rList)
                                                    - self(i, ri) - self(j, rj) - self(k, rk)
                                                    - pair(i, ri, j, rj) - pair(i, ri, k, rk)
                                                    - pair(j, rj, k, rk) - backboneEnergy;
                                            time += System.nanoTime();
                                            logger.info(format(
                                                    " Trimer energy %s %d, %s %d, %s %d: %16.8f at %10.3f Ang in %6.4f.",
                                                    residuei, ri, residuej, rj, residuek, rk,
                                                    threeBodyEnergy[i][ri][j][rj][k][rk], dist, time * 1.0e-9));
                                            if (algorithmListener != null) {
                                                algorithmListener.algorithmUpdate(molecularAssembly);
                                            }
                                        }
                                    } else {
                                        logger.info(format(
                                                " Trimer energy %s %d, %s %d, %s %d: %16.8f at %10.3f Ang.",
                                                residuei, ri, residuej, rj, residuek, rk, 0.0,
                                                threeBodyCutoffDist));
                                    }
                                }
                                // Reset the residue to its default conformation.
                                if (residuek.getResidueType() == NA) {
                                    residuek.revertState(origCoordinates[k]);
                                } else {
                                    RotamerLibrary.applyRotamer(residuek, rotamersk[0]);
                                }
                                turnOffAtoms(residuek);
                            }
                        }
                        // Reset the residue to its default conformation.
                        if (residuej.getResidueType() == NA) {
                            residuej.revertState(origCoordinates[j]);
                        } else {
                            RotamerLibrary.applyRotamer(residuej, rotamersj[0]);
                        }
                        turnOffAtoms(residuej);
                    }
                }
                // Reset the residue to its default conformation.
                if (residuei.getResidueType() == NA) {
                    residuei.revertState(origCoordinates[i]);
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
            try {
                double defaultEnergy = currentEnergy(residues);
                logger.info(format(" Energy of the system with rotamers in their default conformation: %16.8f", defaultEnergy));
            } catch (ArithmeticException ex) {
                logger.severe(String.format(" Exception %s in calculating default energy; FFX shutting down", ex.toString()));
            }
        }

        return backboneEnergy;
    }

    private double computeBackboneEnergy(Residue[] residues) throws ArithmeticException {
        if (residues == null) {
            return 0.0;
        }
        turnOffAllResidues(residues);
        // Checked above this level.
        return currentEnergy(residues);
    }

    private double computeSelfEnergy(Residue[] residues, int i, int ri) {
        if (residues == null) {
            return 0.0;
        }
        int nRes = residues.length;
        if (i < 0 || i >= nRes) {
            return 0.0;
        }
        Residue residue = residues[i];
        Rotamer rotamers[] = residue.getRotamers(library);
        int ni = rotamers.length;
        if (ri < 0 || ri >= ni) {
            return 0.0;
        }
        turnOffAllResidues(residues);
        turnOnAtoms(residue);
        RotamerLibrary.applyRotamer(residue, rotamers[ri]);

        // Unchecked call; would check above this level if this method was used.
        double energy = currentEnergy(residues) - backboneEnergy;

        RotamerLibrary.applyRotamer(residue, rotamers[0]);
        turnOffAtoms(residue);

        return energy;
    }

    private double computePairEnergy(Residue[] residues, int i, int ri, int j, int rj) {
        if (residues == null) {
            return 0.0;
        }
        int nRes = residues.length;
        if (i < 0 || i >= nRes) {
            return 0.0;
        }
        Residue residueI = residues[i];
        Rotamer rotamersI[] = residueI.getRotamers(library);
        int ni = rotamersI.length;
        if (ri < 0 || ri >= ni) {
            return 0.0;
        }
        Residue residueJ = residues[j];
        Rotamer rotamersJ[] = residueJ.getRotamers(library);
        int nj = rotamersJ.length;
        if (rj < 0 || rj >= nj) {
            return 0.0;
        }

        turnOffAllResidues(residues);

        turnOnAtoms(residueI);
        turnOnAtoms(residueJ);
        RotamerLibrary.applyRotamer(residueI, rotamersI[ri]);
        RotamerLibrary.applyRotamer(residueJ, rotamersJ[rj]);

        // Unchecked call; would check above this level if this method was used.
        double energy = currentEnergy(residues) - backboneEnergy - self(i, ri) - self(j, rj);

        turnOffAtoms(residueI);
        turnOffAtoms(residueJ);
        RotamerLibrary.applyRotamer(residueI, rotamersI[0]);
        RotamerLibrary.applyRotamer(residueJ, rotamersJ[0]);

        return energy;
    }

    private double computeTripleEnergy(Residue[] residues, int i, int ri, int j, int rj, int k, int rk) {
        if (residues == null) {
            return 0.0;
        }
        int nRes = residues.length;
        if (i < 0 || i >= nRes) {
            return 0.0;
        }
        Residue residueI = residues[i];
        Rotamer rotamersI[] = residueI.getRotamers(library);
        int ni = rotamersI.length;
        if (ri < 0 || ri >= ni) {
            return 0.0;
        }

        Residue residueJ = residues[j];
        Rotamer rotamersJ[] = residueJ.getRotamers(library);
        int nj = rotamersJ.length;
        if (rj < 0 || rj >= nj) {
            return 0.0;
        }
        Residue residueK = residues[k];
        Rotamer rotamersK[] = residueK.getRotamers(library);
        int nk = rotamersK.length;
        if (rk < 0 || rk >= nk) {
            return 0.0;
        }

        turnOffAllResidues(residues);
        turnOnAtoms(residueI);
        turnOnAtoms(residueJ);
        turnOnAtoms(residueK);
        RotamerLibrary.applyRotamer(residueI, rotamersI[ri]);
        RotamerLibrary.applyRotamer(residueJ, rotamersJ[rj]);
        RotamerLibrary.applyRotamer(residueK, rotamersK[rk]);

        // Unchecked call; would check above this level if this method was used.
        double energy = currentEnergy(residues) - backboneEnergy - self(i, ri) - self(j, rj) - self(k, rk) -
                pair(i, ri, j, rj) - pair(i, ri, k, rk) - pair(j, rj, k, rk);

        turnOffAtoms(residueI);
        turnOffAtoms(residueJ);
        turnOffAtoms(residueK);
        RotamerLibrary.applyRotamer(residueI, rotamersI[0]);
        RotamerLibrary.applyRotamer(residueJ, rotamersJ[0]);
        RotamerLibrary.applyRotamer(residueK, rotamersK[0]);

        return energy;
    }


    private void turnOffAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (int i = 0; i < nRes; i++) {
            turnOffAtoms(residues[i]);
        }
    }

    private void turnOnAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (int i = 0; i < nRes; i++) {
            turnOnAtoms(residues[i]);
        }
    }

    /**
     * Calculates the distance matrix; residue-residue distance is defined as
     * the shortest atom-atom distance in any possible rotamer-rotamer pair if
     * the residues are neighbors (central atom-central atom distances are
     * within a cutoff); else, distances are set to a default of
     * Double.MAX_VALUE.
     * <p>
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
                    logIfMaster(format(" Non-forced Residue i %s has null rotamers.",
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
                            logIfMaster(format(" Residue j %s has null rotamers.", residuej.toString()));
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

        logger.info(format(" Number of pairwise distances: %d", numDistances));

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
            logger.info(format(" Built residue neighbor list:           %8.3f sec", neighborTime * 1.0e-9));

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
            logger.info(format(" Pairwise distance matrix:              %8.3f sec\n", parallelTime * 1.0e-9));

            ResidueState.revertAllCoordinates(allResiduesList, orig);
            try {
                parallelTeam.shutdown();
            } catch (Exception ex) {
                logger.warning(format(" Exception shutting down parallel team for the distance matrix: %s", ex.toString()));
            }
        }
    }

    /**
     * Fills a local distance matrix based on residues passed to the method,
     * using only their non-eliminated Rotamers.
     *
     * @param residues            Residues to create a distance matrix for
     * @param localDistanceMatrix Array to be filled.
     */
    private void localSequentialDistanceMatrix(Residue[] residues, double[][][][] localDistanceMatrix) {
        Crystal crystal = molecularAssembly.getCrystal();
        ArrayList<Residue> residuesList = new ArrayList<>();
        residuesList.addAll(Arrays.asList(residues));
        ResidueState[] originalCoordinates = ResidueState.storeAllCoordinates(residuesList);

        int nSymm = crystal.spaceGroup.getNumberOfSymOps();
        int nResidues = residues.length;
        /**
         * The flag isDistant is used to skip over the ri loop if r-j distance is above the
         * 25 A center-to-center cutoff (ri, rj could not possibly be within 9 A of each other).
         */
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
            /**
             * Loop over residues i.  Testing showed no i-j vs. j-i discrepancies over symmetry operators.
             */
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
                        /**
                         * Use reference atoms for a quick distance check, before looping over rj.
                         */
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

    /**
     * Calculates the minimum distance between two sets of coordinates in a
     * given symmetry operator.
     *
     * @param resi  Coordinates of i by [atom][xyz]
     * @param resj  Coordinates of j by [atom][xyz]
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
            logIfMaster(format(" %3d Residue %c %7s with %2d rotamers.", i + 1, residuei.getChainID(), residuei.toString(), lenri));
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
                minEnergySingles[ri] = self(i, ri);
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
                        logger.info(format(" Inconsistent Pair: %7s %2d, %7s.", residuei, ri, residuej));
                        //eliminateRotamer(residues, i, ri, print);
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
                                logger.info(format(" Inconsistent triple: %s %d, %s %d, %s.",
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
                            logger.info(format(" Eliminating rotamer pair: %s %d, %s %d (%16.8f > %16.8f + %6.6f)",
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

        // Loop over the 2nd residues' rotamers.
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
                    logger.info(format(" Inconsistent Pair: %7s %2d, %7s %2d.", residuei, ri, residuej, rj));
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
//                    logger.info(format(" Eliminating Inconsistent Pair: %s %d, %s %d.",
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
     * @param minMax   The bound on the 3-body energy (minMax[0] = min, minMax[1]
     *                 = max.
     * @param i        Residue i
     * @param ri       Rotamer ri of Residue i
     * @param j        Residue j
     * @param rj       Rotamer rj of Residue j
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
            for (int rk = 0; rk < lenrk; rk++) {
                double current = triple(i, ri, j, rj, k, rk);
                if (current < currentMin) {
                    currentMin = current;
                }
                if (current > currentMax) {
                    currentMax = current;
                }
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
     * @param minMax   The bound on the 3-body energy (minMax[0] = min, minMax[1]
     *                 = max.
     * @param i        Residue i
     * @param ri       Rotamer ri of Residue i
     * @param j        Residue j
     * @param rj       Rotamer rj of Residue j
     * @param k        Residue k
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
                valid = true;

                // (TODO - quad returns 0.0 now)
                // double quadEnergy = quad(i, ri, j, rj, k, rk, l, rl);
                double quadEnergy = 0.0;

                // Collect the 3-body energies and 4-body energy
                double current = triple(i, ri, k, rk, l, rl) + triple(j, rj, k, rk, l, rl) + quadEnergy;

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

    private boolean goldsteinDriver(Residue[] residues) {
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
                // Continue if rotamer (i, riA) is not valid.
                if (check(i, riA)) {
                    continue;
                }
                for (int riB = 0; riB < nri; riB++) {
                    // The eliminating rotamer cannot be riA and must be a valid.
                    if (riA == riB || check(i, riB)) {
                        continue;
                    }
                    if (goldsteinElimination(residues, i, riA, riB)) {
                        eliminated = true;
                        break;
                    }
                }
            }
        }
        return eliminated;
    }

    private boolean goldsteinElimination(Residue residues[], int i, int riA, int riB) {
        int nres = residues.length;
        Residue resi = residues[i];

        // Initialize Goldstein inequality.
        double selfDiff = self(i, riA) - self(i, riB);
        double goldsteinEnergy = selfDiff;

        double sumPairDiff = 0.0;
        double sumTripleDiff = 0.0;

        // Loop over a 2nd residue j.
        for (int j = 0; j < nres; j++) {
            if (j == i) {
                continue;
            }
            Residue resj = residues[j];
            Rotamer rotj[] = resj.getRotamers(library);
            int nrj = rotj.length;
            double minForResJ = Double.MAX_VALUE;
            double minPairDiff = 0.0;
            double minTripleDiff = 0.0;
            int rjEvals = 0;
            // Loop over the rotamers for residue j.
            for (int rj = 0; rj < nrj; rj++) {
                if (check(j, rj)) {
                    continue;
                }

                // This following is NOT checked in Osprey.
                // if (check(i, riA, j, rj) ||  check(i, riB, j, rj)) {
                //    continue;
                // }


                double pairI = pair(i, riA, j, rj);
                double pairJ = pair(i, riB, j, rj);
                double pairDiff = pair(i, riA, j, rj) - pair(i, riB, j, rj);

                /**
                 if (abs(pairI) > 1.0e4 && abs(pairJ) > 1.0e4) {
                 logIfMaster(format(" Goldstein difference [(%7s,%2d),(%7s-%2d)] %12.4f - [(%7s-%2d),(%7s-%2d)] %12.4f = %12.4f",
                 resi, riA, resj, rj, pairI, resi, riB, resj, rj, pairJ, pairDiff));
                 } */

                rjEvals++;

                double tripleDiff = 0.0;
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
                            if (check(k, rk) || check(j, rj, k, rk) || check(i, riA, k, rk)) {
                                continue;
                            }
                            // By anaology, the following is NOT checked in Osprey.
                            // if (check(i, riB, k, rk)) {
                            //     continue;
                            // }

                            rkEvals++;
                            double e = triple(i, riA, j, rj, k, rk) - triple(i, riB, j, rj, k, rk);
                            if (e < minForResK) {
                                minForResK = e;
                            }
                        }
                        /**
                         * If there were no 3-body interactions with residue k, then minForResk is zero.
                         */
                        if (rkEvals == 0) {
                            minForResK = 0.0;
                        }
                        tripleDiff += minForResK;
                    }
                }
                double sumOverK = pairDiff + tripleDiff;
                if (sumOverK < minForResJ) {
                    minForResJ = sumOverK;
                    minPairDiff = pairDiff;
                    minTripleDiff = tripleDiff;
                }
            }
            // If there are no 2-body interactions, then minForResJ is zero.
            if (rjEvals == 0) {
                minForResJ = 0.0;
                minPairDiff = 0.0;
                minTripleDiff = 0.0;
            }
            sumPairDiff += minPairDiff;
            sumTripleDiff += minTripleDiff;
            goldsteinEnergy += minForResJ;
        }

        if (goldsteinEnergy > ensembleBuffer) {
            if (eliminateRotamer(residues, i, riA, print)) {
                logIfMaster(format("  Rotamer elimination of (%7s,%2d) by (%7s,%2d): %12.4f > %6.4f.",
                        resi, riA, resi, riB, goldsteinEnergy, ensembleBuffer));
                logIfMaster(format("   Self: %12.4f, Pairs: %12.4f, Triples: %12.4f.",
                        selfDiff, sumPairDiff, sumTripleDiff));
                return true;
            }
        }
/*        else {
            if (i == 80) {
                logIfMaster(format("  Rotamer (%7s,%2d) not eliminated by (%7s,%2d): %12.4f < %6.4f.",
                        resi, riA, resi, riB, goldsteinEnergy, ensembleBuffer));
                logIfMaster(format("   Self: %12.4f, Pairs: %12.4f, Triples: %12.4f.",
                        selfDiff, sumPairDiff, sumTripleDiff));
            }
        }*/

        return false;
    }

    /**
     * The Goldstein Rotamer Pair driver routine generates pairs of rotamers to
     * be evaluated by the 3-body form of the Goldstein criteria.
     *
     * @param residues The list of residues to be optimized.
     * @return true if a residue is eliminated.
     */
    private boolean goldsteinPairDriver(Residue[] residues) {
        int nres = residues.length;
        // A flag to indicate if any more rotamers or rotamer pairs were eliminated.
        boolean eliminated = false;
        // Loop over residue i.
        for (int i = 0; i < nres; i++) {
            Residue resi = residues[i];
            Rotamer roti[] = resi.getRotamers(library);
            int nri = roti.length;
            // Loop over the set of rotamers for residue i.
            boolean validA = false;
            for (int riA = 0; riA < nri; riA++) {
                if (!validRotamer(residues, i, riA)) {
                    if (riA == nri - 1 && !validA) {
                        logger.info(format(" No valid rotamers remain for %7s.", resi));
                        return false;
                    }
                    continue;
                } else {
                    validA = true;
                }
                // A 2nd loop over the set of rotamers for residue i.
                boolean validB = false;
                for (int riB = 0; riB < nri; riB++) {
                    if (!validRotamer(residues, i, riB)) {
                        if (riB == nri - 1 && !validB) {
                            logger.info(format(" No valid rotamers remain for %7s.", resi));
                            return false;
                        }
                        continue;
                    } else {
                        validB = true;
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
                        boolean validC = false;
                        for (int rjC = 0; rjC < nrj; rjC++) {
                            if (breakOut) {
                                break;
                            }
                            if (!validRotamer(residues, j, rjC)) {
                                if (rjC == nrj - 1 && !validC) {
                                    logger.info(format(" No valid rotamers remain for %7s.", resj));
                                    return false;
                                }
                                continue;
                            } else {
                                validC = true;
                            }
                            // A 2nd loop over the set of rotamers for residue j.
                            boolean validD = false;
                            for (int rjD = 0; rjD < nrj; rjD++) {
                                if (breakOut) {
                                    break;
                                }

                                if (!validRotamer(residues, j, rjD)) {
                                    if (rjD == nrj - 1 && !validD) {
                                        logger.info(format(" No valid rotamers remain for %7s.", resj));
                                        return false;
                                    }
                                    continue;
                                } else {
                                    validD = true;
                                }

                                // At least one rotamer of the eliminating rotamer pair must be different.
                                if (riA == riB && rjC == rjD) {
                                    continue;
                                }
                                if (!validRotamerPair(residues, i, riA, j, rjC)) {
                                    continue;
                                }
                                /**
                                 if (!validRotamerPair(residues, i, riB, j, rjD)) {
                                 continue;
                                 } */
                                // Try to eliminate R_i(riA) & R_j(rjC) using R_i(riB) & R_j(rjD)
                                if (goldsteinPairElimination(residues, i, riA, riB, j, rjC, rjD)) {
                                    eliminated = true;
                                    breakOut = true;
                                }
                            } // End inner loop over residue j's rotamers.
                        } // End outer loop over residue j's rotamers.
                        // Break out here if R_i(riA) & R_j(rjC) are eliminated (no need to eliminate twice).
                    } // End loop over residue j.
                } // End inner loop over residue i's rotamers.
            } // End outer loop over residue i's rotamers.
        } // End loop over residue i.
        return eliminated;
    }

    /**
     * Attempt to eliminate rotamer pair (ResI-RotA, ResJ-RotC) using (ResI-RotB, ResJ-RotD).
     *
     * @param residues
     * @param i        Index of the first residue.
     * @param riA      Index of the first residue's rotamer to eliminate.
     * @param riB      Index of the first residue's rotamer to use for elimination.
     * @param j        Index of the 2nd residue.
     * @param rjC      Index of the 2nd residue's rotamer to eliminate.
     * @param rjD      Index of the 2nd residue's rotamer to use for elimination.
     * @return Return true if eliminated.
     */
    private boolean goldsteinPairElimination(Residue residues[],
                                             int i, int riA, int riB, int j, int rjC, int rjD) {

        ArrayList<Residue> missedResidues = null;
        // Initialize the Goldstein energy.
        double goldsteinEnergy = self(i, riA) + self(j, rjC) + pair(i, riA, j, rjC)
                - self(i, riB) - self(j, rjD) - pair(i, riB, j, rjD);

        try {
            if (parallelTeam == null) {
                parallelTeam = new ParallelTeam();
            }
            if (goldsteinPairRegion == null) {
                goldsteinPairRegion = new GoldsteinPairRegion(parallelTeam.getThreadCount());
            }
            goldsteinPairRegion.init(residues, i, riA, riB, j, rjC, rjD);
            parallelTeam.execute(goldsteinPairRegion);
            goldsteinEnergy += goldsteinPairRegion.getSumOverK();
            missedResidues = goldsteinPairRegion.getMissedResidues();
        } catch (Exception e) {
            logger.log(Level.WARNING, " Exception in GoldsteinPairRegion.", e);
        }
        // goldsteinEnergy += goldsteinPairSumOverK(residues, 0, nres-1, i, riA, riB, j, rjC, rjD);

        if (goldsteinEnergy > ensembleBuffer) {
            if (missedResidues.isEmpty()) {
                if (eliminateRotamerPair(residues, i, riA, j, rjC, print)) {
                    logIfMaster(format("  Pair elimination of [(%7s,%2d),(%7s,%2d)] by [(%7s,%2d),(%7s,%2d)]: %12.4f > %6.4f",
                            residues[i], riA, residues[j], rjC, residues[i], riB, residues[j], rjD, goldsteinEnergy, ensembleBuffer));
                    return true;
                }
            } else {
                logIfMaster(format("  No Pair elimination of [(%7s,%2d),(%7s,%2d)] by [(%7s,%2d),(%7s,%2d)]: %12.4f > %6.4f",
                        residues[i], riA, residues[j], rjC, residues[i], riB, residues[j], rjD, goldsteinEnergy, ensembleBuffer));
                StringBuffer sb = new StringBuffer();
                for (int m = 0; m < missedResidues.size(); m++) {
                    Residue residueM = missedResidues.get(m);
                    sb.append(residueM);
                }
                logIfMaster(format("   due to %s.", sb));
            }
        }
        return false;
    }

    private class GoldsteinPairRegion extends ParallelRegion {

        Residue residues[];
        int i, riA, rjC;
        int j, riB, rjD;
        int nRes;
        GoldsteinRotamerPairLoop goldsteinRotamerPairLoop[];
        SharedDouble sharedSumOverK = new SharedDouble();
        ArrayList<Residue> blockedResidues;

        GoldsteinPairRegion(int nThreads) {
            goldsteinRotamerPairLoop = new GoldsteinRotamerPairLoop[nThreads];
        }

        public void init(Residue residues[], int i, int riA, int riB, int j, int rjC, int rjD) {
            this.residues = residues;
            this.i = i;
            this.riA = riA;
            this.riB = riB;
            this.j = j;
            this.rjC = rjC;
            this.rjD = rjD;
            nRes = residues.length;
        }

        public double getSumOverK() {
            return sharedSumOverK.get();
        }

        public ArrayList<Residue> getMissedResidues() {
            return blockedResidues;
        }

        public void start() {
            sharedSumOverK.set(0.0);
            blockedResidues = new ArrayList<>();
        }

        public void finish() {
            int nThreads = goldsteinRotamerPairLoop.length;
            for (int i = 0; i < nThreads; i++) {
                blockedResidues.addAll(goldsteinRotamerPairLoop[i].blockedResidues);
            }
        }

        @Override
        public void run() {
            int threadID = getThreadIndex();
            if (goldsteinRotamerPairLoop[threadID] == null) {
                goldsteinRotamerPairLoop[threadID] = new GoldsteinRotamerPairLoop();
            }
            try {
                execute(0, nRes - 1, goldsteinRotamerPairLoop[threadID]);
            } catch (Exception e) {
                logger.log(Level.WARNING, " Exception in GoldsteinPairRegion.", e);
            }
        }

        private class GoldsteinRotamerPairLoop extends IntegerForLoop {

            double sumOverK;
            ArrayList<Residue> blockedResidues;

            @Override
            public void start() {
                sumOverK = 0.0;
                blockedResidues = new ArrayList<>();
            }

            @Override
            public void finish() {
                sharedSumOverK.addAndGet(sumOverK);
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                sumOverK += goldsteinPairSumOverK(residues, lb, ub, i, riA, riB, j, rjC, rjD, blockedResidues);
            }
        }
    }

    private double goldsteinPairSumOverK(Residue residues[], int lb, int ub, int i, int riA, int riB,
                                         int j, int rjC, int rjD,
                                         ArrayList<Residue> blockedResidues) {
        double sumOverK = 0.0;
        int nres = residues.length;

        for (int k = lb; k <= ub; k++) {
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
                if (check(k, rk)) {
                    continue;
                }
                // These are not checked in Osprey.
                /**
                 if (check(i, riA, k, rk) || check(j, rjC, k, rk) || check(i, riB, k, rk) || check(j, rjD, k, rk)) {
                 continue;
                 } */
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
                minForResK = 0.0;
                blockedResidues.add(resk);
            }
            sumOverK += minForResK;
        }
        return sumOverK;
    }

    private boolean eliminateRotamer(Residue[] residues, int i, int ri, boolean verbose) {
        // Check if rotamer (i, ri) has already been eliminated.
        if (check(i, ri)) {
            return false;
        }

        // Make sure more than 1 valid rotamers are left.
        int rotCount = rotamerCount(residues, i);
        if (rotCount < 2) {
            return false;
        }

        eliminatedSingles[i][ri] = true;
        rotCount--;

        if (verbose) {
            logIfMaster(format(" Rotamer (%7s,%2d) eliminated (%2d left).", residues[i], ri, rotCount));
        }
        int eliminatedPairs = eliminateRotamerPairs(residues, i, ri, verbose);
        if (eliminatedPairs > 0 && verbose) {
            logIfMaster(format("  Eliminated %2d rotamer pairs.", eliminatedPairs));
        }
        return true;
    }

    private boolean eliminateRotamer(Residue[] residues, int i, int ri, int rj, boolean verbose) {
        // Check the rotamer (i, ri) has not already been eliminated.
        if (check(i, ri) || ri == rj) {
            return false;
        }

        // Check that rotamer rj is valid to eliminate ri.
        if (!validRotamer(residues, i, rj)) {
            return false;
        }

        eliminatedSingles[i][ri] = true;
        int rotCount = rotamerCount(residues, i);
        if (verbose) {
            logIfMaster(format(" %7s rotamer %2d eliminated (%2d left).", residues[i], ri, rotCount));
        }
        int eliminatedPairs = eliminateRotamerPairs(residues, i, ri, verbose);
        if (eliminatedPairs > 0 && verbose) {
            logIfMaster(format("  Eliminated %2d rotamer pairs.", eliminatedPairs));
        }
        return true;
    }

    private int eliminateRotamerPairs(Residue[] residues, int i, int ri, boolean verbose) {
        int nres = residues.length;
        int eliminatedPairs = 0;
        for (int j = 0; j < nres; j++) {
            if (j == i) {
                continue;
            }
            Residue residuej = residues[j];
            Rotamer rotamersj[] = residuej.getRotamers(library);
            int lenrj = rotamersj.length;
            for (int rj = 0; rj < lenrj; rj++) {
                if (eliminateRotamerPair(residues, i, ri, j, rj, verbose)) {
                    eliminatedPairs++;
                }
            }
        }
        return eliminatedPairs;
    }

    private boolean eliminateRotamerPair(Residue[] residues, int i, int ri, int j, int rj, boolean verbose) {
        if (i > j) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        if (!check(i, ri, j, rj)) {
            eliminatedPairs[i][ri][j][rj] = true;
            if (verbose) {
                logIfMaster(format("  Rotamer pair eliminated: [(%7s,%2d) (%7s,%2d)]", residues[i], ri, residues[j], rj));
            }
            return true;
        } else {
            return false;
        }
    }

    /**
     * Count the rotamers remaining for residue i.
     *
     * @param residues Residue array.
     * @param i        The residue number to examine.
     * @return The remaining rotamer count.
     */
    private int rotamerCount(Residue residues[], int i) {
        int nRes = residues.length;
        Rotamer rotI[] = residues[i].getRotamers(library);
        int ni = rotI.length;
        // Initialize the count to 0.
        int rotCount = 0;
        // Loop over the rotamers for residue i.
        for (int ri = 0; ri < ni; ri++) {
            if (check(i, ri)) {
                continue;
            }
            // Check that rotamer ri has valid pairs with all other residues.
            boolean valid = true;
            for (int j = 0; j < nRes; j++) {
                if (i == j) {
                    continue;
                }
                if (rotamerPairCount(residues, i, ri, j) == 0) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                rotCount++;
            }
        }
        return rotCount;
    }

    /**
     * Validate residue i with rotamer ri.
     *
     * @param residues The residues being optimized.
     * @param i        The residue to validate.
     * @param ri       The rotamer to validate.
     * @return The status of this rotamer.
     */
    private boolean validRotamer(Residue residues[], int i, int ri) {
        // Return false if this rotamer has been eliminated.
        if (check(i, ri)) {
            return false;
        }

        // Loop over all residues to check for valid pairs and triples.
        int n = residues.length;
        for (int j = 0; j < n; j++) {
            if (j == i) {
                continue;
            }
            if (rotamerPairCount(residues, i, ri, j) == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Validate rotamer pair (i, ri) and (j, rj).
     *
     * @param residues The residues being optimized.
     * @param i        The first residue to validate.
     * @param ri       The first rotamer to validate.
     * @param j        The 2nd residue to validate.
     * @param rj       The 2nd rotamer to validate.
     * @return The status of this rotamer.
     */
    private boolean validRotamerPair(Residue residues[], int i, int ri, int j, int rj) {
        // Residues i and j must be different.
        if (i == j) {
            return false;
        }

        // Return false if either rotamer is not valid.
        if (!validRotamer(residues, i, ri) || !validRotamer(residues, j, rj)) {
            return false;
        }

        // Return false if the rotamer pair has been eliminated.
        if (check(i, ri, j, rj)) {
            return false;
        }

        // Loop over all residues to check for valid triples.
        int n = residues.length;
        for (int k = 0; k < n; k++) {
            if (k == i || k == j) {
                continue;
            }
            // There must be at least one valid rotamer triple between (i, ri) and (j, rj) with residue k.
            if (rotamerTripleCount(residues, i, ri, j, rj, k) == 0) {
                // Eliminate the pair?
                return false;
            }
        }
        return true;
    }

    /**
     * Count the rotamer pairs remaining for residues i and j.
     *
     * @param residues Residue array.
     * @param i        The first residue to examine.
     * @param j        The second residue to examine.
     * @return The remaining rotamer pair count.
     */
    private int rotamerPairCount(Residue residues[], int i, int j) {
        int pairCount = 0;
        Rotamer rotI[] = residues[i].getRotamers(library);
        int ni = rotI.length;
        Rotamer rotJ[] = residues[j].getRotamers(library);
        int nj = rotJ.length;
        for (int ri = 0; ri < ni; ri++) {
            if (!check(i, ri)) {
                for (int rj = 0; rj < nj; rj++) {
                    if (!check(j, rj) && !check(i, ri, j, rj)) {
                        pairCount++;
                    }
                }
            }
        }
        return pairCount;
    }

    /**
     * Count the rotamer pairs remaining for (residue i, rotamer ri) and residue j.
     *
     * @param residues Residue array.
     * @param i        The first residue to examine.
     * @param ri       The rotamer for the first residue.
     * @param j        The second residue to examine.
     * @return The remaining rotamer pair count.
     */
    private int rotamerPairCount(Residue residues[], int i, int ri, int j) {
        if (i == j || check(i, ri)) {
            return 0;
        }
        int pairCount = 0;
        Rotamer rotJ[] = residues[j].getRotamers(library);
        int nj = rotJ.length;
        // Loop over all rotamers for residue j.
        for (int rj = 0; rj < nj; rj++) {
            // Check for a valid rotamer pair.
            if (!check(j, rj) && !check(i, ri, j, rj)) {
                // Loop over all residues k to check for valid rotamer triples.
                int nRes = residues.length;
                boolean valid = true;
                for (int k = 0; k < nRes; k++) {
                    if (k == i || k == j) {
                        continue;
                    }
                    if (rotamerTripleCount(residues, i, ri, j, rj, k) == 0) {
                        valid = false;
                    }
                }
                if (valid) {
                    pairCount++;
                }
            }
        }
        return pairCount;
    }

    /**
     * Count the rotamer triples remaining for (residue i, rotamer ri) and (residue j, rotamer rj) with residue k.
     *
     * @param residues Residue array.
     * @param i        The first residue to examine.
     * @param ri       The rotamer for the first residue.
     * @param j        The second residue to examine.
     * @param rj       The rotamer for the first residue.
     * @param k        The third residue.
     * @return The remaining rotamer triples count.
     */
    private int rotamerTripleCount(Residue residues[], int i, int ri, int j, int rj, int k) {
        if (i == j || i == k || j == k) {
            return 0;
        }
        int tripleCount = 0;
        Rotamer rotK[] = residues[k].getRotamers(library);
        int nk = rotK.length;
        // Check that each rotamer and their pair have not be eliminated.
        if (!check(i, ri) && !check(j, rj) && !check(i, ri, j, rj)) {
            // Loop over all rotamers for residue k.
            for (int rk = 0; rk < nk; rk++) {
                // Check for a valid rotamer triple.
                if (!check(k, rk) && !check(i, ri, k, rk) && !check(j, rj, k, rk)) {
                    tripleCount++;
                }
            }
        }
        return tripleCount;
    }

    private double self(int i, int ri) {
        try {
            return selfEnergy[i][ri];
        } catch (NullPointerException npe) {
            logger.info(format(" NPE for self energy (%3d,%2d).", i, ri));
            throw npe;
        }
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
        try {
            return twoBodyEnergy[i][ri][j][rj];
        } catch (NullPointerException npe) {
            logger.info(format(" NPE for 2-body energy (%3d,%2d) (%3d,%2d).", i, ri, j, rj));
            throw npe;
        }
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
        try {
            return threeBodyEnergy[i][ri][j][rj][k][rk];
        } catch (NullPointerException npe) {
            logger.info(format(" NPE for 3-body energy (%3d,%2d) (%3d,%2d) (%3d,%2d).", i, ri, j, rj, k, rk));
            throw npe;
        }
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

    public void setSingletonClashThreshold(double singletonClashThreshold) {
        this.clashThreshold = singletonClashThreshold;
    }

    public void setPairClashThreshold(double pairClashThreshold) {
        this.pairClashThreshold = pairClashThreshold;
    }

    public void setWindowSize(int windowSize) {
        this.windowSize = windowSize;
        if (this.increment > windowSize) {
            logger.info(format(" Decreasing increment to match window size %d", windowSize));
            this.increment = windowSize;
        }
    }

    public void setVideoWriter(boolean writeVideo) {
        this.writeVideo = writeVideo;
        if (writeVideo) {
            this.videoWriter = new VideoWriter(molecularAssembly, true, true);
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
     * @param i  Residue i.
     * @param ri Rotamer ri.
     * @return True if rotamer eliminated.
     */
    protected boolean check(int i, int ri) {
        if (eliminatedSingles == null) {
            return false;
        }
        return eliminatedSingles[i][ri];
    }

    /**
     * Check for eliminated rotamer pair; true if eliminated.
     *
     * @param i  Residue i.
     * @param ri Rotamer ri.
     * @param j  Residue j.
     * @param rj Rotamer rj.
     * @return True if eliminated pair.
     */
    protected boolean check(int i, int ri, int j, int rj) {
        if (eliminatedPairs == null) {
            return false;
        }
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
     * Checks to see if any eliminations with j,rj have occurred; assumes i,ri
     * self has already been checked. Checks j,rj self and i,ri,j,rj pair. The
     * intent is to be part of a loop over i,ri,j,rj, and check for eliminations
     * at the j,rj point.
     *
     * @param i  Residue i
     * @param ri Rotamer ri
     * @param j  Residue j
     * @param rj Rotamer rj
     * @return j eliminated with i
     */
    protected boolean checkToJ(int i, int ri, int j, int rj) {
        return (check(j, rj) || check(i, ri, j, rj));
    }

    /**
     * Checks to see if any eliminations with k,rk have occurred; assumes
     * i,ri,j,rj pair has already been checked. Checks the k,rk self, all pairs
     * with k,rk, and the i,ri,j,rj,k,rk triple. The intent is to be part of a
     * loop over i,ri,j,rj,k,rk, and check for eliminations at the k,rk point.
     *
     * @param i  Residue i
     * @param ri Rotamer ri
     * @param j  Residue j
     * @param rj Rotamer rj
     * @param k  Residue k
     * @param rk Rotamer rk
     * @return k eliminated with i,j
     */
    protected boolean checkToK(int i, int ri, int j, int rj, int k, int rk) {
        return (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk));
    }

    /**
     * Checks to see if any eliminations with l,rl have occurred; assumes
     * i,ri,j,rj,k,rk triple has already been checked. Checks the l,rl self, all
     * pairs with l,rl, all triples with l,rl, and the quad. The intent is to be
     * part of a loop over i,ri,j,rj,k,rk,l,rl, and check for eliminations at
     * the l,rl point.
     *
     * @param i  Residue i
     * @param ri Rotamer ri
     * @param j  Residue j
     * @param rj Rotamer rj
     * @param k  Residue k
     * @param rk Rotamer rk
     * @param l  Residue l
     * @param rl Rotamer rl
     * @return l eliminated with i,j,k
     */
    protected boolean checkToL(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        return (check(l, rl) || check(i, ri, l, rl) || check(j, rj, l, rl) || check(k, rk, l, rl));
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
                logger.severe(format(" Coding error: all %d rotamers for residue %s eliminated.", ni, residuei));
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
                    logger.severe(format(" Coding error: all pairs for %s with residue %s eliminated.",
                            residuei, residuej));
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
        StringBuilder sb = new StringBuilder(format(" %d out of %d rotamers eliminated.\n", singles, rotamerCount));
        sb.append(format(" %d out of %d rotamer pairs eliminated.", pairs, pairCount));
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
         *                              indexed by 1-3 min x,y,z, 4-6 max x,y,z
         * @param indices               Index of cell along a, b, and c (x, y, and z).
         * @param linearIndex           Index of box in linear box array.
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
         *                              indexed by 1-3 min x,y,z, 4-6 max x,y,z
         * @param indices               Index of cell along a, b, and c (x, y, and z).
         * @param linearIndex           Index of box in linear box array.
         * @param residuesIn            Array of Residues to initialize the BoxOptCell
         *                              with.
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
         *                                   residue.
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

    private int nameToNumber(String residueString, Residue residues[]) throws NumberFormatException {
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

    private int loadEnergyRestart(File restartFile, Residue residues[]) {
        return loadEnergyRestart(restartFile, residues, -1, null);
    }

    private int loadEnergyRestart(File restartFile, Residue residues[], int boxIteration, int[] cellIndices) {
        try {
            int nResidues = residues.length;
            Path path = Paths.get(restartFile.getCanonicalPath());
            List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
            List<String> linesThisBox = new ArrayList<>();

            try {
                backboneEnergy = computeBackboneEnergy(residues);
            } catch (ArithmeticException ex) {
                logger.severe(String.format(" Exception %s in calculating backbone energy; FFX shutting down.", ex.toString()));
            }

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
                    logIfMaster(format(" Didn't find restart energies for Box %d: %d,%d,%d",
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
                logger.warning(format(" Empty or unreadable energy restart file: %s.", restartFile.getCanonicalPath()));
            }
            if (loaded >= 1) {
                singlesMap.clear();
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
                        singlesMap.put(singleJobIndex, selfJob);
                        String revKey = format("%d %d", i, ri);
                        reverseJobMapSingles.put(revKey, singleJobIndex);
                        singleJobIndex++;
                    }
                }
                // fill in self-energies from file while removing the corresponding jobs from singlesMap
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
                            if (verbose) {
                                logIfMaster(format(" From restart file: Self energy %3d (%7s,%2d): %12.4f", i, residues[i], ri, energy));
                            }
                        } catch (Exception e) {
                            logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d", i, ri);
                        Integer ret[] = singlesMap.remove(reverseJobMapSingles.get(revKey));
                        if (ret == null) {
                            //logIfMaster(format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format(" Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded self energies from restart file.");
                // prune singles
                if (pruneClashes) {
                    pruneSingleClashes(residues);
                }
            }
            if (loaded >= 2) {
                if (singlesMap.size() > 0) {
                    logIfMaster(" Double-check that parameters match original run due to missing self-energies.");
                }
                pairsMap.clear();
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
                                pairsMap.put(pairJobIndex, pairJob);
                                String revKey = format("%d %d %d %d", i, ri, j, rj);
                                reverseJobMapPairs.put(revKey, pairJobIndex);
                                pairJobIndex++;
                            }
                        }
                    }
                }
                // fill in pair-energies from file while removing the corresponding jobs from pairsMap
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
                            if (verbose) {
                                logIfMaster(format(" From restart file: Pair energy [(%7s,%2d),(%7s,%2d)]: %12.4f",
                                        residues[i], ri, residues[j], rj, energy));
                            }
                        } catch (Exception e) {
                            logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d %d %d", i, ri, j, rj);
                        Integer ret[] = pairsMap.remove(reverseJobMapPairs.get(revKey));
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded pair energies from restart file.");
                // prune pairs
                if (prunePairClashes) {
                    prunePairClashes(residues);
                }
            }
            if (loaded >= 3) {
                if (pairsMap.size() > 0) {
                    if (master) {
                        logger.warning("Double-check that parameters match original run!  Found trimers in restart file, but pairs job queue is non-empty.");
                    }
                }
                HashMap<String, Integer> reverseJobMapTrimers = new HashMap<>();
                trimersMap.clear();
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
                                        trimersMap.put(trimerJobIndex, trimerJob);
                                        String revKey = format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                                        reverseJobMapTrimers.put(revKey, trimerJobIndex);
                                        trimerJobIndex++;
                                    }
                                }
                            }
                        }
                    }
                }
                // fill in triple-energies from file while removing the corresponding jobs from trimersMap
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
                            logger.log(Level.SEVERE, format("Restart file contained an out-of-bounds index.  Offending line: %s", line), ex);
                        }
                        if (verbose) {
                            logIfMaster(format(" From restart file: Trimer energy %3d %-2d, %3d %-2d, %3d %-2d: %16.8f", i, ri, j, rj, k, rk, energy));
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                        Integer ret[] = trimersMap.remove(reverseJobMapTrimers.get(revKey));
                        if (ret == null) {
                            //logIfMaster(format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format("Unparsable line in energy restart file: \n%s", line), ex);
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
                // this happens when wrapping globalOptimization with e.g. box optimization
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
            logger.info(format(" Energy restart file: %s", restartFile.getName()));
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
                boxHeader = format("Box %d: %d,%d,%d", iteration, cellIndices[0], cellIndices[1], cellIndices[2]);
            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Couldn't open energy restart file.", ex);
            }
            logger.info(format(" Energy restart file: %s", restartFile.getName()));
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
                                energiesToWrite.add(format("Self %d %d: %16.8f", resi, roti, energy));
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
                                energiesToWrite.add(format("Pair %d %d, %d %d: %16.8f", resi, roti, resj, rotj, energy));
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
                                    energiesToWrite.add(format("Triple %d %d, %d %d, %d %d: %16.8f", resi, roti, resj, rotj, resk, rotk, energy));
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
    private void pruneSingleClashes(Residue residues[]) {
        if (!pruneClashes) {
            return;
        }
        for (int i = 0; i < residues.length; i++) {
            Residue residue = residues[i];
            Rotamer rotamers[] = residue.getRotamers(library);
            int nrot = rotamers.length;
            double minEnergy = Double.MAX_VALUE;
            int minRot = -1;
            for (int ri = 0; ri < nrot; ri++) {
                if (!check(i, ri) && self(i, ri) < minEnergy) {
                    minEnergy = self(i, ri);
                    minRot = ri;
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
                    if (eliminateRotamer(residues, i, ri, print)) {
                        logIfMaster(format("  Rotamer (%7s,%2d) self-energy %12.4f pruned by (%7s,%2d) %12.4f.",
                                residue, ri, self(i, ri), residue, minRot, minEnergy));
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
    private void prunePairClashes(Residue residues[]) {
        if (!prunePairClashes) {
            return;
        }
        int nResidues = residues.length;
        // Loop over first residue.
        for (int i = 0; i < nResidues - 1; i++) {
            Residue resi = residues[i];
            Rotamer[] roti = resi.getRotamers(library);
            int ni = roti.length;
            // Loop over second residue.
            for (int j = i + 1; j < nResidues; j++) {
                Residue resj = residues[j];
                Rotamer[] rotj = resj.getRotamers(library);
                int nj = rotj.length;
                double minPair = Double.MAX_VALUE;
                // Loop over the rotamers for residue i.
                for (int ri = 0; ri < ni; ri++) {
                    if (!validRotamer(residues, i, ri)) {
                        continue;
                    }
                    // Loop over rotamers for residue j.
                    for (int rj = 0; rj < nj; rj++) {
                        if (!validRotamer(residues, j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        double pairEnergy = pair(i, ri, j, rj) + self(i, ri) + self(j, rj);
                        if (minPair > pairEnergy) {
                            minPair = pairEnergy;
                        }
                    }
                }

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
                        break;
                    case 2:
                        threshold *= pruningFactor;
                        break;
                    default:
                        throw new ArithmeticException(" RotamerOptimization.prunePairClashes() has somehow "
                                + "found less than zero or more than two nucleic acid residues in a pair of"
                                + " residues. This result should be impossible.");
                }
                threshold += minPair;

                // Check for elimination of any rotamer ri.
                for (int ri = 0; ri < ni; ri++) {
                    if (!validRotamer(residues, i, ri)) {
                        continue;
                    }

                    double eliminate = Double.MAX_VALUE;
                    for (int rj = 0; rj < nj; rj++) {
                        if (!validRotamer(residues, j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        double pairEnergy = pair(i, ri, j, rj) + self(i, ri) + self(j, rj);
                        if (pairEnergy < eliminate) {
                            eliminate = pairEnergy;
                        }
                    }
                    if (eliminate > threshold) {
                        if (eliminateRotamer(residues, i, ri, print)) {
                            logIfMaster(format(
                                    "  Rotamer (%7s,%2d) pruned by clashes with all %7s rotamers %12.4f >> %12.4f.",
                                    resi, ri, resj, eliminate, minPair));
                        }
                    }
                }

                // Check for elimination of any rotamer rj.
                for (int rj = 0; rj < nj; rj++) {
                    if (!validRotamer(residues, j, rj)) {
                        continue;
                    }

                    double eliminate = Double.MAX_VALUE;
                    for (int ri = 0; ri < ni; ri++) {
                        if (!validRotamer(residues, i, ri) || check(i, ri, j, rj)) {
                            continue;
                        }
                        double pairEnergy = pair(i, ri, j, rj) + self(i, ri) + self(j, rj);
                        if (pairEnergy < eliminate) {
                            eliminate = pairEnergy;
                        }
                    }
                    if (eliminate > threshold) {
                        eliminateRotamer(residues, j, rj, print);
                        logIfMaster(format(
                                "  Rotamer (%7s,%2d) pruned by clashes with all %7s rotamers %12.4f >> %12.4f.",
                                resj, rj, resi, eliminate, minPair));
                    }
                }
            }
        }
    }

    private void revertResidue(Residue res, ResidueState resState, Rotamer rotamer) {
        if (res.getResidueType() == NA) {
            res.revertState(resState);
        } else {
            RotamerLibrary.applyRotamer(res, rotamer);
        }
        turnOffAtoms(res);
    }

    private class SinglesEnergyRegion extends WorkerRegion {

        private final SinglesEnergyLoop energyLoop;
        private final Residue residues[];
        private final int nResidues;
        private final boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        private int singlesEnergyCounter;

        public SinglesEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new SinglesEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
        }

        @Override
        public void start() {
            singlesEnergyCounter = 0;
            /**
             * Setup and compute backbone energy.
             */
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer rotamers[] = residue.getRotamers(library);
                switch (residue.getResidueType()) {
                    case NA:
                        break;
                    case AA:
                    default:
                        RotamerLibrary.applyRotamer(residue, rotamers[0]);
                        break;
                }
                turnOffAtoms(residue);
            }
            try {
                backboneEnergy = currentEnergy(new ArrayList<>());
            } catch (ArithmeticException ex) {
                logger.severe(format(" Error in calculation of backbone energy %s", ex.getMessage()));
            }
            logIfMaster(format(" Backbone energy:  %16.8f\n", backboneEnergy));
        }

        @Override
        public void run() throws Exception {
            if (!singlesMap.isEmpty()) {
                execute(singlesMap.keySet().iterator(), energyLoop);
            }
        }

        @Override
        public void finish() {
            /**
             * Broadcast the "I'm finished" signal.
             */
            double finished[] = new double[3];
            for (int x = 0; x < finished.length; x++) {
                finished[x] = -1.0;
            }
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            /**
             * Wait for everyone else to send in their self energies.
             */
            int waiting = 0;
            while (!selfsDone) {
                try {
                    Thread.sleep(POLLING_FREQUENCY);
                    if (waiting++ == 1000) {
                        logger.warning(format(" Process %d experiencing long wait for others' trimer energies.", Comm.world().rank()));
                    }
                } catch (InterruptedException e) {
                    logger.warning(" Parallel self-energy computation was interrupted.");
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
                        logger.info(format(" Self energy %7s %-2d: %16.8f", residues[i], ri, self(i, ri)));
                    }
                }
            }

            logger.info(format(" Self-energy count: %d.", singlesEnergyCounter));
        }

        private class SinglesEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                /**
                 * Compute the self-energy for one rotamer.
                 */
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!singlesMap.keySet().contains(jobKey)) {
                        continue;
                    }
                    Integer job[] = singlesMap.get(jobKey);
                    int i = job[0];
                    int ri = job[1];
                    if (pruneClashes && check(i, ri)) {
                        continue;
                    }
                    long time = -System.nanoTime();
                    Residue resi = residues[i];
                    Rotamer roti = resi.getRotamers(library)[ri];
                    ResidueState resiOriginal = (resi.getResidueType() == NA ? resi.storeState() : null);
                    turnOnAtoms(resi);
                    RotamerLibrary.applyRotamer(resi, roti);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                    if (writeVideo) {
                        videoWriter.write(format("%03d-%03d", i, ri));
                    }
                    double selfEnergy;
                    if (writeVideo || skipEnergies) {
                        selfEnergy = 0;
                        logger.info(format(" Self %7s %-2d: %16.8f.", resi, ri, selfEnergy));
                    } else {
                        List<Residue> rList = Arrays.asList(resi);
                        try {
                            selfEnergy = currentEnergy(rList) - backboneEnergy;
                            time += System.nanoTime();
                            logger.info(format(" Self %7s %-2d: %16.8f in %6.4f (sec).", resi, ri, selfEnergy, time * 1.0e-9));
                        } catch (ArithmeticException ex) {
                            selfEnergy = 1E100;
                            time += System.nanoTime();
                            logger.info(format(" Self %7s %-2d: set to %10.6g (unreasonable conformation) in %6.4f (sec).", resi, ri, selfEnergy, time * 1.0e-9));
                        }
                        singlesEnergyCounter++;
                    }

                    // Revert residues and turn off atoms.
                    revertResidue(resi, resiOriginal, resi.getRotamers(library)[0]);

                    double anEnergy[] = new double[3];
                    anEnergy[0] = i;
                    anEnergy[1] = ri;
                    anEnergy[2] = selfEnergy;
                    DoubleBuf anEnergyBuf = DoubleBuf.buffer(anEnergy);
                    multicastBuf(anEnergyBuf);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                }
            }
        }
    }

    private class PairsEnergyRegion extends WorkerRegion {

        private final PairsEnergyLoop energyLoop;
        private final Residue residues[];
        private final int nResidues;
        private final boolean useOrigCoordsRot = library.getUsingOrigCoordsRotamer();
        private int pairsEnergyCounter;

        public PairsEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new PairsEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
        }

        @Override
        public void start() {
            pairsEnergyCounter = 0;
        }

        @Override
        public void run() throws Exception {
            if (!pairsMap.isEmpty()) {
                execute(pairsMap.keySet().iterator(), energyLoop);
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
                        logger.warning(format("Process %d experiencing long wait for others' pair energies.", Comm.world().rank()));
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
                                logger.info(format(" (recap) Pair energy %7s %-2d, %7s %-2d: %16.8f", residues[i], ri, residues[j], rj, pair(i, ri, j, rj)));
                            }
                        }
                    }
                }
            }
            logger.info(format(" Pairs Energy Count: %d", pairsEnergyCounter));
        }

        private class PairsEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                /**
                 * Compute the pair-energy for each pair of rotamers using a pair-level job indexing method.
                 */
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!pairsMap.keySet().contains(jobKey)) {
                        continue;
                    }
                    Integer job[] = pairsMap.get(jobKey);
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
                    ResidueState resiOriginal = resi.getResidueType() == NA ? resi.storeState() : null;
                    ResidueState resjOriginal = resj.getResidueType() == NA ? resj.storeState() : null;

                    // Turn on atoms and apply rotamers.
                    long time = -System.nanoTime();
                    turnOnAtoms(resi);
                    turnOnAtoms(resj);
                    RotamerLibrary.applyRotamer(resi, roti);
                    RotamerLibrary.applyRotamer(resj, rotj);
                    if (algorithmListener != null) {
                        algorithmListener.algorithmUpdate(molecularAssembly);
                    }
                    if (writeVideo) {
                        videoWriter.write(format("%03d-%03d_%03d-%03d", i, ri, j, rj));
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
                            logger.info(format(" Pair %7s %-2d, %7s %-2d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang",
                                    resi, ri, resj, rj, dist, superpositionThreshold));
                        } else {
                            String distString = (dist < Double.MAX_VALUE) ? format("%10.3f", dist) : format("     large");
                            try {
                                twoBodyEnergy = currentEnergy(rList) - self(i, ri) - self(j, rj) - backboneEnergy;
                                time += System.nanoTime();
                                logger.info(format(" Pair %7s %-2d, %7s %-2d: %16.8f at %s (Ang) in %6.4f (sec).",
                                        resi, ri, resj, rj, twoBodyEnergy, distString, time * 1.0e-9));
                            } catch (ArithmeticException ex) {
                                twoBodyEnergy = 1E100;
                                time += System.nanoTime();
                                logger.info(format(" Pair %7s %-2d, %7s %-2d: set to %10.6g (unreasonable conformation) at %s (Ang) in %6.4f (sec).",
                                        resi, ri, resj, rj, twoBodyEnergy, distString, time * 1.0e-9));
                            }
                            pairsEnergyCounter++;
                        }
                    } else {
                        try {
                            twoBodyEnergy = currentEnergy(rList) - self(i, ri) - self(j, rj) - backboneEnergy;
                            time += System.nanoTime();
                            logger.info(format(" Pair %7s %-2d, %7s %-2d: %16.8f in %6.4f (sec).",
                                    resi, ri, resj, rj, twoBodyEnergy, time * 1.0e-9));
                        } catch (ArithmeticException ex) {
                            twoBodyEnergy = 1E100;
                            time += System.nanoTime();
                            logger.info(format(" Pair %7s %-2d, %7s %-2d: set to %10.6g (unreasonable conformation) in %6.4f (sec).",
                                    resi, ri, resj, rj, twoBodyEnergy, time * 1.0e-9));
                        }
                        pairsEnergyCounter++;
                    }

                    // Get energy and broadcast it.
                    double commEnergy[] = new double[5];
                    commEnergy[0] = i;
                    commEnergy[1] = ri;
                    commEnergy[2] = j;
                    commEnergy[3] = rj;
                    commEnergy[4] = twoBodyEnergy;
                    DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
                    multicastBuf(commEnergyBuf);

                    // Revert residues and turn off atoms.
                    revertResidue(resi, resiOriginal, resi.getRotamers(library)[0]);
                    revertResidue(resj, resjOriginal, resj.getRotamers(library)[0]);

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
        private double localDistanceMatrix[][][][];
        private int triplesEnergyCounter;

        public TriplesEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new TriplesEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
            localDistanceMatrix = new double[nResidues][][][];
        }

        @Override
        public void start() {
            triplesEnergyCounter = 0;
            if (distance <= 0) {
                logger.info(" Calculating local distance matrix using non-eliminated rotamers.");
                // TODO: check on the location of this call - might need to be done per trimer-job
                localSequentialDistanceMatrix(residues, localDistanceMatrix);
            }
        }

        @Override
        public void run() throws Exception {
            if (!trimersMap.isEmpty()) {
                execute(trimersMap.keySet().iterator(), energyLoop);
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
                        logger.warning(format("Process %d experiencing long wait for others' trimer energies.", Comm.world().rank()));
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
                                        logger.info(format(" (recap) Trimer energy %7s %-2d, %7s %-2d, %7s %-2d: %16.8f", resi, ri, resj, rj, resk, rk, triple(i, ri, j, rj, k, rk)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            logger.info(format(" Triples Energy Count: %d", triplesEnergyCounter));
        }

        private class TriplesEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                // Trimer-level job indexing method
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!trimersMap.keySet().contains(jobKey)) {
                        continue;
                    }
                    Integer job[] = trimersMap.get(jobKey);
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
                            logger.info(format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang.",
                                    resi, ri, resj, rj, resk, rk, dist, superpositionThreshold));
                        } else {
                            long time = -System.nanoTime();
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
                                videoWriter.write(format("%03d-%03d_%03d-%03d_%03d-%03d", i, ri, j, rj, k, rk));
                            }
                            // Compute trimer and broadcast it
                            if (writeVideo || skipEnergies) {
                                threeBodyEnergy = 0;
                            } else {
                                String distString = (dist < Double.MAX_VALUE) ? format("%10.3f", dist) : format("     large");
                                try {
                                    threeBodyEnergy = currentEnergy(rList)
                                            - self(i, ri) - self(j, rj) - self(k, rk)
                                            - pair(i, ri, j, rj) - pair(i, ri, k, rk)
                                            - pair(j, rj, k, rk) - backboneEnergy;
                                    time += System.nanoTime();
                                    logger.info(format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %s (Ang) in %6.4f (sec).",
                                            resi, ri, resj, rj, resk, rk, threeBodyEnergy, distString, time * 1.0e-9));
                                } catch (ArithmeticException ex) {
                                    threeBodyEnergy = 1E100;
                                    time += System.nanoTime();
                                    logger.info(format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d: set to %10.6g (unreasonable conformation) at %s (Ang) in %6.4f (sec).",
                                            resi, ri, resj, rj, resk, rk, threeBodyEnergy, distString, time * 1.0e-9));

                                }
                                triplesEnergyCounter++;
                            }

                            // Revert rotamers and turn off atoms.
                            revertResidue(resi, resiOriginal, resi.getRotamers(library)[0]);
                            revertResidue(resj, resjOriginal, resj.getRotamers(library)[0]);
                            revertResidue(resk, reskOriginal, resk.getRotamers(library)[0]);
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
                        if (logger.isLoggable(Level.FINEST)) {
                            logger.finest(format(" Trimer %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %10.3f Ang.",
                                    resi, ri, resj, rj, resk, rk, 0.0, threeBodyCutoffDist));
                        }
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
        private double localDistanceMatrix[][][][];
        private int quadsEnergyCounter;

        public QuadsEnergyRegion(int nt, Residue residues[]) {
            energyLoop = new QuadsEnergyLoop();
            this.residues = residues;
            nResidues = residues.length;
            localDistanceMatrix = new double[nResidues][][][];
        }

        @Override
        public void start() {
            quadsEnergyCounter = 0;
            if (distance <= 0) {
                logger.info(" Calculating local distance matrix using non-eliminated rotamers.");
                // TODO: check on the location of this call - might need to be done per quad-job
                localSequentialDistanceMatrix(residues, localDistanceMatrix);
            }
        }

        @Override
        public void run() throws Exception {
            if (!quadsMap.isEmpty()) {
                execute(quadsMap.keySet().iterator(), energyLoop);
            }
        }

        @Override
        public void finish() {
            // no "I'm finished" signal for quads
            logger.info(format(" Quads Energy Count: %d", quadsEnergyCounter));
        }

        private class QuadsEnergyLoop extends WorkerIteration<Integer> {

            @Override
            public void run(Integer key) {
                // Quad-level job indexing method
                for (int jobKey = key; jobKey <= key; jobKey++) {
                    if (!quadsMap.keySet().contains(jobKey)) {
                        continue;
                    }
                    Integer job[] = quadsMap.get(jobKey);
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
                            logger.info(format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d:   set to 1.0E100 at %13.6f Ang < %5.3f Ang.",
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
                                String distString = (dist < Double.MAX_VALUE) ? format("%10.3f", dist) : format("     large");
                                try {
                                    quadEnergy = currentEnergy(rList)
                                            - self(i, ri) - self(j, rj) - self(k, rk) - self(l, rl)
                                            - pair(i, ri, j, rj) - pair(i, ri, k, rk) - pair(i, ri, l, rl)
                                            - pair(j, rj, k, rk) - pair(j, rj, l, rl) - pair(k, rk, l, rl)
                                            - triple(i, ri, j, rj, k, rk) - triple(i, ri, j, rj, l, rl) - triple(i, ri, k, rk, l, rl) - triple(j, rj, k, rk, l, rl)
                                            - backboneEnergy;
                                    if (Math.abs(quadEnergy) > 1.0) {
                                        StringBuilder sb = new StringBuilder();
                                        sb.append(format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %s Ang.\n",
                                                resi, ri, resj, rj, resk, rk, resl, rl, quadEnergy, distString));
                                        sb.append(format("   Explain: (ref %d) \n", jobKey));
                                        sb.append(format("     Self %3d %3d:                  %.3f\n", i, ri, self(i, ri)));
                                        sb.append(format("     Self %3d %3d:                  %.3f\n", j, rj, self(j, rj)));
                                        sb.append(format("     Self %3d %3d:                  %.3f\n", k, rk, self(k, rk)));
                                        sb.append(format("     Self %3d %3d:                  %.3f\n", l, rl, self(l, rl)));
                                        sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, j, rj, pair(i, ri, j, rj)));
                                        sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, k, rk, pair(i, ri, k, rk)));
                                        sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, l, rl, pair(i, ri, l, rl)));
                                        sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", j, rj, k, rk, pair(j, rj, k, rk)));
                                        sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", j, rj, l, rl, pair(j, rj, l, rl)));
                                        sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", k, rk, l, rl, pair(k, rk, l, rl)));
                                        sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, j, rj, k, rk, triple(i, ri, j, rj, k, rk)));
                                        sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, j, rj, l, rl, triple(i, ri, j, rj, l, rl)));
                                        sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, k, rk, l, rl, triple(i, ri, k, rk, l, rl)));
                                        sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", j, rj, k, rk, l, rl, triple(j, rj, k, rk, l, rl)));
                                        sb.append(format("     backbone:                      %.3f\n", backboneEnergy));
                                        sb.append(format("     quadEnergy:                 %.3f\n", quadEnergy));
                                        sb.append(format("     --s--\n"));
                                        sb.append(format("     Active residues:\n"));
                                        for (int debug = 0; debug < residues.length; debug++) {
                                            if (residues[debug].getSideChainAtoms().get(0).getUse()) {
                                                sb.append(format("       %s\n", residues[debug].toString()));
                                            }
                                        }
                                        sb.append(format("     --f--\n"));
                                        logger.info(sb.toString());
                                    } else {
                                        logger.info(format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %s Ang.",
                                                resi, ri, resj, rj, resk, rk, resl, rl, quadEnergy, distString));
                                    }
                                } catch (ArithmeticException ex) {
                                    quadEnergy = 1E100;
                                    logger.info(format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: set to %10.6g (unreasonable conformation) at %s Ang.",
                                            resi, ri, resj, rj, resk, rk, resl, rl, quadEnergy, distString));
                                }
                                quadsEnergyCounter++;
                            }

                            // Revert rotamers and turn off atoms.
                            revertResidue(resi, resiOriginalCoordinates, resi.getRotamers(library)[0]);
                            revertResidue(resj, resjOriginalCoordinates, resj.getRotamers(library)[0]);
                            revertResidue(resk, reskOriginalCoordinates, resk.getRotamers(library)[0]);
                            revertResidue(resl, reslOriginalCoordinates, resl.getRotamers(library)[0]);
                        }
                    } else {
                        logger.info(format(" Quad %7s %-2d, %7s %-2d, %7s %-2d, %7s %-2d: %16.8f at %10.3f Ang.",
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
                            logger.warning(format(" Exception in distance loop: %s", ex.toString()));
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
                                    //logger.info(format(" %s %d %s %d %16.8f", residuei.toString(), i, residuej.toString(), j, r));
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

    private void multicastBuf(Buf message) {
        for (int p = 0; p < numProc; p++) {
            try {
                world.send(p, message);
            } catch (IOException ex) {
                logger.log(Level.WARNING, ex.getMessage(), ex);
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

        VideoWriter(MolecularAssembly molecularAssembly, boolean ignoreInactiveAtoms, boolean skipEnergies) {
            this.molecularAssembly = molecularAssembly;
            this.ignoreInactiveAtoms = ignoreInactiveAtoms;
            this.skipEnergies = skipEnergies;
            String molAssFile = molecularAssembly.getFile().getAbsolutePath();
            this.videoDir = FilenameUtils.separatorsToSystem(FilenameUtils.getFullPath(molAssFile) + "video/");
            this.filter = new PDBFilter(new File(videoDir + format("%08d", snapshotNum)), molecularAssembly, null, null);
            logger.info(format(" SDL: video snapshot directory: %s", videoDir));
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
            String name = format("%08d.pdb", snapshotNum);
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
                try {
                    return useFullAMOEBAEnergy ? currentEnergyWrapper(Arrays.asList(residues)) : computeEnergy(residues, currentRots, false);
                } catch (ArithmeticException ex) {
                    return 1E100;
                }
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
         * @param rots        Initial rotamer set
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

    public enum Algorithm {

        ALL, BOX, WINDOW, INDEPENDENT, BRUTE_FORCE
    }

    public enum Direction {

        FORWARD, BACKWARD
    }
}
