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
package ffx.algorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.BiFunction;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.mp.BooleanBuf;
import edu.rit.mp.Buf;
import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import edu.rit.pj.CommStatus;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.WorkerIntegerForLoop;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;
import edu.rit.pj.reduction.SharedDouble;

import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.MCMove;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.Potential;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.ObjectPair;
import static ffx.potential.bonded.Residue.ResidueType.NA;
import static ffx.potential.bonded.RotamerLibrary.applyRotamer;

import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;

// FastUtil libraries for large collections.

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

    /**
     * Logger for this class.
     */
    private static final Logger logger = Logger.getLogger(RotamerOptimization.class.getName());
    /**
     * MolecularAssembly to perform rotamer optimization on.
     */
    protected final MolecularAssembly molecularAssembly;
    /**
     * The Potential to evaluate during rotamer optimization.
     */
    protected Potential potential;
    /**
     * AlgorithmListener who should receive updates as the optimization runs.
     */
    protected AlgorithmListener algorithmListener;
    /**
     * List of Assemblies associated with a multi-topology Potential.
     */
    protected final List<MolecularAssembly> allAssemblies;
    /**
     * World Parallel Java communicator.
     */
    private final Comm world;
    /**
     * Number of Parallel Java processes.
     */
    private final int numProc;
    /**
     * Rank of this process.
     */
    private final int rank;
    /**
     * Flag to indicate if this is the master process.
     */
    private final boolean master;
    /**
     * Flag to indicate verbose logging.
     */
    private boolean verbose = false;
    /**
     * Flag to control the verbosity of printing.
     */
    private boolean print = false;
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
     * The potential energy of the system with all side-chains to be optimized
     * turned off.
     */
    protected double backboneEnergy;
    /**
     * Self-energy of each residue for each rotamer. [residue][rotamer]
     */
    protected double selfEnergy[][];
    /**
     * Two-Body cutoff distance.
     */
    protected double twoBodyCutoffDist = 0;
    // Fallback if there is no vdW node.
    private static final double FALLBACK_TWO_BODY_CUTOFF = 0;
    /**
     * Two-body energies for each pair of residues and pair of rotamers.
     * [residue1][rotamer1][residue2][rotamer2]
     */
    protected double twoBodyEnergy[][][][];
    /**
     * Flag to control use of 3-body terms.
     */
    protected boolean threeBodyTerm = true;
    /**
     * Flag to set three-body energies to zero outside of a cutoff.
     */
    private boolean threeBodyCutoff = true;
    /**
     * Three-body cutoff distance.
     */
    private double threeBodyCutoffDist;
    /**
     * Trimer-energies for each trimer of rotamers.
     * [residue1][rotamer1][residue2][rotamer2][residue3][rotamer3]
     */
    //protected double threeBodyEnergy[][][][][][];
    protected Object2DoubleMap<IntegerKeyset> threeBodyEnergies;
    /**
     * Quad cutoff flag.
     */
    private boolean quadCutoff = true;
    /**
     * Quad cutoff distance.
     */
    private double quadCutoffDist = 5.0;
    /**
     * Flag to prune clashes.
     */
    protected boolean pruneClashes = true;
    /**
     * Flag to prune pair clashes.
     */
    protected boolean prunePairClashes = true;
    /**
     * Eliminated rotamers. [residue][rotamer]
     */
    protected boolean eliminatedSingles[][];
    /**
     * *
     * Pruned rotamers. Only for JUnit testing purposes.
     */
    protected boolean onlyPrunedSingles[][];
    /**
     * Pruned rotamer pairs. Only for JUnit testing purposes.
     */
    protected boolean onlyPrunedPairs[][][][];
    /**
     * RotamerLibrary instance.
     */
    protected RotamerLibrary library = RotamerLibrary.getDefaultLibrary();
    /**
     * Number of permutations whose energy is explicitly evaluated.
     */
    private int evaluatedPermutations = 0;
    /**
     * Permutations are printed when modulo this field is zero.
     */
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
     * ONLY FOR UNIT TEST. DEFAULT VALUE IS TRUE. Turns off the singles
     * elimination criterion.
     */
    private boolean selfEliminationOn = true;
    /**
     * ONLY FOR UNIT TEST. DEFAULT VALUE IS TRUE. Turns off the pairs
     * elimination criterion.
     */
    private boolean pairEliminationOn = true;
    /**
     * The number of most-favorable structures to include as output.
     */
    private int ensembleNumber = 1;
    /**
     * The energy buffer applied to each elimination criteria to affect an
     * ensemble.
     */
    private double ensembleBuffer = 0.5;
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
     * more variation in self and 2-Body energies.
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
     * Eliminated rotamer pairs. [residue1][rotamer1][residue2][rotamer2]
     */
    private boolean eliminatedPairs[][][][];
    /**
     * False unless JUnit testing.
     */
    private boolean testing = false;
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
     * nucleicPairsPruningFactor is the arithmetic mean of 1.0 and the pruning
     * factor, and is applied for AA-NA pairs.
     */
    private double nucleicPruningFactor = 10.0;
    private double nucleicPairsPruningFactor = ((1.0 + nucleicPruningFactor) / 2);

    /**
     * Flag to calculate and print additional energies (mostly for debugging).
     */
    private boolean verboseEnergies = true;
    /**
     * A list of all residues being optimized. Note that Box and Window
     * optimizations operate on subsets of this list.
     */
    private ArrayList<Residue> allResiduesList = null;
    /**
     * An array of all residues being optimized. Note that Box and Window
     * optimizations operate on subsets of this list.
     */
    private Residue allResiduesArray[] = null;
    /**
     * Number of residues being optimized.
     */
    private int numResidues = 0;
    /**
     * The array of optimum rotamers for each residue.
     */
    private int optimum[];
    /**
     * Thread that receives self, 2-body, etc energies computed across nodes.
     */
    private ReceiveThread receiveThread;
    /**
     * Wait time, in milliseconds, for Thread.sleep() loops in rotamerEnergies()
     * and ReceiveThread
     */
    private final int POLLING_FREQUENCY = 500;
    /**
     * Thread that writes out a many-body expansion restart file.
     */
    private EnergyWriterThread energyWriterThread;
    /**
     * List of energy values ready to be written out by the EnergyWriterThread.
     */
    private List<String> energiesToWrite;
    /**
     * If true, write out an energy restart file.
     */
    private boolean writeEnergyRestart = true;
    /**
     * If true, load an energy restart file.
     */
    private boolean loadEnergyRestart = false;
    /**
     * Energy restart File instance.
     */
    private File energyRestartFile;
    /**
     * Partial restart file for box optimization.
     */
    private File partialFile;
    /**
     * ParallelTeam instance.
     */
    private ParallelTeam parallelTeam;
    /**
     * ParallelJava construct used for executing the various EnergyRegions.
     * Declared as a field so as to be reused during box optimization.
     */
    private final WorkerTeam energyWorkerTeam;
    /**
     * Map of self-energy values to compute.
     */
    private final HashMap<Integer, Integer[]> selfEnergyMap = new HashMap<>();
    /**
     * Flag to indicate ready to compute self-energy values.
     */
    private boolean readyForSingles = false;
    /**
     * Flag to indicate self energy computations are done.
     */
    private boolean finishedSelf = false;
    /**
     * Map of 2-body energy values to compute.
     */
    private final HashMap<Integer, Integer[]> twoBodyEnergyMap = new HashMap<>();
    /**
     * Flag to indicate ready to compute 2-body energy values.
     */
    private boolean readyFor2Body = false;
    /**
     * Flag to indicate 2-body energy computations are done.
     */
    private boolean finished2Body = false;
    /**
     * Map of 3-body energy values to compute.
     */
    private final HashMap<Integer, Integer[]> threeBodyEnergyMap = new HashMap<>();
    /**
     * Flag to indicate ready to compute 3-body energy values.
     */
    private boolean readyFor3Body = false;
    /**
     * Flag to indicate 3-body energy computations are done.
     */
    private boolean finished3Body = false;
    /**
     * Map of 4-body energy values to compute.
     */
    private final HashMap<Integer, Integer[]> fourBodyEnergyMap = new HashMap<>();
    /**
     * Flag to indicate computation of 4-body energy values. This is limited to
     * the study 4-body energy magnitudes, but is not included in the rotamer
     * optimization.
     */
    private boolean compute4BodyEnergy = false;
    /**
     * Flag to indicate ready to compute 4-body energy values.
     */
    private boolean readyFor4Body = false;
    /**
     * Maximum number of 4-body energy values to compute.
     */
    private int max4BodyCount = Integer.MAX_VALUE;
    /**
     * Parallel evaluation of quantities used during Goldstain Pair elimination.
     */
    private GoldsteinPairRegion goldsteinPairRegion;
    /**
     * Parallel evaluation of many-body energy sums.
     */
    private EnergyRegion energyRegion;
    /**
     * Flag to indicate use of box optimization.
     */
    private boolean usingBoxOptimization = false;
    /**
     * Number of boxes for box optimization in X, Y, Z.
     */
    private int[] numXYZBoxes = {3, 3, 3};
    /**
     * Box border size.
     */
    private double boxBorderSize = 0;
    /**
     * Approximate box size.
     */
    private double approxBoxLength = 0;
    /**
     * Box optimization inclusion criteria.
     */
    private int boxInclusionCriterion = 1;
    /**
     * Index of the first box to optimize.
     */
    private int boxStart = 0;
    /**
     * Index of the last box to optimize.
     */
    private int boxEnd = -1;
    /**
     * Flag to indicate manual definition of a super box.
     */
    private boolean manualSuperbox = false;
    /**
     * Dimensions of the box.
     */
    private double[] boxDimensions;
    /**
     * Buffer size for the super box.
     */
    private double superboxBuffer = 8.0;
    /**
     * If a pair of residues have two atoms closer together than the
     * superposition threshold, the energy is set to NaN.
     */
    private double superpositionThreshold = 0.1;
    /**
     * Box index loaded during a restart.
     */
    private int boxLoadIndex = -1;
    /**
     * Box indeces loaded during a restart.
     */
    private int[] boxLoadCellIndices;
    /**
     * Flag to indicate use of forced residues during the sliding window
     * algorithm.
     */
    private boolean useForcedResidues = false;
    /**
     * Beginning of the forced residue range.
     */
    private int startForcedResidues = -1;
    /**
     * End of the force residue range.
     */
    private int endForcedResidues = -1;
    /**
     * Flag to indicate computation of a many-body expansion for original
     * coordinates.
     */
    private boolean decomposeOriginal = false;
    /**
     * Use original side-chain coordinates as a rotamer for each residue.
     */
    private boolean addOrigRot = false;
    /**
     * Flag to indicate use of MC optimization.
     */
    private boolean monteCarlo = false;
    /**
     * Number of MC optimization steps.
     */
    private int nMCsteps = 1000000;
    /**
     * MC temperature (K).
     */
    private double mcTemp = 298.15;
    /**
     * Check to see if proposed move has an eliminated 2-body or higher-order
     * term.
     */
    private boolean mcUseAll = false;
    /**
     * Skips brute force enumeration in favor of pure Monte Carlo. Recommended
     * only for testing.
     */
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
    /**
     * Represents the method called to obtain the directory corresponding to the
     * current energy; will be a simple return null for potential energy
     * evaluations. While current energy calls will fill the rotamer list with
     * the current rotamers of the residue, other methods may skip applying the
     * rotamer directly.
     */
    private BiFunction<List<Residue>, List<Rotamer>, File> dirSupplier;
    /**
     * Represents the method called to obtain energy for the current rotamer or
     * state; defaults to the existing potential energy code. May discard the
     * input file.
     */
    private ToDoubleFunction<File> eFunction;
    /**
     * Default value for maxRotCheckDepth.
     */
    private int defaultMaxRotCheckDepth = 2;
    /**
     * Maximum depth to check if a rotamer can be eliminated.
     */
    private int maxRotCheckDepth = defaultMaxRotCheckDepth;
    /**
     * Test Self-Energy Elimination.
     */
    private boolean testSelfEnergyEliminations = false;
    /**
     * Test Pair-Energy Elimination.
     * <p>
     * If greater than or equal to 0, test the specified residue.
     */
    private int testPairEnergyEliminations = -1;
    /**
     * Test Triple-Energy Elimination.
     * <p>
     * If greater than or equal to 0, test the specified residues.
     */
    private int testTripleEnergyEliminations1 = -1;
    /**
     * Test Triple-Energy Elimination.
     * <p>
     * If greater than or equal to 0, test the specified residues.
     */
    private int testTripleEnergyEliminations2 = -1;
    /**
     * Holds onto the "default" low-correction-vector Rotamer for a nucleic acid
     * Residue.
     */
    private Map<Residue, Rotamer> defaultNucleicRotamers;
    /**
     * Writes residue and rotamer numbers of eliminated singles and pairs to a
     * CSV file.
     */
    public boolean outputFile[];
    /**
     * Stores eliminated residue and rotamer singles and pairs.
     */
    public List<Integer> eliminatedResidue = new ArrayList<Integer>();
    public List<Integer> eliminatedRotamer = new ArrayList<Integer>();
    public List<Integer> eliminatedResidue1 = new ArrayList<Integer>();
    public List<Integer> eliminatedRotamer1 = new ArrayList<Integer>();
    public List<Integer> eliminatedResidue2 = new ArrayList<Integer>();
    public List<Integer> eliminatedRotamer2 = new ArrayList<Integer>();
    /**
     * Stores the number of box iterations for single and pair elimination.
     */
    public List<Integer> singlesIterations = new ArrayList<Integer>();
    public List<Integer> pairsIterations = new ArrayList<Integer>();

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
        if (potential instanceof ForceFieldEnergy) {
            ((ForceFieldEnergy) potential).setPrintOnFailure(false, true);
        } else if (potential instanceof DualTopologyEnergy) {
            ((DualTopologyEnergy) potential).setPrintOnFailure(false, true);
        } else if (potential instanceof QuadTopologyEnergy) {
            ((QuadTopologyEnergy) potential).setPrintOnFailure(false, true);
        }
        this.algorithmListener = algorithmListener;
        eFunction = this::currentPE;
        dirSupplier = (List<Residue> resList, List<Rotamer> rotList) -> null;

        world = Comm.world();
        numProc = world.size();
        rank = world.rank();

        if (rank == 0) {
            master = true;
        } else {
            master = false;
        }

        energyWorkerTeam = new WorkerTeam(world);

        if (System.getProperty("verbose") != null) {
            if (System.getProperty("verbose").equalsIgnoreCase("true")) {
                verbose = true;
            }
        }

        /**
         * Set the default 2-body Cutoff to the van der Waals cutoff.
         */
        ForceFieldEnergy forceFieldEnegy = molecularAssembly.getPotentialEnergy();
        VanDerWaals vdW = forceFieldEnegy.getVdwNode();
        if (vdW != null) {
            NonbondedCutoff nonBondedCutoff = vdW.getNonbondedCutoff();
            twoBodyCutoffDist = nonBondedCutoff.off;
        } else {
            twoBodyCutoffDist = FALLBACK_TWO_BODY_CUTOFF;
        }

        // Process relevant system keys.
        String undo = System.getProperty("ro-undo");
        String direction = System.getProperty("ro-direction");
        String increment = System.getProperty("ro-increment");
        String goldstein = System.getProperty("ro-goldstein");
        String superpositionThreshold = System.getProperty("ro-superpositionThreshold");
        String ensembleNumber = System.getProperty("ro-ensembleNumber");
        String ensembleEnergy = System.getProperty("ro-ensembleEnergy");
        String ensembleBuffer = System.getProperty("ro-ensembleBuffer");
        String threeBodyCutoffDist = System.getProperty("ro-threeBodyCutoffDist");
        String nucleicPruningFactor = System.getProperty("ro-nucleicPruningFactor");
        String nucleicCorrectionThreshold = System.getProperty("ro-nucleicCorrectionThreshold");
        String minimumNumberAcceptedNARotamers = System.getProperty("ro-minimumNumberAcceptedNARotamers");
        String singletonClashThreshold = System.getProperty("ro-singletonClashThreshold");
        String multiResClashThreshold = System.getProperty("ro-multiResClashThreshold");
        String pairClashThreshold = System.getProperty("ro-pairClashThreshold");
        String multiResPairClashAddition = System.getProperty("ro-multiResPairClashAddition");
        String boxDimensions = System.getProperty("ro-boxDimensions");
        String computeQuads = System.getProperty("ro-compute4BodyEnergy");
        String quadCutoffDist = System.getProperty("ro-quadCutoffDist");
        String quadMaxout = System.getProperty("ro-max4BodyCount");
        String lazyMatrix = System.getProperty("ro-lazyMatrix");
        String mcTemp = System.getProperty("ro-mcTemp");
        String mcUseAll = System.getProperty("ro-mcUseAll");
        String mcNoEnum = System.getProperty("ro-debug-mcNoEnum");
        String addOrigRotStr = System.getProperty("ro-addOrigRot");
        String origAtEndStr = System.getProperty("ro-origAtEnd");

        if (computeQuads != null) {
            boolean value = Boolean.parseBoolean(computeQuads);
            this.compute4BodyEnergy = value;
            logger.info(format(" (KEY) compute4BodyEnergy: %b", this.compute4BodyEnergy));
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
            this.max4BodyCount = value;
            logger.info(format(" (KEY) max4BodyCount: %d", this.max4BodyCount));
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
        if (nucleicPruningFactor != null) {
            double value = Double.parseDouble(nucleicPruningFactor);
            this.nucleicPruningFactor = (value >= 0 ? value : 1.0);
            this.nucleicPairsPruningFactor = (1.0 + value) / 2;
            logger.info(format(" (KEY) nucleicPruningFactor: %.2f", this.nucleicPruningFactor));
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

        String prop = "ro-maxRotCheckDepth";
        String propStr = System.getProperty(prop);
        if (propStr != null) {
            try {
                maxRotCheckDepth = Integer.parseInt(propStr);
                if (maxRotCheckDepth > 3 || maxRotCheckDepth < 0) {
                    throw new IllegalArgumentException(" ro-maxRotSearchDepth must be between 0-3 inclusive!");
                }
            } catch (Exception ex) {
                maxRotCheckDepth = defaultMaxRotCheckDepth;
                logger.warning(String.format(" Could not parse %s value %s as valid integer; defaulting to %d", prop, propStr, maxRotCheckDepth));
                logger.warning(String.format(" Exception: %s", ex));
            }
        } else {
            maxRotCheckDepth = defaultMaxRotCheckDepth;
        }

        prop = System.getProperty("revertUnfavorable", "false");
        revert = Boolean.parseBoolean(prop);

        allAssemblies = new ArrayList<>();
        allAssemblies.add(molecularAssembly);
    }

    /**
     * Set the "use" flag to true for all variable atoms in a residue.
     *
     * @param residue The Residue to turn on.
     */
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
                break;
        }
    }

    /**
     * Set the "use" flag to true for all variable atoms in a residue.
     *
     * @param residue The residue to turn off.
     */
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
                break;
        }
    }

    /**
     * Control the depth of self-consistency checking with a rotamer is
     * eliminated.
     *
     * @param maxRotCheckDepth
     */
    public void setMaxRotCheckDepth(int maxRotCheckDepth) {
        this.maxRotCheckDepth = maxRotCheckDepth;
    }

    /**
     * Add an additional molecular assembly to the internal list.
     *
     * @param assembly
     */
    public void addAdditionalAssembly(MolecularAssembly assembly) {
        allAssemblies.add(assembly);
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

                double rotEnergy = Double.NaN;
                try {
                    rotEnergy = currentEnergy(resList);
                    logger.info(format(" %d Energy: %s", ++evaluatedPermutations, formatEnergy(rotEnergy)));
                } catch (ArithmeticException ex) {
                    logger.info(String.format(" %d Energy set to NaN (unreasonable conformation)", ++evaluatedPermutations));
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
                if (check(i, ri)) {
                    continue;
                }
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
                if (check(i, ri)) {
                    continue;
                }
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
     * Checks the pair elimination array to see if this permutation has been
     * eliminated.
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
                    double amoebaEnergy = Double.NaN;
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
                    if (threeBodyTerm) {;
                        logIfMaster(format(" %12s %25f %25f", evaluatedPermutations, approximateEnergy, lowEnergy));
                        //logIfMaster(format(" %6e 3-Body: %12.4f (%12.4f)",
                        //       (double) evaluatedPermutations, approximateEnergy, lowEnergy));
                    } else {
                        logIfMaster(format(" %12s %25f %25f", evaluatedPermutations, approximateEnergy, lowEnergy));
                        //logIfMaster(format(" %6e 2-Body: %12.4f (%12.4f)",
                        //        (double) evaluatedPermutations, approximateEnergy, lowEnergy));
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
     * @param i Current depth in residue/rotamer tree.
     * @param currentRotamers Current set of rotamers at this node.
     * @param gmecEnergy Minimum energy for these residues.
     * @param permutationEnergies Energy of all permutations.
     * @param permutations Contains accepted permutations.
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

    /**
     * Flag to control use of 3-body energy terms.
     *
     * @param threeBodyTerm
     */
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
     * Uses existing backbone, self, 2-Body, and 3-body energies from
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
            pairSum = energyRegion.getTwoBody();
            threeBodySum = energyRegion.getThreeBody();
        } catch (Exception e) {
            //logger.log(Level.WARNING, " Exception in EnergyRegion.", e);
            throw new IllegalArgumentException(e);
        }

        double approximateEnergy = backboneEnergy + selfSum + pairSum + threeBodySum;
        if (print) {
            logger.info(format(" Backbone                  %s", formatEnergy(backboneEnergy)));
            logger.info(format(" Self Energy               %s", formatEnergy(selfSum)));
            logger.info(format(" Pair Energy               %s", formatEnergy(pairSum)));
            if (!threeBodyTerm) {
                logger.info(format(" Total Energy up to 2-Body %s", formatEnergy(approximateEnergy)));
            } else {
                logger.info(format(" 3-Body Energy             %s", formatEnergy(threeBodySum)));
                logger.info(format(" Total Energy up to 3-Body %s", formatEnergy(approximateEnergy)));
            }
        }
        return approximateEnergy;
    }

    /**
     * Sets the decompose-original flag.
     *
     * @param decomposeOriginal If true, decompose the energy of the structure
     * without optimizing.
     */
    public void setDecomposeOriginal(boolean decomposeOriginal) {
        this.decomposeOriginal = decomposeOriginal;
    }

    /**
     * Test the self-energy elimination by setting 2-body and 3-body
     * interactions to zero.
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
                                        //threeBodyEnergy[i][ri][j][rj][k][rk] = 0.0;
                                        threeBodyEnergies.put(new IntegerKeyset(i, ri, j, rj, k, rk), 0.0);
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
     * Utility method for formatting energies, using 16 spaces with 8 digits of
     * precision.
     *
     * @param energy Energy to format.
     * @return A string representing the energy.
     */
    private String formatEnergy(double energy) {
        if (abs(energy) < 1.0e6) {
            return format("%16.8f", energy);
        } else {
            return format("*%15.4e", energy);
        }
        // TODO: Possibly replace with %16.8g, which so far as I know, is equivalent.
    }

    /**
     * Test the elimination criteria by setting self and 3-body interactions to
     * zero.
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
                                        //threeBodyEnergy[i][ri][j][rj][k][rk] = 0.0;
                                        threeBodyEnergies.put(new IntegerKeyset(i, ri, j, rj, k, rk), 0);
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
     * Test the elimination criteria by setting self and 2-body interactions to
     * zero. Two residues are at fixed rotamers and all rotamer interactions
     * with those two residues are calculated.
     *
     * @param residues
     * @param resID1 The residue number for one of two fixed residues.
     * @param resID2 The second residue number for one of two fixed residues.
     */
    public void testTripleEnergyElimination(Residue residues[], int resID1, int resID2) {
        int nRes = residues.length;

        if (resID1 >= nRes) {
            return;
        }
        if (resID2 >= nRes) {
            return;
        }
        if (resID1 == resID2) {
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
                        /**
                         * if (i != resID1 && j != resID1) { try {
                         * twoBodyEnergy[i][ri][j][rj] = 0.0; } catch (Exception
                         * e) { // catch NPE. } }
                         */
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

                                    if (i != resID1 && j != resID1 && k != resID1) {
                                        try {
                                            //threeBodyEnergy[i][ri][j][rj][k][rk] = 0.0;
                                            threeBodyEnergies.put(new IntegerKeyset(i, ri, j, rj, k, rk), 0);
                                        } catch (Exception e) {
                                            // catch NPE.
                                        }
                                    }
                                    if (i != resID2 && j != resID2 && k != resID2) {
                                        try {
                                            //threeBodyEnergy[i][ri][j][rj][k][rk] = 0.0;
                                            threeBodyEnergies.put(new IntegerKeyset(i, ri, j, rj, k, rk), 0);
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

    }

    /**
     * Prints based on a molecule-based decomposition of original energies.
     * Would need to be re-implemented.
     *
     * @param residues
     * @param totalEnergy
     * @param backbone
     * @param residueEnergy
     */
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
         * Sum self, 2-Body and 3-Body energies.
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
                    residues[i].toFormattedString(false, true), residueEnergy[0][i], residueEnergy[1][i], residueEnergy[2][i], total));
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
                                residues[i].toFormattedString(false, true), residues[j].toFormattedString(false, true), twoBodyEnergy[i][0][j][0]));
                    }
                }
            }
            logger.info(format("\n Large trimer interactions (>%.2f):", trimerCutoff));
            for (int i = 0; i < nRes; i++) {
                for (int j = i + 1; j < nRes; j++) {
                    for (int k = j + 1; k < nRes; k++) {
                        if (Math.abs(get3Body(residues, i, 0, j, 0, k, 0)) >= trimerCutoff) {
                            logger.info(format(" Large Trimer  %s %s %s:    %16.5f",
                                    residues[i].toFormattedString(false, true), residues[j].toFormattedString(false, true), residues[k].toFormattedString(false, true), get3Body(residues, i, 0, j, 0, k, 0)));
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
                                logger.info(format(" Large Pair %8s %-2d, %8s %-2d: %s",
                                        resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj,
                                        formatEnergy(twoBodyEnergy[i][ri][j][rj])));
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
                                    if (Math.abs(get3Body(residues, i, ri, j, rj, k, rk)) >= trimerCutoff) {
                                        logger.info(format(" Large Trimer %8s %-2d, %8s %-2d, %8s %-2d: %s",
                                                resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk,
                                                formatEnergy(get3Body(residues, i, ri, j, rj, k, rk))));
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
            initResidue(residue, polymer, i);
        }
    }

    /**
     * Initialize a rotamer for optimization: add it to residueList, apply its
     * 0th rotamer, initialize default atomic doordinates, etc.
     *
     * @param residue A Residue to add to optimization.
     * @param polymer The Polymer it belongs to.
     * @param i residues index in polymer.
     */
    private void initResidue(Residue residue, Polymer polymer, int i) {
        if (residue == null) {
            logger.warning(String.format(" Null residue %d for chain %c", i, polymer.getChainID()));
        } else {
            Rotamer[] rotamers = residue.getRotamers(library);
            if (rotamers != null) {
                int lenri = rotamers.length;
                if (lenri > 1 || addOrigRot) {
                    residue.initializeDefaultAtomicCoordinates();
                    RotamerLibrary.applyRotamer(residue, rotamers[0]);
                    residueList.add(residue);
                }
            } else if (useForcedResidues && checkIfForced(i)) {
                residueList.add(residue);
            }
        }
    }
    
    /**
     * Returns the restart file.
     * @return energyRestartFile File with saved side-chain energies.
     */
    public File getRestartFile(){
        return energyRestartFile;
    }
    /**
     * Returns the partial file.
     * @return partialFile File with saved box optimization energies.
     */
    public File getPartial(){
        return partialFile;
    }

    /**
     * Return the residue list.
     *
     * @return
     */
    public ArrayList<Residue> getResidues() {
        return residueList;
    }

    /**
     * Set the residue list.
     *
     * @param residueList
     */
    public void setResidues(ArrayList<Residue> residueList) {
        this.residueList = residueList;
    }

    /**
     * Return an integer array of optimized rotamers following rotamer
     * optimization.
     *
     * @return The optimal rotamer array.
     */
    public int[] getOptimumRotamers() {
        return optimum;
    }

    public void setCoordinatesToEnsemble(int ensnum) {
        if (ensembleStates != null && !ensembleStates.isEmpty()) {
            ensnum %= ensembleStates.size();
            ResidueState.revertAllCoordinates(residueList, ensembleStates.get(ensnum).getVal());
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

    public void setEnsemble(int ensemble) {
        setEnsemble(ensemble, 5.0);
    }

    /**
     * Accepts a list of residues but throws out null residues. Used by the -lR
     * flag.
     *
     * @param residueList
     * @return Added residues.
     */
    public List<Residue> setResiduesIgnoreNull(List<Residue> residueList) {
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

    /**
     * Perform the rotamer optimization using the specified algorithm.
     *
     * @param algorithm
     * @return the lowest energy found.
     */
    public double optimize(Algorithm algorithm) {
        this.algorithm = algorithm;
        return optimize();
    }

    /**
     * Execute the rotamer optimization.
     *
     * @return The lowest energy found.
     */
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

        sortResidues(allResiduesList);
        sortResidues(residueList);

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
         * based on distance cutoff, and for cutoff distances.
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
     * Eliminates NA backbone rotamers with corrections greater than threshold,
     * and sets the "default" rotamer with the minimal correction.
     * <p>
     * A nucleic correction threshold <= 0 will cause
     *
     * @param residues Residues to eliminate bad backbone rotamers over.
     */
    private void eliminateNABackboneRotamers(Residue[] residues) {
        double eliminateBy = Double.MAX_VALUE;
        if (nucleicCorrectionThreshold > 0) {
            eliminateBy = nucleicCorrectionThreshold;
            logIfMaster(format(" Eliminating nucleic acid rotamers with correction vectors larger than %5.3f A", nucleicCorrectionThreshold));
            logIfMaster(format(" A minimum of %d rotamers per NA residue will carry through to energy calculations.", minNumberAcceptedNARotamers));
        }
        ArrayList<Residue> resList = new ArrayList<>();
        resList.addAll(Arrays.asList(residues));
        int numR = residues.length;
        defaultNucleicRotamers = new HashMap<>(numR);

        for (int i = 0; i < numR; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            if (residue.getResidueType() == NA && rotamers != null) {
                int lenri = rotamers.length;
                DoubleIndexPair[] sortedCorrections = new DoubleIndexPair[lenri];

                for (int ri = 0; ri < lenri; ri++) {
                    Rotamer rotamer = rotamers[ri];
                    sortedCorrections[ri] = new DoubleIndexPair(ri, RotamerLibrary.applyNARotamer(residue, rotamer, false));
                }

                Arrays.sort(sortedCorrections);
                int numAccepted = 0;
                for (int indRI = 0; indRI < lenri; indRI++) {
                    DoubleIndexPair correction = sortedCorrections[indRI];
                    int ri = correction.getIndex();
                    double corrMag = correction.getDoubleValue();

                    if (check(i, ri)) {
                        continue;
                    }
                    if (corrMag < eliminateBy || numAccepted < minNumberAcceptedNARotamers) {
                        if (numAccepted++ == 0) {
                            defaultNucleicRotamers.put(residue, rotamers[ri]);
                            logger.info(String.format(" Setting the default rotamer for residue %s to %d",
                                    residue.toFormattedString(true, true), ri));
                            assert ri == 0 || !library.getUsingOrigCoordsRotamer();
                        }
                    } else {
                        eliminateRotamer(residues, i, ri, print);
                        logIfMaster(String.format(" Eliminated nucleic acid rotamer %s-%d that required a "
                                + "linking correction %6.4f A > %6.4f A", residue.toFormattedString(false, true),
                                ri, corrMag, eliminateBy));
                    }
                }
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
     */
    private void reconcileNARotamersWithPriorResidues(Residue[] residues) {
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
                for (int ri = 0; ri < rotamers.length; ri++) {
                    Rotamer roti = rotamers[ri];
                    if (!roti.isState && RotamerLibrary.checkPucker(roti.chi1) != 1) {
                        if (print) {
                            logIfMaster(format(" Rotamer %d of residue %s eliminated "
                                    + "for incompatibility with the sugar pucker of previous "
                                    + "residue %s outside the window.", ri, residuei.toFormattedString(false, true),
                                    prevResidue.toString()));
                        }
                        eliminateRotamer(residues, i, ri, print);
                    }
                }
            } else {
                for (int ri = 0; ri < rotamers.length; ri++) {
                    Rotamer roti = rotamers[ri];
                    if (!roti.isState && RotamerLibrary.checkPucker(roti.chi1) != 2) {
                        if (print) {
                            logIfMaster(format(" Rotamer %d of residue %s eliminated "
                                    + "for incompatibility with the sugar pucker of previous "
                                    + "residue %s outside the window.", ri, residuei.toFormattedString(false, true),
                                    prevResidue.toString()));
                        }
                        eliminateRotamer(residues, i, ri, print);
                    }
                }
            } // TODO: Implement support for the DNA C3'-exo pucker.
        }
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
                                            + "residue %s outside the window.", j, residuei.toFormattedString(false, true),
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
                                            + "residue %s outside the window.", j, residuei.toFormattedString(false, true),
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
     * Specify use of Goldstein optimization.
     *
     * @param useGoldstein
     */
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

    /**
     * Set the starting box index.
     *
     * @param boxStart
     */
    public void setBoxStart(int boxStart) {
        this.boxStart = boxStart;
    }

    /**
     * Set the ending box index.
     *
     * @param boxEnd
     */
    public void setBoxEnd(int boxEnd) {
        // Is -1 if boxes run to completion.
        this.boxEnd = boxEnd;
    }

    /**
     * Set the optimization direction to forward or backward.
     *
     * @param direction
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
    }

    /**
     * Set the cut-off distance for inclusion of residues in sliding box and
     * window methods.
     *
     * @param distance
     */
    public void setDistanceCutoff(double distance) {
        this.distance = distance;
    }

    /**
     * Set the residue increment for sliding window.
     *
     * @param increment
     */
    public void setIncrement(int increment) {
        this.increment = increment;
    }

    /**
     * Set the algorithm to revert to starting coordinates if the energy
     * increases.
     *
     * @param revert
     */
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

    /**
     * The the nucleic acid correction threshold.
     *
     * @param nucleicCorrectionThreshold
     */
    public void setNucleicCorrectionThreshold(double nucleicCorrectionThreshold) {
        if (nucleicCorrectionThreshold >= 0) {
            this.nucleicCorrectionThreshold = nucleicCorrectionThreshold;
        } else {
            logger.warning("\n Correction threshold must be >= 0. Setting to default of 0 (threshold inactive).\n");
            this.nucleicCorrectionThreshold = 0;
        }
    }

    /**
     * Set the minimum number of accepted nucleic acid rotamers.
     *
     * @param minNumberAcceptedNARotamers
     */
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
     * @param nucleicPruningFactor
     */
    public void setNucleicPruningFactor(double nucleicPruningFactor) {
        this.nucleicPruningFactor = nucleicPruningFactor;
        this.nucleicPairsPruningFactor = ((1.0 + nucleicPruningFactor) / 2);
    }

    /**
     * Set the verbose energy flag.
     *
     * @param verboseEnergies
     */
    public void setVerboseEnergies(Boolean verboseEnergies) {
        this.verboseEnergies = verboseEnergies;
    }

    /**
     * Sets whether rotamer optimization should print out any files, or act
     * solely to optimize a structure in memory.
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
                double newE = Double.NaN;
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

        int iterations = 0;
        boolean finalTry = false;
        int bestEnsembleTargetDiffThusFar = Integer.MAX_VALUE;
        double bestBufferThusFar = ensembleBuffer;
        double startingBuffer = ensembleBuffer;

        optimum = new int[nResidues];

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
            double afterPairElim = singletonPermutations - pairTotalElimination;
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
            logIfMaster(format("%30s %35s %35s", "Condition",  "Number of Permutations Left", "Number of Permutations Removed"));
            logIfMaster(format("%30s %35s %35s", "No Eliminations", permutations, ""));
            logIfMaster(format("%30s %35s %35s", "Single Eliminations", singletonPermutations, permutations - singletonPermutations));
            logIfMaster(format("%30s %35s %35s", "Pair Eliminations", afterPairElim, pairTotalElimination));
            logIfMaster(format("%30s %35s %35s", "Single and Pair Eliminations", (double) evaluatedPermutations, pairTotalElimination + (permutations - singletonPermutations)));

            logIfMaster(format("\n Energy of permutations:"));
            logIfMaster(format(" %12s %25s %25s", "Permutation", "Energy", "Lowest Possible Energy"));

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

            logIfMaster(format(" Final rotamers:"));
            logIfMaster(format("%17s %10s %11s %12s %11s", "Residue", "Chi 1", "Chi 2", "Chi 3", "Chi 4"));
            for (int i = 0; i < nResidues; i++) {
                Residue residue = residues[i];
                Rotamer[] rotamers = residue.getRotamers(library);
                int ri = optimum[i];
                Rotamer rotamer = rotamers[ri];
                logIfMaster(format(" %c (%7s,2d) %s", residue.getChainID(), residue, ri, rotamer.toAngleString()));
                RotamerLibrary.applyRotamer(residue, rotamer);
            }
            logIfMaster(format("\n"));

            double sumSelfEnergy = 0;
            double sumPairEnergy = 0;
            double sumTrimerEnergy = 0;
            for (int i = 0; i < nResidues; i++) {
                int ri = optimum[i];
                sumSelfEnergy += getSelf(i, ri);
                logIfMaster(format(" Final self Energy (%8s,%2d): %12.4f", residues[i].toFormattedString(false, true), ri, getSelf(i, ri)));
            }
            for (int i = 0; i < nResidues - 1; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues; j++) {
                    int rj = optimum[j];
                    sumPairEnergy += get2Body(i, ri, j, rj);
                    if (get2Body(i, ri, j, rj) > 10.0) {
                        logIfMaster(format(" Large Final Pair Energy (%8s,%2d) (%8s,%2d): %12.4f", residues[i].toFormattedString(false, true), ri,
                                residues[j].toFormattedString(false, true), rj, get2Body(i, ri, j, rj)));
                    }
                }
            }

            try {
                e = currentEnergy(residueList);
            } catch (ArithmeticException ex) {
                e = Double.NaN;
                logger.severe(String.format(" Exception %s in calculating current energy at the end of triples", ex.toString()));
            }

            logIfMaster(format(" %12s %25s %25s", "Type", "Energy", "Lowest Possible Energy"));
            logIfMaster(format(" %12s %25f %25s", "Self:", sumSelfEnergy, ""));
            logIfMaster(format(" %12s %25f %25s", "Pair:", sumPairEnergy, ""));

            double approximateEnergy = backboneEnergy + sumSelfEnergy + sumPairEnergy;

            if (threeBodyTerm) {
                for (int i = 0; i < nResidues - 2; i++) {
                    int ri = optimum[i];
                    for (int j = i + 1; j < nResidues - 1; j++) {
                        int rj = optimum[j];
                        for (int k = j + 1; k < nResidues; k++) {
                            int rk = optimum[k];
                            try {
                                sumTrimerEnergy += get3Body(residues, i, ri, j, rj, k, rk);
                            } catch (Exception ex) {
                                logger.warning(ex.toString());
                            }
                        }
                    }
                }
                approximateEnergy += sumTrimerEnergy;
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - backboneEnergy;
                logIfMaster(format(" %12s %5s %25f %5s %25s %5s", "Trimer:", sumTrimerEnergy, ""));
                logIfMaster(format(" %12s %5s %25f %5s %25s %5s", "Neglected:", higherOrderEnergy, ""));
            } else {
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - backboneEnergy;
                logIfMaster(format(" %12s %5s %25f %5s %25s %5s", "Neglected:", higherOrderEnergy, ""));
            }
            logIfMaster(format(" %12s %25f %25s", "Approximate:", approximateEnergy, ""));
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

            logIfMaster(format(" Collecting Permutations:"));
            dryRun(residues, 0, currentRotamers);

            double pairTotalElimination = singletonPermutations - (double) evaluatedPermutations;
            double afterPairElim = singletonPermutations - pairTotalElimination;
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
                logIfMaster(format("%30s %35s %35s", "Condition", "Number of Permutations Left", "Number of Permutations Removed"));
                logIfMaster(format("%30s %35s %35s", "No Eliminations", permutations, ""));
                logIfMaster(format("%30s %35s %35s", "Single Eliminations", singletonPermutations, permutations - singletonPermutations));
                logIfMaster(format("%30s %35s %35s", "Pair Eliminations", afterPairElim, pairTotalElimination));
                logIfMaster(format("%30s %35s %35s", "Single and Pair Eliminations", (double) evaluatedPermutations, pairTotalElimination + (permutations - singletonPermutations)));

                logIfMaster(format("\n Energy of permutations:"));
                logIfMaster(format(" %12s %25s %25s", "Permutation", "Energy", "Lowest Possible Energy"));

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

        logIfMaster(format("\n Energy contributions:"));
        logIfMaster(format(" %15s %25s %25s", "Type", "Energy", "Lowest Possible Energy"));

        for (int i = 0; i < nResidues; i++) {
            int ri = optimum[i];
            Residue residue = residues[i];
            Rotamer rotamers[] = residue.getRotamers(library);
            turnOnAtoms(residue);
            RotamerLibrary.applyRotamer(residue, rotamers[ri]);
            double self = getSelf(i, ri);
            residueEnergy[i] = self;
            sumSelfEnergy += self;
            double lowest = lowestSelfEnergy(residues, i);
            sumLowSelfEnergy += lowest;
            if (self - lowest > 10.0) {
                logIfMaster(format(" %15s %25f %25f", "Self (" + residues[i] + "," + ri + "):", self, lowest));
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
                double pair = get2Body(i, ri, j, rj);
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
                    sb.append(format("  Pair Energy (%8s,%2d) (%8s,%2d): %12.4f (Lowest: %12.4f).\n",
                            residues[i].toFormattedString(false, true), ri, residues[j].toFormattedString(false, true), rj, pair, lowest));
                }
            }
            if (sumPairEnergyI - sumLowPairEnergyI > 10.0) {
                logIfMaster(format(" %15s %25f %25f", "Self (" + residues[i] + "," + ri + ")", sumPairEnergyI, sumLowPairEnergyI));
                sb.trimToSize();
                if (!sb.toString().isEmpty()) {
                    logIfMaster(sb.toString());
                }
            }
        }

        double e = Double.NaN;
        try {
            e = currentEnergy(residueList);
        } catch (ArithmeticException ex) {
            logger.severe(String.format(" Exception %s in calculating current energy at the end of self and pairs", ex.toString()));
        }
        logIfMaster(format(" %15s %25f %25s", "Backbone:", backboneEnergy, ""));
        logIfMaster(format(" %15s %25f %25f", "Self:",  sumSelfEnergy, sumLowSelfEnergy));
        logIfMaster(format(" %15s %25f %25f", "Pair:", sumPairEnergy, sumLowPairEnergy));

        double approximateEnergy = backboneEnergy + sumSelfEnergy + sumPairEnergy;

        double sumTrimerEnergy = 0;
        if (threeBodyTerm) {
            for (int i = 0; i < nResidues - 2; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues - 1; j++) {
                    int rj = optimum[j];
                    for (int k = j + 1; k < nResidues; k++) {
                        int rk = optimum[k];
                        try {
                            double triple = get3Body(residues, i, ri, j, rj, k, rk);
                            double thirdTrip = triple / 3.0;
                            residueEnergy[i] += thirdTrip;
                            residueEnergy[j] += thirdTrip;
                            residueEnergy[k] += thirdTrip;
                            sumTrimerEnergy += triple;
                        } catch (Exception ex) {
                            logger.warning(ex.toString());
                        }
                    }
                }
            }
            approximateEnergy += sumTrimerEnergy;
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - backboneEnergy;
            logIfMaster(format(" %15s %25f %25s", "Trimer:", sumTrimerEnergy, ""));
            logIfMaster(format(" %15s %25f %25s", "Neglected:", higherOrderEnergy, ""));
        } else {
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - backboneEnergy;
            logIfMaster(format(" %15s %25f %25s", "Neglected:", higherOrderEnergy, ""));
        }

        logIfMaster(format(" %15s %25f %25s", "Approximate:", approximateEnergy, ""));

        logIfMaster(format("\n Final rotamers:"));
        logIfMaster(format("%17s %10s %11s %12s %11s %14s", "Residue", "Chi 1", "Chi 2", "Chi 3", "Chi 4", "Energy"));
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            int ri = optimum[i];
            Rotamer rotamer = rotamers[ri];
            logIfMaster(format(" %3d %c (%7s,%2d) %s %12.4f ",
                    i + 1, residue.getChainID(), residue, ri, rotamer.toAngleString(), residueEnergy[i]));
            RotamerLibrary.applyRotamer(residue, rotamer);
        }
        logIfMaster(format("\n"));
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
                double e = getSelf(i, ni);
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
                double e = get2Body(i, ri, j, jr);
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
     * of residues.
     *
     * @param residueList Residues to be optimized.
     * @return Global minimum energy conformation energy.
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
         * Compute the number of permutations without eliminating dead-ends.
         */
        double permutations = 1;
        for (int i = 0; i < nResidues; i++) {
            Residue residue = residues[i];
            Rotamer[] rotamers = residue.getRotamers(library);
            permutations *= rotamers.length;
        }

        logger.info(format(" Number of permutations: %16.8e.", permutations));
        double e;
        useFullAMOEBAEnergy = false;

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
     * Default distance method is to find the shortest distance between
     * residues.
     */
    private DistanceMethod distanceMethod = DistanceMethod.RESIDUE;

    /**
     * Allows get2BodyDistance to find the shortest distance between any two
     * rotamers or two residues.
     */
    public enum DistanceMethod {

        ROTAMER, RESIDUE
    };

    /**
     * Checks the distance matrix, finding the shortest distance between two
     * residues' rotamers or any rotamers from two residues under any symmetry
     * operator; will evaluate this if distance matrix not already filled.
     * Default is to find the shortest distance between two residues rather than
     * between two rotamers.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @return Shortest distance
     */
    private double get2BodyDistance(int i, int ri, int j, int rj) {
        if (i > j) {
            int temp = i;
            i = j;
            j = temp;
            temp = ri;
            ri = rj;
            rj = temp;
        }

        if (distanceMethod == DistanceMethod.ROTAMER) {
            return checkDistMatrix(i, ri, j, rj);
        } else {
            return getResidueDistance(i, ri, j, rj);
        }
    }
    
    /** Checks the distance matrix, finding the shortest distance between the closest rotamers of two residues.
     * 
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @return Shortest distance
     */
    private double getResidueDistance(int i, int ri, int j, int rj) {
        if (i > j) {
            int temp = i;
            i = j;
            j = temp;
            temp = ri;
            ri = rj;
            rj = temp;
        }

        double minDist = Double.MAX_VALUE;
        final int lenri = distanceMatrix[i].length;
        final int lenrj = distanceMatrix[i][ri][j].length;
        for (int roti = 0; roti < lenri; roti++) {
            for (int rotj = 0; rotj < lenrj; rotj++) {
                minDist = Math.min(minDist, checkDistMatrix(i, roti, j, rotj));
            }
        }
        return minDist;
    }

    /**
     * Gets the raw distance between two rotamers using lazy loading of the
     * distance matrix. Intended uses: helper method for get2BodyDistance, and
     * for checking for superpositions.
     *
     * @param i A residue index.
     * @param ri A rotamer index for i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for j.
     * @return dist(i, ri, j, rj), ignoring distance matrix method used.
     */
    private double checkDistMatrix(int i, int ri, int j, int rj) {
        if (i > j) {
            int temp = i;
            i = j;
            j = temp;
            temp = ri;
            ri = rj;
            rj = temp;
        }
        double dist = distanceMatrix[i][ri][j][rj];
        if (dist < 0) {
            dist = evaluateDistance(i, ri, j, rj);
            distanceMatrix[i][ri][j][rj] = dist;
        }
        return dist;
    }

    /**
     * Returns the RMS distance between an arbitrary set of rotamers.
     *
     * @param resrot Residue index-rotamer index pairs.
     * @return RMS distance, or Double.MAX_VALUE if ill-defined.
     */
    private double getRawNBodyDistance(int... resrot) {
        int nRes = resrot.length;
        if (nRes % 2 != 0) {
            throw new IllegalArgumentException(" Must have an even number of arguments; res-rot pairs!");
        }
        nRes /= 2;
        if (nRes < 2) {
            throw new IllegalArgumentException(" Must have >= 4 arguments; at least 2 res-rot pairs!");
        }
        // nCr where n is # of residues, choose pairs.
        int numDists = (int) CombinatoricsUtils.binomialCoefficient(nRes, 2);
        double mult = 1.0 / numDists;
        double totDist2 = 0.0;

        for (int i = 0; i < nRes - 1; i++) {
            int i2 = 2 * i;
            for (int j = i + 1; j < nRes; j++) {
                int j2 = 2 * j;
                double rawDist = checkDistMatrix(resrot[i2], resrot[i2 + 1], resrot[j2], resrot[j2 + 1]);
                if (!Double.isFinite(rawDist) || rawDist == Double.MAX_VALUE) {
                    return Double.MAX_VALUE;
                }
                totDist2 += (rawDist * mult);
            }
        }
        return FastMath.sqrt(totDist2);
    }

    /**
     * Returns the RMS separation distance of three 2-body distances. Defaults
     * to Double.MAX_VALUE when there are pair distances outside cutoffs.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @param k Residue k
     * @param rk Rotamer for k
     * @return RMS separation distance
     */
    private double get3BodyDistance(int i, int ri, int j, int rj, int k, int rk) {
        double ij = get2BodyDistance(i, ri, j, rj);
        double ik = get2BodyDistance(i, ri, k, rk);
        double jk = get2BodyDistance(j, rj, k, rk);
        if (areFiniteAndNotMax(ij, ik, jk)) {
            return sqrt((ij * ij + ik * ik + jk * jk) / 3.0);
        } else {
            return Double.MAX_VALUE;
        }
    }

    /**
     * Returns the RMS separation distance for the closest rotamers of three residues' 2-body distances. Defaults
     * to Double.MAX_VALUE when there are pair distances outside cutoffs.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @param k Residue k
     * @param rk Rotamer for k
     * @return RMS separation distance
     */
    private double get3BodyResidueDistance(int i, int ri, int j, int rj, int k, int rk) {
        double ij = getResidueDistance(i, ri, j, rj);
        double ik = getResidueDistance(i, ri, k, rk);
        double jk = getResidueDistance(j, rj, k, rk);
        if (areFiniteAndNotMax(ij, ik, jk)) {
            return sqrt((ij * ij + ik * ik + jk * jk) / 3.0);
        } else {
            return Double.MAX_VALUE;
        }
    }
    
    /**
     * Returns the RMS separation distance of 6 2-body distances. Defaults to
     * Double.MAX_VALUE when there are pair distances outside cutoffs.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @param k Residue k
     * @param rk Rotamer for k
     * @param l Residue l
     * @param rl Rotamer for l
     * @return RMS separation distance
     */
    private double get4BodyDistance(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        double ij = get2BodyDistance(i, ri, j, rj);
        double ik = get2BodyDistance(i, ri, k, rk);
        double il = get2BodyDistance(i, ri, l, rl);
        double jk = get2BodyDistance(j, rj, k, rk);
        double jl = get2BodyDistance(j, rj, l, rl);
        double kl = get2BodyDistance(k, rk, l, rl);
        if (areFiniteAndNotMax(ij, ik, il, jk, jl, kl)) {
            return sqrt((ij * ij + ik * ik + il * il + jk * jk + jl * jl + kl * kl) / 6.0);
        } else {
            return Double.MAX_VALUE;
        }
    }
    
    /**
     * Returns the RMS separation distance for the closest rotamers of 6 2-body distances from four residues. Defaults to
     * Double.MAX_VALUE when there are pair distances outside cutoffs.
     *
     * @param i Residue i
     * @param ri Rotamer for i
     * @param j Residue j
     * @param rj Rotamer for j
     * @param k Residue k
     * @param rk Rotamer for k
     * @param l Residue l
     * @param rl Rotamer for l
     * @return RMS separation distance
     */
    private double get4BodyResidueDistance(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        double ij = getResidueDistance(i, ri, j, rj);
        double ik = getResidueDistance(i, ri, k, rk);
        double il = getResidueDistance(i, ri, l, rl);
        double jk = getResidueDistance(j, rj, k, rk);
        double jl = getResidueDistance(j, rj, l, rl);
        double kl = getResidueDistance(k, rk, l, rl);
        if (areFiniteAndNotMax(ij, ik, il, jk, jl, kl)) {
            return sqrt((ij * ij + ik * ik + il * il + jk * jk + jl * jl + kl * kl) / 6.0);
        } else {
            return Double.MAX_VALUE;
        }
    }

    /**
     * Returns true if all values (usually distances) are both finite and not
     * set to Double.MAX_VALUE.
     *
     * @param values Values to check.
     * @return If all values are regular, finite values.
     */
    private static boolean areFiniteAndNotMax(double... values) {
        return Arrays.stream(values).
                allMatch((double val) -> Double.isFinite(val) && val < Double.MAX_VALUE);
    }

    /**
     * Checks if the i,ri,j,rj pair exceeds the pair distance thresholds.
     *
     * @param i A residue index.
     * @param ri A rotamer index for i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for j.
     * @return If i,ri,j,rj > threshold distance.
     */
    private boolean checkPairDistThreshold(int i, int ri, int j, int rj) {
        if (twoBodyCutoffDist <= 0 || !Double.isFinite(twoBodyCutoffDist)) {
            return false;
        }
        return (get2BodyDistance(i, ri, j, rj) > twoBodyCutoffDist);
    }

    /**
     * Checks if the i,ri,j,rj,k,rk triple exceeds the 3-body threshold, or if
     * any component exceeds the pair distance threshold.
     *
     * @param i A residue index.
     * @param ri A rotamer index for i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for j.
     * @param k A residue index k!=i, k!=j.
     * @param rk A rotamer index for k.
     * @return If i,ri,j,rj,k,rk > threshold distances.
     */
    private boolean checkTriDistThreshold(int i, int ri, int j, int rj, int k, int rk) {
        if (checkPairDistThreshold(i, ri, j, rj) || checkPairDistThreshold(i, ri, k, rk)
                || checkPairDistThreshold(j, rj, k, rk)) {
            return true;
        }
        if (threeBodyCutoffDist <= 0 || !Double.isFinite(threeBodyCutoffDist)) {
            return false;
        }
        return (get3BodyDistance(i, ri, j, rj, k, rk) > threeBodyCutoffDist);
    }

    /**
     * Checks if the i,ri,j,rj,k,rk,l,rl quad exceeds the 3-body threshold, or
     * if any component exceeds the pair/triple distance thresholds.
     *
     * @param i A residue index.
     * @param ri A rotamer index for i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for j.
     * @param k A residue index k!=i, k!=j.
     * @param rk A rotamer index for k.
     * @param l A residue index l!=i, l!=j, l!=k.
     * @param rl A rotamer index for l.
     * @return If i,ri,j,rj,k,rk,l,rl > threshold distances.
     */
    private boolean checkQuadDistThreshold(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        if (checkTriDistThreshold(i, ri, j, rj, k, rk) || checkTriDistThreshold(i, ri, j, rj, l, rl)
                || checkTriDistThreshold(i, ri, k, rk, l, rl) || checkTriDistThreshold(j, rj, k, rk, l, rl)) {
            return true;
        }
        // Use the 3-body cutoff distance for now.
        if (threeBodyCutoffDist <= 0 || !Double.isFinite(threeBodyCutoffDist)) {
            return false;
        }
        return get4BodyDistance(i, ri, j, rj, k, rk, l, rl) > threeBodyCutoffDist;
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
            throw new IllegalStateException(" CheckForcedResidues being called without useForcedResidues.");
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
                                        double rotamerSeparation = get2BodyDistance(indexI, ri, j, rj);
                                        // if (distanceMatrix[indexI][ri][j][rj] <= distance) {
                                        if (rotamerSeparation <= distance) {
                                            if (!currentWindow.contains(residuej)) {
                                                logIfMaster(format(" Adding residue %s at distance %16.8f Ang from %s %d.",
                                                        residuej.toFormattedString(false, true), rotamerSeparation, residuei.toFormattedString(false, true), ri));
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
                        double startingEnergy = Double.NaN;
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
                                double finalEnergy = Double.NaN;
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
                            double finalEnergy = Double.NaN;
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
        for (int i = 0; i < numCells; i++) {
            BoxOptCell celli = cells[i];
            ArrayList<Residue> residuesList = celli.getResiduesAsList();
            int[] cellIndices = celli.getXYZIndex();
            logIfMaster(format("\n Iteration %d of the box optimization.", (i + 1)));
            logIfMaster(format(" Cell index (linear): %d", (celli.getLinearIndex() + 1)));
            logIfMaster(format(" Cell xyz indices: x = %d, y = %d, z = %d", cellIndices[0] + 1, cellIndices[1] + 1, cellIndices[2] + 1));
            int nResidues = residuesList.size();
            if (nResidues > 0) {
                readyForSingles = false;
                finishedSelf = false;
                readyFor2Body = false;
                finished2Body = false;
                readyFor3Body = false;
                finished3Body = false;
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
                    partialFile = new File(filename);
                    if (firstCellSaved) {
                        partialFile.delete();
                    }
                    // Don't write a file if it's the final iteration.
                    if (i == (numCells - 1)) {
                        continue;
                    }
                    PDBFilter windowFilter = new PDBFilter(partialFile, molecularAssembly, null, null);
                    try {
                        windowFilter.writeFile(partialFile, false);
                        if (firstResidue != lastResidue) {
                            logIfMaster(format(" File with residues %s ... %s in window written.", firstResidue.toString(), lastResidue.toString()));
                        } else {
                            logIfMaster(format(" File with residue %s in window written.", firstResidue.toString()));
                        }
                        firstCellSaved = true;
                    } catch (Exception e) {
                        logger.warning(format("Exception writing to file: %s", partialFile.getName()));
                    }
                }
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
     * @param crystal Polymer crystal or dummy crystal
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
            cell.sortBoxResidues();
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
        int nSymm = crystal.spaceGroup.getNumberOfSymOps();

        for (BoxOptCell cell : cells) {
            Set<Residue> toAdd = new HashSet<>();
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                for (Residue residue : residues) {
                    boolean contained;
                    switch (boxInclusionCriterion) {
                        case 2: {
                            contained = cell.residueInsideCell(residue, crystal, symOp, true);
                        }
                        break;
                        case 3: {
                            contained = cell.anyRotamerInsideCell(residue, crystal, symOp, true);
                        }
                        break;
                        case 1:
                        default: {
                            contained = cell.atomInsideCell(residue.getReferenceAtom(), crystal, symOp);
                        }
                    }
                    if (contained) {
                        toAdd.add(residue);
                    }
                }
                // If the identity symop produces nothing, skip checking other symops.
                if (toAdd.isEmpty()) {
                    break;
                }
            }
            toAdd.forEach(cell::addResidue);
        }
    }

    /**
     * Sorts a passed List of Residues by global index.
     *
     * @param residues List of Residues to be sorted.
     */
    @SuppressWarnings("unchecked")
    private void sortResidues(List<Residue> residues) {
        Comparator comparator = Comparator.comparing(Residue::getChainID).thenComparingInt(Residue::getResidueNumber);
        residues.sort(comparator);
    }

    /**
     * ONLY FOR UNIT TESTING. Sets a boolean to turn the self elimination
     * criteria off.
     */
    public void turnRotamerSingleEliminationOff() {
        logger.warning("TURNING SINGLE ELIMINATION OFF.");
        selfEliminationOn = false;
    }

    /**
     * ONLY FOR UNIT TESTING. Sets a boolean to turn the pair elimination
     * criteria off.
     */
    public void turnRotamerPairEliminationOff() {
        logger.warning("TURNING PAIR ELIMINATION OFF.");
        pairEliminationOn = false;
    }

    protected void applyEliminationCriteria(Residue residues[], boolean getEnergies, boolean verbose) {
        Level prevLevel = logger.getLevel();
        if (!verbose) {
            logger.setLevel(Level.WARNING);
        }
        if (getEnergies) {
            applyEliminationCriteria(residues);
        } else {
            //allocateEliminationMemory is now called for all algorithms in rotamerEnergies method.
            //allocateEliminationMemory(residues);
            int i = 0;
            boolean pairEliminated;
            do {
                pairEliminated = false;
                if (useGoldstein) {
                    if (selfEliminationOn) {
                        i++;
                        logIfMaster(format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
                        // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                        while (goldsteinDriver(residues)) {
                            i++;
                            logIfMaster(this.toString());
                            logIfMaster(format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
                        }
                    }
                    if (pairEliminationOn) {
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
                    }
                } else {
                    if (selfEliminationOn) {
                        i++;
                        logIfMaster(format("\n Iteration %d: Applying Single DEE conditions ", i));
                        // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                        while (deeRotamerElimination(residues)) {
                            i++;
                            logIfMaster(toString());
                            logIfMaster(format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
                        }
                    }
                    if (pairEliminationOn) {
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
                }
                validateDEE(residues);
                logIfMaster(toString());
            } while (pairEliminated);
            logIfMaster(" Self-consistent DEE rotamer elimination achieved.\n");
        }
        if (!verbose) {
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
        // allocateEliminationMemory is now called for all algorithms in rotamerEnergies method.
        //allocateEliminationMemory(residues);

        if (verboseEnergies) {
            try {
                logIfMaster(format("\n Beginning Energy %s", formatEnergy(currentEnergy(residues))));
            } catch (ArithmeticException ex) {
                logger.severe(String.format(" Exception %s in calculating beginning energy; FFX shutting down.", ex.toString()));
            }
        }

        rotamerEnergies(residues);

        if (testing) {
            int nres = residues.length;
            onlyPrunedSingles = new boolean[nres][];
            onlyPrunedPairs = new boolean[nres][][][];
            for (int i = 0; i < nres; i++) {
                Residue residuei = residues[i];
                Rotamer rotamersi[] = residuei.getRotamers(library);
                int lenri = rotamersi.length;  // Length rotamers i
                onlyPrunedSingles[i] = new boolean[lenri];
                onlyPrunedSingles[i] = Arrays.copyOf(eliminatedSingles[i], eliminatedSingles[i].length);
                onlyPrunedPairs[i] = new boolean[lenri][][];
                // Loop over the set of rotamers for residue i.
                for (int ri = 0; ri < lenri; ri++) {
                    onlyPrunedPairs[i][ri] = new boolean[nres][];
                    for (int j = i + 1; j < nres; j++) {
                        Residue residuej = residues[j];
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        int lenrj = rotamersj.length;
                        onlyPrunedPairs[i][ri][j] = new boolean[lenrj];
                        onlyPrunedPairs[i][ri][j] = Arrays.copyOf(eliminatedPairs[i][ri][j], eliminatedPairs[i][ri][j].length);
                    }
                }
            }
        }

        if (testSelfEnergyEliminations) {
            testSelfEnergyElimination(residues);
        } else if (testPairEnergyEliminations > -1) {
            testPairEnergyElimination(residues, testPairEnergyEliminations);
        } else if (testTripleEnergyEliminations1 > -1 && testTripleEnergyEliminations2 > -1) {
            testTripleEnergyElimination(residues, testTripleEnergyEliminations1, testTripleEnergyEliminations2);
        }

        // testSelfEnergyElimination(residues);
        // testPairEnergyElimination(residues, 19);
        // Beginning energy
        int currentRotamers[] = new int[residues.length];

        if (pruneClashes) {
            validateDEE(residues);
        }

        int i = 0;
        boolean pairEliminated;
        do {
            pairEliminated = false;
            if (useGoldstein) {
                if (selfEliminationOn) {
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
                    // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                    while (goldsteinDriver(residues)) {
                        i++;
                        logIfMaster(this.toString());
                        logIfMaster(format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
                    }
                }
                if (pairEliminationOn) {
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                    // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                    while (goldsteinPairDriver(residues)) {
                        pairEliminated = true;
                        i++;
                        logIfMaster(this.toString());
                        logIfMaster(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
                    }
                }
            } else {
                if (selfEliminationOn) {
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Single DEE conditions ", i));
                    // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
                    while (deeRotamerElimination(residues)) {
                        i++;
                        logIfMaster(toString());
                        logIfMaster(format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
                    }
                }
                if (pairEliminationOn) {
                    i++;
                    logIfMaster(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                    // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
                    while (deeRotamerPairElimination(residues)) {
                        pairEliminated = true;
                        i++;
                        logIfMaster(toString());
                        logIfMaster(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
                    }
                }
            }
            validateDEE(residues);
            logIfMaster(toString());
        } while (pairEliminated);

        logIfMaster(" Self-consistent DEE rotamer elimination achieved.\n");
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
                filter(res -> res != null).
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
        potential.getCoordinates(x);
        return potential.energy(x);
    }

    /**
     * Wrapper intended for use with RotamerMatrixMC.
     *
     * @param resList
     * @return
     * @throws ArithmeticException
     */
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
        allocateEliminationMemory(residues);

        if (decomposeOriginal) {
            assert library.getUsingOrigCoordsRotamer();
            for (int i = 0; i < nResidues; i++) {
                Residue resi = residues[i];
                Rotamer[] rotsi = resi.getRotamers(library);
                int lenri = rotsi.length;

                // Leave the 0'th original-coordinates rotamer alone.
                for (int ri = 1; ri < lenri; ri++) {
                    eliminateRotamer(residues, i, ri, false);
                }
            }
        }

        /**
         * Initialize all atoms to be used.
         */
        for (int i = 0; i < nAtoms; i++) {
            atoms[i].setUse(true);
        }

        if (!usingBoxOptimization) {
            energiesToWrite = Collections.synchronizedList(new ArrayList<>());
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
        SelfEnergyRegion singlesRegion = new SelfEnergyRegion(residues);
        TwoBodyEnergyRegion pairsRegion = new TwoBodyEnergyRegion(residues);
        ThreeBodyEnergyRegion triplesRegion = new ThreeBodyEnergyRegion(residues);
        FourBodyEnergyRegion quadsRegion = new FourBodyEnergyRegion(residues);

        try {
            if (loaded < 1) {
                selfEnergyMap.clear();
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
                            selfEnergyMap.put(singleJobIndex++, selfJob);
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

            logger.info(String.format(" Number of self energies to calculate: %d", selfEnergyMap.size()));

            energyWorkerTeam.execute(singlesRegion);
            long singlesTime = System.nanoTime() - energyStartTime;
            logIfMaster(format(" Time for single energies: %12.4g", (singlesTime * 1.0E-9)));

            if (loaded < 2) {
                twoBodyEnergyMap.clear();
                // allocate twoBodyEnergy and create jobs
                int pairJobIndex = 0;
                twoBodyEnergy = new double[nResidues][][][];
                for (int i = 0; i < nResidues; i++) {
                    Residue resi = residues[i];
                    int indexI = allResiduesList.indexOf(resi);
                    Rotamer roti[] = resi.getRotamers(library);
                    twoBodyEnergy[i] = new double[roti.length][][];
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (check(i, ri)) {
                            continue;
                        }
                        twoBodyEnergy[i][ri] = new double[nResidues][];
                        for (int j = i + 1; j < nResidues; j++) {
                            Residue resj = residues[j];
                            int indexJ = allResiduesList.indexOf(resj);
                            Rotamer rotj[] = resj.getRotamers(library);
                            twoBodyEnergy[i][ri][j] = new double[rotj.length];
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if (checkToJ(i, ri, j, rj)) {
                                    continue;
                                }

                                // Skip creating a job if the pair is outside pair cut-off.
                                if (checkPairDistThreshold(indexI, ri, indexJ, rj)) {
                                    continue;
                                }

                                Integer pairJob[] = {i, ri, j, rj};
                                if (decomposeOriginal && (ri != 0 || rj != 0)) {
                                    continue;
                                }
                                twoBodyEnergyMap.put(pairJobIndex++, pairJob);
                            }
                        }
                    }
                }
            }

            // broadcast that this proc is done with pruning and allocation; ready for pairs
            multicastBuf(thisProcReadyBuf);
            // launch parallel pair calculation
            while (!readyFor2Body) {
                Thread.sleep(POLLING_FREQUENCY);
            }

            logger.info(String.format(" Number of 2-body energies to calculate: %d", twoBodyEnergyMap.size()));
            energyWorkerTeam.execute(pairsRegion);

            long pairsTime = System.nanoTime() - (singlesTime + energyStartTime);
            long triplesTime = 0;
            long quadsTime = 0;
            logIfMaster(format(" Time for 2-body energies:   %12.4g", (pairsTime * 1.0E-9)));

            if (threeBodyTerm) {
                if (loaded < 3) {
                    alloc3BodyMap(residues);
                    threeBodyEnergyMap.clear();
                    // allocate threeBodyEnergy and create jobs
                    int trimerJobIndex = 0;
                    for (int i = 0; i < nResidues; i++) {
                        Residue resi = residues[i];
                        int indexI = allResiduesList.indexOf(resi);
                        Rotamer roti[] = resi.getRotamers(library);
                        for (int ri = 0; ri < roti.length; ri++) {
                            if (check(i, ri)) {
                                continue;
                            }
                            for (int j = i + 1; j < nResidues; j++) {
                                Residue resj = residues[j];
                                int indexJ = allResiduesList.indexOf(resj);
                                Rotamer rotj[] = resj.getRotamers(library);
                                for (int rj = 0; rj < rotj.length; rj++) {
                                    if (checkToJ(i, ri, j, rj)) {
                                        continue;
                                    }
                                    for (int k = j + 1; k < nResidues; k++) {
                                        Residue resk = residues[k];
                                        int indexK = allResiduesList.indexOf(resk);
                                        Rotamer rotk[] = resk.getRotamers(library);
                                        for (int rk = 0; rk < rotk.length; rk++) {
                                            if (checkToK(i, ri, j, rj, k, rk)) {
                                                continue;
                                            }

                                            // Skip work items outside the 3-body cutoff.
                                            if (checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
                                                continue;
                                            }

                                            Integer trimerJob[] = {i, ri, j, rj, k, rk};
                                            if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0)) {
                                                continue;
                                            }
                                            threeBodyEnergyMap.put(trimerJobIndex++, trimerJob);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // broadcast that this proc is done with pruning and allocation; ready for trimers
                multicastBuf(thisProcReadyBuf);
                // launch parallel 3-Body calculation
                while (!readyFor3Body) {
                    Thread.sleep(POLLING_FREQUENCY);
                }

                logger.info(String.format(" Number of 3-Body energies to calculate: %d", threeBodyEnergyMap.size()));
                energyWorkerTeam.execute(triplesRegion);
                triplesTime = System.nanoTime() - (pairsTime + singlesTime + energyStartTime);
                logIfMaster(format(" Time for 3-Body energies: %12.4g", (triplesTime * 1.0E-9)));
            }

            if (compute4BodyEnergy) {
                logger.info(" Creating 4-Body jobs...");
                fourBodyEnergyMap.clear();
                boolean maxedOut = false;
                // create 4-Body jobs (no memory allocation)
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
                                                fourBodyEnergyMap.put(quadJobIndex++, quadJob);
                                                if (fourBodyEnergyMap.size() > max4BodyCount) {
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
                // logger.info(format(" Proc %d broadcasting ready for quads.", world.rank()));
                multicastBuf(thisProcReadyBuf);
                // launch parallel 3-Body calculation
                while (!readyFor4Body) {
                    Thread.sleep(POLLING_FREQUENCY);
                }

                logger.info(format(" Running quads: %d jobs.", fourBodyEnergyMap.size()));
                energyWorkerTeam.execute(quadsRegion);
                quadsTime = System.nanoTime() - (triplesTime + pairsTime + singlesTime + energyStartTime);
                logIfMaster(format(" Time for 4-Body energies:   %12.4g", quadsTime * 1.0E-9));
            }
            long allTime = singlesTime + pairsTime + triplesTime + quadsTime;
            logIfMaster(format(" Time for all energies:    %12.4g", allTime * 1.0E-9));
        } catch (Exception ex) {
            String message = " Exception computing rotamer energies in parallel.";
            logger.log(Level.SEVERE, message, ex);
        }

        // Turn on all atoms.
        for (int i = 0; i < atoms.length; i++) {
            atoms[i].setUse(true);
        }
        // Print the energy with all rotamers in their default conformation.
        if (verboseEnergies && master) {
            try {
                double defaultEnergy = currentEnergy(residues);
                logger.info(format(" Energy of the system with rotamers in their default conformation: %s",
                        formatEnergy(defaultEnergy)));
            } catch (ArithmeticException ex) {
                logger.severe(String.format(" Exception %s in calculating default energy; FFX shutting down", ex.toString()));
            }
        }
        return backboneEnergy;
    }

    /**
     * Computes the environment/backbone energy, defined as energy with all
     * sidechains under consideration turned off in their 0th rotamer.
     *
     * @param residues Residues under optimization.
     * @return Backbone energy Eenv/bb.
     * @throws ArithmeticException If an exception in calculating energy is
     * found.
     */
    private double computeBackboneEnergy(Residue[] residues) throws ArithmeticException {
        // Set all atoms to be used.
        Atom atoms[] = molecularAssembly.getAtomArray();
        for (Atom atom : atoms) {
            atom.setUse(true);
        }

        // Turn off all Residues and set them to their default conformation.
        turnOffAllResidues(residues);

        // Compute and return the backbone energy.
        return currentEnergy(residues);
    }

    /**
     * Applies the "default" rotamer: currently the 0'th rotamer.
     *
     * @param residue Residue to apply a default rotamer for.
     */
    private void applyDefaultRotamer(Residue residue) {
        RotamerLibrary.applyRotamer(residue, residue.getRotamers(library)[0]);
    }

    /**
     * Computes a self energy, defined as energy with all sidechains but one
     * turned off, minus the backbone energy.
     *
     * @param residues Residues under optimization.
     * @param i A residue index.
     * @param ri A rotamer index for residue i.
     * @return Eself(ri)=E1(ri)-Eenv/bb.
     */
    private double computeSelfEnergy(Residue[] residues, int i, int ri) {
        turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            energy = currentEnergy(residues) - backboneEnergy;
        } finally {
            turnOffResidue(residues[i]);
        }

        return energy;
    }

    /**
     * Computes a pair energy, defined as energy with all sidechains but two
     * turned off, minus the sum of backbone and component self energies.
     *
     * @param residues Residues under optimization.
     * @param i A residue index.
     * @param ri A rotamer index for residue i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for residue j.
     * @return Epair(ri, rj)=E2(ri,rj)-Eself(ri)-Eself(rj)-Eenv/bb.
     */
    private double compute2BodyEnergy(Residue[] residues, int i, int ri, int j, int rj) {
        turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        turnOnResidue(residues[j], rj);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            energy = currentEnergy(residues) - backboneEnergy - getSelf(i, ri) - getSelf(j, rj);
        } finally {
            // Revert if the currentEnergy call throws an exception.
            turnOffResidue(residues[i]);
            turnOffResidue(residues[j]);
        }
        return energy;
    }

    /**
     * Computes a 3-body energy, defined as the energy with all sidechains but
     * three turned off, minus the sum of backbone and component self/2-Body
     * energies.
     *
     * @param residues Residues under optimization.
     * @param i A residue index.
     * @param ri A rotamer index for residue i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for residue j.
     * @param k A residue index k!=j k!=i.
     * @param rk A rotamer index for residue k.
     * @return Etri(ri,
     * rj)=E3(ri,rj,rk)-Epair(ri,rj)-Epair(ri,rk)-Epair(rj,rk)-Eself(ri)-Eself(rj)-Eself(rk)-Eenv/bb.
     */
    private double compute3BodyEnergy(Residue[] residues, int i, int ri, int j, int rj, int k, int rk) {
        turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        turnOnResidue(residues[j], rj);
        turnOnResidue(residues[k], rk);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            energy = currentEnergy(residues) - backboneEnergy - getSelf(i, ri) - getSelf(j, rj) - getSelf(k, rk)
                    - get2Body(i, ri, j, rj) - get2Body(i, ri, k, rk) - get2Body(j, rj, k, rk);
        } finally {
            // Revert if the currentEnergy call throws an exception.
            turnOffResidue(residues[i]);
            turnOffResidue(residues[j]);
            turnOffResidue(residues[k]);
        }
        return energy;
    }

    /**
     * Computes a 4-body energy, defined as the energy with all sidechains but
     * four turned off, minus the sum of backbone and component
     * self/2-Body/3-body energies.
     *
     * @param residues Residues under optimization.
     * @param i A residue index.
     * @param ri A rotamer index for residue i.
     * @param j A residue index j!=i.
     * @param rj A rotamer index for residue j.
     * @param k A residue index k!=j k!=i.
     * @param rk A rotamer index for residue k.
     * @param l A residue index l!=i l!=j l!=k.
     * @param rl A rotamer index for residue l.
     * @return The 4-body energy.
     */
    private double compute4BodyEnergy(Residue[] residues, int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        turnOffAllResidues(residues);
        turnOnResidue(residues[i], ri);
        turnOnResidue(residues[j], rj);
        turnOnResidue(residues[k], rk);
        turnOnResidue(residues[l], rl);
        double energy;
        try {
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            energy = currentEnergy(residues) - backboneEnergy
                    - getSelf(i, ri) - getSelf(j, rj) - getSelf(k, rk) - getSelf(l, rl)
                    - get2Body(i, ri, j, rj) - get2Body(i, ri, k, rk) - get2Body(i, ri, l, rl)
                    - get2Body(j, rj, k, rk) - get2Body(j, rj, l, rl) - get2Body(k, rk, l, rl)
                    - get3Body(residues, i, ri, j, rj, k, rk) - get3Body(residues, i, ri, j, rj, l, rl) - get3Body(residues, i, ri, k, rk, l, rl) - get3Body(residues, j, rj, k, rk, l, rl);

        } finally {
            // Revert if the currentEnergy call throws an exception.
            turnOffResidue(residues[i]);
            turnOffResidue(residues[j]);
            turnOffResidue(residues[k]);
            turnOffResidue(residues[l]);
        }
        return energy;
    }

    private void turnOnResidue(Residue residue, int ri) {
        turnOnAtoms(residue);
        Rotamer rotamers[] = residue.getRotamers(library);
        RotamerLibrary.applyRotamer(residue, rotamers[ri]);
    }

    private void turnOffResidue(Residue residue) {
        turnOffAtoms(residue);
        applyDefaultRotamer(residue);
    }

    private void turnOffAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (int i = 0; i < nRes; i++) {
            turnOffResidue(residues[i]);
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
     * Calculates a residue-residue distance matrix.
     * <p>
     * Residue-residue distance is defined as the shortest atom-atom distance in
     * any possible rotamer-rotamer pair if the residues are neighbors (central
     * atom-central atom distances are within a cutoff). Otherewise, distances
     * are set to a default of Double.MAX_VALUE.
     * </p>
     * <p>
     * The intent of using a neighbor list is to avoid tediously searching
     * rotamer- rotamer pairs when two residues are so far apart we will never
     * need the exact distance. We use the distance matrix for adding residues
     * to the sliding window and determining whether to set 3-body energy to
     * 0.0.
     * </p>
     * <p>
     * If the central atoms are too distant from each other, we can safely
     * assume no atoms will ever be close enough for addition to sliding window
     * or to cause explicit calculation of 3-body energy.
     * </p>
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
                            residuei.toFormattedString(false, true)), Level.WARNING);
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
                            logIfMaster(format(" Residue j %s has null rotamers.", residuej.toFormattedString(false, true)));
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

            double nlistCutoff = Math.max(Math.max(distance, twoBodyCutoffDist), threeBodyCutoffDist);
            /**
             * I think this originated from the fact that side-chain (and later
             * nucleic acid) atoms could be fairly distant from the reference
             * atom.
             */
            double magicNumberBufferOfUnknownOrigin = 25.0;
            nlistCutoff += magicNumberBufferOfUnknownOrigin;

            NeighborList neighborList = new NeighborList(null, crystal,
                    atoms, nlistCutoff, 0.0, parallelTeam);

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
        return sqrt(dist);
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
            logIfMaster(format(" %3d Residue %7s with %2d rotamers.", i + 1, residuei.toFormattedString(false, true), lenri));
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
     * Elimination of rotamers via the original Dead End Elimination algorithm.
     *
     * @param residues Array of residues under consideration.
     * @return True if any rotamers were eliminated.
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
                minEnergySingles[ri] = getSelf(i, ri);
                maxEnergySingles[ri] = minEnergySingles[ri];
                for (int j = 0; j < nres; j++) {
                    if (j == i) {
                        continue;
                    }
                    if (minMaxPairEnergy(residues, minMax, i, ri, j)) {
                        if (Double.isFinite(minMax[0]) && Double.isFinite(minEnergySingles[ri])) {
                            minEnergySingles[ri] += minMax[0];
                        } else {
                            // There is a possibility that ri conflicts with every possible rotamer of some residue j.
                            // In that event, its minimum energy is set NaN and should be easy to eliminate.
                            minEnergySingles[ri] = Double.NaN;
                        }
                        if (Double.isFinite(minMax[0]) && Double.isFinite(maxEnergySingles[ri])) {
                            maxEnergySingles[ri] += minMax[1];
                        } else {
                            // In this branch, ri conflicts with some j,rj and cannot be practically used for elimination.
                            maxEnergySingles[ri] = Double.NaN;
                        }
                    } else {
                        Residue residuej = residues[j];
                        logger.info(format(" Inconsistent Pair: %8s %2d, %8s.", residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true)));
                        //eliminateRotamer(residues, i, ri, print);
                    }
                }
            }

            /**
             * Apply the singles elimination criteria to rotamers of residue i
             * by determining the most favorable maximum energy.
             */
            double eliminationEnergy = Double.MAX_VALUE;
            int eliminatingRotamer = 0;
            for (int ri = 0; ri < lenri; ri++) {
                if (check(i, ri)) {
                    continue;
                }
                if (Double.isFinite(maxEnergySingles[ri]) && maxEnergySingles[ri] < eliminationEnergy) {
                    eliminationEnergy = maxEnergySingles[ri];
                    eliminatingRotamer = ri;
                }
            }

            if (eliminationEnergy == Double.MAX_VALUE) {
                // This branch is taken if every ri conflicts with at least one j,rj. In that case, nothing can be eliminated yet!
                logIfMaster(" Could not eliminate any i,ri because eliminationEnergy was never set!", Level.FINE);
            } else {
                /**
                 * Eliminate rotamers whose minimum energy is greater than the
                 * worst case for another rotamer.
                 */
                for (int ri = 0; ri < lenri; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    // If i,ri has a clash with all of phase space, it can be eliminated by something that doesn't clash with all phase space.
                    if (!Double.isFinite(minEnergySingles[ri])) {
                        if (eliminateRotamer(residues, i, ri, print)) {
                            logIfMaster(format("  Rotamer elimination of (%8s,%2d) that always clashes.", residuei.toFormattedString(false, true), ri));
                            eliminated = true;
                        }
                    }
                    // Otherwise, can eliminate if its best possible energy is still worse than something else's worst possible energy.
                    if (minEnergySingles[ri] > eliminationEnergy + ensembleBuffer) {
                        if (eliminateRotamer(residues, i, ri, print)) {
                            logIfMaster(format("  Rotamer elimination of (%8s,%2d) by (%8s,%2d): %12.4f > %6.4f.",
                                    residuei.toFormattedString(false, true), ri, residuei.toFormattedString(false, true), eliminatingRotamer, minEnergySingles[ri], eliminationEnergy + ensembleBuffer));
                            eliminated = true;
                        }
                    }
                }
            }
        }
        return eliminated;
    }

    /**
     * Rotamer pair elimination driver for many-body Dead End Elimination.
     * Generally less effective than Goldstein.
     *
     * @param residues Residues under consideration.
     * @return If at least one pair eliminated.
     */
    private boolean deeRotamerPairElimination(Residue[] residues) {
        int nres = residues.length;
        boolean eliminated = false;

        for (int i = 0; i < (nres - 1); i++) {
            Residue residuei = residues[i];
            Rotamer rotamersi[] = residuei.getRotamers(library);
            int lenri = rotamersi.length;

            // Minimum and maximum summation found for ri-rj pairs.
            double[][] minPairEnergies = new double[lenri][];
            double[][] maxPairEnergies = new double[lenri][];

            for (int j = i + 1; j < nres; j++) {
                Residue residuej = residues[j];
                Rotamer[] rotamersj = residuej.getRotamers(library);
                int lenrj = rotamersj.length;

                for (int ri = 0; ri < lenri; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    minPairEnergies[ri] = new double[lenrj];
                    maxPairEnergies[ri] = new double[lenrj];

                    for (int rj = 0; rj < lenrj; rj++) {
                        if (check(j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        minPairEnergies[ri][rj] = getSelf(i, ri) + getSelf(j, rj) + get2Body(i, ri, j, rj);
                        maxPairEnergies[ri][rj] = minPairEnergies[ri][rj];

                        // Min and max external summations for ri-rj.
                        double[] minMax = new double[2];

                        // Add contributions from third residues k, and possibly fourth residues l.
                        if (minMaxE2(residues, minMax, i, ri, j, rj)) {
                            if (Double.isFinite(minPairEnergies[ri][rj]) && Double.isFinite(minMax[0])) {
                                minPairEnergies[ri][rj] += minMax[0];
                            } else {
                                logger.severe(String.format(" An ri-rj pair %s-%d %s-%d with NaN minimum was caught incorrectly!", residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj));
                            }
                            if (Double.isFinite(maxPairEnergies[ri][rj]) && Double.isFinite(minMax[1])) {
                                maxPairEnergies[ri][rj] += minMax[1];
                            } else {
                                // ri-rj can clash, and isn't very useful to eliminate by.
                                maxPairEnergies[ri][rj] = Double.NaN;
                            }
                        } else {
                            // A NaN minimum energy for some pair indicates it's definitely not part of the GMEC.
                            minPairEnergies[ri][rj] = Double.NaN;
                            logger.info(String.format(" Eliminating pair %s-%d %s-%d that always clashes.", residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj));
                            eliminateRotamerPair(residues, i, ri, j, rj, print);
                            eliminated = true;
                        }
                    }
                }

                double pairEliminationEnergy = Double.MAX_VALUE;
                for (int ri = 0; ri < lenri; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (check(j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        if (Double.isFinite(maxPairEnergies[ri][rj]) && maxPairEnergies[ri][rj] < pairEliminationEnergy) {
                            pairEliminationEnergy = maxPairEnergies[ri][rj];
                        }
                    }
                }

                if (pairEliminationEnergy == Double.MAX_VALUE) {
                    logIfMaster(String.format(" All rotamer pairs for residues %s and %s have possible conflicts; cannot perform any eliminations!", residuei.toFormattedString(false, true), residuej), Level.FINE);
                } else {
                    double comparisonEnergy = pairEliminationEnergy + ensembleBuffer;
                    for (int ri = 0; ri < lenri; ri++) {
                        if (check(i, ri)) {
                            continue;
                        }
                        for (int rj = 0; rj < lenrj; rj++) {
                            if (check(j, rj) || check(i, ri, j, rj)) {
                                continue;
                            }
                            if (minPairEnergies[ri][rj] > comparisonEnergy) {
                                if (eliminateRotamerPair(residues, i, ri, j, rj, print)) {
                                    eliminated = true;
                                    logIfMaster(format(" Eliminating rotamer pair: %s %d, %s %d (%s > %s + %6.6f)",
                                            residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj,
                                            formatEnergy(minPairEnergies[ri][rj]),
                                            formatEnergy(pairEliminationEnergy), ensembleBuffer), Level.INFO);
                                } else {
                                    // See above check(i, ri, j, rj) for why this should not be taken!
                                    logIfMaster(format(" Already eliminated rotamer pair! %s %d, %s %d (%s > %1s + %6.6f)",
                                            residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj,
                                            formatEnergy(minPairEnergies[ri][rj]),
                                            formatEnergy(pairEliminationEnergy), ensembleBuffer), Level.WARNING);
                                }
                            }
                        }
                    }
                }

                if (pairsToSingleElimination(residues, i, j)) {
                    eliminated = true;
                }
            }
        }

        return eliminated;
    }

    /**
     * Method to check if pairs elimination for some residue pair has enabled a
     * singles rotamer elimination by eliminating all ri-rj for some ri or some
     * rj.
     *
     * @param residues Residues under consideration.
     * @param i A residue index.
     * @param j A residue index j!=i
     * @return If any singletons were eliminated.
     */
    private boolean pairsToSingleElimination(Residue[] residues, int i, int j) {
        assert i != j;
        assert i < residues.length;
        assert j < residues.length;

        Residue residuei = residues[i];
        Residue residuej = residues[j];
        Rotamer[] rotsi = residuei.getRotamers(library);
        Rotamer[] rotsj = residuej.getRotamers(library);
        int lenri = rotsi.length;
        int lenrj = rotsj.length;
        boolean eliminated = false;

        // Now check ris with no remaining pairs to j.
        for (int ri = 0; ri < lenri; ri++) {
            if (check(i, ri)) {
                continue;
            }
            boolean pairRemaining = false;
            for (int rj = 0; rj < lenrj; rj++) {
                if (!check(j, rj) && !check(i, ri, j, rj)) {
                    pairRemaining = true;
                    break;
                }
            }
            if (!pairRemaining) {
                if (eliminateRotamer(residues, i, ri, print)) {
                    eliminated = true;
                    logIfMaster(format(" Eliminating rotamer %s-%d with no remaining pairs to residue %s.", residuei.toFormattedString(false, true), ri, residuej));
                } else {
                    logIfMaster(format(" Already eliminated rotamer %s-%d with no remaining pairs to residue %s.", residuei.toFormattedString(false, true), ri, residuej), Level.WARNING);
                }
            }
        }

        // Check rjs with no remaining pairs to i.
        for (int rj = 0; rj < lenrj; rj++) {
            if (check(j, rj)) {
                continue;
            }
            boolean pairRemaining = false;
            for (int ri = 0; ri < lenri; ri++) {
                if (!check(i, ri) && !check(i, ri, j, rj)) {
                    pairRemaining = true;
                    break;
                }
            }
            if (!pairRemaining) {
                if (eliminateRotamer(residues, j, rj, print)) {
                    eliminated = true;
                    logIfMaster(format(" Eliminating rotamer %s-%d with no remaining pairs to residue %s.", residuej.toFormattedString(false, true), rj, residuei));
                } else {
                    logIfMaster(format(" Already eliminated rotamer J %s-%d with no remaining pairs to residue %s.", residuej.toFormattedString(false, true), rj, residuei), Level.WARNING);
                }
            }
        }

        return eliminated;
    }

    /**
     * Checks all i,j to see if for some ri, all rj are eliminated or
     * vice-versa.
     *
     * @param residues Residues under consideration.
     * @return If any singles wereeliminated by pairs.
     */
    private boolean allPairsToSingleEliminations(Residue[] residues) {
        int nRes = residues.length;
        boolean eliminated = false;
        for (int i = 0; i < nRes - 1; i++) {
            for (int j = i + 1; j < nRes; j++) {
                if (pairsToSingleElimination(residues, i, j)) {
                    eliminated = true;
                }
            }
        }
        return eliminated;
    }

    /**
     * Calculates the minimum and maximum summations over additional residues
     * for some pair ri-rj.
     *
     * @param residues Residues under consideration.
     * @param minMax Result array: 0 is min summation, 1 max summation.
     * @param i Residue i.
     * @param ri Rotamer for residue i.
     * @param j Residue j!=i.
     * @param rj Rotamer for residue j.
     * @return False if ri-rj always clashes with other residues.
     * @throws IllegalArgumentException If ri, rj, or ri-rj eliminated.
     */
    private boolean minMaxE2(Residue[] residues, double[] minMax, int i, int ri, int j, int rj) throws IllegalArgumentException {
        Residue resi = residues[i];
        Residue resj = residues[j];
        if (check(i, ri) || check(j, rj) || check(i, ri, j, rj)) {
            throw new IllegalArgumentException(String.format(" Called for minMaxE2 on an eliminated pair %s-%d %s-%d", resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj));
        }

        // Minimum summation over third residues k.
        minMax[0] = 0;
        // Maximum summation over third residues k.
        minMax[1] = 0;

        int nRes = residues.length;
        for (int k = 0; k < nRes; k++) {
            if (k == i || k == j) {
                continue;
            }
            Residue resk = residues[k];
            Rotamer[] rotsk = resk.getRotamers(library);
            int lenrk = rotsk.length;
            double[] minMaxK = new double[2];
            minMaxK[0] = Double.MAX_VALUE;
            minMaxK[1] = Double.MIN_VALUE;

            for (int rk = 0; rk < lenrk; rk++) {
                if (check(k, rk)) {
                    // Not a valid part of phase space.
                    continue;
                }
                if (check(i, ri, k, rk) || check(j, rj, k, rk)) {
                    // Not implemented: check(i, ri, j, rj, k, rk).

                    // i,ri or j,rj clashes with this rotamer, max will be NaN.
                    // Minimum for this rk will be a clash, which is never a minimum.
                    minMaxK[1] = Double.NaN;
                } else {

                    // Min and max summations over 4th residues l, plus the ri-rk and rj-rk interactions.
                    // If no 3-body term, just the ri-rk and rj-rk interactions.
                    double currentMin = get2Body(i, ri, k, rk) + get2Body(j, rj, k, rk);
                    double currentMax = currentMin;
                    if (threeBodyTerm) {
                        // If the 3-Body eliminated, would fill max to Double.NaN.
                        currentMin += get3Body(residues, i, ri, j, rj, k, rk);
                        currentMax = currentMin;

                        // Obtain min and max summations over l.
                        double[] minMaxTriple = new double[2];
                        if (minMaxE3(residues, minMaxTriple, i, ri, j, rj, k, rk)) {
                            // A non-finite triples minimum should have the code taking the else branch.
                            assert (Double.isFinite(minMaxTriple[0]) && minMaxTriple[0] != Double.MAX_VALUE);

                            // Add the min and max summations over all 4th residues l.
                            currentMin += minMaxTriple[0];

                            if (Double.isFinite(currentMax) && Double.isFinite(minMaxTriple[1])) {
                                currentMax += minMaxTriple[1];
                            } else {
                                currentMax = Double.NaN;
                            }
                        } else {
                            // i, ri, j, rj, k, rk creates an inevitable clash with some residue l.
                            currentMin = Double.NaN;
                            currentMax = Double.NaN;
                        }
                    }

                    assert (threeBodyTerm || currentMax == currentMin);

                    // Now check if rk displaces previously searched rk for min/max over this k.
                    if (Double.isFinite(currentMin) && currentMin < minMaxK[0]) {
                        // rk has a more favorable minimum than previously searched rk.
                        minMaxK[0] = currentMin;
                    } // Else, no new minimum found.

                    if (Double.isFinite(currentMax) && Double.isFinite(minMaxK[1])) {
                        // rk has a less favorable maximum than previously searched rk.
                        minMaxK[1] = (currentMax > minMaxK[1]) ? currentMax : minMaxK[1];
                    } else {
                        // Our maximum is a NaN.
                        minMaxK[1] = Double.NaN;
                    }
                }
            }

            if (Double.isFinite(minMaxK[0])) {
                // Add the minimum contribution from this k to the summation.
                minMax[0] += minMaxK[0];
            } else {
                // Else, ri-rj conflicts with all rk for this k, and can be swiftly eliminated.
                minMax[0] = Double.NaN;
                minMax[1] = Double.NaN;
                return false;
            }
            if (Double.isFinite(minMaxK[1]) && Double.isFinite(minMax[1])) {
                // Add the max contribution from this k to the summation.
                minMax[1] += minMaxK[1];
            } else {
                // Otherwise, the max for ri-rj is a clash.
                minMax[1] = Double.NaN;
            }
        }

        return Double.isFinite(minMax[0]);
    }

    /**
     * Calculates the minimum and maximum summations over additional residues
     * for some rotamer triples ri-rj-rk.
     *
     * @param residues Residues under consideration.
     * @param minMax Result array: 0 is min summation, 1 max summation.
     * @param i Residue i.
     * @param ri Rotamer for residue i.
     * @param j Residue j!=i.
     * @param rj Rotamer for residue j.
     * @param k Residue k!=j and k!=i.
     * @param rk Rotamer for residue k.
     * @return False if ri-rj-rk always clashes with other residues.
     * @throws IllegalArgumentException if there are pre-existing eliminations
     * in ri-rj-rk.
     */
    private boolean minMaxE3(Residue[] residues, double[] minMax, int i, int ri, int j, int rj, int k, int rk) throws IllegalArgumentException {
        Residue resi = residues[i];
        Residue resj = residues[j];
        Residue resk = residues[k];
        if (check(i, ri) || check(j, rj) || check(k, rk) || check(i, ri, j, rj) || check(i, ri, k, rk) || check(j, rj, k, rk)) {
            // Not implemented: check(i, ri, j, rj, k, rk).
            throw new IllegalArgumentException(String.format(" Called for minMaxE2 on an eliminated triple %s-%d %s-%d %s-%d", resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk));
        }

        // These two are a summation of mins/maxes over all fourth residues l.
        minMax[0] = 0;
        minMax[1] = 0;
        int nRes = residues.length;
        for (int l = 0; l < nRes; l++) {
            if (l == i || l == j || l == k) {
                continue;
            }
            Residue resl = residues[l];
            Rotamer[] rotsl = resl.getRotamers(library);
            int lenrl = rotsl.length;

            // Find min/max rl for residue l.
            double currentMax = Double.MIN_VALUE;
            double currentMin = Double.MAX_VALUE;

            for (int rl = 0; rl < lenrl; rl++) {
                if (check(l, rl) || check(k, rk, l, rl)) {
                    // Not valid phase space for anything.
                    continue;
                }

                double current;
                if (check(i, ri, l, rl) || check(j, rj, l, rl)) {
                    // Not implemented: checking ri-rj-rl, ri-rk-rl, rj-rk-rl, or ri-rj-rk-rl.
                    current = Double.NaN;
                } else {
                    // ri-rj-rl is accounted for at a different part of the summation as ri-rj-rk.
                    current = get3Body(residues, i, ri, k, rk, l, rl) + get3Body(residues, j, rj, k, rk, l, rl);
                }

                // TODO: Add quads to the DEE summation.
                // Would have to replace "current" with array "currentQuads".
                //double[] minMaxQuads;
                // minMaxE4(args)
                if (Double.isFinite(current) && current < currentMin) {
                    // rl forms a more favorable 3-body than any prior rl for this residue l.
                    currentMin = current;
                }

                if (Double.isFinite(current) && Double.isFinite(currentMax)) {
                    if (current > currentMax) {
                        currentMax = current;
                    } // Else, no new finite max found.
                } else {
                    currentMax = Double.NaN;
                }
            }

            if (Double.isFinite(currentMin)) {
                minMax[0] += currentMin;
            } else {
                // Else, ri-rj-rk inevitably conflicts with l.
                minMax[0] = Double.NaN;
                minMax[1] = Double.NaN;
                return false;
            }

            if (Double.isFinite(currentMax) && Double.isFinite(minMax[1])) {
                minMax[1] += currentMax;
            } else {
                minMax[1] = Double.NaN;
            }
            // Finished with this residue l.
        }
        return Double.isFinite(minMax[0]);
    }

    /**
     * Computes the maximum and minimum energy i,ri might have with j, and
     * optionally (if three-body energies in use) third residues k.
     * <p>
     * The return value should be redundant with minMax[0] being NaN.
     *
     * @param residues Array of residues under consideration.
     * @param minMax Index 0 to be filled by minimum energy, index 1 filled by
     * maximum energy.
     * @param i Some residue i under consideration.
     * @param ri A rotamer for residue i.
     * @param j Some arbitrary residue i!=j.
     * @return If a valid configuration between i,ri and j could be found.
     */
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

            double currMax = get2Body(i, ri, j, rj);
            double currMin = currMax; // Will remain identical if truncating at 2-body.

            if (threeBodyTerm) {
                double minMaxTriple[] = new double[2];
                // Loop over residue k to find the min/max 3-Body energy.
                boolean validPair = minMax2BodySum(residues, minMaxTriple, i, ri, j, rj);
                if (!validPair) {
                    // Eliminate Rotamer Pair
                    Residue residuei = residues[i];
                    logIfMaster(format(" Inconsistent Pair: %8s %2d, %8s %2d.", residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj), Level.INFO);
                    continue;
                }

                if (Double.isFinite(currMin) && Double.isFinite(minMaxTriple[0])) {
                    currMin += minMaxTriple[0];
                } else {
                    currMin = Double.NaN;
                }

                if (Double.isFinite(currMax) && Double.isFinite(minMaxTriple[1])) {
                    currMax += minMaxTriple[1];
                } else {
                    currMax = Double.NaN;
                }
            }

            valid = true;
            if (Double.isFinite(currMin) && currMin < minMax[0]) {
                minMax[0] = currMin;
            } // Else, we do not have a new minimum.

            if (Double.isFinite(currMax) && Double.isFinite(minMax[1])) {
                if (currMax > minMax[1]) {
                    // We have a new, finite maximum.
                    minMax[1] = currMax;
                }  // Else, if currMax is finite and less than minMax[1], we do not have a new maximum.
            } else {
                // We have a non-finite maximum.
                minMax[1] = Double.NaN;
            }
        }

        // minMax[0] being set to NaN should be redundant with valid being false.
        // It would indicate i,ri clashes with something in every possible configuration.
        minMax[0] = (minMax[0] == Double.MAX_VALUE) ? Double.NaN : minMax[0];
        // minMax[1] always gets set, unless somehow everything turns up as Double.MIN_VALUE.
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
                if (check(k, rk)) {
                    // k,rk is part of no valid phase space, so ignore it.
                    continue;
                }
                if (check(i, ri, k, rk) || check(j, rj, k, rk)) {
                    // Not implemented: check(i, ri, j, rj, k, rk).
                    // k,rk conflicts with i,ri or j,rj, so the max is now Double.NaN. No effect on minimum.
                    currentMax = Double.NaN;
                } else {
                    double current = get3Body(residues, i, ri, j, rj, k, rk);
                    if (Double.isFinite(current) && current < currentMin) {
                        currentMin = current;
                    } // Else, no new minimum found.
                    if (Double.isFinite(current) && Double.isFinite(currentMax)) {
                        if (current > currentMax) {
                            currentMax = current;
                        } // Else, we have failed to find a new finite maximum.
                    } else {
                        // The maximum is NaN.
                        currentMax = Double.NaN;
                    }
                }
            }
            if (currentMin == Double.MAX_VALUE || !Double.isFinite(minSum)) {
                // We have failed to find a viable configuration for i,ri,j,rj, as it conflicts with all rk for this k.
                minMax[0] = Double.NaN;
                minMax[1] = Double.NaN;
                return false;
            } else {
                // Add finite current min to minSum.
                minSum += currentMin;
            }
            if (Double.isFinite(maxSum) && Double.isFinite(currentMax)) {
                maxSum += currentMax;
            } else {
                maxSum = Double.NaN;
            }
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
        if (eliminated == false) {
            logIfMaster(" No more single rotamers to eliminate.");
        }
        return eliminated;
    }

    /**
     * Attemps to eliminate rotamer riA based on riB.
     *
     * @param residues
     * @param i
     * @param riA Rotamer to attempt elimination of.
     * @param riB Rotamer to attempt elimination by.
     * @return If riA was eliminated.
     */
    private boolean goldsteinElimination(Residue residues[], int i, int riA, int riB) {
        int nres = residues.length;
        Residue resi = residues[i];

        // Initialize Goldstein inequality.
        double selfDiff = getSelf(i, riA) - getSelf(i, riB);
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

                if (check(i, riA, j, rj)) {
                    continue; // This is not a part of configuration space accessible to riA.
                }
                if (check(i, riB, j, rj)) {
                    /**
                     * This is a part of configuration space where riA is valid
                     * but not riB. Thus, if j,rj is part of the GMEC, riB is
                     * inconsistent with it. Thus, riB cannot be used to
                     * eliminate riA.
                     */
                    return false;
                }

                double pairI = get2Body(i, riA, j, rj);
                double pairJ = get2Body(i, riB, j, rj);
                double pairDiff = pairI - pairJ;

                rjEvals++;

                // Include three-body interactions.
                double tripleDiff = 0.0;
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
                            /**
                             * If k,rk or j,rj-k,rk are not a part of valid
                             * configuration space, continue. If i,riA-k,rk or
                             * i,riA-j,rj-k,rk are not valid for riA, continue.
                             */
                            if (check(k, rk) || check(j, rj, k, rk) || check(i, riA, k, rk)) {
                                // Not yet implemented: check(i, riA, j, rj, k, rk) because no triples get eliminated.
                                continue;
                            }
                            /**
                             * If i,riB-k,rk or i,riB-j,rj-k,rk are invalid for
                             * riB, there is some part of configuration space
                             * for which riA is valid but not riB.
                             */
                            if (check(i, riB, k, rk)) {
                                // Not yet implemented: check(i, riB, j, rj, k, rk).
                                return false;
                            }

                            rkEvals++;
                            double e = get3Body(residues, i, riA, j, rj, k, rk) - get3Body(residues, i, riB, j, rj, k, rk);
                            if (e < minForResK) {
                                minForResK = e;
                            }
                        }
                        /**
                         * If there were no 3-body interactions with residue k,
                         * then minForResk is zero.
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
                logIfMaster(format("  Rotamer elimination of (%8s,%2d) by (%8s,%2d): %12.4f > %6.4f.",
                        resi.toFormattedString(false, true), riA, resi.toFormattedString(false, true), riB, goldsteinEnergy, ensembleBuffer));
                logIfMaster(format("   Self: %12.4f, Pairs: %12.4f, Triples: %12.4f.",
                        selfDiff, sumPairDiff, sumTripleDiff));
                return true;
            }
        }

        return false;
    }

    /**
     * Finds and eliminates rotamer pairs according to the many-body Goldstein
     * pairs criterion.
     *
     * @param residues Residues under consideration.
     * @return If any rotamer pairs were eliminated.
     */
    private boolean goldsteinPairDriver(Residue[] residues) {
        int nRes = residues.length;
        boolean eliminated = false;

        // First, generate pairs riA-rjC.
        for (int i = 0; i < nRes; i++) {
            Residue resi = residues[i];
            Rotamer[] rotsi = resi.getRotamers(library);
            int lenri = rotsi.length;

            for (int riA = 0; riA < lenri; riA++) {
                // Don't try to eliminate that which is already eliminated.
                if (check(i, riA)) {
                    continue;
                }

                // Residue j can be any other residue, including ones before residue i.
                for (int j = 0; j < nRes; j++) {
                    // Residue j must be distinct from i.
                    if (i == j) {
                        continue;
                    }
                    Residue resj = residues[j];
                    Rotamer[] rotsj = resj.getRotamers(library);
                    int lenrj = rotsj.length;
                    for (int rjC = 0; rjC < lenrj; rjC++) {
                        // Again, no point in eliminating the already-eliminated.
                        if (check(j, rjC) || check(i, riA, j, rjC)) {
                            continue;
                        }
                        boolean breakOut = false;

                        // Now, generate pairs riB-rjD. If any pair riB-rjD eliminates riA-rjC, break out of the loop.
                        for (int riB = 0; riB < lenri; riB++) {
                            if (breakOut) {
                                break;
                            }
                            if (check(i, riB)) {
                                continue;
                            }
                            for (int rjD = 0; rjD < lenrj; rjD++) {
                                if (breakOut) {
                                    break;
                                }
                                // Do not attempt eliminating with an eliminated pair.
                                if (check(j, rjD) || check(i, riB, j, rjD)) {
                                    continue;
                                }
                                // Do not attempt to eliminate a pair with itself.
                                if (riA == riB && rjC == rjD) {
                                    continue;
                                }
                                if (goldsteinPairElimination(residues, i, riA, riB, j, rjC, rjD)) {
                                    breakOut = true;
                                    eliminated = true;
                                }
                            }
                        }
                    }
                    if (pairsToSingleElimination(residues, i, j)) {
                        eliminated = true;
                    }
                }
            }
        }
        return eliminated;
    }

    /**
     * Attempt to eliminate rotamer pair (ResI-RotA, ResJ-RotC) using
     * (ResI-RotB, ResJ-RotD).
     *
     * @param residues
     * @param i Index of the first residue.
     * @param riA Index of the first residue's rotamer to eliminate.
     * @param riB Index of the first residue's rotamer to use for elimination.
     * @param j Index of the 2nd residue.
     * @param rjC Index of the 2nd residue's rotamer to eliminate.
     * @param rjD Index of the 2nd residue's rotamer to use for elimination.
     * @return Return true if eliminated.
     */
    private boolean goldsteinPairElimination(Residue residues[],
            int i, int riA, int riB, int j, int rjC, int rjD) {

        ArrayList<Residue> missedResidues = null;
        // Initialize the Goldstein energy.
        double goldsteinEnergy = getSelf(i, riA) + getSelf(j, rjC) + get2Body(i, riA, j, rjC)
                - getSelf(i, riB) - getSelf(j, rjD) - get2Body(i, riB, j, rjD);

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
        if (missedResidues != null && !missedResidues.isEmpty()) {
            logIfMaster(format(" Skipping energy comparison due to a missed residue: i %d riA %d riB %d j %d rjC %d rjD %d", i, riA, riB, j, rjC, rjD), Level.FINE);
            return false;
        }

        if (goldsteinEnergy > ensembleBuffer) {
            if (missedResidues.isEmpty()) {
                if (eliminateRotamerPair(residues, i, riA, j, rjC, print)) {
                    logIfMaster(format("  Pair elimination of [(%8s,%2d),(%8s,%2d)] by [(%8s,%2d),(%8s,%2d)]: %12.4f > %6.4f",
                            residues[i].toFormattedString(false, true), riA, residues[j].toFormattedString(false, true), rjC, residues[i].toFormattedString(false, true), riB, residues[j].toFormattedString(false, true), rjD, goldsteinEnergy, ensembleBuffer));
                    return true;
                }
            } else {
                logIfMaster(format("  No Pair elimination of [(%8s,%2d),(%8s,%2d)] by [(%8s,%2d),(%8s,%2d)]: %12.4f > %6.4f",
                        residues[i].toFormattedString(false, true), riA, residues[j].toFormattedString(false, true), rjC, residues[i].toFormattedString(false, true), riB, residues[j].toFormattedString(false, true), rjD, goldsteinEnergy, ensembleBuffer));
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
                // Continue if k,rk invalid for riA/rjC.
                if (check(i, riA, k, rk) || check(j, rjC, k, rk)) {
                    // Not implemented: check(i, riA, j, rjC, k, rk).
                    continue;
                }
                // Return false if k,rk invalid for riB/rjD.
                if (check(i, riB, k, rk) || check(j, rjD, k, rk)) {
                    blockedResidues.add(resk);
                    return Double.NaN;
                }

                rkEvals++;
                double currentResK = get2Body(i, riA, k, rk) - get2Body(i, riB, k, rk) + get2Body(j, rjC, k, rk) - get2Body(j, rjD, k, rk);
                // Include 3-body effects.
                if (threeBodyTerm) {
                    double sumOverL = (get3Body(residues, i, riA, j, rjC, k, rk) - get3Body(residues, i, riB, j, rjD, k, rk));
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
                            // If not a part of valid phase space for riA/rjC, continue.
                            if (check(l, rl) || check(k, rk, l, rl) || check(i, riA, l, rl) || check(j, rjC, l, rl)) {
                                // Not implemented: check(i, riA, j, rjC, l, rl) || check(i, riA, k, rk, l, rl) || check(j, rjC, k, rk, l, rl) || check(i, riA, j, rjC, k, rk, l, rl)
                                continue;
                            }
                            if (check(i, riB, l, rl) || check(j, rjD, l, rl)) {
                                // Not implemented: check(i, riB, j, rjD, l, rl) || check(i, riB, k, rk, l, rl) || check(j, rjD, k, rk, l, rl) || check(i, riB, j, rjD, k, rk, l, rl)
                                blockedResidues.add(residuel);
                                return Double.NaN;
                            }
                            rlEvaluations++;
                            double e = get3Body(residues, i, riA, k, rk, l, rl) - get3Body(residues, i, riB, k, rk, l, rl)
                                    + get3Body(residues, j, rjC, k, rk, l, rl) - get3Body(residues, j, rjD, k, rk, l, rl);
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

    /**
     * Safe method to eliminate a rotamer: will not eliminate if there are no
     * alternate rotamers for residue i, or if i-ri is already eliminated.
     *
     * @param residues Residues under consideration.
     * @param i A residue index.
     * @param ri A rotamer to attempt elimination of.
     * @param verbose
     * @return If the rotamer was eliminated.
     */
    private boolean eliminateRotamer(Residue[] residues, int i, int ri, boolean verbose) {
        // Check if rotamer (i, ri) has already been eliminated.
        if (check(i, ri)) {
            return false;
        }

        // Make sure at least one rotamer rii != ri is left.
        int[] validRots = rotamerCount(residues, i);
        int rotCount = 0;
        for (int rii = 0; rii < validRots.length; rii++) {
            if (rii != ri) {
                ++rotCount;
            }
        }

        if (rotCount == 0) {
            // No valid rotamers other than ri are left!
            return false;
        }

        eliminatedSingles[i][ri] = true;

        //adds eliminatedSingles indices to arraylist
        if (eliminatedSingles[i][ri]) {
            eliminatedResidue.add(i);
            eliminatedRotamer.add(ri);
        }

        if (verbose) {
            logIfMaster(format(" Rotamer (%8s,%2d) eliminated (%2d left).", residues[i].toFormattedString(false, true), ri, rotCount));
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
            //adds eliminatedPairs indices to arraylist
            if (eliminatedPairs[i][ri][j][rj]) {
                eliminatedResidue1.add(i);
                eliminatedRotamer1.add(ri);
                eliminatedResidue2.add(j);
                eliminatedRotamer2.add(rj);
            }
            if (verbose) {
                logIfMaster(format("  Rotamer pair eliminated: [(%8s,%2d) (%8s,%2d)]", residues[i].toFormattedString(false, true), ri, residues[j].toFormattedString(false, true), rj));
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
     * @param i The residue number to examine.
     * @return The remaining valid rotamers.
     */
    private int[] rotamerCount(Residue residues[], int i) {
        int nRes = residues.length;
        Rotamer rotI[] = residues[i].getRotamers(library);
        int ni = rotI.length;

        if (maxRotCheckDepth == 0) {
            // Short-circuit on all its rotamers.
            return IntStream.range(0, ni).toArray();
        }

        return IntStream.range(0, ni).filter((int ri) -> {
            if (check(i, ri)) {
                return false;
            }
            if (maxRotCheckDepth > 1) {
                // Check that rotamer ri has valid pairs with all other residues.
                for (int j = 0; j < nRes; j++) {
                    if (i == j) {
                        continue;
                    }
                    if (rotamerPairCount(residues, i, ri, j) == 0) {
                        return false;
                    }
                }
            }
            return true;
        }).toArray();
    }

    /**
     * Validate residue i with rotamer ri.
     *
     * @param residues The residues being optimized.
     * @param i The residue to validate.
     * @param ri The rotamer to validate.
     * @return The status of this rotamer.
     */
    private boolean validRotamer(Residue residues[], int i, int ri) {
        // Return false if this rotamer has been eliminated.
        if (check(i, ri)) {
            return false;
        }

        if (maxRotCheckDepth > 1) {
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
        }
        return true;
    }

    /**
     * Validate rotamer pair (i, ri) and (j, rj).
     *
     * @param residues The residues being optimized.
     * @param i The first residue to validate.
     * @param ri The first rotamer to validate.
     * @param j The 2nd residue to validate.
     * @param rj The 2nd rotamer to validate.
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

        if (maxRotCheckDepth > 1) {
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
        }
        return true;
    }

    /**
     * Count the rotamer pairs remaining for (residue i, rotamer ri) and residue
     * j.
     *
     * @param residues Residue array.
     * @param i The first residue to examine.
     * @param ri The rotamer for the first residue.
     * @param j The second residue to examine.
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
                if (maxRotCheckDepth > 2) {
                    for (int k = 0; k < nRes; k++) {
                        if (k == i || k == j) {
                            continue;
                        }
                        if (rotamerTripleCount(residues, i, ri, j, rj, k) == 0) {
                            valid = false;
                        }
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
     * Count the rotamer triples remaining for (residue i, rotamer ri) and
     * (residue j, rotamer rj) with residue k.
     *
     * @param residues Residue array.
     * @param i The first residue to examine.
     * @param ri The rotamer for the first residue.
     * @param j The second residue to examine.
     * @param rj The rotamer for the first residue.
     * @param k The third residue.
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

    /**
     * Return a previously computed self-energy.
     *
     * @param i Residue i.
     * @param ri Rotamer ri of residue i.
     * @return The self-energy.
     */
    public double getSelf(int i, int ri) {
        try {
            return selfEnergy[i][ri];
        } catch (NullPointerException npe) {
            logger.info(format(" NPE for self energy (%3d,%2d).", i, ri));
            throw npe;
        }
    }

    /**
     * Return a previously computed 2-body energy.
     *
     * @param i Residue i.
     * @param ri Rotamer ri of residue i.
     * @param j Residue j.
     * @param rj Rotamer rj of residue j.
     * @return The 2-Body energy.
     */
    public double get2Body(int i, int ri, int j, int rj) {
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

    /**
     * Return a previously computed 3-body energy.
     *
     * @param i Residue i.
     * @param ri Rotamer ri of residue i.
     * @param j Residue j.
     * @param rj Rotamer rj of residue j.
     * @param k Residue k.
     * @param rk Rotamer rk of residue k.
     * @return The 3-Body energy.
     */
    public double get3Body(Residue[] residues, int i, int ri, int j, int rj, int k, int rk) {
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
        IntegerKeyset ijk = new IntegerKeyset(i, ri, j, rj, k, rk);
        int indexI = allResiduesList.indexOf(residues[i]);
        int indexJ = allResiduesList.indexOf(residues[j]);
        int indexK = allResiduesList.indexOf(residues[k]);
        if (threeBodyEnergies.containsKey(ijk)) {
            return threeBodyEnergies.getOrDefault(ijk, 0);
        } else if (checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
            return 0.00000000;
        } else {
            String message = String.format(" Could not find an energy for 3-body energy (%3d,%2d) (%3d,%2d) (%3d,%2d)", i, ri, j, rj, k, rk);
            logger.info(message);
            throw new IllegalArgumentException(message);
        }
    }

    /**
     * TODO: Implement Quad energy.
     *
     * @param i Residue i.
     * @param ri Rotamer ri of residue i.
     * @param j Residue j.
     * @param rj Rotamer rj of residue j.
     * @param k Residue k.
     * @param rk Rotamer rk of residue k.
     * @param l Residue l.
     * @param rl Rotamer rl of residue l.
     * @return The 4-Body energy.
     */
    private double get4Body(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
        return 0.0;
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
            logger.info(format(" Decreasing increment to match window size %d", windowSize));
            this.increment = windowSize;
        }
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
     * Check for eliminated rotamer; true if eliminated.
     *
     * @param i Residue i.
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
     * @param i Residue i.
     * @param ri Rotamer ri.
     * @param j Residue j.
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
     * Check for pruned rotamer; true if eliminated. Only used during testing.
     *
     * @param i
     * @param ri
     * @return
     */
    protected boolean checkPrunedSingles(int i, int ri) {
        if (onlyPrunedSingles == null) {
            return false;
        }
        return onlyPrunedSingles[i][ri];
    }

    /**
     * Check for pruned rotamer pair; true if eliminated. Only used during
     * testing.
     *
     * @param i Residue i.
     * @param ri Rotamer ri.
     * @param j Residue j.
     * @param rj Rotamer rj.
     * @return
     */
    protected boolean checkPrunedPairs(int i, int ri, int j, int rj) {
        if (onlyPrunedPairs == null) {
            return false;
        }

        if (j < i) {
            int ii = i;
            int iri = ri;
            i = j;
            ri = rj;
            j = ii;
            rj = iri;
        }
        return onlyPrunedPairs[i][ri][j][rj];
    }

    /**
     * Checks to see if any eliminations with j,rj have occurred; assumes i,ri
     * self has already been checked. Checks j,rj self and i,ri,j,rj 2-Body. The
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
        return (check(j, rj) || check(i, ri, j, rj));
    }

    /**
     * Checks to see if any eliminations with k,rk have occurred; assumes
     * i,ri,j,rj 2-Body has already been checked. Checks the k,rk self, all
     * pairs with k,rk, and the i,ri,j,rj,k,rk 3-Body. The intent is to be part
     * of a loop over i,ri,j,rj,k,rk, and check for eliminations at the k,rk
     * point.
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
        return (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk));
    }

    /**
     * Checks to see if any eliminations with l,rl have occurred; assumes
     * i,ri,j,rj,k,rk 3-Body has already been checked. Checks the l,rl self, all
     * pairs with l,rl, all triples with l,rl, and the 4-Body. The intent is to
     * be part of a loop over i,ri,j,rj,k,rk,l,rl, and check for eliminations at
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
                            residuei.toFormattedString(false, true), residuej));
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
            logIfMaster(format(" Backbone energy:  %s\n", formatEnergy(backboneEnergy)));

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
                selfEnergyMap.clear();
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
                        selfEnergyMap.put(singleJobIndex, selfJob);
                        String revKey = format("%d %d", i, ri);
                        reverseJobMapSingles.put(revKey, singleJobIndex);
                        singleJobIndex++;
                    }
                }
                // fill in self-energies from file while removing the corresponding jobs from selfEnergyMap
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
                                logIfMaster(format(" From restart file: Self energy %3d (%8s,%2d): %s", i, residues[i].toFormattedString(false, true), ri,
                                        formatEnergy(energy)));
                            }
                        } catch (Exception e) {
                            if (verbose) {
                                logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                            }
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d", i, ri);
                        Integer ret[] = selfEnergyMap.remove(reverseJobMapSingles.get(revKey));
                        if (ret == null) {
                            //logIfMaster(format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format(" Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded self energies from restart file.");

                //Pre-Prune if self-energy is Double.NaN.
                prePruneSelves(residues);

                // prune singles
                if (pruneClashes) {
                    pruneSingleClashes(residues);
                }
            }

            /**
             * Remap to sequential integer keys.
             */
            condenseEnergyMap(selfEnergyMap);

            if (loaded >= 2) {
                if (selfEnergyMap.size() > 0) {
                    logIfMaster(" Double-check that parameters match original run due to missing self-energies.");
                }
                twoBodyEnergyMap.clear();
                // allocated twoBodyEnergy array and create pair jobs
                HashMap<String, Integer> reverseJobMapPairs = new HashMap<>();
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
                                if (check(j, rj) || check(i, ri, j, rj)) {
                                    continue;
                                }
                                Integer pairJob[] = {i, ri, j, rj};
                                if (decomposeOriginal && (ri != 0 || rj != 0)) {
                                    continue;
                                }
                                twoBodyEnergyMap.put(pairJobIndex, pairJob);
                                String revKey = format("%d %d %d %d", i, ri, j, rj);
                                reverseJobMapPairs.put(revKey, pairJobIndex);
                                pairJobIndex++;
                            }
                        }
                    }
                }
                // fill in pair-energies from file while removing the corresponding jobs from twoBodyEnergyMap
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
                                logIfMaster(format(" From restart file: Pair energy [(%8s,%2d),(%8s,%2d)]: %12.4f",
                                        residues[i].toFormattedString(false, true), ri, residues[j].toFormattedString(false, true), rj, energy));
                            }
                        } catch (Exception e) {
                            if (verbose) {
                                logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                            }
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d %d %d", i, ri, j, rj);
                        Integer ret[] = twoBodyEnergyMap.remove(reverseJobMapPairs.get(revKey));
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded 2-body energies from restart file.");

                // Pre-Prune if pair-energy is Double.NaN.
                prePrunePairs(residues);

                // prune pairs
                if (prunePairClashes) {
                    prunePairClashes(residues);
                }
            }

            /**
             * Remap to sequential integer keys.
             */
            condenseEnergyMap(twoBodyEnergyMap);

            if (loaded >= 3) {
                if (twoBodyEnergyMap.size() > 0) {
                    if (master) {
                        logger.warning("Double-check that parameters match original run!  Found trimers in restart file, but pairs job queue is non-empty.");
                    }
                }
                HashMap<String, Integer> reverseJobMapTrimers = new HashMap<>();
                threeBodyEnergyMap.clear();
                // Allocate 3-body energy Map, fill in 3-Body energies from the restart file.
                alloc3BodyMap(residues);
                int trimerJobIndex = 0;
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
                                if (check(j, rj) || check(i, ri, j, rj)) {
                                    continue;
                                }
                                for (int k = j + 1; k < nResidues; k++) {
                                    Residue resk = residues[k];
                                    Rotamer rotk[] = resk.getRotamers(library);
                                    for (int rk = 0; rk < rotk.length; rk++) {
                                        if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk)) {
                                            // Not implemented: check(i, ri, j, rj, k, rk).
                                            continue;
                                        }
                                        Integer trimerJob[] = {i, ri, j, rj, k, rk};
                                        if (decomposeOriginal && (ri != 0 || rj != 0 || rk != 0)) {
                                            continue;
                                        }
                                        threeBodyEnergyMap.put(trimerJobIndex, trimerJob);
                                        String revKey = format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                                        reverseJobMapTrimers.put(revKey, trimerJobIndex);
                                        trimerJobIndex++;
                                    }
                                }
                            }
                        }
                    }
                }

                // fill in 3-Body energies from file while removing the corresponding jobs from threeBodyEnergyMap
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
                            //threeBodyEnergy[i][ri][j][rj][k][rk] = energy;
                            threeBodyEnergies.put(new IntegerKeyset(i, ri, j, rj, k, rk), energy);
                        } catch (ArrayIndexOutOfBoundsException ex) {
                            if (verbose) {
                                logIfMaster(format(" Restart file out-of-bounds index: %s", line));
                            }
                        } catch (NullPointerException npe) {
                            if (verbose) {
                                logIfMaster(format(" NPE in loading 3-body energies: pruning "
                                        + "likely changed! 3-body %s-%d %s-%d %s-%d", residues[i].toFormattedString(false, true),
                                        ri, residues[j], rj, residues[k], rk));
                            }
                        }
                        if (verbose) {
                            logIfMaster(format(" From restart file: Trimer energy %3d %-2d, %3d %-2d, %3d %-2d: %s",
                                    i, ri, j, rj, k, rk, formatEnergy(energy)));
                        }
                        // remove that job from the pool
                        String revKey = format("%d %d %d %d %d %d", i, ri, j, rj, k, rk);
                        Integer ret[] = threeBodyEnergyMap.remove(reverseJobMapTrimers.get(revKey));
                        if (ret == null) {
                            //logIfMaster(format("(sdl %d) Restart file contained unnecessary value for %s", BOXNUM, revKey));
                        }
                    } catch (NumberFormatException ex) {
                        logger.log(Level.WARNING, format("Unparsable line in energy restart file: \n%s", line), ex);
                    }
                }
                logIfMaster(" Loaded trimer energies from restart file.");
            }

            /**
             * Remap to sequential integer keys.
             */
            condenseEnergyMap(threeBodyEnergyMap);

            return loaded;
        } catch (IOException ex) {
            logger.log(Level.WARNING, "Exception while loading energy restart file.", ex);
        }

        return 0;
    }

    private void condenseEnergyMap(HashMap<Integer, Integer[]> energyMap) {
        Set<Integer> keys = energyMap.keySet();
        HashMap<Integer, Integer[]> tempMap = new HashMap<>();
        int count = 0;
        for (int key : keys) {
            tempMap.put(count, energyMap.get(key));
            count++;
        }
        energyMap.clear();
        energyMap.putAll(tempMap);
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
                if (!check(i, ri) && getSelf(i, ri) < minEnergy) {
                    minEnergy = getSelf(i, ri);
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
            energyToPrune = (residue.getResidueType() == NA) ? energyToPrune * nucleicPruningFactor : energyToPrune;
            energyToPrune += minEnergy;

            for (int ri = 0; ri < nrot; ri++) {
                if (!check(i, ri) && (getSelf(i, ri) > energyToPrune)) {
                    if (eliminateRotamer(residues, i, ri, print)) {
                        logIfMaster(format("  Rotamer (%7s,%2d) self-energy %s pruned by (%7s,%2d) %s.",
                                residue, ri, formatEnergy(getSelf(i, ri)),
                                residue, minRot, formatEnergy(minEnergy)));
                    }
                }
            }
        }
    }

    /**
     * Allocates the 3-body map from i,ri,j,rj,k,rk keys to 3-body energy
     * values.
     *
     * @return If the 3-body map already existed.
     */
    protected boolean alloc3BodyMap(Residue[] residues) {
        boolean retValue = (threeBodyEnergies != null);

        if (retValue) {
            // The clear method only clears the keys, and does not resize.
            // This should help avoid needless upsizing of the map after the first large box/window.
            threeBodyEnergies.clear();
            return true;
        } else {
            // Magic number of 10000 is a random guess at "big enough to cover all small jobs, not too big to be excessive overhead for small jobs".
            // Estimate maybe 500 kB initial size.
            threeBodyEnergies = new Object2DoubleOpenHashMap<>(10000);
            return false;
        }
    }

    /**
     * Pre-prunes any selves that have a self-energy of Double.NaN before pruning and elminations happen.
     * @param residues Array of all residues.
     */
    private void prePruneSelves(Residue residues[]){
        //Pre-Prune if self-energy is Double.NaN.
        for (int i = 0; i < residues.length; i++) {
            Residue residue = residues[i];
            Rotamer rotamers[] = residue.getRotamers(library);
            int nrot = rotamers.length;
            for (int ri = 0; ri < nrot; ri++) {
                if (!check(i, ri) && Double.isNaN(getSelf(i, ri))) {
                    logIfMaster(format(" Rotamer (%7s,%2d) self-energy %12.4f pre-pruned since energy is NaN.", residue, ri, getSelf(i, ri)));
                    eliminateRotamer(residues, i, ri, false);
                }
            }
        }
    }

    /**
     * Pre-prunes any pairs that have a pair-energy of Double.NaN before pruning and eliminations happen.
     * @param residues Array of all residues.
     */
    private void prePrunePairs(Residue residues[]) {
        // Pre-Prune if pair-energy is Double.NaN.
        // Loop over first residue.
        int nResidues = residues.length;
        for (int i = 0; i < nResidues - 1; i++) {
            Residue resi = residues[i];
            Rotamer[] roti = resi.getRotamers(library);
            int ni = roti.length;
            // Loop over second residue.
            for (int j = i + 1; j < nResidues; j++) {
                Residue resj = residues[j];
                Rotamer[] rotj = resj.getRotamers(library);
                int nj = rotj.length;
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
                        if (!check(i, ri, j, rj) && Double.isNaN(get2Body(i, ri, j, rj))) {
                            logIfMaster(format(" Rotamer Pair (%7s,%2d) (%7s,%2d) 2-body energy %12.4f pre-pruned since energy is NaN.", i, ri, j, rj, get2Body(i, ri, j, rj)));
                            eliminateRotamerPair(residues, i, ri, j, rj, print);
                        }
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

        for (int i = 0; i < nResidues - 1; i++) {
            Residue resi = residues[i];
            Rotamer[] rotsi = resi.getRotamers(library);
            int lenri = rotsi.length;
            int indI = allResiduesList.indexOf(resi);
            for (int j = i + 1; j < nResidues; j++) {
                Residue resj = residues[j];
                Rotamer[] rotsj = resj.getRotamers(library);
                int lenrj = rotsj.length;
                int indJ = allResiduesList.indexOf(resj);

                double minPair = Double.MAX_VALUE;
                int minRI = -1;
                int minRJ = -1;

                boolean cutoffPair = true;
                for (int ri = 0; ri < lenri; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (check(j, rj) || check(i, ri, j, rj) || checkPairDistThreshold(indI, ri, indJ, rj)) {
                            continue;
                        }
                        cutoffPair = false;
                        double pairEnergy = get2Body(i, ri, j, rj) + getSelf(i, ri) + getSelf(j, rj);
                        assert Double.isFinite(pairEnergy);
                        if (pairEnergy < minPair) {
                            minPair = pairEnergy;
                            minRI = ri;
                            minRJ = rj;
                        }
                    }
                }
                if (cutoffPair) {
                    // Under this branch: no rotamer pairs that clear the distance threshold.
                    continue;
                }
                assert (minRI >= 0 && minRJ >= 0); // Otherwise, i and j are not on a well-packed backbone.

                // Calculate all the modifiers to the pair clash elimination threshold.
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
                        threshold *= nucleicPairsPruningFactor;
                        break;
                    case 2:
                        threshold *= nucleicPruningFactor;
                        break;
                    default:
                        throw new ArithmeticException(" RotamerOptimization.prunePairClashes() has somehow "
                                + "found less than zero or more than two nucleic acid residues in a pair of"
                                + " residues. This result should be impossible.");
                }
                double toEliminate = threshold + minPair;

                for (int ri = 0; ri < lenri; ri++) {
                    if (check(i, ri)) {
                        continue;
                    }
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (check(j, rj) || check(i, ri, j, rj)) {
                            continue;
                        }
                        double pairEnergy = get2Body(i, ri, j, rj) + getSelf(i, ri) + getSelf(j, rj);
                        assert Double.isFinite(pairEnergy);
                        if (pairEnergy > toEliminate) {
                            logIfMaster(format(" Pruning pair %s-%d %s-%d by %s-%d %s-%d; energy %s > "
                                    + "%s + %s", resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resi.toFormattedString(false, true), minRI, resj.toFormattedString(false, true), minRJ,
                                    formatEnergy(pairEnergy), formatEnergy(threshold), formatEnergy(minPair)));
                        }
                    }
                }
                pairsToSingleElimination(residues, i, j);
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

    public void setTwoBodyCutoff(double twoBodyCutoffDist) {
        this.twoBodyCutoffDist = twoBodyCutoffDist;
    }

    public void setThreeBodyCutoff(double threeBodyCutoffDist) {
        this.threeBodyCutoffDist = threeBodyCutoffDist;
    }

    public void setTestOverallOpt(boolean testing) {
        this.testing = testing;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    public void setTestSelfEnergyEliminations(boolean testSelfEnergyEliminations) {
        this.testing = true;
        this.testSelfEnergyEliminations = testSelfEnergyEliminations;
        testPairEnergyEliminations = -1;
        testTripleEnergyEliminations1 = -1;
        testTripleEnergyEliminations2 = -1;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    public void setTestPairEnergyEliminations(int testPairEnergyEliminations) {
        this.testing = true;
        this.testPairEnergyEliminations = testPairEnergyEliminations;
        testSelfEnergyEliminations = false;
        testTripleEnergyEliminations1 = -1;
        testTripleEnergyEliminations2 = -1;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    public void setTestTripleEnergyEliminations(int testTripleEnergyEliminations1, int testTripleEnergyEliminations2) {
        this.testing = true;
        this.testTripleEnergyEliminations1 = testTripleEnergyEliminations1;
        this.testTripleEnergyEliminations2 = testTripleEnergyEliminations2;
        testSelfEnergyEliminations = false;
        testPairEnergyEliminations = -1;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    public enum Algorithm {

        ALL, BOX, WINDOW, INDEPENDENT, BRUTE_FORCE
    }

    public enum Direction {

        FORWARD, BACKWARD
    }

    private class EnergyRegion extends ParallelRegion {

        private SharedDouble self;
        private SharedDouble twoBody;
        private SharedDouble threeBody;
        private EnergyLoop energyLoops[];
        private Residue[] residues;
        private int rotamers[];
        private int nResidues;

        public EnergyRegion(int nThreads) {
            self = new SharedDouble();
            twoBody = new SharedDouble();
            threeBody = new SharedDouble();
            energyLoops = new EnergyLoop[nThreads];
        }

        public void init(Residue residues[], int rotamers[]) {
            this.rotamers = rotamers;
            this.nResidues = residues.length;
            this.residues = Arrays.copyOf(residues, nResidues);
        }

        public void start() {
            self.set(0.0);
            twoBody.set(0.0);
            threeBody.set(0.0);
        }

        public double getSelf() {
            return self.get();
        }

        public double getTwoBody() {
            return twoBody.get();
        }

        public double getThreeBody() {
            return threeBody.get();
        }

        @Override
        public void run() throws Exception {
            int threadID = getThreadIndex();
            if (energyLoops[threadID] == null) {
                energyLoops[threadID] = new EnergyLoop();
            }
            //try {
            execute(0, nResidues - 1, energyLoops[threadID]);
            /*} catch (Exception e) {
             logger.log(Level.WARNING, " Exception in EnergyLoop.", e);
             }*/
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
                twoBody.addAndGet(pairSum);
                threeBody.addAndGet(threeBodySum);
            }

            @Override
            public void run(int lb, int ub) {
                for (int a = lb; a <= ub; a++) {
                    int ai = rotamers[a];
                    selfSum += RotamerOptimization.this.getSelf(a, ai);
                    for (int b = a + 1; b < nResidues; b++) {
                        int bi = rotamers[b];
                        pairSum += get2Body(a, ai, b, bi);
                        if (threeBodyTerm) {
                            for (int c = b + 1; c < nResidues; c++) {
                                int ci = rotamers[c];
                                threeBodySum += get3Body(residues, a, ai, b, bi, c, ci);
                            }
                        }
                    }
                }
            }
        }
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

        /**
         * Initializes a ParallelRegion to attempt the elimination of riA,rjC by
         * riB,rjD.
         *
         * @param residues
         * @param i First residue of the pair.
         * @param riA First member of the pair to attempt eliminating.
         * @param riB First member of the pair to try eliminating by.
         * @param j Second residue of the pair.
         * @param rjC Second member of the pair to attempt eliminating.
         * @param rjD Second member of the pair to try eliminating by.
         */
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
                if (blockedResidues.isEmpty()) {
                    //sumOverK += goldsteinPairSumOverK(residues, lb, ub, i, riA, riB, j, rjC, rjD, blockedResidues);
                    double locSumOverK = goldsteinPairSumOverK(residues, lb, ub, i, riA, riB, j, rjC, rjD, blockedResidues);

                    // Should be redundant checks.
                    if (Double.isFinite(locSumOverK) && blockedResidues.isEmpty()) {
                        sumOverK += locSumOverK;
                    } else {
                        sumOverK = 0;
                    }
                } else {
                    logIfMaster(format(" Skipping %d to %d because we cannot eliminate", lb, ub), Level.FINE);
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

        // fracCoords indexed by 1-3 min x,y,z, 4-6 are max x,y,z
        private final double[] fracCoords = new double[6];
        private final int[] indexXYZ = new int[3];
        private final int linearIndex;
        private ArrayList<Residue> residues = new ArrayList<>();

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
         * Sorts residues in the box.
         */
        public void sortBoxResidues() {
            sortResidues(residues);
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
         * Checks if any rotamer of a Residue is inside this BoxOptCell.
         *
         * @param residue Residue to check.
         * @param crystal A Crystal.
         * @param symOp A symmetry operator to apply.
         * @param variableOnly If using only variable (protein side-chain, nucleic acid backbone) atoms.
         * @return If contained inside this BoxOptCell.
         */
        public boolean anyRotamerInsideCell(Residue residue, Crystal crystal, SymOp symOp, boolean variableOnly) {
            ResidueState incomingState = residue.storeState();
            Rotamer[] rotamers = residue.getRotamers(library);
            boolean inside = Arrays.stream(rotamers).anyMatch((Rotamer r) -> {
                RotamerLibrary.applyRotamer(residue, r);
                return residueInsideCell(residue, crystal, symOp, variableOnly);
            });
            residue.revertState(incomingState);
            return inside;
        }

        /**
         * Checks if a Residue is inside this BoxOptCell.
         *
         * @param residue Residue to check.
         * @param crystal A Crystal.
         * @param symOp A symmetry operator to apply.
         * @param variableOnly If using only variable (protein side-chain, nucleic acid backbone) atoms.
         * @return If contained inside this BoxOptCell.
         */
        public boolean residueInsideCell(Residue residue, Crystal crystal, SymOp symOp, boolean variableOnly) {
            List<Atom> atoms = variableOnly ? residue.getVariableAtoms() : residue.getAtomList();
            return atoms.stream().anyMatch(a -> atomInsideCell(a, crystal, symOp));
        }

        /**
         * Checks if an Atom would be contained inside this cell.
         * @param atom Atom to check.
         * @param crystal A Crystal.
         * @param symOp A symmetry operator to apply.
         * @return If contained.
         */
        public boolean atomInsideCell(Atom atom, Crystal crystal, SymOp symOp) {
            double[] atXYZ = new double[3];
            atXYZ = atom.getXYZ(atXYZ);
            crystal.toFractionalCoordinates(atXYZ, atXYZ);
            NeighborList.moveValuesBetweenZeroAndOne(atXYZ);
            crystal.applyFracSymOp(atXYZ, atXYZ, symOp);
            return checkIfContained(atXYZ);
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

    /**
     * Writes files that allow restarting from partially-completed call to
     * rotamerEnergies(). Spawned only in a parallel environment and only by the
     * master process.
     */
    private class EnergyWriterThread extends Thread {

        private final int writeFrequency = 100;
        private ReceiveThread receiveThread;
        private File restartFile;
        private BufferedWriter bw;
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
                energyRestartFile = restartFile;
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
                if (receiveThread.getState() == State.TERMINATED) {
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
                        logger.log(Level.SEVERE, " Exception writing energy restart file.", ex);
                    }
                    writing.clear();
                }
            }
            try {
                bw.close();
            } catch (IOException ex) {
                logger.log(Level.SEVERE, " Exception while closing energy restart file.", ex);
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
            Residue[] resArr = new Residue[residueList.size()];
            resArr = residueList.toArray(resArr);
            while (alive) {
                try {
                    CommStatus cs = world.receive(null, incSelfBuf);
                    if (cs.length > 0) {
                        int resi = (int) incSelf[0];
                        int roti = (int) incSelf[1];
                        double energy = incSelf[2];
                        if (energy == Double.NaN) {
                            logger.info("Rotamer eliminated.");
                            eliminateRotamer(resArr, resi, roti, false);
                        }
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
                    finishedSelf = true;
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
                    readyFor2Body = true;
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
                        if (energy == Double.NaN) {
                            logger.info("Rotamer pair eliminated.");
                            eliminateRotamerPair(resArr, resi, roti, resj, rotj, false);
                        }
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
                    finished2Body = true;
                    break;
                }
            }

            if (threeBodyTerm) {
                // Barrier; wait for everyone to be done pruning and allocating threeBodyEnergy arrays before computing trimers.
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
                        readyFor3Body = true;
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
                                //threeBodyEnergy[resi][roti][resj][rotj][resk][rotk] = energy;
                                threeBodyEnergies.put(new IntegerKeyset(resi, roti, resj, rotj, resk, rotk), energy);
                                if (writeEnergyRestart && printFiles) {
                                    energiesToWrite.add(format("Triple %d %d, %d %d, %d %d: %16.8f", resi, roti, resj, rotj, resk, rotk, energy));
                                }
                            }
                        }
                    } catch (IOException ex) {
                        logger.log(Level.WARNING, ex.getMessage(), ex);
                    }

                    if (procsDone >= numProc) {
                        finished3Body = true;
                        break;
                    }
                }
            }

            if (compute4BodyEnergy) {
                // Barrier; wait for everyone to be done pruning before starting 4-Body energies.
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
                        readyFor4Body = true;
                        break;
                    }
                }
                // No receive mechanism for 4-Body energies since we only print them.
            }
            if (verbose) {
                logger.info(" Receive thread shutting down.");
            }
        }
    }

    /**
     * Compute residue self-energy values in parallel across nodes.
     */
    private class SelfEnergyRegion extends WorkerRegion {

        private final SelfEnergyLoop energyLoop;
        private final Residue residues[];
        Set<Integer> keySet;

        public SelfEnergyRegion(Residue residues[]) {
            energyLoop = new SelfEnergyLoop();
            this.residues = residues;
        }

        @Override
        public void start() {
            /**
             * Load the keySet of self energies.
             */
            keySet = selfEnergyMap.keySet();

            /**
             * Compute backbone energy.
             */
            try {
                backboneEnergy = computeBackboneEnergy(residues);
            } catch (ArithmeticException ex) {
                logger.severe(format(" Error in calculation of backbone energy %s", ex.getMessage()));
            }
            logIfMaster(format(" Backbone energy:  %s\n", formatEnergy(backboneEnergy)));
        }

        @Override
        public void run() throws Exception {
            if (!keySet.isEmpty()) {
                execute(0, keySet.size() - 1, energyLoop);
            }
        }

        @Override
        public void finish() {
            /**
             * Broadcast completion signal.
             */
            double finished[] = {-1.0, -1.0, -1.0};
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            /**
             * Wait for everyone else to send in their self energies.
             */
            int waiting = 0;
            while (!finishedSelf) {
                try {
                    Thread.sleep(POLLING_FREQUENCY);
                    if (waiting++ == 1000) {
                        logger.warning(format(" Process %d experiencing long wait for others' self energies.", Comm.world().rank()));
                    }
                } catch (InterruptedException e) {
                    logger.warning(" Parallel self-energy computation was interrupted.");
                }
            }

            //Pre-Prune if self-energy is Double.NaN.
            prePruneSelves(residues);

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
                        logger.info(format(" Self energy %8s %-2d: %s", residues[i].toFormattedString(false, true), ri, formatEnergy(getSelf(i, ri))));
                    }
                }
            }
        }

        private class SelfEnergyLoop extends WorkerIntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.guided();
                // /return IntegerSchedule.dynamic(1);
            }

            @Override
            public void run(int lb, int ub) {
                for (int key = lb; key <= ub; key++) {
                    Integer job[] = selfEnergyMap.get(key);
                    int i = job[0];
                    int ri = job[1];

                    if (check(i, ri)) {
                        continue;
                    }

                    long time = -System.nanoTime();
                    double selfEnergy = 0.0;
                    try {
                        selfEnergy = computeSelfEnergy(residues, i, ri);
                        time += System.nanoTime();
                        logger.info(format(" Self %8s %-2d: %s in %6.4f (sec).", residues[i].toFormattedString(false, true), ri,
                                formatEnergy(selfEnergy), time * 1.0e-9));
                    } catch (ArithmeticException ex) {
                        eliminateRotamer(residues, i, ri, false);
                        selfEnergy = Double.NaN;
                        time += System.nanoTime();
                        logger.info(format(" Self %8s %-2d:\t    pruned in %6.4f (sec).", residues[i].toFormattedString(false, true), ri, time * 1.0e-9));
                    }

                    // Broadcast the result.
                    double anEnergy[] = {i, ri, selfEnergy};
                    DoubleBuf anEnergyBuf = DoubleBuf.buffer(anEnergy);
                    multicastBuf(anEnergyBuf);
                }
            }
        }
    }

    /**
     * Compute 2-Body energy values in parallel across nodes.
     */
    private class TwoBodyEnergyRegion extends WorkerRegion {

        private final TwoBodyEnergyLoop energyLoop;
        private final Residue residues[];
        Set<Integer> keySet;

        public TwoBodyEnergyRegion(Residue residues[]) {
            energyLoop = new TwoBodyEnergyLoop();
            this.residues = residues;
        }

        @Override
        public void start() {
            keySet = twoBodyEnergyMap.keySet();
        }

        @Override
        public void run() throws Exception {
            if (!keySet.isEmpty()) {
                execute(0, keySet.size() - 1, energyLoop);
            }
        }

        @Override
        public void finish() {
            // Broadcast the completion signal.
            double finished[] = {-1.0, -1.0, -1.0, -1.0, -1.0};
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            // Wait for other nodes to finish.
            int waiting = 0;
            while (!finished2Body) {
                try {
                    Thread.sleep(POLLING_FREQUENCY);
                    if (waiting++ == 1000) {
                        logger.warning(format("Process %d experiencing long wait for 2-body energies from other processes.", Comm.world().rank()));
                    }
                } catch (InterruptedException ex) {
                }
            }

            // Pre-Prune if pair-energy is Double.NaN.
            prePrunePairs(residues);

            // Prune each rotamer that clashes with all rotamers from a 2nd residue.
            if (prunePairClashes) {
                prunePairClashes(residues);
            }

            // Print what we've got so far.
            if (master && verbose) {
                for (int i = 0; i < residues.length; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (check(i, ri)) {
                            continue;
                        }
                        for (int j = i + 1; j < residues.length; j++) {
                            Residue resj = residues[j];
                            Rotamer rotj[] = resj.getRotamers(library);
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if (check(j, rj) || check(i, ri, j, rj)) {
                                    continue;
                                }
                                logger.info(format(" Pair energy %8s %-2d, %8s %-2d: %s",
                                        residues[i].toFormattedString(false, true), ri, residues[j].toFormattedString(false, true), rj, formatEnergy(get2Body(i, ri, j, rj))));
                            }
                        }
                    }
                }
            }
        }

        private class TwoBodyEnergyLoop extends WorkerIntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.guided();
            }

            @Override
            public void run(int lb, int ub) {
                for (int key = lb; key <= ub; key++) {
                    long time = -System.nanoTime();
                    Integer job[] = twoBodyEnergyMap.get(key);
                    int i = job[0];
                    int ri = job[1];
                    int j = job[2];
                    int rj = job[3];

                    if (check(i, ri) || check(j, rj) || check(i, ri, j, rj)) {
                        logger.severe(String.format(" Invalid pair energy job for %d-%d %d-%d: eliminated i %b j %b i-j %b", i, ri, j, rj, check(i, ri), check(j, rj), check(i, ri, j, rj)));
                        continue;
                    }

                    Residue residueI = residues[i];
                    Residue residueJ = residues[j];
                    int indexI = allResiduesList.indexOf(residueI);
                    int indexJ = allResiduesList.indexOf(residueJ);
                    double resDist = getResidueDistance(indexI, ri, indexJ, rj);
                    String resDistString = format("large");
                    if (resDist < Double.MAX_VALUE) {
                        resDistString = format("%5.3f", resDist);
                    }
                    
                    double dist = checkDistMatrix(indexI, ri, indexJ, rj);

                    String distString = format("     large");
                    if (dist < Double.MAX_VALUE) {
                        distString = format("%10.3f", dist);
                    }

                    double twoBodyEnergy = 0.0;
                    if (dist < superpositionThreshold) {
                        // Set the energy to NaN for superposed atoms.
                        twoBodyEnergy = Double.NaN;
                        logger.info(format(" Pair %8s %-2d, %8s %-2d:\t    NaN at %13.6f Ang (%s Ang by residue) < %5.3f Ang",
                                residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, dist, resDist, superpositionThreshold));
                    } else if (checkPairDistThreshold(indexI, ri, indexJ, rj)) {
                        // Set the two-body energy to 0.0 for separation distances larger than the two-body cutoff.
                        twoBodyEnergy = 0.0;
                        time += System.nanoTime();
                        logger.info(format(" Pair %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, formatEnergy(twoBodyEnergy), distString, resDistString, time * 1.0e-9));
                    } else {
                        try {
                            twoBodyEnergy = compute2BodyEnergy(residues, i, ri, j, rj);
                            time += System.nanoTime();
                            logger.info(format(" Pair %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                    residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, formatEnergy(twoBodyEnergy), distString, resDistString, time * 1.0e-9));
                        } catch (ArithmeticException ex) {
                            twoBodyEnergy = Double.NaN;
                            time += System.nanoTime();
                            logger.info(format(" Pair %8s %-2d, %8s %-2d:              NaN at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                    residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, distString, resDistString, time * 1.0e-9));
                        }
                    }

                    if (!Double.isFinite(twoBodyEnergy)) {
                        eliminateRotamerPair(residues, i, ri, j, rj, false);
                    }

                    // Broadcast the energy.
                    double commEnergy[] = {i, ri, j, rj, twoBodyEnergy};
                    DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
                    multicastBuf(commEnergyBuf);
                }
            }
        }
    }

    /**
     * Compute 3-Body energy values in parallel across nodes.
     */
    private class ThreeBodyEnergyRegion extends WorkerRegion {

        private final ThreeBodyEnergyLoop threeBodyEnergyLoop;
        private final Residue residues[];
        Set<Integer> keySet;

        public ThreeBodyEnergyRegion(Residue residues[]) {
            threeBodyEnergyLoop = new ThreeBodyEnergyLoop();
            this.residues = residues;
        }

        @Override
        public void start() {
            keySet = threeBodyEnergyMap.keySet();
        }

        @Override
        public void run() throws Exception {
            if (!keySet.isEmpty()) {
                execute(0, keySet.size() - 1, threeBodyEnergyLoop);
            }
        }

        @Override
        public void finish() {
            // Broadcast completion signal.
            double finished[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
            DoubleBuf finishedBuf = DoubleBuf.buffer(finished);
            multicastBuf(finishedBuf);

            // wait for everyone else to send in their trimer energies
            int waiting = 0;
            while (!finished3Body) {
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
                for (int i = 0; i < residues.length; i++) {
                    Residue resi = residues[i];
                    Rotamer roti[] = resi.getRotamers(library);
                    for (int ri = 0; ri < roti.length; ri++) {
                        if (check(i, ri)) {
                            continue;
                        }
                        for (int j = i + 1; j < residues.length; j++) {
                            Residue resj = residues[j];
                            Rotamer rotj[] = resj.getRotamers(library);
                            for (int rj = 0; rj < rotj.length; rj++) {
                                if (check(j, rj) || check(i, ri, j, rj)) {
                                    continue;
                                }
                                for (int k = j + 1; k < residues.length; k++) {
                                    Residue resk = residues[k];
                                    Rotamer rotk[] = resk.getRotamers(library);
                                    for (int rk = 0; rk < rotk.length; rk++) {
                                        if (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk)) {
                                            continue;
                                        }
                                        logger.info(format(" 3-Body energy %8s %-2d, %8s %-2d, %8s %-2d: %s",
                                                resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk, formatEnergy(get3Body(residues, i, ri, j, rj, k, rk))));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        private class ThreeBodyEnergyLoop extends WorkerIntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.guided();
            }

            @Override
            public void run(int lb, int ub) {
                for (int key = lb; key <= ub; key++) {
                    long time = -System.nanoTime();
                    Integer job[] = threeBodyEnergyMap.get(key);
                    int i = job[0];
                    int ri = job[1];
                    int j = job[2];
                    int rj = job[3];
                    int k = job[4];
                    int rk = job[5];
                    if ((check(i, ri) || check(j, rj) || check(k, rk)
                            || check(i, ri, j, rj) || check(i, ri, k, rk) || check(j, rj, k, rk))) {
                        continue;
                    }

                    Residue residueI = residues[i];
                    Residue residueJ = residues[j];
                    Residue residueK = residues[k];

                    int indexI = allResiduesList.indexOf(residueI);
                    int indexJ = allResiduesList.indexOf(residueJ);
                    int indexK = allResiduesList.indexOf(residueK);

                    double rawDist = getRawNBodyDistance(indexI, ri, indexJ, rj, indexK, rk);
                    double dIJ = checkDistMatrix(indexI, ri, indexJ, rj);
                    double dIK = checkDistMatrix(indexI, ri, indexK, rk);
                    double dJK = checkDistMatrix(indexJ, rj, indexK, rk);
                    double minDist = Math.min(Math.min(dIJ, dIK), dJK);
                    
                    double resDist = get3BodyResidueDistance(indexI, ri, indexJ, rj, indexK, rk);
                    String resDistString = format("     large");
                    if (resDist < Double.MAX_VALUE){
                        resDistString = format("%5.3f", resDist);
                    }
                    
                    String distString = format("     large");
                    if (rawDist < Double.MAX_VALUE) {
                        distString = format("%10.3f", rawDist);
                    }

                    double threeBodyEnergy = 0.0;
                    if (minDist < superpositionThreshold) {
                        threeBodyEnergy = Double.NaN;
                        logger.info(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d:\t    NaN      at %13.6f Ang (%s Ang by residue) < %5.3f Ang.",
                                residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, residueK.toFormattedString(false, true), rk, minDist, resDistString, superpositionThreshold));
                    } else if (checkTriDistThreshold(indexI, ri, indexJ, rj, indexK, rk)) {
                        // Set the two-body energy to 0.0 for separation distances larger than the two-body cutoff.
                        threeBodyEnergy = 0.0;
                        time += System.nanoTime();
                        logger.fine(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, residueK.toFormattedString(false, true), rk,
                                formatEnergy(threeBodyEnergy), distString, resDistString, time * 1.0e-9));
                    } else {
                        try {
                            threeBodyEnergy = compute3BodyEnergy(residues, i, ri, j, rj, k, rk);
                            time += System.nanoTime();
                            logger.info(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                    residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, residueK.toFormattedString(false, true), rk,
                                    formatEnergy(threeBodyEnergy), distString, resDistString, time * 1.0e-9));
                        } catch (ArithmeticException ex) {
                            threeBodyEnergy = Double.NaN;
                            time += System.nanoTime();
                            logger.info(format(" 3-Body %8s %-2d, %8s %-2d, %8s %-2d:\t    NaN      at %s Ang (%s Ang by residue) in %6.4f (sec).",
                                    residueI.toFormattedString(false, true), ri, residueJ.toFormattedString(false, true), rj, residueK.toFormattedString(false, true), rk, distString, resDistString, time * 1.0e-9));
                        }
                    }

                    double commEnergy[] = {i, ri, j, rj, k, rk, threeBodyEnergy};
                    DoubleBuf commEnergyBuf = DoubleBuf.buffer(commEnergy);
                    multicastBuf(commEnergyBuf);
                }
            }
        }
    }

    /**
     * Compute 4-Body energies. This code is experimental.
     */
    private class FourBodyEnergyRegion extends WorkerRegion {

        private final QuadsEnergyLoop energyLoop;
        private final Residue residues[];
        Set<Integer> keySet;

        public FourBodyEnergyRegion(Residue residues[]) {
            energyLoop = new QuadsEnergyLoop();
            this.residues = residues;
        }

        @Override
        public void start() {
            keySet = fourBodyEnergyMap.keySet();
        }

        @Override
        public void run() throws Exception {
            if (!keySet.isEmpty()) {
                execute(0, keySet.size() - 1, energyLoop);
            }
        }

        private class QuadsEnergyLoop extends WorkerIntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.guided();
            }

            @Override
            public void run(int lb, int ub) {
                for (int key = lb; key <= ub; key++) {
                    long time = -System.nanoTime();
                    if (!fourBodyEnergyMap.keySet().contains(key)) {
                        continue;
                    }

                    Integer job[] = fourBodyEnergyMap.get(key);
                    int i = job[0];
                    int ri = job[1];
                    int j = job[2];
                    int rj = job[3];
                    int k = job[4];
                    int rk = job[5];
                    int l = job[6];
                    int rl = job[7];

                    if (check(i, ri) || check(j, rj) || check(k, rk) || check(l, rl)
                            || check(i, ri, j, rj) || check(i, ri, k, rk)
                            || check(i, ri, l, rl) || check(j, rj, k, rk)
                            || check(j, rj, l, rl) || check(k, rk, l, rl)) {
                        // Not implemented: 3-body or 4-body checks.
                        continue;
                    }

                    Residue resi = residues[i];
                    Residue resj = residues[j];
                    Residue resk = residues[k];
                    Residue resl = residues[l];

                    int indexI = allResiduesList.indexOf(residues[i]);
                    int indexJ = allResiduesList.indexOf(residues[j]);
                    int indexK = allResiduesList.indexOf(residues[k]);
                    int indexL = allResiduesList.indexOf(residues[l]);

                    double rawDist = getRawNBodyDistance(indexI, ri, indexJ, rj, indexK, rk, indexL, rl);
                    double dIJ = checkDistMatrix(indexI, ri, indexJ, rj);
                    double dIK = checkDistMatrix(indexI, ri, indexK, rk);
                    double dIL = checkDistMatrix(indexI, ri, indexL, rl);
                    double dJK = checkDistMatrix(indexJ, rj, indexK, rk);
                    double dJL = checkDistMatrix(indexJ, rj, indexL, rl);
                    double dKL = checkDistMatrix(indexK, rk, indexL, rl);

                    double minDist = DoubleStream.of(dIJ, dIK, dIL, dJK, dJL, dKL).min().getAsDouble();

                    String distString = format("     large");
                    if (rawDist < Double.MAX_VALUE) {
                        distString = format("%10.3f", rawDist);
                    }
                    
                    double resDist = get4BodyResidueDistance(indexI, ri, indexJ, rj, indexK, rk, indexL, rl);
                    String resDistString = format("     large");
                    if (resDist < Double.MAX_VALUE){
                        resDistString = format("%5.3f", resDist);
                    }

                    double fourBodyEnergy = 0.0;
                    if (minDist < superpositionThreshold) {
                        fourBodyEnergy = Double.NaN;
                        logger.info(format(" Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d:   set to NaN at %13.6f Ang (%s Ang by residue)  < %5.3f Ang.",
                                residues[i], ri, residues[j].toFormattedString(false, true), rj, residues[k].toFormattedString(false, true), rk, residues[l].toFormattedString(false, true), rl, minDist, resDistString, superpositionThreshold));
                    } //else if (quadCutoff && (dist < quadCutoffDist)) {
                    else if (checkQuadDistThreshold(indexI, ri, indexJ, rj, indexK, rk, indexL, rl)) {
                        // Set the 4-body energy to 0.0 for separation distances larger than the 4-body cutoff.
                        fourBodyEnergy = 0.0;
                        time += System.nanoTime();
                        logger.info(format(" Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).",
                                resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk, resl.toFormattedString(false, true), rl, formatEnergy(fourBodyEnergy), distString, resDistString));
                    } else {
                        try {
                            fourBodyEnergy = compute4BodyEnergy(residues, i, ri, j, rj, k, rk, l, rl);
                            time += System.nanoTime();
                            logger.info(format(" Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).",
                                    resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk, resl.toFormattedString(false, true), rl, formatEnergy(fourBodyEnergy), distString, resDistString));
                            if (abs(fourBodyEnergy) > 1.0) {
                                StringBuilder sb = new StringBuilder();
                                sb.append(format(" Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).\n",
                                        resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk, resl.toFormattedString(false, true), rl, formatEnergy(fourBodyEnergy), distString, resDistString));
                                sb.append(format("   Explain: (ref %d) \n", key));
                                sb.append(format("     Self %3d %3d:                  %.3f\n", i, ri, getSelf(i, ri)));
                                sb.append(format("     Self %3d %3d:                  %.3f\n", j, rj, getSelf(j, rj)));
                                sb.append(format("     Self %3d %3d:                  %.3f\n", k, rk, getSelf(k, rk)));
                                sb.append(format("     Self %3d %3d:                  %.3f\n", l, rl, getSelf(l, rl)));
                                sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, j, rj, get2Body(i, ri, j, rj)));
                                sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, k, rk, get2Body(i, ri, k, rk)));
                                sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", i, ri, l, rl, get2Body(i, ri, l, rl)));
                                sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", j, rj, k, rk, get2Body(j, rj, k, rk)));
                                sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", j, rj, l, rl, get2Body(j, rj, l, rl)));
                                sb.append(format("     Pair %3d %3d %3d %3d:          %.3f\n", k, rk, l, rl, get2Body(k, rk, l, rl)));
                                sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, j, rj, k, rk, get3Body(residues, i, ri, j, rj, k, rk)));
                                sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, j, rj, l, rl, get3Body(residues, i, ri, j, rj, l, rl)));
                                sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", i, ri, k, rk, l, rl, get3Body(residues, i, ri, k, rk, l, rl)));
                                sb.append(format("     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n", j, rj, k, rk, l, rl, get3Body(residues, j, rj, k, rk, l, rl)));
                                sb.append(format("     backbone:                      %.3f\n", backboneEnergy));
                                sb.append(format("     quadEnergy:                 %.3f\n", fourBodyEnergy));
                                sb.append(format("     --s--\n"));
                                sb.append(format("     Active residues:\n"));
                                for (int debug = 0; debug < residues.length; debug++) {
                                    if (residues[debug].getSideChainAtoms().get(0).getUse()) {
                                        sb.append(format("       %s\n", residues[debug].toString()));
                                    }
                                }
                                sb.append(format("     --f--\n"));
                                logger.info(sb.toString());
                            }
                        } catch (ArithmeticException ex) {
                            fourBodyEnergy = Double.NaN;
                            time += System.nanoTime();
                            logger.info(format(" Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: NaN at %s Ang (%s Ang by residue).",
                                    resi.toFormattedString(false, true), ri, resj.toFormattedString(false, true), rj, resk.toFormattedString(false, true), rk, resl.toFormattedString(false, true), rl, distString, resDistString));
                        }
                    }
                }
            }
        }
    }

    /**
     * Compute the minimum distance between each pair of residues for all
     * rotamer permutations.
     */
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

        private double[][] getCoordinates(Residue residue, Rotamer rotamer, boolean forced) {
            synchronized (residue) {
                if (!forced) {
                    RotamerLibrary.applyRotamer(residue, rotamer);
                    return residue.storeCoordinateArray();
                } else {
                    return residue.storeCoordinateArray();
                }
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
                            double xi[][] = getCoordinates(residuei, rotamersi[ri], forcedResidueI);

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
                                    double xj[][] = getCoordinates(residuej, rotamersj[rj], forcedResidueJ);

                                    if (getThreadIndex() == 0 && algorithmListener != null) {
                                        algorithmListener.algorithmUpdate(molecularAssembly);
                                    }

                                    double r = interResidueDistance(xi, xj, symOp);
                                    if (i < j) {
                                        if (r < distanceMatrix[i][ri][j][rj]) {
                                            distanceMatrix[i][ri][j][rj] = r;
                                        }
                                    } else if (r < distanceMatrix[j][rj][i][ri]) {
                                        distanceMatrix[j][rj][i][ri] = r;
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
            return "Rotamer moves utilizing a rotamer energy matrix";
        }
    }

    /**
     * Class contains an integer array, comparing and hashing based on the
     * integer values.
     */
    private class IntegerKeyset {

        private final int[] keys;

        public IntegerKeyset(int... keySet) {
            this.keys = Arrays.copyOf(keySet, keySet.length);
        }

        public int getSize() {
            return keys.length;
        }

        public int[] getKeys() {
            return Arrays.copyOf(keys, keys.length);
        }

        /**
         * Equality based on shallow value equality.
         *
         * @param other
         * @return
         */
        @Override
        public boolean equals(Object other) {
            if (other instanceof IntegerKeyset) {
                return Arrays.equals(keys, ((IntegerKeyset) other).getKeys());
            }
            return false;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(keys);
        }
    }

    /**
     * Writes eliminated singles and pairs to a CSV file. Reads in .log file.
     */
    public void outputEliminated() throws FileNotFoundException, IOException {
        /**
         * eliminated.csv stores the eliminated singles and pairs.
         */
        File fileName = new File("eliminated.csv");
        String path = fileName.getCanonicalPath();
        File outputTxt = new File(path);
        Path currentRelativePath = Paths.get("eliminated.csv");
        String logPath = currentRelativePath.toAbsolutePath().toString();
        logPath = logPath.replaceAll("/eliminated.csv", "");
        File logDirectory = new File(logPath);
        /**
         * Searches for .log file within the directory.
         */
        File[] logArray = logDirectory.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String filename) {
                return filename.endsWith(".log");
            }
        });
        File logFile = logArray[0];
        PrintWriter out = new PrintWriter(new FileWriter(outputTxt, true));
        BufferedReader br = new BufferedReader(new FileReader(logFile));
        String line;
        String line3;
        List<String> boxes = new ArrayList<String>();
        while ((line3 = br.readLine()) != null) {
            if (line3.contains("-a, 5")) {
                /**
                 * Stores the number of box optimization iterations and the
                 * coordinates of each box.
                 */
                while ((line = br.readLine()) != null) {
                    if (line.contains("xyz indices")) {
                        String coordinates = line.replaceAll(" Cell xyz indices:", "");
                        String nextLine = br.readLine();
                        if (!nextLine.contains("Empty box")) {
                            boxes.add(coordinates);
                            String nextLine2;
                            String nextLine4;
                            String nextLine5;
                            while (!(nextLine2 = br.readLine()).contains("Collecting Permutations")) {
                                if (nextLine2.contains("Applying Rotamer Pair")) {
                                    int total = 0;
                                    int total2 = 0;
                                    nextLine2 = br.readLine();
                                    nextLine4 = br.readLine();
                                    nextLine5 = br.readLine();
                                    if (nextLine5.contains("Self-consistent")) {
                                        String[] split = nextLine2.split(" ");
                                        int singlesAmount = Integer.parseInt(split[1]);
                                        if (singlesIterations.size() > 0) {
                                            total = singlesIterations.get(singlesIterations.size() - 1);
                                        }
                                        singlesIterations.add(singlesAmount + total);
                                        System.out.println(nextLine4);
                                        String[] split2 = nextLine4.split(" ");
                                        int pairsAmount = Integer.parseInt(split2[1]);
                                        if (pairsIterations.size() > 0) {
                                            total2 = pairsIterations.get(pairsIterations.size() - 1);
                                        }
                                        pairsIterations.add(pairsAmount + total2);
                                    }
                                }
                            }
                        }
                    }
                }
                /**
                 * Writes eliminated singles and pairs to eliminated.csv, sorted
                 * by box iteration.
                 */
                out.append("Eliminated Singles\n");
                out.append("Residue,Rotamer\n");
                out.append(boxes.get(0) + "\n");
                int stopper = 0;
                for (int i = 0; i < eliminatedResidue.size(); i++) {
                    while (stopper == 0) {
                        for (int x = 0; x < (singlesIterations.size() - 1); x++) {
                            if (singlesIterations.get(x) == 0) {
                                String boxHeader = (boxes.get(x + 1));
                                out.append(boxHeader + System.lineSeparator());
                            }
                            stopper++;
                        }
                    }
                    String first = eliminatedResidue.get(i).toString();
                    String second = eliminatedRotamer.get(i).toString();
                    String resrot = first + "," + second;
                    out.append(resrot + System.lineSeparator());
                    for (int x = 0; x < (singlesIterations.size() - 1); x++) {
                        if ((i + 1) == singlesIterations.get(x)) {
                            String boxHeader = (boxes.get(x + 1));
                            out.append(boxHeader + System.lineSeparator());
                        }
                    }
                }
                out.append("Eliminated Pairs\n");
                out.append("Residue1,Rotamer1,Residue2,Rotamer2" + System.lineSeparator());
                out.append(boxes.get(0) + "\n");
                int stopper2 = 0;
                for (int i = 0; i < eliminatedResidue1.size(); i++) {
                    while (stopper2 == 0) {
                        for (int x = 0; x < (pairsIterations.size() - 1); x++) {
                            if (pairsIterations.get(x) == 0) {
                                String boxHeader = (boxes.get(x + 1));
                                out.append(boxHeader + System.lineSeparator());
                            }
                            stopper2++;
                        }
                    }
                    String first = eliminatedResidue1.get(i).toString();
                    String second = eliminatedRotamer1.get(i).toString();
                    String third = eliminatedResidue2.get(i).toString();
                    String fourth = eliminatedRotamer2.get(i).toString();
                    String resrot = first + "," + second + "," + third + "," + fourth;
                    out.append(resrot + System.lineSeparator());
                    for (int x = 0; x < (pairsIterations.size() - 1); x++) {
                        if ((i + 1) == pairsIterations.get(x)) {
                            String boxHeader = (boxes.get(x + 1));
                            out.append(boxHeader + System.lineSeparator());
                        }
                    }
                }
                /**
                 * Writes eliminated singles and pairs to eliminated.csv for
                 * global optimization.
                 */
            } else if ((line3.contains("-a, 1")) || (line3.contains("-a, 2")) || (line3.contains("-a, 3")) || (line3.contains("-a, 4"))) {
                out.append("Eliminated Singles\n");
                out.append("Residue,Rotamer\n");
                for (int i = 0; i < eliminatedResidue.size(); i++) {
                    String first = eliminatedResidue.get(i).toString();
                    String second = eliminatedRotamer.get(i).toString();
                    String resrot = first + "," + second;
                    out.append(resrot + System.lineSeparator());
                }
                out.append("Eliminated Pairs\n");
                out.append("Residue1,Rotamer1,Residue2,Rotamer2" + System.lineSeparator());
                for (int i = 0; i < eliminatedResidue1.size(); i++) {
                    String first = eliminatedResidue1.get(i).toString();
                    String second = eliminatedRotamer1.get(i).toString();
                    String third = eliminatedResidue2.get(i).toString();
                    String fourth = eliminatedRotamer2.get(i).toString();
                    String resrot = first + "," + second + "," + third + "," + fourth;
                    out.append(resrot + System.lineSeparator());
                }
            }
        }
        out.close();
        br.close();
    }
}
