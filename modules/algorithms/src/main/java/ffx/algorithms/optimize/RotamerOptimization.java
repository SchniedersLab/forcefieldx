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
package ffx.algorithms.optimize;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.abs;

import edu.rit.pj.Comm;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.WorkerTeam;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.algorithms.mc.MCMove;
import ffx.algorithms.optimize.manybody.BoxOptCell;
import ffx.algorithms.optimize.manybody.DistanceMatrix;
import ffx.algorithms.optimize.manybody.EliminatedRotamers;
import ffx.algorithms.optimize.manybody.EnergyExpansion;
import ffx.algorithms.optimize.manybody.EnergyRegion;
import ffx.algorithms.optimize.manybody.FourBodyEnergyRegion;
import ffx.algorithms.optimize.manybody.GoldsteinPairRegion;
import ffx.algorithms.optimize.manybody.RotamerMatrixMC;
import ffx.algorithms.optimize.manybody.RotamerMatrixMove;
import ffx.algorithms.optimize.manybody.SelfEnergyRegion;
import ffx.algorithms.optimize.manybody.ThreeBodyEnergyRegion;
import ffx.algorithms.optimize.manybody.TwoBodyEnergyRegion;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.Potential;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parsers.PDBFilter;
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
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class RotamerOptimization implements Terminatable {

    /**
     * Logger for this class.
     */
    private static final Logger logger = Logger.getLogger(RotamerOptimization.class.getName());

    public enum Algorithm {

        ALL, BOX, WINDOW, INDEPENDENT, BRUTE_FORCE
    }

    /**
     * Allows get2BodyDistance to find the shortest distance between any two
     * rotamers or two residues.
     */
    public enum DistanceMethod {

        ROTAMER, RESIDUE
    }

    public enum Direction {

        FORWARD, BACKWARD
    }

    /**
     * RotamerLibrary instance.
     */
    protected RotamerLibrary library = RotamerLibrary.getDefaultLibrary();
    private DistanceMatrix dM;
    private EnergyExpansion eE;
    private EliminatedRotamers eR;

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
    private final List<MolecularAssembly> allAssemblies;
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
     * Two-Body cutoff distance.
     */
    private double twoBodyCutoffDist;
    /**
     * Fallback if there is no vdW node.
     */
    private static final double FALLBACK_TWO_BODY_CUTOFF = 0;
    /**
     * Flag to control use of 3-body terms.
     */
    private boolean threeBodyTerm = false;
    /**
     * Three-body cutoff distance.
     */
    private double threeBodyCutoffDist;
    /**
     * Interaction partners of a Residue that come after it.
     */
    private int[][] resNeighbors;
    /**
     * All interaction partners of a Residue, including prior residues.
     */
    private int[][] bidiResNeighbors;
    /**
     * Flag to prune clashes.
     */
    private boolean pruneClashes = true;
    /**
     * Flag to prune pair clashes.
     */
    private boolean prunePairClashes = true;
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
    private boolean revert;
    /**
     * The distance that the distance matrix checks for.
     */
    private double distance = 2.0;
    /**
     * Default distance method is to find the shortest distance between
     * residues.
     */
    private DistanceMethod distanceMethod = DistanceMethod.RESIDUE;
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
     * Flag to load the distance matrix as needed; if false, matrix is prefilled
     * at the beginning of rotamer optimization.
     */
    private boolean lazyMatrix = false;
    /**
     * Clash energy threshold (kcal/mole).
     */
    private double clashThreshold = 25.0;
    /**
     * Clash energy threshold (kcal/mol) for MultiResidues, which can have much
     * more variation in self and 2-Body energies.
     */
    private double multiResClashThreshold = 80.0;
    /**
     * Clash energy threshold (kcal/mole).
     */
    private double pairClashThreshold = 25.0;
    /**
     * Pair clash energy threshold (kcal/mol) for MultiResidues.
     */
    private double multiResPairClashAddn = 80.0;
    /**
     * False unless ManyBodyTest is occurring.
     */
    private boolean monteCarloTesting = false;
    /**
     * An array of atomic coordinates (length 3 * the number of atoms).
     */
    private double[] x = null;
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
     * The approximate energy from summing the backbone energy, self-energies, pair-energies, etc.
     */
    private double approximateEnergy = 0;
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
    private Residue[] allResiduesArray = null;
    /**
     * Number of residues being optimized.
     */
    private int numResidues = 0;
    /**
     * The array of optimum rotamers for each residue.
     */
    private int[] optimum;
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
     * ParallelTeam instance.
     */
    private ParallelTeam parallelTeam;
    /**
     * Flag to indicate computation of 4-body energy values. This is limited to
     * the study 4-body energy magnitudes, but is not included in the rotamer
     * optimization.
     */
    private boolean compute4BodyEnergy = false;
    /**
     * Parallel evaluation of quantities used during Goldstein Pair elimination.
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
    private double superpositionThreshold = 0.25;
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
    private boolean addOrigRot = true;
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
     * Maximum depth to check if a rotamer can be eliminated.
     */
    private int maxRotCheckDepth;
    /**
     * Writes energies to restart file.
     */
    private BufferedWriter energyWriter;

    /**
     * RotamerOptimization constructor.
     *
     * @param molecularAssembly The MolecularAssembly to search rotamers for.
     * @param potential         a {@link ffx.numerics.Potential} object.
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
                logger.info(format("Warning: threeBodyCuoffDist should not be less than 0."));
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
        int defaultMaxRotCheckDepth = 1;
        if (propStr != null) {
            try {
                maxRotCheckDepth = Integer.parseInt(propStr);
                if (maxRotCheckDepth > 3 || maxRotCheckDepth < 0) {
                    throw new IllegalArgumentException(" ro-maxRotSearchDepth must be between 0-3 inclusive!");
                }
            } catch (Exception ex) {
                maxRotCheckDepth = defaultMaxRotCheckDepth;
                logger.warning(format(" Could not parse %s value %s as valid integer; defaulting to %d", prop, propStr, maxRotCheckDepth));
                logger.warning(format(" Exception: %s", ex));
            }
        } else {
            maxRotCheckDepth = defaultMaxRotCheckDepth;
        }

        prop = System.getProperty("revertUnfavorable", "false");
        revert = Boolean.parseBoolean(prop);

        allAssemblies = new ArrayList<>();
        allAssemblies.add(molecularAssembly);
        setUpRestart();
    }

    public EnergyExpansion getEnergyExpansion() {
        return eE;
    }

    public EliminatedRotamers getEliminatedRotamers() {
        return eR;
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
     * @param maxRotCheckDepth a int.
     */
    public void setMaxRotCheckDepth(int maxRotCheckDepth) {
        this.maxRotCheckDepth = maxRotCheckDepth;
    }

    /**
     * A brute-force global optimization over side-chain rotamers using a
     * recursive algorithm.
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     *                          object.
     * @param residues          an array of {@link ffx.potential.bonded.Residue} objects.
     * @param i                 a int.
     * @param lowEnergy         a double.
     * @param optimum           an array of {@link int} objects.
     * @return the current energy.
     */
    public double rotamerOptimization(MolecularAssembly molecularAssembly, Residue[] residues, int i,
                                      double lowEnergy, int[] optimum) {

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
                    logger.info(format(" %d Energy set to NaN (unreasonable conformation)", ++evaluatedPermutations));
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
    private double decomposedRotamerOptimization(MolecularAssembly molecularAssembly, Residue[] residues, int i,
                                                 double lowEnergy, int[] optimum, int[] currentRotamers) {

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
                if (eR.check(i, ri)) {
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
                if (eR.check(i, ri)) {
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
                    arraycopy(currentRotamers, 0, optimum, 0, optimum.length);
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
        arraycopy(initialRots, 0, optimum, 0, nRes);
        assert optimum.length == nRes;
        assert initialRots.length == nRes;

        RotamerMatrixMC rmc = new RotamerMatrixMC(initialRots, residues, useFullAMOEBAEnergy, this);
        rmc.setTemperature(mcTemp);
        RotamerMatrixMove rmove = new RotamerMatrixMove(useAllElims, initialRots, residues,
                library, this, eR, monteCarloTesting);
        List<MCMove> rmList = new ArrayList<>(1);
        rmList.add(rmove);

        double initialEnergy = computeEnergy(residues, initialRots, false);
        double optimumEnergy = initialEnergy;
        double currentEnergy = initialEnergy;

        int nAccept = 0;
        logIfMaster(format(" Beginning %d iterations of Monte Carlo search "
                + "starting from energy %10.6f", maxIters, initialEnergy));

        for (int i = 0; i < maxIters; i++) {
            if (rmc.mcStep(rmList, currentEnergy)) {
                currentEnergy = rmc.lastEnergy();
                ++nAccept;
                if (currentEnergy < optimumEnergy) {
                    optimumEnergy = currentEnergy;
                    arraycopy(initialRots, 0, optimum, 0, nRes);
                }
            } else {
                currentEnergy = rmc.lastEnergy();
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
                if (!eR.check(i, ri)) {
                    allowedRots.add(ri);
                }
            }

            Random rand = new Random();
            int nRots = allowedRots.size();
            if (monteCarloTesting) {
                rand.setSeed(nRots);
            }
            if (nRots > 1) {
                boolean validMove = !useAllElims;
                int indexRI;
                do {
                    int ri = rand.nextInt(nRots);
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
     * @param i  Residue number
     * @param ri Rotamer number
     * @return If valid
     */
    public boolean checkValidMove(int i, int ri, int[] currentRots) {
        int nRes = currentRots.length;
        for (int j = 0; j < nRes; j++) {
            if (j == i) {
                continue;
            }
            if (eR.check(j, currentRots[j], i, ri)) {
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
    private double rotamerOptimizationDEE(MolecularAssembly molecularAssembly, Residue[] residues, int i,
                                          int[] currentRotamers, double lowEnergy, int[] optimum,
                                          double[] permutationEnergies) {
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
                if (eR.check(i, ri)) {
                    continue;
                }
                /**
                 * Check if rotamer ri has been eliminated by an upstream
                 * rotamer (any residue's rotamer from j = 0 .. i-1).
                 */
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
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
                if (eR.check(i, ri)) {
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
                    deadEnd = eR.check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                applyRotamer(residuei, rotamersi[ri]);
                // Compute the energy based on a 3-body approximation
                approximateEnergy = computeEnergy(residues, currentRotamers, false);
                double comparisonEnergy = approximateEnergy;
                evaluatedPermutations++;
                // Compute the AMOEBA energy
                if (useFullAMOEBAEnergy) {
                    double amoebaEnergy = Double.NaN;
                    try {
                        amoebaEnergy = currentEnergy(resList);
                    } catch (ArithmeticException ex) {
                        logger.warning(format(" Exception %s in calculating full AMOEBA energy for permutation %d", ex.toString(), evaluatedPermutations));
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
                        ;
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
    private double dryRun(Residue[] residues, int i, int[] currentRotamers) {
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
                if (eR.check(i, ri)) {
                    continue;
                }
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
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
                if (eR.check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
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
    private double dryRunForEnsemble(Residue[] residues, int i, int[] currentRotamers,
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
                if (eR.check(i, ri)) {
                    continue;
                }
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
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
                if (eR.check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (deadEnd) {
                    continue;
                }
                if (permutationEnergies[evaluatedPermutations] - gmecEnergy < ensembleEnergy) {
                    permutations[evaluatedPermutations] = new int[nResidues];
                    arraycopy(currentRotamers, 0, permutations[evaluatedPermutations], 0, nResidues);
                }
                evaluatedPermutations++;
            }
        }
        return 0.0;
    }

    /**
     * Flag to control use of 3-body energy terms.
     *
     * @param threeBodyTerm a boolean.
     */
    public void setThreeBodyEnergy(boolean threeBodyTerm) {
        this.threeBodyTerm = threeBodyTerm;
    }

    /**
     * Uses existing backbone, self, 2-Body, and 3-body energies from
     * rotamerEnergies() to calculate an approximate energy for a rotamer
     * permutation.
     *
     * @param residues Current window of optimization.
     * @param rotamers Set of rotamers to calculate an approximate energy for.
     * @param print    Verbosity flag
     * @return Approximate permutation energy (backbone + selfs + pairs +
     * trimers).
     */
    public double computeEnergy(Residue[] residues, int[] rotamers, boolean print) {
        double selfSum;
        double pairSum;
        double threeBodySum;
        try {
            if (parallelTeam == null) {
                parallelTeam = new ParallelTeam();
            }
            if (energyRegion == null) {
                energyRegion = new EnergyRegion(parallelTeam.getThreadCount());
            }
            energyRegion.init(eE, residues, rotamers, threeBodyTerm);
            parallelTeam.execute(energyRegion);
            selfSum = energyRegion.getSelf();
            pairSum = energyRegion.getTwoBody();
            threeBodySum = energyRegion.getThreeBody();
        } catch (Exception e) {
            throw new IllegalArgumentException(e);
        }

        approximateEnergy = eE.getBackboneEnergy() + selfSum + pairSum + threeBodySum;
        if (print) {
            logger.info(format(" Backbone                  %s", formatEnergy(eE.getBackboneEnergy())));
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
     *                          without optimizing.
     */
    public void setDecomposeOriginal(boolean decomposeOriginal) {
        this.decomposeOriginal = decomposeOriginal;
    }

    /**
     * Utility method for formatting energies, using 16 spaces with 8 digits of
     * precision.
     *
     * @param energy Energy to format.
     * @return A string representing the energy.
     */
    public String formatEnergy(double energy) {
        if (abs(energy) < 1.0e6) {
            return format("%16.8f", energy);
        } else {
            return format("*%15.4e", energy);
        }
        // TODO: Possibly replace with %16.8g, which so far as I know, is equivalent.
    }

    /**
     * Set a contiguous block of residues to optimize in a specific chain.
     *
     * @param chain      a {@link java.lang.String} object.
     * @param startResID a int.
     * @param finalResID a int.
     */
    public void setResidues(String chain, int startResID, int finalResID) {
        this.chain = chain;
        this.setResidues(startResID, finalResID);
    }

    /**
     * Set a contiguous block of residues to optimize.
     *
     * @param startResID a int.
     * @param finalResID a int.
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
     * @param i       residues index in polymer.
     */
    private void initResidue(Residue residue, Polymer polymer, int i) {
        if (residue == null) {
            logger.warning(format(" Null residue %d for chain %c", i, polymer.getChainID()));
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
     *
     * @return energyRestartFile File with saved side-chain energies.
     */
    public File getRestartFile() {
        return energyRestartFile;
    }

    /**
     * Return the residue list.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Residue> getResidues() {
        return residueList;
    }

    /**
     * Set the residue list.
     *
     * @param residueList a {@link java.util.ArrayList} object.
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
    int[] getOptimumRotamers() {
        return optimum;
    }

    public double getApproximate(){
        return approximateEnergy;
    }
    /**
     * <p>
     * setCoordinatesToEnsemble.</p>
     *
     * @param ensnum a int.
     */
    public void setCoordinatesToEnsemble(int ensnum) {
        if (ensembleStates != null && !ensembleStates.isEmpty()) {
            ensnum %= ensembleStates.size();
            ResidueState.revertAllCoordinates(residueList, ensembleStates.get(ensnum).getVal());
        } else {
            throw new IllegalArgumentException(" Ensemble states not initialized!");
        }
    }

    /**
     * <p>
     * getEnsemble.</p>
     *
     * @return a {@link java.util.List} object.
     */
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
     * <p>
     * setEnsemble.</p>
     *
     * @param ensemble a int.
     */
    public void setEnsemble(int ensemble) {
        setEnsemble(ensemble, 5.0);
    }

    /**
     * Accepts a list of residues but throws out null residues. Used by the -lR
     * flag.
     *
     * @param residueList a {@link java.util.List} object.
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
     * @param algorithm a {@link RotamerOptimization.Algorithm}
     *                  object.
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
        double e = 0;
        try {
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
                dM = new DistanceMatrix(this, molecularAssembly, algorithmListener, allResiduesArray, allResiduesList,
                        library, distanceMethod, distance, twoBodyCutoffDist, threeBodyCutoffDist, lazyMatrix, useForcedResidues);
            }

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
        } finally {
            try {
                if (energyWriter != null) {
                    energyWriter.close();
                }
            } catch (IOException ex) {
                logger.severe(format(" Exception in closing buffered energy writer."));
            }
        }

        return e;
    }

    /**
     * Specify use of Goldstein optimization.
     *
     * @param useGoldstein a boolean.
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
        arraycopy(numXYZBoxes, 0, this.numXYZBoxes, 0, this.numXYZBoxes.length);
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
     * @param boxStart a int.
     */
    public void setBoxStart(int boxStart) {
        this.boxStart = boxStart;
    }

    /**
     * Set the ending box index.
     *
     * @param boxEnd a int.
     */
    public void setBoxEnd(int boxEnd) {
        // Is -1 if boxes run to completion.
        this.boxEnd = boxEnd;
    }

    /**
     * Set the optimization direction to forward or backward.
     *
     * @param direction a {@link RotamerOptimization.Direction}
     *                  object.
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
    }

    /**
     * Set the cut-off distance for inclusion of residues in sliding box and
     * window methods.
     *
     * @param distance a double.
     */
    public void setDistanceCutoff(double distance) {
        this.distance = distance;
    }

    /**
     * Set the residue increment for sliding window.
     *
     * @param increment a int.
     */
    public void setIncrement(int increment) {
        this.increment = increment;
    }

    /**
     * Set the algorithm to revert to starting coordinates if the energy
     * increases.
     *
     * @param revert a boolean.
     */
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

    /**
     * The the nucleic acid correction threshold.
     *
     * @param nucleicCorrectionThreshold a double.
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
     * @param minNumberAcceptedNARotamers a int.
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
     * @param nucleicPruningFactor a double.
     */
    public void setNucleicPruningFactor(double nucleicPruningFactor) {
        this.nucleicPruningFactor = nucleicPruningFactor;
        this.nucleicPairsPruningFactor = ((1.0 + nucleicPruningFactor) / 2);
    }

    /**
     * Sets whether rotamer optimization should print out any files, or act
     * solely to optimize a structure in memory.
     *
     * @param printFiles a boolean.
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
                double newE = Double.NaN;
                try {
                    newE = currentEnergy(rList);
                } catch (ArithmeticException ex) {
                    logger.fine(format(" Exception %s in energy calculations during independent for %s-%d", ex.toString(), residue, j));
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
    private boolean firstValidPerm(Residue[] residues, int i, int[] currentRotamers) {
        int nResidues = residues.length;
        Residue residuei = residues[i];
        Rotamer[] rotamersi = residuei.getRotamers(library);
        int lenri = rotamersi.length;
        if (i < nResidues - 1) {
            for (int ri = 0; ri < lenri; ri++) {
                if (eR.check(i, ri)) {
                    continue;
                }
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
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
                if (eR.check(i, ri)) {
                    continue;
                }
                currentRotamers[i] = ri;
                boolean deadEnd = false;
                for (int j = 0; j < i; j++) {
                    int rj = currentRotamers[j];
                    deadEnd = eR.check(j, rj, i, ri);
                    if (deadEnd) {
                        break;
                    }
                }
                if (!deadEnd) {
                    break;
                }
            }
            return true;
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
        Residue[] residues = residueList.toArray(new Residue[residueList.size()]);
        int nResidues = residues.length;
        int[] currentRotamers = new int[nResidues];

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
                        if (!eR.eliminatedSingles[i][ri]) {
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
            logIfMaster(format("%30s %35s %35s", "Condition", "Number of Permutations Left", "Number of Permutations Removed"));
            logIfMaster(format("%30s %35s %35s", "No Eliminations", permutations, ""));
            logIfMaster(format("%30s %35s %35s", "Single Eliminations", singletonPermutations, permutations - singletonPermutations));
            logIfMaster(format("%30s %35s %35s", "Pair Eliminations", afterPairElim, pairTotalElimination));
            logIfMaster(format("%30s %35s %35s", "Single and Pair Eliminations", (double) evaluatedPermutations, pairTotalElimination + (permutations - singletonPermutations)));

            logIfMaster(format("\n Energy of permutations:"));
            logIfMaster(format(" %12s %25s %25s", "Permutation", "Energy", "Lowest Possible Energy"));

            double e;
            if (useMonteCarlo()) {
                firstValidPerm(residues, 0, currentRotamers);
                arraycopy(currentRotamers, 0, optimum, 0, nResidues);
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
                sumSelfEnergy += eE.getSelf(i, ri);
                logIfMaster(format(" Final self Energy (%8s,%2d): %12.4f",
                        residues[i].toFormattedString(false, true), ri, eE.getSelf(i, ri)));
            }
            for (int i = 0; i < nResidues - 1; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues; j++) {
                    int rj = optimum[j];
                    sumPairEnergy += eE.get2Body(i, ri, j, rj);
                    if (eE.get2Body(i, ri, j, rj) > 10.0) {
                        logIfMaster(format(" Large Final Pair Energy (%8s,%2d) (%8s,%2d): %12.4f",
                                residues[i].toFormattedString(false, true), ri,
                                residues[j].toFormattedString(false, true), rj,
                                eE.get2Body(i, ri, j, rj)));
                    }
                }
            }

            try {
                e = currentEnergy(residueList);
            } catch (ArithmeticException ex) {
                e = Double.NaN;
                logger.severe(format(" Exception %s in calculating current energy at the end of triples", ex.toString()));
            }

            logIfMaster(format(" %12s %25s %25s", "Type", "Energy", "Lowest Possible Energy"));
            logIfMaster(format(" %12s %25f %25s", "Self:", sumSelfEnergy, ""));
            logIfMaster(format(" %12s %25f %25s", "Pair:", sumPairEnergy, ""));

            approximateEnergy = eE.getBackboneEnergy() + sumSelfEnergy + sumPairEnergy;

            if (threeBodyTerm) {
                for (int i = 0; i < nResidues - 2; i++) {
                    int ri = optimum[i];
                    for (int j = i + 1; j < nResidues - 1; j++) {
                        int rj = optimum[j];
                        for (int k = j + 1; k < nResidues; k++) {
                            int rk = optimum[k];
                            try {
                                sumTrimerEnergy += eE.get3Body(residues, i, ri, j, rj, k, rk);
                            } catch (Exception ex) {
                                logger.warning(ex.toString());
                            }
                        }
                    }
                }
                approximateEnergy += sumTrimerEnergy;
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - eE.getBackboneEnergy();
                logIfMaster(format(" %12s %5s %25f %5s %25s %5s", "Trimer:", sumTrimerEnergy, ""));
                logIfMaster(format(" %12s %5s %25f %5s %25s %5s", "Neglected:", higherOrderEnergy, ""));
            } else {
                double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - eE.getBackboneEnergy();
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
                        if (!eR.eliminatedSingles[i][ri]) {
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
            double self = eE.getSelf(i, ri);
            residueEnergy[i] = self;
            sumSelfEnergy += self;
            double lowest = eE.lowestSelfEnergy(residues, i);
            sumLowSelfEnergy += lowest;
            if (self - lowest > 10.0) {
                logIfMaster(format(" %15s %25f %25f", "Self (" + residues[i] + "," + ri + "):", self, lowest));
            }
        }

        double sumPairEnergy = 0.0;
        double sumLowPairEnergy = 0.0;
        double[] resPairEnergy = new double[nResidues];
        double[] lowPairEnergy = new double[nResidues];
        for (int i = 0; i < nResidues - 1; i++) {
            StringBuilder sb = new StringBuilder();
            int ri = optimum[i];
            double sumPairEnergyI = 0;
            double sumLowPairEnergyI = 0;
            for (int j = i + 1; j < nResidues; j++) {
                int rj = optimum[j];
                double pair = eE.get2Body(i, ri, j, rj);
                residueEnergy[i] += 0.5 * pair;
                residueEnergy[j] += 0.5 * pair;
                sumPairEnergy += pair;
                sumPairEnergyI += pair;
                double lowest = eE.lowestPairEnergy(residues, i, ri, j);
                sumLowPairEnergy += lowest;
                sumLowPairEnergyI += lowest;
                resPairEnergy[i] = 0.5 * pair;
                resPairEnergy[j] = 0.5 * pair;
                lowPairEnergy[i] = 0.5 * lowest;
                lowPairEnergy[j] = 0.5 * lowest;
                if (resPairEnergy[i] - lowPairEnergy[i] > 10.0) {
                    sb.append(format("  Pair Energy (%8s,%2d) (%8s,%2d): %12.4f (Lowest: %12.4f).\n",
                            residues[i].toFormattedString(false, true), ri,
                            residues[j].toFormattedString(false, true), rj, pair, lowest));
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
            logger.severe(format(" Exception %s in calculating current energy at the end of self and pairs", ex.toString()));
        }
        logIfMaster(format(" %15s %25f %25s", "Backbone:", eE.getBackboneEnergy(), ""));
        logIfMaster(format(" %15s %25f %25f", "Self:", sumSelfEnergy, sumLowSelfEnergy));
        logIfMaster(format(" %15s %25f %25f", "Pair:", sumPairEnergy, sumLowPairEnergy));

        approximateEnergy = eE.getBackboneEnergy() + sumSelfEnergy + sumPairEnergy;

        double sumTrimerEnergy = 0;
        if (threeBodyTerm) {
            for (int i = 0; i < nResidues - 2; i++) {
                int ri = optimum[i];
                for (int j = i + 1; j < nResidues - 1; j++) {
                    int rj = optimum[j];
                    for (int k = j + 1; k < nResidues; k++) {
                        int rk = optimum[k];
                        try {
                            double triple = eE.get3Body(residues, i, ri, j, rj, k, rk);
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
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - eE.getBackboneEnergy();
            logIfMaster(format(" %15s %25f %25s", "Trimer:", sumTrimerEnergy, ""));
            logIfMaster(format(" %15s %25f %25s", "Neglected:", higherOrderEnergy, ""));
        } else {
            double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - eE.getBackboneEnergy();
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
                logger.severe(format(" Exception %s in calculating full energy; FFX shutting down", ex.toString()));
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
                    logger.fine(format(" Exception %s in calculating full AMOEBA energy at the end of brute force", ex.toString()));
                }
            }
        }
        return e;
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
    public boolean checkIfForced(Residue residue) throws NullPointerException {
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
                                        double rotamerSeparation = dM.get2BodyDistance(indexI, ri, j, rj);
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
                            logger.severe(format(" Exception %s in calculating starting energy of a window; FFX shutting down", ex.toString()));
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
                                    logger.severe(format(" Exception %s in calculating final energy of a window; FFX shutting down", ex.toString()));
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
                                logger.severe(format(" Exception %s in calculating final energy of a window; FFX shutting down", ex.toString()));
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
                    maxXYZ[i] = Math.max(refAtomCoords[i], maxXYZ[i]);
                    minXYZ[i] = Math.min(refAtomCoords[i], minXYZ[i]);
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
                if (master && writeEnergyRestart && printFiles) {
                    String boxHeader = format(" Box %d: %d,%d,%d", i + 1, cellIndices[0], cellIndices[1], cellIndices[2]);
                    try {
                        energyWriter.append(boxHeader);
                        energyWriter.newLine();
                        boxHeader = null;
                    } catch (IOException ex) {
                        logger.log(Level.SEVERE, " Exception writing box header to energy restart file.", ex);
                    }
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
                        logger.severe(format(" Exception %s in calculating starting energy of a box; FFX shutting down", ex.toString()));
                    }
                    globalOptimization(residuesList);
                    try {
                        finalEnergy = currentEnergy(residuesList);
                    } catch (ArithmeticException ex) {
                        logger.severe(format(" Exception %s in calculating starting energy of a box; FFX shutting down", ex.toString()));
                    }

                    if (startingEnergy <= finalEnergy) {
                        logger.warning("Optimization did not yield a better energy. Reverting to original coordinates.");
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
                    // Don't write a file if it's the final iteration.
                    if (i == (numCells - 1)) {
                        continue;
                    }
                    try {
                        if (firstResidue != lastResidue) {
                            logIfMaster(format(" File with residues %s ... %s in window written.", firstResidue.toString(), lastResidue.toString()));
                        } else {
                            logIfMaster(format(" File with residue %s in window written.", firstResidue.toString()));
                        }
                        firstCellSaved = true;
                    } catch (Exception e) {
                        logger.warning(format("Exception writing to file."));
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
     * @param crystal  Crystal group.
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
                            contained = cell.anyRotamerInsideCell(residue, crystal, symOp, true, library);
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
     * <p>
     * applyEliminationCriteria.</p>
     *
     * @param residues    an array of {@link ffx.potential.bonded.Residue} objects.
     * @param getEnergies a boolean.
     * @param verbose     a boolean.
     */
    private void applyEliminationCriteria(Residue[] residues, boolean getEnergies, boolean verbose) {
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
                eR.validateDEE(residues);
                logIfMaster(toString());
            } while (pairEliminated);
            logIfMaster(" Self-consistent DEE rotamer elimination achieved.\n");
        }
        if (!verbose) {
            logger.setLevel(prevLevel);
        }
    }

    public void logIfMaster(String msg) {
        if (master) {
            logger.info(msg);
        }
    }

    public void logIfMaster(String msg, Level level) {
        if (master) {
            logger.log(level, msg);
        }
    }

    private void applyEliminationCriteria(Residue[] residues) {
        // allocateEliminationMemory is now called for all algorithms in rotamerEnergies method.
        //allocateEliminationMemory(residues);

        if (verboseEnergies) {
            try {
                logIfMaster(format("\n Beginning Energy %s", formatEnergy(currentEnergy(residues))));
            } catch (ArithmeticException ex) {
                logger.severe(format(" Exception %s in calculating beginning energy; FFX shutting down.", ex.toString()));
            }
        }

        rotamerEnergies(residues);

        if (testing) {
            int nres = residues.length;
            eR.onlyPrunedSingles = new boolean[nres][];
            eR.onlyPrunedPairs = new boolean[nres][][][];
            for (int i = 0; i < nres; i++) {
                Residue residuei = residues[i];
                Rotamer rotamersi[] = residuei.getRotamers(library);
                int lenri = rotamersi.length;  // Length rotamers i
                eR.onlyPrunedSingles[i] = new boolean[lenri];
                eR.onlyPrunedSingles[i] = Arrays.copyOf(eR.eliminatedSingles[i], eR.eliminatedSingles[i].length);
                eR.onlyPrunedPairs[i] = new boolean[lenri][][];
                // Loop over the set of rotamers for residue i.
                for (int ri = 0; ri < lenri; ri++) {
                    eR.onlyPrunedPairs[i][ri] = new boolean[nres][];
                    for (int j = i + 1; j < nres; j++) {
                        Residue residuej = residues[j];
                        Rotamer rotamersj[] = residuej.getRotamers(library);
                        int lenrj = rotamersj.length;
                        eR.onlyPrunedPairs[i][ri][j] = new boolean[lenrj];
                        eR.onlyPrunedPairs[i][ri][j] = Arrays.copyOf(eR.eliminatedPairs[i][ri][j], eR.eliminatedPairs[i][ri][j].length);
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
            eR.validateDEE(residues);
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
            eR.validateDEE(residues);
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
    public double currentEnergy(Residue[] resArray) throws ArithmeticException {
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
     * Forces the use of a ForceFieldEnergyOpenMM's underlying ForceFieldEnergy.
     *
     * @return Current potential energy as calculated by FFX reference platform.
     */
    public double currentFFXPE() {
        potential.getCoordinates(x);
        return ((ForceFieldEnergyOpenMM) potential).ffxEnergy(x, false);
    }

    /**
     * Wrapper intended for use with RotamerMatrixMC.
     *
     * @param resList
     * @return
     * @throws ArithmeticException
     */
    public double currentEnergyWrapper(List<Residue> resList) throws ArithmeticException {
        return currentEnergy(resList);
    }

    /**
     * Generates list of subsequent, neighboring Residues.
     *
     * @param residues To generate neighbor lists for.
     */
    private void generateResidueNeighbors(Residue[] residues) {
        int nRes = residues.length;
        resNeighbors = new int[nRes][];
        bidiResNeighbors = new int[nRes][];
        for (int i = 0; i < nRes; i++) {
            Set<Integer> nearby = new HashSet<>();
            Residue resi = residues[i];
            Rotamer[] rotsi = resi.getRotamers(library);
            int lenri = rotsi.length;
            int indexI = allResiduesList.indexOf(resi);

            for (int j = 0; j < nRes; j++) {
                if (i == j) {
                    continue;
                }
                Residue resj = residues[j];
                Rotamer[] rotsj = resj.getRotamers(library);
                int lenrj = rotsj.length;
                int indexJ = allResiduesList.indexOf(resj);

                boolean foundClose = false;
                for (int ri = 0; ri < lenri; ri++) {
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (!dM.checkPairDistThreshold(indexI, ri, indexJ, rj)) {
                            foundClose = true;
                            break;
                        }
                    }
                    if (foundClose) {
                        break;
                    }
                }
                if (foundClose) {
                    nearby.add(j);
                }
            }

            // Collect all neighbors.
            int[] nI = nearby.stream().mapToInt(Integer::intValue).toArray();
            bidiResNeighbors[i] = nI;

            // Collect only subsequent neighbors.
            final int fi = i; // Final copy of i.
            nI = nearby.stream().
                    mapToInt(Integer::intValue).
                    filter(j -> j > fi).
                    toArray();
            resNeighbors[i] = nI;
        }
    }

    private void setUpRestart() {
        File restartFile;
        if (loadEnergyRestart) {
            restartFile = energyRestartFile;
        } else {
            File file = molecularAssembly.getFile();
            String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
            Path restartPath = Paths.get(filename + ".restart");
            restartFile = restartPath.toFile();
            energyRestartFile = restartFile;
        }
        try {
            energyWriter = new BufferedWriter(new FileWriter(restartFile, true));
        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Couldn't open energy restart file.", ex);
        }
        logger.info(format(" Energy restart file: %s", restartFile.getName()));
    }

    /**
     * Turn off non-bonded contributions from all residues except for one.
     * Compute the self-energy for each residue relative to the backbone
     * contribution.
     *
     * @param residues A list of residues that we undergo rotamer optimization.
     * @return template energy
     */
    private double rotamerEnergies(Residue[] residues) {

        if (residues == null) {
            logger.warning(" Attempt to compute rotamer energies for an empty array of residues.");
            return 0.0;
        }

        int nResidues = residues.length;
        Atom[] atoms = molecularAssembly.getAtomArray();
        generateResidueNeighbors(residues);

        eR = new EliminatedRotamers(this, dM, library, allResiduesList, maxRotCheckDepth, clashThreshold,
                pairClashThreshold, multiResClashThreshold, nucleicPruningFactor, nucleicPairsPruningFactor,
                multiResPairClashAddn, pruneClashes, prunePairClashes, print, residues);

        if (decomposeOriginal) {
            assert library.getUsingOrigCoordsRotamer();
            for (int i = 0; i < nResidues; i++) {
                Residue resi = residues[i];
                Rotamer[] rotsi = resi.getRotamers(library);
                int lenri = rotsi.length;

                // Leave the 0'th original-coordinates rotamer alone.
                for (int ri = 1; ri < lenri; ri++) {
                    eR.eliminateRotamer(residues, i, ri, false);
                }
            }
        }

        // Initialize all atoms to be used.
        for (Atom atom : atoms) {
            atom.setUse(true);
        }

        eE = new EnergyExpansion(this, dM, eR, molecularAssembly, potential, library, algorithmListener,
                allResiduesList, resNeighbors, threeBodyTerm, decomposeOriginal, usingBoxOptimization,
                verbose, pruneClashes, prunePairClashes, master);

        // Update the EliminatedRotamers instance with the EnergyExpansion instance.
        eR.setEnergyExpansion(eE);

        int loaded = 0;
        if (loadEnergyRestart) {
            if (usingBoxOptimization) {
                loaded = eE.loadEnergyRestart(energyRestartFile, residues, boxLoadIndex, boxLoadCellIndices);
            } else {
                loaded = eE.loadEnergyRestart(energyRestartFile, residues);
            }
        }

        long energyStartTime = System.nanoTime();
        WorkerTeam energyWorkerTeam = new WorkerTeam(world);

        try {
            if (loaded < 1) {
                eE.allocateSelfJobMap(residues, nResidues, false);
            }

            SelfEnergyRegion selfEnergyRegion = new SelfEnergyRegion(this, eE, eR, residues, library,
                    energyWriter, world, numProc, pruneClashes, master,
                    rank, verbose, writeEnergyRestart, printFiles);
            energyWorkerTeam.execute(selfEnergyRegion);
            long singlesTime = System.nanoTime() - energyStartTime;
            logIfMaster(format(" Time for single energies: %12.4g", (singlesTime * 1.0E-9)));

            if (loaded < 2) {
                eE.allocate2BodyJobMap(residues, nResidues, false);
            }

            TwoBodyEnergyRegion twoBodyEnergyRegion = new TwoBodyEnergyRegion(
                    this, dM, eE, eR, residues, allResiduesList,
                    library, energyWriter, world, numProc, prunePairClashes, superpositionThreshold,
                    master, rank, verbose, writeEnergyRestart, printFiles);
            energyWorkerTeam.execute(twoBodyEnergyRegion);
            long pairsTime = System.nanoTime() - (singlesTime + energyStartTime);

            long triplesTime = 0;
            long quadsTime = 0;
            logIfMaster(format(" Time for 2-body energies:   %12.4g", (pairsTime * 1.0E-9)));

            if (threeBodyTerm) {
                if (loaded < 3) {
                    eE.allocate3BodyJobMap(residues, nResidues, false);
                }

                ThreeBodyEnergyRegion threeBodyEnergyRegion = new ThreeBodyEnergyRegion(
                        this, dM, eE, eR, residues, allResiduesList, library, energyWriter, world, numProc,
                        superpositionThreshold, master, rank, verbose, writeEnergyRestart, printFiles);
                energyWorkerTeam.execute(threeBodyEnergyRegion);
                triplesTime = System.nanoTime() - (pairsTime + singlesTime + energyStartTime);
                logIfMaster(format(" Time for 3-Body energies: %12.4g", (triplesTime * 1.0E-9)));
            }

            if (compute4BodyEnergy) {
                eE.allocate4BodyJobMap(residues, nResidues);

                FourBodyEnergyRegion fourBodyEnergyRegion = new FourBodyEnergyRegion(
                        this, dM, eE, eR, residues, allResiduesList, superpositionThreshold);
                energyWorkerTeam.execute(fourBodyEnergyRegion);
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
        for (Atom atom : atoms) {
            atom.setUse(true);
        }
        // Print the energy with all rotamers in their default conformation.
        if (verboseEnergies && master) {
            try {
                double defaultEnergy = currentEnergy(residues);
                logger.info(format(" Energy of the system with rotamers in their default conformation: %s",
                        formatEnergy(defaultEnergy)));
            } catch (ArithmeticException ex) {
                logger.severe(format(" Exception %s in calculating default energy; FFX shutting down", ex.toString()));
            }
        }
        return eE.getBackboneEnergy();
    }

    /**
     * Computes the environment/backbone energy, defined as energy with all
     * sidechains under consideration turned off in their 0th rotamer.
     *
     * @param residues Residues under optimization.
     * @return Backbone energy Eenv/bb.
     * @throws ArithmeticException If an exception in calculating energy is
     *                             found.
     */
    public double computeBackboneEnergy(Residue[] residues) throws ArithmeticException {
        // Set all atoms to be used.
        Atom[] atoms = molecularAssembly.getAtomArray();
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

    public double getBackboneEnergy() {
        return eE.getBackboneEnergy();
    }

    public void turnOnResidue(Residue residue, int ri) {
        turnOnAtoms(residue);
        Rotamer[] rotamers = residue.getRotamers(library);
        RotamerLibrary.applyRotamer(residue, rotamers[ri]);
    }

    public void turnOffResidue(Residue residue) {
        turnOffAtoms(residue);
        applyDefaultRotamer(residue);
    }

    public void turnOffAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (Residue residue : residues) {
            turnOffResidue(residue);
        }
    }

    public void turnOnAllResidues(Residue[] residues) {
        if (residues == null) {
            return;
        }
        int nRes = residues.length;
        for (Residue residue : residues) {
            turnOnAtoms(residue);
        }
    }

    /**
     * Find clashes between side-chain rotamers and environment atoms.
     *
     * @param atoms    All atoms in the system.
     * @param crystal  The crystal contain PBC and SymOp info.
     * @param residues The residue being optimized.
     */
    private void computeBackboneRotamerClashes(Atom[] atoms, Crystal crystal, Residue[] residues) {
        // Create a NeighborList with a short cut-off
        double cutoff = superpositionThreshold;
        double buffer = 0.0;
        parallelTeam = (parallelTeam == null) ? new ParallelTeam() : parallelTeam;
        NeighborList backboneNeighborList = new NeighborList(null, crystal,
                atoms, cutoff, buffer, parallelTeam);

        int nAtoms = atoms.length;
        boolean use[] = new boolean[nAtoms];

        int nSymm = crystal.spaceGroup.getNumberOfSymOps();
        double[][] xyz = new double[nSymm][nAtoms * 3];
        int[][][] lists = new int[nSymm][nAtoms][];

        // Turn off all residues.
        turnOffAllResidues(residues);

        // Fill in the coordinate array for SymOp 0 (i.e. the identity operator).
        int index = 0;
        int atomIndex = 1;
        for (Atom atom : atoms) {
            int xyzIndex = atom.getXyzIndex();
            if (atomIndex != xyzIndex) {
                // ToDo: make this compatible with MultiResidues, or make MultiResidues behave better.
                logger.severe(
                        format(" Unexpected atom ordering in RotamerOptimization (Expected: %d, Actual: %d).",
                                xyzIndex, atomIndex));
            }
            atomIndex++;

            xyz[0][index++] = atom.getX();
            xyz[0][index++] = atom.getY();
            xyz[0][index++] = atom.getZ();
        }

        // Loop over all residues
        for (int resIndex = 0; resIndex < residues.length; resIndex++) {
            Residue residue = residues[resIndex];
            // Collect rotamers for this residue.
            Rotamer[] rotamers = residue.getRotamers(library);
            int nrot = rotamers.length;

            // Turn on this residue to its 0th rotamer.
            turnOnResidue(residue, 0);

            List<Atom> resAtoms = residue.getSideChainAtoms();
            Set<Integer> sideChainAtomIndices = new HashSet<>();

            // Configure the per atom "use" flag. Atoms not being used will not be included in the NeighborList.
            for (Atom resAtom : resAtoms) {
                // Stored XYZ index is 1+ its index in the standard Atom arrays.
                int xyzIndex = resAtom.getXyzIndex() - 1;
                use[xyzIndex] = true;
                sideChainAtomIndices.add(xyzIndex);
            }

            // Loop over all rotamers
            for (int ri = 0; ri < nrot; ri++) {
                // Apply the rotamer.
                RotamerLibrary.applyRotamer(residue, rotamers[ri]);

                // Update the coordinate array for SymOp 0 (i.e. the identity operator).
                for (Atom resAtom : resAtoms) {
                    index = (resAtom.getXyzIndex() - 1) * 3;
                    xyz[0][index++] = resAtom.getX();
                    xyz[0][index++] = resAtom.getY();
                    xyz[0][index] = resAtom.getZ();
                }

                // Expand the coordinate array.
                if (nSymm > 1) {
                    // ToDo: refactor coordinate expansion.
                    logger.severe(
                            format(" Backbone clashes with rotamers is not yet supported for symmetry operators."));
                }

                // Build the neighbor list.
                backboneNeighborList.buildList(xyz, lists, use, true, false);

                // Search for a rotamer -> environment clash.
                boolean clash = false;
                clashBreak:
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        // Case 1: check the neighbor list for a side chain atom.
                        if (sideChainAtomIndices.contains(i)) {
                            if (lists[iSymm][i] != null && lists[iSymm][i].length > 0) {
                                clash = true;
                                break clashBreak;
                            }
                        } else if (lists[iSymm][i] != null && lists[iSymm][i].length > 0) {
                            int n = lists[iSymm][i].length;
                            for (int k = 0; k < n; k++) {
                                if (sideChainAtomIndices.contains(lists[iSymm][i][k])) {
                                    clash = true;
                                    break clashBreak;
                                }
                            }
                        }
                    }
                }

                if (clash) {
                    eR.eliminateRotamer(residues, resIndex, ri, verbose);
                    logger.info(format("Eliminated rotamer %s-%d due to backbone clash.", residue, ri));
                }
            }

            // Turn off per atom "use" flags for the current residue.
            for (Integer integer : sideChainAtomIndices) {
                use[integer] = false;
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
            Rotamer[] rotamersi = residuei.getRotamers(library);
            int lenri = rotamersi.length;
            if (minEnergySingles == null || minEnergySingles.length < lenri) {
                minEnergySingles = new double[lenri];
                maxEnergySingles = new double[lenri];
            }
            // Loop over the set of rotamers for residue i.
            for (int ri = 0; ri < lenri; ri++) {
                // Check for an eliminated single.
                if (eR.check(i, ri)) {
                    continue;
                }
                // Start the min/max summation with the self-energy.
                minEnergySingles[ri] = eE.getSelf(i, ri);
                maxEnergySingles[ri] = minEnergySingles[ri];
                for (int j = 0; j < nres; j++) {
                    if (j == i) {
                        continue;
                    }
                    if (eE.minMaxPairEnergy(residues, minMax, i, ri, j)) {
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
                if (eR.check(i, ri)) {
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
                    if (eR.check(i, ri)) {
                        continue;
                    }
                    // If i,ri has a clash with all of phase space, it can be eliminated by something that doesn't clash with all phase space.
                    if (!Double.isFinite(minEnergySingles[ri])) {
                        if (eR.eliminateRotamer(residues, i, ri, print)) {
                            logIfMaster(format("  Rotamer elimination of (%8s,%2d) that always clashes.", residuei.toFormattedString(false, true), ri));
                            eliminated = true;
                        }
                    }
                    // Otherwise, can eliminate if its best possible energy is still worse than something else's worst possible energy.
                    if (minEnergySingles[ri] > eliminationEnergy + ensembleBuffer) {
                        if (eR.eliminateRotamer(residues, i, ri, print)) {
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
                    if (eR.check(i, ri)) {
                        continue;
                    }
                    minPairEnergies[ri] = new double[lenrj];
                    maxPairEnergies[ri] = new double[lenrj];

                    for (int rj = 0; rj < lenrj; rj++) {
                        if (eR.check(j, rj) || eR.check(i, ri, j, rj)) {
                            continue;
                        }
                        minPairEnergies[ri][rj] = eE.getSelf(i, ri) + eE.getSelf(j, rj) + eE.get2Body(i, ri, j, rj);
                        maxPairEnergies[ri][rj] = minPairEnergies[ri][rj];

                        // Min and max external summations for ri-rj.
                        double[] minMax = new double[2];

                        // Add contributions from third residues k, and possibly fourth residues l.
                        if (eE.minMaxE2(residues, minMax, i, ri, j, rj)) {
                            if (Double.isFinite(minPairEnergies[ri][rj]) && Double.isFinite(minMax[0])) {
                                minPairEnergies[ri][rj] += minMax[0];
                            } else {
                                logger.severe(format(" An ri-rj pair %s-%d %s-%d with NaN minimum was caught incorrectly!", residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj));
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
                            logger.info(format(" Eliminating pair %s-%d %s-%d that always clashes.", residuei.toFormattedString(false, true), ri, residuej.toFormattedString(false, true), rj));
                            eR.eliminateRotamerPair(residues, i, ri, j, rj, print);
                            eliminated = true;
                        }
                    }
                }

                double pairEliminationEnergy = Double.MAX_VALUE;
                for (int ri = 0; ri < lenri; ri++) {
                    if (eR.check(i, ri)) {
                        continue;
                    }
                    for (int rj = 0; rj < lenrj; rj++) {
                        if (eR.check(j, rj) || eR.check(i, ri, j, rj)) {
                            continue;
                        }
                        if (Double.isFinite(maxPairEnergies[ri][rj]) && maxPairEnergies[ri][rj] < pairEliminationEnergy) {
                            pairEliminationEnergy = maxPairEnergies[ri][rj];
                        }
                    }
                }

                if (pairEliminationEnergy == Double.MAX_VALUE) {
                    logIfMaster(format(" All rotamer pairs for residues %s and %s have possible conflicts; cannot perform any eliminations!", residuei.toFormattedString(false, true), residuej), Level.FINE);
                } else {
                    double comparisonEnergy = pairEliminationEnergy + ensembleBuffer;
                    for (int ri = 0; ri < lenri; ri++) {
                        if (eR.check(i, ri)) {
                            continue;
                        }
                        for (int rj = 0; rj < lenrj; rj++) {
                            if (eR.check(j, rj) || eR.check(i, ri, j, rj)) {
                                continue;
                            }
                            if (minPairEnergies[ri][rj] > comparisonEnergy) {
                                if (eR.eliminateRotamerPair(residues, i, ri, j, rj, print)) {
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

                if (eR.pairsToSingleElimination(residues, i, j)) {
                    eliminated = true;
                }
            }
        }

        return eliminated;
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
                if (eR.check(i, riA)) {
                    continue;
                }
                for (int riB = 0; riB < nri; riB++) {
                    // The eliminating rotamer cannot be riA and must be a valid.
                    if (riA == riB || eR.check(i, riB)) {
                        continue;
                    }
                    if (goldsteinElimination(residues, i, riA, riB)) {
                        eliminated = true;
                        break;
                    }
                }
            }
        }
        if (!eliminated) {
            logIfMaster(" No more single rotamers to eliminate.");
        }
        return eliminated;
    }

    /**
     * Checks if residue i is considered to be interacting with residue j, and
     * thus has non-null elements in the pair energies matrix.
     *
     * @param i A residue index.
     * @param j A residue index j != i.
     * @return If i and j interact.
     */
    public boolean checkNeighboringPair(int i, int j) {
        assert i != j;
        final int first;
        final int second;
        if (i > j) {
            first = j;
            second = i;
        } else {
            first = i;
            second = j;
        }
        return Arrays.stream(resNeighbors[first]).anyMatch(l -> l == second);
    }

    /**
     * Checks if residue i is considered to be interacting with residue j, that
     * residue k is interacting with either i or j, and thus i-j-k has non-null
     * elements in the triple energies matrix.
     *
     * @param i A residue index.
     * @param j A residue index j != i.
     * @param k A residue index k != i, k != j.
     * @return If i, j, and k form an interacting triple.
     */
    public boolean checkNeighboringTriple(int i, int j, int k) {
        assert i != j && i != k && j != k;
        int[] vals = new int[]{i, j, k};
        Arrays.sort(vals);
        final int first = vals[0];
        final int second = vals[1];
        final int third = vals[2];
        if (!checkNeighboringPair(first, second)) {
            return false;
        }
        return (checkNeighboringPair(first, third) || checkNeighboringPair(second, third));
    }

    /**
     * Attemps to eliminate rotamer riA based on riB.
     *
     * @param residues
     * @param i
     * @param riA      Rotamer to attempt elimination of.
     * @param riB      Rotamer to attempt elimination by.
     * @return If riA was eliminated.
     */
    private boolean goldsteinElimination(Residue[] residues, int i, int riA, int riB) {
        Residue resi = residues[i];

        // Initialize Goldstein inequality.
        double selfDiff = eE.getSelf(i, riA) - eE.getSelf(i, riB);
        double goldsteinEnergy = selfDiff;

        double sumPairDiff = 0.0;
        double sumTripleDiff = 0.0;

        // Loop over a 2nd residue j.
        for (int j : bidiResNeighbors[i]) {
            Residue resj = residues[j];
            Rotamer rotj[] = resj.getRotamers(library);
            int nrj = rotj.length;
            double minForResJ = Double.MAX_VALUE;
            double minPairDiff = 0.0;
            double minTripleDiff = 0.0;
            int rjEvals = 0;

            // Loop over the rotamers for residue j.
            for (int rj = 0; rj < nrj; rj++) {
                if (eR.check(j, rj)) {
                    continue;
                }

                if (eR.check(i, riA, j, rj)) {
                    continue; // This is not a part of configuration space accessible to riA.
                }
                if (eR.check(i, riB, j, rj)) {
                    /**
                     * This is a part of configuration space where riA is valid
                     * but not riB. Thus, if j,rj is part of the GMEC, riB is
                     * inconsistent with it. Thus, riB cannot be used to
                     * eliminate riA.
                     */
                    return false;
                }

                double pairI = eE.get2Body(i, riA, j, rj);
                double pairJ = eE.get2Body(i, riB, j, rj);
                double pairDiff = pairI - pairJ;

                rjEvals++;

                // Include three-body interactions.
                double tripleDiff = 0.0;
                if (threeBodyTerm) {
                    IntStream kStream = IntStream.concat(Arrays.stream(bidiResNeighbors[i]), Arrays.stream(bidiResNeighbors[j]));
                    int[] possKs = kStream.distinct().sorted().toArray();
                    for (int k : possKs) {
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
                            if (eR.check(k, rk) || eR.check(j, rj, k, rk) || eR.check(i, riA, k, rk)) {
                                //Not yet implemented: check(i, riA, j, rj, k, rk) because no triples get eliminated.
                                continue;
                            }
                            /**
                             * If i,riB-k,rk or i,riB-j,rj-k,rk are invalid for
                             * riB, there is some part of configuration space
                             * for which riA is valid but not riB.
                             */
                            if (eR.check(i, riB, k, rk)) {
                                // Not yet implemented: check(i, riB, j, rj, k, rk).
                                return false;
                            }

                            rkEvals++;
                            double e = eE.get3Body(residues, i, riA, j, rj, k, rk) - eE.get3Body(residues, i, riB, j, rj, k, rk);
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
            if (eR.eliminateRotamer(residues, i, riA, print)) {
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
                if (eR.check(i, riA)) {
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
                        if (eR.check(j, rjC) || eR.check(i, riA, j, rjC)) {
                            continue;
                        }
                        boolean breakOut = false;

                        // Now, generate pairs riB-rjD. If any pair riB-rjD eliminates riA-rjC, break out of the loop.
                        for (int riB = 0; riB < lenri; riB++) {
                            if (breakOut) {
                                break;
                            }
                            if (eR.check(i, riB)) {
                                continue;
                            }
                            for (int rjD = 0; rjD < lenrj; rjD++) {
                                if (breakOut) {
                                    break;
                                }
                                // Do not attempt eliminating with an eliminated pair.
                                if (eR.check(j, rjD) || eR.check(i, riB, j, rjD)) {
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
                    if (eR.pairsToSingleElimination(residues, i, j)) {
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
     * @param i        Index of the first residue.
     * @param riA      Index of the first residue's rotamer to eliminate.
     * @param riB      Index of the first residue's rotamer to use for elimination.
     * @param j        Index of the 2nd residue.
     * @param rjC      Index of the 2nd residue's rotamer to eliminate.
     * @param rjD      Index of the 2nd residue's rotamer to use for elimination.
     * @return Return true if eliminated.
     */
    private boolean goldsteinPairElimination(Residue[] residues,
                                             int i, int riA, int riB, int j, int rjC, int rjD) {

        ArrayList<Residue> missedResidues = null;

        // Initialize the Goldstein energy.
        double goldsteinEnergy = eE.getSelf(i, riA) + eE.getSelf(j, rjC) + eE.get2Body(i, riA, j, rjC)
                - eE.getSelf(i, riB) - eE.getSelf(j, rjD) - eE.get2Body(i, riB, j, rjD);

        try {
            if (parallelTeam == null) {
                parallelTeam = new ParallelTeam();
            }
            if (goldsteinPairRegion == null) {
                goldsteinPairRegion = new GoldsteinPairRegion(parallelTeam.getThreadCount());
            }
            goldsteinPairRegion.init(residues, i, riA, riB, j, rjC, rjD, bidiResNeighbors, this);
            parallelTeam.execute(goldsteinPairRegion);
            goldsteinEnergy += goldsteinPairRegion.getSumOverK();
            missedResidues = goldsteinPairRegion.getMissedResidues();
        } catch (Exception e) {
            logger.log(Level.WARNING, " Exception in GoldsteinPairRegion.", e);
        }
        if (missedResidues != null && !missedResidues.isEmpty()) {
            logIfMaster(format(" Skipping energy comparison due to a missed residue: i %d riA %d riB %d j %d rjC %d rjD %d", i, riA, riB, j, rjC, rjD), Level.FINE);
            return false;
        }

        if (goldsteinEnergy > ensembleBuffer) {
            if (missedResidues.isEmpty()) {
                if (eR.eliminateRotamerPair(residues, i, riA, j, rjC, print)) {
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

    public double goldsteinPairSumOverK(Residue[] residues, int lb, int ub, int i, int riA, int riB,
                                        int j, int rjC, int rjD,
                                        ArrayList<Residue> blockedResidues,
                                        int[] possK) {
        double sumOverK = 0.0;

        for (int indK = lb; indK <= ub; indK++) {
            int k = possK[indK];
            double minForResK = Double.MAX_VALUE;
            Residue resk = residues[k];
            Rotamer rotk[] = resk.getRotamers(library);
            int nrk = rotk.length;
            int rkEvals = 0;
            // Loop over residue k's rotamers.
            for (int rk = 0; rk < nrk; rk++) {
                if (eR.check(k, rk)) {
                    continue;
                }
                // Continue if k,rk invalid for riA/rjC.
                if (eR.check(i, riA, k, rk) || eR.check(j, rjC, k, rk)) {
                    // Not implemented: check(i, riA, j, rjC, k, rk).
                    continue;
                }
                // Return false if k,rk invalid for riB/rjD.
                if (eR.check(i, riB, k, rk) || eR.check(j, rjD, k, rk)) {
                    blockedResidues.add(resk);
                    return Double.NaN;
                }

                rkEvals++;
                double currentResK = eE.get2Body(i, riA, k, rk) - eE.get2Body(i, riB, k, rk)
                        + eE.get2Body(j, rjC, k, rk) - eE.get2Body(j, rjD, k, rk);
                // Include 3-body effects.
                if (threeBodyTerm) {
                    double sumOverL = (eE.get3Body(residues, i, riA, j, rjC, k, rk)
                            - eE.get3Body(residues, i, riB, j, rjD, k, rk));
                    // Loop over a 4th residue l.
                    int[] nK = bidiResNeighbors[k];
                    IntStream lStream = IntStream.concat(Arrays.stream(possK), Arrays.stream(nK));
                    int[] possL = lStream.filter(l -> (l != i && l != j && l != k)).
                            sorted().
                            distinct().
                            toArray();

                    for (int l : possL) {
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
                            if (eR.check(l, rl) || eR.check(k, rk, l, rl) || eR.check(i, riA, l, rl) || eR.check(j, rjC, l, rl)) {
                                // Not implemented: check(i, riA, j, rjC, l, rl) || check(i, riA, k, rk, l, rl) || check(j, rjC, k, rk, l, rl) || check(i, riA, j, rjC, k, rk, l, rl)
                                continue;
                            }
                            if (eR.check(i, riB, l, rl) || eR.check(j, rjD, l, rl)) {
                                // Not implemented: check(i, riB, j, rjD, l, rl) || check(i, riB, k, rk, l, rl) || check(j, rjD, k, rk, l, rl) || check(i, riB, j, rjD, k, rk, l, rl)
                                blockedResidues.add(residuel);
                                return Double.NaN;
                            }
                            rlEvaluations++;
                            double e = eE.get3Body(residues, i, riA, k, rk, l, rl) - eE.get3Body(residues, i, riB, k, rk, l, rl)
                                    + eE.get3Body(residues, j, rjC, k, rk, l, rl) - eE.get3Body(residues, j, rjD, k, rk, l, rl);
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
     * <p>
     * Setter for the field <code>threeBodyCutoffDist</code>.</p>
     *
     * @param dist a double.
     */
    public void setThreeBodyCutoffDist(double dist) {
        this.threeBodyCutoffDist = dist;
        if (threeBodyCutoffDist < 0) {
            logger.info(format("Warning: threeBodyCutoffDist should not be less than 0."));
        }
    }

    /**
     * <p>
     * Setter for the field <code>superpositionThreshold</code>.</p>
     *
     * @param superpositionThreshold a double.
     */
    public void setSuperpositionThreshold(double superpositionThreshold) {
        this.superpositionThreshold = superpositionThreshold;
    }

    /**
     * <p>
     * setGoldstein.</p>
     *
     * @param set a boolean.
     */
    public void setGoldstein(boolean set) {
        this.useGoldstein = set;
    }

    public void setRotamerLibrary(RotamerLibrary lib) {
        library = lib;
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

    /**
     * <p>
     * setSingletonClashThreshold.</p>
     *
     * @param singletonClashThreshold a double.
     */
    public void setSingletonClashThreshold(double singletonClashThreshold) {
        this.clashThreshold = singletonClashThreshold;
    }

    /**
     * <p>
     * Setter for the field <code>pairClashThreshold</code>.</p>
     *
     * @param pairClashThreshold a double.
     */
    public void setPairClashThreshold(double pairClashThreshold) {
        this.pairClashThreshold = pairClashThreshold;
    }

    /**
     * <p>
     * Setter for the field <code>windowSize</code>.</p>
     *
     * @param windowSize a int.
     */
    public void setWindowSize(int windowSize) {
        this.windowSize = windowSize;
        if (this.increment > windowSize) {
            logger.info(format(" Decreasing increment to match window size %d", windowSize));
            this.increment = windowSize;
        }
    }

    /**
     * <p>
     * setEnsemble.</p>
     *
     * @param ensemble       a int.
     * @param ensembleBuffer a double.
     */
    public void setEnsemble(int ensemble, double ensembleBuffer) {
        this.ensembleNumber = ensemble;
        this.ensembleBuffer = ensembleBuffer;
        this.ensembleBufferStep = 0.1 * ensembleBuffer;
        if (ensemble > 1) {
            setPruning(0);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (eR != null) {
            return eR.toString();
        } else return null;
    }

    /**
     * {@inheritDoc}
     */
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
     * <p>
     * Setter for the field <code>energyRestartFile</code>.</p>
     *
     * @param file a {@link java.io.File} object.
     */
    public void setEnergyRestartFile(File file) {
        loadEnergyRestart = true;
        energyRestartFile = file;
    }

    /**
     * Sets the twoBodyCutoffDist. All two-body energies where the rotamers have
     * a separation distance larger than the cutoff are set to 0.
     *
     * @param twoBodyCutoffDist Separation distance at which the interaction of
     *                          two side-chains is assumed to have an energy of 0.
     */
    public void setTwoBodyCutoff(double twoBodyCutoffDist) {
        this.twoBodyCutoffDist = twoBodyCutoffDist;
        if (this.twoBodyCutoffDist < 0) {
            logger.info(format("Warning: threeBodyCutoffDist should not be less than 0."));
        }
    }

    /**
     * Sets the threeBodyCutoffDist. All three-body energies where the rotamers
     * have a separation distance larger than the cutoff are set to 0.
     *
     * @param threeBodyCutoffDist Separation distance at which the interaction
     *                            of three side-chains is assumed to have an energy of 0.
     */
    public void setThreeBodyCutoff(double threeBodyCutoffDist) {
        this.threeBodyCutoffDist = threeBodyCutoffDist;
        if (this.threeBodyCutoffDist < 0) {
            logger.info(format("Warning: threeBodyCutoffDist should not be less than 0."));
        }
    }

    /**
     * Sets the monteCarloTesting boolean in RotamerOptimization.java to true or
     * false. This should only be set to true when monte carlo is being tested
     * through the ManyBodyTest.java script. When true, the method sets a seed
     * for the pseudo-random number generator and allows the monte carlo rotamer
     * optimization to be deterministic.
     *
     * @param bool True or false.
     */
    public void setMonteCarloTesting(boolean bool) {
        this.monteCarloTesting = bool;

    }


    /**
     * False unless JUnit testing.
     */
    private boolean testing = false;
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

    /**
     * Method allows for testing of the elimination criteria by setting
     * parameters appropriately.
     *
     * @param testing True only during RotamerOptimizationTest.
     */
    void setTestOverallOpt(boolean testing) {
        this.testing = testing;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    /**
     * Method allows for testing of the elimination criteria by setting
     * parameters appropriately for a specific test case.
     *
     * @param testSelfEnergyEliminations True when only self energies are
     *                                   calculated; pairs, triples, etc., are assumed to be 0.
     */
    void setTestSelfEnergyEliminations(boolean testSelfEnergyEliminations) {
        this.testing = true;
        this.testSelfEnergyEliminations = testSelfEnergyEliminations;
        testPairEnergyEliminations = -1;
        testTripleEnergyEliminations1 = -1;
        testTripleEnergyEliminations2 = -1;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    /**
     * Method allows for testing of the elimination criteria by setting
     * parameters appropriately for a specific test case.
     *
     * @param testPairEnergyEliminations True when only pair energies are
     *                                   calculated; selves, triples, etc., are assumed to be 0.
     */
    void setTestPairEnergyEliminations(int testPairEnergyEliminations) {
        this.testing = true;
        this.testPairEnergyEliminations = testPairEnergyEliminations;
        testSelfEnergyEliminations = false;
        testTripleEnergyEliminations1 = -1;
        testTripleEnergyEliminations2 = -1;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    /**
     * Method allows for testing of the elimination criteria by setting
     * parameters appropriately for a specific test case.
     *
     * @param testTripleEnergyEliminations1 True when only triple energies are
     *                                      calculated; selves, pairs, etc., are assumed to be 0.
     * @param testTripleEnergyEliminations2 True when only triple energies are
     *                                      calculated; selves, pairs, etc., are assumed to be 0.
     */
    void setTestTripleEnergyEliminations(int testTripleEnergyEliminations1, int testTripleEnergyEliminations2) {
        this.testing = true;
        this.testTripleEnergyEliminations1 = testTripleEnergyEliminations1;
        this.testTripleEnergyEliminations2 = testTripleEnergyEliminations2;
        testSelfEnergyEliminations = false;
        testPairEnergyEliminations = -1;
        distanceMethod = DistanceMethod.ROTAMER;
        setTwoBodyCutoff(Double.MAX_VALUE);
        setThreeBodyCutoff(9.0);
    }

    /**
     * Test the self-energy elimination by setting 2-body and 3-body
     * interactions to zero.
     *
     * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
     */
    private void testSelfEnergyElimination(Residue[] residues) {
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
                            eE.set2Body(i, ri, j, rj, 0, true);
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
                                        eE.set3Body(residues, i, ri, j, rj, k, rk, 0, true);
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
     * Test the elimination criteria by setting self and 3-body interactions to
     * zero.
     *
     * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
     * @param resID    a int.
     */
    private void testPairEnergyElimination(Residue[] residues, int resID) {
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
                    eE.setSelf(i, ri, 0, true);
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
                                eE.set2Body(i, ri, j, rj, 0, true);
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
                                        eE.set3Body(residues, i, ri, j, rj, k, rk, 0, true);
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
     * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
     * @param resID1   The residue number for one of two fixed residues.
     * @param resID2   The second residue number for one of two fixed residues.
     */
    private void testTripleEnergyElimination(Residue[] residues, int resID1, int resID2) {
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
                    eE.setSelf(i, ri, 0, true);
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
                            eE.set2Body(i, ri, j, rj, 0, true);
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
                                            eE.set3Body(residues, i, ri, j, rj, k, rk, 0, true);
                                        } catch (Exception e) {
                                            // catch NPE.
                                        }
                                    }
                                    if (i != resID2 && j != resID2 && k != resID2) {
                                        try {
                                            eE.set3Body(residues, i, ri, j, rj, k, rk, 0, true);
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

    //    /**
//     * Writes eliminated singles and pairs to a CSV file. Reads in .log file.
//     *
//     * @throws java.io.FileNotFoundException if any.
//     * @throws java.io.IOException if any.
//     */
//    public void outputEliminated() throws FileNotFoundException, IOException {
//        /**
//         * eliminated.csv stores the eliminated singles and pairs.
//         */
//        File fileName = new File("eliminated.csv");
//        String path = fileName.getCanonicalPath();
//        File outputTxt = new File(path);
//        Path currentRelativePath = Paths.get("eliminated.csv");
//        String logPath = currentRelativePath.toAbsolutePath().toString();
//        logPath = logPath.replaceAll("/eliminated.csv", "");
//        File logDirectory = new File(logPath);
//        /**
//         * Searches for .log file within the directory.
//         */
//        File[] logArray = logDirectory.listFiles(new FilenameFilter() {
//            public boolean accept(File dir, String filename) {
//                return filename.endsWith(".log");
//            }
//        });
//        File logFile = logArray[0];
//        PrintWriter out = new PrintWriter(new FileWriter(outputTxt, true));
//        BufferedReader br = new BufferedReader(new FileReader(logFile));
//        String line;
//        String line3;
//        List<String> boxes = new ArrayList<String>();
//        while ((line3 = br.readLine()) != null) {
//            if (line3.contains("-a, 5")) {
//                /**
//                 * Stores the number of box optimization iterations and the
//                 * coordinates of each box.
//                 */
//                while ((line = br.readLine()) != null) {
//                    if (line.contains("xyz indices")) {
//                        String coordinates = line.replaceAll(" Cell xyz indices:", "");
//                        String nextLine = br.readLine();
//                        if (!nextLine.contains("Empty box")) {
//                            boxes.add(coordinates);
//                            String nextLine2;
//                            String nextLine4;
//                            String nextLine5;
//                            while (!(nextLine2 = br.readLine()).contains("Collecting Permutations")) {
//                                if (nextLine2.contains("Applying Rotamer Pair")) {
//                                    int total = 0;
//                                    int total2 = 0;
//                                    nextLine2 = br.readLine();
//                                    nextLine4 = br.readLine();
//                                    nextLine5 = br.readLine();
//                                    if (nextLine5.contains("Self-consistent")) {
//                                        String[] split = nextLine2.split(" ");
//                                        int singlesAmount = Integer.parseInt(split[1]);
//                                        if (singlesIterations.size() > 0) {
//                                            total = singlesIterations.get(singlesIterations.size() - 1);
//                                        }
//                                        singlesIterations.add(singlesAmount + total);
//                                        System.out.println(nextLine4);
//                                        String[] split2 = nextLine4.split(" ");
//                                        int pairsAmount = Integer.parseInt(split2[1]);
//                                        if (pairsIterations.size() > 0) {
//                                            total2 = pairsIterations.get(pairsIterations.size() - 1);
//                                        }
//                                        pairsIterations.add(pairsAmount + total2);
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                /**
//                 * Writes eliminated singles and pairs to eliminated.csv, sorted
//                 * by box iteration.
//                 */
//                out.append("Eliminated Singles\n");
//                out.append("Residue,Rotamer\n");
//                out.append(boxes.get(0) + "\n");
//                int stopper = 0;
//                for (int i = 0; i < eliminatedResidue.size(); i++) {
//                    while (stopper == 0) {
//                        for (int x = 0; x < (singlesIterations.size() - 1); x++) {
//                            if (singlesIterations.get(x) == 0) {
//                                String boxHeader = (boxes.get(x + 1));
//                                out.append(boxHeader + System.lineSeparator());
//                            }
//                            stopper++;
//                        }
//                    }
//                    String first = eliminatedResidue.get(i).toString();
//                    String second = eliminatedRotamer.get(i).toString();
//                    String resrot = first + "," + second;
//                    out.append(resrot + System.lineSeparator());
//                    for (int x = 0; x < (singlesIterations.size() - 1); x++) {
//                        if ((i + 1) == singlesIterations.get(x)) {
//                            String boxHeader = (boxes.get(x + 1));
//                            out.append(boxHeader + System.lineSeparator());
//                        }
//                    }
//                }
//                out.append("Eliminated Pairs\n");
//                out.append("Residue1,Rotamer1,Residue2,Rotamer2" + System.lineSeparator());
//                out.append(boxes.get(0) + "\n");
//                int stopper2 = 0;
//                for (int i = 0; i < eliminatedResidue1.size(); i++) {
//                    while (stopper2 == 0) {
//                        for (int x = 0; x < (pairsIterations.size() - 1); x++) {
//                            if (pairsIterations.get(x) == 0) {
//                                String boxHeader = (boxes.get(x + 1));
//                                out.append(boxHeader + System.lineSeparator());
//                            }
//                            stopper2++;
//                        }
//                    }
//                    String first = eliminatedResidue1.get(i).toString();
//                    String second = eliminatedRotamer1.get(i).toString();
//                    String third = eliminatedResidue2.get(i).toString();
//                    String fourth = eliminatedRotamer2.get(i).toString();
//                    String resrot = first + "," + second + "," + third + "," + fourth;
//                    out.append(resrot + System.lineSeparator());
//                    for (int x = 0; x < (pairsIterations.size() - 1); x++) {
//                        if ((i + 1) == pairsIterations.get(x)) {
//                            String boxHeader = (boxes.get(x + 1));
//                            out.append(boxHeader + System.lineSeparator());
//                        }
//                    }
//                }
//                /**
//                 * Writes eliminated singles and pairs to eliminated.csv for
//                 * global optimization.
//                 */
//            } else if ((line3.contains("-a, 1")) || (line3.contains("-a, 2")) || (line3.contains("-a, 3")) || (line3.contains("-a, 4"))) {
//                out.append("Eliminated Singles\n");
//                out.append("Residue,Rotamer\n");
//                for (int i = 0; i < eliminatedResidue.size(); i++) {
//                    String first = eliminatedResidue.get(i).toString();
//                    String second = eliminatedRotamer.get(i).toString();
//                    String resrot = first + "," + second;
//                    out.append(resrot + System.lineSeparator());
//                }
//                out.append("Eliminated Pairs\n");
//                out.append("Residue1,Rotamer1,Residue2,Rotamer2" + System.lineSeparator());
//                for (int i = 0; i < eliminatedResidue1.size(); i++) {
//                    String first = eliminatedResidue1.get(i).toString();
//                    String second = eliminatedRotamer1.get(i).toString();
//                    String third = eliminatedResidue2.get(i).toString();
//                    String fourth = eliminatedRotamer2.get(i).toString();
//                    String resrot = first + "," + second + "," + third + "," + fourth;
//                    out.append(resrot + System.lineSeparator());
//                }
//            }
//        }
//        out.close();
//        br.close();
//    } 
}
