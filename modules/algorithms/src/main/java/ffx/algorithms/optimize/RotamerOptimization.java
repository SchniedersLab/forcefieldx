// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
// ******************************************************************************
package ffx.algorithms.optimize;

import edu.rit.pj.Comm;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.WorkerTeam;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.algorithms.mc.MCMove;
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
import ffx.numerics.Potential;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.openmm.OpenMMEnergy;
import ffx.potential.parameters.TitrationUtils;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Constants;
import ffx.utilities.ObjectPair;
import ffx.utilities.Resources;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

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
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static ffx.potential.bonded.Residue.ResidueType.NA;
import static ffx.potential.bonded.RotamerLibrary.applyRotamer;
import static java.lang.Boolean.parseBoolean;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;

/**
 * Optimize protein side-chain conformations and nucleic acid backbone conformations using rotamers.
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
  /**
   * Fallback if there is no vdW node.
   */
  private static final double FALLBACK_TWO_BODY_CUTOFF = 0;
  /**
   * MolecularAssembly to perform rotamer optimization on.
   */
  protected MolecularAssembly molecularAssembly;
  /**
   * The Potential to evaluate during rotamer optimization.
   */
  protected final Potential potential;
  /**
   * AlgorithmListener who should receive updates as the optimization runs.
   */
  protected final AlgorithmListener algorithmListener;
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
  protected final boolean rank0;
  /**
   * Flag to control the verbosity of printing.
   */
  private final boolean print = false;
  /**
   * Flag to calculate and print additional energies (mostly for debugging).
   */
  private final boolean verboseEnergies = true;

  /**
   * If true, write out an energy restart file.
   */
  protected boolean writeEnergyRestart = true;
  /**
   * Parameters for box optimization are stored here.
   */
  BoxOptimization boxOpt;
  /**
   * Represents the method called to obtain the directory corresponding to the current energy; will
   * be a simple return null for potential energy evaluations. While current energy calls will fill
   * the rotamer list with the current rotamers of the residue, other methods may skip applying the
   * rotamer directly.
   */
  private final BiFunction<List<Residue>, List<Rotamer>, File> dirSupplier;
  /**
   * Represents the method called to obtain energy for the current rotamer or state; defaults to the
   * existing potential energy code. May discard the input file.
   */
  private final ToDoubleFunction<File> eFunction;
  /**
   * Flag to indicate verbose logging.
   */
  private final boolean verbose;
  /**
   * The DistanceMatrix class handles calculating distances between residues.
   */
  private DistanceMatrix dM;
  /**
   * The EnergyExpansion class compute terms in the many-body energy expansion.
   */
  private EnergyExpansion eE;
  /**
   * The EliminatedRotamers class tracks eliminated rotamers and rotamer paris.
   */
  private EliminatedRotamers eR;
  /**
   * RotamerLibrary instance.
   */
  protected RotamerLibrary library = RotamerLibrary.getDefaultLibrary();
  /**
   * Parallel evaluation of quantities used during Goldstein Pair elimination.
   */
  private GoldsteinPairRegion goldsteinPairRegion;
  /**
   * Parallel evaluation of many-body energy sums.
   */
  private EnergyRegion energyRegion;
  /**
   * Flag to indicate a request to terminate the optimization.
   */
  private boolean terminate = false;
  /**
   * Flag to indicate if the algorithm is running (done == false) or completed (done == true).
   */
  private boolean done = true;
  /**
   * Two-Body cutoff distance.
   */
  private double twoBodyCutoffDist;
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
   * The temperature used during either Monte Carlo sampling,
   * calculation of Boltzmann weights and the partition function.
   * The default is 298.15 Kelvin.
   */
  private double temperature = 298.15;
  /**
   * Total boltzmann is the value of the partition function.
   */
  private double totalBoltzmann = 0;
  /**
   * Reference energy for calculating boltzmanns in the partition function
   */
  private double refEnergy = 0;
  /**
   * Rotamer populations from the partition function
   */
  private double[][] fraction;
  /**
   * Botlzmann weights of every rotamer
   */
  private double[][] populationBoltzmann;
  /**
   * pH during titration rotamer optimization
   */
  private double pH;
  /**
   * Energy restraint on titration
   */
  private double pHRestraint = 0;
  /**
   * Recompute the self energies from restart when changing the pH
   */
  private boolean recomputeSelf = false;
  /**
   * True when running the GenZ algorithm
   */
  boolean genZ = false;
  /**
   * List of residues to optimize; they may not be contiguous or all members of the same chain.
   */
  private List<Residue> residueList;
  /**
   * This is the optimal rotamers corresponding to residueList.
   */
  protected int[] optimum;

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
  protected boolean revert;
  /**
   * The sliding window direction.
   */
  protected Direction direction = Direction.FORWARD;

  /**
   * The distance that the distance matrix checks for.
   */
  private double distance = 2.0;
  /**
   * Default distance method is to find the shortest distance between residues.
   */
  private DistanceMethod distanceMethod = DistanceMethod.RESIDUE;
  /**
   * The algorithm to use for rotamer optimization.
   */
  private Algorithm algorithm = null;

  /**
   * Flag to indicate use of the Goldstein criteria instead of the less stringent Dead-End
   * Elimination criteria.
   */
  private boolean useGoldstein = true;
  /**
   * The number of most-favorable structures to include as output.
   */
  private int ensembleNumber = 1;
  /**
   * The energy buffer applied to each elimination criteria to affect an ensemble.
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
   * Clash energy threshold (kcal/mole).
   */
  private double clashThreshold = 25.0;
  /**
   * Clash energy threshold (kcal/mol) for MultiResidues, which can have much more variation in self
   * and 2-Body energies.
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
   * An array of atomic coordinates (length 3 * the number of atoms).
   */
  private double[] x = null;
  /**
   * A flag to indicate use of the full N-Body AMOEBA potential energy during the rotamer
   * optimization.
   */
  private boolean useForceFieldEnergy = false;
  /**
   * Threshold to eliminate nucleic acid Rotamers based on excessive correction distances; 0
   * indicates the threshold is not being implemented.
   */
  private double nucleicCorrectionThreshold = 0;
  /**
   * The approximate energy from summing the backbone energy, self-energies, pair-energies, etc.
   */
  private double approximateEnergy = 0;
  /**
   * Minimum number of nucleic acid Rotamers to use for rotamer optimization should some be
   * eliminated by the nucleic correction threshold.
   */
  private int minNumberAcceptedNARotamers = 10;
  /**
   * Factor by which to multiply the pruning constraints for nucleic acids.
   */
  private double nucleicPruningFactor = 10.0;
  /**
   * The arithmetic mean of 1.0 and the pruning factor, and is applied for AA-NA pairs.
   */
  private double nucleicPairsPruningFactor = ((1.0 + nucleicPruningFactor) / 2);
  /**
   * A list of all residues in the system, which is used to compute a distance matrix.
   */
  private List<Residue> allResiduesList = null;
  /**
   * An array of all residues in the system, which is used to compute a distance matrix.
   */
  private Residue[] allResiduesArray = null;
  /**
   * Number of residues being optimized.
   */
  private int nAllResidues = 0;
  /**
   * The array of optimum rotamers for the subset of residues being optimized during box or window
   * optimization.
   */
  protected int[] optimumSubset;
  /**
   * If true, load an energy restart file.
   */
  protected boolean loadEnergyRestart = false;
  /**
   * Energy restart File instance.
   */
  private File energyRestartFile;
  /**
   * ParallelTeam instance.
   */
  private ParallelTeam parallelTeam;
  /**
   * Flag to indicate computation of 4-body energy values. This is limited to the study 4-body energy
   * magnitudes, but is not included in the rotamer optimization.
   */
  private boolean compute4BodyEnergy = false;
  /**
   * Flag to indicate use of box optimization.
   */
  protected boolean usingBoxOptimization = false;
  /**
   * If a pair of residues have two atoms closer together than the superposition threshold, the
   * energy is set to NaN.
   */
  private double superpositionThreshold = 0.25;
  /**
   * Flag to indicate computation of a many-body expansion for original coordinates.
   */
  private boolean decomposeOriginal = false;
  /**
   * Flag to indicate use of MC optimization.
   */
  private boolean monteCarlo = false;
  /**
   * Number of MC optimization steps.
   */
  private int nMCSteps = 1000000;
  /**
   * Check to see if proposed move has an eliminated 2-body or higher-order term.
   */
  private boolean mcUseAll = false;
  /**
   * Skips brute force enumeration in favor of pure Monte Carlo. Recommended only for testing.
   */
  private boolean mcNoEnum = false;
  /**
   * Sets whether files should be printed; true for standalone applications, false for some
   * applications which use rotamer optimization as part of a larger process.
   */
  protected boolean printFiles = true;
  /**
   * Stores states of each ensemble if printFiles is false.
   */
  private List<ObjectPair<ResidueState[], Double>> ensembleStates;
  /**
   * Maximum depth to check if a rotamer can be eliminated.
   */
  private int maxRotCheckDepth;
  /**
   * Writes energies to restart file.
   */
  protected BufferedWriter energyWriter;

  /**
   * False unless JUnit testing.
   */
  private boolean testing = false;
  /**
   * False unless ManyBodyTest is occurring.
   */
  private boolean monteCarloTesting = false;
  /**
   * Test Self-Energy Elimination.
   */
  private boolean testSelfEnergyEliminations = false;
  /**
   * Test Pair-Energy Elimination.
   *
   * <p>If greater than or equal to 0, test the specified residue.
   */
  private int testPairEnergyEliminations = -1;
  /**
   * Test Triple-Energy Elimination.
   *
   * <p>If greater than or equal to 0, test the specified residues.
   */
  private int testTripleEnergyEliminations1 = -1;
  /**
   * Test Triple-Energy Elimination.
   *
   * <p>If greater than or equal to 0, test the specified residues.
   */
  private int testTripleEnergyEliminations2 = -1;
  /**
   * Only for unit testing; turns off the singles elimination criterion.
   */
  private boolean selfEliminationOn = true;
  /**
   * Only for unit testing; turns off the pairs elimination criterion.
   */
  private boolean pairEliminationOn = true;

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

    rank0 = rank == 0;

    boxOpt = new BoxOptimization(this);

    CompositeConfiguration properties = molecularAssembly.getProperties();
    verbose = properties.getBoolean("verbose", false);

    // Monte Carlo Options.
    String temp = properties.getString("temperature", "298.15");
    if (temp != null) {
      this.temperature = parseDouble(temp);
      logIfRank0(format(" Temperature: %10.6f", this.temperature));
    }

    // Set the default 2-body Cutoff to the van der Waals cutoff.
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    VanDerWaals vdW = forceFieldEnergy.getVdwNode();
    if (vdW != null) {
      NonbondedCutoff nonBondedCutoff = vdW.getNonbondedCutoff();
      twoBodyCutoffDist = nonBondedCutoff.off;
    } else {
      twoBodyCutoffDist = FALLBACK_TWO_BODY_CUTOFF;
    }

    // Process properties; most are seldom used and not available in ManyBodyOptions.
    // Handle testing flags.
    boolean testing = properties.getBoolean("manybody-testing", false);
    if (testing) {
      turnRotamerSingleEliminationOff();
      turnRotamerPairEliminationOff();
    }
    boolean monteCarloTesting = properties.getBoolean("manybody-testing-mc", false);
    if (monteCarloTesting) {
      setMonteCarloTesting(true);
    }

    // Compute 4-body energies.
    String computeQuads = properties.getString("ro-compute4BodyEnergy");
    if (computeQuads != null) {
      this.compute4BodyEnergy = parseBoolean(computeQuads);
      logger.info(format(" (KEY) compute4BodyEnergy: %b", this.compute4BodyEnergy));
    }

    // Box / sliding window options.
    String direction = properties.getString("ro-direction");
    String boxDimensions = properties.getString("ro-boxDimensions");
    if (direction != null) {
      this.direction = Direction.valueOf(direction);
      logger.info(format(" (KEY) direction: %s", this.direction));
    }
    if (boxDimensions != null) {
      boxOpt.update(boxDimensions);
    }

    // Ensemble options.
    String ensembleNumber = properties.getString("ro-ensembleNumber");
    String ensembleEnergy = properties.getString("ro-ensembleEnergy");
    String ensembleBuffer = properties.getString("ro-ensembleBuffer");
    if (ensembleNumber != null) {
      this.ensembleNumber = parseInt(ensembleNumber);
      this.ensembleBuffer = 5.0;
      this.ensembleBufferStep = 0.1 * this.ensembleBuffer;
      logger.info(format(" (KEY) ensembleNumber: %d", this.ensembleNumber));
    }
    if (ensembleBuffer != null) {
      this.ensembleBuffer = parseDouble(ensembleBuffer);
      this.ensembleBufferStep = 0.1 * this.ensembleBuffer;
      logger.info(format(" (KEY) ensembleBuffer: %.2f", this.ensembleBuffer));
    }
    if (ensembleEnergy != null) {
      this.ensembleEnergy = parseDouble(ensembleEnergy);
      logger.info(format(" (KEY) ensembleEnergy: %.2f", this.ensembleEnergy));
    }

    // Nucleic Acid Options.
    String nucleicPruningFactor = properties.getString("ro-nucleicPruningFactor");
    String nucleicCorrectionThreshold = properties.getString("ro-nucleicCorrectionThreshold");
    String minimumNumberAcceptedNARotamers = properties.getString(
        "ro-minimumNumberAcceptedNARotamers");
    if (nucleicPruningFactor != null) {
      double value = parseDouble(nucleicPruningFactor);
      this.nucleicPruningFactor = (value >= 0 ? value : 1.0);
      this.nucleicPairsPruningFactor = (1.0 + value) / 2;
      logger.info(format(" (KEY) nucleicPruningFactor: %.2f", this.nucleicPruningFactor));
    }
    if (nucleicCorrectionThreshold != null) {
      double value = parseDouble(nucleicCorrectionThreshold);
      this.nucleicCorrectionThreshold = (value >= 0 ? value : 0);
      logger.info(
          format(" (KEY) nucleicCorrectionThreshold: %.2f", this.nucleicCorrectionThreshold));
    }
    if (minimumNumberAcceptedNARotamers != null) {
      int value = parseInt(minimumNumberAcceptedNARotamers);
      this.minNumberAcceptedNARotamers = (value > 0 ? value : 10);
      logger.info(
          format(" (KEY) minimumNumberAcceptedNARotamers: %d", this.minNumberAcceptedNARotamers));
    }

    // Superposition
    String superpositionThreshold = properties.getString("ro-superpositionThreshold");
    if (superpositionThreshold != null) {
      this.superpositionThreshold = parseDouble(superpositionThreshold);
      logger.info(format(" (KEY) superpositionThreshold: %.2f", this.superpositionThreshold));
    }

    // Multi-residue clash thresholds
    String multiResClashThreshold = properties.getString("ro-multiResClashThreshold");
    String multiResPairClashAddition = properties.getString("ro-multiResPairClashAddition");
    if (multiResClashThreshold != null) {
      this.multiResClashThreshold = parseDouble(multiResClashThreshold);
      logger.info(format(" (KEY) multiResClashThreshold: %.2f", this.multiResClashThreshold));
    }
    if (multiResPairClashAddition != null) {
      this.multiResPairClashAddn = parseDouble(multiResPairClashAddition);
      logger.info(format(" (KEY) multiResPairClashAddition: %.2f", this.multiResPairClashAddn));
    }

    String mcUseAll = properties.getString("ro-mcUseAll");
    String mcNoEnum = properties.getString("ro-debug-mcNoEnum");
    if (mcUseAll != null) {
      this.mcUseAll = parseBoolean(mcUseAll);
      logIfRank0(format(" (KEY) mcUseAll: %b", this.mcUseAll));
    }
    if (mcNoEnum != null) {
      this.mcNoEnum = parseBoolean(mcNoEnum);
      logIfRank0(format(" (KEY) debug-mcNoEnum: %b", this.mcNoEnum));
    }

    String propStr = properties.getString("ro-maxRotCheckDepth");
    int defaultMaxRotCheckDepth = 1;
    if (propStr != null) {
      try {
        maxRotCheckDepth = parseInt(propStr);
        if (maxRotCheckDepth > 3 || maxRotCheckDepth < 0) {
          throw new IllegalArgumentException(" ro-maxRotCheckDepth must be between 0-3 inclusive!");
        }
      } catch (Exception ex) {
        maxRotCheckDepth = defaultMaxRotCheckDepth;
        logger.warning(format(" Could not parse %s value %s as valid integer; defaulting to %d",
            "ro-maxRotCheckDepth", propStr, maxRotCheckDepth));
        logger.warning(format(" Exception: %s", ex));
      }
    } else {
      maxRotCheckDepth = defaultMaxRotCheckDepth;
    }

    setUpRestart();
  }

  /**
   * Set the writeEnergyRestart flag.
   *
   * @param writeEnergyRestart If true, write out an energy restart file.
   */
  public void setWriteEnergyRestart(boolean writeEnergyRestart) {
    this.writeEnergyRestart = writeEnergyRestart;
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
   * Checks if residue i is considered to be interacting with residue j, and thus has non-null
   * elements in the pair energies matrix.
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
   * Checks if residue i is considered to be interacting with residue j, that residue k is
   * interacting with either i or j, and thus i-j-k has non-null elements in the triple energy
   * matrix.
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
   * Checks the pair elimination array to see if this permutation has been eliminated.
   *
   * @param i           Residue number
   * @param ri          Rotamer number
   * @param currentRots Array of current rotamers indeces.
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
   * Computes the environment/backbone energy, defined as energy with all sidechains under
   * consideration turned off in their 0th rotamer.
   *
   * @param residues Residues under optimization.
   * @return Backbone energy Eenv/bb.
   * @throws ArithmeticException If an exception in calculating energy is found.
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
   * Uses existing backbone, self, 2-Body, and 3-body energies from rotamerEnergies() to calculate an
   * approximate energy for a rotamer permutation.
   *
   * @param residues Current window of optimization.
   * @param rotamers Set of rotamers to calculate an approximate energy for.
   * @param print    Verbosity flag
   * @return Approximate permutation energy (backbone + selfs + pairs + trimers).
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
   * Calculates the energy at the current state.
   *
   * @param resArray Array of residues in current energy term.
   * @return Energy of the current state.
   */
  public double currentEnergy(Residue[] resArray) throws ArithmeticException {
    return currentEnergy(Arrays.asList(resArray));
  }

  /**
   * Wrapper intended for use with RotamerMatrixMC.
   *
   * @param resList Reside list.
   * @return Returns the current energy.
   * @throws ArithmeticException Thrown if there is an arithmetic exception computing the
   *                             energy.
   */
  public double currentEnergyWrapper(List<Residue> resList) throws ArithmeticException {
    return currentEnergy(resList);
  }

  /**
   * Forces the use of a ForceFieldEnergyOpenMM's underlying ForceFieldEnergy.
   *
   * @return Current potential energy as calculated by FFX reference platform.
   */
  public double currentFFXPE() {
    if (x == null) {
      int n = potential.getNumberOfVariables();
      x = new double[n];
    }
    potential.getCoordinates(x);
    return ((OpenMMEnergy) potential).energyFFX(x, false);
  }

  /**
   * Utility method for formatting energies, using 16 spaces with 8 digits of precision.
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
  }

  public double getApproximate() {
    return approximateEnergy;
  }

  public double getBackboneEnergy() {
    return eE.getBackboneEnergy();
  }

  public EliminatedRotamers getEliminatedRotamers() {
    return eR;
  }

  public EnergyExpansion getEnergyExpansion() {
    return eE;
  }

  /**
   * getEnsemble.
   *
   * @return a {@link java.util.List} object.
   */
  public List<ResidueState[]> getEnsemble() {
    if (ensembleStates == null) {
      return null;
    } else {
      List<ResidueState[]> states = new ArrayList<>(ensembleStates.size());
      ensembleStates.forEach((es) -> states.add(es.val()));
      return states;
    }
  }

  /**
   * setEnsemble.
   *
   * @param ensemble a int.
   */
  public void setEnsemble(int ensemble) {
    setEnsemble(ensemble, 5.0);
  }

  /**
   * Return the residue list.
   *
   * @return a {@link java.util.List} object.
   */
  public List<Residue> getResidues() {
    return residueList;
  }

  /**
   * Set the residue list.
   *
   * @param residueList a {@link java.util.List} object.
   */
  public void setResidues(List<Residue> residueList) {
    this.residueList = residueList;
  }

  /**
   * Init fraction array
   *
   * @param residueList a {@link java.util.List} object.
   */
  public void initFraction(List<Residue> residueList) {
    fraction = new double[residueList.size()][56];
    genZ = true;
  }

  /**
   * Returns the restart file.
   *
   * @return energyRestartFile File with saved side-chain energies.
   */
  public File getRestartFile() {
    return energyRestartFile;
  }

  public double goldsteinPairSumOverK(Residue[] residues, int lb, int ub, int i, int riA, int riB,
                                      int j, int rjC, int rjD, List<Residue> blockedResidues, int[] possK) {
    double sumOverK = 0.0;

    for (int indK = lb; indK <= ub; indK++) {
      int k = possK[indK];
      double minForResK = Double.MAX_VALUE;
      Residue resk = residues[k];
      Rotamer[] rotk = resk.getRotamers();
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
        double currentResK =
            eE.get2Body(i, riA, k, rk) - eE.get2Body(i, riB, k, rk) + eE.get2Body(j, rjC, k, rk)
                - eE.get2Body(j, rjD, k, rk);
        // Include 3-body effects.
        if (threeBodyTerm) {
          double sumOverL = (eE.get3Body(residues, i, riA, j, rjC, k, rk) - eE.get3Body(residues, i,
              riB, j, rjD, k, rk));
          // Loop over a 4th residue l.
          int[] nK = bidiResNeighbors[k];
          IntStream lStream = IntStream.concat(Arrays.stream(possK), Arrays.stream(nK));
          int[] possL = lStream.filter(l -> (l != i && l != j && l != k)).sorted().distinct()
              .toArray();

          for (int l : possL) {
            if (l == k || l == i || l == j) {
              continue;
            }
            Residue residuel = residues[l];
            Rotamer[] rotamersl = residuel.getRotamers();
            int nrl = rotamersl.length;
            int rlEvaluations = 0;
            double minForResL = Double.MAX_VALUE;
            // Loop over rotamers for residue l.
            for (int rl = 0; rl < nrl; rl++) {
              // If not a part of valid phase space for riA/rjC, continue.
              if (eR.check(l, rl) || eR.check(k, rk, l, rl) || eR.check(i, riA, l, rl) || eR.check(j,
                  rjC, l, rl)) {
                // Not implemented: check(i, riA, j, rjC, l, rl) || check(i, riA, k, rk, l, rl) ||
                // check(j, rjC, k, rk, l, rl) || check(i, riA, j, rjC, k, rk, l, rl)
                continue;
              }
              if (eR.check(i, riB, l, rl) || eR.check(j, rjD, l, rl)) {
                // Not implemented: check(i, riB, j, rjD, l, rl) || check(i, riB, k, rk, l, rl) ||
                // check(j, rjD, k, rk, l, rl) || check(i, riB, j, rjD, k, rk, l, rl)
                blockedResidues.add(residuel);
                return Double.NaN;
              }
              rlEvaluations++;
              double e =
                  eE.get3Body(residues, i, riA, k, rk, l, rl) - eE.get3Body(residues, i, riB, k, rk,
                      l, rl) + eE.get3Body(residues, j, rjC, k, rk, l, rl) - eE.get3Body(residues, j,
                      rjD, k, rk, l, rl);
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

  public void logIfRank0(String msg) {
    if (rank0) {
      logger.info(msg);
    }
  }

  public void logIfRank0(String msg, Level level) {
    if (rank0) {
      logger.log(level, msg);
    }
  }

  /**
   * Perform the rotamer optimization using the specified algorithm.
   *
   * @param algorithm a {@link RotamerOptimization.Algorithm} object.
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

      CompositeConfiguration properties = molecularAssembly.getProperties();
      boolean ignoreNA = properties.getBoolean("ignoreNA", false);
      if (ignoreNA) {
        logger.info(" Ignoring nucleic acids.");
      }

      logger.info(format("\n Rotamer Library:     %s", library.getLibrary()));
      logger.info(format(" Algorithm:           %s", algorithm));
      logger.info(format(" Goldstein Criteria:  %b", useGoldstein));
      logger.info(format(" Three-Body Energies: %b\n", threeBodyTerm));

      /*
       * Collect all residues in the MolecularAssembly. Use all Residues with Rotamers.
       */
      allResiduesList = new ArrayList<>();
      // An array of polymers from the MolecularAssembly.
      Polymer[] polymers = molecularAssembly.getChains();
      for (Polymer polymer : polymers) {
        List<Residue> current = polymer.getResidues();
        for (Residue residuej : current) {
          residuej.setRotamers(library);
          if (residuej.getRotamers() != null) {
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
        List<Residue> onlyAA = new ArrayList<>();
        for (Residue res : residueList) {
          if (res.getResidueType() != Residue.ResidueType.NA) {
            onlyAA.add(res);
          }
        }
        residueList = onlyAA;
      }

      RotamerLibrary.initializeDefaultAtomicCoordinates(
          molecularAssembly.getChains()); // for NA only

      nAllResidues = allResiduesList.size();
      allResiduesArray = allResiduesList.toArray(new Residue[nAllResidues]);

      /*
       * Distance matrix is used to add residues to the sliding window based on distance cutoff,
       * and for cutoff distances.
       *
       * The memory and compute overhead can be a problem for some very large structures.
       */
      if (distance > 0) {
        dM = new DistanceMatrix(molecularAssembly, algorithmListener, allResiduesArray,
            allResiduesList, distanceMethod, distance, twoBodyCutoffDist, threeBodyCutoffDist);
      }

      if (residueList != null) {

        // Allocate memory for storing optimal rotamers.
        optimum = new int[residueList.size()];

        done = false;
        terminate = false;

        switch (algorithm) {
          case INDEPENDENT -> e = independent(residueList);
          case BRUTE_FORCE -> e = bruteForce(residueList);
          case ALL -> {
            e = globalOptimization(residueList);
            arraycopy(optimumSubset, 0, optimum, 0, residueList.size());
          }
          case WINDOW -> {
            if (genZ) {
              e = slidingWindowCentered(residueList);
            } else {
              e = slidingWindowOptimization(residueList, windowSize, increment, revert, distance,
                  direction, -1);
            }
          }
          case BOX -> e = boxOpt.boxOptimization(residueList);
          default -> {
          }
        }
        terminate = false;
        done = true;
      }
    } catch (Exception exception) {
      exception.printStackTrace();
    } finally {
      try {
        if (energyWriter != null) {
          energyWriter.close();
        }
      } catch (IOException ex) {
        logger.severe(" Exception in closing buffered energy writer.");
      }
    }

    return e;
  }

  /**
   * A brute-force global optimization over side-chain rotamers using a recursive algorithm.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param residues          an array of {@link ffx.potential.bonded.Residue} objects.
   * @param i                 a int.
   * @param lowEnergy         a double.
   * @param optimum           the minimum energy rotamer indices for each residue.
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
    Rotamer[] rotamers = current.getRotamers();
    int lenri = rotamers.length;
    double currentEnergy = Double.MAX_VALUE;
    List<Residue> resList = Arrays.asList(residues);
    if (i < nResidues - 1) {
      // As long as there are more residues, continue the recursion for each rotamer of the current
      // residue.
      int minRot = -1;
      for (int ri = 0; ri < lenri; ri++) {
        applyRotamer(current, rotamers[ri]);
        double rotEnergy = rotamerOptimization(molecularAssembly, residues, i + 1, lowEnergy,
            optimum);
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
      /*
       At the end of the recursion, compute the potential energy for each rotamer of the final
       residue. If a lower potential energy is discovered, the rotamers of each residue will be
       collected as the recursion returns up the chain.
      */
      int[] rotArray = copyOf(optimum, nResidues);
      for (int ri = 0; ri < lenri; ri++) {
        applyRotamer(current, rotamers[ri]);
        double rotEnergy = Double.NaN;
        try {
          rotArray[nResidues - 1] = ri;
          rotEnergy = currentEnergy(resList) + eE.getTotalRotamerPhBias(resList, rotArray, pH, pHRestraint);
          logger.info(format(" %d Energy: %s", ++evaluatedPermutations, formatEnergy(rotEnergy)));
        } catch (ArithmeticException ex) {
          logger.info(
              format(" %d Energy set to NaN (unreasonable conformation)", ++evaluatedPermutations));
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
   * Set the K for the harmonic pH restraint
   *
   * @param pHRestraint KpH
   */
  public void setPHRestraint(double pHRestraint) {
    this.pHRestraint = pHRestraint;
  }

  /**
   * Set the environment pH
   *
   * @param pH
   */
  public void setpH(double pH) {
    this.pH = pH;
  }

  /**
   * Sets to recompute self energies at a different pH using an energy restart file
   *
   * @param recomputeSelf
   */
  public void setRecomputeSelf(boolean recomputeSelf) {
    this.recomputeSelf = recomputeSelf;
  }

  /**
   * Return the K in the harmonic pH restraint
   *
   * @return double KpH
   */
  public double getPHRestraint() {
    return pHRestraint;
  }

  /**
   * Get the enviroment pH
   *
   * @return double pH
   */
  public double getPH() {
    return pH;
  }

  /**
   * Sets the approximate dimensions of boxes, over-riding numXYZBoxes in determining box size.
   * Rounds box size up and number of boxes down to get a whole number of boxes along each axis.
   *
   * @param approxBoxLength Optional box dimensions parameter (Angstroms).
   */
  public void setApproxBoxLength(double approxBoxLength) {
    boxOpt.approxBoxLength = approxBoxLength;
  }

  /**
   * Sets the amount of overlap between adjacent boxes for box optimization.
   *
   * @param boxBorderSize Box overlap in Angstroms.
   */
  public void setBoxBorderSize(double boxBorderSize) {
    boxOpt.cellBorderSize = boxBorderSize;
  }

  /**
   * Set the ending box index.
   *
   * @param boxEnd a int.
   */
  public void setBoxEnd(int boxEnd) {
    // Is -1 if boxes run to completion.
    boxOpt.cellEnd = boxEnd;
  }

  /**
   * Sets behavior for how Residues are added to boxOptCells; 1 uses just reference atom (C alpha for
   * protein, N1/9 for nucleic acids), 2 uses any atom, 3 uses any atom in any rotamer.
   *
   * @param boxInclusionCriterion Criterion to use
   */
  public void setBoxInclusionCriterion(int boxInclusionCriterion) {
    boxOpt.boxInclusionCriterion = boxInclusionCriterion;
  }

  /**
   * Set the starting box index.
   *
   * @param boxStart a int.
   */
  public void setBoxStart(int boxStart) {
    boxOpt.cellStart = boxStart;
  }

  /**
   * Sets titratble residues to be the center of the box.
   *
   * @param titrationBoxes a boolean.
   */
  public void setTitrationBoxes(boolean titrationBoxes) {
    boxOpt.titrationBoxes = titrationBoxes;
  }

  /**
   * Sets the size around the titratable residues.
   *
   * @param titrationBoxSize double size of the titration box.
   */
  public void setTitrationBoxSize(double titrationBoxSize) {
    boxOpt.titrationBoxSize = titrationBoxSize;
  }

  /**
   * setCoordinatesToEnsemble.
   *
   * @param ensnum a int.
   */
  public void setCoordinatesToEnsemble(int ensnum) {
    if (ensembleStates != null && !ensembleStates.isEmpty()) {
      ensnum %= ensembleStates.size();
      ResidueState.revertAllCoordinates(residueList, ensembleStates.get(ensnum).val());
    } else {
      throw new IllegalArgumentException(" Ensemble states not initialized!");
    }
  }

  /**
   * Sets the decompose-original flag.
   *
   * @param decomposeOriginal If true, decompose the energy of the structure without optimizing.
   */
  public void setDecomposeOriginal(boolean decomposeOriginal) {
    this.decomposeOriginal = decomposeOriginal;
  }

  /**
   * Set the optimization direction to forward or backward.
   *
   * @param direction a {@link RotamerOptimization.Direction} object.
   */
  public void setDirection(Direction direction) {
    this.direction = direction;
  }

  /**
   * Set the cut-off distance for inclusion of residues in sliding box and window methods.
   *
   * @param distance a double.
   */
  public void setDistanceCutoff(double distance) {
    this.distance = distance;
  }

  /**
   * Setter for the field <code>energyRestartFile</code>.
   *
   * @param file a {@link java.io.File} object.
   */
  public void setEnergyRestartFile(File file) {
    loadEnergyRestart = true;
    energyRestartFile = file;
  }

  /**
   * setEnsemble.
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
   * Set the residue increment for sliding window.
   *
   * @param increment a int.
   */
  public void setIncrement(int increment) {
    this.increment = increment;
  }

  /**
   * Control the depth of self-consistency checking with a rotamer is eliminated.
   *
   * @param maxRotCheckDepth a int.
   */
  public void setMaxRotCheckDepth(int maxRotCheckDepth) {
    this.maxRotCheckDepth = maxRotCheckDepth;
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
      logger.warning(
          "\n Minimum number of accepted NA rotamers must be a positive integer.\n Setting to default value 10.\n");
      this.minNumberAcceptedNARotamers = 10;
    }
  }

  /**
   * Sets the option to use a number of Monte Carlo steps for final optimization.
   *
   * @param monteCarlo If Monte Carlo is to be considered
   * @param nMCsteps   Number of steps to be taken
   */
  public void setMonteCarlo(boolean monteCarlo, int nMCsteps) {
    this.monteCarlo = monteCarlo;
    this.nMCSteps = nMCsteps;
  }

  /**
   * Sets the monteCarloTesting boolean in RotamerOptimization.java to true or false. This should
   * only be set to true when monte carlo is being tested through the ManyBodyTest script. When true,
   * the method sets a seed for the pseudo-random number generator and allows the monte carlo rotamer
   * optimization to be deterministic.
   *
   * @param bool True or false.
   */
  public void setMonteCarloTesting(boolean bool) {
    this.monteCarloTesting = bool;
  }

  /**
   * The nucleic acid correction threshold.
   *
   * @param nucleicCorrectionThreshold a double.
   */
  public void setNucleicCorrectionThreshold(double nucleicCorrectionThreshold) {
    if (nucleicCorrectionThreshold >= 0) {
      this.nucleicCorrectionThreshold = nucleicCorrectionThreshold;
    } else {
      logger.warning(
          "\n Correction threshold must be >= 0. Setting to default of 0 (threshold inactive).\n");
      this.nucleicCorrectionThreshold = 0;
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
   * Sets the number of boxes in the x, y, and z axes if the box optimization is to be carried out.
   *
   * @param numXYZBoxes Int[3] of number of boxes in x, y, z.
   */
  public void setNumXYZBoxes(int[] numXYZBoxes) {
    arraycopy(numXYZBoxes, 0, boxOpt.numXYZCells, 0, boxOpt.numXYZCells.length);
  }

  /**
   * Setter for the field <code>pairClashThreshold</code>.
   *
   * @param pairClashThreshold a double.
   */
  public void setPairClashThreshold(double pairClashThreshold) {
    this.pairClashThreshold = pairClashThreshold;
  }

  /**
   * Sets whether rotamer optimization should print out any files, or act solely to optimize a
   * structure in memory.
   *
   * @param printFiles a boolean.
   */
  public void setPrintFiles(boolean printFiles) {
    this.printFiles = printFiles;
  }

  /**
   * Sets level of pruning: 0 for fully off, 1 for only singles, 2 for single and pair pruning.
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
   * Accepts a list of residues but throws out null residues. Used by the -lR flag.
   *
   * @param residues a {@link java.util.List} object.
   */
  public void setResiduesIgnoreNull(List<Residue> residues) {
    residueList = new ArrayList<>();
    logger.fine(" Optimizing these residues: ");
    for (Residue r : residues) {
      if (r.getRotamers() != null) {
        residueList.add(r);
        logger.fine(format("\t%s", r));
      } else {
        logger.fine(format(" not \t%s", r));
      }
    }
  }

  /**
   * Set the algorithm to revert to starting coordinates if the energy increases.
   *
   * @param revert a boolean.
   */
  public void setRevert(boolean revert) {
    this.revert = revert;
  }

  public void setRotamerLibrary(RotamerLibrary lib) {
    library = lib;
  }

  /**
   * setSingletonClashThreshold.
   *
   * @param singletonClashThreshold a double.
   */
  public void setSingletonClashThreshold(double singletonClashThreshold) {
    this.clashThreshold = singletonClashThreshold;
  }

  /**
   * Setter for the field <code>superpositionThreshold</code>.
   *
   * @param superpositionThreshold a double.
   */
  public void setSuperpositionThreshold(double superpositionThreshold) {
    this.superpositionThreshold = superpositionThreshold;
  }

  /**
   * Sets the threeBodyCutoffDist. All three-body energies where the rotamers have a separation
   * distance larger than the cutoff are set to 0.
   *
   * @param threeBodyCutoffDist Separation distance at which the interaction of three side-chains
   *                            is assumed to have an energy of 0.
   */
  public void setThreeBodyCutoff(double threeBodyCutoffDist) {
    this.threeBodyCutoffDist = threeBodyCutoffDist;
    if (this.threeBodyCutoffDist < 0) {
      logger.info("Warning: threeBodyCutoffDist should not be less than 0.");
    }
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
   * Sets the twoBodyCutoffDist. All two-body energies where the rotamers have a separation distance
   * larger than the cutoff are set to 0.
   *
   * @param twoBodyCutoffDist Separation distance at which the interaction of two side-chains is
   *                          assumed to have an energy of 0.
   */
  public void setTwoBodyCutoff(double twoBodyCutoffDist) {
    this.twoBodyCutoffDist = twoBodyCutoffDist;
    if (this.twoBodyCutoffDist < 0) {
      logger.info("Warning: threeBodyCutoffDist should not be less than 0.");
    }
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
   * Setter for the field <code>windowSize</code>.
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
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    if (eR != null) {
      return eR.toString();
    } else {
      return null;
    }
  }

  public void turnOffAllResidues(Residue[] residues) {
    if (residues == null) {
      return;
    }
    for (Residue residue : residues) {
      turnOffResidue(residue);
    }
  }

  public void turnOffResidue(Residue residue) {
    turnOffAtoms(residue);
    applyDefaultRotamer(residue);
  }

  public void turnOnAllResidues(Residue[] residues) {
    if (residues == null) {
      return;
    }
    for (Residue residue : residues) {
      turnOnAtoms(residue);
    }
  }

  public void turnOnResidue(Residue residue, int ri) {
    turnOnAtoms(residue);
    Rotamer[] rotamers = residue.getRotamers();
    applyRotamer(residue, rotamers[ri]);
  }

  /**
   * ONLY FOR UNIT TESTING. Sets a boolean to turn the pair elimination criteria off.
   */
  public void turnRotamerPairEliminationOff() {
    logger.info(" Turning off pair eliminations.");
    pairEliminationOn = false;
  }

  /**
   * ONLY FOR UNIT TESTING. Sets a boolean to turn the self elimination criteria off.
   */
  public void turnRotamerSingleEliminationOff() {
    logger.info(" Turning off single eliminations.");
    selfEliminationOn = false;
  }

  /**
   * Recursive brute-force method which uses single, pair, and potentially trimer energies to
   * calculate an optimum set of rotamers.
   *
   * @param residues        Optimization window
   * @param i               Current residue in the recursion.
   * @param lowEnergy       Minimum energy yet found by the recursion.
   * @param optimum         Optimum rotamer set yet found by the recursion.
   * @param currentRotamers Rotamer permutation under investigation.
   * @return Minimum energy found under this node in the recursion.
   */
  private double decomposedRotamerOptimization(Residue[] residues, int i, double lowEnergy,
                                               int[] optimum, int[] currentRotamers) {

    // This is the initialization condition.
    if (i == 0) {
      evaluatedPermutations = 0;
    }

    int nResidues = residues.length;
    Residue current = residues[i];
    Rotamer[] rotamers = current.getRotamers();
    int lenri = rotamers.length;
    double currentEnergy = Double.MAX_VALUE;
    if (i < nResidues - 1) {
      // As long as there are more residues, continue the recursion for each rotamer of the current
      // residue.
      for (int ri = 0; ri < lenri; ri++) {
        if (eR.check(i, ri)) {
          continue;
        }
        currentRotamers[i] = ri;
        double rotEnergy = decomposedRotamerOptimization(residues, i + 1, lowEnergy, optimum,
            currentRotamers);
        if (rotEnergy < lowEnergy) {
          lowEnergy = rotEnergy;
        }
        if (rotEnergy < currentEnergy) {
          currentEnergy = rotEnergy;
        }
      }
    } else {
      // At the end of the recursion, compute the potential energy for each rotamer of the final
      // residue and update optimum[].
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
          /*
           Because we print the rotamer set immediately on finding a
           more optimal structure, we have to reset the entire length
           of optimum instead of lazily doing it on the way down.
          */
          arraycopy(currentRotamers, 0, optimum, 0, optimum.length);
          if (evaluatedPermutations > 1) {
            logger.info(
                format(" Minimum energy update: %f < %f, permutation %d", rotEnergy, lowEnergy,
                    evaluatedPermutations));
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
   * Runs Monte Carlo side chain optimization using the rotamer energy matrix and potentially some
   * information from dead-end or Goldstein elimination. The useAllElims variable should be set false
   * if detailed balance is to be maintained. At present, no support for ensembles.
   *
   * @param residues      Optimization window
   * @param optimum       Array to store optimum rotamers
   * @param initialRots   Array with starting rotamers
   * @param maxIters      Number of MC steps to run
   * @param randomizeRots Scramble initialRots
   * @param useAllElims   Use pair/triple elimination information
   * @return Lowest energy found
   */
  private double rotamerOptimizationMC(Residue[] residues, int[] optimum, int[] initialRots,
                                       int maxIters, boolean randomizeRots, boolean useAllElims) {

    long initTime = -System.nanoTime();
    if (randomizeRots) {
      randomizeRotamers(initialRots, residues, true);
    }

    int nRes = residues.length;
    arraycopy(initialRots, 0, optimum, 0, nRes);
    assert optimum.length == nRes;
    assert initialRots.length == nRes;

    RotamerMatrixMC rmc = new RotamerMatrixMC(initialRots, residues, useForceFieldEnergy, this);
    rmc.setTemperature(temperature);
    RotamerMatrixMove rmove = new RotamerMatrixMove(useAllElims, initialRots, residues, this, eR,
        monteCarloTesting);
    List<MCMove> rmList = new ArrayList<>(1);
    rmList.add(rmove);

    double initialEnergy = computeEnergy(residues, initialRots, false);
    double optimumEnergy = initialEnergy;
    double currentEnergy = initialEnergy;

    int nAccept = 0;
    logIfRank0(
        format(" Beginning %d iterations of Monte Carlo search " + "starting from energy %10.6f",
            maxIters, initialEnergy));

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
    logIfRank0(
        format(" %d steps of DEE-MC completed in %10.6f seconds", maxIters, (initTime * 1.0E-9)));
    logIfRank0(format(" Number of steps accepted: %d for %10.6f of total", nAccept, fractAccept));
    logIfRank0(format(" Lowest energy found: %10.6f kcal/mol", optimumEnergy));
    logIfRank0(format(" Final energy found: %10.6f kcal/mol", currentEnergy));

    return optimumEnergy;
  }

  /**
   * Scrambles an array of rotamers.
   *
   * @param rotamers    Array of rotamers.
   * @param residues    Array of residues.
   * @param useAllElims If trye, use all elliminations.
   */
  private void randomizeRotamers(int[] rotamers, Residue[] residues, boolean useAllElims) {
    int nRes = rotamers.length;
    for (int i = 0; i < nRes; i++) {
      Rotamer[] rotsi = residues[i].getRotamers();
      int lenri = rotsi.length;
      List<Integer> allowedRots = new ArrayList<>(lenri);

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
   * A global optimization over side-chain rotamers using a recursive algorithm and information about
   * eliminated rotamers, rotamer pairs and rotamer triples
   *
   * @param molecularAssembly   MolecularAssmebly to use.
   * @param residues            Array of residues.
   * @param i                   Current permutation.
   * @param currentRotamers     Current array of rotamers.
   * @param lowEnergy           Lowest energy found.
   * @param optimum             Optimum set of rotamers.
   * @param permutationEnergies Energies of visited permutations or null.
   * @return current energy.
   */
  private double rotamerOptimizationDEE(MolecularAssembly molecularAssembly, Residue[] residues,
                                        int i, int[] currentRotamers, double lowEnergy, int[] optimum, double[] permutationEnergies) {
    // This is the initialization condition.
    if (i == 0) {
      evaluatedPermutations = 0;
    }

    int nResidues = residues.length;
    Residue residuei = residues[i];
    Rotamer[] rotamersi = residuei.getRotamers();
    int lenri = rotamersi.length;
    double currentEnergy = Double.MAX_VALUE;
    List<Residue> resList = Arrays.asList(residues);

    // As long as there are more residues, continue the recursion for each rotamer of the current
    // residue.
    if (i < nResidues - 1) {
      // Loop over rotamers of residue i.
      for (int ri = 0; ri < lenri; ri++) {
        // Check if rotamer ri has been eliminated by DEE.
        if (eR.check(i, ri)) {
          continue;
        }
        // Check if rotamer ri has been eliminated by an upstream rotamer (any residue's rotamer
        // from j = 0 .. i-1).
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

      /*
       At the end of the recursion, compute the potential energy for each rotamer of the final
       residue. If a lower potential energy is discovered, the rotamers of each residue will be
       collected as the recursion returns up the chain.
      */
      for (int ri = 0; ri < lenri; ri++) {
        // Check if rotamer ri has been eliminated by DEE.
        if (eR.check(i, ri)) {
          continue;
        }
        currentRotamers[i] = ri;
        // Check if rotamer ri has been eliminated by an upstream rotamer (any residue's rotamer
        // from 0 .. i-1).
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
        if (useForceFieldEnergy) {
          double amoebaEnergy = Double.NaN;
          try {
            // Add the rotamer pH bias to the force field energy.
            amoebaEnergy =
                currentEnergy(resList) + eE.getTotalRotamerPhBias(resList, currentRotamers, pH, pHRestraint);
          } catch (ArithmeticException ex) {
            logger.warning(
                format(" Exception %s in calculating full AMOEBA energy for permutation %d", ex,
                    evaluatedPermutations));
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
          if (rank0 && printFiles) {
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
              bw.write("ENDMDL");
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

        if (useForceFieldEnergy) {
          // Log current results
          logIfRank0(format(" %6e AMOEBA: %12.4f 3-Body: %12.4f Neglected: %12.4f (%12.4f)",
              (double) evaluatedPermutations, comparisonEnergy, approximateEnergy,
              comparisonEnergy - approximateEnergy, lowEnergy));
        } else {
          logIfRank0(
              format(" %12s %25f %25f", evaluatedPermutations, approximateEnergy, lowEnergy));
        }
      }

      ensembleStates.sort(null);
    }
    return currentEnergy;
  }

  /**
   * A global optimization over side-chain rotamers using a recursive algorithm and information about
   * eliminated rotamers, rotamer pairs and rotamer triples.
   *
   * @param residues        Residue array.
   * @param i               Current number of permutations.
   * @param currentRotamers Current rotamer list.
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
        logIfRank0(
            format(" The permutations have reached %10.4e.", (double) evaluatedPermutationsPrint));
        evaluatedPermutationsPrint *= 10;
      }
    }

    int nResidues = residues.length;
    Residue residuei = residues[i];
    Rotamer[] rotamersi = residuei.getRotamers();
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
      // At the end of the recursion, check each rotamer of the final residue.
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
   * @param residues            Array of residues.
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
    Rotamer[] rotamersi = residuei.getRotamers();
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
        dryRunForEnsemble(residues, i + 1, currentRotamers, gmecEnergy, permutationEnergies,
            permutations);
      }
    } else {
      // At the end of the recursion, check each rotamer of the final residue.
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
   * A global optimization over side-chain rotamers using a recursive algorithm and information about
   * eliminated rotamers, rotamer pairs and rotamer triples.
   *
   * @param residues        Residue array.
   * @param i               Current number of permutations.
   * @param currentRotamers Current rotamer list.
   */
  public void partitionFunction(Residue[] residues, int i, int[] currentRotamers) throws Exception {
    // This is the initialization condition.
    final double LOG10 = log(10.0);
    final double beta = 1.0 / (Constants.kB * temperature);

    if (i == 0) {
      totalBoltzmann = 0;
      evaluatedPermutations = 0;
      evaluatedPermutationsPrint = 1000;
    }

    if (evaluatedPermutations >= evaluatedPermutationsPrint) {
      if (evaluatedPermutations % evaluatedPermutationsPrint == 0) {
        logIfRank0(format(" The permutations have reached %10.4e.", (double) evaluatedPermutationsPrint));
        evaluatedPermutationsPrint *= 10;
      }
    }

    int nResidues = residues.length;
    Residue residuei = residues[i];
    Rotamer[] rotamersi = residuei.getRotamers();
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
        partitionFunction(residues, i + 1, currentRotamers);
      }
    } else {
      // At the end of the recursion, check each rotamer of the final residue.
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

          energyRegion.init(eE, residues, currentRotamers, threeBodyTerm);
          parallelTeam.execute(energyRegion);
          double selfEnergy = energyRegion.getSelf();
          // Recompute the self energy from a restart file run at pH 7.0
          if (recomputeSelf) {
            int count = 0;
            for (Residue residue : residues) {
              double bias7 = 0;
              double biasCurrent = 0;
              Rotamer[] rotamers = residue.getRotamers();
              int currentRotamer = currentRotamers[count];
              switch (rotamers[currentRotamer].getName()) {
                case "HIE" -> {
                  bias7 = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.HIStoHIE.pKa - 7)) -
                      TitrationUtils.Titration.HIStoHIE.freeEnergyDiff;
                  biasCurrent = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.HIStoHIE.pKa - pH)) -
                      TitrationUtils.Titration.HIStoHIE.freeEnergyDiff;
                }
                case "HID" -> {
                  bias7 = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.HIStoHID.pKa - 7)) -
                      TitrationUtils.Titration.HIStoHID.freeEnergyDiff;
                  biasCurrent = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.HIStoHID.pKa - pH)) -
                      TitrationUtils.Titration.HIStoHID.freeEnergyDiff;
                }
                case "ASP" -> {
                  bias7 = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.ASHtoASP.pKa - 7)) -
                      TitrationUtils.Titration.ASHtoASP.freeEnergyDiff;
                  biasCurrent = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.ASHtoASP.pKa - pH)) -
                      TitrationUtils.Titration.ASHtoASP.freeEnergyDiff;
                }
                case "GLU" -> {
                  bias7 = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.GLHtoGLU.pKa - 7)) -
                      TitrationUtils.Titration.GLHtoGLU.freeEnergyDiff;
                  biasCurrent = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.GLHtoGLU.pKa - pH)) -
                      TitrationUtils.Titration.GLHtoGLU.freeEnergyDiff;
                }
                case "LYD" -> {
                  bias7 = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.LYStoLYD.pKa - 7)) -
                      TitrationUtils.Titration.LYStoLYD.freeEnergyDiff;
                  biasCurrent = (LOG10 * Constants.R * temperature * (TitrationUtils.Titration.LYStoLYD.pKa - pH)) -
                      TitrationUtils.Titration.LYStoLYD.freeEnergyDiff;
                }
                default -> {
                }
              }
              selfEnergy = selfEnergy - bias7 + biasCurrent;
              count += 1;
            }
          }

          // Calculate the total energy of a permutation/conformation
          double totalEnergy = eE.getBackboneEnergy() + selfEnergy +
              energyRegion.getTwoBody() + energyRegion.getThreeBody();

          // Set a reference energy to evaluate all follow energies against for the Boltzmann calculations to avoid Nan/Inf errors
          if (evaluatedPermutations == 1) {
            refEnergy = totalEnergy;
          }
          double boltzmannWeight = exp(-beta * (totalEnergy - refEnergy));

          // Collect Boltzmann weight for every rotamer for residues included in the optimization
          for (int res = 0; res < residues.length; res++) {
            int currentRotamer = currentRotamers[res];
            populationBoltzmann[res][currentRotamer] += boltzmannWeight;
          }

          // Sum Boltzmann of all permutations
          totalBoltzmann += boltzmannWeight;
        }
      }
    }
  }

  /**
   * Return reference energy for partition function boltzmann weights
   *
   * @return ref energy
   */
  public double getRefEnergy() {
    return refEnergy;
  }

  /**
   * Return the total boltzmann weight for an ensemble
   *
   * @return total boltzmann
   */
  public double getTotalBoltzmann() {
    return totalBoltzmann;
  }

  /**
   * Return the ensemble average of protonated rotamers for all titratable sites
   *
   * @return fraction of protonated residues
   */
  public double[][] getFraction() {
    return fraction;
  }

  /**
   * Return the Protonated Boltzmann for all sites
   *
   * @return fraction of protonated residues
   */
  public double[][] getPopulationBoltzmann() {
    return populationBoltzmann;
  }

  /**
   * Calculate population of each rotamer for residues in the system
   *
   * @param residues        residue array
   * @param i               int
   * @param currentRotamers empty array
   * @throws Exception too many permutations to continue
   */
  public void getFractions(Residue[] residues, int i, int[] currentRotamers) throws Exception {
    populationBoltzmann = new double[residues.length][56];
    if (!usingBoxOptimization) {
      partitionFunction(residues, i, currentRotamers);
      optimum = new int[residues.length];
      for (int m = 0; m < fraction.length; m++) {
        for (int n = 0; n < 56; n++) {
          fraction[m][n] = populationBoltzmann[m][n] / totalBoltzmann;
          if (n > 0 && fraction[m][n] > fraction[m][n - 1]) {
            optimum[m] = n;
          } else if (n == 0) {
            optimum[m] = n;
          }
        }
        Rotamer highestPopRot = residues[m].getRotamers()[optimum[m]];
        RotamerLibrary.applyRotamer(residues[m], highestPopRot);
      }
    }
    logger.info("\n   Total permutations evaluated: " + evaluatedPermutations + "\n");
  }

  /**
   * Return population of each rotamer for residues in the system
   *
   * @param residues        residue array
   * @param i               int
   * @param currentRotamers empty array
   * @throws Exception too many permutations to continue
   */
  public void getFractions(Residue[] residues, int i, int[] currentRotamers, boolean usingBoxOptimization) throws Exception {
    double[][] fractionSubset = new double[residues.length][56];
    populationBoltzmann = new double[residues.length][56];
    partitionFunction(residues, i, currentRotamers);
    optimumSubset = new int[residues.length];
    for (int m = 0; m < fractionSubset.length; m++) {
      for (int n = 0; n < 56; n++) {
        fractionSubset[m][n] = populationBoltzmann[m][n] / totalBoltzmann;
        if (n > 0 && fractionSubset[m][n] > fractionSubset[m][n - 1]) {
          optimumSubset[m] = n;
        } else if (n == 0) {
          optimumSubset[m] = n;
        }
      }
      Rotamer highestPopRot = residues[m].getRotamers()[optimumSubset[m]];
      RotamerLibrary.applyRotamer(residues[m], highestPopRot);
      Residue residue = residues[m];
      int index = residueList.indexOf(residue);
      fraction[index] = fractionSubset[m];
      optimum[index] = optimumSubset[m];
    }
    logger.info("\n   Total permutations evaluated: " + evaluatedPermutations + "\n");
  }

  /**
   * Return the populations of the titratable residue states and print
   *
   * @param residues residues in the system
   * @return double array of populations
   */
  public double[][] getProtonationPopulations(Residue[] residues) {
    double[][] populations = new double[residues.length][3];
    int residueIndex = 0;
    for (Residue residue : residues) {
      int index = residueList.indexOf(residue);
      // Set sums for to protonated, deprotonated, and tautomer states of titratable residues
      double protSum = 0;
      double deprotSum = 0;
      double tautomerSum = 0;
      Rotamer[] rotamers = residue.getRotamers();
      for (int rotIndex = 0; rotIndex < rotamers.length; rotIndex++) {
        switch (rotamers[rotIndex].getName()) {
          case "HIS":
          case "LYS":
          case "GLH":
          case "ASH":
          case "CYS":
            protSum += fraction[index][rotIndex];
            populations[residueIndex][0] = protSum;
            break;
          case "HIE":
          case "LYD":
          case "GLU":
          case "ASP":
          case "CYD":
            deprotSum += fraction[index][rotIndex];
            populations[residueIndex][1] = deprotSum;
            break;
          case "HID":
            tautomerSum += fraction[index][rotIndex];
            populations[residueIndex][2] = tautomerSum;
            break;
          default:
            break;
        }
      }
      String formatedProtSum = format("%.6f", protSum);
      String formatedDeprotSum = format("%.6f", deprotSum);
      String formatedTautomerSum = format("%.6f", tautomerSum);
      switch (residue.getName()) {
        case "HIS":
        case "HIE":
        case "HID":
          logger.info(residue.getResidueNumber() + "\tHIS" + "\t" + formatedProtSum + "\t" +
              "HIE" + "\t" + formatedDeprotSum + "\t" +
              "HID" + "\t" + formatedTautomerSum);
          break;
        case "LYS":
        case "LYD":
          logger.info(residue.getResidueNumber() + "\tLYS" + "\t" + formatedProtSum + "\t" +
              "LYD" + "\t" + formatedDeprotSum);
          break;
        case "ASH":
        case "ASP":
          logger.info(residue.getResidueNumber() + "\tASP" + "\t" + formatedDeprotSum + "\t" +
              "ASH" + "\t" + formatedProtSum);
          break;
        case "GLH":
        case "GLU":
          logger.info(residue.getResidueNumber() + "\tGLU" + "\t" + formatedDeprotSum + "\t" +
              "GLH" + "\t" + formatedProtSum);
          break;
        case "CYS":
        case "CYD":
          logger.info(residue.getResidueNumber() + "\tCYS" + "\t" + formatedProtSum + "\t" +
              "CYD" + "\t" + formatedDeprotSum);
          break;
        default:
          break;
      }
      residueIndex++;
    }


    return populations;
  }

  public double slidingWindowCentered(List<Residue> residueList) throws Exception {
    String[] titratableResidues = new String[]{"HIS", "HIE", "HID", "GLU", "GLH", "ASP", "ASH", "LYS", "LYD", "CYS", "CYD"};
    List<String> titratableResiudesList = Arrays.asList(titratableResidues);
    // TO DO make generic for a list of given residue centers, make array that is filled with TR or take array of centers
    double e = 0.0;
    for (Residue titrationResidue : residueList) {
      if (titratableResiudesList.contains(titrationResidue.getName())) {
        e = slidingWindowOptimization(residueList, windowSize, increment, revert, distance, Direction.FORWARD, residueList.indexOf(titrationResidue));
      }
    }
    return e;
  }

  /**
   * Return the rotamer index for each conformer (A,B,C) in xray and realspace genZ
   *
   * @return int array of rotamer indexes for each conformer (A,B,C)
   * @throws Exception
   */
  public int[][] getConformers() throws Exception {
    int[][] conformers = new int[fraction.length][3];
    for (int i = 0; i < fraction.length; i++) {
      double[] tempArray = new double[fraction[0].length];
      java.lang.System.arraycopy(fraction[i], 0, tempArray, 0, fraction[0].length);
      List<Double> elements = new ArrayList<Double>();
      for (int j = 0; j < tempArray.length; j++) {
        elements.add(tempArray[j]);
      }
      Arrays.sort(tempArray);
      int count = -1;
      for (int k = tempArray.length - 3; k < tempArray.length; k++) {
        count++;
        conformers[i][count] = elements.indexOf(tempArray[k]);
      }
    }
    return conformers;
  }


  /**
   * Return an integer array of optimized rotamers following rotamer optimization.
   *
   * @return The optimal rotamer array.
   */
  public int[] getOptimumRotamers() {
    return optimum;
  }

  /**
   * Independent optimization treats each residue sequentially.
   *
   * @param residues The list of residues to be optimized.
   * @return The final energy.
   */
  private double independent(List<Residue> residues) {
    double e = 0.0;
    List<Residue> singletonResidue = new ArrayList<>(Collections.nCopies(1, null));
    for (int i = 0; i < residues.size(); i++) {
      Residue residue = residues.get(i);
      singletonResidue.set(0, residue);
      logger.info(format(" Optimizing %s side-chain.", residue));
      Rotamer[] rotamers = residue.getRotamers();
      e = Double.MAX_VALUE;
      int bestRotamer = 0;
      double startingEnergy = 0.0;
      for (int j = 0; j < rotamers.length; j++) {
        Rotamer rotamer = rotamers[j];
        RotamerLibrary.applyRotamer(residue, rotamer);
        if (algorithmListener != null) {
          algorithmListener.algorithmUpdate(molecularAssembly);
        }
        double newE = Double.NaN;
        try {
          if (rotamer.isTitrating) {
            newE = currentEnergy(singletonResidue) + rotamer.getRotamerPhBias();
          } else {
            newE = currentEnergy(singletonResidue);
          }
          // Report energies relative to the first rotamer.
          if (j == 0) {
            startingEnergy = newE;
            newE = 0.0;
          } else {
            newE -= startingEnergy;
          }
          logger.info(format("  Energy %8s %-2d: %s", residue.toString(rotamers[j]), j, formatEnergy(newE)));
          double singularityThreshold = -100000;
          if (newE < singularityThreshold) {
            String message = format("   Rejecting as energy (%s << %s) is likely an error.", formatEnergy(newE), formatEnergy(singularityThreshold));
            logger.info(message);
            newE = Double.MAX_VALUE;
          }
        } catch (ArithmeticException ex) {
          logger.info(format(" Exception %s in energy calculations during independent for %s-%d", ex, residue, j));
        }
        if (newE < e) {
          e = newE;
          bestRotamer = j;
        }
      }
      Rotamer rotamer = rotamers[bestRotamer];
      RotamerLibrary.applyRotamer(residue, rotamer);
      optimum[i] = bestRotamer;
      logger.info(format(" Best Energy %8s %-2d: %s", residue.toString(rotamer), bestRotamer, formatEnergy(e)));

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
   * @param residues        Array of residues.
   * @param i               Current permutation.
   * @param currentRotamers Current list of rotamers.
   * @return If valid permutation found.
   */
  private boolean firstValidPerm(Residue[] residues, int i, int[] currentRotamers) {
    int nResidues = residues.length;
    Residue residuei = residues[i];
    Rotamer[] rotamersi = residuei.getRotamers();
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
      // At the end of the recursion, check each rotamer of the final residue.
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
  protected double globalOptimization(List<Residue> residueList) {
    int currentEnsemble = Integer.MAX_VALUE;
    Residue[] residues = residueList.toArray(new Residue[0]);
    int nResidues = residues.length;
    int[] currentRotamers = new int[nResidues];
    optimumSubset = new int[nResidues];

    int iterations = 0;
    boolean finalTry = false;
    int bestEnsembleTargetDiffThusFar = Integer.MAX_VALUE;
    double bestBufferThusFar = ensembleBuffer;
    double startingBuffer = ensembleBuffer;

    if (ensembleEnergy > 0.0) {
      ensembleBuffer = ensembleEnergy;
      applyEliminationCriteria(residues, true, true);
      // Compute the number of permutations without eliminating dead-ends and compute the number of
      // permutations using singleton elimination.
      double permutations = 1;
      double singletonPermutations = 1;
      for (int i = 0; i < nResidues; i++) {
        Residue residue = residues[i];
        Rotamer[] rotamers = residue.getRotamers();
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
        logger.severe(
            " No valid path through rotamer space found; try recomputing without pruning or using ensemble.");
      }
      if (rank0 && printFiles && ensembleFile == null) {
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
            logger.warning(
                format(" Versioning failed: appending to end of file %s", ensembleFile.getName()));
          }
        }
        ensembleFilter = new PDBFilter(new File(ensembleFile.getName()), molecularAssembly, null,
            null);
        logger.info(format(" Ensemble file: %s", ensembleFile.getName()));
      }
      logIfRank0(format("%30s %35s %35s", "Condition", "Number of Permutations Left",
          "Number of Permutations Removed"));
      logIfRank0(format("%30s %35s %35s", "No Eliminations", permutations, ""));
      logIfRank0(format("%30s %35s %35s", "Single Eliminations", singletonPermutations,
          permutations - singletonPermutations));
      logIfRank0(
          format("%30s %35s %35s", "Pair Eliminations", afterPairElim, pairTotalElimination));
      logIfRank0(
          format("%30s %35s %35s", "Single and Pair Eliminations", (double) evaluatedPermutations,
              pairTotalElimination + (permutations - singletonPermutations)));

      logIfRank0("\n Energy of permutations:");
      logIfRank0(format(" %12s %25s %25s", "Permutation", "Energy", "Lowest Possible Energy"));

      double e;
      if (useMonteCarlo()) {
        firstValidPerm(residues, 0, currentRotamers);
        arraycopy(currentRotamers, 0, optimumSubset, 0, nResidues);
        rotamerOptimizationMC(residues, optimumSubset, currentRotamers, nMCSteps, false, mcUseAll);

        logIfRank0(" Ensembles not currently compatible with Monte Carlo search");
        // Not currently compatible with ensembles.
      } else {
        double[] permutationEnergies = new double[evaluatedPermutations];
        ensembleStates = new ArrayList<>();

        e = rotamerOptimizationDEE(molecularAssembly, residues, 0, currentRotamers, Double.MAX_VALUE,
            optimumSubset, permutationEnergies);
        int[][] acceptedPermutations = new int[evaluatedPermutations][];
        logIfRank0(format(
            "\n Checking permutations for distance < %5.3f kcal/mol from GMEC energy %10.8f kcal/mol",
            ensembleEnergy, e));
        dryRunForEnsemble(residues, 0, currentRotamers, e, permutationEnergies,
            acceptedPermutations);
        int numAcceptedPermutations = 0;

        for (int i = 0; i < acceptedPermutations.length; i++) {
          if (acceptedPermutations[i] != null) {
            ++numAcceptedPermutations;
            logIfRank0(
                format(" Accepting permutation %d at %8.6f < %8.6f", i, permutationEnergies[i] - e,
                    ensembleEnergy));
            for (int j = 0; j < nResidues; j++) {
              Residue residuej = residues[j];
              Rotamer[] rotamersj = residuej.getRotamers();
              RotamerLibrary.applyRotamer(residuej, rotamersj[acceptedPermutations[i][j]]);
            }

            ResidueState[] states = ResidueState.storeAllCoordinates(residues);
            ensembleStates.add(new ObjectPair<>(states, permutationEnergies[i]));

            if (printFiles && rank0) {
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
                bw.write("ENDMDL");
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
        logIfRank0(format(" Number of permutations within %5.3f kcal/mol of GMEC energy: %6.4e",
            ensembleEnergy, (double) numAcceptedPermutations));
        ensembleStates.sort(null);
      }

      logIfRank0(" Final rotamers:");
      logIfRank0(format("%17s %10s %11s %12s %11s", "Residue", "Chi 1", "Chi 2", "Chi 3", "Chi 4"));
      for (int i = 0; i < nResidues; i++) {
        Residue residue = residues[i];
        Rotamer[] rotamers = residue.getRotamers();
        int ri = optimumSubset[i];
        Rotamer rotamer = rotamers[ri];
        logIfRank0(format(" %c (%7s,%2d) %s", residue.getChainID(), residue.toString(rotamer), ri,
            rotamer.toAngleString()));
        RotamerLibrary.applyRotamer(residue, rotamer);
      }
      logIfRank0("\n");

      double sumSelfEnergy = 0;
      double sumPairEnergy = 0;
      double sumTrimerEnergy = 0;
      for (int i = 0; i < nResidues; i++) {
        Residue residue = residues[i];
        Rotamer[] rotamers = residue.getRotamers();
        int ri = optimumSubset[i];
        sumSelfEnergy += eE.getSelf(i, ri);
        logIfRank0(
            format(" Final self Energy (%8s,%2d): %12.4f", residue.toString(rotamers[ri]), ri,
                eE.getSelf(i, ri)));
      }
      for (int i = 0; i < nResidues - 1; i++) {
        Residue residueI = residues[i];
        Rotamer[] rotI = residueI.getRotamers();
        int ri = optimumSubset[i];
        for (int j = i + 1; j < nResidues; j++) {
          Residue residueJ = residues[j];
          Rotamer[] rotJ = residueJ.getRotamers();
          int rj = optimumSubset[j];
          sumPairEnergy += eE.get2Body(i, ri, j, rj);
          if (eE.get2Body(i, ri, j, rj) > 10.0) {
            logIfRank0(format(" Large Final Pair Energy (%8s,%2d) (%8s,%2d): %12.4f",
                residueI.toString(rotI[ri]), ri, residueJ.toString(rotJ[rj]), rj,
                eE.get2Body(i, ri, j, rj)));
          }
        }
      }

      try {
        // Add the force field energy to the pH bias.
        e = currentEnergy(residueList) + eE.getTotalRotamerPhBias(residueList, optimumSubset, pH, pHRestraint);
      } catch (ArithmeticException ex) {
        e = Double.NaN;
        logger.severe(
            format(" Exception %s in calculating current energy at the end of triples", ex));
      }

      logIfRank0(format(" %12s %25s %25s", "Type", "Energy", "Lowest Possible Energy"));
      logIfRank0(format(" %12s %25f %25s", "Self:", sumSelfEnergy, ""));
      logIfRank0(format(" %12s %25f %25s", "Pair:", sumPairEnergy, ""));

      approximateEnergy = eE.getBackboneEnergy() + sumSelfEnergy + sumPairEnergy;

      if (threeBodyTerm) {
        for (int i = 0; i < nResidues - 2; i++) {
          int ri = optimumSubset[i];
          for (int j = i + 1; j < nResidues - 1; j++) {
            int rj = optimumSubset[j];
            for (int k = j + 1; k < nResidues; k++) {
              int rk = optimumSubset[k];
              try {
                sumTrimerEnergy += eE.get3Body(residues, i, ri, j, rj, k, rk);
              } catch (Exception ex) {
                logger.warning(ex.toString());
              }
            }
          }
        }
        approximateEnergy += sumTrimerEnergy;
        double higherOrderEnergy =
            e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - eE.getBackboneEnergy();
        logIfRank0(format(" %12s %25f %25s", "Trimer:", sumTrimerEnergy, ""));
        logIfRank0(format(" %12s %25f %25s", "Neglected:", higherOrderEnergy, ""));
      } else {
        double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - eE.getBackboneEnergy();
        logIfRank0(format(" %12s %25f %25s", "Neglected:", higherOrderEnergy, ""));
      }
      logIfRank0(format(" %12s %25f %25s", "Approximate:", approximateEnergy, ""));
      return e;
    }

    // Permutations used only to set maximum bound on ensembleNumber, thus it is safe here to put
    // that value in a 32-bit int.
    int nPerms = 1;
    for (Residue residue : residues) {
      Rotamer[] rotamers = residue.getRotamers();
      int nr = rotamers.length;
      if (nr > 1) {
        nPerms *= rotamers.length;
      }
      if (nPerms > ensembleNumber) {
        break;
      }
    }

    if (nPerms < ensembleNumber) {
      logger.info(format(
          " Requested an ensemble of %d, but only %d permutations exist; returning full ensemble",
          ensembleNumber, nPerms));
      ensembleNumber = nPerms;
    }

    while (currentEnsemble != ensembleNumber) {
      if (monteCarlo) {
        logIfRank0(" Ensemble search not currently compatible with Monte Carlo");
        ensembleNumber = 1;
      }

      if (iterations == 0) {
        applyEliminationCriteria(residues, true, true);
      } else {
        applyEliminationCriteria(residues, false, false);
      }

      // Compute the number of permutations without eliminating dead-ends and compute the number of
      // permutations using singleton elimination.
      double permutations = 1;
      double singletonPermutations = 1;
      for (int i = 0; i < nResidues; i++) {
        Residue residue = residues[i];
        Rotamer[] rotamers = residue.getRotamers();
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

      logIfRank0(" Collecting Permutations:");
      dryRun(residues, 0, currentRotamers);

      double pairTotalElimination = singletonPermutations - (double) evaluatedPermutations;
      double afterPairElim = singletonPermutations - pairTotalElimination;
      currentEnsemble = evaluatedPermutations;
      if (ensembleNumber == 1 && currentEnsemble == 0) {
        logger.severe(
            " No valid path through rotamer space found; try recomputing without pruning or using ensemble.");
      }
      if (ensembleNumber > 1) {
        if (rank0 && printFiles && ensembleFile == null) {
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
              logger.warning(
                  format(" Versioning failed: appending to end of file %s", ensembleFile.getName()));
            }
          }
          ensembleFilter = new PDBFilter(new File(ensembleFile.getName()), molecularAssembly, null,
              null);
          logger.info(format(" Ensemble file: %s", ensembleFile.getName()));
        }
        logIfRank0(format(" Ensemble Search Stats: (buffer: %5.3f, current: %d, target: %d)",
            ensembleBuffer, currentEnsemble, ensembleNumber));
      }
      if (ensembleNumber == 1 || finalTry) {
        logIfRank0(format("%30s %35s %35s", "Condition", "Number of Permutations Left",
            "Number of Permutations Removed"));
        logIfRank0(format("%30s %35s %35s", "No Eliminations", permutations, ""));
        logIfRank0(format("%30s %35s %35s", "Single Eliminations", singletonPermutations,
            permutations - singletonPermutations));
        logIfRank0(
            format("%30s %35s %35s", "Pair Eliminations", afterPairElim, pairTotalElimination));
        logIfRank0(
            format("%30s %35s %35s", "Single and Pair Eliminations", (double) evaluatedPermutations,
                pairTotalElimination + (permutations - singletonPermutations)));

        logIfRank0("\n Energy of permutations:");
        logIfRank0(format(" %12s %25s %25s", "Permutation", "Energy", "Lowest Possible Energy"));

        break;
      }
      if (abs(currentEnsemble - ensembleNumber) < bestEnsembleTargetDiffThusFar) {
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
          // TODO: Decide whether we like these next four lines.  Has the potential to produce a
          // crazy amount of permutations.
          logIfRank0(" Ensemble still empty; increasing buffer energy.");
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
      logger.warning(
          " No valid rotamer permutations found; results will be unreliable.  Try increasing the starting ensemble buffer.");
    }
    double[] permutationEnergyStub = null;

    if (useMonteCarlo()) {
      firstValidPerm(residues, 0, currentRotamers);
      rotamerOptimizationMC(residues, optimumSubset, currentRotamers, nMCSteps, false, mcUseAll);
    } else {
      rotamerOptimizationDEE(molecularAssembly, residues, 0, currentRotamers, Double.MAX_VALUE,
          optimumSubset, permutationEnergyStub);
    }

    double[] residueEnergy = new double[nResidues];

    double sumSelfEnergy = 0;
    double sumLowSelfEnergy = 0;

    logIfRank0("\n Energy contributions:");
    logIfRank0(format(" %15s %25s %25s", "Type", "Energy", "Lowest Possible Energy"));

    for (int i = 0; i < nResidues; i++) {
      int ri = optimumSubset[i];
      Residue residue = residues[i];
      Rotamer[] rotamers = residue.getRotamers();
      turnOnAtoms(residue);
      RotamerLibrary.applyRotamer(residue, rotamers[ri]);
      double self = eE.getSelf(i, ri);
      residueEnergy[i] = self;
      sumSelfEnergy += self;
      double lowest = eE.lowestSelfEnergy(residues, i);
      sumLowSelfEnergy += lowest;
      if (self - lowest > 10.0) {
        logIfRank0(
            format(" %15s %25f %25f", "Self (" + residues[i] + "," + ri + "):", self, lowest));
      }
    }

    double sumPairEnergy = 0.0;
    double sumLowPairEnergy = 0.0;
    double[] resPairEnergy = new double[nResidues];
    double[] lowPairEnergy = new double[nResidues];
    for (int i = 0; i < nResidues - 1; i++) {
      StringBuilder sb = new StringBuilder();
      int ri = optimumSubset[i];
      double sumPairEnergyI = 0;
      double sumLowPairEnergyI = 0;
      for (int j = i + 1; j < nResidues; j++) {
        int rj = optimumSubset[j];
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
        logIfRank0(
            format(" %15s %25f %25f", "Self (" + residues[i] + "," + ri + ")", sumPairEnergyI,
                sumLowPairEnergyI));
        sb.trimToSize();
        if (!sb.toString().isEmpty()) {
          logIfRank0(sb.toString());
        }
      }
    }

    double e = Double.NaN;
    try {
      // Add the force field energy to the pH bias.
      e = currentEnergy(residueList) + eE.getTotalRotamerPhBias(residueList, optimumSubset, pH, pHRestraint);
    } catch (ArithmeticException ex) {
      logger.severe(
          format(" Exception %s in calculating current energy at the end of self and pairs", ex));
    }
    logIfRank0(format(" %15s %25f %25s", "Backbone:", eE.getBackboneEnergy(), ""));
    logIfRank0(format(" %15s %25f %25f", "Self:", sumSelfEnergy, sumLowSelfEnergy));
    logIfRank0(format(" %15s %25f %25f", "Pair:", sumPairEnergy, sumLowPairEnergy));

    approximateEnergy = eE.getBackboneEnergy() + sumSelfEnergy + sumPairEnergy;

    double sumTrimerEnergy = 0;
    if (threeBodyTerm) {
      for (int i = 0; i < nResidues - 2; i++) {
        int ri = optimumSubset[i];
        for (int j = i + 1; j < nResidues - 1; j++) {
          int rj = optimumSubset[j];
          for (int k = j + 1; k < nResidues; k++) {
            int rk = optimumSubset[k];
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
      double higherOrderEnergy =
          e - sumSelfEnergy - sumPairEnergy - sumTrimerEnergy - eE.getBackboneEnergy();
      logIfRank0(format(" %15s %25f %25s", "Trimer:", sumTrimerEnergy, ""));
      logIfRank0(format(" %15s %25f %25s", "Neglected:", higherOrderEnergy, ""));
    } else {
      double higherOrderEnergy = e - sumSelfEnergy - sumPairEnergy - eE.getBackboneEnergy();
      logIfRank0(format(" %15s %25f %25s", "Neglected:", higherOrderEnergy, ""));
    }

    logIfRank0(format(" %15s %25f %25s", "Approximate:", approximateEnergy, ""));

    logIfRank0("\n Final rotamers:");
    logIfRank0(
        format("%17s %10s %11s %12s %11s %14s", "Residue", "Chi 1", "Chi 2", "Chi 3", "Chi 4",
            "Energy"));
    for (int i = 0; i < nResidues; i++) {
      Residue residue = residues[i];
      Rotamer[] rotamers = residue.getRotamers();
      int ri = optimumSubset[i];
      Rotamer rotamer = rotamers[ri];
      logIfRank0(format(" %3d %c (%7s,%2d) %s %12.4f ", i + 1, residue.getChainID(),
          residue.toString(rotamers[ri]), ri, rotamer.toAngleString(), residueEnergy[i]));
      RotamerLibrary.applyRotamer(residue, rotamer);
    }
    logIfRank0("\n");
    return e;
  }

  /**
   * Use Monte Carlo if monteCarlo specified, and either skipDEE specified or nMCsteps is smaller
   * than the remaining permutation size.
   *
   * @return Finish DEE search with Monte Carlo.
   */
  private boolean useMonteCarlo() {
    return monteCarlo && (mcNoEnum || (nMCSteps < evaluatedPermutations));
  }

  /**
   * Performs a recursive brute-force rotamer optimization over a passed list of residues.
   *
   * @param residueList Residues to be optimized.
   * @return Global minimum energy conformation energy.
   */
  private double bruteForce(List<Residue> residueList) {

    Residue[] residues = residueList.toArray(new Residue[0]);
    int nResidues = residues.length;

    // Compute the number of permutations without eliminating dead-ends.
    double permutations = 1;
    for (Residue residue : residues) {
      Rotamer[] rotamers = residue.getRotamers();
      permutations *= rotamers.length;
    }

    logger.info(format(" Number of permutations: %16.8e.", permutations));

    double e;
    useForceFieldEnergy = molecularAssembly.getProperties()
        .getBoolean("ro-use-force-field-energy", false);
    if (!useForceFieldEnergy) {
      // Use a many-body energy sum to evaluate permutations.
      setPruning(0);
      rotamerEnergies(residues);
      int[] currentRotamers = new int[nResidues];
      e = decomposedRotamerOptimization(residues, 0, Double.MAX_VALUE, optimum, currentRotamers);

      // Set each residue to their optimal rotamer.
      for (int i = 0; i < nResidues; i++) {
        Residue residue = residues[i];
        Rotamer[] rotamers = residue.getRotamers();
        RotamerLibrary.applyRotamer(residue, rotamers[optimum[i]]);
        turnOnAtoms(residue);
      }

      // Compute the full energy for comparison (i.e., without many-body truncation).
      double fullEnergy = 0;
      try {
        // Add the force field energy to the pH bias.
        fullEnergy = currentEnergy(residueList) + eE.getTotalRotamerPhBias(residueList, optimum, pH, pHRestraint);
      } catch (Exception ex) {
        logger.severe(format(" Exception %s in calculating full energy; FFX shutting down", ex));
      }

      logger.info(format(" Final summation of energies:    %16.5f", e));
      logger.info(format(" Final energy of optimized structure:    %16.5f", fullEnergy));
      logger.info(format(" Neglected:    %16.5f", fullEnergy - e));
    } else {
      e = rotamerOptimization(molecularAssembly, residues, 0, Double.MAX_VALUE, optimum);
    }

    for (int i = 0; i < nResidues; i++) {
      Residue residue = residues[i];
      Rotamer[] rotamers = residue.getRotamers();
      int ri = optimum[i];
      Rotamer rotamer = rotamers[ri];
      logger.info(format(" %s %s (%d)", residue.getResidueNumber(), rotamer.toString(), ri));
      RotamerLibrary.applyRotamer(residue, rotamer);
      if (useForceFieldEnergy) {
        try {
          e = currentEnergy(residueList) + eE.getTotalRotamerPhBias(residueList, optimum, pH, pHRestraint);
        } catch (ArithmeticException ex) {
          logger.fine(
              format(" Exception %s in calculating full AMOEBA energy at the end of brute force",
                  ex));
        }
      }
    }
    return e;
  }

  @SuppressWarnings("fallthrough")
  private double slidingWindowOptimization(List<Residue> residueList, int windowSize, int increment,
                                           boolean revert, double distance, Direction direction, int windowCenter) throws Exception {

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
        List<Residue> temp = new ArrayList<>();
        for (int i = nOptimize - 1; i >= 0; i--) {
          temp.add(residueList.get(i));
        }
        residueList = temp;
        // Fall through into the FORWARD case.
      case FORWARD:
        for (int windowStart = 0; windowStart + (windowSize - 1) < nOptimize;
             windowStart += increment) {
          long windowTime = -System.nanoTime();
// Set the start at the residue based on my center
          if (windowCenter > -1) {
            windowStart = windowCenter - windowSize / 2;
            if (windowStart < 0) {
              windowStart = 0;
            }
            windowEnd = windowCenter + windowSize / 2;
            if (windowEnd >= residueList.size()) {
              windowEnd = residueList.size() - 1;
            }
          } else {
            windowEnd = windowStart + (windowSize - 1);
          }

          if (windowCenter > -1) {
            logIfRank0(format("\n Center window at residue %d.\n", residueList.get(windowCenter).getResidueNumber()));
          } else {
            logIfRank0(format("\n Iteration %d of the sliding window.\n", counter++));
          }

          Residue firstResidue = residueList.get(windowStart);
          Residue lastResidue = residueList.get(windowEnd);
          if (firstResidue != lastResidue) {
            logIfRank0(
                format(" Residues %s ... %s", firstResidue.toString(), lastResidue.toString()));
          } else {
            logIfRank0(format(" Residue %s", firstResidue.toString()));
          }
          List<Residue> currentWindow = new ArrayList<>();
          for (int i = windowStart; i <= windowEnd; i++) {
            Residue residue = residueList.get(i);
            currentWindow.add(residue);
          }

          if (distance > 0) {
            for (int i = windowStart; i <= windowEnd; i++) {
              if (windowCenter > -1) {
                i = windowCenter;
              }
              Residue residuei = residueList.get(i);
              int indexI = allResiduesList.indexOf(residuei);
              int lengthRi;
              lengthRi = residuei.getRotamers().length;
              for (int ri = 0; ri < lengthRi; ri++) {
                for (int j = 0; j < nAllResidues; j++) {
                  Residue residuej = allResiduesArray[j];
                  Rotamer[] rotamersj = residuej.getRotamers();
                  if (currentWindow.contains(residuej) || rotamersj == null) {
                    continue;
                  }
                  int lengthRj = rotamersj.length;
                  for (int rj = 0; rj < lengthRj; rj++) {
                    double rotamerSeparation = dM.get2BodyDistance(indexI, ri, j, rj);
                    if (rotamerSeparation <= distance) {
                      if (!currentWindow.contains(residuej)) {
                        logIfRank0(format(" Adding residue %s at distance %16.8f Ang from %s %d.",
                            residuej.toFormattedString(false, true), rotamerSeparation,
                            residuei.toFormattedString(false, true), ri));
                        currentWindow.add(residuej);
                      }
                      break;
                    }
                  }
                }
              }
              if (windowCenter > -1) {
                break;
              }
            }
          }

          /*
           If the window starts with a nucleic acid, and there is a 5' NA residue, ensure that 5'
           NA residue has been included in the window. Otherwise, that previous residue may not
           have had a chance to be flexible about its sugar pucker.
           <p>
           If window size is greater than increment, however, this has already been handled.
           Additionally, do not perform this for the first window (counter is already incremented
           by the time this check is performed, so first window's counter will be 2). Furthermore,
           do not include Residues with null Rotamer lists (this breaks things).
           <p>
           The issue: if window size = increment, the last NA residue in each window will not
           have flexibility about its sugar pucker, because its self-energy includes the O3' (i)
           to P (i+1) bond, so it must remain in the original sugar pucker to meet the i+1
           residue. However, this problem can be solved by ensuring that final residue is included
           in the next window, where it will have flexibility about its sugar pucker.
           <p>
           If you are running successive sliding window jobs on the same file, I would suggest
           starting the next job on the last residue of the previous job, unless you know your
           settings will include it.
          */
          if (counter > 2 && windowSize <= increment && firstResidue.getResidueType() == NA) {
            Residue prevResidue = firstResidue.getPreviousResidue();
            if (prevResidue != null && prevResidue.getResidueType() == NA
                && !currentWindow.contains(prevResidue) && prevResidue.getRotamers() != null) {
              logIfRank0(format(
                  " Adding nucleic acid residue 5' of window start %s to give it flexibility about its sugar pucker.",
                  prevResidue));
              currentWindow.add(prevResidue);
            }
          }
          sortResidues(currentWindow);

          if (revert) {
            ResidueState[] coordinates = ResidueState.storeAllCoordinates(currentWindow);
            double startingEnergy = Double.NaN;
            try {
              startingEnergy = currentEnergy(currentWindow);
            } catch (ArithmeticException ex) {
              logger.severe(format(
                  " Exception %s in calculating starting energy of a window; FFX shutting down",
                  ex));
            }
            globalOptimization(currentWindow);
            double finalEnergy = Double.NaN;
            try {
              finalEnergy = currentEnergy(currentWindow);
            } catch (ArithmeticException ex) {
              logger.severe(
                  format(" Exception %s in calculating final energy of a window; FFX shutting down",
                      ex));
            }
            if (startingEnergy <= finalEnergy) {
              logger.info(
                  "Optimization did not yield a better energy. Reverting to original coordinates.");
              ResidueState.revertAllCoordinates(currentWindow, coordinates);
            } else {
              // Copy sliding window optimal rotamers into the overall optimum array.
              int i = 0;
              for (Residue residue : currentWindow) {
                int index = residueList.indexOf(residue);
                optimum[index] = optimumSubset[i++];
              }
            }
          } else {
            globalOptimization(currentWindow);
            // Copy sliding window optimal rotamers into the overall optimum array.
            int i = 0;
            for (Residue residue : currentWindow) {
              int index = residueList.indexOf(residue);
              optimum[index] = optimumSubset[i++];
            }
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

          if (rank0 && printFiles) {
            File file = molecularAssembly.getFile();
            if (firstWindowSaved) {
              file.delete();
            }
            // Don't write a file if it's the final iteration.
            if (windowStart + windowSize == nOptimize) {
              continue;
            }
            PDBFilter windowFilter = new PDBFilter(file, molecularAssembly, null, null);
            //   StringBuilder header = new StringBuilder(format("Iteration %d of the sliding
            // window\n", counter - 1));
            try {
              windowFilter.writeFile(file, false);
              if (firstResidue != lastResidue) {
                logger.info(
                    format(" File with residues %s ... %s in window written to.", firstResidue,
                        lastResidue));
              } else {
                logger.info(format(" File with residue %s in window written to.", firstResidue));
              }
            } catch (Exception e) {
              logger.warning(format("Exception writing to file: %s", file.getName()));
            }
            firstWindowSaved = true;
          }
          long currentTime = System.nanoTime();
          windowTime += currentTime;
          logIfRank0(format(" Time elapsed for this iteration: %11.3f sec", windowTime * 1.0E-9));
          logIfRank0(
              format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
          if (genZ) {
            int[] currentRotamers = new int[optimumSubset.length];
            usingBoxOptimization = true;
            Residue[] residueSubsetArray = currentWindow.toArray(new Residue[currentWindow.size()]);
            getFractions(residueSubsetArray, 0, currentRotamers, true);
            Residue[] titrationResidueArray = new Residue[]{residueList.get(windowCenter)};
            getProtonationPopulations(titrationResidueArray);
          }
          if (windowCenter > -1) {
            break;
          }
        }


        break;
      default:
        // No default case.
        break;
    }
    return 0.0;
  }

  /**
   * Sorts a passed List of Residues by global index.
   *
   * @param residues List of Residues to be sorted.
   */
  private void sortResidues(List<Residue> residues) {
    Comparator<Residue> comparator = Comparator.comparing(Residue::getChainID)
        .thenComparingInt(Residue::getResidueNumber);
    residues.sort(comparator);
  }

  /**
   * applyEliminationCriteria.
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
      int i = 0;
      boolean pairEliminated;
      do {
        pairEliminated = false;
        if (useGoldstein) {
          if (selfEliminationOn) {
            i++;
            logIfRank0(format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
            // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
            while (goldsteinDriver(residues)) {
              i++;
              logIfRank0(this.toString());
              logIfRank0(
                  format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
            }
          }
          if (pairEliminationOn) {
            i++;
            logIfRank0(
                format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
            // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
            pairEliminated = false;
            while (goldsteinPairDriver(residues)) {
              pairEliminated = true;
              i++;
              logIfRank0(this.toString());
              logIfRank0(
                  format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
            }
          }
        } else {
          if (selfEliminationOn) {
            i++;
            logIfRank0(format("\n Iteration %d: Applying Single DEE conditions ", i));
            // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
            while (deeRotamerElimination(residues)) {
              i++;
              logIfRank0(toString());
              logIfRank0(format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
            }
          }
          if (pairEliminationOn) {
            i++;
            logIfRank0(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
            // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
            pairEliminated = false;
            while (deeRotamerPairElimination(residues)) {
              pairEliminated = true;
              i++;
              logIfRank0(toString());
              logIfRank0(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
            }
          }
        }
        eR.validateDEE(residues);
        logIfRank0(toString());
      } while (pairEliminated);
      logIfRank0(" Self-consistent DEE rotamer elimination achieved.\n");
    }
    if (!verbose) {
      logger.setLevel(prevLevel);
    }
  }

  private void applyEliminationCriteria(Residue[] residues) {
    if (verboseEnergies) {
      try {
        logIfRank0(format("\n Beginning Energy %s", formatEnergy(currentEnergy(residues))));
      } catch (ArithmeticException ex) {
        logger.severe(format(" Exception %s in calculating beginning energy; FFX shutting down.", ex));
      }
    }

    rotamerEnergies(residues);

    if (testing) {
      int nres = residues.length;
      eR.onlyPrunedSingles = new boolean[nres][];
      eR.onlyPrunedPairs = new boolean[nres][][][];
      for (int i = 0; i < nres; i++) {
        Residue residuei = residues[i];
        Rotamer[] rotamersi = residuei.getRotamers();
        int lenri = rotamersi.length; // Length rotamers i
        eR.onlyPrunedSingles[i] = new boolean[lenri];
        eR.onlyPrunedSingles[i] = copyOf(eR.eliminatedSingles[i], eR.eliminatedSingles[i].length);
        eR.onlyPrunedPairs[i] = new boolean[lenri][][];
        // Loop over the set of rotamers for residue i.
        for (int ri = 0; ri < lenri; ri++) {
          eR.onlyPrunedPairs[i][ri] = new boolean[nres][];
          for (int j = i + 1; j < nres; j++) {
            Residue residuej = residues[j];
            Rotamer[] rotamersj = residuej.getRotamers();
            int lenrj = rotamersj.length;
            eR.onlyPrunedPairs[i][ri][j] = new boolean[lenrj];
            eR.onlyPrunedPairs[i][ri][j] = copyOf(eR.eliminatedPairs[i][ri][j],
                eR.eliminatedPairs[i][ri][j].length);
          }
        }
      }
    }

    if (testSelfEnergyEliminations) {
      testSelfEnergyElimination(residues);
    } else if (testPairEnergyEliminations > -1) {
      testPairEnergyElimination(residues, testPairEnergyEliminations);
    } else if (testTripleEnergyEliminations1 > -1 && testTripleEnergyEliminations2 > -1) {
      testTripleEnergyElimination(residues, testTripleEnergyEliminations1,
          testTripleEnergyEliminations2);
    }

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
          logIfRank0(format("\n Iteration %d: Applying Single Goldstein DEE conditions ", i));
          // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
          while (goldsteinDriver(residues)) {
            i++;
            logIfRank0(this.toString());
            logIfRank0(
                format("\n Iteration %d: Applying Single Rotamer Goldstein DEE conditions ", i));
          }
        }
        if (pairEliminationOn) {
          i++;
          logIfRank0(format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
          // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
          while (goldsteinPairDriver(residues)) {
            pairEliminated = true;
            i++;
            logIfRank0(this.toString());
            logIfRank0(
                format("\n Iteration %d: Applying Rotamer Pair Goldstein DEE conditions ", i));
          }
        }
      } else {
        if (selfEliminationOn) {
          i++;
          logIfRank0(format("\n Iteration %d: Applying Single DEE conditions ", i));
          // While there are eliminated rotamers, repeatedly apply single rotamer elimination.
          while (deeRotamerElimination(residues)) {
            i++;
            logIfRank0(toString());
            logIfRank0(format("\n Iteration %d: Applying Single Rotamer DEE conditions ", i));
          }
        }
        if (pairEliminationOn) {
          i++;
          logIfRank0(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
          // While there are eliminated rotamer pairs, repeatedly apply rotamer pair elimination.
          while (deeRotamerPairElimination(residues)) {
            pairEliminated = true;
            i++;
            logIfRank0(toString());
            logIfRank0(format("\n Iteration %d: Applying Rotamer Pair DEE conditions ", i));
          }
        }
      }
      eR.validateDEE(residues);
      logIfRank0(toString());
    } while (pairEliminated);

    logIfRank0(" Self-consistent DEE rotamer elimination achieved.\n");
  }

  /**
   * Calculates the energy at the current state.
   *
   * @param resList List of residues in current energy term.
   * @return Energy of the current state.
   */
  private double currentEnergy(List<Residue> resList) throws ArithmeticException {
    List<Rotamer> rots = resList.stream().filter(Objects::nonNull).map(Residue::getRotamer)
        .collect(Collectors.toList());
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
      Rotamer[] rotsi = resi.getRotamers();
      int lenri = rotsi.length;
      int indexI = allResiduesList.indexOf(resi);

      for (int j = 0; j < nRes; j++) {
        if (i == j) {
          continue;
        }
        Residue resj = residues[j];
        Rotamer[] rotsj = resj.getRotamers();
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
      nI = nearby.stream().mapToInt(Integer::intValue).filter(j -> j > fi).toArray();
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
    logger.info(format("\n Energy restart file: %s", restartFile.getName()));
  }

  /**
   * Turn off non-bonded contributions from all residues except for one. Compute the self-energy for
   * each residue relative to the backbone contribution.
   *
   * @param residues A list of residues that will undergo rotamer optimization.
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

    eR = new EliminatedRotamers(this, dM, allResiduesList, maxRotCheckDepth, clashThreshold,
        pairClashThreshold, multiResClashThreshold, nucleicPruningFactor, nucleicPairsPruningFactor,
        multiResPairClashAddn, pruneClashes, prunePairClashes, print, residues);

    if (decomposeOriginal) {
      assert library.getUsingOrigCoordsRotamer();
      for (int i = 0; i < nResidues; i++) {
        Residue resi = residues[i];
        Rotamer[] rotsi = resi.getRotamers();
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

    eE = new EnergyExpansion(this, dM, eR, molecularAssembly, potential, algorithmListener,
        allResiduesList, resNeighbors, threeBodyTerm, decomposeOriginal, usingBoxOptimization,
        verbose, pruneClashes, prunePairClashes, rank0);

    // Update the EliminatedRotamers instance with the EnergyExpansion instance.
    eR.setEnergyExpansion(eE);

    int loaded = 0;
    if (loadEnergyRestart) {
      if (usingBoxOptimization) {
        loaded = eE.loadEnergyRestart(energyRestartFile, residues, boxOpt.boxLoadIndex, boxOpt.boxLoadCellIndices);
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

      SelfEnergyRegion selfEnergyRegion = new SelfEnergyRegion(this, eE, eR, residues, energyWriter,
          world, numProc, pruneClashes, rank0, rank, verbose, writeEnergyRestart, printFiles);
      energyWorkerTeam.execute(selfEnergyRegion);
      long singlesTime = System.nanoTime() - energyStartTime;
      logIfRank0(format(" Time for single energies: %12.4g", (singlesTime * 1.0E-9)));
      if (logger.isLoggable(Level.FINE)) {
        Resources.logResources();
      }

      if (loaded < 2) {
        eE.allocate2BodyJobMap(residues, nResidues, false);
      }

      TwoBodyEnergyRegion twoBodyEnergyRegion = new TwoBodyEnergyRegion(this, dM, eE, eR, residues,
          allResiduesList, energyWriter, world, numProc, prunePairClashes, superpositionThreshold,
          rank0, rank, verbose, writeEnergyRestart, printFiles);
      energyWorkerTeam.execute(twoBodyEnergyRegion);
      long pairsTime = System.nanoTime() - (singlesTime + energyStartTime);

      long triplesTime = 0;
      long quadsTime = 0;
      logIfRank0(format(" Time for 2-body energies:   %12.4g", (pairsTime * 1.0E-9)));
      if (logger.isLoggable(Level.FINE)) {
        Resources.logResources();
      }

      if (threeBodyTerm) {
        if (loaded < 3) {
          eE.allocate3BodyJobMap(residues, nResidues, false);
        }
        ThreeBodyEnergyRegion threeBodyEnergyRegion = new ThreeBodyEnergyRegion(this, dM, eE, eR,
            residues, allResiduesList, energyWriter, world, numProc, superpositionThreshold, rank0,
            rank, verbose, writeEnergyRestart, printFiles);
        energyWorkerTeam.execute(threeBodyEnergyRegion);
        triplesTime = System.nanoTime() - (pairsTime + singlesTime + energyStartTime);
        logIfRank0(format(" Time for 3-Body energies: %12.4g", (triplesTime * 1.0E-9)));
        if (logger.isLoggable(Level.FINE)) {
          Resources.logResources();
        }
      }

      if (compute4BodyEnergy) {
        eE.allocate4BodyJobMap(residues, nResidues);

        FourBodyEnergyRegion fourBodyEnergyRegion = new FourBodyEnergyRegion(this, dM, eE, eR,
            residues, allResiduesList, superpositionThreshold);
        energyWorkerTeam.execute(fourBodyEnergyRegion);
        quadsTime = System.nanoTime() - (triplesTime + pairsTime + singlesTime + energyStartTime);
        logIfRank0(format(" Time for 4-Body energies:   %12.4g", quadsTime * 1.0E-9));
      }

      long allTime = singlesTime + pairsTime + triplesTime + quadsTime;
      logIfRank0(format(" Time for all energies:    %12.4g", allTime * 1.0E-9));
    } catch (Exception ex) {
      String message = " Exception computing rotamer energies in parallel.";
      logger.log(Level.SEVERE, message, ex);
    }

    // Turn on all atoms.
    for (Atom atom : atoms) {
      atom.setUse(true);
    }
    // Print the energy with all rotamers in their default conformation.
    if (verboseEnergies && rank0) {
      try {
        double defaultEnergy = currentEnergy(residues);
        logger.info(format(" Energy of the system with rotamers in their default conformation: %s",
            formatEnergy(defaultEnergy)));
      } catch (ArithmeticException ex) {
        logger.severe(format(" Exception %s in calculating default energy; FFX shutting down", ex));
      }
    }
    return eE.getBackboneEnergy();
  }

  /**
   * Applies the "default" rotamer: currently the 0'th rotamer.
   *
   * @param residue Residue to apply a default rotamer for.
   */
  private void applyDefaultRotamer(Residue residue) {
    RotamerLibrary.applyRotamer(residue, residue.getRotamers()[0]);
  }

  /**
   * Elimination of rotamers via the original Dead End Elimination algorithm.
   *
   * @param residues Array of residues under consideration.
   * @return True if any rotamers were eliminated.
   */
  private boolean deeRotamerElimination(Residue[] residues) {
    int nres = residues.length;
    // A flag to indicate if more rotamers or rotamer pairs were eliminated.
    boolean eliminated = false;
    // Loop over residues.
    double[] minMax = new double[2];
    double[] minEnergySingles = null;
    double[] maxEnergySingles = null;
    for (int i = 0; i < nres; i++) {
      Residue residuei = residues[i];
      Rotamer[] rotamersi = residuei.getRotamers();
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
              // There is a chance that ri conflicts with every possible rotamer of some
              // residue j.
              // In that event, its minimum energy is set NaN and should be easy to eliminate.
              minEnergySingles[ri] = Double.NaN;
            }
            if (Double.isFinite(minMax[0]) && Double.isFinite(maxEnergySingles[ri])) {
              maxEnergySingles[ri] += minMax[1];
            } else {
              // In this branch, ri conflicts with some j,rj and cannot be practically used for
              // elimination.
              maxEnergySingles[ri] = Double.NaN;
            }
          } else {
            Residue residuej = residues[j];
            logger.info(
                format(" Inconsistent Pair: %8s %2d, %8s.", residuei.toFormattedString(false, true),
                    ri, residuej.toFormattedString(false, true)));
            // eliminateRotamer(residues, i, ri, print);
          }
        }
      }

      /*
       Apply the singles elimination criteria to rotamers of residue i by determining the most
       favorable maximum energy.
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
        // This branch is taken if every ri conflicts with at least one j,rj. In that case, nothing
        // can be eliminated yet!
        logIfRank0(" Could not eliminate any i,ri because eliminationEnergy was never set!",
            Level.FINE);
      } else {
        /*
         Eliminate rotamers whose minimum energy is greater than the worst case for another
         rotamer.
        */
        for (int ri = 0; ri < lenri; ri++) {
          if (eR.check(i, ri)) {
            continue;
          }
          // If i,ri has a clash with all phase space, it can be eliminated by something that
          // doesn't clash with all phase space.
          if (!Double.isFinite(minEnergySingles[ri])) {
            if (eR.eliminateRotamer(residues, i, ri, print)) {
              logIfRank0(format("  Rotamer elimination of (%8s,%2d) that always clashes.",
                  residuei.toFormattedString(false, true), ri));
              eliminated = true;
            }
          }
          // Otherwise, can eliminate if its best possible energy is still worse than something
          // else's worst possible energy.
          if (minEnergySingles[ri] > eliminationEnergy + ensembleBuffer) {
            if (eR.eliminateRotamer(residues, i, ri, print)) {
              logIfRank0(format("  Rotamer elimination of (%8s,%2d) by (%8s,%2d): %12.4f > %6.4f.",
                  residuei.toFormattedString(false, true), ri,
                  residuei.toFormattedString(false, true), eliminatingRotamer, minEnergySingles[ri],
                  eliminationEnergy + ensembleBuffer));
              eliminated = true;
            }
          }
        }
      }
    }
    return eliminated;
  }

  /**
   * Rotamer pair elimination driver for many-body Dead End Elimination. Generally less effective
   * than Goldstein.
   *
   * @param residues Residues under consideration.
   * @return If at least one pair eliminated.
   */
  private boolean deeRotamerPairElimination(Residue[] residues) {
    int nres = residues.length;
    boolean eliminated = false;

    for (int i = 0; i < (nres - 1); i++) {
      Residue residuei = residues[i];
      Rotamer[] rotamersi = residuei.getRotamers();
      int lenri = rotamersi.length;

      // Minimum and maximum summation found for ri-rj pairs.
      double[][] minPairEnergies = new double[lenri][];
      double[][] maxPairEnergies = new double[lenri][];

      for (int j = i + 1; j < nres; j++) {
        Residue residuej = residues[j];
        Rotamer[] rotamersj = residuej.getRotamers();
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
            minPairEnergies[ri][rj] =
                eE.getSelf(i, ri) + eE.getSelf(j, rj) + eE.get2Body(i, ri, j, rj);
            maxPairEnergies[ri][rj] = minPairEnergies[ri][rj];

            // Min and max external summations for ri-rj.
            double[] minMax = new double[2];

            // Add contributions from third residues k, and possibly fourth residues l.
            if (eE.minMaxE2(residues, minMax, i, ri, j, rj)) {
              if (Double.isFinite(minPairEnergies[ri][rj]) && Double.isFinite(minMax[0])) {
                minPairEnergies[ri][rj] += minMax[0];
              } else {
                logger.severe(
                    format(" An ri-rj pair %s-%d %s-%d with NaN minimum was caught incorrectly!",
                        residuei.toFormattedString(false, true), ri,
                        residuej.toFormattedString(false, true), rj));
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
              logger.info(format(" Eliminating pair %s-%d %s-%d that always clashes.",
                  residuei.toFormattedString(false, true), ri,
                  residuej.toFormattedString(false, true), rj));
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
            if (Double.isFinite(maxPairEnergies[ri][rj])
                && maxPairEnergies[ri][rj] < pairEliminationEnergy) {
              pairEliminationEnergy = maxPairEnergies[ri][rj];
            }
          }
        }

        if (pairEliminationEnergy == Double.MAX_VALUE) {
          logIfRank0(format(
              " All rotamer pairs for residues %s and %s have possible conflicts; cannot perform any eliminations!",
              residuei.toFormattedString(false, true), residuej), Level.FINE);
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
                  logIfRank0(format(" Eliminating rotamer pair: %s %d, %s %d (%s > %s + %6.6f)",
                      residuei.toFormattedString(false, true), ri,
                      residuej.toFormattedString(false, true), rj,
                      formatEnergy(minPairEnergies[ri][rj]), formatEnergy(pairEliminationEnergy),
                      ensembleBuffer), Level.INFO);
                } else {
                  // See above check(i, ri, j, rj) for why this should not be taken!
                  logIfRank0(
                      format(" Already eliminated rotamer pair! %s %d, %s %d (%s > %1s + %6.6f)",
                          residuei.toFormattedString(false, true), ri,
                          residuej.toFormattedString(false, true), rj,
                          formatEnergy(minPairEnergies[ri][rj]), formatEnergy(pairEliminationEnergy),
                          ensembleBuffer), Level.WARNING);
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
      Rotamer[] roti = resi.getRotamers();
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
      logIfRank0(" No more single rotamers to eliminate.");
    }
    return eliminated;
  }

  /**
   * Attemps to eliminate rotamer riA based on riB.
   *
   * @param residues Array of Residues.
   * @param i        Residue i index.
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
      Rotamer[] rotj = resj.getRotamers();
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
          /*
           This is a part of configuration space where riA is valid but not riB. Thus, if j,rj is
           part of the GMEC, riB is inconsistent with it. Thus, riB cannot be used to eliminate
           riA.
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
          IntStream kStream = IntStream.concat(Arrays.stream(bidiResNeighbors[i]),
              Arrays.stream(bidiResNeighbors[j]));
          int[] possKs = kStream.distinct().sorted().toArray();
          for (int k : possKs) {
            if (k == i || k == j) {
              continue;
            }
            Residue resk = residues[k];
            Rotamer[] rotk = resk.getRotamers();
            int nrk = rotk.length;
            int rkEvals = 0;
            double minForResK = Double.MAX_VALUE;
            for (int rk = 0; rk < nrk; rk++) {
              /*
               If k,rk or j,rj-k,rk are not a part of valid configuration space, continue. If
               i,riA-k,rk or i,riA-j,rj-k,rk are not valid for riA, continue.
              */
              if (eR.check(k, rk) || eR.check(j, rj, k, rk) || eR.check(i, riA, k, rk)) {
                // Not yet implemented: check(i, riA, j, rj, k, rk) because no triples get
                // eliminated.
                continue;
              }
              /*
               If i,riB-k,rk or i,riB-j,rj-k,rk are invalid for riB, there is some part of
               configuration space for which riA is valid but not riB.
              */
              if (eR.check(i, riB, k, rk)) {
                // Not yet implemented: check(i, riB, j, rj, k, rk).
                return false;
              }

              rkEvals++;
              double e =
                  eE.get3Body(residues, i, riA, j, rj, k, rk) - eE.get3Body(residues, i, riB, j, rj,
                      k, rk);
              if (e < minForResK) {
                minForResK = e;
              }
            }
            /* If there were no 3-body interactions with residue k, then minForResk is zero. */
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
        logIfRank0(format("  Rotamer elimination of (%8s,%2d) by (%8s,%2d): %12.4f > %6.4f.",
            resi.toFormattedString(false, true), riA, resi.toFormattedString(false, true), riB,
            goldsteinEnergy, ensembleBuffer));
        logIfRank0(format("   Self: %12.4f, Pairs: %12.4f, Triples: %12.4f.", selfDiff, sumPairDiff,
            sumTripleDiff));
        return true;
      }
    }
    return false;
  }

  /**
   * Finds and eliminates rotamer pairs according to the many-body Goldstein pairs criterion.
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
      Rotamer[] rotsi = resi.getRotamers();
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
          Rotamer[] rotsj = resj.getRotamers();
          int lenrj = rotsj.length;
          for (int rjC = 0; rjC < lenrj; rjC++) {
            // Again, no point in eliminating the already-eliminated.
            if (eR.check(j, rjC) || eR.check(i, riA, j, rjC)) {
              continue;
            }
            boolean breakOut = false;

            // Now, generate pairs riB-rjD. If any pair riB-rjD eliminates riA-rjC, break out of the
            // loop.
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
   * Attempt to eliminate rotamer pair (ResI-RotA, ResJ-RotC) using (ResI-RotB, ResJ-RotD).
   *
   * @param residues Array of residues.
   * @param i        Index of the first residue.
   * @param riA      Index of the first residue's rotamer to eliminate.
   * @param riB      Index of the first residue's rotamer to use for elimination.
   * @param j        Index of the 2nd residue.
   * @param rjC      Index of the 2nd residue's rotamer to eliminate.
   * @param rjD      Index of the 2nd residue's rotamer to use for elimination.
   * @return Return true if eliminated.
   */
  private boolean goldsteinPairElimination(Residue[] residues, int i, int riA, int riB, int j,
                                           int rjC, int rjD) {

    List<Residue> missedResidues = null;

    // Initialize the Goldstein energy.
    double goldsteinEnergy =
        eE.getSelf(i, riA) + eE.getSelf(j, rjC) + eE.get2Body(i, riA, j, rjC) - eE.getSelf(i, riB)
            - eE.getSelf(j, rjD) - eE.get2Body(i, riB, j, rjD);

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
      logIfRank0(format(
          " Skipping energy comparison due to a missed residue: i %d riA %d riB %d j %d rjC %d rjD %d",
          i, riA, riB, j, rjC, rjD), Level.FINE);
      return false;
    }

    if (goldsteinEnergy > ensembleBuffer) {
      if (missedResidues.isEmpty()) {
        if (eR.eliminateRotamerPair(residues, i, riA, j, rjC, print)) {
          logIfRank0(format(
              "  Pair elimination of [(%8s,%2d),(%8s,%2d)] by [(%8s,%2d),(%8s,%2d)]: %12.4f > %6.4f",
              residues[i].toFormattedString(false, true), riA,
              residues[j].toFormattedString(false, true), rjC,
              residues[i].toFormattedString(false, true), riB,
              residues[j].toFormattedString(false, true), rjD, goldsteinEnergy, ensembleBuffer));
          return true;
        }
      } else {
        logIfRank0(format(
            "  No Pair elimination of [(%8s,%2d),(%8s,%2d)] by [(%8s,%2d),(%8s,%2d)]: %12.4f > %6.4f",
            residues[i].toFormattedString(false, true), riA,
            residues[j].toFormattedString(false, true), rjC,
            residues[i].toFormattedString(false, true), riB,
            residues[j].toFormattedString(false, true), rjD, goldsteinEnergy, ensembleBuffer));
        StringBuilder sb = new StringBuilder();
        for (Residue residueM : missedResidues) {
          sb.append(residueM);
        }
        logIfRank0(format("   due to %s.", sb));
      }
    }
    return false;
  }

  /**
   * Method allows for testing of the elimination criteria by setting parameters appropriately.
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
   * Method allows for testing of the elimination criteria by setting parameters appropriately for a
   * specific test case.
   *
   * @param testSelfEnergyEliminations True when only self energies are calculated; pairs,
   *                                   triples, etc., are assumed to be 0.
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
   * Method allows for testing of the elimination criteria by setting parameters appropriately for a
   * specific test case.
   *
   * @param testPairEnergyEliminations True when only pair energies are calculated; selves,
   *                                   triples, etc., are assumed to be 0.
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
   * Method allows for testing of the elimination criteria by setting parameters appropriately for a
   * specific test case.
   *
   * @param testTripleEnergyEliminations1 True when only triple energies are calculated; selves,
   *                                      pairs, etc., are assumed to be 0.
   * @param testTripleEnergyEliminations2 True when only triple energies are calculated; selves,
   *                                      pairs, etc., are assumed to be 0.
   */
  void setTestTripleEnergyEliminations(int testTripleEnergyEliminations1,
                                       int testTripleEnergyEliminations2) {
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
   * Test the self-energy elimination by setting 2-body and 3-body interactions to zero.
   *
   * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
   */
  private void testSelfEnergyElimination(Residue[] residues) {
    int nRes = residues.length;
    for (int i = 0; i < nRes; i++) {
      Residue resI = residues[i];
      Rotamer[] rotI = resI.getRotamers();
      int nI = rotI.length;
      for (int ri = 0; ri < nI; ri++) {
        for (int j = i + 1; j < nRes; j++) {
          Residue resJ = residues[j];
          Rotamer[] rotJ = resJ.getRotamers();
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
                Rotamer[] rotK = resK.getRotamers();
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
   * Test the elimination criteria by setting self and 3-body interactions to zero.
   *
   * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
   * @param resID    The residue to test.
   */
  private void testPairEnergyElimination(Residue[] residues, int resID) {
    int nRes = residues.length;

    if (resID >= nRes) {
      return;
    }

    for (int i = 0; i < nRes; i++) {
      Residue resI = residues[i];
      Rotamer[] rotI = resI.getRotamers();
      int nI = rotI.length;
      for (int ri = 0; ri < nI; ri++) {
        try {
          eE.setSelf(i, ri, 0, true);
        } catch (Exception e) {
          // catch NPE.
        }
        for (int j = i + 1; j < nRes; j++) {
          Residue resJ = residues[j];
          Rotamer[] rotJ = resJ.getRotamers();
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
                Rotamer[] rotK = resK.getRotamers();
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
   * Test the elimination criteria by setting self and 2-body interactions to zero. Two residues are
   * at fixed rotamers and all rotamer interactions with those two residues are calculated.
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
      Rotamer[] rotI = resI.getRotamers();
      int nI = rotI.length;
      for (int ri = 0; ri < nI; ri++) {
        try {
          eE.setSelf(i, ri, 0, true);
        } catch (Exception e) {
          // catch NPE.
        }
        for (int j = i + 1; j < nRes; j++) {
          Residue resJ = residues[j];
          Rotamer[] rotJ = resJ.getRotamers();
          int nJ = rotJ.length;
          for (int rj = 0; rj < nJ; rj++) {
            /*
             if (i != resID1 && j != resID1) { try { twoBodyEnergy[i][ri][j][rj] = 0.0; } catch
             (Exception e) { // catch NPE. } }
            */
            try {
              eE.set2Body(i, ri, j, rj, 0, true);
            } catch (Exception e) {
              // catch NPE.
            }
            if (threeBodyTerm) {
              for (int k = j + 1; k < nRes; k++) {
                Residue resK = residues[k];
                Rotamer[] rotK = resK.getRotamers();
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

  /**
   * Rotamer Optimization Methods.
   */
  public enum Algorithm {
    INDEPENDENT,  // 1
    ALL,          // 2
    BRUTE_FORCE,  // 3
    WINDOW,       // 4
    BOX;          // 5

    public static Algorithm getAlgorithm(int algorithm) {
      switch (algorithm) {
        case 1:
          return INDEPENDENT;
        case 2:
          return ALL;
        case 3:
          return BRUTE_FORCE;
        case 4:
          return WINDOW;
        case 5:
          return BOX;
        default:
          throw new IllegalArgumentException(
              format(" Algorithm choice was %d, not in range 1-5!", algorithm));
      }
    }
  }

  /**
   * Allows get2BodyDistance to find the shortest distance between any two rotamers or two residues.
   */
  public enum DistanceMethod {
    ROTAMER, RESIDUE
  }

  public enum Direction {
    FORWARD, BACKWARD
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
  //                            while (!(nextLine2 = br.readLine()).contains("Collecting
  // Permutations")) {
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
  //                                            total =
  // singlesIterations.get(singlesIterations.size() - 1);
  //                                        }
  //                                        singlesIterations.add(singlesAmount + total);
  //                                        System.out.println(nextLine4);
  //                                        String[] split2 = nextLine4.split(" ");
  //                                        int pairsAmount = Integer.parseInt(split2[1]);
  //                                        if (pairsIterations.size() > 0) {
  //                                            total2 = pairsIterations.get(pairsIterations.size()
  // - 1);
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
  //            } else if ((line3.contains("-a, 1")) || (line3.contains("-a, 2")) ||
  // (line3.contains("-a, 3")) || (line3.contains("-a, 4"))) {
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
