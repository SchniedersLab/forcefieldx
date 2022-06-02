// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.potential.utils;

import static org.apache.commons.math3.util.FastMath.ceil;
import static ffx.crystal.Crystal.mod;
import static ffx.crystal.SymOp.Tr_0_0_0;
import static ffx.crystal.SymOp.ZERO_ROTATION;
import static ffx.crystal.SymOp.applyCartesianSymOp;
import static ffx.crystal.SymOp.invertSymOp;
import static ffx.numerics.math.DoubleMath.dist2;
import static ffx.numerics.math.MatrixMath.mat4Mat4;
import static ffx.potential.parsers.DistanceMatrixFilter.writeDistanceMatrixRow;
import static ffx.potential.utils.StructureMetrics.momentsOfInertia;
import static ffx.potential.utils.StructureMetrics.radiusOfGyration;
import static ffx.potential.utils.StructureMetrics.radiusOfGyrationComponents;
import static ffx.potential.utils.Superpose.applyRotation;
import static ffx.potential.utils.Superpose.applyTranslation;
import static ffx.potential.utils.Superpose.calculateRotation;
import static ffx.potential.utils.Superpose.calculateTranslation;
import static ffx.potential.utils.Superpose.rmsd;
import static ffx.potential.utils.Superpose.rotate;
import static ffx.potential.utils.Superpose.translate;
import static ffx.utilities.Resources.logResources;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static java.util.Arrays.sort;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.numerics.math.RunningStatistics;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.DistanceMatrixFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.DoubleIndexPair;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

/**
 * Class ProgressiveAlignmentOfCrystals holds the majority of the functionality necessary to quantify
 * crystal similarity following the PAC method.
 *
 * @author Okimasa OKADA, Aaron J. Nessler,  and Michael J. Schnieders
 * @since 1.0
 */
public class ProgressiveAlignmentOfCrystals {

  /**
   * Logger for the ProgressiveAlignmentOfCrystals Class.
   */
  private static final Logger logger = Logger
      .getLogger(ProgressiveAlignmentOfCrystals.class.getName());
  /**
   * SystemFilter containing structures for crystal 1.
   */
  SystemFilter baseFilter;
  /**
   * Number of structures stored in SystemFilter for crystal 1.
   */
  private final int baseSize;
  /**
   * Label for the first crystal.
   */
  private final String baseLabel;
  /**
   * SystemFilter containing structures for crystal 2.
   */
  SystemFilter targetFilter;
  /**
   * Number of structures stored in SystemFilter for crystal 2.
   */
  private final int targetSize;
  /**
   * The amount of work per row for each process.
   */
  private final int numWorkItems;
  /**
   * If this flag is true, then the RMSD matrix is symmetric (e.g., comparing an archive of
   * structures to itself).
   */
  private final boolean isSymmetric;
  /**
   * Label for the second crystal.
   */
  private final String targetLabel;
  /**
   * Label to use for the RSMD logging
   */
  private String rmsdLabel;
  /**
   * Row of RMSD values (length = targetSize).
   */
  public final double[] distRow;
  /**
   * The default restart row is 0. A larger value may be set by the "readMatrix" method if a restart
   * is requested.
   */
  private int restartRow = 0;
  /**
   * The default restart column is 0. Larger values being set by the "readMatrix" method during
   * restarts are not currently supported. A restart must begin from the beginning of a row for
   * simplicity.
   */
  private int restartColumn = 0;
  /**
   * Parallel Java world communicator.
   */
  private final Comm world;
  /**
   * If false, do not use MPI communication.
   */
  private final boolean useMPI;
  /**
   * Number of processes.
   */
  private final int numProc;
  /**
   * Rank of this process.
   */
  private final int rank;
  /**
   * The distances matrix stores a single RSMD value from each process. The array is of size
   * [numProc][1].
   */
  private final double[][] distances;
  /**
   * The assemblies that each process needs for the inner loop are cached.
   */
  private File[] fileCache;
  /**
   * The assemblies that each process needs for the inner loop are cached.
   */
  private String[] nameCache;
  /**
   * The assemblies that each process needs for the inner loop are cached.
   */
  private ArrayList<Bond>[] bondCache;
  /**
   * The assemblies that each process needs for the inner loop are cached.
   */
  private Atom[][] atomCache;
  /**
   * The assemblies that each process needs for the inner loop are cached.
   */
  private ForceField[] forceFieldCache;
  /**
   * The assemblies that each process needs for the inner loop are cached.
   */
  private Crystal[] crystalCache;
  /**
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private final DoubleBuf[] buffers;
  /**
   * Convenience reference to the RMSD value for this process.
   */
  private final double[] myDistances;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private final DoubleBuf myBuffer;
  /**
   * Mass of N asymmetric units.
   */
  private double[] massN;
  /**
   * Sum of masses for one asymmetric unit.
   */
  private double massSum = 0;
  /**
   * Indices of atoms from first crystal to be used in comparison (Z' = 1).
   */
  private int[] comparisonAtoms;
  /**
   * List of each AU index's distance to center of expanded crystal for first crystal.
   */
  private DoubleIndexPair[] molDist1;
  /**
   * Working copy of molDist1, updated when treating a new AU as center.
   */
  private DoubleIndexPair[] molDist1_2;
  /**
   * List of each AU index's distance to center of expanded crystal for second crystal.
   */
  private DoubleIndexPair[] molDist2;
  /**
   * Working copy of molDist2, updated when treating a new AU as center.
   */
  private DoubleIndexPair[] molDist2_2;
  /**
   * List containing coordinates for each unique AU in first crystal.
   */
  private final ArrayList<double[]> tempListXYZ = new ArrayList<>();
  /**
   * Array containing coordinates for each unique AU in first crystal.
   */
  private double[][] baseXYZs;
  /**
   * Array containing coordinates for each unique AU in second crystal.
   */
  private double[][] targetXYZs;
  /**
   * List containing indices for each unique AU in first crystal.
   */
  private final List<Integer> tempListIndices = new ArrayList<>();
  /**
   * Array containing indices for each unique AU in first crystal.
   */
  private Integer[] baseIndices;
  /**
   * Array containing indices for each unique AU in second crystal.
   */
  private Integer[] targetIndices;
  /**
   * Array containing XYZ coordinates for first crystal.
   */
  private double[] baseXYZ;
  /**
   * Array containing XYZ coordinates for first crystal.
   */
  private double[] baseXYZoriginal;
  /**
   * Array containing XYZ coordinates for second crystal.
   */
  private double[] targetXYZ;
  /**
   * Array containing XYZ coordinates for second crystal.
   */
  private double[] targetXYZoriginal;
  /**
   * Coordinates of N AUs from first crystal of the closest match (lowest RMSD).
   */
  private double[] bestBaseNAUs;
  /**
   * Coordinates of N AUs from second crystal of the closest match (lowest RMSD).
   */
  private double[] bestTargetNAUs;
  /**
   * Center of masses (or geometric center) for every AU in first replicates crystal
   */
  private double[][] baseCoM;
  /**
   * Coordinates for central AU of first crystal.
   */
  private double[] baseAU;
  /**
   * Original coordinates for central AU of second crystal.
   */
  private double[] baseAUoriginal;
  /**
   * Coordinates for central 3 AUs of first crystal.
   */
  private double[] base3AUs;
  /**
   * Coordinates for central N AUs of first crystal.
   */
  private double[] baseNAUs;
  /**
   * Center of masses (or geometric center) for every AU in second replicates crystal
   */
  private double[][] targetCoM;
  /**
   * Coordinates for central AU of second crystal.
   */
  private double[] targetAU;
  /**
   * Original coordinates for central AU of second crystal.
   */
  private double[] targetAUoriginal;
  /**
   * Coordinates for central 3 AUs of second crystal.
   */
  private double[] target3AUs;
  /**
   * Coordinates for central N AUs of second crystal.
   */
  private double[] targetNAUs;
  /**
   * Indices and distances of matching AUs between first and second crystal.
   */
  private DoubleIndexPair[] pairedAUs;
  /**
   * Radius of gyration for the best matching clusters (position 0 for first and 1 for second
   * crystal).
   */
  double[] gyrations = new double[2];
  /**
   * Moments of inertia and vectors of application for first crystal.
   */
  double[][] bestBaseMandV = new double[3][4];
  /**
   * Moments of inertia and vectors of application for second crystal.
   */
  double[][] bestTargetMandV = new double[3][4];
  /**
   * Radius of gyration components for first crystal.
   */
  double[] bestBaseRg = new double[3];
  /**
   * Radius of gyration components for second crystal.
   */
  double[] bestTargetRg = new double[3];
  /**
   * Replicates Sym Op applied to asymmetric unit of second system.
   */
  SymOp baseSymOp;
  /**
   * Replicates Sym Op applied to asymmetric unit of first system for best comparison.
   */
  SymOp bestBaseSymOp;
  /**
   * Replicates Sym Op applied to asymmetric unit of second system.
   */
  SymOp targetSymOp;
  /**
   * Replicates Sym Op applied to asymmetric unit of second system for best comparison.
   */
  SymOp bestTargetSymOp;
  /**
   * List of symmetry operators used to create replicates crystal for 1st structure.
   */
  ArrayList<SymOp> baseSymOps = new ArrayList<>();
  /**
   * List of indices created through expansion of replicates crystal for 1st structures.
   */
  ArrayList<Integer> centerB = new ArrayList<>();
  /**
   * Crystal object for first crystal.
   */
  Crystal baseCrystal;
  /**
   * List of symmetry operators used to create replicates crystal for 2nd structure.
   */
  ArrayList<SymOp> targetSymOps = new ArrayList<>();
  /**
   * List of indices created through expansion of replicates crystal for 2nd structures.
   */
  ArrayList<Integer> centerT = new ArrayList<>();
  /**
   * Crystal object for 2nd crystal.
   */
  Crystal targetCrystal;
  /**
   * Boolean to determine if comparison should print symmetry operators used.
   */
  private double printSym;
  /**
   * String builder to attempt to off load work from logger.
   */
  private final StringBuilder stringBuilder = new StringBuilder();
  /**
   * Radius determined based on the volume of one asymmetric unit.
   */
  private double radius;
  /**
   * Symmetry operator used to generate central base molecule updated each iteration.
   */
  private SymOp baseTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
  /**
   * Symmetry operator used to generate central target molecules updated each iteration.
   */
  private SymOp targetTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
  /**
   * Best symmetry operator used to generate central base molecule.
   */
  private SymOp bestBaseTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
  /**
   * Best symmetry operator used to generate central target molecules.
   */
  private SymOp bestTargetTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
  /**
   * AU used from base system.
   */
  int baseSymIndex = -2;
  /**
   * AU used from target system.
   */
  int targetSymIndex = -2;
  /**
   * If molecules between two crystals differ below this tolerance, they are treated as equivalent.
   */
  private static final double MATCH_TOLERANCE = 1.0E-12;

  /**
   * Constructor for the ProgressiveAlignmentOfCrystals class.
   *
   * @param baseFilter SystemFilter containing a set of crystal structures to compare.
   * @param targetFilter SystemFilter containing the other set of crystals to compare.
   */
  public ProgressiveAlignmentOfCrystals(SystemFilter baseFilter,
      SystemFilter targetFilter, boolean isSymmetric) {
    this.baseFilter = baseFilter;
    this.targetFilter = targetFilter;
    this.isSymmetric = isSymmetric;

    // Number of models to be evaluated.
    baseSize = baseFilter.countNumModels();
    baseLabel = getName(baseFilter.getFile().getAbsolutePath());
    targetSize = targetFilter.countNumModels();
    targetLabel = getName(targetFilter.getFile().getAbsolutePath());

    assert !isSymmetric || (baseSize == targetSize);

    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format("\n Conformations for %s: %d", baseLabel, baseSize));
      logger.finer(format(" Conformations for %s: %d", targetLabel, targetSize));
    }

    // Distance matrix to store compared values (dimensions are "human readable" [m x n]).
    distRow = new double[targetSize];

    CompositeConfiguration properties = baseFilter.getActiveMolecularSystem().getProperties();
    useMPI = properties.getBoolean("pj.use.mpi", true);

    if (useMPI) {
      world = Comm.world();
      // Number of processes is equal to world size (often called size).
      numProc = world.size();
      // Each processor gets its own rank (ID of sorts).
      rank = world.rank();
    } else {
      world = null;
      numProc = 1;
      rank = 0;
    }

    // Padding of the target array size (inner loop limit) is for parallelization.
    // Target conformations are parallelized over available nodes.
    // For example, if numProc = 8 and targetSize = 12, then paddedTargetSize = 16.
    int extra = targetSize % numProc;
    int paddedTargetSize = targetSize;
    if (extra != 0) {
      paddedTargetSize = targetSize - extra + numProc;
    }
    numWorkItems = paddedTargetSize / numProc;

    if (numProc > 1) {
      logger.info(format(" Number of MPI Processes:  %d", numProc));
      logger.info(format(" Rank of this MPI Process: %d", rank));
      logger.info(format(" Work per process per row: %d", numWorkItems));
    }

    // Initialize array to -1.0 as -1.0 is not a viable RMSD.
    fill(distRow, -1.0);

    // Each process will complete the following amount of work per row.
    distances = new double[numProc][numWorkItems];

    // Initialize each distance as -2.0.
    for (int i = 0; i < numProc; i++) {
      fill(distances[i], -2.0);
    }

    // DoubleBuf is a wrapper used by MPI Comm methods to transfer data between processors.
    buffers = new DoubleBuf[numProc];
    for (int i = 0; i < numProc; i++) {
      buffers[i] = DoubleBuf.buffer(distances[i]);
    }

    // Convenience reference to the storage for each process.
    myDistances = distances[rank];
    myBuffer = buffers[rank];
  }

  /**
   * Perform single comparison between two crystals (staticAssembly and mobileAssembly).
   *
   * @param name1 Crystal structure (1st or base) that will remain relatively static
   *     (only translations).
   * @param name2 Crystal strucutre (2nd or target) that will rotate to match static
   *     assembly.
   * @param compareAtomsSize Number of active atoms in asymmetric unit of first crystal.
   * @param nAU Number of asymmetric units to compare.
   * @param baseSearchValue Number of anticipated unique entities in 1st system.
   * @param targetSearchValue Number of anticipated unique entities in 2nd system.
   * @param matchTol Tolerance to determine whether two AUs are the same.
   * @param compNum Comparison number based on all file submitted (logging).
   * @param permute Compare all unique AUs between crystals.
   * @param save Save out files of compared crystals.
   * @param machineLearning Save out PDBs and CSVs of compared crystals.
   * @param linkage Criteria to select nearest AUs (0=single, 1=average, 2=complete linkage).
   * @return the computed RMSD.
   */
  private double compare(final File file1, final String name1, final List<Bond> bondList1, final Atom[] atoms1,
      final ForceField forceField1, final File file2, final String name2, final List<Bond> bondList2,
      final Atom[] atoms2, final ForceField forceField2, final int z1, final int z2, final int compareAtomsSize,
      final int nAU, final int baseSearchValue, final int targetSearchValue, final double matchTol, final int compNum,
      final boolean permute, final int save, final boolean machineLearning, final int linkage, final boolean inertia,
      final boolean gyrationComponents) {
    // TODO: Does PAC work for a combination of molecules and polymers?
    // It does not compare them on an individual basis, but can compare AUs as a whole (set zp/zp2 to 1).

    int nCoords = compareAtomsSize * 3;
    // Number of species in expanded crystals.
    int nBaseMols = baseXYZ.length / nCoords;
    int nTargetMols = targetXYZ.length / nCoords;

    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(
          format("\n Base: %s Target: %s", name1, name2));
      stringBuilder.append(format("\n Comparing %3d of %5d in base sphere.\n" +
          " Comparing %3d of %5d in target sphere.\n", nAU, nBaseMols, nAU, nTargetMols));
    }

    arraycopy(baseXYZ, 0, baseAU, 0, nCoords);

    if (molDist1 == null || molDist1.length != nBaseMols) {
      molDist1 = new DoubleIndexPair[nBaseMols];
    }
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(" Prioritize Base System.\n");
    }
    // Center of Masses for crystal 1.
    if (baseCoM == null || baseCoM.length != nBaseMols) {
      baseCoM = new double[nBaseMols][3];
    }
    //Update coordinates for each new comparison.
    centerOfMass(baseCoM, baseXYZ, massN, massSum, compareAtomsSize);
    prioritizeReplicates(baseXYZ, massN, massSum, baseCoM, compareAtomsSize, molDist1, 0, linkage);

    //Used for debugging. Can be removed.
    if (logger.isLoggable(Level.FINEST)) {
      stringBuilder.append(" System 1 distances to center of sphere:\n");
      for (int i = 0; i < nAU; i++) {
        stringBuilder.append(
            format(" %d\t%16.8f\n", molDist1[i].getIndex(), sqrt(molDist1[i].getDoubleValue())));
      }
    }
    // Translate system to the origin.

    boolean useSym = printSym >= 0.0;
    boolean useSave = save > 0;
    arraycopy(targetXYZ, 0, targetAU, 0, nCoords);

    if (molDist2 == null || molDist2.length != nTargetMols) {
      molDist2 = new DoubleIndexPair[nTargetMols];
    }

    // Reorder molDist2 as we shift a different molecule (m) to the center each loop.
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(" Prioritize target system.\n");
    }

    if (targetCoM == null || targetCoM.length != nTargetMols) {
      targetCoM = new double[nTargetMols][3];
    }

    centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
    prioritizeReplicates(targetXYZ, massN, massSum, targetCoM, compareAtomsSize, molDist2, 0,
        linkage);

    if (logger.isLoggable(Level.FINEST)) {
      stringBuilder.append(" System 2 distances to center of sphere:\n");
      for (int i = 0; i < nAU; i++) {
        stringBuilder.append(
            format(" %d\t%16.8f\n", molDist2[i].getIndex(), sqrt(molDist2[i].getDoubleValue())));
      }
    }

    //Determine if AUs in first system are same hand as center most in first (stereoisomer handling).
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(" Search conformations of each crystal:\n");
    }
    tempListXYZ.clear();
    tempListIndices.clear();
    numberUniqueAUs(baseXYZ, molDist1, baseAU.clone(), nCoords, baseSearchValue, permute,
        nBaseMols, massN, matchTol);
    if (baseXYZs == null || baseXYZs.length != tempListXYZ.size()) {
      baseXYZs = new double[tempListXYZ.size()][nCoords];
    }
    tempListXYZ.toArray(baseXYZs);
    int baseLength = tempListIndices.size();
    if (baseIndices == null || baseIndices.length != tempListIndices.size()) {
      baseIndices = new Integer[baseLength];
    }
    tempListIndices.toArray(baseIndices);
    tempListXYZ.clear();
    tempListIndices.clear();
    numberUniqueAUs(targetXYZ, molDist2, targetAU.clone(), nCoords, targetSearchValue, permute,
        nTargetMols, massN, matchTol);
    if (targetXYZs == null || targetXYZs.length != tempListXYZ.size()) {
      targetXYZs = new double[tempListXYZ.size()][nCoords];
    }
    tempListXYZ.toArray(targetXYZs);
    int targetLength = tempListIndices.size();
    if (targetIndices == null || targetIndices.length != targetLength) {
      targetIndices = new Integer[targetLength];
    }
    tempListIndices.toArray(targetIndices);
    int maxLength = max(baseLength, targetLength);
    final boolean reverse = baseLength > targetLength;
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(format("\n Reverse: %b\n" +
              "  %d conformations detected out of %d in base crystal.\n " +
              " %d conformations detected out of %d in target crystal.\n",
          reverse, baseLength, baseSearchValue, targetLength, targetSearchValue));
    }

    // Determine which unique AUs are most similar between crystals
    // Minimum difference between each unique target AU and closest matching base AU.
    // If RMSD_1 is symmetric, could save time by using shorter of two.
    double[] targetBaseDiff = new double[maxLength];
    // Index of the closest matching base AU to each target AU.
    int[] targetBaseIndices = new int[maxLength];
    int[] baseTargetIndices = new int[maxLength];
    double min1 = Double.MAX_VALUE;
    int minBIndex = -1;
    int minTIndex = -1;
    double[][] rotationT = new double[3][3];
    double[] translationB = new double[3];
    double[] translationT = new double[3];

    if (reverse) {
      for (int i = 0; i < baseLength; i++) {
        int minIndex = -1;
        double minDiff = Double.MAX_VALUE;
        int bIndex = baseIndices[i];
        double[] baseXYZ = new double[nCoords];
        for (int j = 0; j < targetLength; j++) {
          // Need to remove translation/rotation from previous check.
          arraycopy(this.baseXYZ, molDist1[bIndex].getIndex() * nCoords, baseXYZ, 0, nCoords);
          int tIndex = targetIndices[j];
          double[] targetXYZ = new double[nCoords];
          arraycopy(this.targetXYZ, molDist2[tIndex].getIndex() * nCoords, targetXYZ, 0, nCoords);
          double[] tempTranB = calculateTranslation(baseXYZ, massN);
          applyTranslation(baseXYZ, tempTranB);
          double[] tempTranT = calculateTranslation(targetXYZ, massN);
          applyTranslation(targetXYZ, tempTranT);
          double[][] tempRotT = calculateRotation(baseXYZ, targetXYZ, massN);
          applyRotation(targetXYZ, tempRotT);
          double value = rmsd(baseXYZ, targetXYZ, massN);
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(
                format(" BaseInd: %3d TargetInd: %3d Diff: %8.4f\n", bIndex, tIndex, value));
          }
          if (value < minDiff) {
            minDiff = value;
            minIndex = j;
            if (minDiff < min1) {
              min1 = minDiff;
              minBIndex = bIndex;
              minTIndex = tIndex;
              baseAU = baseXYZ;
              targetAU = targetXYZ;
              if (nAU == 1 && !permute) {
                // Record values for short circuit.
                translationB = tempTranB;
                translationT = tempTranT;
                rotationT = tempRotT;
                baseNAUs = baseXYZ;
                targetNAUs = targetXYZ;
              }
            }
          }
        }
        targetBaseDiff[i] = minDiff;
        baseTargetIndices[i] = targetIndices[minIndex];
      }
    } else {
      for (int i = 0; i < targetLength; i++) {
        int minIndex = -1;
        double minDiff = Double.MAX_VALUE;
        int tIndex = targetIndices[i];
        double[] targetXYZ = new double[nCoords];
        for (int j = 0; j < baseLength; j++) {
          // Need to remove translation/rotation from previous check.
          arraycopy(this.targetXYZ, molDist2[tIndex].getIndex() * nCoords, targetXYZ, 0, nCoords);
          int bIndex = baseIndices[j];
          double[] baseXYZ = new double[nCoords];
          arraycopy(this.baseXYZ, molDist1[bIndex].getIndex() * nCoords, baseXYZ, 0, nCoords);
          double[] tempTranB = calculateTranslation(baseXYZ, massN);
          applyTranslation(baseXYZ, tempTranB);
          double[] tempTranT = calculateTranslation(targetXYZ, massN);
          applyTranslation(targetXYZ, tempTranT);
          double[][] tempRotT = calculateRotation(baseXYZ, targetXYZ, massN);
          applyRotation(targetXYZ, tempRotT);
          double value = rmsd(baseXYZ, targetXYZ, massN);
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(
                format(" BaseInd: %3d TargetInd: %3d Diff: %8.4f\n", bIndex, tIndex, value));
          }
          if (value < minDiff) {
            minDiff = value;
            minIndex = j;
            if (minDiff < min1) {
              min1 = minDiff;
              minBIndex = bIndex;
              minTIndex = tIndex;
              arraycopy(baseXYZ, 0, baseAU, 0, nCoords);
              arraycopy(targetXYZ, 0, targetAU, 0, nCoords);
              if (nAU == 1 && !permute) {
                // Record values for short circuit.
                translationB = tempTranB;
                translationT = tempTranT;
                rotationT = tempRotT;
                arraycopy(baseAU, 0, baseNAUs, 0, nCoords);
                arraycopy(targetAU, 0, targetNAUs, 0, nCoords);
              }
            }
          }
        }
        targetBaseDiff[i] = minDiff;
        targetBaseIndices[i] = baseIndices[minIndex];
      }
    }
    baseSymIndex = molDist1[minBIndex].getIndex() % z1;
    targetSymIndex = molDist2[minTIndex].getIndex() % z2;
    if (useSym) {
      arraycopy(baseXYZoriginal, centerB.indexOf(baseSymIndex) * nCoords, baseAUoriginal, 0, nCoords);
      arraycopy(targetXYZoriginal, centerT.indexOf(targetSymIndex) * nCoords, targetAUoriginal, 0, nCoords);
      if (logger.isLoggable(Level.FINER)) {
        stringBuilder.append(" \tInitial AU Coords     vs    Center AU Coords");
        for (int i = 0; i < nAU * compareAtomsSize; i++) {
          int k = i * 3;
          stringBuilder.append(format("\n %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f",
                  targetAUoriginal[k], targetAUoriginal[k + 1], targetAUoriginal[k + 2], targetAU[k],
                  targetAU[k + 1], targetAU[k + 2]));
        }
        stringBuilder.append("\n");
      }
    }
    if (nAU == 1 && !permute) {
      // If only comparing one entity we are done. Return values.
      // Rg ranges from 1-->sqrt(3). Normalize to get from 0-->1
      double baseGyration = radiusOfGyration(baseNAUs.clone(), massN);
      gyrations[0] = baseGyration;
      double targetGyration = radiusOfGyration(targetNAUs.clone(), massN);
      gyrations[1] = targetGyration;
      arraycopy(baseNAUs, 0, bestBaseNAUs, 0, nAU * nCoords);
      arraycopy(targetNAUs, 0, bestTargetNAUs, 0, nAU * nCoords);
      if (useSym) {
        baseTransformSymOp = new SymOp(ZERO_ROTATION, translationB);
        bestBaseTransformSymOp = baseTransformSymOp;
        // SymOp.append defaults to rotation then translation. Break apart to have translation then rotation.
        targetTransformSymOp = new SymOp(ZERO_ROTATION, translationT);
        targetTransformSymOp = targetTransformSymOp.append(new SymOp(rotationT, Tr_0_0_0));
        bestTargetTransformSymOp = targetTransformSymOp;
        baseSymOp = baseSymOps.get(centerB.get(molDist1[minBIndex].getIndex()));
        targetSymOp = targetSymOps.get(centerT.get(molDist2[minTIndex].getIndex()));
        if(z1 > 1 || z2 > 1) {
          stringBuilder.append(format("\n Utilizing Sym Ops when Z'>1 requires use of the proper molecule in asymmetric unit.\n" +
                          " Base molecule: %2d of %2d\n" +
                          " Target molecule: %2d of %2d\n",
                  baseSymIndex + 1, z1, targetSymIndex + 1, z2));
        }
        bestBaseSymOp = baseSymOp;
        bestTargetSymOp = targetSymOp;
        StringBuilder alchemicalAtoms = null;
        StringBuilder separation = null;
        for (int i = 0; i < nAU * compareAtomsSize; i++) {
          int index = i * 3;
          double value = rmsd(
              new double[] {bestBaseNAUs[index], bestBaseNAUs[index + 1], bestBaseNAUs[index + 2]},
              new double[] {bestTargetNAUs[index], bestTargetNAUs[index + 1],
                  bestTargetNAUs[index + 2]}, massN);
          if (logger.isLoggable(Level.INFO)) {
            if (printSym < value) {
              // Collect atomic separation distances
              if (separation == null) {
                separation = new StringBuilder("\n Atom  Separation (A)  Description\n");
              }
              separation.append(format(" %4d  %14.6f\n", comparisonAtoms[i] + 1, value));
              // Collect alchemical atoms
              if (alchemicalAtoms == null) {
                alchemicalAtoms = new StringBuilder(" ").append(comparisonAtoms[i] + 1);
              } else {
                alchemicalAtoms.append(",").append(comparisonAtoms[i] + 1);
              }
            }
          }
        }
        if (separation != null) {
          stringBuilder.append(separation);
          stringBuilder.append("\n Alchemical Atoms:\n").append(alchemicalAtoms).append("\n");
        }
        if (useSave) {
          if (machineLearning) {
            saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs, comparisonAtoms, "_c1", 0.000, compNum, save);
            saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs, comparisonAtoms, "_c2", min1, compNum, save);
          } else {
            saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs, comparisonAtoms, "_c1", compNum, save);
            saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs, comparisonAtoms, "_c2", compNum, save);
          }
        }
      }
      if(logger.isLoggable(Level.FINER) && useSym) {
        printSym(compareAtomsSize, file2, name2, bondList2, atoms2, forceField2, save);
      }
      return min1;
    }
    if (logger.isLoggable(Level.FINER)) {
      if (reverse) {
        stringBuilder.append(
            " Minimum RMSD_1 Between Unique Base and Target AUs:\n i bInd RMSD_1    tInd\n");
        int btIndicesLength = baseTargetIndices.length;
        for (int i = 0; i < btIndicesLength; i++) {
          stringBuilder.append(
              format(" %d %4d %4.4f %4d\n", i, molDist1[baseIndices[i]].getIndex(),
                  targetBaseDiff[i], molDist2[baseTargetIndices[i]].getIndex()));
        }
      } else {
        stringBuilder.append(
            " Minimum RMSD_1 Between Unique Target and Base AUs:\n i tInd RMSD_1    bInd\n");
        int tbIndicesLength = targetBaseIndices.length;
        for (int i = 0; i < tbIndicesLength; i++) {
          stringBuilder.append(
              format(" %d %4d %4.4f %4d\n", i, molDist2[targetIndices[i]].getIndex(),
                  targetBaseDiff[i], molDist1[targetBaseIndices[i]].getIndex()));
        }
      }
    }

    // Coordinate arrays to save out structures at the end.
    double bestRMSD = Double.MAX_VALUE;
    if (logger.isLoggable(Level.FINE)) {
      stringBuilder.append(format("\n  Trial     RMSD_1 (%8s)  RMSD_3 (%8s)  %8s  G(r1)    G(r2)\n",
          rmsdLabel, rmsdLabel, rmsdLabel));
    }
    molDist1_2 = new DoubleIndexPair[nBaseMols];
    molDist2_2 = new DoubleIndexPair[nTargetMols];
    // Begin comparison
    // Integer used only for user display logging.
    int currentComparison = 1;
    for (int l = 0; l < baseLength; l++) {
      int index = baseIndices[l];
      final int centerB = molDist1[index].getIndex();
      if (useSym){
        if(logger.isLoggable(Level.FINER)) {
          stringBuilder.append(format("\n centerB: %d", centerB));
          stringBuilder.append(format("\n This centerB: %d", this.centerB.get(centerB)));
          stringBuilder.append(format("\n RemainderB: %d\n", this.centerB.get(centerB) % z1));
          stringBuilder.append(format("\n Base Index: %d", index))
                  .append(format("\n Center: %d", centerB))
                  .append(format("\n SymOp: %d\n", this.centerB.get(centerB)));
        }
        baseSymOp = baseSymOps.get(this.centerB.get(centerB));
        // Reset values for printing symmetry operations
        arraycopy(baseXYZoriginal, 0, baseXYZ, 0, baseXYZoriginal.length);
        centerOfMass(baseCoM, baseXYZoriginal, massN, massSum, compareAtomsSize);
        baseTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
        if (logger.isLoggable(Level.FINER)) {
          stringBuilder.append(matrixToString(baseSymOp.rot, index, "Base Sym Rot"));
          stringBuilder.append(format("\n Base %d Sym Op Translation: ", index))
                  .append(Arrays.toString(baseSymOp.tr)).append("\n");
        }
      }

      //Re-prioritize based on center-most molecule
      prioritizeReplicates(baseXYZ, massN, massSum, baseCoM, compareAtomsSize, molDist1_2, centerB,
          linkage);

      //Translate base system based on center-most molecule
      int baseAUIndex = molDist1_2[0].getIndex() * nCoords;
      arraycopy(baseXYZ, baseAUIndex, baseAU, 0, nCoords);

      // Acquire coordinates based on center 3 molecules
      if (logger.isLoggable(Level.FINEST)) {
        stringBuilder.append(" Base 3 Conformations\n");
      }
      for (int i = 0; i < 3; i++) {
        int baseIndex = molDist1_2[i].getIndex() * nCoords;
        arraycopy(baseXYZ, baseIndex, base3AUs, i * nCoords, nCoords);
      }

      // Acquire coordinates for final comparison
      double maxDist = molDist1[molDist1.length - 1].getDoubleValue();
      if (logger.isLoggable(Level.FINEST)) {
        stringBuilder.append(format("\n System 1 Max Dist: %4.8f", sqrt(maxDist)));
      }
      for (int i = 0; i < nAU; i++) {
        // Only checks for spherical replicates crystal. If one axis is large check will miss selected molecules near short borders.
        if (maxDist - molDist1_2[i].getDoubleValue() < radius) {
          logger.warning(
              " Selected base entity is near the replicates border. Please increase inflation factor (--ni).");
        }
        int molIndex = molDist1_2[i].getIndex() * nCoords;
        arraycopy(baseXYZ, molIndex, baseNAUs, i * nCoords, nCoords);
      }
      for (int m = 0; m < targetLength; m++) {
        if (permute || !reverse && Objects.equals(targetBaseIndices[m], index)
            || reverse && Objects.equals(baseTargetIndices[l], targetIndices[m])) {
          if (useSym) {
            // Reset values for printing symmetry operations
            arraycopy(targetXYZoriginal, 0, targetXYZ, 0, targetXYZoriginal.length);
            targetTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
          }
          // Switch m center most molecules (looking for stereoisomers)
          final int centerT = molDist2[targetIndices[m]].getIndex();
          targetSymOp = targetSymOps.get(this.centerT.get(centerT));

          //Re-prioritize based on central AU if different from first prioritization.
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(" Re-prioritize target system.\n");
          }
          prioritizeReplicates(targetXYZ, massN, massSum, targetCoM, compareAtomsSize, molDist2_2,
              centerT, linkage);
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(" Rotation 1:\n");
          }
          // Isolate the AU closest to center in the second crystal (targetAU).
          arraycopy(targetXYZ, molDist2_2[0].getIndex() * nCoords, targetAU, 0, nCoords);
//          if (useSym) {
//            baseSymOp = baseSymOps.get(this.centerB.get(centerB));
//            targetSymOp = targetSymOps.get(this.centerT.get(centerT));
//            // TODO: finish making printSym compatible for Z'>1.
//              int value = centerT / z2 * z2;
//              stringBuilder.append("\n Sym Op rounded number: ").append(value).append("\n");
//              int newCenter = 0;
//              for(int i = 0; i < z2; i++){
//                newCenter = value + i;
//                targetSymOp = targetSymOps.get(this.centerT.get(newCenter));
//                stringBuilder.append(" Sym Op number: ").append(newCenter).append("\n");
//                if(false){ // Replace with boolean
//                  break;
//                }else{
//
//                }
//              }
//            if (logger.isLoggable(Level.FINER)) {
//              stringBuilder.append(format("\n\n centerT: %d", centerT));
//              stringBuilder.append(format("\n This centerT: %d", this.centerT.get(centerT)));
//              stringBuilder.append(format("\n RemainderT: %d\n", this.centerT.get(centerT) % z2));
//              stringBuilder.append(format("\n Target Index: %d", targetIndices[m]))
//                  .append(format("\n center: %d", centerT))
//                  .append(format("\n targetBaseInd: %d", targetBaseIndices[m]))
//                  .append(format("\n m counter: %d", m))
//                  .append(
//                      format("\n Index of target Index: %d", this.centerT.indexOf(targetIndices[m])))
//                  .append(format("\n Index of center: %d", this.centerT.indexOf(centerT)))
//                  .append(format("\n SymOp: %d\n", this.centerT.get(centerT)));
//              stringBuilder.append(matrixToString(targetSymOp.rot, index, "Target Sym Op Rotation"));
//              stringBuilder.append(format("\n Target %d Sym Op Translation: ", index))
//                  .append(Arrays.toString(targetSymOp.tr)).append("\n");
//              printSym(compareAtomsSize, mobileAssembly, save);
//            }
//          }
          // Translate each system to the origin
          double[] translation = calculateTranslation(baseAU, massN);
          if (useSym) {
            baseTransformSymOp = baseTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(format("\n Base %d to Origin Translation: ", index))
                  .append(Arrays.toString(translation)).append("\n");
            }
          }
          applyTranslation(baseAU, translation);
          applyTranslation(base3AUs, translation);
          applyTranslation(baseNAUs, translation);
          applyTranslation(baseXYZ, translation);

          // Translate to origin
          translation = calculateTranslation(targetAU, massN);
          // Move the center current AU
          applyTranslation(targetAU, translation);
          // Move the entire system
          applyTranslation(targetXYZ, translation);

          // Rotate the target molecule onto the base molecule.
          double[][] rotation = calculateRotation(baseAU, targetAU, massN);
          if (useSym) {
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(rotation, Tr_0_0_0));
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(format("\n Trans target %d sys to origin: \n\t", index))
                  .append(Arrays.toString(translation)).append("\n");
              stringBuilder.append(matrixToString(rotation, index, "1st Rot"));
            }
          }
          applyRotation(targetAU, rotation);
          applyRotation(targetXYZ, rotation);

          // At this point both systems have completed first rotation/translation
          //  Therefore both center-most molecules should be overlapped.
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append("\n Match molecules between systems.\n");
          }

          //Update CoMs with translation
          centerOfMass(baseCoM, baseXYZ, massN, massSum, compareAtomsSize);
          //Update center of masses with the first trans/rot
          centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
          // Need 3 for progressive alignment even if RMSD_1 or RMSD_2 is desired
          pairEntities(Math.max(3, nAU));

          double checkRMSD1 = -3.0;
          double n1RMSD = -4.0;
          if (logger.isLoggable(Level.FINE)) {
            int targetCenterMol = molDist2_2[0].getIndex() * nCoords;
            arraycopy(targetXYZ, targetCenterMol, targetAU, 0, nCoords);
            checkRMSD1 = rmsd(baseAU, targetAU, massN);
            for (int i = 0; i < nAU; i++) {
              int offset = i * nCoords;
              int molIndex = pairedAUs[i].getIndex() * nCoords;
              arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
            }
            n1RMSD = rmsd(baseNAUs, targetNAUs, massN);
            if (logger.isLoggable(Level.FINEST) && useSave) {
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseAU, comparisonAtoms, "_c1_1", compNum, save);
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs, comparisonAtoms, "_c1_1N", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetAU, comparisonAtoms, "_c2_1", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, "_c2_1N", compNum, save);
            }
          }

          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(" Rotation 2:\n");
          }

          translation = calculateTranslation(base3AUs, massN);
          applyTranslation(baseAU, translation);
          applyTranslation(base3AUs, translation);
          applyTranslation(baseNAUs, translation);
          applyTranslation(baseXYZ, translation);
          if (useSym) {
            baseTransformSymOp = baseTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(format("\n Base %d 2nd Translation: ", index))
                  .append(Arrays.toString(translation)).append("\n");
            }
          }

          // Load coordinates for 3 molecules for the target systems
          for (int i = 0; i < 3; i++) {
            int targetIndex = pairedAUs[i].getIndex() * nCoords;
            int coordIndex = i * nCoords;
            arraycopy(targetXYZ, targetIndex, target3AUs, coordIndex, nCoords);
          }

          translation = calculateTranslation(target3AUs, massN);
          applyTranslation(target3AUs, translation);
          applyTranslation(targetNAUs, translation);
          applyTranslation(targetXYZ, translation);

          // Perform 2nd Rotation
          // Calculate the rotation matrix and apply it to the target system.
          rotation = calculateRotation(base3AUs, target3AUs, massN);
          applyRotation(target3AUs, rotation);
          applyRotation(targetNAUs, rotation);
          applyRotation(targetXYZ, rotation);
          if (useSym) {
            applyTranslation(targetAU, translation);
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(rotation, Tr_0_0_0));
            applyRotation(targetAU, rotation);
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(format("\n Target %d 2nd Translation: ", index))
                  .append(Arrays.toString(translation)).append("\n");
              printSym(compareAtomsSize, file2, name2, bondList2, atoms2, forceField2, save);
            }
          }

          double checkRMSD2 = -5.0;
          double n3RMSD = -6.0;
          if (logger.isLoggable(Level.FINE)) {
            checkRMSD2 = rmsd(base3AUs, target3AUs, massN);
            n3RMSD = rmsd(baseNAUs, targetNAUs, massN);
            if (logger.isLoggable(Level.FINEST) && useSave) {
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, base3AUs, comparisonAtoms, "_c1_3", compNum, save);
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs, comparisonAtoms, "_c1_3N", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, target3AUs, comparisonAtoms, "_c2_3", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, "_c2_3N", compNum, save);
            }
          }

          // Update center of masses with the second trans/rot.
          centerOfMass(baseCoM, baseXYZ, massN, massSum, compareAtomsSize);
          centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
          // Rotations 1 and 2 have been completed and both systems should be overlapped
          //  Isolate center most nAU from System 1 and matching molecules from System 2
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(" Match Molecules:\n");
          }
          pairEntities(nAU);

          for (int i = 0; i < nAU; i++) {
            int offset = i * nCoords;
            int molIndex = pairedAUs[i].getIndex() * nCoords;
            arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
          }

          translation = calculateTranslation(baseNAUs, massN);
          if (useSym) {
            baseTransformSymOp = baseTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(format("\n Base %d Final Translation: ", index))
                  .append(Arrays.toString(translation)).append("\n");
            }
          }
          applyTranslation(baseAU, translation);
          applyTranslation(base3AUs, translation);
          applyTranslation(baseNAUs, translation);
          applyTranslation(baseXYZ, translation);

          translation = calculateTranslation(targetNAUs, massN);
          applyTranslation(targetNAUs, translation);
          applyTranslation(targetXYZ, translation);
          rotation = calculateRotation(baseNAUs, targetNAUs, massN);
          if (useSym) {
            applyTranslation(targetAU, translation);
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(rotation, Tr_0_0_0));
            applyRotation(targetAU, rotation);
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(matrixToString(rotation, index, "Target System Final Rotation"));
              stringBuilder.append(format("\n Target %d Final Translation: ", index))
                  .append(Arrays.toString(translation)).append("\n");
              printSym(compareAtomsSize, file2, name2, bondList2, atoms2, forceField2, save);
            }
          }
          if (logger.isLoggable(Level.FINEST) && useSave) {
            saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs, comparisonAtoms, "_c1_N", compNum, save);
            saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, "_c2_N", compNum, save);
          }
          applyRotation(targetNAUs, rotation);
          applyRotation(targetXYZ, rotation);
          double rmsdSymOp = rmsd(baseNAUs, targetNAUs, massN);

          if (logger.isLoggable(Level.FINEST) && useSave) {
            saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs, comparisonAtoms, "_c1_N2", compNum, save);
            saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, "_c2_N2", compNum, save);
          }

          double baseGyration = radiusOfGyration(baseNAUs, massN);
          double targetGyration = radiusOfGyration(targetNAUs, massN);

          if (logger.isLoggable(Level.FINE)) {
            int totalComparisons = (permute) ? baseLength * targetLength
                : maxLength;
            String output = format(
                " %2d of %2d: %7.4f (%8.4f) %7.4f (%8.4f) %8.4f %8.4f %8.4f",
                currentComparison, totalComparisons, checkRMSD1, n1RMSD, checkRMSD2, n3RMSD,
                rmsdSymOp,
                baseGyration, targetGyration);

            if (logger.isLoggable(Level.FINER)) {
              if (reverse) {
                output += format(" b: %2d t: %2d bt: %2d", molDist1[index].getIndex(),
                    molDist2[targetIndices[m]].getIndex(),
                    molDist2[baseTargetIndices[m]].getIndex());
              } else {
                output += format(" b: %2d t: %2d tb: %2d", molDist1[index].getIndex(),
                    molDist2[targetIndices[m]].getIndex(),
                    molDist1[targetBaseIndices[m]].getIndex());
              }
            }
            stringBuilder.append(output).append("\n");
          }

          if (rmsdSymOp < bestRMSD) {
            // Rg ranges from 1-->sqrt(3). Normalize to get from 0-->1
            gyrations[0] = baseGyration;
            gyrations[1] = targetGyration;
            bestRMSD = rmsdSymOp;
            arraycopy(baseNAUs, 0, bestBaseNAUs, 0, nAU * nCoords);
            arraycopy(targetNAUs, 0, bestTargetNAUs, 0, nAU * nCoords);
            if (useSym) {
              bestTargetTransformSymOp = new SymOp(targetTransformSymOp.asMatrix());
              bestTargetSymOp = targetSymOp;
              bestBaseSymOp = baseSymOp;
              bestBaseTransformSymOp = new SymOp(baseTransformSymOp.asMatrix());
            }
          }
          currentComparison++;
        }
      }
    }

    double finalRMSD;
    if (bestRMSD < Double.MAX_VALUE) {
      finalRMSD = bestRMSD;
    } else {
      stringBuilder.append(" This RMSD was filtered out! Try the -p flag or increasing --if.\n");
      // TODO: Double.NaN causes an error in RunningStatistics... Set to -7.0 for now...
      finalRMSD = -7.0;
    }
    if (nAU == 1 && abs(min1 - finalRMSD) > MATCH_TOLERANCE) {
      stringBuilder.append(
          format(" WARNING: Single molecule overlay (%8.4f) does not match PAC RMSD_1 (%8.4f)", min1,
              finalRMSD));
    }
    if (useSave) {
      if (machineLearning) {
        saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs, comparisonAtoms, "_c1", 0.000, compNum, save);
        saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs, comparisonAtoms, "_c2", finalRMSD, compNum,
            save);
      } else {
        saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs, comparisonAtoms, "_c1", compNum, save);
        saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs, comparisonAtoms, "_c2", compNum, save);
      }
    }

    if (inertia) {
      bestBaseMandV = momentsOfInertia(bestBaseNAUs, massN, false, false, true);
      bestTargetMandV = momentsOfInertia(bestTargetNAUs, massN, false, false, true);
    }

    if (gyrationComponents) {
      bestBaseRg = radiusOfGyrationComponents(bestBaseNAUs, massN, true);
      bestTargetRg = radiusOfGyrationComponents(bestTargetNAUs, massN, true);
    }

    return finalRMSD;
  }

  /**
   * Perform default comparison
   *
   * @return RunningStatistics of results.
   */
  public RunningStatistics comparisons() {
    return comparisons(20, 4, 0.1, -1, -1,
        false, false, false, 0, false, 0, false,
        false, false, false, false, 0, -1.0, false, "default");
  }

  /**
   * Compare the crystals within the SystemFilters that were inputted into the constructor of this
   * class.
   *
   * @param nAU Number of asymmetric units to compare.
   * @param inflationFactor Specify safety factor when generating replicates crystal.
   * @param matchTol Tolerance to determine whether two AUs are the same (increases efficiency).
   * @param zPrime Number of asymmetric units in first crystal (-1 attempts detection).
   * @param zPrime2 Number of asymmetric units in second crystal (-1 attempts detection).
   * @param alphaCarbons Perform comparisons on only alpha carbons.
   * @param includeHydrogen Perform comparisons without hydrogen atoms.
   * @param massWeighted Perform comparisons with mass weighted coordinates (center of mass
   *     instead of geometric center).
   * @param crystalPriority Prioritize most dense (0), least dense (1), or first inputted file
   *     (2).
   * @param permute Compare all unique AUs between crystals.
   * @param save Save out files of the resulting superposition.
   * @param restart Try to restart from a previous job.
   * @param write Save out a PAC RMSD file.
   * @param machineLearning Save out CSV files for machine learning input (saves PDBs as well).
   * @param inertia Compute moments of inertia for final clusters.
   * @param gyrationComponents Compute axial components for radius of gyration of final
   *     clusters.
   * @param linkage Prioritize entities based on single, average, or complete linkage.
   * @param printSym Print final symmetry operator used to superimpose mobile assembly onto
   *     static assembly.
   * @param pacFileName The filename to use.
   * @return RunningStatistics Statistics for comparisons performed.
   */
  public RunningStatistics comparisons(int nAU, double inflationFactor, double matchTol, int zPrime,
      int zPrime2, boolean alphaCarbons, boolean includeHydrogen, boolean massWeighted,
      int crystalPriority, boolean permute, int save, boolean restart, boolean write,
      boolean machineLearning, boolean inertia, boolean gyrationComponents,
      int linkage, double printSym, boolean lowMemory, String pacFileName) {
    this.printSym = printSym;
    //TODO: Incorporate graphic user interface (gui: ffx)
    //TODO: Save systems out as original molecule regardless of selection
    //TODO: Handle ring flipping or atoms in equivalent positions (mislabeled atoms)
    //TODO: Handle Z' > 1 for heterogeneous/co-crystals.
    RunningStatistics runningStatistics;
    if (restart) {
      runningStatistics = readMatrix(pacFileName);
      if (runningStatistics == null) {
        runningStatistics = new RunningStatistics();
      }
    } else {
      runningStatistics = new RunningStatistics();
      File file = new File(pacFileName);
      if (file.exists() && file.delete()) {
        logger.info(format(" PAC RMSD file (%s) was deleted.", pacFileName));
        logger.info(" To restart from a previous run, use the '-r' flag.");
      }
    }

    // restartRow and restartColumn are initialized to zero when this class was constructed.
    // They are updated by the "readMatrix" method if a restart is requested.

    // Read ahead to the base starting conformation.
    for (int row = 0; row < restartRow; row++) {
      baseFilter.readNext(false, false, row + 1 >= restartRow);
    }

    MolecularAssembly baseAssembly = baseFilter.getActiveMolecularSystem();

    // Atom arrays from the 1st assembly.
    Atom[] atoms1 = baseAssembly.getAtomArray();
    int nAtoms = atoms1.length;

    // Collect selected atoms.
    ArrayList<Integer> activeIndices = new ArrayList<>();
    determineActiveAtoms(atoms1, activeIndices, alphaCarbons, includeHydrogen);

    if (activeIndices.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD in first crystal.");
      return null;
    }

    int[] comparisonAtoms = activeIndices.stream().mapToInt(i -> i).toArray();

    activeIndices.clear();
    // Atom arrays from the 2nd assembly.
    MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
    Atom[] atoms2 = targetAssembly.getAtomArray();
    int nAtoms2 = atoms2.length;

    // Collect selected atoms.
    determineActiveAtoms(atoms2, activeIndices, alphaCarbons, includeHydrogen);

    if (activeIndices.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD in second crystal.");
      return null;
    }
    int[] comparisonAtoms2 = activeIndices.stream().mapToInt(i -> i).toArray();

    int compareAtomsSize = comparisonAtoms.length;
    int compareAtomsSize2 = comparisonAtoms2.length;

    //Determine number of species within asymmetric unit.
    //TODO: Handle Z' > 1 for heterogeneous/co-crystals.
    int z1 = guessZPrime(zPrime, baseAssembly.getMolecules().size());
    int z2 = guessZPrime(zPrime2, targetAssembly.getMolecules().size());
    // Each ASU contains z * comparisonAtoms species so treat each species individually.
    if (z1 > 1) {
      compareAtomsSize /= z1;
    }
    if (z2 > 1) {
      compareAtomsSize2 /= z2;
    }

    if (compareAtomsSize != compareAtomsSize2) {
      logger.warning(" Selected atom sizes differ between crystals.");
    }
    // When printing Sym Ops it is imperative to find all molecules (matchTol --> 0.0).
    if(z1 > 1 && printSym >= 0.0 || z2 > 1 && printSym >= 0.0){
      matchTol = 0.0;
    }

    baseCrystal = baseAssembly.getCrystal().getUnitCell();
    targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
    // Sohncke groups are non-enantiogenic, so only 1 copy, 2 if mirror exists.
    int baseSearchValue = (baseCrystal.spaceGroup.respectsChirality()) ? z1 : 2 * z1;
    int targetSearchValue = (targetCrystal.spaceGroup.respectsChirality()) ? z2 : 2 * z2;

    // Number of used coordinates for atoms in one AU.
    int nCoords = compareAtomsSize * 3;
    // Remove duplicated atoms from Z' > 1.
    if (this.comparisonAtoms == null) {
      this.comparisonAtoms = new int[compareAtomsSize];
      arraycopy(comparisonAtoms, 0, this.comparisonAtoms, 0, compareAtomsSize);
    }
    if (bestBaseNAUs == null) {
      bestBaseNAUs = new double[nAU * nCoords];
    }
    if (bestTargetNAUs == null) {
      bestTargetNAUs = new double[nAU * nCoords];
    }
    if (baseAU == null) {
      baseAU = new double[nCoords];
    }
    if (targetAU == null) {
      targetAU = new double[nCoords];
    }
    if (baseAUoriginal == null) {
      baseAUoriginal = new double[nCoords];
    }
    if (targetAUoriginal == null) {
      targetAUoriginal = new double[nCoords];
    }
    if (pairedAUs == null || pairedAUs.length != Math.max(3, nAU)) {
      pairedAUs = new DoubleIndexPair[Math.max(3, nAU)];
    }
    if (base3AUs == null) {
      base3AUs = new double[nCoords * 3];
    }
    if (target3AUs == null) {
      target3AUs = new double[nCoords * 3];
    }
    if (baseNAUs == null) {
      baseNAUs = new double[nAU * nCoords];
    }
    if (targetNAUs == null) {
      targetNAUs = new double[nAU * nCoords];
    }
    if (baseSymOps == null) {
      baseSymOps = new ArrayList<>();
    }
    if (targetSymOps == null) {
      targetSymOps = new ArrayList<>();
    }

    int massIndex = 0;
    double[] mass = new double[compareAtomsSize];
    for (Integer value : this.comparisonAtoms) {
      Atom atom = atoms1[value];
      double m = atom.getMass();
      mass[massIndex++] = (massWeighted) ? m : 1.0;
    }

    massIndex = 0;
    double[] mass2 = new double[compareAtomsSize2];
    for (Integer value : comparisonAtoms2) {
      if (massIndex == compareAtomsSize2) {
        //ComparisonAtoms2 contains all indices for the unit cell.
        break;
      }
      Atom atom = atoms2[value];
      double m = atom.getMass();
      mass2[massIndex++] = (massWeighted) ? m : 1.0;
    }

    if (!Arrays.equals(mass, mass2)) {
      logger.warning(" Mass arrays do not match. Check atom ordering.");
      if (logger.isLoggable(Level.FINER)) {
        for (int i = 0; i < compareAtomsSize; i++) {
          if (mass[i] != mass2[i]) {
            logger.finer(format(" Index: %d Mass 1: %4.4f Mass 2: %4.4f", i, mass[i], mass2[i]));
          }
        }
      }
    }

    //Mass of N species
    if (massN == null) {
      massN = (nAU > 3) ? new double[compareAtomsSize * nAU] : new double[compareAtomsSize * 3];
      for (int i = 0; i < nAU; i++) {
        arraycopy(mass, 0, massN, i * compareAtomsSize, compareAtomsSize);
      }
    }
    // Sum of all masses for one species.
    if (massSum == 0) {
      for (int i = 0; i < compareAtomsSize; i++) {
        massSum += massN[i];
      }
    }

    if (machineLearning) {
      save = 1;
    }
    if (save > 2) {
      save = 0;
      logger.info(" Save flag specified incorrectly (1:PDB; 2:XYZ). Not saving files.");
    }

    if(logger.isLoggable(Level.FINER)) {
      if (linkage == 0) {
        logger.finer(" Single linkage will be used.");
      } else if (linkage == 2) {
        logger.finer(" Complete linkage will be used.");
      } else if (linkage == 1) {
        logger.finer(" Average linkage will be used.");
      }
    }
    if(linkage != 1 && linkage !=0 && linkage !=2) {
      logger.warning(
              "Prioritization method specified incorrectly (--pm {0, 1, 2}). Using default of average linkage.");
      linkage = 1;
    }

    // Number of atoms included in the PAC RMSD.
    if (compareAtomsSize == 0 || compareAtomsSize2 == 0) {
      logger.severe(" No atoms were selected for comparison.");
      return runningStatistics;
    }
    logger.info(format("\n %d atoms will be used for the PAC RMSD out of %d in first crystal.",
        compareAtomsSize * z1, nAtoms));
    logger.info(format(" %d atoms will be used for the PAC RMSD out of %d in second crystal.\n",
        compareAtomsSize2 * z2, nAtoms2));

    // Label for logging.
    rmsdLabel = format("RMSD_%d", nAU);

    // Minimum amount of time for a single comparison.
    double minTime = Double.MAX_VALUE;
    double maxTime = -Double.MIN_VALUE;
    // Assemblies used for comparisons. Allows ordering.
    if(!lowMemory) {
      // Store necessary information in arrays for faster access.
      int size = (int)ceil((double)targetSize/(double)numProc);
      atomCache = new Atom[size][nAtoms2];
      fileCache = new File[size];
      nameCache = new String[size];
      bondCache = new ArrayList[size];
      for (int i = 0; i < size; i++) {
        bondCache[i] = new ArrayList<>();
      }
      forceFieldCache = new ForceField[size];
      crystalCache = new Crystal[size];

      // Cache assemblies needed for inner loop.
      for (int column = restartColumn; column < targetSize; column++) {
        int targetRank = column % numProc;
        if (targetRank == rank) {
          int assemblyNum = column / numProc;
          MolecularAssembly currentAssembly = targetFilter.getActiveMolecularSystem();
          Atom[] arrayAtom = currentAssembly.getAtomArray();
          for (int i = 0; i < arrayAtom.length; i++) {
            double[] xyz = new double[3];
            arrayAtom[i].getXYZ(xyz);
            atomCache[assemblyNum][i] = new Atom(i, arrayAtom[i].getName(), arrayAtom[i].getAtomType(), xyz);
          }
          List<Bond> currentBonds = currentAssembly.getBondList();
          for (Bond b : currentBonds) {
            bondCache[assemblyNum].add(new Bond(b.getAtom(0), b.getAtom(1)));
          }
          ForceField currentForcefield = currentAssembly.getForceField();
          fileCache[assemblyNum] = new File(currentAssembly.getFile().getName());
          forceFieldCache[assemblyNum] = new ForceField(currentForcefield.getProperties());
          nameCache[assemblyNum] = currentAssembly.getName();
          Crystal cXtal = currentAssembly.getCrystal().getUnitCell();
          crystalCache[assemblyNum] = new Crystal(cXtal.a, cXtal.b, cXtal.c, cXtal.alpha, cXtal.beta, cXtal.gamma, cXtal.spaceGroup.pdbName);
        }
        targetFilter.readNext(false, false, (column + 1) % numProc == rank);
      }
    }

    if(logger.isLoggable(Level.FINER)) {
      logResources();
    }

    // Loop over conformations in the base assembly.
    for (int row = restartRow; row < baseSize; row++) {
      // Initialize the distance this rank is responsible for to zero.
      fill(myDistances, -8.0);
      int myIndex = 0;
      // Base unit cell for logging.
      baseAssembly = baseFilter.getActiveMolecularSystem();
      baseCrystal = baseAssembly.getCrystal().getUnitCell();
      double baseDensity = baseCrystal.getDensity(baseAssembly.getMass());
      if (baseCrystal == null || baseCrystal.aperiodic()) {
        logger.warning(" " + baseAssembly.getName() + " does not have a crystal.\n");
        continue;
      }
      for (int column = restartColumn; column < targetSize; column++) {
        int targetRank = column % numProc;
        if (targetRank == rank) {
          stringBuilder.setLength(0);
          long time = -System.nanoTime();
          int assemblyNum = column / numProc;
          if (!lowMemory) {
            targetCrystal = crystalCache[assemblyNum];
          } else {
            targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
          }
          if (targetCrystal == null || targetCrystal.aperiodic()) {
            logger.warning(" " + nameCache[assemblyNum] + " does not have a crystal.\n");
            continue;
          }
          double rmsd = -9.0;
          if (isSymmetric && row == column) {
            stringBuilder.append(
                    format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s\n",
                            row + 1, baseCrystal.toShortString(), baseLabel,
                            column + 1, targetCrystal.toShortString(), targetLabel));
            // Fill the diagonal.
            rmsd = 0.0;
            // Log the final result.
            stringBuilder.append(format("\n PAC %s: %12s %7.4f A\n", rmsdLabel, "", rmsd));
          } else if (isSymmetric && row > column) {
            // Do not compute lower triangle values.
            rmsd = -10.0;
          } else {
            stringBuilder.append(
                    format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s\n",
                            row + 1, baseCrystal.toShortString(), baseLabel,
                            column + 1, targetCrystal.toShortString(), targetLabel));

            // Prioritize crystal order based on user specification (High/low density or file order).
            double densityMass = 0;
            if (!lowMemory) {
              for (Atom atom : atomCache[assemblyNum]) {
                densityMass += atom.getMass();
              }
            } else {
              densityMass = targetFilter.getActiveMolecularSystem().getMass();
            }
            double targetDensity = targetCrystal.getDensity(densityMass);
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(
                      format("\n Base Density: %4.4f Target Density: %4.4f\n", baseDensity,
                              targetDensity));
            }
            boolean densityCheck =
                    (crystalPriority == 1) ? baseDensity < targetDensity : baseDensity > targetDensity;
            // Flip system order based on crystalPriority if needed.
            File file1;
            File file2;
            String name1;
            String name2;
            List<Bond> bondList1;
            List<Bond> bondList2;
            ForceField forceField1;
            ForceField forceField2;

            if (densityCheck || crystalPriority == 2) {
              atoms1 = baseAssembly.getAtomArray();
              baseCrystal = baseAssembly.getCrystal();
              file1 = baseAssembly.getFile();
              name1 = baseAssembly.getName();
              bondList1 = baseAssembly.getBondList();
              forceField1 = baseAssembly.getForceField();
              if (!lowMemory) {
                atoms2 = atomCache[assemblyNum];
                targetCrystal = crystalCache[assemblyNum];
                file2 = fileCache[assemblyNum];
                name2 = nameCache[assemblyNum];
                bondList2 = bondCache[assemblyNum];
                forceField2 = forceFieldCache[assemblyNum];
              } else {
                MolecularAssembly mobileAssembly = targetFilter.getActiveMolecularSystem();
                atoms2 = mobileAssembly.getAtomArray();
                targetCrystal = mobileAssembly.getCrystal().getUnitCell();
                file2 = mobileAssembly.getFile();
                name2 = mobileAssembly.getName();
                bondList2 = mobileAssembly.getBondList();
                forceField2 = mobileAssembly.getForceField();
              }
            } else {
              atoms2 = baseAssembly.getAtomArray();
              targetCrystal = baseAssembly.getCrystal();
              file2 = baseAssembly.getFile();
              name2 = baseAssembly.getName();
              bondList2 = baseAssembly.getBondList();
              forceField2 = baseAssembly.getForceField();
              if (!lowMemory) {
                atoms1 = atomCache[assemblyNum];
                baseCrystal = crystalCache[assemblyNum];
                file1 = fileCache[assemblyNum];
                name1 = nameCache[assemblyNum];
                bondList1 = bondCache[assemblyNum];
                forceField1 = forceFieldCache[assemblyNum];
              } else {
                MolecularAssembly staticAssembly = targetFilter.getActiveMolecularSystem();
                atoms1 = staticAssembly.getAtomArray();
                baseCrystal = staticAssembly.getCrystal().getUnitCell();
                file1 = staticAssembly.getFile();
                name1 = staticAssembly.getName();
                bondList1 = staticAssembly.getBondList();
                forceField1 = staticAssembly.getForceField();
              }

              int[] tempAtoms = comparisonAtoms.clone();
              comparisonAtoms = comparisonAtoms2.clone();
              comparisonAtoms2 = tempAtoms.clone();
              int temp = z1;
              z1 = z2;
              z2 = temp;
              temp = baseSearchValue;
              baseSearchValue = targetSearchValue;
              targetSearchValue = temp;
            }

            //Setup for comparison with crystal specific information.
            // Density changes based on mass weighted flag, therefore use volume.
            double baseVolume = baseCrystal.volume / baseCrystal.getNumSymOps() / z1;
            double targetVolume = targetCrystal.volume / targetCrystal.getNumSymOps() / z2;
            double asymmetricUnitVolume = max(abs(baseVolume), abs(targetVolume));
            radius = cbrt(0.75 / PI * asymmetricUnitVolume * nAU * inflationFactor);

            // Estimate a radius that will include desired number of asymmetric units (inflationFactor).
            if (logger.isLoggable(Level.FINER)) {
              logger.finer(
                      format(" Unit Cell Volume:  (Base) %4.2f (Target) %4.2f", baseCrystal.volume,
                              targetCrystal.volume));
              logger.finer(
                      format(" Unit Cell Symm Ops: (Base) %d (Target) %d", baseCrystal.getNumSymOps(),
                              targetCrystal.getNumSymOps()));
              logger.finer(format(" Z': (Base) %d (Target) %d", z1, z2));
              logger.finer(format(" Larger Asymmetric Unit Volume:  %4.2f", asymmetricUnitVolume));
              logger.finer(
                      format(" Larger N Asymmetric Units Volume:  %4.2f", asymmetricUnitVolume * nAU));
              logger.finer(format(" Replicates Cutoff Radius:  %4.2f", radius));
            }

            // When the system was read in, a replicates crystal may have been created to satisfy the cutoff.
            // Retrieve a reference to the unit cell (not the replicates crystal).
            // Here we will use the unit cell, to create a new replicates crystal that may be
            // a different size (i.e. larger).
            //Remove atoms not used in comparisons from the original molecular assembly (crystal 1).
            double[] reducedBaseCoords = reduceSystem(atoms1,
                    comparisonAtoms);

            baseSymOps.clear();
            centerB.clear();
            baseXYZoriginal = generateInflatedSphere(baseCrystal.getUnitCell(), radius,
                    reducedBaseCoords, z1, mass,
                    baseSymOps, centerB);
            int nBaseCoords = baseXYZoriginal.length;
            baseXYZ = new double[nBaseCoords];
            arraycopy(baseXYZoriginal, 0, baseXYZ, 0, baseXYZoriginal.length);
            // Center of Masses for crystal 1.

            //Remove atoms not used in comparisons from the original molecular assembly (crystal 2).
            double[] reducedTargetCoords = reduceSystem(atoms2,
                    comparisonAtoms2);

            targetSymOps.clear();
            centerT.clear();
            targetXYZoriginal = generateInflatedSphere(targetCrystal.getUnitCell(), radius,
                    reducedTargetCoords, z2, mass2,
                    targetSymOps, centerT);
            int nTargetCoords = targetXYZoriginal.length;
            targetXYZ = new double[nTargetCoords];
            arraycopy(targetXYZoriginal, 0, targetXYZ, 0, nTargetCoords);

            if (logger.isLoggable(Level.FINEST) && save > 0) {
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseXYZoriginal, comparisonAtoms, "_c1_rep", -1, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetXYZoriginal, comparisonAtoms2, "_c2_rep", -1, save);
            }
            // Compute the PAC RMSD.
            rmsd = compare(file1, name1, bondList1, atoms1, forceField1, file2, name2, bondList2, atoms2, forceField2, z1, z2, compareAtomsSize, nAU,
                    baseSearchValue, targetSearchValue, matchTol, row * targetSize + column,
                    permute, save, machineLearning, linkage, inertia, gyrationComponents);
            time += System.nanoTime();
            double timeSec = time * 1.0e-9;
            // Record the fastest comparison.
            if (timeSec < minTime) {
              minTime = timeSec;
            }
            if (timeSec > maxTime) {
              maxTime = timeSec;
            }
            // Log the final result.
            double avgGyration = (gyrations[0] + gyrations[1]) / 2;
            double diff = max(gyrations[0], gyrations[1]) - avgGyration;
            stringBuilder.append(
                    format("\n PAC %7s: %12s %7.4f A (%5.3f sec) G(r) %7.4f A +/- %7.4f\n", rmsdLabel, "",
                            rmsd,
                            timeSec, avgGyration, diff));
            if (inertia) {
              stringBuilder.append("\n Moments of Inertia and Principle Axes for first crystal:\n" +
                      "  Moments (amu Ang^2): \t X-, Y-, and Z-Components of Axes:\n");
              for (int i = 1; i < 4; i++) {
                stringBuilder.append(
                        format("  %16.3f %12.6f %12.6f %12.6f\n", bestBaseMandV[i - 1][0],
                                bestBaseMandV[0][i], bestBaseMandV[1][i], bestBaseMandV[2][i]));
              }
              stringBuilder.append(" Moments of Inertia and Principle Axes for second crystal:\n" +
                      "  Moments (amu Ang^2): \t X-, Y-, and Z-Components of Axes:\n");
              for (int i = 1; i < 4; i++) {
                stringBuilder.append(
                        format("  %16.3f %12.6f %12.6f %12.6f\n", bestTargetMandV[i - 1][0],
                                bestTargetMandV[0][i], bestTargetMandV[1][i], bestTargetMandV[2][i]));
              }
            }
            if (gyrationComponents) {
              // Crystal 1 metrics
              double Rmax = max(max(bestBaseRg[0], bestBaseRg[1]), bestBaseRg[2]);
              Rmax *= Rmax;
              double Rmed = max(min(bestBaseRg[0], bestBaseRg[1]), bestBaseRg[2]);
              Rmed *= Rmed;
              double Rmin = min(min(bestBaseRg[0], bestBaseRg[1]), bestBaseRg[2]);
              Rmin *= Rmin;
              double Rg = Rmax + Rmed + Rmin;
              // 0.0 indicates a perfect sphere, and 1.0 indicates all atoms along a single axis.
              double asphericity = Rmax - 0.5 * (Rmin + Rmed);
              double acylindricity = Rmed - Rmin;
              double anisotropy = (asphericity * asphericity + 0.75 * acylindricity * acylindricity);
              double totalMass = Arrays.stream(mass).sum();
              // Density changes based on mass weighted flag.
              double volume = totalMass / max(baseCrystal.getDensity(massSum),
                      targetCrystal.getDensity(massSum));
              double targetRadius = cbrt(0.75 / PI * volume);
              stringBuilder.append(format(" Base Radius Metric: %7.4f", avgGyration / targetRadius));
              if (logger.isLoggable(Level.FINER)) {
                stringBuilder.append(
                        format("\n Target Base Radius:           %9.3f Å\n", targetRadius));
                stringBuilder.append(
                        format(" Asphericity:   (%6.3f Å^2) %9.3f\n", asphericity, asphericity / Rmax));
                stringBuilder.append(format(" Acylindricity: (%6.3f Å^2) %9.3f\n", acylindricity,
                        acylindricity / Rmax));
                stringBuilder.append(
                        format(" Anisotropy:    (%6.3f Å) %9.3f\n", sqrt(sqrt(anisotropy)),
                                anisotropy / (Rg * Rg)));
                stringBuilder.append(format("  Ryz %9.3fÅ  Rxz %9.3fÅ  Rxy %9.3fÅ\n",
                        bestBaseRg[0], bestBaseRg[1], bestBaseRg[2]));
              }
              // Crystal 2 metrics
              Rmax = max(max(bestTargetRg[0], bestTargetRg[1]), bestTargetRg[2]);
              Rmax *= Rmax;
              Rmed = max(min(bestTargetRg[0], bestTargetRg[1]), bestTargetRg[2]);
              Rmed *= Rmed;
              Rmin = min(min(bestTargetRg[0], bestTargetRg[1]), bestTargetRg[2]);
              Rmin *= Rmin;
              Rg = Rmax + Rmed + Rmin;
              // 0.0 indicates a perfect sphere, and 1.0 indicates all atoms along a single axis.
              asphericity = Rmax - 0.5 * (Rmin + Rmed);
              acylindricity = Rmed - Rmin;
              anisotropy = (asphericity * asphericity + 0.75 * acylindricity * acylindricity);
              totalMass = Arrays.stream(mass).sum();
              volume = totalMass / targetDensity;
              targetRadius = cbrt(0.75 / PI * volume);
              if (logger.isLoggable(Level.FINER)) {
                stringBuilder.append(
                        format("\n Target Target Radius:           %9.3f Å\n", targetRadius));
                stringBuilder.append(
                        format(" Asphericity:   (%6.3f Å) %9.3f\n", asphericity, asphericity / Rmax));
                stringBuilder.append(format(" Acylindricity: (%6.3f Å) %9.3f\n", acylindricity,
                        acylindricity / Rmax));
                stringBuilder.append(
                        format(" Anisotropy:    (%6.3f Å) %9.3f\n", sqrt(sqrt(anisotropy)),
                                anisotropy / (Rg * Rg)));
                stringBuilder.append(format("  Ryz %9.3fÅ  Rxz %9.3fÅ  Rxy %9.3fÅ\n",
                        bestTargetRg[0], bestTargetRg[1], bestTargetRg[2]));
              }
            }
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(
                      format("\n Gyration Crystal 1 (%s): %7.4f Crystal 2 (%s): %7.4f.\n",
                              name1,
                              gyrations[0], name2, gyrations[1]));
            }
            if (printSym >= 0.0) {
              bestTargetTransformSymOp = bestTargetTransformSymOp.prepend(bestTargetSymOp);
              // Apply inverse of base operations:
              // For orthogonal matrices the inverse matrix = the transpose. True iff det(A)== +/- 1.
              // No inverse if det(A)==0.
              bestBaseTransformSymOp = bestBaseTransformSymOp.prepend(bestBaseSymOp);
              bestTargetTransformSymOp = bestTargetTransformSymOp.append(
                      invertSymOp(bestBaseTransformSymOp));

              stringBuilder.append(
                      format("\n Sym Op to move %s onto %s:\n", name2,
                              name1)).append(bestTargetTransformSymOp).append("\n");
              stringBuilder.append(asMatrixString(bestTargetTransformSymOp));
              SymOp inverted = invertSymOp(bestTargetTransformSymOp);
              stringBuilder.append(
                      format("\n\n Inverted Sym Op to move %s onto %s:\n", name1,
                              name2)).append(inverted).append("\n");
              stringBuilder.append(asMatrixString(inverted));
              if (logger.isLoggable(Level.FINE)) {
                stringBuilder.append("\n Sym Op Application from Scratch:");
                double[] symMol = new double[nCoords];
                for (int i = 0; i < compareAtomsSize; i++) {
                  int k = i * 3;
                  double[] xyz = new double[]{targetAUoriginal[k], targetAUoriginal[k + 1],
                          targetAUoriginal[k + 2]};
                  applyCartesianSymOp(xyz, xyz, bestTargetTransformSymOp);
                  symMol[k] = xyz[0];
                  symMol[k + 1] = xyz[1];
                  symMol[k + 2] = xyz[2];
                  stringBuilder.append(format("\n %8.3f %8.3f %8.3f moved to %8.3f %8.3f %8.3f compared to %8.3f %8.3f %8.3f",
                          targetAUoriginal[k], targetAUoriginal[k + 1], targetAUoriginal[k + 2],
                          symMol[k], symMol[k + 1], symMol[k + 2],
                          baseAUoriginal[k], baseAUoriginal[k + 1], baseAUoriginal[k + 2]));
                }
                if (save > 0) {
                  saveAssembly(file2, name2, bondList2, atoms2, forceField2, symMol, comparisonAtoms, "_moved", 0, save);
                }

                double value = rmsd(baseAUoriginal, symMol, massN);
                stringBuilder.append(
                        format("\n\n Sym Op RMSD: %8.8f ", value));
                double[] symMol2 = new double[nCoords];
                stringBuilder.append("\n\n Sym Op Inverse Application:");
                for (int i = 0; i < compareAtomsSize; i++) {
                  int k = i * 3;
                  double[] xyz = new double[]{symMol[k], symMol[k + 1],
                          symMol[k + 2]};
                  applyCartesianSymOp(xyz, xyz, inverted);
                  symMol2[k] = xyz[0];
                  symMol2[k + 1] = xyz[1];
                  symMol2[k + 2] = xyz[2];
                  stringBuilder.append(format("\n %8.3f %8.3f %8.3f moved to %8.3f %8.3f %8.3f compared to %8.3f %8.3f %8.3f",
                          symMol[k], symMol[k + 1], symMol[k + 2],
                          symMol2[k], symMol2[k + 1], symMol2[k + 2],
                          targetAUoriginal[k], targetAUoriginal[k + 1], targetAUoriginal[k + 2]));
                }
              }
            }
          }
          myDistances[myIndex] = rmsd;
          myIndex++;
          if(stringBuilder.length() > 0) {
            logger.info(stringBuilder.toString());
          }
        }
        if (lowMemory) {
          targetFilter.readNext(false, false, (column + 1) % numProc == rank);
        }
      }
      restartColumn = 0;
      if(lowMemory) {
        targetFilter.readNext(true, false, 0 % numProc == rank);
      }
      baseFilter.readNext(false, false);

      // Gather RMSDs for this row.
      gatherRMSDs(row, runningStatistics);

      // Write out this row.
      if (rank == 0 && write) {
        int firstColumn = (isSymmetric) ? row : 0;
        writeDistanceMatrixRow(pacFileName, distRow, firstColumn);
      }
    }

    if (minTime < Double.MAX_VALUE) {
      logger.info(format("\n Minimum PAC Time: %7.4f", minTime));
    }
    if (maxTime > Double.MIN_VALUE && maxTime != minTime) {
      logger.info(format(" Maximum PAC Time: %7.4f", maxTime));
    }

    baseFilter.closeReader();
    targetFilter.closeReader();

    if (rank == 0 || logger.isLoggable(Level.FINER)) {
      logger.info(format(" RMSD Minimum:  %8.6f", runningStatistics.getMin()));
      logger.info(format(" RMSD Maximum:  %8.6f", runningStatistics.getMax()));
      logger.info(format(" RMSD Mean:     %8.6f", runningStatistics.getMean()));
      double variance = runningStatistics.getVariance();
      if (!Double.isNaN(variance)) {
        logger.info(format(" RMSD Variance: %8.6f", variance));
      }
    }

    // Return distMatrix for validation if this is for the test script
    return runningStatistics;
  }

  /**
   * Print a Sym Op matrix as a continued line string.
   *
   * @param symOp Symmetry operation to print.
   * @return Continued line string.
   */
  private static String asMatrixString(SymOp symOp) {
    double[][] values = symOp.asMatrix();
    return format(
        " %14.8f %14.8f %14.8f \\\n" +
            " %14.8f %14.8f %14.8f \\\n" +
            " %14.8f %14.8f %14.8f \\\n" +
            " %14.8f %14.8f %14.8f ",
        values[0][0], values[0][1], values[0][2],
        values[1][0], values[1][1], values[1][2],
        values[2][0], values[2][1], values[2][2],
        values[0][3], values[1][3], values[2][3]);
  }

  /**
   * Add a value to a list of doubles if its difference to all listed values is greater than the
   * tolerance.
   *
   * @param values List of values already found.
   * @param value Potential new value if it is not already in list.
   * @param tol Tolerance that determine whether values are equal.
   */
  private static boolean addLooseUnequal(List<Double> values, double value, double tol) {
    boolean found = false;
    for (Double dbl : values) {
      if (abs(dbl - value) < tol) {
        found = true;
        break;
      }
    }
    if (!found) {
      values.add(value);
    }
    // If value is not found it was added. Otherwise, not added.
    return !found;
  }

  /**
   * Accumulate rotations (matrix multiplication)
   *
   * @param rot Rotation matrix to add.
   * @param totalTransform Array to be updated with rotation (4x4).
   * @param prepend If true prepend translation, false append to end.
   */
  public void addRotation(double[][] rot, double[][] totalTransform, boolean prepend) {
    double[][] transform = new double[][] {
        {rot[0][0], rot[0][1], rot[0][2], 0.0},
        {rot[1][0], rot[1][1], rot[1][2], 0.0},
        {rot[2][0], rot[2][1], rot[2][2], 0.0},
        {0.0, 0.0, 0.0, 1.0}};
    transform =
        (prepend) ? mat4Mat4(totalTransform, transform) : mat4Mat4(transform, totalTransform);
    for (int i = 0; i < totalTransform.length; i++) {
      arraycopy(transform[i], 0, totalTransform[i], 0, totalTransform[i].length);
    }
  }

  /**
   * Accumulate translations (matrix multiplication)
   *
   * @param translation Translation matrix to add.
   * @param totalTransform Array to be updated with translation (4x4).
   * @param prepend If true prepend translation, false append to end.
   */
  public static void addTranslation(double[] translation, double[][] totalTransform,
      boolean prepend) {
    double[][] transform = new double[][] {
        {1.0, 0.0, 0.0, translation[0]},
        {0.0, 1.0, 0.0, translation[1]},
        {0.0, 0.0, 1.0, translation[2]},
        {0.0, 0.0, 0.0, 1.0}};
    transform =
        (prepend) ? mat4Mat4(totalTransform, transform) : mat4Mat4(transform, totalTransform);
    for (int i = 0; i < totalTransform.length; i++) {
      arraycopy(transform[i], 0, totalTransform[i], 0, totalTransform[i].length);
    }
  }

  /**
   * Calculate the center of mass for a given set of masses for the asymmetric unit and coordinates
   * (xyz)
   *
   * @param centersOfMass Returned center of mass for each asymmetric unit
   * @param coords Coordinates of every atom in system.
   * @param mass Masses of each atom in asymmetric unit.
   * @param massSum Sum of masses within asymmetric unit.
   * @param nAtoms Number of coordinates in an entity.
   */
  private static void centerOfMass(double[][] centersOfMass, double[] coords, double[] mass,
      double massSum, int nAtoms) {
    int size = centersOfMass.length;
    for (int i = 0; i < size; i++) {
      int molIndex = i * nAtoms * 3;
      double[] comI = new double[3];
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = j * 3;
        double m = mass[j];
        for (int k = 0; k < 3; k++) {
          comI[k] += coords[molIndex + atomIndex + k] * m;
        }
      }
      for (int j = 0; j < 3; j++) {
        comI[j] /= massSum;
      }
      centersOfMass[i] = comI;
    }
  }

  /**
   * Deep copy of values from one array to another.
   *
   * @param newArray Destination of values.
   * @param oldArray Current location of values.
   */
  private static void copyArrayValues(double[][] newArray, double[][] oldArray) {
    for (int i = 0; i < newArray.length; i++) {
      arraycopy(oldArray[i], 0, newArray[i], 0, newArray[i].length);
    }
  }

  /**
   * Determine the indices of the atoms from the assembly that are active for this comparison.
   *
   * @param atoms Atoms we potentially wish to use in comparison.
   * @param indices Array list containing atom indices that will be used for this comparison.
   * @param alphaCarbons Boolean whether to include only alpha carbons/nitrogens.
   * @param includeHydrogen Boolean whether to include hydrogens.
   */
  private static void determineActiveAtoms(Atom[] atoms, ArrayList<Integer> indices,
      boolean alphaCarbons,
      boolean includeHydrogen) {
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        if (alphaCarbons) {
          String atomName = atom.getName();
          int atomAtNum = atom.getAtomicNumber();
          boolean proteinCheck = atomName.equals("CA") && atomAtNum == 6;
          boolean aminoCheck = (atomName.equals("N1") || atomName.equals("N9")) && atomAtNum == 7;
          if (proteinCheck || aminoCheck) {
            indices.add(i);
          }
        } else if (includeHydrogen || !atom.isHydrogen()) {
          indices.add(i);
        }
      }
      // Reset all atoms to active once the selection is recorded.
      atom.setActive(true);
    }
  }

  /**
   * Determine translation needed to move center from x2 to x1.
   *
   * @param x1 Coordinates to whose center we wish to move
   * @param x2 Coordinates whose center we wish to move
   * @param mass Mass of the system.
   * @return double[] translation needed to move center of x2 to x1.
   */
  public static double[] determineTranslation(double[] x1, double[] x2, final double[] mass) {
    assert (x1.length == x2.length);
    assert (x1.length % 3 == 0);
    assert (x1.length / 3 == mass.length);

    double x1mid = 0.0;
    double y1mid = 0.0;
    double z1mid = 0.0;
    double x2mid = 0.0;
    double y2mid = 0.0;
    double z2mid = 0.0;
    double norm = 0.0;
    int n = x1.length / 3;

    for (int i = 0; i < n; i++) {
      int k = 3 * i;
      double weigh = mass[i];
      x1mid += x1[k] * weigh;
      y1mid += x1[k + 1] * weigh;
      z1mid += x1[k + 2] * weigh;
      x2mid += x2[k] * weigh;
      y2mid += x2[k + 1] * weigh;
      z2mid += x2[k + 2] * weigh;
      norm += weigh;
    }
    x1mid /= norm;
    y1mid /= norm;
    z1mid /= norm;
    x2mid /= norm;
    y2mid /= norm;
    z2mid /= norm;

    return new double[] {x1mid - x2mid, y1mid - y2mid, z1mid - z2mid};
  }

  /**
   * This method calls <code>world.gather</code> to collect numProc PAC RMSD values.
   *
   * @param row Current row of the PAC RMSD matrix.
   * @param runningStatistics Stats for the RMSDs.
   */
  private void gatherRMSDs(int row, RunningStatistics runningStatistics) {
    if (useMPI) {
      try {
        if (logger.isLoggable(Level.FINER)) {
          logger.finer(" Receiving MPI results.");
        }
        // TODO: Node 0 is the only process that updates. Remove from other nodes.
        world.gather(0, myBuffer, buffers);
        if (rank == 0) {
          for (int workItem = 0; workItem < numWorkItems; workItem++) {
            for (int proc = 0; proc < numProc; proc++) {
              int column = numProc * workItem + proc;
              // Do not include padded results.
              if (column < targetSize) {
                distRow[column] = distances[proc][workItem];
                if (!isSymmetric) {
                  runningStatistics.addValue(distRow[column]);
                } else if (column >= row) {
                  // Only collect stats for the upper triangle.
                  runningStatistics.addValue(distRow[column]);
                }
                if (logger.isLoggable(Level.FINER)) {
                  logger.finer(format(" %d %d %14.8f", row, column, distances[proc][workItem]));
                }
              }
            }
          }
        } else {
          for (int workItem = 0; workItem < numWorkItems; workItem++) {
            int column = numProc * workItem + rank;
            // Do not include padded results.
            if (column < targetSize) {
              distRow[column] = distances[rank][workItem];
              if (!isSymmetric) {
                runningStatistics.addValue(distRow[column]);
              } else if (column >= row) {
                // Only collect stats for the upper triangle.
                runningStatistics.addValue(distRow[column]);
              }
              if (logger.isLoggable(Level.FINER)) {
                logger.finer(format(" %d %d %14.8f", row, column, distances[rank][workItem]));
              }
            }
          }
        }
      } catch (Exception e) {
        logger.severe(" Exception collecting distance values." + e);
      }
    } else {
      for (int i = 0; i < targetSize; i++) {
        distRow[i] = myDistances[i];
        if (!isSymmetric) {
          runningStatistics.addValue(distRow[i]);
        } else if (i >= row) {
          // Only collect stats for the upper triangle.
          runningStatistics.addValue(distRow[i]);
        }
      }
    }
  }

  /**
   * Generate and expanded sphere of asymmetric unit with the intention of observing a crystals'
   * distribution of replicates rather to facilitate comparisons that go beyond lattice parameters.
   *
   * @param unitCell Crystal to define expansion.
   * @param radius Minimum radius desired for replicates crystal.
   * @param reducedCoords Coordinates of asymmetric unit we wish to expand.
   * @param mass Masses for atoms within reduced asymmetric unit.
   * @return double[] containing the coordinates for the expanded crystal.
   */
  private static double[] generateInflatedSphere(Crystal unitCell, double radius,
      double[] reducedCoords, int zPrime,
      double[] mass, ArrayList<SymOp> symOps, ArrayList<Integer> center) {
    int nAtoms = reducedCoords.length / 3;
    int zAtoms = nAtoms / zPrime;
    // Collect asymmetric unit atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];
    for (int i = 0; i < nAtoms; i++) {
      int atomIndex = i * 3;
      x[i] = reducedCoords[atomIndex];
      y[i] = reducedCoords[atomIndex + 1];
      z[i] = reducedCoords[atomIndex + 2];
    }

    Crystal replicatesCrystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, radius * 2.0);

    // Symmetry coordinates for each molecule in replicates crystal
    int nSymm = replicatesCrystal.getNumSymOps();

    int numEntities = nSymm * zPrime;

    double[][] xS = new double[nSymm][nAtoms];
    double[][] yS = new double[nSymm][nAtoms];
    double[][] zS = new double[nSymm][nAtoms];
    // Cartesian center of each molecule
    double[][] centerMolsCart = new double[numEntities][3];

    // Loop over replicate crystal SymOps
    List<SymOp> inflatedSymOps = replicatesCrystal.spaceGroup.symOps;
    for (int iSym = 0; iSym < nSymm; iSym++) {
      SymOp symOp = inflatedSymOps.get(iSym);
      //Convert sym op into cartesian for later use.
      double[][] rot = new double[3][3];
      replicatesCrystal.getTransformationOperator(symOp, rot);
      // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
      replicatesCrystal.applySymOp(nAtoms, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
      for (int zp = 0; zp < zPrime; zp++) {
        int symIndex = zp * zAtoms;
        // Compute center-of-mass (CoM) for Cartesian coordinates
        double[] centerOfMass = new double[3];
        double totalMass = 0.0;
        for (int i = 0; i < zAtoms; i++) {
          double m = mass[i];
          centerOfMass[0] += xS[iSym][symIndex + i] * m;
          centerOfMass[1] += yS[iSym][symIndex + i] * m;
          centerOfMass[2] += zS[iSym][symIndex + i] * m;
          totalMass += m;
        }
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;

        double[] translate = moveIntoCrystal(replicatesCrystal, centerOfMass);
        //Convert sym op into cartesian for later use.
        double[] trans = symOp.tr.clone();
        replicatesCrystal.toCartesianCoordinates(trans, trans);
        for (int i = 0; i < translate.length; i++) {
          trans[i] += translate[i];
        }
        symOps.add(new SymOp(rot, trans));
        for (int i = 0; i < zAtoms; i++) {
          xS[iSym][symIndex + i] += translate[0];
          yS[iSym][symIndex + i] += translate[1];
          zS[iSym][symIndex + i] += translate[2];
        }

        // Save CoM cartesian coordinates
        centerMolsCart[iSym * zPrime + zp] = centerOfMass;
      }
    }

    //Determine molecular distances to "center" of sphere.
    //  In PACCOM the center is the geometric average of coordinates.
    //  In FFX the center is the center of the replicates crystal.
    DoubleIndexPair[] molsDists = new DoubleIndexPair[numEntities];
    double[] cartCenter = new double[3];

    // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
    // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
    double[] fracCenter = {0.5, 0.5, 0.5};
    replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

    if (logger.isLoggable(Level.FINER)) {
      if (logger.isLoggable(Level.FINEST)) {
        logger.finer(" Replicates crystal " + replicatesCrystal);
      }
      logger.finer(format(" Replicates Volume: %8.4f", replicatesCrystal.volume));
      logger.finer(format(" Expanded Crystal Center: %14.8f %14.8f %14.8f",
          cartCenter[0], cartCenter[1], cartCenter[2]));
    }

    for (int i = 0; i < numEntities; i++) {
      // Then compute Euclidean distance from Cartesian center of the replicates cell
      molsDists[i] = new DoubleIndexPair(i, dist2(cartCenter, centerMolsCart[i]));
    }

    // Sort the molecules by their distance from the center.
    // Note that the smallest distances are first in the array after the sort.
    sort(molsDists);
    for (DoubleIndexPair molsDist : molsDists) {
      center.add(molsDist.getIndex());
    }

    if (logger.isLoggable(Level.FINEST)) {
      logger.finest("\n Copy  SymOp        Distance");
    }
    double[] systemCoords = new double[numEntities * zAtoms * 3];
    for (int n = 0; n < nSymm; n++) {
      for (int zp = 0; zp < zPrime; zp++) {
        int index = n * zPrime + zp;
        // Current molecule
        int iSym = molsDists[index].getIndex();
        double distance = molsDists[index].getDoubleValue();
        if (logger.isLoggable(Level.FINEST) && n < 50) {
          logger.finest(format(" %4d  %5d  %14.8f", index, iSym, sqrt(distance)));
        }

        // Create a new set of Atoms for each SymOp molecule
        for (int i = 0; i < zAtoms; i++) {
          int symIndex = index * zAtoms * 3;
          int atomIndex = i * 3;
          int conversion = iSym / zPrime;
          int shift = (iSym % zPrime) * zAtoms;
          systemCoords[symIndex + atomIndex] = xS[conversion][i + shift];
          systemCoords[symIndex + atomIndex + 1] = yS[conversion][i + shift];
          systemCoords[symIndex + atomIndex + 2] = zS[conversion][i + shift];
        }
      }
    }
    return systemCoords;
  }

  /**
   * Try to automatically determine number of species in asymmetric unit (only works for molecules).
   *
   * @param zPrime User input overrides detection method.
   * @param numEntities Number of species detected.
   * @return Number of expected species in asymmetric unit.
   */
  private static int guessZPrime(int zPrime, int numEntities) {
    int z = (zPrime > 0) ? zPrime : Math.max(numEntities, 1);
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" Number of species in asymmetric unit (Z Prime): %d", z));
    }
    return z;
  }

  /**
   * Parse values of a matrix into a string.
   *
   * @param matrix Values to present
   * @param index Identifier
   * @param description Identifier
   * @return String of values.
   */
  public static String matrixToString(double[][] matrix, int index, String description) {
    StringBuilder value = new StringBuilder(format("\n %s %d: ", description, index));
    for (double[] doubles : matrix) {
      value.append(format("\n\t%s", Arrays.toString(doubles)));
    }
    value.append("\n");
    return value.toString();
  }

  /**
   * Produce a translation vector necessary to move an object with the current center of mass (com)
   * into the provided crystal.
   *
   * @param crystal Replicates crystal within whom coordinates should be moved.
   * @param com Center of mass (x, y, z) for the object of concern
   * @return double[] translation vector to move the object within the provided crystal.
   */
  private static double[] moveIntoCrystal(Crystal crystal, double[] com) {

    double[] translate = new double[3];
    double[] currentCoM = new double[3];
    currentCoM[0] = com[0];
    currentCoM[1] = com[1];
    currentCoM[2] = com[2];

    // Move the COM to the Replicates Crystal.
    crystal.toFractionalCoordinates(com, translate);
    translate[0] = mod(translate[0], 1.0);
    translate[1] = mod(translate[1], 1.0);
    translate[2] = mod(translate[2], 1.0);
    crystal.toCartesianCoordinates(translate, translate);

    // Correct center of mass.
    com[0] = translate[0];
    com[1] = translate[1];
    com[2] = translate[2];

    // The translation vector is difference between the new location and the current COM.
    translate[0] -= currentCoM[0];
    translate[1] -= currentCoM[1];
    translate[2] -= currentCoM[2];

    return translate;
  }

  /**
   * Determine the number of unique AUs within the replicates crystal to a tolerance.
   *
   * @param allCoords Coordinates for every atom in replicates crystal ([x1, y1, z1, x2, y2,
   *     z2...].
   * @param molIndices Prioritization of molecules.
   * @param auCoords Coordinates for single AU.
   * @param nCoords Number of coordinates in an AU (number of atoms * 3).
   * @param upperLimit The largest number of unique AUs (0 for no upper limit).
   * @param permute Search entire replicates crystal if true, otherwise only the expected.
   * @param nAUinReplicates Number of AUs in replicates crystal.
   * @param massN Array containing masses for each atom in AU.
   * @param matchTol Tolerance to determine whether two AUs are the same.
   */
  private void numberUniqueAUs(double[] allCoords, DoubleIndexPair[] molIndices, double[] auCoords,
      int nCoords, int upperLimit, boolean permute,
      int nAUinReplicates, double[] massN, double matchTol) {
    // uniqueDiffs is only recorded for logging... could remove later.
    // List of differences (RMSD_1) for AUs in replicates crystal.
    List<Double> uniqueDiffs = new ArrayList<>();
    tempListIndices.add(0);
    uniqueDiffs.add(rmsd(auCoords, auCoords, massN));
    tempListXYZ.add(auCoords);
    // Start from 1 as zero is automatically added.
    int index = 1;
    // Determine number of conformations in first crystal
    int numConfCheck = (permute || upperLimit <= 0) ? nAUinReplicates : upperLimit;
    if (logger.isLoggable(Level.FINEST)) {
      logger.finest("RMSD Differences in Replicates Crystal:");
    }
    double[] target = new double[nCoords];
    arraycopy(tempListXYZ.get(0), 0, target, 0, nCoords);
    translate(target, massN);
    while (uniqueDiffs.size() < numConfCheck) {
      double[] baseCheckMol = new double[nCoords];
      arraycopy(allCoords, molIndices[index].getIndex() * nCoords, baseCheckMol, 0,
          nCoords);
      translate(baseCheckMol, massN);
      rotate(target, baseCheckMol, massN);
      double value = rmsd(target, baseCheckMol, massN);
      if (logger.isLoggable(Level.FINEST)) {
        logger.finest(format("%4d %4.4f", molIndices[index].getIndex(), value));
      }
      if (addLooseUnequal(uniqueDiffs, value, matchTol)) {
        tempListIndices.add(index);
        tempListXYZ.add(baseCheckMol);
      }
      index++;
      if (index >= molIndices.length) {
        break;
      }
    }
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(" RMSD_1 from 1st AU:\n i RMSD     AU Index");
      int numUniqueIndices = tempListIndices.size();
      for (int i = 0; i < numUniqueIndices; i++) {
        logger.finer(format(" %d %4.4f    %4d", i, uniqueDiffs.get(i),
            molIndices[tempListIndices.get(i)].getIndex()));
      }
    }
  }

  /**
   * Pair species between two crystals based on distances between centers.
   *
   * @param desiredAUs Number of pairs to determine.
   */
  private void pairEntities(int desiredAUs) {
    tempListIndices.clear();
    // List of indexes for second system.
    for (DoubleIndexPair doubleIndexPair : molDist2_2) {
      // Only search molecules within range of the desired number of molecules.
      // Must have enough molecules for matching (using exhaustive till better heuristic is determined)
      tempListIndices.add(doubleIndexPair.getIndex());
    }

    // Compare distances between center of masses from system 1 and 2.
    for (int i = 0; i < desiredAUs; i++) {
      double minDist = Double.MAX_VALUE;
      Integer minIndex = -1;
      double[] baseCoMCurr = baseCoM[molDist1_2[i].getIndex()];
      for (Integer target : tempListIndices) {
        double dist = dist2(baseCoMCurr, targetCoM[target]);
        if (dist < minDist) {
          minIndex = target;
          minDist = dist;
        }
        if (abs(minDist) < MATCH_TOLERANCE) {
          // Distance between center of masses is ~0 is the best scenario assuming no coordinate overlaps.
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(" \tExit out of match loop.\n");
          }
          break;
        }
      }
      pairedAUs[i] = new DoubleIndexPair(minIndex, minDist);
      if (!tempListIndices.remove(minIndex)) {
        logger.warning(format(" Index value of %d was not found (%4.4f).", minIndex, minDist));
      }
      if (logger.isLoggable(Level.FINEST)) {
        stringBuilder.append(
            format("\n Base position:   %d: %8.4f %8.4f %8.4f\n", i,
                baseCoM[molDist1_2[i].getIndex()][0],
                baseCoM[molDist1_2[i].getIndex()][1],
                baseCoM[molDist1_2[i].getIndex()][2]));
        stringBuilder.append(
            format(" Match Distance:  %d: %8.4f\n", i, sqrt(pairedAUs[i].getDoubleValue())));
        stringBuilder.append(
            format(" Target position: %d: %8.4f %8.4f %8.4f\n", i,
                targetCoM[pairedAUs[i].getIndex()][0],
                targetCoM[pairedAUs[i].getIndex()][1],
                targetCoM[pairedAUs[i].getIndex()][2]));
      }
    }
    if (logger.isLoggable(Level.FINEST)) {
      stringBuilder.append("  Distance between pairs after rot 2:\n" +
          " Index  Base Index  Target Index    Match Index    Distance\n");
      for (int i = 0; i < desiredAUs; i++) {
        stringBuilder.append(format(" %5d %10d %14d %14d %10.4f\n", i,
            molDist1_2[i].getIndex(), molDist2_2[i].getIndex(), pairedAUs[i].getIndex(),
            sqrt(pairedAUs[i].getDoubleValue())));
      }
    }
  }

  /**
   * Print values for the symmetry operations performed so far (useful for debugging printSym flag).
   *
   * @param compareAtomsSize Number of atoms being compared from each AU.
   * @param save Save integer switch to determine if/how to save file.
   */
  private void printSym(int compareAtomsSize, File file, String name, List<Bond> bondList, Atom[] atoms,
                        ForceField forceField, int save) {
    // Apply inverse of base operations:
    // For orthogonal matrices the inverse matrix = the transpose. True iff det(A)== +/- 1.
    // No inverse if det(A)==0.
    double[][] tempTransform = new double[4][4];
    copyArrayValues(tempTransform, targetTransformSymOp.asMatrix());
    addTranslation(targetSymOp.tr, tempTransform, true);
    addRotation(targetSymOp.rot, tempTransform, true);
    double[] bestTranslation = new double[] {
        tempTransform[0][3] / tempTransform[3][3],
        tempTransform[1][3] / tempTransform[3][3],
        tempTransform[2][3] / tempTransform[3][3]};
    SymOp symOp = new SymOp(Arrays.copyOf(tempTransform, 3), bestTranslation);
    double[] newMol = new double[compareAtomsSize * 3];
    stringBuilder.append(
        "\n Print Current Symmetry Operator:\n \tOriginal Coords \t\t Current Coords \t ==\t Applied Sym Op Coords");
    for (int i = 0; i < compareAtomsSize; i++) {
      int k = i * 3;
      double[] xyz = new double[] {targetAUoriginal[k], targetAUoriginal[k + 1],
          targetAUoriginal[k + 2]};
      applyCartesianSymOp(xyz, xyz, symOp);
      newMol[k] = xyz[0];
      newMol[k + 1] = xyz[1];
      newMol[k + 2] = xyz[2];
      stringBuilder.append(format("\n %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f",
              targetAUoriginal[k], targetAUoriginal[k + 1], targetAUoriginal[k + 2],
              targetAU[k], targetAU[k + 1], targetAU[k + 2],
              newMol[k], newMol[k + 1], newMol[k + 2]));
    }
    saveAssembly(file, name, bondList, atoms, forceField, newMol, comparisonAtoms, "_moveMethod", 0, save);
    stringBuilder.append("\n").append(symOp).append("\n");
  }

  /**
   * Prioritize asymmetric units within the system based on distance to specified index.
   *
   * @param coordsXYZ Coordinates for expanded crystal (should contain 3 * nAtoms * nMols
   *     entries).
   * @param mass Mass of atoms within asymmetric unit (should contain one mass per atom in asym
   *     unit).
   * @param massSum Sum of atomic masses within asymmetric unit.
   * @param centerOfMasses Center of masses for each replicate within inflated crystal.
   * @param nAtoms Number of coordinates in an entity.
   * @param molDists Prioritization of AUs in expanded system based on linkage criteria
   * @param index Index of molecules to be center.
   * @param linkage User specified criteria to determine prioritization.
   */
  private static void prioritizeReplicates(double[] coordsXYZ, double[] mass,
      double massSum, double[][] centerOfMasses, int nAtoms,
      DoubleIndexPair[] molDists, int index, int linkage) {
    // Find AU to be treated as the new center.
    // AUs added to system based on distance to center of all atoms. Index = 0 AU should be closest to all atom center.
    int nCoords = nAtoms * 3;
    int nMols = coordsXYZ.length / nCoords;
    if (linkage == 0) {
      // Prioritize based on closest atomic distances.
      int centerIndex = index * nCoords;
      for (int i = 0; i < nMols; i++) {
        double tempDist = Double.MAX_VALUE;
        int molIndex = i * nCoords;
        for (int j = 0; j < nAtoms; j++) {
          int centerAtomIndex = j * 3;
          double[] centerXYZ = {coordsXYZ[centerIndex + centerAtomIndex],
              coordsXYZ[centerIndex + centerAtomIndex + 1],
              coordsXYZ[centerIndex + centerAtomIndex + 2]};
          for (int k = 0; k < nAtoms; k++) {
            int atomIndex = k * 3;
            double[] xyz = {coordsXYZ[molIndex + atomIndex],
                coordsXYZ[molIndex + atomIndex + 1],
                coordsXYZ[molIndex + atomIndex + 2]};
            double currDist = dist2(centerXYZ, xyz);
            if (currDist < tempDist) {
              tempDist = currDist;
            }
            if (abs(tempDist) < MATCH_TOLERANCE) {
              break;
            }
          }
          if (abs(tempDist) < MATCH_TOLERANCE) {
            break;
          }
        }
        molDists[i] = new DoubleIndexPair(i, tempDist);
      }
      // Sort so the smallest distance is at position 0.
      sort(molDists);

      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists.length);
        double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        logger.finer(" Prioritize replicates:");
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
              format("  AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i, index,
                  molDists[i].getIndex(),
                  sqrt(molDists[i].getDoubleValue()), targetMol[molDists[i].getIndex()][0],
                  targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
        }
        logger.finest("\n");
      }
    } else if (linkage == 1) {
      // Prioritize based on geometric center/center of mass
      double[] coordCenter = centerOfMasses[index];
      for (int i = 0; i < nMols; i++) {
        double[] moleculeCenter = centerOfMasses[i];
        molDists[i] = new DoubleIndexPair(i, dist2(coordCenter, moleculeCenter));
      }
      // Reorder based on distance to AU closest to Index.
      sort(molDists);

      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists.length);
        double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
              format("  1AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i, index,
                  molDists[i].getIndex(),
                  sqrt(molDists[i].getDoubleValue()), targetMol[molDists[i].getIndex()][0],
                  targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
        }
      }

      // Molecules in crystal sorted based on distance to center of mass of center most molecule
      // Want the first two molecules chosen in this manner, but third molecule to be closest to both
      // Assigning distance from molDists ensures correct ordering of center most molecule.
      DoubleIndexPair[] molDists2 = new DoubleIndexPair[nMols];
      molDists2[0] = new DoubleIndexPair(molDists[0].getIndex(), molDists[0].getDoubleValue());

      double[] avgCenter = new double[3];
      avgCenter[0] =
          (centerOfMasses[molDists[0].getIndex()][0] + centerOfMasses[molDists[1].getIndex()][0])
              / 2;
      avgCenter[1] =
          (centerOfMasses[molDists[0].getIndex()][1] + centerOfMasses[molDists[1].getIndex()][1])
              / 2;
      avgCenter[2] =
          (centerOfMasses[molDists[0].getIndex()][2] + centerOfMasses[molDists[1].getIndex()][2])
              / 2;

      for (int i = 1; i < nMols; i++) {
        double[] moleculeCenter = centerOfMasses[molDists[i].getIndex()];
        molDists2[i] = new DoubleIndexPair(molDists[i].getIndex(), dist2(avgCenter, moleculeCenter));
      }
      sort(molDists2);
      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists2.length);
        double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
              format("  1,2AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i,
                  index,
                  molDists2[i].getIndex(),
                  sqrt(molDists2[i].getDoubleValue()), targetMol[molDists2[i].getIndex()][0],
                  targetMol[molDists2[i].getIndex()][1], targetMol[molDists2[i].getIndex()][2]));
        }
      }
      // Molecules in crystal sorted based on distance to center of mass of center 2 most molecule
      // Want the first three molecules chosen in this manner, but rest to be closest to all three
      // Assigning distance from molDists2 ensures correct ordering of center most molecule.
      DoubleIndexPair[] molDists3 = new DoubleIndexPair[nMols];
      molDists3[0] = new DoubleIndexPair(molDists2[0].getIndex(), molDists2[0].getDoubleValue());
      molDists3[1] = new DoubleIndexPair(molDists2[1].getIndex(), molDists2[1].getDoubleValue());
      avgCenter = new double[3];
      for (int i = 0; i < 3; i++) {
        avgCenter[0] += centerOfMasses[molDists2[i].getIndex()][0];
        avgCenter[1] += centerOfMasses[molDists2[i].getIndex()][1];
        avgCenter[2] += centerOfMasses[molDists2[i].getIndex()][2];
      }
      for (int i = 0; i < 3; i++) {
        avgCenter[i] /= 3;
      }

      for (int i = 2; i < nMols; i++) {
        double[] moleculeCenter = centerOfMasses[molDists2[i].getIndex()];
        molDists3[i] = new DoubleIndexPair(molDists2[i].getIndex(),
            dist2(avgCenter, moleculeCenter));
      }
      //Reorder based on center point between center-most AU to all atom center and closest AU to center-most AU.
      sort(molDists3);
      arraycopy(molDists3, 0, molDists, 0, nMols);
      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists3.length);
        double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
              format("  1,2,3AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i,
                  index,
                  molDists3[i].getIndex(),
                  sqrt(molDists3[i].getDoubleValue()), targetMol[molDists3[i].getIndex()][0],
                  targetMol[molDists3[i].getIndex()][1], targetMol[molDists3[i].getIndex()][2]));
        }
      }
    } else if (linkage == 2) {
      // Prioritize based on minimum distance between the farthest atom.
      int centerIndex = index * nCoords;
      for (int i = 0; i < nMols; i++) {
        double tempDist = 0.0;
        int molIndex = i * nCoords;
        for (int j = 0; j < nAtoms; j++) {
          int centerAtomIndex = j * 3;
          double[] centerXYZ = {coordsXYZ[centerIndex + centerAtomIndex],
              coordsXYZ[centerIndex + centerAtomIndex + 1],
              coordsXYZ[centerIndex + centerAtomIndex + 2]};
          for (int k = 0; k < nAtoms; k++) {
            int atomIndex = k * 3;
            double[] xyz = {coordsXYZ[molIndex + atomIndex],
                coordsXYZ[molIndex + atomIndex + 1],
                coordsXYZ[molIndex + atomIndex + 2]};
            double currDist = dist2(centerXYZ, xyz);
            if (currDist > tempDist) {
              tempDist = currDist;
            }
            if (abs(tempDist) < MATCH_TOLERANCE) {
              break;
            }
          }
          if (abs(tempDist) < MATCH_TOLERANCE) {
            break;
          }
        }
        molDists[i] = new DoubleIndexPair(i, tempDist);
      }
      // Sort so the smallest distance is at position 0.
      sort(molDists);

      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists.length);
        double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
              format("  AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i, index,
                  molDists[i].getIndex(),
                  sqrt(molDists[i].getDoubleValue()), targetMol[molDists[i].getIndex()][0],
                  targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
        }
      }
    } else {
      logger.severe(format(" Linkage value of %d was unrecognized. Try 0 (, 1, or 2.", linkage));
    }
  }

  /**
   * Read in the distance matrix.
   *
   * @param filename The PAC RMSD matrix file to read from.
   * @return Stats for all read in distance matrix values.
   */
  private RunningStatistics readMatrix(String filename) {
    restartRow = 0;
    restartColumn = 0;

    DistanceMatrixFilter distanceMatrixFilter = new DistanceMatrixFilter();
    RunningStatistics runningStatistics = distanceMatrixFilter.readDistanceMatrix(
        filename, baseSize, targetSize);

    if (runningStatistics != null && runningStatistics.getCount() > 0) {
      restartRow = distanceMatrixFilter.getRestartRow();
      restartColumn = distanceMatrixFilter.getRestartColumn();

      if (isSymmetric) {
        // Only the diagonal entry (0.0) is on the last row for a symmetric matrix.
        if (restartRow == baseSize && restartColumn == 1) {
          logger.info(format(" Complete symmetric distance matrix found (%d x %d).", restartRow,
              restartRow));
        } else {
          restartColumn = 0;
          logger.info(format(
              " Incomplete symmetric distance matrix found.\n Restarting at row %d, column %d.",
              restartRow + 1, restartColumn + 1));
        }
      } else if (restartRow == baseSize && restartColumn == targetSize) {
        logger.info(format(" Complete distance matrix found (%d x %d).", restartRow, restartColumn));
      } else {
        restartColumn = 0;
        logger.info(format(" Incomplete distance matrix found.\n Restarting at row %d, column %d.",
            restartRow + 1, restartColumn + 1));
      }
    }

    return runningStatistics;
  }

  /**
   * Reduce asymmetric unit to atoms that are going to be used in final RMSD.
   *
   * @param atoms Atoms we wish to reduce.
   * @param comparisonAtoms Atoms of interest within asymmetric unit.
   * @return Linear coordinates for only atoms of interest.
   */
  private static double[] reduceSystem(Atom[] atoms, int[] comparisonAtoms) {
    // Collect asymmetric unit atomic coordinates.
    double[] reducedCoords = new double[comparisonAtoms.length * 3];
    int coordIndex = 0;
    for (Integer value : comparisonAtoms) {
      Atom atom = atoms[value];
      reducedCoords[coordIndex++] = atom.getX();
      reducedCoords[coordIndex++] = atom.getY();
      reducedCoords[coordIndex++] = atom.getZ();
    }
    return reducedCoords;
  }

  /**
   * Save the provided coordinates as a PDB file.
   *
   * @param coords Coordinates to be saved within the PDB.
   * @param comparisonAtoms Atoms of interest within the initial asymmetric unit.
   * @param description Unique identifier that will be added to PDB file name.
   * @param compNum Unique number for the current comparison
   * @param save Type of file to save (0=none, 1=PDB, 2=XYZ)
   */
  private void saveAssembly(File file, String name, List<Bond> bondList, Atom[] atoms, ForceField forceField0, final double[] coords,
      final int[] comparisonAtoms, String description, int compNum, int save) {
    //TODO: Save systems out as original molecule regardless of selection
    String fileName = FilenameUtils.removeExtension(file.getName());
    File saveLocation;
    if (save == 2) {
      saveLocation = new File(fileName + description + ".xyz");
    } else if (save == 1) {
      saveLocation = new File(fileName + description + ".pdb");
    } else {
      return;
    }
    // Save aperiodic system of n_mol closest atoms for visualization
    MolecularAssembly currentAssembly = new MolecularAssembly(name);
    ArrayList<Atom> newAtomList = new ArrayList<>();
    int atomIndex = 1;
    int compareAtomsSize = comparisonAtoms.length;
    int nCoords = compareAtomsSize * 3;
    // Reset atom Index for indexing new atoms
    int numMols = coords.length / (nCoords);
    for (int n = 0; n < numMols; n++) {
      // Obtain atoms from moved AU (create copy to move to origin)
      // move original and copy to origin
      // rotate original to match copy
      // translate back to moved location
      // Add atoms from moved original to atom list.
      ArrayList<Atom> atomList = new ArrayList<>();
      // Create a new set of Atoms for each SymOp molecule
      int atomValue = 0;
      //Add atoms from comparison to output assembly.
      for (Integer i : comparisonAtoms) {
        Atom a = atoms[i];
        double[] xyz = new double[3];
        xyz[0] = coords[n * nCoords + atomValue];
        xyz[1] = coords[n * nCoords + atomValue + 1];
        xyz[2] = coords[n * nCoords + atomValue + 2];
        Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
        atomList.add(atom);
        atomValue += 3;
      }
      // Setup bonds for AUs in comparison.
      if (save == 2) {
        //TODO make more robust. Currently only handles all atom selection. Not subsets.
        for (Bond bond : bondList) {
          Atom a1 = bond.getAtom(0);
          Atom a2 = bond.getAtom(1);
          //Indices stored as human-readable.
          int a1Ind = a1.getIndex() - 1;
          int a2Ind = a2.getIndex() - 1;
          if (IntStream.of(comparisonAtoms).anyMatch(x -> x == a1Ind) && IntStream.of(
              comparisonAtoms).anyMatch(x -> x == a2Ind)) {
            Atom newA1 = atomList.get(a1Ind);
            Atom newA2 = atomList.get(a2Ind);
            Bond b = new Bond(newA1, newA2);
            b.setBondType(bond.getBondType());
          }
        }
      }
      newAtomList.addAll(atomList);
    }

    // Construct the force field for the expanded set of molecules
    ForceField forceField = new ForceField(forceField0.getProperties());

    // Clear all periodic boundary keywords to achieve an aperiodic system.
    forceField.clearProperty("a-axis");
    forceField.clearProperty("b-axis");
    forceField.clearProperty("c-axis");
    forceField.clearProperty("alpha");
    forceField.clearProperty("beta");
    forceField.clearProperty("gamma");
    forceField.clearProperty("spacegroup");
    currentAssembly.setForceField(forceField);

    // The biochemistry method is designed to load chemical entities into the
    // Polymer, Molecule, Water and Ion data structure.
    Utilities.biochemistry(currentAssembly, newAtomList);

    currentAssembly.setFile(file);
    if (save == 2) {
      File key = new File(fileName + ".key");
      File properties = new File(fileName + ".properties");
      if (key.exists()) {
        File keyComparison = new File(fileName + description + ".key");
        try {
          if (keyComparison.createNewFile()) {
            FileUtils.copyFile(key, keyComparison);
          } else {
            stringBuilder.append("\n Could not create properties file.\n");
          }
        } catch (Exception ex) {
          // Likely using properties file.
          if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(ex);
          }
        }
      } else if (properties.exists()) {
        File propertiesComparison = new File(fileName + description + ".properties");
        try {
          if (propertiesComparison.createNewFile()) {
            FileUtils.copyFile(properties, propertiesComparison);
          } else {
            stringBuilder.append("\n Could not create properties file.\n");
          }
        } catch (Exception ex) {
          // Likely not using a key/properties file... so PDB?
          stringBuilder.append(
              "\n Neither key nor properties file detected therefore only creating XYZ.\n");
          if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(ex);
          }
        }
      }
      XYZFilter xyzFilter = new XYZFilter(saveLocation, currentAssembly, forceField,
          currentAssembly.getProperties());
      xyzFilter.writeFile(saveLocation, true);
    } else {
      PDBFilter pdbfilter = new PDBFilter(saveLocation, currentAssembly, forceField,
          currentAssembly.getProperties());
      pdbfilter.setModelNumbering(compNum);
      pdbfilter.writeFile(saveLocation, true, false, false);
    }
    currentAssembly.destroy();
  }

  /**
   * Save the provided coordinates as a PDB file with accompanying CSV containing RMSD.
   *
   * @param coords Coordinates to be saved within the PDB.
   * @param comparisonAtoms Atoms of interest within the initial asymmetric unit.
   * @param description Unique identifier that will be added to PDB file name.
   * @param finalRMSD RMSD to be saved to CSV file.
   * @param compNum Unique number for the current comparison
   * @param save Type of file to save (0=none, 1=PDB, 2=XYZ)
   */
  private void saveAssembly(File file, String name,List<Bond> bondList, Atom[] atoms, ForceField forceField, double[] coords,
      int[] comparisonAtoms, String description, double finalRMSD, int compNum, int save) {
    saveAssembly(file, name, bondList, atoms, forceField, coords, comparisonAtoms, description, compNum, save);
    String fileName = FilenameUtils.removeExtension(file.getName());
    try {
      // Needs same name as PDB so follow PDB format.
      File csv = PDBFilter.version(new File(fileName + description + ".csv"));
      if (csv.createNewFile()) {
        BufferedWriter bw = new BufferedWriter(new FileWriter(csv, false));
        bw.append("rms\n");
        bw.append(format("%4.4f\n", finalRMSD));
        bw.close();
      } else {
        logger.warning(" Could not create CSV file \"" + csv.getName() + "\"");
      }
    } catch (Exception ex) {
      logger.info(ex.toString());
    }
  }
}

//        // Used in QEtoXYZ.groovy which is not ready for git which is why this method appears unused.
//    /**
//     * Orient coordinates so that the second index is on the x axis, and the third index is on the X-Y
//     * plane. First index should be at the origin (0, 0, 0).
//     *
//     * @param coordsXYZ   An array of XYZ positions (e.g. [x0, y0, z0, x1, y1, z1, x2, y2, z2]
//     * @param atomIndices Indices for three desired sets from the XYZ list (e.g. [0, 1, 2]). Index
//     *                    0 should be at origin.
//     */
//    public static void standardOrientation(double[] coordsXYZ, int[] atomIndices) {
//        int numCoords = coordsXYZ.length / 3;
//        double[] atomCoords = new double[3 * 3];
//        atomCoords[0] = coordsXYZ[atomIndices[0]];
//        atomCoords[1] = coordsXYZ[atomIndices[0] + 1];
//        atomCoords[2] = coordsXYZ[atomIndices[0] + 2];
//        atomCoords[3] = coordsXYZ[atomIndices[1]];
//        atomCoords[4] = coordsXYZ[atomIndices[1] + 1];
//        atomCoords[5] = coordsXYZ[atomIndices[1] + 2];
//        atomCoords[6] = coordsXYZ[atomIndices[2]];
//        atomCoords[7] = coordsXYZ[atomIndices[2] + 1];
//        atomCoords[8] = coordsXYZ[atomIndices[2] + 2];
//
//        // TODO: Delete coordsXYZOrig?
//        double[] coordsXYZOrig = new double[numCoords * 3];
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            coordsXYZOrig[atomIndex] = coordsXYZ[atomIndex];
//            coordsXYZOrig[atomIndex + 1] = coordsXYZ[atomIndex + 1];
//            coordsXYZOrig[atomIndex + 2] = coordsXYZ[atomIndex + 2];
//        }
//
//        // TODO: Delete atomsCoordsOrig?
//        double[] atomsCoordsOrig = new double[3 * 3];
//        arraycopy(atomCoords, 0, atomsCoordsOrig, 0, 9);
//
//        logger.fine(
//                format(" START: N1:\t%16.15f %16.15f %16.15f", atomCoords[0], atomCoords[1], atomCoords[2]));
//        logger.fine(
//                format(" START: N2:\t%16.15f %16.15f %16.15f", atomCoords[3], atomCoords[4], atomCoords[5]));
//        logger.fine(
//                format(" START: N3:\t%16.15f %16.15f %16.15f", atomCoords[6], atomCoords[7], atomCoords[8]));
//
//        double p1n2 = coordsXYZ[atomIndices[1]];
//        double q1n2 = coordsXYZ[atomIndices[1] + 1];
//        double r1n2 = coordsXYZ[atomIndices[1] + 2];
//
//        // Calculation of sigma, phai, and cita angles needed to get specified atoms to desired loci
//        double cita0 = acos(p1n2 / sqrt(p1n2 * p1n2 + q1n2 * q1n2));
//        double phai0 = acos(sqrt(p1n2 * p1n2 + q1n2 * q1n2) /
//                sqrt(p1n2 * p1n2 + q1n2 * q1n2 + r1n2 * r1n2));
//        if (q1n2 < 0.0) {
//            cita0 = -cita0;
//        }
//
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            double ptmp = coordsXYZ[atomIndex] * cos(cita0) + coordsXYZ[atomIndex + 1] * sin(cita0);
//            double qtmp = -coordsXYZ[atomIndex] * sin(cita0) + coordsXYZ[atomIndex + 1] * cos(cita0);
//            coordsXYZ[atomIndex] = ptmp;
//            coordsXYZ[atomIndex + 1] = qtmp;
//        }
//
//        p1n2 = coordsXYZ[atomIndices[1]];
//        q1n2 = coordsXYZ[atomIndices[1] + 1];
//
//        if (r1n2 > 0.0) {
//            phai0 = -phai0;
//        }
//
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            double ptmp = coordsXYZ[atomIndex] * cos(phai0) - coordsXYZ[atomIndex + 2] * sin(phai0);
//            double rtmp = coordsXYZ[atomIndex] * sin(phai0) + coordsXYZ[atomIndex + 2] * cos(phai0);
//            coordsXYZ[atomIndex] = ptmp;
//            coordsXYZ[atomIndex + 2] = rtmp;
//        }
//
//        p1n2 = coordsXYZ[atomIndices[1]];
//        r1n2 = coordsXYZ[atomIndices[1] + 2];
//        double p1n3 = coordsXYZ[atomIndices[2]];
//        double q1n3 = coordsXYZ[atomIndices[2] + 1];
//        double r1n3 = coordsXYZ[atomIndices[2] + 2];
//
//        double sigma0 = acos(q1n3 / sqrt(q1n3 * q1n3 + r1n3 * r1n3));
//        if (r1n3 < 0.0) {
//            sigma0 = -sigma0;
//        }
//
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            double qtmp = coordsXYZ[atomIndex + 1] * cos(sigma0) + coordsXYZ[atomIndex + 2] * sin(sigma0);
//            double rtmp = -coordsXYZ[atomIndex + 1] * sin(sigma0) + coordsXYZ[atomIndex + 2] * cos(sigma0);
//            coordsXYZ[atomIndex + 1] = qtmp;
//            coordsXYZ[atomIndex + 2] = rtmp;
//        }
//
//        q1n2 = coordsXYZ[atomIndices[1] + 1];
//        r1n2 = coordsXYZ[atomIndices[1] + 2];
//        q1n3 = coordsXYZ[atomIndices[2] + 1];
//        r1n3 = coordsXYZ[atomIndices[2] + 2];
//
//        if (logger.isLoggable(Level.FINER)) {
//            logger.finer(
//                    format(" DONE N1: %16.15f %16.15f %16.15f", atomCoords[0], atomCoords[1], atomCoords[2]));
//            logger.finer(format(" DONE N2: %16.15f %16.15f %16.15f", p1n2, q1n2, r1n2));
//            logger.finer(format(" DONE N3: %16.15f %16.15f %16.15f", p1n3, q1n3, r1n3));
//        }
//    }

//  Frank-Kasper phases of metallic ions can reach a coordination number of 16...