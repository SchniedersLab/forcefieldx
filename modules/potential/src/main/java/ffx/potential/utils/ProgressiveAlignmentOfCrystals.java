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

import static ffx.crystal.SymOp.Tr_0_0_0;
import static ffx.crystal.SymOp.ZERO_ROTATION;
import static ffx.crystal.SymOp.applyCartesianSymOp;
import static ffx.crystal.SymOp.invertSymOp;
import static ffx.numerics.math.ScalarMath.mod;
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
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static java.util.Arrays.sort;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.math3.util.FastMath.*;

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
import ffx.potential.cli.TopologyOptions;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.DistanceMatrixFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.IndexIndexPair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
  private static final Logger logger = Logger.getLogger(ProgressiveAlignmentOfCrystals.class.getName());
  /**
   * SystemFilter containing structures for first file.
   */
  private final SystemFilter baseFilter;
  /**
   * Number of structures stored in SystemFilter of first file.
   */
  private final int baseSize;
  /**
   * Label for the first file.
   */
  private final String baseLabel;
  /**
   * SystemFilter containing structures for second file.
   */
  private final SystemFilter targetFilter;
  /**
   * Number of structures stored in SystemFilter for second file.
   */
  private final int targetSize;
  /**
   * Label for the second file.
   */
  private final String targetLabel;
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
   * Masses for N asymmetric units (nAU * nCoords).
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
   * List of each AU index's distance to central AU for first crystal (Maps baseXYZ).
   */
  private DoubleIndexPair[] baseAUDist;
  /**
   * Working copy of baseAUDist, updated when treating a new AU as center.
   */
  private DoubleIndexPair[] baseAUDist_2;
  /**
   * List of each AU index's distance to central AU for second crystal (Maps targetXYZ).
   */
  private DoubleIndexPair[] targetAUDist;
  /**
   * Working copy of targetAUDist, updated when treating a new AU as center.
   */
  private DoubleIndexPair[] targetAUDist_2;
  /**
   * List containing indices for each unique AU in current crystal.
   */
  private final List<Integer> tempListIndices = new ArrayList<>();
  /**
   * Array containing indices for each unique AU in first crystal (Maps baseAUDist).
   */
  private Integer[] uniqueBase;
  /**
   * Array containing indices for each unique AU in second crystal (Maps targetAUDist).
   */
  private Integer[] uniqueTarget;
  /**
   * Array containing XYZ coordinates for first crystal sorted based on distance to replicates center.
   */
  private double[] baseXYZoriginal;
  /**
   * Array containing XYZ coordinates for first crystal sorted based on distance to replicates center.
   * Working copy of baseXYZoriginal that undergoes rotations/translations during comparisons.
   */
  private double[] baseXYZ;
  /**
   * Array containing XYZ coordinates for second crystal sorted based on distance to replicates center.
   */
  private double[] targetXYZoriginal;
  /**
   * Array containing XYZ coordinates for second crystal sorted based on distance to replicates center.
   * Working copy of baseXYZoriginal that undergoes rotations/translations during comparisons.
   */
  private double[] targetXYZ;
  /**
   * Center of masses (or geometric center) for every AU in first replicates crystal (comparison, AU num, XYZ).
   */
  private double[][][] baseCoM;
  /**
   * Original coordinates for central AU of second crystal (first index Z', second index X, Y, Z coords).
   * These coordinates should match coordinates in input file.
   */
  private double[][] baseAUoriginal;
  /**
   * Coordinates for central AU of first crystal (X, Y, Z coords).
   */
  private double[] baseAU;
  /**
   * Coordinates for central 3 AUs of first crystal  ( X, Y, Z coords).
   */
  private double[] base3AUs;
  /**
   * Coordinates for central N AUs of first crystal  (comparison value, second index X, Y, Z coords).
   */
  private double[][] baseNAUs;
  /**
   * Coordinates of N AUs from first crystal of the closest match (lowest RMSD).
   */
  private double[][][] bestBaseNAUs;
  /**
   * Center of masses (or geometric center) for every AU in second replicates crystal
   */
  private double[][] targetCoM;
  /**
   * Original coordinates for central AU of second crystal (first index Z', second index X, Y, Z coords).
   * These coordinates should match coordinates in input file.
   */
  private double[][] targetAUoriginal;
  /**
   * Coordinates for central AU of second crystal (first index Z', second index X, Y, Z coords).
   */
  private double[] targetAU;
  /**
   * Coordinates for central 3 AUs of second crystal (first index Z', second index X, Y, Z coords).
   */
  private double[] target3AUs;
  /**
   * Coordinates for central N AUs of second crystal (first index Z', second index X, Y, Z coords).
   */
  private double[] targetNAUs;
  /**
   * Coordinates of N AUs from second crystal of the closest match (lowest RMSD) (first index Z', second index X, Y, Z coords).
   */
  private double[][][] bestTargetNAUs;
  /**
   * Indices and distances of matching AUs between first and second crystal.
   */
  private DoubleIndexPair[] pairedAUs;
  /**
   * Default PAC caches files to promote efficiency (otherwise use low memory flag).
   */
  private File[] fileCache;
  /**
   * Default PAC caches file names to promote efficiency (otherwise use low memory flag).
   */
  private String[] nameCache;
  /**
   * Default PAC caches bonds to promote efficiency (otherwise use low memory flag).
   */
  private ArrayList<Bond>[] bondCache;
  /**
   * Default PAC caches atoms to promote efficiency (otherwise use low memory flag).
   */
  private Atom[][] atomCache;
  /**
   * Default PAC caches force field properties to promote efficiency (otherwise use low memory flag).
   */
  private ForceField[] forceFieldCache;
  /**
   * Default PAC caches crystals to promote efficiency (otherwise use low memory flag).
   */
  private Crystal[] crystalCache;
  /**
   * Radius of gyration for the best matching clusters (position 0 for first and 1 for second
   * crystal).
   */
  private final double[] gyrations = new double[2];
  /**
   * Moments of inertia and vectors of application for first crystal.
   */
  private double[][] bestBaseMandV = new double[3][4];
  /**
   * Moments of inertia and vectors of application for second crystal.
   */
  private double[][] bestTargetMandV = new double[3][4];
  /**
   * Radius of gyration components for first crystal.
   */
  private double[] bestBaseRg = new double[3];
  /**
   * Radius of gyration components for second crystal.
   */
  private double[] bestTargetRg = new double[3];
  /**
   * Working copy of replicates Sym Op applied to asymmetric unit of first system.
   */
  private SymOp baseSymOp;
  /**
   * Replicates Sym Op applied to asymmetric unit of first system for best comparison.
   */
  private SymOp[][] bestBaseSymOp;
  /**
   * Working copy of replicates Sym Op applied to asymmetric unit of second system.
   */
  private SymOp targetSymOp;
  /**
   * Replicates Sym Op applied to asymmetric unit of second system for best comparison.
   */
  private SymOp[][] bestTargetSymOp;
  /**
   * List of symmetry operators used to create replicates crystal for 1st structure.
   */
  private ArrayList<SymOp> baseSymOps = new ArrayList<>();
  /**
   * Index order closest to replicates (get(0) closest to rep center; indexOf(i) ith Sym Op).
   */
  private ArrayList<Integer> baseDistMap = new ArrayList<>();
  /**
   * Index order closest to replicates (get(0) closest to rep center; indexOf(i) ith Sym Op).
   */
  private ArrayList<SymOp> targetSymOps = new ArrayList<>();
  /**
   * List of indices created through expansion of replicates crystal for 2nd structures.
   */
  private ArrayList<Integer> targetDistMap = new ArrayList<>();
  /**
   * Boolean to determine if comparison should print symmetry operators used.
   */
  private double printSym;
  /**
   * String builder to attempt to off load work from logger.
   */
  private final StringBuilder stringBuilder = new StringBuilder();
  /**
   * Working copy of symmetry operator used to generate central base molecule updated each iteration.
   */
  private SymOp baseTransformSymOp;
  /**
   * Working copy of symmetry operator used to generate central target molecules updated each iteration.
   */
  private SymOp targetTransformSymOp;
  /**
   * Best symmetry operator used to generate central base molecule.
   */
  private SymOp[][] bestBaseTransformSymOp;
  /**
   * Best symmetry operator used to generate central target molecules.
   */
  private SymOp[][] bestTargetTransformSymOp;
  /**
   * If molecules between two crystals differ below this tolerance, they are treated as equivalent.
   */
  private static final double MATCH_TOLERANCE = 1.0E-12;

  /**
   * Constructor for the ProgressiveAlignmentOfCrystals class.
   *
   * @param baseFilter   SystemFilter containing a set of crystal structures to compare.
   * @param targetFilter SystemFilter containing the other set of crystals to compare.
   * @param isSymmetric  Whether the comparison can be limited to the upper triangle.
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
   * @param file1              First file object containing crystals to compare
   * @param name1              Crystal structure (1st or base) that will remain relatively static (only translations).
   * @param bondList1          List of bonds in crystal used to save XYZ files.
   * @param atoms1             Atoms in first crystal used to save XYZ files.
   * @param forceField1        Force field values used to save files/print symmetry operators
   * @param file2              Second file object containing crystals to compare
   * @param name2              Crystal structure (2nd or target) that will rotate to match static assembly.
   * @param bondList2          List of bonds in crystal used to save XYZ files.
   * @param atoms2             Atoms in first crystal used to save XYZ files.
   * @param forceField2        Force field values used to save files/print symmetry operators
   * @param z1                 Number of molecules in asymmetric unit of first crystal.
   * @param z2                 Number of molecules in asymmetric unit of second crystal.
   * @param compareAtomsSize   Number of active atoms in asymmetric unit of first crystal.
   * @param nAU                Number of asymmetric units to compare.
   * @param baseSearchValue    Number of anticipated unique entities in 1st system.
   * @param targetSearchValue  Number of anticipated unique entities in 2nd system.
   * @param matchTol           Tolerance to determine whether two AUs are the same.
   * @param compNum            Comparison number based on all file submitted (logging).
   * @param strict            Compare all unique AUs between crystals.
   * @param save               Save out files of compared crystals.
   * @param machineLearning    Save out PDBs and CSVs of compared crystals.
   * @param linkage            Criteria to select nearest AUs (0=single, 1=average, 2=complete linkage).
   * @param inertia            Compute and display components of inertia tensor.
   * @param gyrationComponents Compute and display gyration components (asphericity, acylindricity, anisotropy, etc.)
   * @return the computed RMSD.
   */
  private double compare(final File file1, final String name1, final List<Bond> bondList1, final Atom[] atoms1,
                         final ForceField forceField1, final File file2, final String name2, final List<Bond> bondList2,
                         final Atom[] atoms2, final ForceField forceField2, final int z1, final int z2, final int compareAtomsSize,
                         final int nAU, final int baseSearchValue, final int targetSearchValue, final double matchTol, final int compNum,
                         final boolean strict, final int save, final boolean machineLearning, final int linkage, final boolean inertia,
                         final boolean gyrationComponents) {
    // TODO: Does PAC work for a combination of molecules and polymers?
    // It does not compare them on an individual basis, but can compare AUs as a whole (set zp/zp2 to 1).
    boolean useSave = save > 0;
    boolean useSym = printSym >= 0.0;
    int nCoords = compareAtomsSize * 3;
    // Number of species in expanded crystals.
    int nBaseMols = baseXYZ.length / nCoords;
    int nTargetMols = targetXYZ.length / nCoords;

    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(format("\n Base: %s Target: %s", name1, name2));
      stringBuilder.append(format("\n Comparing %3d of %5d in base sphere.\n" +
              " Comparing %3d of %5d in target sphere.\n", nAU, nBaseMols, nAU, nTargetMols));
    }

    if (useSym) {
      // Allocate variables used for Sym Op
      if (bestBaseSymOp == null || bestBaseSymOp.length != z1 || bestBaseSymOp[0].length != z2) {
        bestBaseSymOp = new SymOp[z1][z2];
      }
      if (bestTargetSymOp == null || bestTargetSymOp.length != z1 || bestTargetSymOp[0].length != z2) {
        bestTargetSymOp = new SymOp[z1][z2];
      }
      if (bestBaseTransformSymOp == null || bestBaseTransformSymOp.length != z1 || bestBaseTransformSymOp[0].length != z2) {
        bestBaseTransformSymOp = new SymOp[z1][z2];
      }
      if (bestTargetTransformSymOp == null || bestTargetTransformSymOp.length != z1 || bestTargetTransformSymOp[0].length != z2) {
        bestTargetTransformSymOp = new SymOp[z1][z2];
      }
    }

    // Reorder molDist2 as we shift a different molecule (m) to the center each loop.
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(" Prioritize target system.\n");
    }

    //Used for debugging. Can be removed.
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(" System 1 distances to center AU:\n");
      for (int i = 0; i < 10; i++) {
        stringBuilder.append(
                format(" %2d %2d\t%16.8f\n", i, baseAUDist[i].getIndex(), sqrt(baseAUDist[i].getDoubleValue())));
      }
      stringBuilder.append(" System 2 distances to center AU:\n");
      for (int i = 0; i < 10; i++) {
        stringBuilder.append(
                format(" %2d %2d\t%16.8f\n", i, targetAUDist[i].getIndex(), sqrt(targetAUDist[i].getDoubleValue())));
      }
    }

    //Determine if AUs in first system are same hand as center most in first (stereoisomer handling).
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(" Search conformations of each crystal:\n");
    }
    // NOTE not Z' representatives, just first Z closest AUs.
    arraycopy(baseXYZ, baseAUDist[0].getIndex() * nCoords, baseAU, 0, nCoords);

    tempListIndices.clear();
    // Compare replicates to initial molecule to determine unique conformations.
    numberUniqueAUs(baseXYZ, baseAUDist, baseAU.clone(), nCoords, baseSearchValue, strict,
            nBaseMols, massN, matchTol);
    int baseLength = tempListIndices.size();
    if (uniqueBase == null || uniqueBase.length != baseLength) {
      uniqueBase = new Integer[baseLength];
    }
    tempListIndices.toArray(uniqueBase);

    // NOTE not Z' representatives, just first Z closest AUs.
    arraycopy(targetXYZ, targetAUDist[0].getIndex() * nCoords, targetAU, 0, nCoords);
    tempListIndices.clear();
    // Compare replicates to initial molecule to determine unique conformations.
    numberUniqueAUs(targetXYZ, targetAUDist, targetAU.clone(), nCoords, targetSearchValue, strict,
            nTargetMols, massN, matchTol);
    int targetLength = tempListIndices.size();
    if (uniqueTarget == null || uniqueTarget.length != targetLength) {
      uniqueTarget = new Integer[targetLength];
    }
    tempListIndices.toArray(uniqueTarget);
    final boolean reverse = baseLength > targetLength;
    if (logger.isLoggable(Level.FINE)) {
      stringBuilder.append(format("\n  %d conformations detected out of %d in base crystal.\n " +
                      " %d conformations detected out of %d in target crystal.\n",
              baseLength, baseSearchValue, targetLength, targetSearchValue));
    }

    if (useSym) {
      // Obtain coordinates contained in file (assumes MoveIntoUnitCell was performed before PAC).
      if (z1 > 1 || z2 > 1) {
        stringBuilder.append("\n Utilizing Sym Ops with Z'>1. " +
                "Attempting to map all unique structures.\n");
      }
    }

    // Determine which unique AUs are most similar between crystals
    // Minimum difference between each unique target AU and closest matching base AU.
    double min1 = Double.MAX_VALUE;
    int minBIndex = -1;
    int minTIndex = -1;

    // Maintains index of minimum Z' differences for Sym Op prep
    IndexIndexPair[][] minZinds = new IndexIndexPair[z1][z2];
    for (IndexIndexPair[] row : minZinds) {
      fill(row, new IndexIndexPair(-1, -1));
    }

    // Maintains the minimum difference between Z' AUs
    double[][] minZdiffs = new double[z1][z2];
    for (double[] row : minZdiffs) {
      fill(row, Double.MAX_VALUE);
    }
    // Array to store all RMSD_1 comparisons
    double[][] minDiffs = new double[baseLength][targetLength];
    for (double[] row : minDiffs) {
      fill(row, Double.MAX_VALUE);
    }
    for (int i = 0; i < baseLength; i++) {
      int bIndex = uniqueBase[i];
      int baseInd = baseAUDist[bIndex].getIndex();
      int currZ1 = baseDistMap.get(baseInd) % z1;

      double[] baseXYZ = new double[nCoords];
      for (int j = 0; j < targetLength; j++) {
        // Need to remove translation/rotation from previous check.
        arraycopy(this.baseXYZ, baseInd * nCoords, baseXYZ, 0, nCoords);
        int tIndex = uniqueTarget[j];
        int targetInd = targetAUDist[tIndex].getIndex();
        int currZ2 = targetDistMap.get(targetInd) % z2;

        double[] targetXYZ = new double[nCoords];
        arraycopy(this.targetXYZ, targetInd * nCoords, targetXYZ, 0, nCoords);
        double[] tempTranB = calculateTranslation(baseXYZ, massN);
        applyTranslation(baseXYZ, tempTranB);
        double[] tempTranT = calculateTranslation(targetXYZ, massN);
        applyTranslation(targetXYZ, tempTranT);
        double[][] tempRotT = calculateRotation(baseXYZ, targetXYZ, massN);
        applyRotation(targetXYZ, tempRotT);
        double value = rmsd(baseXYZ, targetXYZ, massN);

        if (logger.isLoggable(Level.FINEST)) {
          stringBuilder.append(
                  format("\n Comp %3d: Base %2d IND %3d (Z'=%2d) \t Target %2d IND %3d (Z'=%2d) Diff: %8.4f\n",
                          j + i * targetLength, bIndex, baseInd, currZ1, tIndex, targetInd, currZ2, value));
          if (useSym) {
            stringBuilder.append(" Base Temp Translation: ").append(Arrays.toString(tempTranB)).append("\n");
            stringBuilder.append(" Target Temp Translation: ").append(Arrays.toString(tempTranT));
            stringBuilder.append(matrixToString(tempRotT, j + i * targetLength, "Temp Rot")).append("\n");
          }
        }

        minDiffs[i][j] = value;
        if (nAU == 1 && !strict) {
          // Record values for short circuit
          arraycopy(baseXYZ, 0, baseNAUs[0], 0, nAU * nCoords);
          arraycopy(targetXYZ, 0, targetNAUs, 0, nAU * nCoords);
          arraycopy(baseXYZ, 0, bestBaseNAUs[currZ1][currZ2], 0, nAU * nCoords);
          arraycopy(targetXYZ, 0, bestTargetNAUs[currZ1][currZ2], 0, nAU * nCoords);
        }
        if (useSym && value < minZdiffs[currZ1][currZ2]) {
          minZdiffs[currZ1][currZ2] = value;
          minZinds[currZ1][currZ2] = new IndexIndexPair(baseInd, targetInd);
          arraycopy(baseXYZ, 0, baseAU, 0, nCoords);
          arraycopy(targetXYZ, 0, targetAU, 0, nCoords);
          if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(format("\n Base AU: %2d of %2d" +
                            "\t\tvs\t\t Target AU: %2d of %2d \t (%9.4f)",
                    currZ1 + 1, z1, currZ2 + 1, z2, minZdiffs[currZ1][currZ2]));
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append("\n Base Translation: ").append(Arrays.toString(tempTranB)).append("\n");
              stringBuilder.append(" Target Translation: ").append(Arrays.toString(tempTranT));
              stringBuilder.append(matrixToString(tempRotT, currZ1, "Target Rotation"));
            }
          }
          baseTransformSymOp = new SymOp(ZERO_ROTATION, tempTranB);
          bestBaseTransformSymOp[currZ1][currZ2] = baseTransformSymOp;
          // SymOp.append defaults to rotation then translation. Break apart to have translation then rotation.
          targetTransformSymOp = new SymOp(ZERO_ROTATION, tempTranT);
          targetTransformSymOp = targetTransformSymOp.append(new SymOp(tempRotT, Tr_0_0_0));
          bestTargetTransformSymOp[currZ1][currZ2] = targetTransformSymOp;
          SymOp base = baseSymOps.get(baseDistMap.get(minZinds[currZ1][currZ2].getSortedIndex()));
          baseSymOp = new SymOp(base.rot, base.tr);
          SymOp target = targetSymOps.get(targetDistMap.get(minZinds[currZ1][currZ2].getReferenceIndex()));
          targetSymOp = new SymOp(target.rot, target.tr);
          bestBaseSymOp[currZ1][currZ2] = baseSymOp;
          bestTargetSymOp[currZ1][currZ2] = targetSymOp;

          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append(format("\n \t Base sym op %2d (%3d):  \t\t\t  Target sym op %2d (%3d):",
                    currZ1, minZinds[currZ1][currZ2].getSortedIndex(),
                    currZ2, minZinds[currZ1][currZ2].getReferenceIndex()));
            for (int k = 0; k < nAU * compareAtomsSize; k++) {
              int l = k * 3;
              double[] bXYZ = new double[nCoords];
              double[] tXYZ = new double[nCoords];
              double[] bXYZsym = new double[nCoords];
              double[] tXYZsym = new double[nCoords];
              arraycopy(baseAUoriginal[currZ1], 0, bXYZ, 0, nCoords);
              arraycopy(targetAUoriginal[currZ2], 0, tXYZ, 0, nCoords);
              arraycopy(baseAUoriginal[currZ1], 0, bXYZsym, 0, nCoords);
              arraycopy(targetAUoriginal[currZ2], 0, tXYZsym, 0, nCoords);
              applyCartesianSymOp(bXYZsym, bXYZsym, bestBaseSymOp[currZ1][currZ2]);
              applyCartesianSymOp(tXYZsym, tXYZsym, bestTargetSymOp[currZ1][currZ2]);

              stringBuilder.append(format("\n from %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f compare %9.3f %9.3f %9.3f from %9.3f %9.3f %9.3f",
                      bXYZ[l], bXYZ[l + 1], bXYZ[l + 2],
                      bXYZsym[l], bXYZsym[l + 1], bXYZsym[l + 2],
                      tXYZsym[l], tXYZsym[l + 1], tXYZsym[l + 2],
                      tXYZ[l], tXYZ[l + 1], tXYZ[l + 2]));
            }
          }
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append("\n Base Sym Op Translation: ").append(Arrays.toString(bestBaseSymOp[currZ1][currZ2].tr));
            stringBuilder.append(matrixToString(bestBaseSymOp[currZ1][currZ2].rot, currZ2 + currZ1 * z2,
                    "Base Sym Op Rot")).append("\n");
            stringBuilder.append(" Base Transform Translation: ").append(Arrays.toString(bestBaseTransformSymOp[currZ1][currZ2].tr));
            stringBuilder.append(matrixToString(bestBaseTransformSymOp[currZ1][currZ2].rot, currZ2 + currZ1 * z2,
                    "Base Transform Rot")).append("\n");
            stringBuilder.append(" Target Sym Op Translation: ").append(Arrays.toString(bestTargetSymOp[currZ1][currZ2].tr));
            stringBuilder.append(matrixToString(bestTargetSymOp[currZ1][currZ2].rot, currZ2 + currZ1 * z2,
                    "Target Sym Op Rot")).append("\n");
            stringBuilder.append(" Target Transform Translation: ").append(Arrays.toString(bestTargetTransformSymOp[currZ1][currZ2].tr));
            stringBuilder.append(matrixToString(bestTargetTransformSymOp[currZ1][currZ2].rot, currZ2 + currZ1 * z2,
                    "Target Transform Rot")).append("\n");
          }
          if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(
                    format(" \n Saved %3d: Base %2d IND %3d (Z'=%2d) \t Target %2d IND %3d (Z'=%2d) Diff: %8.4f\n",
                            j + i * targetLength, bIndex, baseInd, currZ1, tIndex, targetInd, currZ2, value));
            if (logger.isLoggable(Level.FINEST)) {
              stringBuilder.append(format(" %35s to %30s compare %30s", "Base transform", "overlap Coords:",
                      "Target overlap Coords from transform:"));
              for (int k = 0; k < compareAtomsSize; k++) {
                int l = k * 3;
                double[] bXYZ = new double[nCoords];
                double[] tXYZ = new double[nCoords];
                arraycopy(this.baseXYZ, baseInd * nCoords, bXYZ, 0, nCoords);
                arraycopy(this.targetXYZ, targetInd * nCoords, tXYZ, 0, nCoords);
                stringBuilder.append(format("\n from %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f compare %9.3f %9.3f %9.3f from %9.3f %9.3f %9.3f",
                        bXYZ[l], bXYZ[l + 1], bXYZ[l + 2],
                        baseXYZ[l], baseXYZ[l + 1], baseXYZ[l + 2],
                        targetXYZ[l], targetXYZ[l + 1], targetXYZ[l + 2],
                        tXYZ[l], tXYZ[l + 1], tXYZ[l + 2]));
              }
            }
          }
        }
        if (minDiffs[i][j] < min1) {
          min1 = minDiffs[i][j];
          minBIndex = bIndex;
          minTIndex = tIndex;
        }
      }
    }
    // Index of the most similar uniques... need to expand to Z'>1
    int baseBestZ = baseDistMap.get(baseAUDist[minBIndex].getIndex()) % z1;
    int targetBestZ = targetDistMap.get(targetAUDist[minTIndex].getIndex()) % z2;
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(format("\n Base min index: %d (Z'=%d) Target min Index %d (Z'=%d)", minBIndex, baseBestZ,
              minTIndex, targetBestZ));
    }
    // Identify which AUs are most similar between crystals.
    int maxLength = max(baseLength, targetLength);
    int[] baseTargetMap = new int[maxLength];
    // Ensure comparison selection is correct even when order is reversed.
    if (reverse) {
      for (int i = 0; i < baseLength; i++) {
        double minValue = Double.MAX_VALUE;
        for (int j = 0; j < targetLength; j++) {
          if (minValue > minDiffs[i][j]) {
            minValue = minDiffs[i][j];
            baseTargetMap[i] = uniqueTarget[j];
          }
        }
      }
    } else {
      for (int i = 0; i < targetLength; i++) {
        double minValue = Double.MAX_VALUE;
        for (int j = 0; j < baseLength; j++) {
          if (minValue > minDiffs[j][i]) {
            minValue = minDiffs[j][i];
            baseTargetMap[i] = uniqueBase[j];
          }
        }
      }
    }
    if (nAU == 1 && !strict) {
      // If only comparing one entity we are done. Return values.
      // Rg ranges from 1-->sqrt(3). Normalize to get from 0-->1
      double baseGyration = radiusOfGyration(bestBaseNAUs[baseBestZ][targetBestZ].clone(), massN);
      gyrations[0] = baseGyration;
      double targetGyration = radiusOfGyration(bestTargetNAUs[baseBestZ][targetBestZ].clone(), massN);
      gyrations[1] = targetGyration;
      if (useSym) {
        StringBuilder alchemicalAtoms = null;
        StringBuilder alchemicalAtoms2 = null;
        StringBuilder separation = new StringBuilder();
        for (int i = 0; i < z1; i++) {
          for (int j = 0; j < z2; j++) {
            separation.append("\n Atom Separation (A)  Description");
            if (z1 > 1 || z2 > 1) {
              separation.append(format(" (z1=%2d z2=%2d)\n  Base: Target:  Distance:\n", i, j));
            } else {
              separation.append(" \n");
            }
            for (int k = 0; k < compareAtomsSize; k++) {
              int index = k * 3;
              double value = rmsd(
                      new double[]{bestBaseNAUs[i][j][index], bestBaseNAUs[i][j][index + 1], bestBaseNAUs[i][j][index + 2]},
                      new double[]{bestTargetNAUs[i][j][index], bestTargetNAUs[i][j][index + 1],
                              bestTargetNAUs[i][j][index + 2]}, massN);
              if (logger.isLoggable(Level.INFO)) {
                if (printSym < value) {
                  // Collect atomic separation distances
                  separation.append(format(" %4d %4d  %14.6f\n", i * compareAtomsSize + comparisonAtoms[k] + 1,
                          j * compareAtomsSize + comparisonAtoms[k] + 1, value));
                  // Collect alchemical atoms

                  if (alchemicalAtoms == null) {
                    alchemicalAtoms = new StringBuilder(format(" %3d %3d  ", i, j)).append(i * compareAtomsSize + comparisonAtoms[k] + 1);
                    alchemicalAtoms2 = new StringBuilder(format(" %3d %3d  ", i, j)).append(j * compareAtomsSize + comparisonAtoms[k] + 1);
                  } else if (alchemicalAtoms.charAt(alchemicalAtoms.length() - 1) == '\n') {
                    alchemicalAtoms.append(format(" %3d %3d  ", i, j)).append(i * compareAtomsSize + comparisonAtoms[k] + 1);
                    alchemicalAtoms2.append(format(" %3d %3d  ", i, j)).append(j * compareAtomsSize + comparisonAtoms[k] + 1);
                  } else {
                    alchemicalAtoms.append(",").append(i * compareAtomsSize + comparisonAtoms[k] + 1);
                    alchemicalAtoms2.append(",").append(j * compareAtomsSize + comparisonAtoms[k] + 1);
                  }
                }
              }
            }
            if (useSave) {
              if (machineLearning) {
                saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs[i][j], comparisonAtoms, nAU,
                        "_c1", 0.000, compNum, save);
                saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs[i][j], comparisonAtoms, nAU,
                        "_c2", min1, compNum, save);
              } else {
                saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs[i][j], comparisonAtoms, nAU,
                        "_c1", compNum, save);
                saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs[i][j], comparisonAtoms, nAU,
                        "_c2", compNum, save);
              }
            }
            if (alchemicalAtoms != null) {
              alchemicalAtoms.append("\n");
              alchemicalAtoms2.append("\n");
            }
          }
        }
        if (alchemicalAtoms != null) {
          stringBuilder.append(separation);
          stringBuilder.append("\n Alchemical Atoms:\n");
          // If more than one entity being compared, specify which atoms belong to which molecule of each crystal.
          if (z1 > 1 || z2 > 1) {
            stringBuilder.append(" Base: (").append(baseLabel).append(")\n");
          }
          stringBuilder.append(alchemicalAtoms).append("\n");
          if (z1 > 1 || z2 > 1) {
            stringBuilder.append(" Target: (").append(targetLabel).append(")\n").append(alchemicalAtoms2).append("\n");
          }
        }
      }
      if (useSave && !useSym) {
        if (machineLearning) {
          saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs[baseBestZ][targetBestZ],
                  comparisonAtoms, nAU, "_c1", 0.000, compNum, save);
          saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs[baseBestZ][targetBestZ],
                  comparisonAtoms, nAU, "_c2", min1, compNum, save);
        } else {
          saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs[baseBestZ][targetBestZ],
                  comparisonAtoms, nAU, "_c1", compNum, save);
          saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs[baseBestZ][targetBestZ],
                  comparisonAtoms, nAU, "_c2", compNum, save);
        }
      }
      if (logger.isLoggable(Level.FINER) && useSym) {
        for (int i = 0; i < z1; i++) {
          for (int j = 0; j < z2; j++) {
            printSym(compareAtomsSize, file2, name2, bondList2, atoms2, forceField2, save, j);
          }
        }
      }
      return min1;
    }
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append(
              " Minimum RMSD_1 Between Unique Base and Target AUs:\n i  j (bInd tInd) RMSD_1\n");
      for (int i = 0; i < baseLength; i++) {
        for (int j = 0; j < targetLength; j++) {
          stringBuilder.append(
                  format(" %d %d (%4d %4d) %4.4f\n", i, j,
                          baseAUDist[uniqueBase[i]].getIndex(), targetAUDist[uniqueTarget[j]].getIndex(), minDiffs[i][j]));
        }
      }
    }

    // Coordinate arrays to save out structures at the end.
    double bestRMSD = Double.MAX_VALUE;
    // Reset values from minZdiffs... or allocate another array.
    for (double[] row : minZdiffs) {
      fill(row, Double.MAX_VALUE);
    }
    if (logger.isLoggable(Level.FINE)) {
      stringBuilder.append(format("\n  Trial     RMSD_1 (%8s)  RMSD_3 (%8s)  %8s  G(r1)    G(r2)\n",
              rmsdLabel, rmsdLabel, rmsdLabel));
    }
    baseAUDist_2 = new DoubleIndexPair[nBaseMols];
    targetAUDist_2 = new DoubleIndexPair[nTargetMols];
    // Begin comparison
    // Integer used only for user display logging.
    int currentComparison = 1;
    for (int l = 0; l < baseLength; l++) {
      final int bIndex = uniqueBase[l];
      final int baseInd = baseAUDist[bIndex].getIndex();
      final int currZ1 = baseDistMap.get(baseInd) % z1;

      if (l > 0) {
        // Reset coords for next comparison.
        arraycopy(baseXYZoriginal, 0, baseXYZ, 0, baseXYZoriginal.length);
      }
      //Re-prioritize based on center-most molecule
      prioritizeReplicates(baseXYZ, massN, massSum, baseCoM[0], compareAtomsSize, baseAUDist_2, baseInd, linkage);
      //Translate base system based on center-most molecule
      int baseAUIndex = baseAUDist_2[0].getIndex() * nCoords;
      arraycopy(baseXYZ, baseAUIndex, baseAU, 0, nCoords);
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < nAU; j++) {
          int molIndex = baseAUDist_2[j].getIndex() * nCoords;
          arraycopy(baseXYZ, molIndex, baseNAUs[i], j * nCoords, nCoords);
        }
      }
      // Translate base to the origin
      double[] translation = calculateTranslation(baseAU, massN);
      applyTranslation(baseAU, translation);
      applyTranslation(baseNAUs[1], translation);
      applyTranslation(baseNAUs[2], translation);
      applyTranslation(baseNAUs[3], translation);
      applyTranslation(baseXYZ, translation);
      //Update CoMs with translation
      centerOfMass(baseCoM[1], baseXYZ, massN, massSum, compareAtomsSize);
      if (useSym) {
        // Create new identity SymOp for comparison.
        baseTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
        // All target comparisons are performed based on the same base AU.
        SymOp base = baseSymOps.get(baseInd);
        baseSymOp = new SymOp(base.rot, base.tr);
        baseTransformSymOp = baseTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
        if (logger.isLoggable(Level.FINEST)) {
          stringBuilder.append(format("\n Base %d (%2d) to Origin Translation: ", l, bIndex))
                  .append(Arrays.toString(translation)).append("\n");
        }
      }

      // Acquire coordinates based on center 3 molecules
      if (logger.isLoggable(Level.FINEST)) {
        stringBuilder.append("\n Base 3 Conformations:\n");
      }
      for (int i = 0; i < 3; i++) {
        int baseIndex = baseAUDist_2[i].getIndex() * nCoords;
        arraycopy(baseXYZ, baseIndex, base3AUs, i * nCoords, nCoords);
      }
      translation = calculateTranslation(base3AUs, massN);
      applyTranslation(base3AUs, translation);
      applyTranslation(baseNAUs[2], translation);
      applyTranslation(baseNAUs[3], translation);
      applyTranslation(baseXYZ, translation);
      //Update CoMs with translation
      centerOfMass(baseCoM[2], baseXYZ, massN, massSum, compareAtomsSize);
      if (useSym) {
        baseTransformSymOp = baseTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
        if (logger.isLoggable(Level.FINEST)) {
          stringBuilder.append(format("\n Base %d (%2d) 2nd Translation: ", l, bIndex))
                  .append(Arrays.toString(translation)).append("\n");
        }
      }

      // Acquire coordinates for final comparison
      double maxDist = baseAUDist[baseAUDist.length - 1].getDoubleValue();
      if (logger.isLoggable(Level.FINEST)) {
        stringBuilder.append(format("\n System 1 farthest distance: %4.8f", sqrt(maxDist)));
      }
      for (int i = 0; i < nAU; i++) {
        int molIndex = baseAUDist_2[i].getIndex() * nCoords;
        arraycopy(baseXYZ, molIndex, baseNAUs[3], i * nCoords, nCoords);
      }
      translation = calculateTranslation(baseNAUs[3], massN);
      // array copy to baseNAUs
      applyTranslation(baseNAUs[3], translation);
      applyTranslation(baseXYZ, translation);
      if (useSym) {
        baseTransformSymOp = baseTransformSymOp.append(new SymOp(ZERO_ROTATION, translation));
        if (logger.isLoggable(Level.FINEST)) {
          stringBuilder.append(format("\n Base %d (%2d) Final Translation: ", l, bIndex))
                  .append(Arrays.toString(translation)).append("\n");
        }
      }

      for (int m = 0; m < targetLength; m++) {
        final int tIndex = uniqueTarget[m];
        final int targetInd = targetAUDist[tIndex].getIndex();
        final int currZ2 = targetDistMap.get(targetInd) % z2;
        if (strict || !reverse && baseTargetMap[m] == bIndex
                || reverse && baseTargetMap[l] == tIndex) {

          arraycopy(targetXYZoriginal, 0, targetXYZ, 0, targetXYZoriginal.length);

          //Re-prioritize based on central AU if different from first prioritization.
          prioritizeReplicates(targetXYZ, massN, massSum, targetCoM, compareAtomsSize, targetAUDist_2, targetInd, linkage);

          // Isolate the AU closest to center in the second crystal (targetAU).
          arraycopy(targetXYZ, targetAUDist_2[0].getIndex() * nCoords, targetAU, 0, nCoords);

          // Translate target to origin
          translation = calculateTranslation(targetAU, massN);
          // Move the center current AU
          applyTranslation(targetAU, translation);
          // Move the entire system
          applyTranslation(targetXYZ, translation);

          // Rotate the target molecule onto the base molecule.
          double[][] rotation = calculateRotation(baseAU, targetAU, massN);
          applyRotation(targetAU, rotation);
          applyRotation(targetXYZ, rotation);
          //Update center of masses with the first trans/rot
          centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);

          if (useSym) {
            // Reset values for printing symmetry operations
            targetTransformSymOp = new SymOp(ZERO_ROTATION, Tr_0_0_0);
            // Switch m center most molecules (looking for stereoisomers)
            SymOp target = targetSymOps.get(targetInd);
            targetSymOp = new SymOp(target.rot, target.tr);
            //Matrix order matters. Perform translation then rotation. Standard append applies rotation then translation.
//            targetTransformSymOp[currZ1][currZ2] = targetTransformSymOp[currZ1][currZ2].append(new SymOp(rotation, Tr_0_0_0));
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(ZERO_ROTATION, translation)).append(new SymOp(rotation, Tr_0_0_0));

            if (logger.isLoggable(Level.FINEST)) {
              stringBuilder.append(matrixToString(baseSymOp.rot, bIndex, "Base Sym Rot"));
              stringBuilder.append(format("\n Base %d Sym Op Translation: ", bIndex))
                      .append(Arrays.toString(baseSymOp.tr)).append("\n");
              stringBuilder.append(matrixToString(targetSymOp.rot, tIndex, "Target Sym Rot"));
              stringBuilder.append(format("\n Target %d Sym Op Translation: ", tIndex))
                      .append(Arrays.toString(targetSymOp.tr)).append("\n");
              stringBuilder.append(format("\n Trans target %d (%2d) sys to origin: \n\t", m, tIndex))
                      .append(Arrays.toString(translation)).append("\n");
              stringBuilder.append(matrixToString(rotation, tIndex, "1st Rot"));
            }
          }

          // At this point both systems have completed first rotation/translation
          //  Therefore both center-most molecules should be overlapped.
          // Need 3 for progressive alignment even if RMSD_1 or RMSD_2 is desired
          pairEntities(Math.max(3, nAU), 1);

          double checkRMSD1 = -3.0;
          double n1RMSD = -4.0;
          if (logger.isLoggable(Level.FINE)) {
            checkRMSD1 = rmsd(baseAU, targetAU, massN);
            for (int i = 0; i < nAU; i++) {
              int offset = i * nCoords;
              int molIndex = pairedAUs[i].getIndex() * nCoords;
              arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
            }
            n1RMSD = rmsd(baseNAUs[1], targetNAUs, massN);
            if (logger.isLoggable(Level.FINEST) && useSave) {
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseAU, comparisonAtoms, 1,
                      "_c1_1", compNum, save);
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs[1], comparisonAtoms, nAU,
                      "_c1_1N", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetAU, comparisonAtoms, 1,
                      "_c2_1", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, nAU,
                      "_c2_1N", compNum, save);
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
          // Update center of masses with the second trans/rot.
          centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
          if (useSym) {
            applyTranslation(targetAU, translation);
            applyRotation(targetAU, rotation);
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(ZERO_ROTATION, translation)).append(new SymOp(rotation, Tr_0_0_0));
            if (logger.isLoggable(Level.FINEST)) {
              stringBuilder.append(format("\n Target %d (%2d) 2nd Translation/Rotation: ", m, tIndex))
                      .append(Arrays.toString(translation)).append("\n");
              printSym(compareAtomsSize, file2, name2, bondList2, atoms2, forceField2, save, currZ2);
            }
          }

          double checkRMSD2 = -5.0;
          double n3RMSD = -6.0;
          if (logger.isLoggable(Level.FINE)) {
            checkRMSD2 = rmsd(base3AUs, target3AUs, massN);
            n3RMSD = rmsd(baseNAUs[2], targetNAUs, massN);
            if (logger.isLoggable(Level.FINEST) && useSave) {
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, base3AUs, comparisonAtoms, 3,
                      "_c1_3", compNum, save);
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs[2], comparisonAtoms, nAU,
                      "_c1_3N", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, target3AUs, comparisonAtoms, 3,
                      "_c2_3", compNum, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, nAU,
                      "_c2_3N", compNum, save);
            }
          }

          // Rotations 1 and 2 have been completed and both systems should be overlapped
          //  Isolate center most nAU from System 1 and matching molecules from System 2
          pairEntities(nAU, 2);

          for (int i = 0; i < nAU; i++) {
            int offset = i * nCoords;
            int molIndex = pairedAUs[i].getIndex() * nCoords;
            arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
          }

          translation = calculateTranslation(targetNAUs, massN);
          applyTranslation(targetNAUs, translation);
          applyTranslation(targetXYZ, translation);
          rotation = calculateRotation(baseNAUs[3], targetNAUs, massN);
          if (useSym) {
            applyTranslation(targetAU, translation);
            applyRotation(targetAU, rotation);
            targetTransformSymOp = targetTransformSymOp.append(new SymOp(ZERO_ROTATION, translation)).append(new SymOp(rotation, Tr_0_0_0));
//            targetTransformSymOp[currZ1][currZ2] = targetTransformSymOp[currZ1][currZ2].append(new SymOp(rotation, Tr_0_0_0));
            if (logger.isLoggable(Level.FINEST)) {
              stringBuilder.append(matrixToString(rotation, bIndex, "Target System Final Rotation"));
              stringBuilder.append(format("\n Target %d Final Translation: ", bIndex))
                      .append(Arrays.toString(translation)).append("\n");
              printSym(compareAtomsSize, file2, name2, bondList2, atoms2, forceField2, save, currZ2);
            }
          }
          if (logger.isLoggable(Level.FINEST) && useSave) {
            saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs[3], comparisonAtoms, nAU,
                    "_c1_N", compNum, save);
            saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, nAU,
                    "_c2_N", compNum, save);
          }
          applyRotation(targetNAUs, rotation);
          applyRotation(targetXYZ, rotation);
          double rmsdSymOp = rmsd(baseNAUs[3], targetNAUs, massN);

          if (logger.isLoggable(Level.FINEST) && useSave) {
            saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseNAUs[3], comparisonAtoms, nAU,
                    "_c1_N2", compNum, save);
            saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetNAUs, comparisonAtoms, nAU,
                    "_c2_N2", compNum, save);
          }

          double baseGyration = radiusOfGyration(baseNAUs[3], massN);
          double targetGyration = radiusOfGyration(targetNAUs, massN);

          if (logger.isLoggable(Level.FINE)) {
            int totalComparisons = (strict) ? baseLength * targetLength : maxLength;
            String output = format(
                    " %2d of %2d: %7.4f (%8.4f) %7.4f (%8.4f) %8.4f %8.4f %8.4f",
                    currentComparison, totalComparisons, checkRMSD1, n1RMSD, checkRMSD2, n3RMSD,
                    rmsdSymOp,
                    baseGyration, targetGyration);

            if (logger.isLoggable(Level.FINER)) {
              if (reverse) {
                output += format(" b: %2d t: %2d bt: %2d", baseAUDist[bIndex].getIndex(),
                        targetAUDist[uniqueTarget[m]].getIndex(),
                        targetAUDist[baseTargetMap[m]].getIndex());
              } else {
                output += format(" b: %2d t: %2d tb: %2d", baseAUDist[bIndex].getIndex(),
                        targetAUDist[uniqueTarget[m]].getIndex(),
                        baseAUDist[baseTargetMap[m]].getIndex());
              }
            }
            stringBuilder.append(output).append("\n");
          }

          if (useSym && rmsdSymOp < minZdiffs[currZ1][currZ2]) {
            minZdiffs[currZ1][currZ2] = rmsdSymOp;
            bestTargetTransformSymOp[currZ1][currZ2] = new SymOp(targetTransformSymOp.rot, targetTransformSymOp.tr);
            bestTargetSymOp[currZ1][currZ2] = new SymOp(targetSymOp.rot, targetSymOp.tr);
            bestBaseSymOp[currZ1][currZ2] = new SymOp(baseSymOp.rot, baseSymOp.tr);
            bestBaseTransformSymOp[currZ1][currZ2] = new SymOp(baseTransformSymOp.rot, baseTransformSymOp.tr);
          }
          if (rmsdSymOp < bestRMSD) {
            // Rg ranges from 1-->sqrt(3). Normalize to get from 0-->1
            gyrations[0] = baseGyration;
            gyrations[1] = targetGyration;
            bestRMSD = rmsdSymOp;
            baseBestZ = currZ1;
            targetBestZ = currZ2;
            arraycopy(baseNAUs[3], 0, bestBaseNAUs[currZ1][currZ2], 0, nAU * nCoords);
            arraycopy(targetNAUs, 0, bestTargetNAUs[currZ1][currZ2], 0, nAU * nCoords);
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
        saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs[baseBestZ][targetBestZ],
                comparisonAtoms, nAU, "_c1", 0.000, compNum, save);
        saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs[baseBestZ][targetBestZ],
                comparisonAtoms, nAU, "_c2", finalRMSD, compNum, save);
      } else {
        saveAssembly(file1, name1, bondList1, atoms1, forceField1, bestBaseNAUs[baseBestZ][targetBestZ],
                comparisonAtoms, nAU, "_c1", compNum, save);
        saveAssembly(file2, name2, bondList2, atoms2, forceField2, bestTargetNAUs[baseBestZ][targetBestZ],
                comparisonAtoms, nAU, "_c2", compNum, save);
      }
    }

    if (inertia) {
      bestBaseMandV = momentsOfInertia(bestBaseNAUs[baseBestZ][targetBestZ], massN, false, false, true);
      bestTargetMandV = momentsOfInertia(bestTargetNAUs[baseBestZ][targetBestZ], massN, false, false, true);
    }

    if (gyrationComponents) {
      bestBaseRg = radiusOfGyrationComponents(bestBaseNAUs[baseBestZ][targetBestZ], massN, true);
      bestTargetRg = radiusOfGyrationComponents(bestTargetNAUs[baseBestZ][targetBestZ], massN, true);
    }

    return finalRMSD;
  }

  /**
   * Compare the crystals within the SystemFilters that were inputted into the constructor of this
   * class.
   *
   * @param nAU                Number of asymmetric units to compare.
   * @param inflationFactor    Specify safety factor when generating replicates crystal.
   * @param matchTol           Tolerance to determine whether two AUs are the same (increases efficiency).
   * @param zPrime             Number of asymmetric units in first crystal (default attempts detection).
   * @param zPrime2            Number of asymmetric units in second crystal (default attempts detection).
   * @param unsharedA          List of atoms specific to first crystal.
   * @param unsharedB          List of atoms specific to second crystal.
   * @param alphaCarbons       Perform comparisons on only alpha carbons.
   * @param includeHydrogen    Perform comparisons without hydrogen atoms.
   * @param massWeighted       Perform comparisons with mass weighted coordinates (center of mass instead of geometric center).
   * @param crystalPriority    Prioritize most dense (0), least dense (1), or first inputted file (2).
   * @param strict             More intensive, but less efficient version of PAC.
   * @param save               Save out files of the resulting superposition.
   * @param restart            Try to restart from a previous job.
   * @param write              Save out a PAC RMSD file.
   * @param machineLearning    Save out CSV files for machine learning input (saves PDBs as well).
   * @param inertia            Compute moments of inertia for final clusters.
   * @param gyrationComponents Compute axial components for radius of gyration of final clusters.
   * @param linkage            Prioritize entities based on single, average, or complete linkage.
   * @param printSym           Print final symmetry operator used to superimpose mobile assembly onto static assembly.
   * @param lowMemory          Crystals will be read in as needed (slower performance, but less memory intensive)
   * @param pacFileName        The filename to use.
   * @return RunningStatistics Statistics for comparisons performed.
   */
  public RunningStatistics comparisons(int nAU, double inflationFactor, double matchTol, int zPrime,
                                       int zPrime2, String unsharedA, String unsharedB, boolean alphaCarbons, boolean includeHydrogen, boolean massWeighted,
                                       int crystalPriority, boolean strict, int save, boolean restart, boolean write,
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

    // Atom arrays from the 2nd assembly.
    MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
    Atom[] atoms2 = targetAssembly.getAtomArray();
    int nAtoms2 = atoms2.length;

    // Collect selected atoms.
    ArrayList<Integer> unique = new ArrayList<>(TopologyOptions.getUniqueAtoms(baseAssembly, "base", unsharedA));
    // TopologyOptions assumes we want indices for a concatenated XYZ object... recreate atom index.
    ArrayList<Integer> uniqueB = new ArrayList<>();
    for (Integer value : unique) {
      if (value % 3 == 0) {
        int index = (value / 3);
        if (!uniqueB.contains(index)) {
          uniqueB.add(index);
        }
      }
    }
    unique.clear();
    determineComparableAtoms(atoms1, unique, uniqueB, alphaCarbons, includeHydrogen);
    if (unique.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD in first crystal.");
      return null;
    }
    int[] comparisonAtoms = unique.stream().mapToInt(i -> i).toArray();
    unique.clear();
    // Collect selected atoms.
    unique = new ArrayList<>(TopologyOptions.getUniqueAtoms(targetAssembly, "target", unsharedB));
    ArrayList<Integer> uniqueT = new ArrayList<>();
    for (Integer value : unique) {
      if (value % 3 == 0) {
        int index = (value / 3);
        if (!uniqueT.contains(index)) {
          uniqueT.add(index);
        }
      }
    }
    unique.clear();
    determineComparableAtoms(atoms2, unique, uniqueT, alphaCarbons, includeHydrogen);

    if (unique.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD in second crystal.");
      return null;
    }
    int[] comparisonAtoms2 = unique.stream().mapToInt(i -> i).toArray();

    // Each ASU contains z * comparisonAtoms species so treat each species individually.
    int compareAtomsSize = comparisonAtoms.length;
    int compareAtomsSize2 = comparisonAtoms2.length;
    //Determine number of species within asymmetric unit.
    //TODO: Handle Z' > 1 for heterogeneous/co-crystals.
    int z1 = guessZPrime(zPrime, baseAssembly.getMolecules().size());
    int z2 = guessZPrime(zPrime2, targetAssembly.getMolecules().size());
    if (z1 > 1) {
      compareAtomsSize /= z1;
    }
    if (z2 > 1) {
      compareAtomsSize2 /= z2;
    }
    if (compareAtomsSize != compareAtomsSize2) {
      logger.warning(" Selected atom sizes differ between crystals.");
    }

    // When printing Sym Ops it is imperative to find all molecules (decrease matchTol).
    boolean useSym = printSym >= 0.0;
    if (useSym && (z1 > 1 || z2 > 1)) {
      matchTol = 0.00000002;
    }
    if(useSym && z1 != z2){
      logger.warning(" Comparisons to determine symmetry should have the same number of molecules in the asymmetric unit.");
    }

    Crystal baseCrystal = baseAssembly.getCrystal().getUnitCell();
    Crystal targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
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

    if (logger.isLoggable(Level.FINER)) {
      if (linkage == 0) {
        logger.finer(" Single linkage will be used.");
      } else if (linkage == 2) {
        logger.finer(" Complete linkage will be used.");
      } else if (linkage == 1) {
        logger.finer(" Average linkage will be used.");
      }
    }
    if (linkage != 1 && linkage != 0 && linkage != 2) {
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
    if (!lowMemory) {
      // Store necessary information in arrays for faster access.
      int size = (int) ceil((double) targetSize / (double) numProc);
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
          int atomSize = arrayAtom.length;
          for (int i = 0; i < atomSize; i++) {
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

    if (logger.isLoggable(Level.FINER)) {
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
      atoms1 = baseAssembly.getAtomArray();
      //Remove atoms not used in comparisons from the original molecular assembly (crystal 1).
      double[] reducedBaseCoords = reduceSystem(atoms1, comparisonAtoms);
      double baseDensity = baseCrystal.getUnitCell().getDensity(baseAssembly.getMass());
      if (baseCrystal == null || baseCrystal.aperiodic()) {
        logger.warning(" " + baseAssembly.getName() + " does not have a crystal.\n");
        continue;
      }

      //Setup for comparison with crystal specific information.
      // Density changes based on mass weighted flag, therefore use volume.

      // Estimate a radius that will include desired number of asymmetric units (inflationFactor).
      if (logger.isLoggable(Level.FINER)) {
        logger.finer(format(" Unit Cell Volume:  (Base) %4.2f (Target) %4.2f", baseCrystal.volume, targetCrystal.volume));
        logger.finer(format(" Unit Cell Symm Ops: (Base) %d (Target) %d", baseCrystal.getNumSymOps(), targetCrystal.getNumSymOps()));
        logger.finer(format(" Z': (Base) %d (Target) %d", z1, z2));
      }

      // When the system was read in, a replicates crystal may have been created to satisfy the cutoff.
      // Retrieve a reference to the unit cell (not the replicates crystal).
      // Here we will use the unit cell, to create a new replicates crystal that may be
      // a different size (i.e. larger).
      // Used for LMN specific replicates expansion (add LMN as input to generateInflatedSphere).
//      int[] baseLMN = determineExpansion(baseCrystal.getUnitCell(), reducedBaseCoords, comparisonAtoms, mass, nAU, z1, inflationFactor);
      int[] baseLMN = new int[3];
      double inflationFactorOrig = inflationFactor;
      baseXYZoriginal = generateInflatedSphere(baseCrystal.getUnitCell(),
              reducedBaseCoords, z1, nAU, linkage, mass, baseLMN, strict, baseSymOps, baseDistMap, inflationFactor);
      if(inflationFactorOrig != inflationFactor){
        logger.warning(format(" Default replicates crystal was too small for comparison. Increasing inflation factor from %9.3f to %9.3f",
            inflationFactor, inflationFactorOrig));
      }
      int nBaseCoords = baseXYZoriginal.length;
      baseXYZ = new double[nBaseCoords];
      arraycopy(baseXYZoriginal, 0, baseXYZ, 0, nBaseCoords);

      // Center of Masses for crystal 1.
      int nBaseMols = nBaseCoords / (compareAtomsSize * 3);
      if (baseAUDist == null || baseAUDist.length != nBaseMols) {
        baseAUDist = new DoubleIndexPair[nBaseMols];
      }
      if (baseCoM == null || baseCoM.length != nBaseMols) {
        baseCoM = new double[3][nBaseMols][3];
      }
      centerOfMass(baseCoM[0], baseXYZ, massN, massSum, compareAtomsSize);
      prioritizeReplicates(baseXYZ, massN, massSum, baseCoM[0], compareAtomsSize, baseAUDist, 0, linkage);

      for (int column = restartColumn; column < targetSize; column++) {
        int targetRank = column % numProc;
        if (targetRank == rank) {
          stringBuilder.setLength(0);
          long time = -System.nanoTime();
          int assemblyNum = column / numProc;
          if (!lowMemory) {
            targetCrystal = crystalCache[assemblyNum];
            atoms2 = atomCache[assemblyNum];
          } else {
            targetAssembly = targetFilter.getActiveMolecularSystem();
            targetCrystal = targetAssembly.getCrystal().getUnitCell();
            atoms2 = targetAssembly.getAtomArray();
          }
          if (targetCrystal == null || targetCrystal.aperiodic()) {
            if (!lowMemory) {
              logger.warning(" " + nameCache[assemblyNum] + " does not have a crystal.\n");
            } else {
              logger.warning(" " + targetAssembly.getName() + " does not have a crystal.\n");
            }
            continue;
          }
          double rmsd = -9.0;
          // Used to determine static/mobile crystal. Needed to switch back after comparison.
          boolean densityCheck;
          if (isSymmetric && row == column) {
            stringBuilder.append(
                    format("\n Comparing Model %4d (%s) of %s\n with      Model %4d (%s) of %s\n",
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
                    format("\n Comparing Model %4d (%s) of %s\n with      Model %4d (%s) of %s\n",
                            row + 1, baseCrystal.toShortString(), baseLabel,
                            column + 1, targetCrystal.toShortString(), targetLabel));

            //Remove atoms not used in comparisons from the original molecular assembly (crystal 2).
            double[] reducedTargetCoords = reduceSystem(atoms2, comparisonAtoms2);
            // Used for LMN specific replicates expansion (add LMN as input to generateInflatedSphere).
//            int[] targetLMN = determineExpansion(targetCrystal.getUnitCell(), reducedTargetCoords, comparisonAtoms2,
//                mass2, nAU, z2, inflationFactor);
            int[] targetLMN = new int[3];
            inflationFactorOrig = inflationFactor;
            targetXYZoriginal = generateInflatedSphere(targetCrystal.getUnitCell(), reducedTargetCoords, z2, nAU,
                linkage, mass2, targetLMN, strict, targetSymOps, targetDistMap, inflationFactor);
            if(inflationFactorOrig != inflationFactor){
              logger.warning(format(" Default replicates crystal was too small for comparison. Increasing inflation factor from %9.3f to %9.3f",
                  inflationFactor, inflationFactorOrig));
            }
            int nTargetCoords = targetXYZoriginal.length;
            targetXYZ = new double[nTargetCoords];
            arraycopy(targetXYZoriginal, 0, targetXYZ, 0, nTargetCoords);

            // Prioritize crystal order based on user specification (High/low density or file order).
            double densityMass = 0;
            if (!lowMemory) {
              for (Atom atom : atomCache[assemblyNum]) {
                densityMass += atom.getMass();
              }
            } else {
              densityMass = targetFilter.getActiveMolecularSystem().getMass();
            }
            double targetDensity = targetCrystal.getUnitCell().getDensity(densityMass);
            if (logger.isLoggable(Level.FINER)) {
              stringBuilder.append(
                      format("\n Base Density: %4.4f Target Density: %4.4f\n", baseDensity,
                              targetDensity));
            }
            densityCheck =
                    (crystalPriority == 1) ? baseDensity < targetDensity : baseDensity > targetDensity;

            if (logger.isLoggable(Level.FINE) || strict) {
              int z_1 = guessZPrime(zPrime, baseAssembly.getMolecules().size());
              int z_2 = guessZPrime(zPrime2, targetAssembly.getMolecules().size());
              if (z1 != z_1 || z2 != z_2) {
                stringBuilder.append(format("\n Different Z' value detected within file: \n" +
                    " WARNING: Z' assumed to be constant for all structures in input.\n " +
                        "       z1: was %2d (now %2d) z2: was %2d (now %2d)\n", z1, z_1, z2, z_2));
              }
            }

            // Center of Masses for crystal 1.
            int nTargetMols = targetXYZ.length / (compareAtomsSize2 * 3);
            if (targetAUDist == null || targetAUDist.length != nTargetMols) {
              targetAUDist = new DoubleIndexPair[nTargetMols];
            }
            if (targetCoM == null || targetCoM.length != nTargetMols) {
              targetCoM = new double[nTargetMols][3];
            }
            //Update center of mass and priority for each target crystal..
            centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
            prioritizeReplicates(targetXYZ, massN, massSum, targetCoM, compareAtomsSize, targetAUDist, 0, linkage);

            // Flip system order based on crystalPriority if needed.
            File file1;
            File file2;
            String name1;
            String name2;
            List<Bond> bondList1;
            List<Bond> bondList2;
            ForceField forceField1;
            ForceField forceField2;
            arraycopy(baseXYZoriginal, 0, baseXYZ, 0, nBaseCoords);
            if (densityCheck || crystalPriority == 2) {
              atoms1 = baseAssembly.getAtomArray();
              file1 = baseAssembly.getFile();
              name1 = baseAssembly.getName();
              bondList1 = baseAssembly.getBondList();
              forceField1 = baseAssembly.getForceField();
              if (!lowMemory) {
                atoms2 = atomCache[assemblyNum];
                file2 = fileCache[assemblyNum];
                name2 = nameCache[assemblyNum];
                bondList2 = bondCache[assemblyNum];
                forceField2 = forceFieldCache[assemblyNum];
              } else {
                MolecularAssembly mobileAssembly = targetFilter.getActiveMolecularSystem();
                atoms2 = mobileAssembly.getAtomArray();
                file2 = mobileAssembly.getFile();
                name2 = mobileAssembly.getName();
                bondList2 = mobileAssembly.getBondList();
                forceField2 = mobileAssembly.getForceField();
              }
            } else {
              double[] tempXYZ = new double[nBaseCoords];
              arraycopy(baseXYZ, 0, tempXYZ, 0, nBaseCoords);
              baseXYZ = new double[nTargetCoords];
              arraycopy(targetXYZ, 0, baseXYZ, 0, nTargetCoords);
              targetXYZ = new double[nBaseCoords];
              arraycopy(tempXYZ, 0, targetXYZ, 0, nBaseCoords);
              arraycopy(baseXYZoriginal, 0, tempXYZ, 0, nBaseCoords);
              baseXYZoriginal = new double[nTargetCoords];
              arraycopy(targetXYZoriginal, 0, baseXYZoriginal, 0, nTargetCoords);
              targetXYZoriginal = new double[nBaseCoords];
              arraycopy(tempXYZ, 0, targetXYZoriginal, 0, nBaseCoords);
              tempXYZ = new double[nCoords];
              arraycopy(reducedBaseCoords, 0, tempXYZ, 0, nCoords);
              arraycopy(reducedTargetCoords, 0, reducedBaseCoords, 0, nCoords);
              arraycopy(tempXYZ, 0, reducedTargetCoords, 0, nCoords);
              ArrayList<SymOp> also = new ArrayList<>(baseSymOps);
              baseSymOps = new ArrayList<>(targetSymOps);
              targetSymOps = new ArrayList<>(also);
              ArrayList<Integer> ali = new ArrayList<>(baseDistMap);
              baseDistMap = new ArrayList<>(targetDistMap);
              targetDistMap = new ArrayList<>(ali);
              DoubleIndexPair[] tempDIP = new DoubleIndexPair[nBaseMols];
              arraycopy(baseAUDist, 0, tempDIP, 0, nBaseMols);
              baseAUDist = new DoubleIndexPair[nTargetMols];
              arraycopy(targetAUDist, 0, baseAUDist, 0, nTargetMols);
              targetAUDist = new DoubleIndexPair[nBaseMols];
              arraycopy(tempDIP, 0, targetAUDist, 0, nBaseMols);
              double[][] tempDD = new double[nBaseMols][3];
              arraycopy(baseCoM[0], 0, tempDD, 0, nBaseMols);
              baseCoM = new double[3][nTargetMols][3];
              arraycopy(targetCoM, 0, baseCoM[0], 0, nTargetMols);
              targetCoM = new double[nBaseMols][3];
              arraycopy(tempDD, 0, targetCoM, 0, nBaseMols);
              atoms2 = baseAssembly.getAtomArray();
              file2 = baseAssembly.getFile();
              name2 = baseAssembly.getName();
              bondList2 = baseAssembly.getBondList();
              forceField2 = baseAssembly.getForceField();
              int[] tempAtoms = comparisonAtoms.clone();
              comparisonAtoms = comparisonAtoms2.clone();
              comparisonAtoms2 = tempAtoms.clone();
              int temp = z1;
              z1 = z2;
              z2 = temp;
              temp = baseSearchValue;
              baseSearchValue = targetSearchValue;
              targetSearchValue = temp;
              temp = nBaseCoords;
              nBaseCoords = nTargetCoords;
              nTargetCoords = temp;
              temp = nBaseMols;
              nBaseMols = nTargetMols;
              nTargetMols = temp;
              if (!lowMemory) {
                atoms1 = atomCache[assemblyNum];
                file1 = fileCache[assemblyNum];
                name1 = nameCache[assemblyNum];
                bondList1 = bondCache[assemblyNum];
                forceField1 = forceFieldCache[assemblyNum];
              } else {
                MolecularAssembly staticAssembly = targetFilter.getActiveMolecularSystem();
                atoms1 = staticAssembly.getAtomArray();
                file1 = staticAssembly.getFile();
                name1 = staticAssembly.getName();
                bondList1 = staticAssembly.getBondList();
                forceField1 = staticAssembly.getForceField();
              }
            }
            // Now that Base and Target structures are finalized, allocate variables.
            if (baseAU == null) {
              baseAU = new double[nCoords];
            }
            if (targetAU == null) {
              targetAU = new double[nCoords];
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
              baseNAUs = new double[4][nAU * nCoords];
            }
            if (targetNAUs == null) {
              targetNAUs = new double[nAU * nCoords];
            }
            if (bestBaseNAUs == null) {
              bestBaseNAUs = new double[z1][z2][nAU * nCoords];
            }
            if (bestTargetNAUs == null) {
              bestTargetNAUs = new double[z1][z2][nAU * nCoords];
            }
            if (useSym) {
              if (baseSymOps == null) {
                baseSymOps = new ArrayList<>();
              }
              if (targetSymOps == null) {
                targetSymOps = new ArrayList<>();
              }
              if (baseAUoriginal == null || baseAUoriginal.length != z1) {
                baseAUoriginal = new double[z1][nCoords];
              }
              if (targetAUoriginal == null || targetAUoriginal.length != z2) {
                targetAUoriginal = new double[z2][nCoords];
              }
            }

            if (useSym) {
              for (int i = 0; i < z1; i++) {
                arraycopy(reducedBaseCoords, i * compareAtomsSize * 3, baseAUoriginal[i], 0, compareAtomsSize * 3);
              }
              for (int i = 0; i < z2; i++) {
                arraycopy(reducedTargetCoords, i * compareAtomsSize2 * 3, targetAUoriginal[i], 0, compareAtomsSize2 * 3);
              }
            }

            if (logger.isLoggable(Level.FINER)) {
              logger.finer(format(" %3d Molecules in Base Crystal \t %3d Molecules in Target Crystal",
                      (nBaseCoords / 3) / compareAtomsSize, (nTargetCoords / 3) / compareAtomsSize));
            }

            if (logger.isLoggable(Level.FINEST) && save > 0) {
              saveAssembly(file1, name1, bondList1, atoms1, forceField1, baseXYZoriginal, comparisonAtoms,
                      nBaseCoords / 3, "_c1_rep", -1, save);
              saveAssembly(file2, name2, bondList2, atoms2, forceField2, targetXYZoriginal, comparisonAtoms2,
                      nTargetCoords / 3, "_c2_rep", -1, save);
            }
            // Compute the PAC RMSD.
            rmsd = compare(file1, name1, bondList1, atoms1, forceField1, file2, name2, bondList2, atoms2, forceField2,
                    z1, z2, compareAtomsSize, nAU,
                    baseSearchValue, targetSearchValue, matchTol, row * targetSize + column,
                    strict, save, machineLearning, linkage, inertia, gyrationComponents);
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
            if (useSym) {
              StringBuilder sbSO = new StringBuilder(format("\n Sym Op to move %s onto %s:\nsymop ", name2,
                      name1));
              StringBuilder sbInv = new StringBuilder(format("\n Inverted Sym Op to move %s onto %s:\nsymop ", name1, name2));
              for (int i = 0; i < z1; i++) {
                for (int j = 0; j < z2; j++) {
                  if (bestTargetTransformSymOp[i][j] != null) {
                    bestTargetTransformSymOp[i][j] = bestTargetTransformSymOp[i][j].prepend(bestTargetSymOp[i][j]);
                    // Apply inverse of base operations:
                    // For orthogonal matrices the inverse matrix = the transpose. True iff det(A)== +/- 1.
                    // No inverse if det(A)==0.
                    bestBaseTransformSymOp[i][j] = bestBaseTransformSymOp[i][j].prepend(bestBaseSymOp[i][j]);
                    bestTargetTransformSymOp[i][j] = bestTargetTransformSymOp[i][j].append(
                            invertSymOp(bestBaseTransformSymOp[i][j]));
                    // Default to only printing out one sym Op per pair.
                    if (strict || logger.isLoggable(Level.FINEST) || z1 != z2 || i == j) {
                      SymOp inverted = invertSymOp(bestTargetTransformSymOp[i][j]);
                      if (logger.isLoggable(Level.FINEST)) {
                        stringBuilder.append(
                                format("\n Sym Op to move %s (%2d) onto %s (%2d):\n", name2, j,
                                        name1, i)).append(bestTargetTransformSymOp[i][j]).append("\n");
                        stringBuilder.append(
                                format("\n\n Inverted Sym Op to move %s (%2d) onto %s (%2d):\n", name1, i,
                                        name2, j)).append(inverted).append("\n");
                      }
                      sbSO.append(format("    %2d     %2d", j, i)).append(asMatrixString(bestTargetTransformSymOp[i][j]));
                      sbInv.append(format("    %2d     %2d", j, i)).append(asMatrixString(inverted));
                      if (i + 1 < z1 || j + 1 < z2) {
                        sbSO.append("\\\n");
                        sbInv.append("\\\n");
                      }

                      if (logger.isLoggable(Level.FINER)) {
                        stringBuilder.append(format("\n Sym Op Application from Scratch: Target: %s (%2d) onto Base: %s (%2d)",
                                name2, j, name1, i));
                        double[] symMol = new double[nCoords];
                        for (int k = 0; k < compareAtomsSize; k++) {
                          int l = k * 3;
                          double[] xyz = new double[]{targetAUoriginal[j][l], targetAUoriginal[j][l + 1],
                                  targetAUoriginal[j][l + 2]};
                          applyCartesianSymOp(xyz, xyz, bestTargetTransformSymOp[i][j]);
                          symMol[l] = xyz[0];
                          symMol[l + 1] = xyz[1];
                          symMol[l + 2] = xyz[2];
                          stringBuilder.append(format("\n %8.3f %8.3f %8.3f moved to %8.3f %8.3f %8.3f compared to %8.3f %8.3f %8.3f",
                                  targetAUoriginal[j][l], targetAUoriginal[j][l + 1], targetAUoriginal[j][l + 2],
                                  symMol[l], symMol[l + 1], symMol[l + 2],
                                  baseAUoriginal[i][l], baseAUoriginal[i][l + 1], baseAUoriginal[i][l + 2]));
                        }
                        if (save > 0) {
                          saveAssembly(file2, name2, bondList2, atoms2, forceField2, symMol, comparisonAtoms, 1,
                                  "_moved", 0, save);
                        }

                        double value = rmsd(baseAUoriginal[i], symMol, massN);
                        double[] symMol2 = new double[nCoords];
                        stringBuilder.append("\n\n Sym Op Inverse Application:");
                        for (int k = 0; k < compareAtomsSize; k++) {
                          int l = k * 3;
                          double[] xyz = new double[]{symMol[l], symMol[l + 1],
                                  symMol[l + 2]};
                          applyCartesianSymOp(xyz, xyz, inverted);
                          symMol2[l] = xyz[0];
                          symMol2[l + 1] = xyz[1];
                          symMol2[l + 2] = xyz[2];
                          stringBuilder.append(format("\n %8.3f %8.3f %8.3f moved to %8.3f %8.3f %8.3f compared to %8.3f %8.3f %8.3f",
                                  symMol[l], symMol[l + 1], symMol[l + 2],
                                  symMol2[l], symMol2[l + 1], symMol2[l + 2],
                                  targetAUoriginal[j][l], targetAUoriginal[j][l + 1], targetAUoriginal[j][l + 2]));
                        }
                        stringBuilder.append(
                                format("\n\n Sym Op RMSD: %8.8f \n\n", value));
                      }
                    }
                  } else {
                    if (logger.isLoggable(Level.INFO)) {
                      logger.info(format(" Symmetry Operator %2d %2d was filtered out. " +
                              "Try using a single asymmetric unit (--na 1).", i, j));
                    }
                  }
                }
              }
              stringBuilder.append(sbSO.toString()).append("\n").append(sbInv.toString()).append("\n");
            }
            // Reset values for base system if they were swapped.
            if (!densityCheck && crystalPriority != 2) {
              baseXYZ = new double[nTargetCoords];
              arraycopy(targetXYZ, 0, baseXYZ, 0, nTargetCoords);
              baseXYZoriginal = new double[nTargetCoords];
              arraycopy(targetXYZoriginal, 0, baseXYZoriginal, 0, nTargetCoords);
              arraycopy(reducedTargetCoords, 0, reducedBaseCoords, 0, nCoords);
              baseSymOps = new ArrayList<>(targetSymOps);
              baseDistMap = new ArrayList<>(targetDistMap);
              int[] tempAtoms = comparisonAtoms.clone();
              comparisonAtoms = comparisonAtoms2.clone();
              comparisonAtoms2 = tempAtoms.clone();
              DoubleIndexPair[] tempDIP = new DoubleIndexPair[nBaseMols];
              arraycopy(baseAUDist, 0, tempDIP, 0, nBaseMols);
              baseAUDist = new DoubleIndexPair[nTargetMols];
              arraycopy(targetAUDist, 0, baseAUDist, 0, nTargetMols);
              targetAUDist = new DoubleIndexPair[nBaseMols];
              arraycopy(tempDIP, 0, targetAUDist, 0, nBaseMols);
              double[][] tempDD = new double[nBaseMols][3];
              arraycopy(baseCoM[0], 0, tempDD, 0, nBaseMols);
              baseCoM = new double[3][nTargetMols][3];
              arraycopy(targetCoM, 0, baseCoM[0], 0, nTargetMols);
              targetCoM = new double[nBaseMols][3];
              arraycopy(tempDD, 0, targetCoM, 0, nBaseMols);
              int temp = z1;
              z1 = z2;
              z2 = temp;
              temp = baseSearchValue;
              baseSearchValue = targetSearchValue;
              targetSearchValue = temp;
              nBaseCoords = nTargetCoords;
              nBaseMols = nTargetMols;
            }
          }
          myDistances[myIndex] = rmsd;
          myIndex++;
          if (stringBuilder.length() > 0) {
            logger.info(stringBuilder.toString());
          }
        }
        if (lowMemory) {
          targetFilter.readNext(false, false, (column + 1) % numProc == rank);
        }
      }
      restartColumn = 0;
      if (lowMemory) {
        targetFilter.readNext(true, false, 0 == rank);
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
    if (logger.isLoggable(Level.FINE) && maxTime > Double.MIN_VALUE && maxTime != minTime) {
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
                    "           %14.8f %14.8f %14.8f \\\n" +
                    "           %14.8f %14.8f %14.8f \\\n" +
                    "           %14.8f %14.8f %14.8f ",
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
   * @param value  Potential new value if it is not already in list.
   * @param tol    Tolerance that determine whether values are equal.
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
   * @param rot            Rotation matrix to add.
   * @param totalTransform Array to be updated with rotation (4x4).
   * @param prepend        If true prepend translation, false append to end.
   */
  public void addRotation(double[][] rot, double[][] totalTransform, boolean prepend) {
    double[][] transform = new double[][]{
            {rot[0][0], rot[0][1], rot[0][2], 0.0},
            {rot[1][0], rot[1][1], rot[1][2], 0.0},
            {rot[2][0], rot[2][1], rot[2][2], 0.0},
            {0.0, 0.0, 0.0, 1.0}};
    transform = (prepend) ? mat4Mat4(totalTransform, transform) : mat4Mat4(transform, totalTransform);
    int length = totalTransform.length;
    for (int i = 0; i < length; i++) {
      arraycopy(transform[i], 0, totalTransform[i], 0, totalTransform[i].length);
    }
  }

  /**
   * Accumulate translations (matrix multiplication)
   *
   * @param translation    Translation matrix to add.
   * @param totalTransform Array to be updated with translation (4x4).
   * @param prepend        If true prepend translation, false append to end.
   */
  public static void addTranslation(double[] translation, double[][] totalTransform,
                                    boolean prepend) {
    double[][] transform = new double[][]{
            {1.0, 0.0, 0.0, translation[0]},
            {0.0, 1.0, 0.0, translation[1]},
            {0.0, 0.0, 1.0, translation[2]},
            {0.0, 0.0, 0.0, 1.0}};
    transform = (prepend) ? mat4Mat4(totalTransform, transform) : mat4Mat4(transform, totalTransform);
    int length = totalTransform.length;
    for (int i = 0; i < length; i++) {
      arraycopy(transform[i], 0, totalTransform[i], 0, totalTransform[i].length);
    }
  }

  /**
   * Calculate the center of mass for a given set of masses for the asymmetric unit and coordinates
   * (xyz)
   *
   * @param centersOfMass Returned center of mass for each asymmetric unit
   * @param coords        Coordinates of every atom in system.
   * @param mass          Masses of each atom in asymmetric unit.
   * @param massSum       Sum of masses within asymmetric unit.
   * @param nAtoms        Number of atoms in an entity.
   */
  private static void centerOfMass(double[][] centersOfMass, double[] coords, double[] mass,
                                   double massSum, int nAtoms) {
    int nAU = coords.length / (nAtoms * 3);
    for (int i = 0; i < nAU; i++) {
      int auIndex = i * nAtoms * 3;
      double[] comI = new double[3];
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = auIndex + j * 3;
        double m = mass[j];
        for (int k = 0; k < 3; k++) {
          comI[k] += coords[atomIndex + k] * m;
        }
      }
      for (int j = 0; j < 3; j++) {
        comI[j] /= massSum;
      }
      centersOfMass[i] = comI;
    }
  }

  /**
   * Determine if replicates crystal is large enough for approximate cluster to be used in PAC comparison.
   *
   * @param xyz              XYZ coordinates of structures.
   * @param massN            Masses per atom.
   * @param massSum          Sum of masses.
   * @param compareAtomsSize Number of atoms being compared from each structure.
   * @param nAU              Number of structures to compare between crystals.
   * @param linkage          Linkage method to be used.
   * @param startLMN         Current values for L x M x N (overwritten with new values)
   * @return Whether the LMN values have changed. If yes, recalculate replicates.
   */
  private static boolean checkInflatedSphere(double[] xyz, double[] massN,
         double massSum, int compareAtomsSize, int nAU, int linkage, int[] startLMN,
         ArrayList<SymOp> symOps, ArrayList<Integer> indexOrder) {
    // Recalculate center of mass based on XYZ order ([0] is the closest to center).
    double[][] centerOfMass = new double[xyz.length/3][3];
    centerOfMass(centerOfMass, xyz, massN, massSum, compareAtomsSize);
    DoubleIndexPair[] auDist = new DoubleIndexPair[xyz.length/3 / compareAtomsSize];
    prioritizeReplicates(xyz, massN, massSum, centerOfMass, compareAtomsSize, auDist, 0, linkage);

    boolean redo = false;
    for(int i = 0; i < nAU; i++) {
      if(logger.isLoggable(Level.FINE)) {
        logger.fine(format(" i: %3d\t auDist Index: %3d \t indexOrder: %3d\t indget(au): %3d\t ReplicatesVector:  %2d (%2d) %2d (%2d) %2d (%2d)",
            i, auDist[i].getIndex(), indexOrder.get(i), indexOrder.get(auDist[i].getIndex()),
            symOps.get(indexOrder.get(auDist[i].getIndex())).replicatesVector[0] + 1, startLMN[0],
            symOps.get(indexOrder.get(auDist[i].getIndex())).replicatesVector[1] + 1, startLMN[1],
            symOps.get(indexOrder.get(auDist[i].getIndex())).replicatesVector[2] + 1, startLMN[2]));
      }
      int[] lmn = symOps.get(indexOrder.get(auDist[i].getIndex())).replicatesVector;
      if(lmn[0] == 0 || lmn[0] == startLMN[0] - 1){
        redo = true;
      }
      if(lmn[1] == 0 || lmn[1] == startLMN[1] - 1){
        redo = true;
      }
      if(lmn[2] == 0 || lmn[2] == startLMN[2] - 1){
        redo = true;
      }
      if(redo){
        break;
      }
    }
    if(redo && logger.isLoggable(Level.FINE)){
      logger.fine(" REDO EXPANSION.");
    }
    return redo;
  }

  /**
   * Deep copy of values from one array to another.
   *
   * @param newArray Destination of values.
   * @param oldArray Current location of values.
   */
  private static void copyArrayValues(double[][] newArray, double[][] oldArray) {
    int arrayLength = newArray.length;
    for (int i = 0; i < arrayLength; i++) {
      arraycopy(oldArray[i], 0, newArray[i], 0, newArray[i].length);
    }
  }

  /**
   * Determine the indices of the atoms from the assembly that are active for this comparison.
   *
   * @param atoms           Atoms we potentially wish to use in comparison.
   * @param indices         Array list containing atom indices that will be used for this comparison.
   * @param unique          Indices that are unique to the system and should not be included in comparison.
   * @param alphaCarbons    Boolean whether to include only alpha carbons/nitrogens.
   * @param includeHydrogen Boolean whether to include hydrogens.
   */
  private static void determineComparableAtoms(Atom[] atoms, ArrayList<Integer> indices, ArrayList<Integer> unique,
                                               boolean alphaCarbons, boolean includeHydrogen) {
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      if (!unique.contains(i)) {
        Atom atom = atoms[i];
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
        // Reset all atoms to active once the selection is recorded.
        atom.setActive(true);
      }
    }
  }

  /**
   * This method calls <code>world.gather</code> to collect numProc PAC RMSD values.
   *
   * @param row               Current row of the PAC RMSD matrix.
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
   * distribution of replicates to facilitate comparisons that go beyond lattice parameters.
   *
   * @param unitCell      Crystal to define expansion.
   * @param reducedCoords Coordinates of asymmetric unit we wish to expand.
   * @param zPrime        Number of molecules in asymmetric unit.
   * @param nAU           Number of asymmetric units to compare in final shell (used to determine expansion guess).
   * @param linkage       Type of linkage to be used in comparison (single, average, or complete).
   * @param mass          Masses for atoms within reduced asymmetric unit.
   * @param lmn           Replicates lattice vecor lengths (L x M x N).
   * @param symOps        List of symmetry operators used to expand to replicates crystal.
   * @param indexOrder    Sorted list of index values (e.g. index 0 is contained in input file but may not be first).
   * @param inflationFactor Scalar to over expand replicates crystal to obtain all necessary orientations.
   * @return double[] containing the coordinates for the expanded crystal.
   */
  private static double[] generateInflatedSphere(final Crystal unitCell, final double[] reducedCoords, final int zPrime,
                   final int nAU, final int linkage, final double[] mass, int[] lmn, boolean strict,
                                                 ArrayList<SymOp> symOps, ArrayList<Integer> indexOrder, double inflationFactor) {
    symOps.clear();
    indexOrder.clear();
    // Num coords in asymmetric unit
    int nCoords = reducedCoords.length;
    // Num atoms in asymmetric unit
    int nAtoms = nCoords / 3;
    // Num atoms per molecule
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

    double volume = unitCell.volume / unitCell.getNumSymOps() / zPrime;
    double radius = cbrt(0.75 / PI * volume * max(3, nAU) * inflationFactor);
    Crystal replicatesCrystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, radius * 2.0, lmn);
    // Used for LMN specific replicates expansion.
//    Crystal replicatesCrystal = new ReplicatesCrystal(unitCell, lmn[0], lmn[1], lmn[2], radius * 2.0);
//    startLMN = lmn;
    // Symmetry coordinates for each molecule in replicates crystal
    int nSymm = replicatesCrystal.getNumSymOps();
    double[][] xS = new double[nSymm][nAtoms];
    double[][] yS = new double[nSymm][nAtoms];
    double[][] zS = new double[nSymm][nAtoms];

    int numEntities = nSymm * zPrime;
    // Cartesian center of each molecule
    double[][] cartCenterOfMass = new double[numEntities][3];

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
        double zMass = 0.0;
        for (int i = 0; i < zAtoms; i++) {
          double m = mass[i];
          int coordIndex = symIndex + i;
          centerOfMass[0] += xS[iSym][coordIndex] * m;
          centerOfMass[1] += yS[iSym][coordIndex] * m;
          centerOfMass[2] += zS[iSym][coordIndex] * m;
          zMass += m;
        }
        centerOfMass[0] /= zMass;
        centerOfMass[1] /= zMass;
        centerOfMass[2] /= zMass;

        double[] translate = moveIntoCrystal(replicatesCrystal, centerOfMass);
        //Convert sym op into cartesian for later use.
        double[] trans = symOp.tr.clone();
        replicatesCrystal.toCartesianCoordinates(trans, trans);
        int translateLength = translate.length;
        for (int i = 0; i < translateLength; i++) {
          trans[i] += translate[i];
        }
        symOps.add(new SymOp(rot, trans, symOp.replicatesVector));
        for (int i = 0; i < zAtoms; i++) {
          int coordIndex = symIndex + i;
          xS[iSym][coordIndex] += translate[0];
          centerOfMass[0] += translate[0];
          yS[iSym][coordIndex] += translate[1];
          centerOfMass[1] += translate[1];
          zS[iSym][coordIndex] += translate[2];
          centerOfMass[2] += translate[2];
        }

        // Save CoM cartesian coordinates
        cartCenterOfMass[iSym * zPrime + zp] = centerOfMass;
      }
    }

    //Determine molecular distances to "center" of sphere.
    //  In PACCOM the center is the geometric average of coordinates.
    //  In FFX the center is the center of the replicates crystal.
    DoubleIndexPair[] auDist = new DoubleIndexPair[numEntities];
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
    logger.fine(format(" Expanded Crystal Center: %14.8f %14.8f %14.8f",
        cartCenter[0], cartCenter[1], cartCenter[2]));

    for (int i = 0; i < numEntities; i++) {
      // Then compute Euclidean distance from Cartesian center of the replicates cell
      auDist[i] = new DoubleIndexPair(i, dist2(cartCenter, cartCenterOfMass[i]));
    }

    // Sort the molecules by their distance from the center.
    // Note that the smallest distances are first in the array after the sort.
    sort(auDist);
    for (DoubleIndexPair molsDist : auDist) {
      indexOrder.add(molsDist.getIndex());
    }

    if (logger.isLoggable(Level.FINEST)) {
      logger.finest("\n Copy  SymOp        Distance");
    }
    double[] systemCoords = new double[numEntities * zAtoms * 3];
    for (int n = 0; n < nSymm; n++) {
      for (int zp = 0; zp < zPrime; zp++) {
        int index = n * zPrime + zp;
        // Current molecule
        int iSym = auDist[index].getIndex();
        double distance = auDist[index].getDoubleValue();
        if (logger.isLoggable(Level.FINEST) && n < 30) {
          logger.finest(format(" %4d  %5d  %14.8f", index, iSym, sqrt(distance)));
        }

        // Create a new set of Atoms for each SymOp molecule
        for (int i = 0; i < zAtoms; i++) {
          int symIndex = index * zAtoms * 3;
          int atomIndex = symIndex + i * 3;
          int conversion = iSym / zPrime;
          int shift = i + (iSym % zPrime) * zAtoms;
          systemCoords[atomIndex] = xS[conversion][shift];
          systemCoords[atomIndex + 1] = yS[conversion][shift];
          systemCoords[atomIndex + 2] = zS[conversion][shift];
        }
      }
    }
    double massSum = Arrays.stream(mass).sum();
    if(logger.isLoggable(Level.FINE)) {
      logger.fine(" Checking replicates crystal.");
      logger.fine(format(" SysCoords: %3d (%3d)(%3d), reducedCoords Size: %3d (%3d)\n mass Size %3d, massSum: %6.3f, nAU: %2d, nCoords: %3d\n " +
              "auDist Size: %3d, l: %2d, m: %2d, n: %2d, numEntities: %3d, nAtoms: %3d, zPrime: %3d, zAtoms: %3d",
          systemCoords.length, systemCoords.length / 3, systemCoords.length / 3 / zAtoms, reducedCoords.length,
          reducedCoords.length / 3, mass.length, massSum, nAU, nCoords, auDist.length, lmn[0], lmn[1], lmn[2],
           numEntities, nAtoms, zPrime, zAtoms));
    }
    if (strict && checkInflatedSphere(systemCoords, mass, massSum, zAtoms, nAU, linkage, lmn, symOps, indexOrder)) {
        // Used for LMN specific replicates expansion.
//        startLMN = determineExpansion(unitCell, reducedCoords, comparisonAtoms, mass, nAU, zPrime, inflationFactor);
      inflationFactor += 5;
        systemCoords = generateInflatedSphere(unitCell, reducedCoords, zPrime, nAU, linkage, mass, lmn, strict,
            symOps, indexOrder, inflationFactor);
    }
    return systemCoords;
  }

  /**
   * Try to automatically determine number of species in asymmetric unit (only works for molecules).
   *
   * @param zPrime      User input overrides detection method.
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
   * @param matrix      Values to present
   * @param index       Identifier
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
   * @param com     Center of mass (x, y, z) for the object of concern
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
   * @param allCoords       Coordinates for every atom in replicates crystal ([x1, y1, z1, x2, y2,
   *                        z2...].
   * @param auDist          Prioritization of molecules.
   * @param auCoords        Coordinates for single AU.
   * @param nCoords         Number of coordinates in an AU (number of atoms * 3).
   * @param upperLimit      The largest number of unique AUs (0 for no upper limit).
   * @param strict         Search entire replicates crystal if true, otherwise only the expected.
   * @param nAUinReplicates Number of AUs in replicates crystal.
   * @param massN           Array containing masses for each atom in AU.
   * @param matchTol        Tolerance to determine whether two AUs are the same.
   */
  private void numberUniqueAUs(double[] allCoords, DoubleIndexPair[] auDist, double[] auCoords,
                               int nCoords, int upperLimit, boolean strict,
                               int nAUinReplicates, double[] massN, double matchTol) {
    // uniqueDiffs is only recorded for logging... could remove later.
    // List of differences (RMSD_1) for AUs in replicates crystal.
    List<Double> uniqueDiffs = new ArrayList<>();
    ArrayList<double[]> tempListXYZ = new ArrayList<>();
    tempListIndices.add(0);
    uniqueDiffs.add(rmsd(auCoords, auCoords, massN));
    tempListXYZ.add(auCoords);
    // Start from 1 as zero is automatically added.
    int index = 1;
    // Determine number of conformations in first crystal
    boolean useSym = printSym >= 0.0;
    int numConfCheck = (strict || upperLimit <= 0 || useSym) ? nAUinReplicates : upperLimit;
    if (logger.isLoggable(Level.FINEST)) {
      logger.finest("RMSD Differences in Replicates Crystal:");
    }
    double[] target = new double[nCoords];
    arraycopy(tempListXYZ.get(0), 0, target, 0, nCoords);
    translate(target, massN);
    while ((uniqueDiffs.size() < numConfCheck)) {
      double[] baseCheckMol = new double[nCoords];
      int auIndex = auDist[index].getIndex();
      arraycopy(allCoords, auIndex * nCoords, baseCheckMol, 0, nCoords);
      translate(baseCheckMol, massN);
      rotate(target, baseCheckMol, massN);
      double value = rmsd(target, baseCheckMol, massN);
      if (addLooseUnequal(uniqueDiffs, value, matchTol) || useSym) {
        if (logger.isLoggable(Level.FINEST)) {
          logger.finest(format(" Sorted Index: %4d  Dist: %18.16f", auIndex, value));
        }
        tempListIndices.add(index);
        tempListXYZ.add(baseCheckMol);
      }
      index++;
      if (index >= auDist.length) {
        break;
      }
    }
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(" RMSD_1 from 1st AU:\n i RMSD     AU Index");
      int numUnique = uniqueDiffs.size();
      for (int i = 0; i < numUnique; i++) {
        logger.finer(format(" %d %4.4f    %4d", i, uniqueDiffs.get(i),
                auDist[tempListIndices.get(i)].getIndex()));
      }
    }
  }

  /**
   * Pair species between two crystals based on distances between centers.
   *
   * @param desiredAUs    Number of pairs to determine.
   * @param comparisonNum Determine which centers of mass to use for first crystal.
   */
  private void pairEntities(int desiredAUs, int comparisonNum) {
    tempListIndices.clear();
    // List of indexes for second system.
    for (DoubleIndexPair doubleIndexPair : targetAUDist_2) {
      // Only search molecules within range of the desired number of molecules.
      // Must have enough molecules for matching (using exhaustive till better heuristic is determined)
      tempListIndices.add(doubleIndexPair.getIndex());
    }

    // Compare distances between center of masses from system 1 and 2.
    for (int i = 0; i < desiredAUs; i++) {
      double minDist = Double.MAX_VALUE;
      Integer minIndex = -1;
      double[] baseCoMCurr = baseCoM[comparisonNum][baseAUDist_2[i].getIndex()];
      for (Integer target : tempListIndices) {
        double dist = dist2(baseCoMCurr, targetCoM[target]);
        if (dist < minDist) {
          minIndex = target;
          minDist = dist;
        }
        if (abs(minDist) < MATCH_TOLERANCE) {
          // Distance between center of masses is ~0 is the best scenario assuming no coordinate overlaps.
          if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append("\n \tExit out of match loop.\n");
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
                        baseCoM[comparisonNum][baseAUDist_2[i].getIndex()][0],
                        baseCoM[comparisonNum][baseAUDist_2[i].getIndex()][1],
                        baseCoM[comparisonNum][baseAUDist_2[i].getIndex()][2]));
        stringBuilder.append(
                format(" Match Distance:  %d: %8.4f\n", i, sqrt(pairedAUs[i].getDoubleValue())));
        int pairedIndex = pairedAUs[i].getIndex();
        stringBuilder.append(
                format(" Target position: %d: %8.4f %8.4f %8.4f\n", i,
                        targetCoM[pairedIndex][0],
                        targetCoM[pairedIndex][1],
                        targetCoM[pairedIndex][2]));
      }
    }
    if (logger.isLoggable(Level.FINER)) {
      stringBuilder.append("  Distance between pairs:\n" +
              " Index  Base Index  Target Index    Match Index    Distance\n");
      for (int i = 0; i < desiredAUs; i++) {
        stringBuilder.append(format(" %5d %10d %14d %14d %10.4f\n", i,
                baseAUDist_2[i].getIndex(), targetAUDist_2[i].getIndex(), pairedAUs[i].getIndex(),
                sqrt(pairedAUs[i].getDoubleValue())));
      }
    }
  }

  /**
   * Print values for the symmetry operations performed so far (useful for debugging printSym flag).
   *
   * @param compareAtomsSize Number of atoms being compared from each AU.
   * @param file             File to save current coordinates
   * @param name             Description for the current symmetry operator
   * @param bondList         List of bonds for saving in XYZ format
   * @param atoms            List of atoms for saving XYZ format
   * @param forceField       Force field for saving in XYZ format
   * @param save             Save integer switch to determine if/how to save file.
   * @param currZ2           Documentation for which molecule from asymmetric unit is being used
   */
  private void printSym(int compareAtomsSize, File file, String name, List<Bond> bondList, Atom[] atoms,
                        ForceField forceField, int save, int currZ2) {
    // Apply inverse of base operations:
    // For orthogonal matrices the inverse matrix = the transpose. True iff det(A)== +/- 1.
    // No inverse if det(A)==0.
    double[][] tempTransform = new double[4][4];
    copyArrayValues(tempTransform, targetTransformSymOp.asMatrix());
    addTranslation(targetSymOp.tr, tempTransform, true);
    addRotation(targetSymOp.rot, tempTransform, true);
    double[] bestTranslation = new double[]{
            tempTransform[0][3] / tempTransform[3][3],
            tempTransform[1][3] / tempTransform[3][3],
            tempTransform[2][3] / tempTransform[3][3]};
    SymOp symOp = new SymOp(copyOf(tempTransform, 3), bestTranslation);
    double[] newMol = new double[compareAtomsSize * 3];
    StringBuilder printString = new StringBuilder();
    String title = format("\n Print Current Symmetry Operator (Z'=%2d):\n" +
            " \t Original Coords \t\t  Current Coords \t==\t Applied Sym Op Coords", currZ2);
    double[] originalAU = targetAUoriginal[currZ2];
    for (int i = 0; i < compareAtomsSize; i++) {
      int k = i * 3;
      double[] xyz = new double[]{originalAU[k], originalAU[k + 1], originalAU[k + 2]};
      applyCartesianSymOp(xyz, xyz, symOp);
      newMol[k] = xyz[0];
      newMol[k + 1] = xyz[1];
      newMol[k + 2] = xyz[2];
      printString.append(format("\n %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f to %9.3f %9.3f %9.3f",
              originalAU[k], originalAU[k + 1], originalAU[k + 2],
              targetAU[k], targetAU[k + 1], targetAU[k + 2],
              newMol[k], newMol[k + 1], newMol[k + 2]));
    }
    title += format("    RMSD: %9.3f (should be 0.000)", rmsd(newMol, targetAU, massN));
    stringBuilder.append(title).append(printString);
    saveAssembly(file, name, bondList, atoms, forceField, newMol, comparisonAtoms, 1, "_moveMethod",
            0, save);
    stringBuilder.append("\n").append(symOp).append("\n");
  }

  /**
   * Prioritize asymmetric units within the system based on distance to specified index.
   *
   * @param coordsXYZ      Coordinates for expanded crystal (should contain 3 * nAtoms * nMols
   *                       entries).
   * @param mass           Mass of atoms within asymmetric unit (should contain one mass per atom in asym
   *                       unit).
   * @param massSum        Sum of atomic masses within asymmetric unit.
   * @param centerOfMasses Center of masses for each replicate within inflated crystal.
   * @param nAtoms         Number of coordinates in an entity.
   * @param auDist         Prioritization of AUs in expanded system based on linkage criteria
   * @param index          Index of AU to be center.
   * @param linkage        User specified criteria to determine prioritization.
   */
  private static void prioritizeReplicates(double[] coordsXYZ, double[] mass,
                                           double massSum, double[][] centerOfMasses, int nAtoms,
                                           DoubleIndexPair[] auDist, int index, int linkage) {
    // Find AU to be treated as the new center.
    // AUs added to system based on distance to center of all atoms. Index = 0 AU should be closest to all atom center.
    int nCoords = nAtoms * 3;
    int length = coordsXYZ.length;
    int nMols = length / nCoords;
    if(auDist.length != nMols){
      logger.warning(" Number of molecules does not match distance sort array length.");
    }
    if (linkage == 0) {
      // Prioritize based on closest atomic distances.
      int centerIndex = index * nCoords;
      for (int i = 0; i < nMols; i++) {
        double tempDist = Double.MAX_VALUE;
        int molIndex = i * nCoords;
        for (int j = 0; j < nAtoms; j++) {
          int centerAtomIndex = centerIndex + j * 3;
          double[] centerXYZ = {coordsXYZ[centerAtomIndex],
                  coordsXYZ[centerAtomIndex + 1],
                  coordsXYZ[centerAtomIndex + 2]};
          for (int k = 0; k < nAtoms; k++) {
            int atomIndex = molIndex + k * 3;
            double[] xyz = {coordsXYZ[atomIndex],
                    coordsXYZ[atomIndex + 1],
                    coordsXYZ[atomIndex + 2]};
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
        auDist[i] = new DoubleIndexPair(i, tempDist);
      }
      // Sort so the smallest distance is at position 0.
      sort(auDist);

      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, auDist.length);
        double[][] targetMol = new double[length / nCoords][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        logger.finer(" Prioritize replicates:");
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
                  format("  AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i, index,
                          auDist[i].getIndex(),
                          sqrt(auDist[i].getDoubleValue()), targetMol[auDist[i].getIndex()][0],
                          targetMol[auDist[i].getIndex()][1], targetMol[auDist[i].getIndex()][2]));
        }
        logger.finest("\n");
      }
    } else if (linkage == 1) {
      // Prioritize based on geometric center/center of mass
      double[] coordCenter = centerOfMasses[index];
      for (int i = 0; i < nMols; i++) {
        double[] moleculeCenter = centerOfMasses[i];
        auDist[i] = new DoubleIndexPair(i, dist2(coordCenter, moleculeCenter));
      }
      // Reorder based on distance to AU closest to Index.
      sort(auDist);

      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, auDist.length);
        double[][] targetMol = new double[length / nCoords][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          logger.finest(
                  format("  1AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i, index,
                          auDist[i].getIndex(),
                          sqrt(auDist[i].getDoubleValue()), targetMol[auDist[i].getIndex()][0],
                          targetMol[auDist[i].getIndex()][1], targetMol[auDist[i].getIndex()][2]));
        }
      }

      // Molecules in crystal sorted based on distance to center of mass of center most molecule
      // Want the first two molecules chosen in this manner, but third molecule to be closest to both
      // Assigning distance from auDist ensures correct ordering of center most molecule.
      DoubleIndexPair[] molDists2 = new DoubleIndexPair[nMols];
      int auIndex0 = auDist[0].getIndex();
      int auIndex1 = auDist[1].getIndex();
      molDists2[0] = new DoubleIndexPair(auIndex0, auDist[0].getDoubleValue());

      double[] avgCenter = new double[3];
      avgCenter[0] =
              (centerOfMasses[auIndex0][0] + centerOfMasses[auIndex1][0]) / 2;
      avgCenter[1] =
              (centerOfMasses[auIndex0][1] + centerOfMasses[auIndex1][1]) / 2;
      avgCenter[2] =
              (centerOfMasses[auIndex0][2] + centerOfMasses[auIndex1][2]) / 2;

      for (int i = 1; i < nMols; i++) {
        auIndex0 = auDist[i].getIndex();
        double[] moleculeCenter = centerOfMasses[auIndex0];
        molDists2[i] = new DoubleIndexPair(auIndex0, dist2(avgCenter, moleculeCenter));
      }
      sort(molDists2);
      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists2.length);
        double[][] targetMol = new double[length / nCoords][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          logger.finest(format("  1,2AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i,
                          index, molDists2[i].getIndex(),
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
      int bound = min(3, molDists2.length);
      for (int i = 0; i < bound; i++) {
        int molIndex = molDists2[i].getIndex();
        avgCenter[0] += centerOfMasses[molIndex][0];
        avgCenter[1] += centerOfMasses[molIndex][1];
        avgCenter[2] += centerOfMasses[molIndex][2];
      }
      for (int i = 0; i < 3; i++) {
        avgCenter[i] /= 3;
      }

      for (int i = 2; i < nMols; i++) {
        auIndex0 = molDists2[i].getIndex();
        double[] moleculeCenter = centerOfMasses[auIndex0];
        molDists3[i] = new DoubleIndexPair(auIndex0,
                dist2(avgCenter, moleculeCenter));
      }
      //Reorder based on center point between center-most AU to all atom center and closest AU to center-most AU.
      sort(molDists3);
      arraycopy(molDists3, 0, auDist, 0, nMols);
      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, molDists3.length);
        double[][] targetMol = new double[length / nCoords][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          auIndex0 = molDists3[i].getIndex();
          logger.finest(format("  1,2,3AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i,
                          index, auIndex0,
                          sqrt(molDists3[i].getDoubleValue()), targetMol[auIndex0][0],
                          targetMol[auIndex0][1], targetMol[auIndex0][2]));
        }
      }
    } else if (linkage == 2) {
      // Prioritize based on minimum distance between the farthest atom.
      int centerIndex = index * nCoords;
      for (int i = 0; i < nMols; i++) {
        double tempDist = 0.0;
        int molIndex = i * nCoords;
        for (int j = 0; j < nAtoms; j++) {
          int centerAtomIndex = centerIndex + j * 3;
          double[] centerXYZ = {coordsXYZ[centerAtomIndex],
                  coordsXYZ[centerAtomIndex + 1],
                  coordsXYZ[centerAtomIndex + 2]};
          for (int k = 0; k < nAtoms; k++) {
            int atomIndex = molIndex + k * 3;
            double[] xyz = {coordsXYZ[atomIndex], coordsXYZ[atomIndex + 1],
                    coordsXYZ[atomIndex + 2]};
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
        auDist[i] = new DoubleIndexPair(i, tempDist);
      }
      // Sort so the smallest distance is at position 0.
      sort(auDist);

      if (logger.isLoggable(Level.FINEST)) {
        int numCheck = Math.min(5, auDist.length);
        double[][] targetMol = new double[length / nCoords][3];
        centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
        for (int i = 0; i < numCheck; i++) {
          int auIndex = auDist[i].getIndex();
          logger.finest(format("  AU %d Center: %3d Index: %3d Dist %4.4f CoM:(%6.4f %6.4f %6.4f)", i, index,
                          auIndex, sqrt(auDist[i].getDoubleValue()), targetMol[auIndex][0],
                          targetMol[auIndex][1], targetMol[auIndex][2]));
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
    RunningStatistics runningStatistics = distanceMatrixFilter.readDistanceMatrix(filename, baseSize, targetSize);

    if (runningStatistics != null && runningStatistics.getCount() > 0) {
      restartRow = distanceMatrixFilter.getRestartRow();
      restartColumn = distanceMatrixFilter.getRestartColumn();

      if (isSymmetric) {
        // Only the diagonal entry (0.0) is on the last row for a symmetric matrix.
        if (restartRow == baseSize && restartColumn == 1) {
          logger.info(format(" Complete symmetric distance matrix found (%d x %d).",
                  restartRow, restartRow));
        } else {
          restartColumn = 0;
          logger.info(format(" Incomplete symmetric distance matrix found.\n Restarting at row %d, column %d.",
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
   * @param atoms           Atoms we wish to reduce.
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
   * @param file            File to save coordinates
   * @param name            Desired name for file
   * @param bondList        List of bonds for saving in XYZ format
   * @param forceField0     Force field for saving current assembly
   * @param coords          Coordinates to be saved within the PDB
   * @param comparisonAtoms Atoms of interest within the initial asymmetric unit.
   * @param nAU             Number of desired asymmetric units in comparison
   * @param description     Unique identifier that will be added to PDB file name.
   * @param compNum         Unique number for the current comparison
   * @param save            Type of file to save (0=none, 1=PDB, 2=XYZ)
   */
  private void saveAssembly(File file, String name, List<Bond> bondList, Atom[] atoms, ForceField forceField0, final double[] coords,
                            final int[] comparisonAtoms, int nAU, String description, int compNum, int save) {
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
    int nCoords = coords.length / (nAU);
    // Reset atom Index for indexing new atoms
    try { // Assumes same atom selection per AU in Z'>1. Catch exceptions for now.
      for (int n = 0; n < nAU; n++) {
        // Obtain atoms from moved AU (create copy to move to origin)
        // move original and copy to origin
        // rotate original to match copy
        // translate back to moved location
        // Add atoms from moved original to atom list.
        ArrayList<Atom> atomList = new ArrayList<>();
        // Create a new set of Atoms for each SymOp molecule
        int atomValue = 0;
        //Add atoms from comparison to output assembly.
        // Assumes same selection between au when Z'>1...
        // TODO make work for co-crystals.
        int[] atomSelection = new int[nCoords / 3];
        arraycopy(comparisonAtoms, 0, atomSelection, 0, nCoords / 3);
        for (Integer i : atomSelection) {
          Atom a = atoms[i];

          double[] xyz = new double[3];
          int coordIndex = n * nCoords + atomValue;
          xyz[0] = coords[coordIndex];
          xyz[1] = coords[coordIndex + 1];
          xyz[2] = coords[coordIndex + 2];
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
    } catch (Exception exception) {
      stringBuilder.append("\n Error saving moved coordinates to PDB.\n").append(exception).append("\n");
      exception.printStackTrace();
    }

    if (logger.isLoggable(Level.FINEST)) {
      int newSize = newAtomList.size();
      stringBuilder.append(format("\n Save Num AU: %3d", nAU));
      stringBuilder.append(format("\n Save Num Active: %3d", nCoords));
      stringBuilder.append(format("\n Save Num Coords: %3d", coords.length));
      stringBuilder.append(format("\n Save Assembly Length: %3d", newSize));
      if (newSize > 0) {
        Atom zero = newAtomList.get(0);
        stringBuilder.append(format("\n Save Assembly First Atom: %3s %9.3f %9.3f %9.3f\n",
                zero.getName(), zero.getX(), zero.getY(), zero.getZ()));
      } else {
        stringBuilder.append("\n WARNING: New list containing new atoms was empty.\n");
      }
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
          stringBuilder.append("\n Neither key nor properties file detected therefore only creating XYZ.\n");
          if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(ex);
          }
        }
      }
      XYZFilter xyzFilter = new XYZFilter(saveLocation, currentAssembly, forceField, currentAssembly.getProperties());
      xyzFilter.writeFile(saveLocation, true);
    } else {
      PDBFilter pdbfilter = new PDBFilter(saveLocation, currentAssembly, forceField, currentAssembly.getProperties());
      pdbfilter.setModelNumbering(compNum);
      pdbfilter.writeFile(saveLocation, true, false, false);
    }
    currentAssembly.destroy();
  }

  /**
   * Save the provided coordinates as a PDB file with accompanying CSV containing RMSD.
   *
   * @param coords          Coordinates to be saved within the PDB.
   * @param comparisonAtoms Atoms of interest within the initial asymmetric unit.
   * @param description     Unique identifier that will be added to PDB file name.
   * @param finalRMSD       RMSD to be saved to CSV file.
   * @param compNum         Unique number for the current comparison
   * @param save            Type of file to save (0=none, 1=PDB, 2=XYZ)
   */
  private void saveAssembly(File file, String name, List<Bond> bondList, Atom[] atoms, ForceField forceField, double[] coords,
                            int[] comparisonAtoms, int nAU, String description, double finalRMSD, int compNum, int save) {
    saveAssembly(file, name, bondList, atoms, forceField, coords, comparisonAtoms, nAU, description, compNum, save);
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

// Used for LMN expansion. Computes ratio between a, b, c axis lengths to guess final L x M x N distribution.
//  /**
//   * FFX Replicates crystal has to obey lattice requirements. Update LMN accordingly.
//   * @param lattice Lattice system of the current crystal.
//   * @param lmn     Proposed LMN values (will be updated to satisfy lattice requirements).
//   */
//  private static void checkLattice(LatticeSystem lattice, int[] lmn){
//    int restriction;
//    switch(lattice) {
//      case HEXAGONAL_LATTICE:
//      case TETRAGONAL_LATTICE:
//        //a==b
//        restriction = 0;
//        break;
//      case CUBIC_LATTICE:
//      case RHOMBOHEDRAL_LATTICE:
//        //a==b==c
//        restriction = 1;
//        break;
//      default:
//        //no restriction
//        restriction = -1;
//        break;
//    }
//    if(restriction == 0){
//      if(lmn[0] != lmn[1]){
//        int max = max(lmn[0],lmn[1]);
//        lmn[0] = max;
//        lmn[1] = max;
//      }
//    }else if(restriction == 1){
//      if(lmn[0] != lmn[2] || lmn[1] != lmn[2]){
//        int max = max(lmn[0],max(lmn[1],lmn[2]));
//        lmn[0] = max;
//        lmn[1] = max;
//        lmn[2] = max;
//      }
//    }
//  }
//    // Used in LMN expansion as well.
//  /**
//   * Determine replicates expansion based on the size of molecules contained in unit cell.
//   *
//   * @param unitCell        Unit cell to be expanded.
//   * @param reducedCoords           Atoms in system.
//   * @param comparisonAtoms Atoms active in comparison.
//   * @param nAU             Number of asymmetric units to be compared.
//   * @param scaleFactor     Scaling factor to determine final size (default nAU * 10).
//   * @return LMN values for replicates crystal.
//   */
//  private static int[] determineExpansion(final Crystal unitCell, final double[] reducedCoords, final int[] comparisonAtoms, final double[] mass,
//                                          final int nAU, final int zPrime, final double scaleFactor) {
//    // TODO Create better expansion criteria (hopefully can eliminate checkInflatedSphere method some day...).
//    final int XX = 0;
//    final int YY = 1;
//    final int ZZ = 2;
//    // Determine number of symmetry operators in unit cell.
//    final int nSymm = unitCell.spaceGroup.getNumberOfSymOps();
//    final int numEntities = nSymm * zPrime;
//    final int compareAtomsSize = comparisonAtoms.length;
//
//    int zAtoms = compareAtomsSize / zPrime;
//    // Collect asymmetric unit atomic coordinates.
//    double[] x = new double[compareAtomsSize];
//    double[] y = new double[compareAtomsSize];
//    double[] z = new double[compareAtomsSize];
//    for (int i = 0; i < compareAtomsSize; i++) {
//      int atomIndex = i * 3;
//      x[i] = reducedCoords[atomIndex];
//      y[i] = reducedCoords[atomIndex + 1];
//      z[i] = reducedCoords[atomIndex + 2];
//    }
//
//    // Symmetry coordinates for each molecule in replicates crystal
//    double[][] xS = new double[nSymm][compareAtomsSize];
//    double[][] yS = new double[nSymm][compareAtomsSize];
//    double[][] zS = new double[nSymm][compareAtomsSize];
//
//    // Loop over replicate crystal SymOps
//    List<SymOp> inflatedSymOps = unitCell.spaceGroup.symOps;
//    for (int iSym = 0; iSym < nSymm; iSym++) {
//      SymOp symOp = inflatedSymOps.get(iSym);
//      //Convert sym op into cartesian for later use.
//      double[][] rot = new double[3][3];
//      unitCell.getTransformationOperator(symOp, rot);
//      // Apply SymOp to the asymmetric unit reducedCoords Cartesian Coordinates.
//      unitCell.applySymOp(compareAtomsSize, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
//      for (int zp = 0; zp < zPrime; zp++) {
//        int symIndex = zp * zAtoms;
//        // Compute center-of-mass (CoM) for Cartesian coordinates
//        double[] centerOfMass = new double[3];
//        double totalMass = 0.0;
//        for (int i = 0; i < zAtoms; i++) {
//          double m = mass[i];
//          int coordIndex = symIndex + i;
//          centerOfMass[0] += xS[iSym][coordIndex] * m;
//          centerOfMass[1] += yS[iSym][coordIndex] * m;
//          centerOfMass[2] += zS[iSym][coordIndex] * m;
//          totalMass += m;
//        }
//        centerOfMass[0] /= totalMass;
//        centerOfMass[1] /= totalMass;
//        centerOfMass[2] /= totalMass;
//
//        double[] translate = moveIntoCrystal(unitCell, centerOfMass);
//        //Convert sym op into cartesian for later use.
//        double[] trans = symOp.tr.clone();
//        unitCell.toCartesianCoordinates(trans, trans);
//        int translateLength = translate.length;
//        for (int i = 0; i < translateLength; i++) {
//          trans[i] += translate[i];
//        }
//        for (int i = 0; i < zAtoms; i++) {
//          int coordIndex = symIndex + i;
//          xS[iSym][coordIndex] += translate[0];
//          yS[iSym][coordIndex] += translate[1];
//          zS[iSym][coordIndex] += translate[2];
//        }
//      }
//    }
//
//    // Loop over molecule(s) in unit cell to determine maximum dimensions of unit cell.
//    double maxDistX = 0.0;
//    double maxDistY = 0.0;
//    double maxDistZ = 0.0;
//    for (int i = 0; i < nSymm; i++) {
//      double maxSymDistX = 0.0;
//      double maxSymDistY = 0.0;
//      double maxSymDistZ = 0.0;
//      for (int j = 0; j < compareAtomsSize; j++) {
//        double x1 = xS[i][j];
//        double y1 = yS[i][j];
//        double z1 = zS[i][j];
//        for (int m = j + 1; m < compareAtomsSize; m++) {
//          double x2 = xS[i][m];
//          double y2 = yS[i][m];
//          double z2 = zS[i][m];
//          double distX = abs(x2 - x1);
//          double distY = abs(y2 - y1);
//          double distZ = abs(z2 - z1);
//          if (distX > maxSymDistX) {
//            maxSymDistX = distX;
//          }
//          if (distY > maxSymDistY) {
//            maxSymDistY = distY;
//          }
//          if (distZ > maxSymDistZ) {
//            maxSymDistZ = distZ;
//          }
//        }
//      }
//      if (maxSymDistX > maxDistX) {
//        maxDistX = maxSymDistX;
//      }
//      if (maxSymDistY > maxDistY) {
//        maxDistY = maxSymDistY;
//      }
//      if (maxSymDistZ > maxDistZ) {
//        maxDistZ = maxSymDistZ;
//      }
//    }
//    // We have calculated max distance between center of atoms. Add wiggle room to prevent overlap.
//    maxDistX += 4;
//    maxDistY += 4;
//    maxDistZ += 4;
//    // Convert X,Y,Z differences to % of a,b,c lattice lengths.
//    double[] output = new double[3];
//    unitCell.toFractionalCoordinates(new double[]{maxDistX, maxDistY, maxDistZ}, output);
//    // Initial guesses for l, m, n.
//    double target = nAU * scaleFactor;
//    output[XX] = ceil(output[XX]);
//    output[YY] = ceil(output[YY]);
//    output[ZZ] = ceil(output[ZZ]);
//    double l = cbrt((target * output[XX] * output[XX])/(numEntities * output[YY] * output[ZZ]));
//    double m = cbrt((target * output[YY] * output[YY])/(numEntities * output[XX] * output[ZZ]));
//    double n = cbrt((target * output[ZZ] * output[ZZ])/(numEntities * output[XX] * output[YY]));
//    int[] lmn = new int[]{(int) ceil(l), (int) ceil(m), (int) ceil(n)};
//    // Minimum number of asymmetric units desired in replicates crystal.
//    int total = (lmn[0] * lmn[1] * lmn[2]) * numEntities;
//    int[] adjust = new int[]{0, 0, 0};
//    if (logger.isLoggable(Level.FINER)) {
//      logger.info(format(" Number Sym Ops: %3d", nSymm));
//      logger.info(format(" Fractional Distances: %9.3f %9.3f %9.3f", output[XX], output[YY], output[ZZ]));
//      logger.info(format(" LMN Guess: %3d %3d %3d (Total: %3d > %9.3f)", lmn[0], lmn[1], lmn[2], total, target));
//    }
//    // Increase l,m,n while maintaining ratio between them until target number of structures have been obtained.
//    while (total < target) {
//      // Maintain approximate ratio by enforcing each value has to increase before any increase twice.
//      boolean A = adjust[1] == 1 && adjust[2] == 1;
//      boolean B = adjust[0] == 1 && adjust[2] == 1;
//      boolean C = adjust[0] == 1 && adjust[1] == 1;
//      if (A || B || C) {
//        if (A) {
//          adjust[1] = 0;
//          adjust[2] = 0;
//          lmn[0]++;
//        }
//        if (B) {
//          adjust[0] = 0;
//          adjust[2] = 0;
//          lmn[1]++;
//        }
//        if (C) {
//          adjust[0] = 0;
//          adjust[1] = 0;
//          lmn[2]++;
//        }
//      } else {
//        if (adjust[0] == 0 && lmn[0] < lmn[1]) {
//          if (lmn[0] < lmn[2]) {
//            lmn[0]++;
//            adjust[0]++;
//          } else if (lmn[0] == lmn[2]) {
//            if (output[XX] > output[ZZ]) {
//              lmn[0]++;
//              adjust[0]++;
//            } else if (adjust[2] == 0) {
//              lmn[2]++;
//              adjust[2]++;
//            } else {
//              lmn[0]++;
//              adjust[0]++;
//            }
//          } else {
//            lmn[1]++;
//            adjust[1]++;
//          }
//        } else if (adjust[1] == 0 && lmn[1] < lmn[0]) {
//          if (lmn[1] < lmn[2]) {
//            lmn[1]++;
//            adjust[1]++;
//          } else if (lmn[1] == lmn[2]) {
//            if (output[YY] > output[ZZ]) {
//              lmn[1]++;
//              adjust[1]++;
//            } else if (adjust[2] == 0) {
//              lmn[2]++;
//              adjust[2]++;
//            } else {
//              lmn[1]++;
//              adjust[1]++;
//            }
//          } else {
//            lmn[2]++;
//            adjust[2]++;
//          }
//        } else if (adjust[2] == 0 && lmn[2] < lmn[0]) {
//          lmn[2]++;
//          adjust[2]++;
//        } else {
//          if (adjust[0] == 1) {
//            if (output[YY] >= output[ZZ]) {
//              lmn[1]++;
//              adjust[1]++;
//            } else {
//              lmn[2]++;
//              adjust[2]++;
//            }
//          } else if (adjust[1] == 1) {
//            if (output[XX] >= output[ZZ]) {
//              lmn[0]++;
//              adjust[0]++;
//            } else {
//              lmn[2]++;
//              adjust[2]++;
//            }
//          } else if (adjust[2] == 1) {
//            if (output[XX] >= output[YY]) {
//              lmn[0]++;
//              adjust[0]++;
//            } else {
//              lmn[1]++;
//              adjust[1]++;
//            }
//          } else {
//            if (output[XX] >= output[YY] && output[XX] >= output[ZZ]) {
//              lmn[0]++;
//              adjust[0]++;
//            } else if (output[YY] >= output[XX] && output[YY] >= output[ZZ]) {
//              lmn[1]++;
//              adjust[1]++;
//            } else {
//              lmn[2]++;
//              adjust[2]++;
//            }
//          }
//        }
//      }
//      // TODO determine if below steps are necessary (change replicates expansion method)
//      LatticeSystem lattice = unitCell.spaceGroup.latticeSystem;
//      checkLattice(lattice, lmn);
//      total = (lmn[0] * lmn[1] * lmn[2]) * numEntities;
//      if (logger.isLoggable(Level.FINEST)) {
//        logger.finest(format(" \tLMN Guess: %3d %3d %3d (Total: %3d > %9.3f)", lmn[0], lmn[1], lmn[2], total, target));
//      }
//    }
//    // Want a 1 unit cell buffer on all sides. Must be at least 3 wide in each direction.
//    if(lmn[0] < 3){
//      lmn[0] = 3;
//    }
//    if(lmn[1] < 3){
//      lmn[1] = 3;
//    }
//    if(lmn[2] < 3){
//      lmn[2] = 3;
//    }
//    if (logger.isLoggable(Level.FINE)) {
//      logger.fine(format(" Final LMN: %3d %3d %3d (Total: %3d > %9.3f)", lmn[0], lmn[1], lmn[2], total, target));
//    }
//
//    //Return l,m,n for desired replicates crystal.
//    return new int[]{lmn[0], lmn[1], lmn[2]};
//  }

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