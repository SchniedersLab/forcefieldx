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

import static ffx.crystal.Crystal.mod;
import static ffx.numerics.math.DoubleMath.dist;
import static ffx.potential.parsers.DistanceMatrixFilter.writeDistanceMatrixRow;
import static ffx.potential.utils.Superpose.applyRotation;
import static ffx.potential.utils.Superpose.applyTranslation;
import static ffx.potential.utils.Superpose.calculateRotation;
import static ffx.potential.utils.Superpose.calculateTranslation;
import static ffx.potential.utils.Superpose.rmsd;
import static ffx.potential.utils.Superpose.rotate;
import static ffx.potential.utils.Superpose.translate;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cbrt;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.numerics.math.RunningStatistics;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.DistanceMatrixFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.DoubleIndexPair;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * Class ProgressiveAlignmentOfCrystals holds the majority of the functionality necessary to quantify
 * crystal similarity following the PACCOM method.
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
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private final DoubleBuf[] buffers;
  /**
   * Convenience reference to the RMSD value for this process.
   */
  private final double[] myDistance;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private final DoubleBuf myBuffer;
  /**
   * If molecules between two crystals match below this tolerance, it is assumed they are paired.
   */
  private static final double MATCH_TOLERANCE = 1.0E-10;

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

    logger.info(format(" %s conformations: %d", baseLabel, baseSize));
    logger.info(format(" %s conformations: %d", targetLabel, targetSize));

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

    if (numProc > 1) {
      logger.info(format(" Number of MPI Processes:  %d", numProc));
      logger.info(format(" Rank of this MPI Process: %d", rank));
    }

    // Initialize array as -1.0 as -1.0 is not a viable RMSD.
    fill(distRow, -1.0);

    distances = new double[numProc][1];

    // Initialize each distance as -1.0.
    for (int i = 0; i < numProc; i++) {
      fill(distances[i], -1.0);
    }

    // DoubleBuf is a wrapper used by Comm to transfer data between processors.
    buffers = new DoubleBuf[numProc];
    for (int i = 0; i < numProc; i++) {
      buffers[i] = DoubleBuf.buffer(distances[i]);
    }

    // Reference to each processors individual task (for convenience).
    myDistance = distances[rank];
    myBuffer = buffers[rank];
  }

  /**
   * Compare the crystals within the SystemFilters that were inputted into the constructor of this
   * class.
   *
   * @param comparisonAtoms Number of atoms for comparison (cA <= nA).
   * @param nAU Number of asymmetric units to compare.
   * @param inflatedAU Minimum number of asymmetric units in inflated crystal
   * @param numSearch Number of loops to search for mirrors in first system.
   * @param numSearch2 Number of loops to search for mirrors in second system.
   * @param full Perform full number of comparisons (takes longer, but more information
   *     returned).
   * @param save Save out PDBs of the resulting superposition.
   * @param restart Try to restart from a previous job.
   * @param write Save out a PAC RMSD file.
   * @param pacFileName The filename to use.
   * @return Statistics for the computed RMSD values.
   */
  public RunningStatistics comparisons(List<Integer> comparisonAtoms, int nAU, int inflatedAU,
      int numSearch, int numSearch2, int zPrime, boolean full, boolean save, boolean restart,
      boolean write, boolean ares, String pacFileName) {

    RunningStatistics runningStatistics;
    if (restart) {
      runningStatistics = readMatrix(pacFileName, isSymmetric, baseSize, targetSize);
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

    // If shockne group must respect chirality so numSearch = 1
    Crystal baseCrystal = baseFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
    if (!full && baseCrystal.spaceGroup.respectsChirality()) {
      numSearch = 1;
    }
    Crystal targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
    if (!full && targetCrystal.spaceGroup.respectsChirality()) {
      numSearch2 = 1;
    }

    // Minimum amount of time for a single comparison.
    double minTime = Double.MAX_VALUE;

    // restartRow and restartColumn are initialized to zero when this class was constructed.
    // They are updated by the "readMatrix" method if a restart is requested.

    // Read ahead to the base starting conformation.
    for (int row = 0; row < restartRow; row++) {
      baseFilter.readNext(false, false);
    }

    // Padding of the target array size (inner loop limit) is for parallelization.
    // Target conformations are parallelized over available nodes.
    // For example, if numProc = 8 and targetSize = 12, then paddedTargetSize = 16.
    int paddedTargetSize = targetSize;
    int extra = targetSize % numProc;
    if (extra != 0) {
      paddedTargetSize = targetSize - extra + numProc;
      logger.fine(format(" Target size %d vs. Padded size %d", targetSize, paddedTargetSize));
    }

    // Label for logging.
    rmsdLabel = format("RMSD_%d", nAU);

    // Loop over conformations in the base assembly.
    for (int row = restartRow; row < baseSize; row++) {
      // Initialize the distance this rank is responsible for to zero.
      myDistance[0] = -1.0;
      // Base unit cell for logging.
      baseCrystal = baseFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
      for (int column = restartColumn; column < paddedTargetSize; column++) {
        // Make sure this is not a padded value of column.
        if (column < targetSize) {
          int targetRank = column % numProc;
          if (targetRank == rank) {
            long time = -System.nanoTime();
            targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();

            double rmsd;
            if (isSymmetric && row == column) {
              logger.info(format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s",
                  row + 1, baseCrystal.toShortString(), baseLabel,
                  column + 1, targetCrystal.toShortString(), targetLabel));
              // Fill the diagonal.
              rmsd = 0.0;
              // Log the final result.
              logger.info(format(" PAC %s: %12s %7.4f A", rmsdLabel, "", rmsd));
            } else if (isSymmetric && row > column) {
              // Do not compute lower triangle values.
              rmsd = -1.0;
            } else {
              logger.info(format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s",
                  row + 1, baseCrystal.toShortString(), baseLabel,
                  column + 1, targetCrystal.toShortString(), targetLabel));
              // Compute the PAC RMSD.
              rmsd = compare(comparisonAtoms, nAU, inflatedAU, numSearch, numSearch2, zPrime,
                  full, save, ares);
              time += System.nanoTime();
              double timeSec = time * 1.0e-9;
              // Record the fastest comparison.
              if (timeSec < minTime) {
                minTime = timeSec;
              }
              // Log the final result.
              logger.info(format(" PAC %s: %12s %7.4f A (%5.3f sec)", rmsdLabel, "", rmsd, timeSec));
            }
            myDistance[0] = rmsd;

          }
          targetFilter.readNext(false, false);
        }

        // Every numProc iterations, send the results of each rank.
        if ((column + 1) % numProc == 0) {
          gatherRMSDs(row, column, runningStatistics);
        }
      }

      restartColumn = 0;
      targetFilter.readNext(true, false);
      baseFilter.readNext(false, false);

      // Write out this row.
      if (rank == 0 && write) {
        int firstColumn = 0;
        if (isSymmetric) {
          firstColumn = row;
        }
        writeDistanceMatrixRow(pacFileName, distRow, firstColumn);
      }
    }

    if (minTime < Double.MAX_VALUE) {
      logger.info(format("\n Minimum PAC time: %7.4f", minTime));
    }

    baseFilter.closeReader();
    targetFilter.closeReader();

    logger.info(format(" RMSD Minimum:  %8.6f", runningStatistics.getMin()));
    logger.info(format(" RMSD Maximum:  %8.6f", runningStatistics.getMax()));
    logger.info(format(" RMSD Mean:     %8.6f", runningStatistics.getMean()));
    double variance = runningStatistics.getVariance();
    if (!Double.isNaN(variance)) {
      logger.info(format(" RMSD Variance: %8.6f", variance));
    }

    // Return distMatrix for validation if this is for the test script
    return runningStatistics;
  }

  /**
   * Perform single comparison between two crystals.
   *
   * @param comparisonAtoms Number of atoms for comparison (cA <= nA).
   * @param nAU Number of asymmetric units to compare.
   * @param inflatedAU Number of asymmetric units in expanded system.
   * @param numSearch Number of molecules to check for mirror inversions in the first system.
   * @param numSearch2 Number of molecules to check for mirror inversions in the scecond system.
   * @param full Perform full number of comparisons.
   * @param save Save out PDBs of compared crystals.
   * @return the computed RMSD.
   */
  private double compare(List<Integer> comparisonAtoms, int nAU, int inflatedAU,
      int numSearch, int numSearch2, int zPrime, boolean full, boolean save, boolean ares) {
    // TODO: Does PAC work for a combination of molecules and polymers?
    // Number of atoms being used in comparisons
    int compareAtomsSize = comparisonAtoms.size();

    MolecularAssembly baseAssembly = baseFilter.getActiveMolecularSystem();
    baseAssembly.moveAllIntoUnitCell();
    int z;
    if (zPrime > 0) {
      z = zPrime;
    } else if (baseAssembly.getMolecules() != null && baseAssembly.getChainNames() != null) {
      z = baseAssembly.getMolecules().size() + baseAssembly.getChainNames().length;
    } else if (baseAssembly.getChainNames() != null) {
      z = baseAssembly.getChainNames().length;
    } else if (baseAssembly.getMolecules() != null) {
      z = baseAssembly.getMolecules().size();
    } else {
      logger.warning(
          " Number of species in asymmetric unit was not determined. Setting Z'=1. Use --zp flag to set manually.");
      z = 1;
    }
    logger.finer(format(" Number of species in asymmetric unit (Z'): %d", z));

    double[] massStart = new double[compareAtomsSize];
    double[] reducedBaseCoords = reduceSystem(baseAssembly, comparisonAtoms, massStart);
    Crystal baseXtal = baseAssembly.getCrystal();
    double[] baseXYZOrig = generateInflatedSphere(baseXtal, reducedBaseCoords, massStart,
        inflatedAU);

    MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
    targetAssembly.moveAllIntoUnitCell();
    double[] reducedTargetCoords = reduceSystem(targetAssembly, comparisonAtoms, massStart);
    Crystal targetXtal = targetAssembly.getCrystal();
    double[] targetXYZOrig = generateInflatedSphere(targetXtal, reducedTargetCoords, massStart,
        inflatedAU);

    if (z > 1) {
      compareAtomsSize /= z;
    }
    double[] mass = new double[compareAtomsSize];
    arraycopy(massStart, 0, mass, 0, compareAtomsSize);

    // Number of used coordinates for atoms in one AU.
    int nCoords = compareAtomsSize * 3;

    int nBaseMols = baseXYZOrig.length / nCoords;
    int nTargetMols = targetXYZOrig.length / nCoords;
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" Number of copies to compare:    %4d", nAU));
      logger.finer(format(" Number entities in base sphere: %4d", nBaseMols));
      logger.finer(format(" Number entities in target sphere: %d", nTargetMols));
    }

    double[] mass3 = new double[nCoords];
    for (int i = 0; i < 3; i++) {
      arraycopy(mass, 0, mass3, i * compareAtomsSize, compareAtomsSize);
    }

    double[] massN = new double[compareAtomsSize * nAU];
    for (int i = 0; i < nAU; i++) {
      arraycopy(mass, 0, massN, i * compareAtomsSize, compareAtomsSize);
    }

    // Translate asymmetric unit of 0th index (closest to all atom center) to the origin.
    translateAUtoOrigin(baseXYZOrig, mass);

    DoubleIndexPair[] molDist1 = new DoubleIndexPair[nBaseMols];
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(" Prioritize Base System.");
    }

    double massSum = 0;
    for (double value : mass) {
      massSum += value;
    }

    double[][] baseCoMOrig = new double[nBaseMols][3];
    centerOfMass(baseCoMOrig, baseXYZOrig, mass, massSum);
    prioritizeReplicates(compareAtomsSize, baseXYZOrig, mass, massSum, nBaseMols, baseCoMOrig,
        molDist1);

    //Used for debugging. can be removed.
    if (logger.isLoggable(Level.FINER)) {
      int printSize = 20;
      logger.finer(" System 1 distances to center of sphere:");
      for (int i = 0; i < printSize; i++) {
        logger.finer(format(" %d\t%16.8f", molDist1[i].getIndex(), molDist1[i].getDoubleValue()));
      }
    }

    // Translate system to the origin.
    translateAUtoOrigin(targetXYZOrig, mass);

    DoubleIndexPair[] molDist2 = new DoubleIndexPair[nTargetMols];
    // Reorder molDist2 as we shift a different molecule (m) to the center each loop.
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(" Prioritize target system.");
    }
    double[][] targetCoMOrig = new double[nTargetMols][3];
    centerOfMass(targetCoMOrig, targetXYZOrig, mass, massSum);
    prioritizeReplicates(compareAtomsSize, targetXYZOrig, mass, massSum, nTargetMols, targetCoMOrig,
        molDist2);

    if (logger.isLoggable(Level.FINER)) {
      int printSize = 20;
      logger.finer(" System 2 distances to center of sphere:");
      for (int i = 0; i < printSize; i++) {
        logger.finer(format(" %d\t%16.8f", molDist2[i].getIndex(), molDist2[i].getDoubleValue()));
      }
    }

    // TODO determine better numSearch/numSearch2 (exhaustive not excessive?)
    //  Frank-Kasper phases of metallic ions can reach a coordination number of 16...

    // TODO if not using loop for first crystal (l=0 < numSearch):
    //  1. Use molDists1 and molDists12 (not molDists1_2 and molDists12_2).
    //  2. Remove re-prioritization.

    //Determine if AUs in first system are same hand as center most in first.
    logger.finer(" Base Search Handedness:");
    List<Boolean> match1b = new ArrayList<>();
    double[] base0XYZ = new double[nCoords];
    arraycopy(baseXYZOrig, molDist1[0].getIndex() * nCoords, base0XYZ, 0, nCoords);
    int index = 0;
    int match1count = 0;
    List<Integer> baseindices = new ArrayList<>();

    while (Collections.frequency(match1b, true) < numSearch
        || Collections.frequency(match1b, false) < numSearch) {
      double[] baseCheckMol = new double[nCoords];
      arraycopy(baseXYZOrig, molDist1[match1count].getIndex() * nCoords, baseCheckMol, 0, nCoords);
      double value = doMoleculesMatch(base0XYZ, baseCheckMol, mass);
      boolean check = value < MATCH_TOLERANCE;
      if ((Collections.frequency(match1b, true) < numSearch && check) || (
          Collections.frequency(match1b, false) < numSearch && !check)) {
        match1b.add(check);
        logger.finer(format(" %d %4.4f %b", index, value, match1b.get(index)));
        baseindices.add(match1count);
      }
      if (baseXtal.spaceGroup.respectsChirality() && match1b.size() == numSearch) {
        break;
      }
      match1count++;
      if (!(match1count < nBaseMols)) {
        logger.warning(" Only one handedness detected in base crystal.");
        match1b.add(check);
        baseindices.add(match1count);
        break;
      }
    }
    Boolean[] match1 = new Boolean[match1b.size()];
    match1b.toArray(match1);
    Integer[] baseIndices = new Integer[baseindices.size()];
    baseindices.toArray(baseIndices);

    // Determine if AUs in second system are same hand as center most in first
    double highOrLow = -1;
    double high = -1;
    logger.finer(" Target Search Handedness:");
    int matchTcount = 0;
    List<Boolean> match1t = new ArrayList<>();
    List<Integer> targetindices = new ArrayList<>();

    boolean handDone = false;
    while (Collections.frequency(match1t, true) < numSearch2 || (
        Collections.frequency(match1t, false) < numSearch2
            || targetXtal.spaceGroup.respectsChirality())) {
      double[] targetCheckMol = new double[nCoords];
      arraycopy(targetXYZOrig, molDist2[matchTcount].getIndex() * nCoords, targetCheckMol, 0,
          nCoords);
      double value = doMoleculesMatch(base0XYZ, targetCheckMol, mass);
      logger.finer(format(" %d %4.17f", matchTcount, value));
      if (matchTcount == 0) {
        // Starting value
        highOrLow = value;
      }
      boolean check = abs(highOrLow - value) < MATCH_TOLERANCE;
      // value is same as starting value... not sure if High or Low.

      // Second value has been encountered. Determine if High or Low.
      if (handDone) {
        if (value > high) {
          high = value;
        }
        if (high - MATCH_TOLERANCE > value && Collections.frequency(match1t, true) < numSearch2) {
          match1t.add(true);
          targetindices.add(matchTcount);
        } else if (Collections.frequency(match1t, false) < numSearch2) {
          match1t.add(false);
          targetindices.add(matchTcount);
        }
      } else if (!check) {
        handDone = true;
        if (highOrLow > value) {
          for (int i = 0; i < Math.min(matchTcount, numSearch2); i++) {
            match1t.add(false);
            targetindices.add(i);
          }
          match1t.add(true);
          targetindices.add(matchTcount);
          high = highOrLow;
        } else {
          for (int i = 0; i < Math.min(matchTcount, numSearch2); i++) {
            match1t.add(true);
            targetindices.add(i);
          }
          match1t.add(false);
          targetindices.add(matchTcount);
          // Want highOrLow to be the High value in the end.
          high = value;
        }
      }
      if (targetXtal.spaceGroup.respectsChirality()) {
        match1t.add(true);
        targetindices.add(matchTcount);
        if (match1t.size() == numSearch2) {
          break;
        }
      }
      matchTcount++;
      if (!(matchTcount < nTargetMols)) {
        logger.warning(" Only one handedness detected in target crystal.");
        match1t.add(true);
        targetindices.add(matchTcount);
        break;
      }
    }
    Boolean[] matchT1 = new Boolean[match1t.size()];
    match1t.toArray(matchT1);
    Integer[] targetIndices = new Integer[targetindices.size()];
    targetindices.toArray(targetIndices);

    for (int i = 0; i < numSearch2; i++) {
      logger.finer(format(" %d %b", i, matchT1[i]));
    }

    boolean foundt = false;
    boolean foundf = false;
    double[] bestBaseNMols = new double[nAU * nCoords];
    double[] bestTargetNMols = new double[nAU * nCoords];
    double bestRMSD = Double.MAX_VALUE;
    List<Double> crystalCheckRMSDs = new ArrayList<>();

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format("\n Trial      RMSD_1  RMSD_3 %7s", rmsdLabel));
    }

    int baseBooleanIndex = 0;
    for (int l : baseIndices) {
      if (full || (!foundt && match1[baseBooleanIndex] || !foundf && !match1[baseBooleanIndex])) {
        // Place rest of comparison code here.
        double[] baseXYZ = new double[nBaseMols * nCoords];
        double[][] baseCoM = new double[nBaseMols][3];
        arraycopy(baseXYZOrig, 0, baseXYZ, 0, baseXYZOrig.length);
        arraycopy(baseCoMOrig, 0, baseCoM, 0, baseCoMOrig.length);
        int center = molDist1[l].getIndex();
        DoubleIndexPair[] molDist1_2 = new DoubleIndexPair[nBaseMols];
        //Re-prioritize based on center-most molecule.
        prioritizeReplicates(compareAtomsSize, baseXYZ, mass, massSum, nBaseMols, baseCoM,
            molDist1_2, center);

        // Determine densest base3Mols. Use that going forward.
        //Translate base system based on center-most molecule
        double[] baseMol = new double[nCoords];
        int baseAUIndex = molDist1_2[0].getIndex() * nCoords;
        arraycopy(baseXYZ, baseAUIndex, baseMol, 0, nCoords);
        double[] translation = calculateTranslation(baseMol, mass);
        applyTranslation(baseMol, translation);
        applyTranslation(baseXYZ, translation);

        //Update CoMs with translation
        centerOfMass(baseCoM, baseXYZ, mass, massSum);

        // Acquire coordinates based on center 3 molecules
        logger.finer("Base 3 Mol Handedness");
        double[] base3Mols = new double[nCoords * 3];
        for (int i = 0; i < 3; i++) {
          int baseIndex = molDist1_2[i].getIndex() * nCoords;
          arraycopy(baseXYZ, baseIndex, base3Mols, i * nCoords, nCoords);
        }
        logger.finer(" Base 3 Handedness");
        boolean[] match3 = new boolean[3];
        for (int i = 0; i < 3; i++) {
          double[] base1of3Mols = new double[nCoords];
          arraycopy(base3Mols, i * nCoords, base1of3Mols, 0, nCoords);
          double value = doMoleculesMatch(base0XYZ, base1of3Mols, mass);
          match3[i] = value < MATCH_TOLERANCE;
          logger.finer(format(" %d %4.4f %b", i, value, match3[i]));
        }

        reprioritizeReplicates(compareAtomsSize, baseXYZ, mass, massSum, nBaseMols, baseCoM,
            molDist1_2, center);
        // Acquire coordinates for final comparison
        double[] baseNMols = new double[nAU * nCoords];
        for (int i = 0; i < nAU; i++) {
          int molIndex = molDist1_2[i].getIndex() * nCoords;
          arraycopy(baseXYZ, molIndex, baseNMols, i * nCoords, nCoords);
        }
        boolean[] matchN = new boolean[nAU];
        logger.finer(" Base N Handedness");
        for (int i = 0; i < nAU; i++) {
          double[] base1ofNMols = new double[nCoords];
          arraycopy(baseNMols, i * nCoords, base1ofNMols, 0, nCoords);
          double value = doMoleculesMatch(base0XYZ, base1ofNMols, mass);
          matchN[i] = value < MATCH_TOLERANCE;
          logger.finer(format(" %d %4.4f %b", i, value, matchN[i]));
        }
        int targetBooleanIndex = 0;
        for (int m : targetIndices) {
          if (full || match1[baseBooleanIndex] == matchT1[targetBooleanIndex] && (
              !foundt && match1[baseBooleanIndex] || !foundf && !match1[baseBooleanIndex])) {
            double[] targetXYZ = new double[nTargetMols * nCoords];
            double[][] targetCoM = new double[nTargetMols][3];
            arraycopy(targetXYZOrig, 0, targetXYZ, 0, targetXYZOrig.length);
            arraycopy(targetCoMOrig, 0, targetCoM, 0, targetCoMOrig.length);
            if (logger.isLoggable(Level.FINER)) {
              logger.finer("\n");
            }
            // Switch m center most molecules (looking for stereoisomers)
            center = molDist2[m].getIndex();
            DoubleIndexPair[] molDist2_2 = new DoubleIndexPair[nTargetMols];

            //Re-prioritize based on center most molecule.
            if (logger.isLoggable(Level.FINER)) {
              logger.finer(" Re-prioritize target system.");
            }
            prioritizeReplicates(compareAtomsSize, targetXYZ, mass, massSum, nTargetMols, targetCoM,
                molDist2_2, center);

            if (logger.isLoggable(Level.FINER)) {
              logger.finer(" Rotation 1:");
            }
            firstRotation(targetXYZ, mass, baseMol, molDist2_2[0].getIndex());

            double[] targetMol = new double[nCoords];
            int targetCenterMol = molDist2_2[0].getIndex() * nCoords;
            arraycopy(targetXYZ, targetCenterMol, targetMol, 0, nCoords);
            double checkRMSD1 = rmsd(baseMol, targetMol, mass);

            if (logger.isLoggable(Level.FINER)) {
              logger.finer(format(" Center Molecule RMSD after rot 1: %16.8f", checkRMSD1));
            }

            // Check center of mass for center numCheck most entities and distance to center-most.
            if (logger.isLoggable(Level.FINER)) {
              int numCheck = 7;
              double[][] baseCenters = new double[numCheck][3];
              double[][] targetCenters = new double[numCheck][3];
              for (int j = 0; j < numCheck; j++) {
                int baseCenterMols = molDist1_2[j].getIndex() * nCoords;
                int targetCenterMols = molDist2_2[j].getIndex() * nCoords;
                for (int i = 0; i < compareAtomsSize; i++) {
                  int atomIndex = i * 3;
                  int baseIndex = baseCenterMols + atomIndex;
                  int targetIndex = targetCenterMols + atomIndex;
                  double mi = mass[i];
                  baseCenters[j][0] += baseXYZ[baseIndex] * mi;
                  baseCenters[j][1] += baseXYZ[baseIndex + 1] * mi;
                  baseCenters[j][2] += baseXYZ[baseIndex + 2] * mi;
                  targetCenters[j][0] += targetXYZ[targetIndex] * mi;
                  targetCenters[j][1] += targetXYZ[targetIndex + 1] * mi;
                  targetCenters[j][2] += targetXYZ[targetIndex + 2] * mi;
                }
                for (int i = 0; i < 3; i++) {
                  baseCenters[j][i] /= massSum;
                  targetCenters[j][i] /= massSum;
                }
                logger.finer(format(
                    " Position: %d Base(%4.4f): %8.4f %8.4f %8.4f Target(%4.4f): %8.4f %8.4f %8.4f",
                    j,
                    molDist1_2[j].getDoubleValue(),
                    baseCenters[j][0], baseCenters[j][1], baseCenters[j][2],
                    molDist2_2[j].getDoubleValue(),
                    targetCenters[j][0], targetCenters[j][1], targetCenters[j][2]));
              }
            }

            // At this point both systems have completed first rotation/translation
            //  Therefore both center-most molecules should be overlapped.
            // TODO could have prioritization favor closest distance from molDist12 or molDist22.
            if (logger.isLoggable(Level.FINER)) {
              logger.finer(" Match molecules between systems");
            }

            //Update center of masses with the first trans/rot
            centerOfMass(targetCoM, targetXYZ, mass, massSum);

            DoubleIndexPair[] matchMols = new DoubleIndexPair[3];
            matchMolecules(matchMols, baseCoM, targetCoM, molDist1_2, molDist2_2);

            if (logger.isLoggable(Level.FINER)) {
              logger.finer("  Distance between pairs after rot 1:");
              for (DoubleIndexPair matchMol : matchMols) {
                logger.finer(
                    format(" %2d %16.8f", matchMol.getIndex(), matchMol.getDoubleValue()));
              }
              logger.finer(" Index  MolDist1    MolDist2    MatchMols");
              for (int i = 0; i < matchMols.length; i++) {
                logger.finer(format(" %2d base: %2d target: %2d Match: %2d", i,
                    molDist1_2[i].getIndex(), molDist2_2[i].getIndex(), matchMols[i].getIndex()));
              }

              logger.finer(" Rotation 2:");
            }

            secondRotation(base3Mols, targetXYZ, mass3, matchMols);
            double[] target3Mols = new double[nCoords * 3];
            for (int i = 0; i < 3; i++) {
              int molIndex = i * nCoords;
              targetCenterMol = matchMols[i].getIndex() * nCoords;
              arraycopy(targetXYZ, targetCenterMol, target3Mols, molIndex, nCoords);
            }
            double checkRMSD2 = rmsd(base3Mols, target3Mols, mass3);
            logger.finer(" Target 3 Handedness");
            boolean[] matchT3 = new boolean[3];
            boolean hand3Match = true;
            for (int i = 0; i < 3; i++) {
              double[] target1of3Mols = new double[nCoords];
              arraycopy(target3Mols, i * nCoords, target1of3Mols, 0, nCoords);
              double value = doMoleculesMatch(base0XYZ, target1of3Mols, mass);
              matchT3[i] = high - MATCH_TOLERANCE > value;
              logger.finer(format(" %d %4.4f %b", i, value, matchT3[i]));
              if (match3[i] != matchT3[i]) {
                hand3Match = false;
              }
            }

            //Update center of masses with the second rot (only one crystal moves).
            centerOfMass(targetCoM, targetXYZ, mass, massSum);
            reprioritizeReplicates(compareAtomsSize, targetXYZ, mass, massSum, nTargetMols,
                targetCoM,
                molDist2_2, center);
            // Rotations 1 and 2 have been completed and both systems should be overlapped
            //  Isolate center most nAU from System 1 and matching molecules from System 2
            if (logger.isLoggable(Level.FINER)) {
              logger.finer(" Final rotation:");
            }

            if (logger.isLoggable(Level.FINER)) {
              logger.finer(" Match Molecules:");
            }

            matchMols = new DoubleIndexPair[nAU];
            matchMolecules(matchMols, baseCoM, targetCoM, molDist1_2, molDist2_2);

            if (logger.isLoggable(Level.FINER)) {
              int printSize = Math.min(nAU, 10);
              logger.finer(format("  Distance between %d pairs after rot 2:", printSize));
              for (int i = 0; i < printSize; i++) {
                logger.finer(
                    format(" %2d %16.8f", matchMols[i].getIndex(), matchMols[i].getDoubleValue()));
              }
            }

            double[] targetNMols = new double[nAU * nCoords];
            for (int i = 0; i < nAU; i++) {
              int offset = i * nCoords;
              int molIndex = matchMols[i].getIndex() * nCoords;
              arraycopy(targetXYZ, molIndex, targetNMols, offset, nCoords);
            }

            boolean[] matchTN = new boolean[nAU];
            logger.finer(" Target N Handedness");
            boolean handNMatch = true;
            for (int i = 0; i < nAU; i++) {
              double[] target1ofNMols = new double[nCoords];
              arraycopy(targetNMols, i * nCoords, target1ofNMols, 0, nCoords);
              double value = doMoleculesMatch(base0XYZ, target1ofNMols, mass);
              matchTN[i] = high - MATCH_TOLERANCE > value;
              logger.finer(format(" %d %4.4f %b", i, value, matchTN[i]));
              if (matchN[i] != matchTN[i]) {
                handNMatch = false;
              }
            }

            translate(baseNMols, massN, targetNMols, massN);
            rotate(baseNMols, targetNMols, massN);

            // Check center most molecules RMSD between both systems after trans/rot2.
            // Debug check can be removed later.
            target3Mols = new double[nCoords];
            arraycopy(targetNMols, 0, target3Mols, 0, nCoords);
            double checkRMSD3 = rmsd(baseMol, target3Mols, mass);

            if (logger.isLoggable(Level.FINER)) {
              logger.finer(format(" Center Molecule RMSD final: %16.8f", checkRMSD3));
            }

            double rmsdSymOp = rmsd(baseNMols, targetNMols, massN);
            if (logger.isLoggable(Level.FINE)) {
              boolean handednessMatches =
                  match1[baseBooleanIndex] == matchT1[targetBooleanIndex] && hand3Match
                      && handNMatch;
              String output = format(" %2d of %2d: %7.4f %7.4f %7.4f",
                  baseBooleanIndex * numSearch2 + targetBooleanIndex + 1,
                  numSearch * numSearch2, checkRMSD1, checkRMSD2, rmsdSymOp);

              if (logger.isLoggable(Level.FINER)) {
                int countTrue3 = 0;
                for (boolean value : match3) {
                  if (value) {
                    countTrue3++;
                  }
                }
                int countTrueTarget3 = 0;
                for (boolean value : matchT3) {
                  if (value) {
                    countTrueTarget3++;
                  }
                }
                int countTrueN = 0;
                for (boolean value : matchN) {
                  if (value) {
                    countTrueN++;
                  }
                }
                int countTrueTargetN = 0;
                for (boolean value : matchTN) {
                  if (value) {
                    countTrueTargetN++;
                  }
                }
                char b1;
                if (match1[baseBooleanIndex]) {
                  b1 = 't';
                } else {
                  b1 = 'f';
                }
                char t1;
                if (matchT1[targetBooleanIndex]) {
                  t1 = 't';
                } else {
                  t1 = 'f';
                }
                output += format(" %c%c b%dt%d b%dt%d", b1, t1, countTrue3, countTrueTarget3,
                    countTrueN, countTrueTargetN);
                ;
              }

              if (!foundt && match1[baseBooleanIndex] && handednessMatches) {
                foundt = true;
              }
              if (!foundf && !match1[baseBooleanIndex] && handednessMatches) {
                foundf = true;
              }
              logger.fine(output);
            }

//                        if (save) {
//                            int loop = l *numSearch2 + m;
//                            saveAssemblyPDB(baseAssembly, bestBaseNMols, compareAtomsSize, "_c1", loop, false);
//                            saveAssemblyPDB(targetAssembly, bestTargetNMols, compareAtomsSize, "_c2", loop, false);
//                        }

            if (rmsdSymOp < bestRMSD) {
              bestRMSD = rmsdSymOp;
              bestBaseNMols = baseNMols;
              bestTargetNMols = targetNMols;
            }
            addLooseUnequal(crystalCheckRMSDs, rmsdSymOp);
          }
          targetBooleanIndex++;
        }
      }
      baseBooleanIndex++;
    }

    double finalRMSD = Double.NaN;
    if (bestRMSD < Double.MAX_VALUE) {
      finalRMSD = bestRMSD;
    } else {
      logger.warning(" This RMSD was filtered out! Try the -f flag." +
          "\nAlternatively increase --ns and/or --ns2.");
      // TODO: Double.NaN causes an error in RunningStatistics... Set to -2.0 for now...
      finalRMSD = -2.0;
    }
    if (save) {
      saveAssemblyPDB(baseAssembly, bestBaseNMols, comparisonAtoms, "_r" + rank + "_c1",
          0.000, ares, baseXtal.getDensity(massSum));
      saveAssemblyPDB(targetAssembly, bestTargetNMols, comparisonAtoms, "_r" + rank + "_c2",
          finalRMSD, ares, targetXtal.getDensity(massSum));
    }

    // Logging to check number of RMSD values determined.
    StringBuilder dblOut = new StringBuilder();
    Double[] uniqueRMSDs = new Double[crystalCheckRMSDs.size()];
    crystalCheckRMSDs.toArray(uniqueRMSDs);
    Arrays.sort(uniqueRMSDs);
    for (double dbl : uniqueRMSDs) {
      dblOut.append(" ").append(format("%4.4f", dbl));
    }
    String message = format(" Unique %s Values: %s", rmsdLabel, dblOut);

    int numUnique = crystalCheckRMSDs.size();
    if (numUnique > 2) {
      logger.fine(format(
          " PAC determined %2d unique values. Consider increasing the number of inflated molecules.\n %s",
          numUnique, message));
    } else {
      logger.fine(message);
    }

    logger.info(" ");

    return finalRMSD;
  }

  /**
   * Read in the distance matrix.
   *
   * @param filename The PAC RMSD matrix file to read from.
   * @param isSymmetric Is the distance matrix symmetric.
   * @param expectedRows The expected number of rows.
   * @param expectedColumns The expected number of columns.
   * @return Stats for all read in distance matrix values.
   */
  private RunningStatistics readMatrix(String filename, boolean isSymmetric, int expectedRows,
      int expectedColumns) {
    restartRow = 0;
    restartColumn = 0;

    DistanceMatrixFilter distanceMatrixFilter = new DistanceMatrixFilter();
    RunningStatistics runningStatistics = distanceMatrixFilter.readDistanceMatrix(
        filename, expectedRows, expectedColumns);

    if (runningStatistics != null && runningStatistics.getCount() > 0) {
      restartRow = distanceMatrixFilter.getRestartRow();
      restartColumn = distanceMatrixFilter.getRestartColumn();

      if (isSymmetric) {
        // Only the diagonal entry (0.0) is on the last row for a symmetric matrix.
        if (restartRow == expectedRows && restartColumn == 1) {
          logger.info(format(" Complete symmetric distance matrix found (%d x %d).", restartRow,
              restartRow));
        } else {
          restartColumn = 0;
          logger.info(format(
              " Incomplete symmetric distance matrix found.\n Restarting at row %d, column %d.",
              restartRow + 1, restartColumn + 1));
        }
      } else if (restartRow == expectedRows && restartColumn == expectedColumns) {
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
   * This method calls <code>world.AllGather</code> to collect numProc PAC RMSD values.
   *
   * @param row Current row of the PAC RMSD matrix.
   * @param column Current column of the PAC RMSD matrix.
   */
  private void gatherRMSDs(int row, int column, RunningStatistics runningStatistics) {
    if (useMPI) {
      try {
        if (logger.isLoggable(Level.FINER)) {
          logger.finer(" Receiving results.");
        }
        world.allGather(myBuffer, buffers);
        for (int i = 0; i < numProc; i++) {
          int c = (column + 1) - numProc + i;
          if (c < targetSize) {
            distRow[c] = distances[i][0];
            if (!isSymmetric) {
              runningStatistics.addValue(distRow[c]);
            } else if (c >= row) {
              // Only collect stats for the upper triangle.
              runningStatistics.addValue(distRow[c]);
            }
            if (logger.isLoggable(Level.FINER)) {
              logger.finer(format(" %d %d %16.8f", row, c, distances[i][0]));
            }
          }
        }
      } catch (Exception e) {
        logger.severe(" Exception collecting distance values." + e);
      }
    } else {
      distRow[column] = myDistance[0];
      if (!isSymmetric) {
        runningStatistics.addValue(distRow[column]);
      } else if (column >= row) {
        // Only collect stats for the upper triangle.
        runningStatistics.addValue(distRow[column]);
      }
    }
  }

  /**
   * Perform first rotation to match center most molecules.
   *
   * @param targetXYZ Coordinates for system 2
   * @param mass Masses for atoms in a molecule
   * @param index2 Index of center most molecule in coords
   */
  private static void firstRotation(double[] targetXYZ, double[] mass,
      double[] baseMol, int index2) {
    int nAtoms = mass.length;
    int nCoords = nAtoms * 3;

    // Copy base and target coordinates for the center molecule.
    double[] targetMol = new double[nCoords];
    int targetIndex = index2 * nCoords;
    arraycopy(targetXYZ, targetIndex, targetMol, 0, nCoords);

    // Translate the center of mass to the origin.
    double[] translation = calculateTranslation(targetMol, mass);
    applyTranslation(targetMol, translation);
    applyTranslation(targetXYZ, translation);

    // Rotate the target molecule onto the base molecule.
    double[][] rotation = calculateRotation(baseMol, targetMol, mass);
    applyRotation(targetXYZ, rotation);
  }

  /**
   * Translate asymmetric unit to the origin (assumes ordered by distance to center).
   *
   * @param baseXYZ Coordinates of the system
   * @param mass Masses of atoms in one asymmetric unit.
   */
  private static void translateAUtoOrigin(double[] baseXYZ, double[] mass) {
    translateAUtoOrigin(baseXYZ, mass, 0);
  }

  /**
   * Translate an asymmetric unit to the origin.
   *
   * @param baseXYZ Coordinates of the system
   * @param mass Masses of atoms in one asymmetric unit.
   * @param index Index of the asymmetric unit to move.
   */
  private static void translateAUtoOrigin(double[] baseXYZ, double[] mass, int index) {
    int nAtoms = mass.length;
    int nCoords = nAtoms * 3;

    // Load the coordinates.
    double[] baseMol = new double[nCoords];
    int baseIndex = index * nCoords;
    arraycopy(baseXYZ, baseIndex, baseMol, 0, nCoords);

    double[] translateBase = calculateTranslation(baseMol, mass);
    applyTranslation(baseXYZ, translateBase);
  }

  /**
   * Perform second rotation to better match crystal systems.
   *
   * @param base3Mols Coordinates for system 1
   * @param targetXYZ Coordinates for system 2
   * @param mass Masses for atoms in one molecule
   * @param matchMols Indices for system 2 that match molecules in system 1
   */
  private static void secondRotation(double[] base3Mols, double[] targetXYZ, double[] mass,
      DoubleIndexPair[] matchMols) {
    int nCoords = mass.length;

    // Load coordinates for 3 molecules for the base and target systems
    double[] target3Mols = new double[nCoords * 3];
    for (int i = 0; i < 3; i++) {
      int index = i * nCoords;
      int targetIndex = matchMols[i].getIndex() * nCoords;
      arraycopy(targetXYZ, targetIndex, target3Mols, index, nCoords);
    }

    // Calculate the rotation matrix and apply it to the target system.
    applyRotation(targetXYZ, calculateRotation(base3Mols, target3Mols, mass));
  }

  private static void matchMolecules(DoubleIndexPair[] matchMols, double[][] baseCoM,
      double[][] targetCoM,
      DoubleIndexPair[] molDist12, DoubleIndexPair[] molDist22) {
    int desiredMols = matchMols.length;
    int nTargetMols = targetCoM.length;
    // List of indexes for second system.
    int counter = 0;
    List<Integer> targetIndex = new ArrayList<>(nTargetMols);
    for (DoubleIndexPair doubleIndexPair : molDist22) {
      // Only search molecules within range of the desired number of molecules.
      // Must have enough molecules for matching
      targetIndex.add(doubleIndexPair.getIndex());
      counter++;
    }
    logger.finer(format(" Evaluated: %d of %d", counter, nTargetMols));

    // Compare distances between center of masses from system 1 and 2.
    for (int i = 0; i < desiredMols; i++) {
      double minDist = Double.MAX_VALUE;
      Integer minIndex = -1;
      for (Integer target : targetIndex) {
        double dist = dist(baseCoM[molDist12[i].getIndex()], targetCoM[target]);
        if (dist < minDist) {
          minDist = dist;
          minIndex = target;
        }
        if (abs(minDist) < MATCH_TOLERANCE) {
          // Distance between center of masses is ~0 is the best scenario assuming no coordinate overlaps.
          break;
        }
      }
      matchMols[i] = new DoubleIndexPair(minIndex, minDist);
      if (!targetIndex.remove(minIndex)) {
        logger.warning(format(" Index value of %d was not found (%4.4f).", minIndex, minDist));
      }
      if (logger.isLoggable(Level.FINER)) {
        logger.finer(
            format(" Base position:   %d: %8.4f %8.4f %8.4f", i, baseCoM[molDist12[i].getIndex()][0],
                baseCoM[molDist12[i].getIndex()][1],
                baseCoM[molDist12[i].getIndex()][2]));
        logger.finer(format(" Match Distance:  %d: %8.4f", i, matchMols[i].getDoubleValue()));
        logger.finer(format(" Target position: %d: %8.4f %8.4f %8.4f", i,
            targetCoM[molDist22[i].getIndex()][0],
            targetCoM[molDist22[i].getIndex()][1], targetCoM[molDist22[i].getIndex()][2]));
      }
    }
  }

  /**
   * Determine if two molecules are of the same handedness.
   */
  private static double doMoleculesMatch(double[] baseXYZ, double[] targetXYZ,
      double[] mass) {

    translate(baseXYZ, mass, targetXYZ, mass);
    rotate(baseXYZ, targetXYZ, mass);
    return rmsd(baseXYZ, targetXYZ, mass);
  }

  /**
   * Calculate the center of mass for a given set of masses for the asymmetric unit and coordinates
   * (xyz)
   *
   * @param centersOfMass Returned center of mass for each asymmetric unit
   * @param coords Coordinates of every atom in system.
   * @param mass Masses of each atom in asymetric unit.
   */
  private static void centerOfMass(double[][] centersOfMass, double[] coords, double[] mass,
      double massSum) {
    int size = centersOfMass.length;
    int nAtoms = mass.length;
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

  private static double[] reduceSystem(MolecularAssembly assembly, List<Integer> comparisonAtoms,
      double[] mass) {
    Atom[] atoms = assembly.getAtomArray();
    // Collect asymmetric unit atomic coordinates.
    double[] reducedCoords = new double[comparisonAtoms.size() * 3];
    int coordIndex = 0;
    int massIndex = 0;
    for (Integer value : comparisonAtoms) {
      Atom atom = atoms[value];
      double m = atom.getMass();
      mass[massIndex++] = m;
      reducedCoords[coordIndex++] = atom.getX();
      reducedCoords[coordIndex++] = atom.getY();
      reducedCoords[coordIndex++] = atom.getZ();
    }
    return reducedCoords;
  }

  /**
   * Generate and expanded sphere of asymetric unit with the intention of observing a crystals
   * distribution of replicates rather to facilitate comparisons that go beyond lattice parameters.
   *
   * @param crystal Crystal to define expansion.
   * @param inflatedAU Minimum number of replicates desired in resultant sphere.
   * @return double[] containing the coordinates for the expanded crystal.
   */
  private static double[] generateInflatedSphere(Crystal crystal, double[] reducedCoords,
      double[] mass,
      int inflatedAU) {
    int nAtoms = mass.length;
    // Collect asymmetric unit atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];
    double[] xf = new double[nAtoms];
    double[] yf = new double[nAtoms];
    double[] zf = new double[nAtoms];
    for (int i = 0; i < nAtoms; i++) {
      int atomIndex = i * 3;
      x[i] = reducedCoords[atomIndex];
      y[i] = reducedCoords[atomIndex + 1];
      z[i] = reducedCoords[atomIndex + 2];
    }

    // When the system was read in, a replicates crystal may have been created to satisfy the cutoff.
    // Retrieve a reference to the unit cell (not the replicates crystal).
    // Here we will use the unit cell, to create a new replicates crystal that may be
    // a different size (i.e. larger).
    Crystal unitCell = crystal.getUnitCell();

    unitCell.toFractionalCoordinates(nAtoms, x, y, z, xf, yf, zf);

    double asymmetricUnitVolume = unitCell.volume / unitCell.getNumSymOps();

    // Estimate a radius that will include desired number of asymmetric units (inflatedAU).
    double radius = cbrt(inflatedAU * asymmetricUnitVolume) / 2;
    Crystal replicatesCrystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, radius * 2.0);

    // Symmetry coordinates for each molecule in replicates crystal
    int nSymm = replicatesCrystal.getNumSymOps();

    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" Desired copies in target sphere:     %3d", inflatedAU));
      logger.finer(format(" Asymmetric Unit Volume:  %4.2f", asymmetricUnitVolume));
      logger.finer(format(" Estimated spherical radius:  %4.2f", radius));
      logger.finer(" Replicates crystal " + replicatesCrystal);
      logger.finer(format(" Number of replicates: %3d", nSymm));
    }

    double[][] xS = new double[nSymm][nAtoms];
    double[][] yS = new double[nSymm][nAtoms];
    double[][] zS = new double[nSymm][nAtoms];
    // Cartesian center of each molecule
    double[][] centerMolsCart = new double[nSymm][3];

    // Loop over replicate crystal SymOps
    List<SymOp> inflatedSymOps = replicatesCrystal.spaceGroup.symOps;
    for (int iSym = 0; iSym < nSymm; iSym++) {
      SymOp symOp = inflatedSymOps.get(iSym);
      // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
      replicatesCrystal.applySymOp(nAtoms, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
      // Compute center-of-mass (CoM) for Cartesian coordinates
      double[] centerOfMass = new double[3];
      int index = 0;
      double totalMass = 0.0;
      for (double m : mass) {
        centerOfMass[0] += xS[iSym][index] * m;
        centerOfMass[1] += yS[iSym][index] * m;
        centerOfMass[2] += zS[iSym][index++] * m;
        totalMass += m;
      }
      centerOfMass[0] /= totalMass;
      centerOfMass[1] /= totalMass;
      centerOfMass[2] /= totalMass;

      double[] translate = moveIntoCrystal(replicatesCrystal, centerOfMass);
      for (int i = 0; i < nAtoms; i++) {
        xS[iSym][i] += translate[0];
        yS[iSym][i] += translate[1];
        zS[iSym][i] += translate[2];
      }

      // Save CoM cartesian coordinates
      centerMolsCart[iSym] = centerOfMass;
    }

    //Determine molecular distances to "center" of sphere.
    //  In PACCOM the center is the geometric average of coordinates.
    //  In FFX the center is the center of the replicates crystal.
    DoubleIndexPair[] molsDists = new DoubleIndexPair[nSymm];
    double[] cartCenter = new double[3];

    // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
    // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
    double[] fracCenter = {0.5, 0.5, 0.5};
    replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" Expanded Crystal Center: %16.8f %16.8f %16.8f",
          cartCenter[0], cartCenter[1], cartCenter[2]));
    }

    for (int i = 0; i < nSymm; i++) {
      // Then compute Euclidean distance from Cartesian center of the replicates cell
      molsDists[i] = new DoubleIndexPair(i, dist(cartCenter, centerMolsCart[i]));
    }

    // Sort the molecules by their distance from the center.
    // Note that the smallest distances are first in the array after the sort.

    // The problem with parallelSort is that we don't have control over how many threads it
    // will use. If we MPI parallelize PAC over many processes on one node,
    // the parallelSort of each process cannot try to use all threads of the node.
    // Arrays.parallelSort(molsDists);
    Arrays.sort(molsDists);

    if (logger.isLoggable(Level.FINEST)) {
      logger.finest("\n Copy  SymOp        Distance");
    }
    double[] systemCoords = new double[nSymm * nAtoms * 3];
    for (int n = 0; n < nSymm; n++) {
      // Current molecule
      int iSym = molsDists[n].getIndex();
      double distance = molsDists[n].getDoubleValue();
      if (logger.isLoggable(Level.FINEST)) {
        logger.finest(format(" %4d  %5d  %16.8f", n, iSym, distance));
      }

      // Create a new set of Atoms for each SymOp molecule
      for (int i = 0; i < nAtoms; i++) {
        int symIndex = n * nAtoms * 3;
        int atomIndex = i * 3;
        systemCoords[symIndex + atomIndex] = xS[iSym][i];
        systemCoords[symIndex + atomIndex + 1] = yS[iSym][i];
        systemCoords[symIndex + atomIndex + 2] = zS[iSym][i];
      }
    }

    return systemCoords;
  }

  /**
   * Produce a translation vector necessary to move an object with the current center of mass (com)
   * into the provided crystal.
   *
   * @param crystal Replicates crystal within whom coordinates should be moved.
   * @param com Center of mass (x, y, z) for the object of concern
   * @return double[] translation vector to move the object to within the provided crystal.
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

    if (logger.isLoggable(Level.FINEST)) {
      logger.finest(format(" Center of Mass Prior: %16.8f %16.8f %16.8f",
          currentCoM[0], currentCoM[1], currentCoM[2]));
      logger.finest(format(" Center of Mass Post: %16.8f %16.8f %16.8f", com[0], com[1], com[2]));
    }

    return translate;
  }

  /**
   * Prioritize asymmetric units within the system based on distance to center-most (0th index).
   *
   * @param nAtoms Number of atoms being compared in asymmetric unit.
   * @param coordsXYZ Coordinates for expanded crystal.
   * @param mass Mass of atoms within asymmetric unit.
   * @param nMols Number of asymmetric units in expanded crystal.
   * @param molDists DoubleIndexPair of distances from center most asymmetric unit to asymmetric
   *     unit at index.
   */
  private static void prioritizeReplicates(int nAtoms, double[] coordsXYZ, double[] mass,
      double massSum, int nMols,
      double[][] centerOfMasses, DoubleIndexPair[] molDists) {
    // AUs added to system based on distance to center of all atoms. 0th AU should be closest to all atom center.
    prioritizeReplicates(nAtoms, coordsXYZ, mass, massSum, nMols, centerOfMasses, molDists, 0);
  }

  /**
   * Prioritize asymmetric units within the system based on distance to specified index.
   *
   * @param nAtoms Number of atoms being compared in asymmetric unit.
   * @param coordsXYZ Coordinates for expanded crystal.
   * @param mass Mass of atoms within asymmetric unit.
   * @param nMols Number of molecules in expanded crystal.
   * @param molDists Prioritization of molecules in expanded system on distance to center-most
   *     molecule.
   * @param index Index of molecules to be center.
   */
  private static void prioritizeReplicates(int nAtoms, double[] coordsXYZ, double[] mass,
      double massSum, int nMols,
      double[][] centerOfMasses, DoubleIndexPair[] molDists,
      int index) {
    // Find AU to be treated as the new center.
    double[] coordCenter = centerOfMasses[index];
    for (int i = 0; i < nMols; i++) {
      double[] moleculeCenter = centerOfMasses[i];
      molDists[i] = new DoubleIndexPair(i, dist(coordCenter, moleculeCenter));
    }
    // Reorder based on distance to AU closest to Index.
    Arrays.sort(molDists);

    if (logger.isLoggable(Level.FINER)) {
      int numCheck = Math.min(7, molDists.length);
      double[][] targetMol = new double[coordsXYZ.length / (nAtoms * 3)][3];
      centerOfMass(targetMol, coordsXYZ, mass, massSum);
      for (int i = 0; i < numCheck; i++) {
        logger.finer(format("Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
            molDists[i].getIndex(),
            molDists[i].getDoubleValue(), targetMol[molDists[i].getIndex()][0],
            targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
      }
    }

//        if(lineCheck) {
//            //Determine if the three points are in a line.
//            double xcoord = 0;
//            double ycoord = 0;
//            double zcoord = 0;
//            int counter = 3;
//            while (xcoord == ycoord && xcoord == zcoord) {
//                xcoord = centerOfMasses[molDists[2].getIndex()][0] - centerOfMasses[molDists[0].getIndex()][0] /
//                        centerOfMasses[molDists[1].getIndex()][0] - centerOfMasses[molDists[0].getIndex()][0];
//                ycoord = centerOfMasses[molDists[2].getIndex()][1] - centerOfMasses[molDists[0].getIndex()][1] /
//                        centerOfMasses[molDists[1].getIndex()][1] - centerOfMasses[molDists[0].getIndex()][1];
//                zcoord = centerOfMasses[molDists[2].getIndex()][2] - centerOfMasses[molDists[0].getIndex()][2] /
//                        centerOfMasses[molDists[1].getIndex()][2] - centerOfMasses[molDists[0].getIndex()][2];
//                if (xcoord == ycoord && xcoord == zcoord) {
//                    if (logger.isLoggable(Level.FINER)) {
//                        logger.warning(" Closest 3 molecules in line.");
//                    }
//                    DoubleIndexPair temp = molDists[2];
//                    molDists[2] = molDists[counter];
//                    molDists[counter++] = temp;
//                }
//            }
//        }else {
    //First attempt to avoid CoMs from being in a straight line.
    // Molecules in crystal sorted based on distance to center of mass of center most molecule
    // Want the first two molecules chosen in this manner, but third molecule to be closest to both
    // Assigning distance from molDists ensures correct ordering of center most molecule.
    DoubleIndexPair[] molDists2 = new DoubleIndexPair[nMols];
    molDists2[0] = new DoubleIndexPair(molDists[0].getIndex(), molDists[0].getDoubleValue());

    double[] avgCenter = new double[3];
    avgCenter[0] =
        (centerOfMasses[molDists[0].getIndex()][0] + centerOfMasses[molDists[1].getIndex()][0]) / 2;
    avgCenter[1] =
        (centerOfMasses[molDists[0].getIndex()][1] + centerOfMasses[molDists[1].getIndex()][1]) / 2;
    avgCenter[2] =
        (centerOfMasses[molDists[0].getIndex()][2] + centerOfMasses[molDists[1].getIndex()][2]) / 2;

    for (int i = 1; i < nMols; i++) {
      double[] moleculeCenter = centerOfMasses[molDists[i].getIndex()];
      molDists2[i] = new DoubleIndexPair(molDists[i].getIndex(), dist(avgCenter, moleculeCenter));
    }
    //Reorder based on center point between center-most AU to all atom center and closest AU to center-most AU.
    Arrays.sort(molDists2);
    arraycopy(molDists2, 0, molDists, 0, nMols);
//        }
  }

  /**
   * Reprioritize asymmetric units within the system based on distance to specified index.
   *
   * @param nAtoms Number of atoms being compared in asymmetric unit.
   * @param coordsXYZ Coordinates for expanded crystal.
   * @param mass Mass of atoms within asymmetric unit.
   * @param nMols Number of molecules in expanded crystal.
   * @param molDists Prioritization of molecules in expanded system on distance to center-most
   *     molecule.
   * @param index Index of molecules to be center.
   */
  private static void reprioritizeReplicates(int nAtoms, double[] coordsXYZ, double[] mass,
      double massSum, int nMols,
      double[][] centerOfMasses, DoubleIndexPair[] molDists,
      int index) {
    //First attempt to avoid CoMs from being in a straight line.
    // Molecules in crystal sorted based on distance to center of mass of center most molecule
    // Want the first two molecules chosen in this manner, but third molecule to be closest to both
    // Assigning distance from molDists ensures correct ordering of center most molecule.
    DoubleIndexPair[] molDists2 = new DoubleIndexPair[nMols];
    molDists2[0] = new DoubleIndexPair(molDists[0].getIndex(), molDists[0].getDoubleValue());
    molDists2[1] = new DoubleIndexPair(molDists[1].getIndex(), molDists[1].getDoubleValue());
    double[] avgCenter = new double[3];
    for (int i = 0; i < 3; i++) {
      avgCenter[0] += centerOfMasses[molDists[i].getIndex()][0];
      avgCenter[1] += centerOfMasses[molDists[i].getIndex()][1];
      avgCenter[2] += centerOfMasses[molDists[i].getIndex()][2];
    }
    for (int i = 0; i < 3; i++) {
      avgCenter[i] /= 3;
    }

    for (int i = 2; i < nMols; i++) {
      double[] moleculeCenter = centerOfMasses[molDists[i].getIndex()];
      molDists2[i] = new DoubleIndexPair(molDists[i].getIndex(), dist(avgCenter, moleculeCenter));
    }
    //Reorder based on center point between center-most AU to all atom center and closest AU to center-most AU.
    Arrays.sort(molDists2);
    arraycopy(molDists2, 0, molDists, 0, nMols);

    if (logger.isLoggable(Level.FINER)) {
      int numCheck = 7;
      double[][] targetMol = new double[coordsXYZ.length / (nAtoms * 3)][3];
      centerOfMass(targetMol, coordsXYZ, mass, massSum);
      for (int i = 0; i < numCheck; i++) {
        logger.finer(format("Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
            molDists[i].getIndex(),
            molDists[i].getDoubleValue(), targetMol[molDists[i].getIndex()][0],
            targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
      }
    }
  }

  /**
   * Add a value to a list of long if its difference to all listed values is greater than the
   * tolerance.
   *
   * @param values List of values already found.
   * @param value Potential new value if it is not already in list.
   */
  private static void addLooseUnequal(List<Double> values, double value) {
    boolean found = false;
    for (Double dbl : values) {
      if (abs(dbl - value) < MATCH_TOLERANCE) {
        found = true;
        break;
      }
    }
    if (!found) {
      values.add(value);
    }
  }

  /**
   * Save the current assembly as a PDB file.
   */
  private static void saveAssemblyPDB(MolecularAssembly molecularAssembly, double[] coords,
      List<Integer> comparisonAtoms, String description,
      double finalRMSD, boolean ares, double density) {
    String fileName = FilenameUtils.removeExtension(molecularAssembly.getName());
    File saveLocationPDB = new File(fileName + description + ".pdb");
    int fileNumber = 0;
    while (saveLocationPDB.exists()) {
      saveLocationPDB = new File(fileName + description + "_" + fileNumber++ + ".pdb");
    }
    // Save aperiodic system of n_mol closest atoms for visualization
    MolecularAssembly currentAssembly = new MolecularAssembly(molecularAssembly.getName());
//        List<Bond> bondList = molecularAssembly.getBondList();
    ArrayList<Atom> newAtomList = new ArrayList<>();
    Atom[] atoms = molecularAssembly.getAtomArray();
    int atomIndex = 0;
    int compareAtomsSize = comparisonAtoms.size();
    int numMols = coords.length / (3 * compareAtomsSize);
    for (int n = 0; n < numMols; n++) {
      ArrayList<Atom> atomList = new ArrayList<>();
      // Create a new set of Atoms for each SymOp molecule
      int atomValue = 0;
      for (Integer i : comparisonAtoms) {
        Atom a = atoms[i];
        double[] xyz = new double[3];
        xyz[0] = coords[n * compareAtomsSize * 3 + atomValue];
        xyz[1] = coords[n * compareAtomsSize * 3 + atomValue + 1];
        xyz[2] = coords[n * compareAtomsSize * 3 + atomValue + 2];
        Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
        atomList.add(atom);
        atomValue += 3;
      }
      // Create a new set of Bonds for each SymOp molecule
//            for (Bond bond : bondList) {
//                Atom a1 = bond.getAtom(0);
//                Atom a2 = bond.getAtom(1);
//                Atom newA1 = atomList.get(a1.getIndex() - 1);
//                Atom newA2 = atomList.get(a2.getIndex() - 1);
//                Bond b = new Bond(newA1, newA2);
//                b.setBondType(bond.getBondType());
//            }
      newAtomList.addAll(atomList);
    }

    // Construct the force field for the expanded set of molecules
    ForceField forceField = molecularAssembly.getForceField();

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

    currentAssembly.setFile(molecularAssembly.getFile());
    PDBFilter pdbfilter = new PDBFilter(saveLocationPDB, currentAssembly, forceField,
        currentAssembly.getProperties());
    pdbfilter.writeFile(saveLocationPDB, false);
    if (ares) {
      try {
        BufferedWriter bw = new BufferedWriter(new FileWriter(saveLocationPDB, true));
        bw.append(format("rms %4.4f\n", finalRMSD));
        //May not be needed. Energy would likely be a better value, but might take more time...
        bw.append(format("density %4.4f\n", density));
        bw.append(format("energy %4.4f\n", molecularAssembly.getPotentialEnergy().energy()));
        bw.close();
      } catch (Exception ex) {
        logger.info(ex.toString());
      }
    }
  }

//    /**
//     * Orient coordinates so that the second index is on the x axis, and the third index is on the X-Y
//     * plane. First index should be at the origin (0, 0, 0).
//     *
//     * @param coordsXYZ   An array of XYZ positions (e.g. [x0, y0, z0, x1, y1, z1, x2, y2, z2]
//     * @param atomIndices Indices for three desired sets from the XYZ list (e.g. [0, 1, 2]). Index
//     *                    0 should be at origin.
//     */
//    public static void standardOrientation(double[] coordsXYZ, int[] atomIndices) {
//        // Used in QEtoXYZ.groovy which is not ready for git which is why this method appears unused.
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
}
