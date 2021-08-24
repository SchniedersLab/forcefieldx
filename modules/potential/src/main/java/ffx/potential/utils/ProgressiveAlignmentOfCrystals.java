// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import static ffx.potential.parsers.DistanceMatrixFilter.toDistanceMatrixString;
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
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.DistanceMatrixFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.DoubleIndexPair;
import java.io.File;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FilenameUtils;

/**
 * Class ProgressiveAlignmentOfCrystals holds the majority of the functionality necessary to quantify
 * crystal similarity following the PACCOM method.
 *
 * @author Aaron J. Nessler, Okimasa OKADA and Michael J. Schnieders
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
   * Label for the second crystal.
   */
  private final String targetLabel;
  /** Label to use for the RSMD logging */
  private String rmsdLabel;
  /**
   * Matrix of RMSD values [baseSize][targetSize].
   */
  public final double[][] distMatrix;
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
   * Number of processes.
   */
  private final int numProc;
  /**
   * Rank of this process.
   */
  private final int rank;
  /**
   * The distances matrix stores a single RSMD value from each process.
   * The array is of size [numProc][1].
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
  private static final double MATCH_TOLERANCE = 1.0E-4;

  /**
   * Constructor for the ProgressiveAlignmentOfCrystals class.
   *
   * @param baseFilter SystemFilter containing a set of crystal structures to compare.
   * @param targetFilter SystemFilter containing the other set of crystals to compare.
   */
  public ProgressiveAlignmentOfCrystals(SystemFilter baseFilter, SystemFilter targetFilter) {
    this.baseFilter = baseFilter;
    this.targetFilter = targetFilter;

    // Number of models to be evaluated.
    baseSize = baseFilter.countNumModels();
    baseLabel = getName(baseFilter.getFile().getAbsolutePath());
    targetSize = targetFilter.countNumModels();
    targetLabel = getName(targetFilter.getFile().getAbsolutePath());

    logger.info(format(" %s conformations: %d", baseLabel, baseSize));
    logger.info(format(" %s conformations: %d", targetLabel, targetSize));

    world = Comm.world();
    // Number of processes is equal to world size (often called size).
    numProc = world.size();
    // Each processor gets its own rank (ID of sorts).
    rank = world.rank();
    if (numProc > 1) {
      logger.info(format(" Number of MPI Processes:  %d", numProc));
      logger.info(format(" Rank of this MPI Process: %d", rank));
    }

    // Distance matrix to store compared values (dimensions are "human readable" [m x n]).
    distMatrix = new double[baseSize][targetSize];
    // Initialize array as -1.0 as -1.0 is not a viable RMSD.
    for (int i = 0; i < baseSize; i++) {
      fill(distMatrix[i], -1.0);
    }

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
   * @param nAtoms Number of atoms in asymmetric unit.
   * @param comparisonAtoms Number of atoms for comparison (cA <= nA).
   * @param nAU Number of asymmetric units to compare.
   * @param inflatedAU Minimum number of asymmetric units in inflated crystal
   * @param numSearch Number of loops to search for mirrors in first system.
   * @param numSearch2 Number of loops to search for mirrors in second system.
   * @param force Compare regardless of single molecule RMSD.
   * @param symmetric Perform 2 sided calculation (A vs B; B vs A)
   * @param save Save out PDBs of the resulting superposition.
   * @param restart Try to restart from a previous job.
   * @param write Save out a PAC RMSD file.
   * @param pacFileName The filename to use.
   */
  public double[][] comparisons(int nAtoms, List<Integer> comparisonAtoms, int nAU, int inflatedAU,
      int numSearch, int numSearch2, boolean force, boolean symmetric, boolean save, boolean restart,
      boolean write, String pacFileName) {

    if (restart) {
      readMatrix(pacFileName);
    } else {
      File file = new File(pacFileName);
      if (file.exists() && file.delete()) {
        logger.info(format(" PAC RMSD file (%s) was deleted.", pacFileName));
        logger.info(" To restart from a previous run, use the '-r' flag.");
      }
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
      Crystal baseCrystal = baseFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
      for (int column = restartColumn; column < paddedTargetSize; column++) {
        // Make sure this is not a padded value of column.
        if (column < targetSize) {
          int targetRank = column % numProc;
          if (targetRank == rank) {
            long time = -System.nanoTime();
            Crystal targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
            logger.info(format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s",
                    row + 1, baseCrystal.toShortString(), baseLabel,
                    column + 1, targetCrystal.toShortString(), targetLabel));

            double rmsd;
            if (row == column) {
              // Fill the diagonal.
              rmsd = 0.0;
              // Log the final result.
              logger.info(format(" PAC %s: %12s %7.4f A", rmsdLabel, "", rmsd));
            } else if (row >= column) {
              // Fill the lower triangle from the upper triangle.
              rmsd = distMatrix[column][row];
              // Log the final result.
              logger.info(format(" PAC %s: %12s %7.4f A", rmsdLabel, "", rmsd));
            } else {
              // Compute the PAC RMSD.
              rmsd = compare(nAtoms, comparisonAtoms, nAU, inflatedAU, numSearch, numSearch2,
                  force, symmetric, save);
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
          gatherRMSDs(row, column);
        }
      }

      restartColumn = 0;
      targetFilter.readNext(true, false);
      baseFilter.readNext(false, false);

      // Write out this row.
      if (rank == 0 && write) {
        writeDistanceMatrixRow(pacFileName, distMatrix[row]);
      }
    }

    if (minTime < Double.MAX_VALUE) {
      logger.info(format("\n Minimum PAC time: %7.4f", minTime));
    }

    baseFilter.closeReader();
    targetFilter.closeReader();

    // Print the PAC RMSD matrix.
    if (baseSize > 1 || targetSize > 1) {
      logger.info("\n" + toDistanceMatrixString(distMatrix));
    }

    // Return distMatrix for validation if this is for the test script
    return distMatrix;
  }

  /**
   * Perform single comparison between two crystals.
   *
   * @param nAtoms Number of atoms in asymmetric unit.
   * @param comparisonAtoms Number of atoms for comparison (cA <= nA).
   * @param nAU Number of asymmetric units to compare.
   * @param inflatedAU Number of asymmetric units in expanded system.
   * @param numSearch Number of molecules to check for mirror inversions in the first system.
   * @param numSearch2 Number of molecules to check for mirror inversions in the scecond system.
   * @param force Compare regardless of single molecule RMSD.
   * @param symmetric Perform a symmetric PAC calculation (A vs B; B vs A).
   * @param save Save out PDBs of compared crystals.
   * @return the computed RMSD.
   */
  private double compare(int nAtoms, List<Integer> comparisonAtoms, int nAU, int inflatedAU,
      int numSearch, int numSearch2, boolean force, boolean symmetric, boolean save) {

    logger.finer(format(" Number of copies to compare: %4d", nAU));
    MolecularAssembly baseAssembly = baseFilter.getActiveMolecularSystem();
    MolecularAssembly baseSphere = generateInflatedSphere(baseAssembly, inflatedAU);

    int nBaseMols;
    if (baseSphere.getMolecules() != null && baseSphere.getChainNames() != null) {
      nBaseMols = baseSphere.getMolecules().size() + baseSphere.getChainNames().length;
    } else if (baseSphere.getChainNames() != null) {
      nBaseMols = baseSphere.getChainNames().length;
    } else if (baseSphere.getMolecules() != null) {
      nBaseMols = baseSphere.getMolecules().size();
    } else {
      logger.severe(" Could not find crystal entities for base crystal.");
      return Double.NaN;
    }

    logger.finer(format(" Number entities in base sphere: %3d", nBaseMols));
    // Save coordinates from expanded sphere and save out central molecule coordinates.
    Atom[] baseAtoms = baseSphere.getAtomArray();

    int compareAtomsSize = comparisonAtoms.size();
    double[] baseXYZ = new double[compareAtomsSize * nBaseMols * 3];
    double[] mass = new double[compareAtomsSize];
    double massSum = 0;
    int index = 0;
    int massIndex = 0;
    for (int i = 0; i < baseAtoms.length; i++) {
      if (comparisonAtoms.contains(i % nAtoms)) {
        Atom atom = baseAtoms[i];
        baseXYZ[index++] = atom.getX();
        baseXYZ[index++] = atom.getY();
        baseXYZ[index++] = atom.getZ();
        if (i < nAtoms) {
          mass[massIndex++] = atom.getMass();
          massSum += atom.getMass();
        }
      }
    }

    // Translate asymmetric unit of 0th index (closest to center) to the origin.
    translateAUtoOrigin(baseXYZ, mass);

    DoubleIndexPair[] molDist1 = new DoubleIndexPair[nBaseMols];
    DoubleIndexPair[] molDist12 = new DoubleIndexPair[nBaseMols];
    logger.finer(" Prioritize Base System.");
    prioritizeReplicates(compareAtomsSize, baseXYZ, mass, massSum, nBaseMols, molDist1, molDist12);

    //Used for debugging.. can be removed.
    int printSize = 20;
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(" System 1 distances to center of sphere:");
      for (int i = 0; i < printSize; i++) {
        logger.finer(format(" %d\t%16.8f", molDist12[i].getIndex(), molDist12[i].getDoubleValue()));
      }
    }

    MolecularAssembly molecularAssembly = targetFilter.getActiveMolecularSystem();

    MolecularAssembly targetSphere = generateInflatedSphere(molecularAssembly, inflatedAU);

    logger.fine(format("\n Trial      RMSD_1  RMSD_3 %7s", rmsdLabel));
    int nTargetMols;
    if (targetSphere.getMolecules() != null && targetSphere.getChainNames() != null) {
      nTargetMols = targetSphere.getMolecules().size() + targetSphere.getChainNames().length;
    } else if (targetSphere.getChainNames() != null) {
      nTargetMols = targetSphere.getChainNames().length;
    } else if (targetSphere.getMolecules() != null) {
      nTargetMols = targetSphere.getMolecules().size();
    } else {
      logger.severe(" Could not find crystal entities for target sphere.");
      return Double.NaN;
    }

    logger.finer(format(" Number entities in target sphere: %d", nTargetMols));
    // Save coordinates from expanded sphere and save out central molecule coordinates
    Atom[] targetAtoms = targetSphere.getAtomArray();
    double[] targetXYZ = new double[compareAtomsSize * nTargetMols * 3];

    //Assumes masses are same as first system.
    index = 0;
    for (int i = 0; i < targetAtoms.length; i++) {
      if (comparisonAtoms.contains(i % nAtoms)) {
        Atom atom = targetAtoms[i];
        targetXYZ[index++] = atom.getX();
        targetXYZ[index++] = atom.getY();
        targetXYZ[index++] = atom.getZ();
      }
    }

    //Translate system to the origin.
    translateAUtoOrigin(targetXYZ, mass);

    DoubleIndexPair[] molDist2 = new DoubleIndexPair[nTargetMols];
    DoubleIndexPair[] molDist22 = new DoubleIndexPair[nTargetMols];
    // Reorder molDist2 as we shift a different molecule (m) to the center each loop.
    logger.finer(" Prioritize target system.");
    prioritizeReplicates(compareAtomsSize, targetXYZ, mass, massSum, nTargetMols, molDist2,
        molDist22);

    if (logger.isLoggable(Level.FINER)) {
      logger.finer(" System 2 distances to center of sphere:");
      for (int i = 0; i < printSize; i++) {
        logger.finer(format(" %d\t%16.8f", molDist22[i].getIndex(), molDist22[i].getDoubleValue()));
      }
    }

    double bestRMSD = Double.MAX_VALUE;
    List<Double> crystalCheckRMSDs = new ArrayList<>();

    // TODO determine better numSearch/numSearch2 (exhaustive not excessive?)
    //  Frank-Kasper phases of metallic ions can reach a coordination number of 16...

    // TODO if not using loop for first crystal (l=0 < numSearch):
    //  1. Use molDists1 and molDists12 (not molDists1_2 and molDists12_2).
    //  2. Remove re-prioritization.
    for (int l = 0; l < numSearch; l++) {
      int center = molDist12[l].getIndex();
      DoubleIndexPair[] molDist1_2 = new DoubleIndexPair[nBaseMols];
      DoubleIndexPair[] molDist12_2 = new DoubleIndexPair[nBaseMols];
      //Re-prioritize based on center most molecule.
      logger.finer(" Re-prioritize base system.");
      prioritizeReplicates(compareAtomsSize, baseXYZ, mass, massSum, nBaseMols, molDist1_2, center,
          molDist12_2);

      for (int m = 0; m < numSearch2; m++) {
        logger.finer("\n");
        // Switch m center most molecules (looking for chirality)
        center = molDist22[m].getIndex();
        DoubleIndexPair[] molDist2_2 = new DoubleIndexPair[nTargetMols];
        DoubleIndexPair[] molDist22_2 = new DoubleIndexPair[nTargetMols];

        //Re-prioritize based on center most molecule.
        logger.finer(" Re-prioritize target system.");
        prioritizeReplicates(compareAtomsSize, targetXYZ, mass, massSum, nTargetMols, molDist2_2,
            center, molDist22_2);

        // Check center most molecules RMSD between both systems.
        double[] baseMolCoords = new double[compareAtomsSize * 3];
        double[] targetMolCoords = new double[compareAtomsSize * 3];
        int baseCenterMol = molDist12_2[0].getIndex() * 3 * compareAtomsSize;
        int targetCenterMol = molDist22_2[0].getIndex() * 3 * compareAtomsSize;
        for (int i = 0; i < compareAtomsSize; i++) {
          int atomIndex = i * 3;
          baseMolCoords[atomIndex] = baseXYZ[baseCenterMol + atomIndex];
          baseMolCoords[atomIndex + 1] = baseXYZ[baseCenterMol + atomIndex + 1];
          baseMolCoords[atomIndex + 2] = baseXYZ[baseCenterMol + atomIndex + 2];
          targetMolCoords[atomIndex] = targetXYZ[targetCenterMol + atomIndex];
          targetMolCoords[atomIndex + 1] = targetXYZ[targetCenterMol + atomIndex + 1];
          targetMolCoords[atomIndex + 2] = targetXYZ[targetCenterMol + atomIndex + 2];
        }
        translate(baseMolCoords, mass, targetMolCoords, mass);
        rotate(baseMolCoords, targetMolCoords, mass);
        double checkRMSD0 = rmsd(baseMolCoords, targetMolCoords, mass);
        logger.finer(format(" Center Molecule RMSD: %16.8f", checkRMSD0));
        // If center most molecule's RMSD is large any crystal matches may not be accurate.

        if (checkRMSD0 < 5.0 || force) {
          logger.finer(" Rotation 1:");
          //Note "mass" uses masses "massCheck" is not weighted.
          firstRotation(baseXYZ, targetXYZ, mass, molDist12_2[0].getIndex(),
              molDist22_2[0].getIndex());

          //Check center most molecules RMSD between both systems after trans/rot1. Should be same as prior.
//                    if (logger.isLoggable(Level.FINE)) {
          baseMolCoords = new double[compareAtomsSize * 3];
          targetMolCoords = new double[compareAtomsSize * 3];
          baseCenterMol = molDist12_2[0].getIndex() * 3 * compareAtomsSize;
          targetCenterMol = molDist22_2[0].getIndex() * 3 * compareAtomsSize;
          for (int i = 0; i < compareAtomsSize; i++) {
            int atomIndex = i * 3;
            baseMolCoords[atomIndex] = baseXYZ[baseCenterMol + atomIndex];
            baseMolCoords[atomIndex + 1] = baseXYZ[baseCenterMol + atomIndex + 1];
            baseMolCoords[atomIndex + 2] = baseXYZ[baseCenterMol + atomIndex + 2];
            targetMolCoords[atomIndex] = targetXYZ[targetCenterMol + atomIndex];
            targetMolCoords[atomIndex + 1] = targetXYZ[targetCenterMol + atomIndex + 1];
            targetMolCoords[atomIndex + 2] = targetXYZ[targetCenterMol + atomIndex + 2];
          }
          double checkRMSD1 = rmsd(baseMolCoords, targetMolCoords, mass);
          logger.finer(format(" Center Molecule RMSD after rot 1: %16.8f", checkRMSD1));
//                    }

          // Check center of mass for center numCheck most entities and distance to center-most.
          if (logger.isLoggable(Level.FINER)) {
            int numCheck = 7;
            double[][] baseCenters = new double[numCheck][3];
            double[][] targetCenters = new double[numCheck][3];
            for (int j = 0; j < numCheck; j++) {
              int baseCenterMols = molDist12_2[j].getIndex() * 3 * compareAtomsSize;
              int targetCenterMols = molDist22_2[j].getIndex() * 3 * compareAtomsSize;
              for (int i = 0; i < compareAtomsSize; i++) {
                int atomIndex = i * 3;
                baseCenters[j][0] += baseXYZ[baseCenterMols + atomIndex] * mass[i];
                baseCenters[j][1] += baseXYZ[baseCenterMols + atomIndex + 1] * mass[i];
                baseCenters[j][2] += baseXYZ[baseCenterMols + atomIndex + 2] * mass[i];
                targetCenters[j][0] += targetXYZ[targetCenterMols + atomIndex] * mass[i];
                targetCenters[j][1] += targetXYZ[targetCenterMols + atomIndex + 1] * mass[i];
                targetCenters[j][2] += targetXYZ[targetCenterMols + atomIndex + 2] * mass[i];
              }
              for (int i = 0; i < 3; i++) {
                baseCenters[j][i] /= massSum;
                targetCenters[j][i] /= massSum;
              }
              logger.finer(format(
                  " Position: %d Base(%4.4f): %8.4f %8.4f %8.4f Target(%4.4f): %8.4f %8.4f %8.4f", j,
                  molDist12_2[j].getDoubleValue(), baseCenters[j][0], baseCenters[j][1],
                  baseCenters[j][2],
                  molDist22_2[j].getDoubleValue(), targetCenters[j][0], targetCenters[j][1],
                  targetCenters[j][2]));
            }
          }

          // At this point both systems have completed first rotation/translation
          //  Therefore both center-most molecules should be overlapped.
          // TODO could have prioritization favor closest distance from molDist12 or molDist22.
          logger.finer(" Match molecules between systems");

          int nMols = Math.min(nBaseMols, nTargetMols);

          DoubleIndexPair[] matchMols = new DoubleIndexPair[nMols];

          boolean reverse = false;
          if (symmetric) {
            double dist1 = molDist12_2[1].getDoubleValue();
            double dist2 = molDist22_2[1].getDoubleValue();
            if (dist1 < dist2) {
              // Use dist2
              reverse = true;
              matchMolecules(matchMols, targetXYZ, baseXYZ, mass, molDist22_2, molDist12_2);
            } else {
              matchMolecules(matchMols, baseXYZ, targetXYZ, mass, molDist12_2, molDist22_2);
            }
          } else {
            matchMolecules(matchMols, baseXYZ, targetXYZ, mass, molDist12_2, molDist22_2);
          }

          int numMolMatch = matchMols.length;

          if (printSize > numMolMatch) {
            printSize = numMolMatch;
          }

          if (logger.isLoggable(Level.FINER)) {
            logger.finer("  Distance between pairs after rot 1:");
            for (int i = 0; i < printSize; i++) {
              logger.finer(
                  format(" %2d %16.8f", matchMols[i].getIndex(), matchMols[i].getDoubleValue()));
            }
            logger.finer(" Index  MolDist1    MolDist2    MatchMols");
            for (int i = 0; i < printSize; i++) {
              logger.finer(format(" %2d base: %2d target: %2d Match: %2d", i,
                  molDist12_2[i].getIndex(), molDist22_2[i].getIndex(), matchMols[i].getIndex()));
            }
          }

          logger.finer(" Rotation 2:");
          //Note "mass" uses masses "massCheck" is not weighted.
          if (reverse) {
            secondRotation(targetXYZ, baseXYZ, mass, molDist22_2, matchMols);
          } else {
            secondRotation(baseXYZ, targetXYZ, mass, molDist12_2, matchMols);
          }

          baseMolCoords = new double[compareAtomsSize * 3 * 3];
          targetMolCoords = new double[compareAtomsSize * 3 * 3];
          for (int i = 0; i < 3; i++) {
            int molIndex = i * compareAtomsSize;
            baseCenterMol = molDist12_2[i].getIndex() * 3 * compareAtomsSize;
            targetCenterMol = matchMols[i].getIndex() * 3 * compareAtomsSize;
            for (int j = 0; j < compareAtomsSize; j++) {
              int atomIndex = j * 3;
              baseMolCoords[molIndex + atomIndex] = baseXYZ[baseCenterMol + atomIndex];
              baseMolCoords[molIndex + atomIndex + 1] = baseXYZ[baseCenterMol + atomIndex + 1];
              baseMolCoords[molIndex + atomIndex + 2] = baseXYZ[baseCenterMol + atomIndex + 2];
              targetMolCoords[molIndex + atomIndex] = targetXYZ[targetCenterMol + atomIndex];
              targetMolCoords[molIndex + atomIndex + 1] = targetXYZ[targetCenterMol + atomIndex + 1];
              targetMolCoords[molIndex + atomIndex + 2] = targetXYZ[targetCenterMol + atomIndex + 2];
            }
          }

          double[] mass3 = new double[compareAtomsSize * 3];
          for (int i = 0; i < 3; i++) {
            arraycopy(mass, 0, mass3, i * compareAtomsSize, compareAtomsSize);
          }
          double checkRMSD2 = rmsd(baseMolCoords, targetMolCoords, mass3);

          // Rotations 1 and 2 have been completed and both systems should be overlapped
          //  Isolate center most nAU from System 1 and matching molecules from System 2
          logger.finer(" Final rotation:");
          double[] baseNMols = new double[nAU * compareAtomsSize * 3];
          if (reverse) {
            for (int i = 0; i < nAU; i++) {
              int offset = i * compareAtomsSize * 3;
              int molIndex = molDist22_2[i].getIndex() * compareAtomsSize * 3;
              for (int j = 0; j < compareAtomsSize; j++) {
                int atomIndex = j * 3;
                baseNMols[offset + atomIndex] = targetXYZ[molIndex + atomIndex];
                baseNMols[offset + atomIndex + 1] = targetXYZ[molIndex + atomIndex + 1];
                baseNMols[offset + atomIndex + 2] = targetXYZ[molIndex + atomIndex + 2];
              }
            }
          } else {
            for (int i = 0; i < nAU; i++) {
              int offset = i * compareAtomsSize * 3;
              int molIndex = molDist12_2[i].getIndex() * compareAtomsSize * 3;
              for (int j = 0; j < compareAtomsSize; j++) {
                int atomIndex = j * 3;
                baseNMols[offset + atomIndex] = baseXYZ[molIndex + atomIndex];
                baseNMols[offset + atomIndex + 1] = baseXYZ[molIndex + atomIndex + 1];
                baseNMols[offset + atomIndex + 2] = baseXYZ[molIndex + atomIndex + 2];
              }
            }
          }

          logger.finer(" Match Molecules:");
          matchMols = new DoubleIndexPair[nAU];
          if (reverse) {
            matchMolecules(matchMols, targetXYZ, baseXYZ, mass, molDist22_2, molDist12_2);
          } else {
            matchMolecules(matchMols, baseXYZ, targetXYZ, mass, molDist12_2, molDist22_2);
          }

          if (printSize > nAU) {
            printSize = nAU;
          }

          if (logger.isLoggable(Level.FINER)) {
            logger.finer("  Distance between pairs after rot 2:");
            for (int i = 0; i < printSize; i++) {
              logger.finer(
                  format(" %2d %16.8f", matchMols[i].getIndex(), matchMols[i].getDoubleValue()));
            }
          }

          double[] targetNMols = new double[nAU * compareAtomsSize * 3];
          if (reverse) {
            for (int i = 0; i < nAU; i++) {
              int offset = i * compareAtomsSize * 3;
              int molIndex = matchMols[i].getIndex() * compareAtomsSize * 3;
              for (int j = 0; j < compareAtomsSize; j++) {
                int atomIndex = j * 3;
                targetNMols[offset + atomIndex] = baseXYZ[molIndex + atomIndex];
                targetNMols[offset + atomIndex + 1] = baseXYZ[molIndex + atomIndex + 1];
                targetNMols[offset + atomIndex + 2] = baseXYZ[molIndex + atomIndex + 2];
              }
            }
          } else {
            for (int i = 0; i < nAU; i++) {
              int offset = i * compareAtomsSize * 3;
              int molIndex = matchMols[i].getIndex() * compareAtomsSize * 3;
              for (int j = 0; j < compareAtomsSize; j++) {
                int atomIndex = j * 3;
                targetNMols[offset + atomIndex] = targetXYZ[molIndex + atomIndex];
                targetNMols[offset + atomIndex + 1] = targetXYZ[molIndex + atomIndex + 1];
                targetNMols[offset + atomIndex + 2] = targetXYZ[molIndex + atomIndex + 2];
              }
            }
          }

          double[] massN = new double[compareAtomsSize * nAU];
          for (int i = 0; i < nAU; i++) {
            arraycopy(mass, 0, massN, i * compareAtomsSize, compareAtomsSize);
          }

          translate(baseNMols, massN, targetNMols, massN);
          rotate(baseNMols, targetNMols, massN);

          if (save) {
            saveAssemblyPDB(baseAssembly, baseNMols, "_" + l + "_" + m + "_c1");
            saveAssemblyPDB(molecularAssembly, targetNMols, "_" + l + "_" + m + "_c2");
          }

          //Check center most molecules RMSD between both systems after trans/rot2.
          //   Debug check can be removed later.
          baseMolCoords = new double[compareAtomsSize * 3];
          targetMolCoords = new double[compareAtomsSize * 3];
          for (int i = 0; i < compareAtomsSize; i++) {
            int atomIndex = i * 3;
            baseMolCoords[atomIndex] = baseNMols[atomIndex];
            baseMolCoords[atomIndex + 1] = baseNMols[atomIndex + 1];
            baseMolCoords[atomIndex + 2] = baseNMols[atomIndex + 2];
            targetMolCoords[atomIndex] = targetNMols[atomIndex];
            targetMolCoords[atomIndex + 1] = targetNMols[atomIndex + 1];
            targetMolCoords[atomIndex + 2] = targetNMols[atomIndex + 2];
          }
          double checkRMSD3 = rmsd(baseMolCoords, targetMolCoords, mass);
          logger.finer(format(" Center Molecule RMSD final: %16.8f", checkRMSD3));

          double rmsdSymOp = rmsd(baseNMols, targetNMols, massN);
          logger.fine(format(" %2d of %2d: %7.4f %7.4f %7.4f", l * numSearch2 + m + 1,
              numSearch * numSearch2, checkRMSD1, checkRMSD2, rmsdSymOp));
          if (rmsdSymOp < bestRMSD) {
            bestRMSD = rmsdSymOp;
          }
          addLooseUnequal(crystalCheckRMSDs, round(rmsdSymOp, 6));
        } else {
          logger.finer(
              format(" Center molecule %2d has an RMSD of %7.4f and was skipped.", m, checkRMSD0));
        }
      }
    }

    double finalRMSD = Double.NaN;
    if (bestRMSD < Double.MAX_VALUE) {
      finalRMSD = bestRMSD;
    } else {
      logger.warning(" RMSD was filtered out! " +
          "Check if molecules are the same and have the same atom ordering." +
          "\n If you really need a value try the -f flag.");
    }

    // Logging to check number of RMSD values determined.
    StringBuilder dblOut = new StringBuilder();
    for (Double dbl : crystalCheckRMSDs) {
      dblOut.append(" ").append(format("%4.4f", dbl));
    }
    String message = format(" Unique %s Values: %s", rmsdLabel, dblOut);

    int numUnique = crystalCheckRMSDs.size();
    if (numUnique > 4) {
      logger.finer(format(
          " PAC determined %2d unique values. Consider increasing the number of inflated molecules.\n %s",
          numUnique, message));
    } else {
      logger.finer(message);
    }

    return finalRMSD;
  }

  /**
   * Read in the distance matrix.
   *
   * @param filename The PAC RMSD matrix file to read from.
   */
  private void readMatrix(String filename) {
    restartRow = 0;
    restartColumn = 0;

    DistanceMatrixFilter distanceMatrixFilter = new DistanceMatrixFilter();
    if (distanceMatrixFilter.readDistanceMatrix(filename, distMatrix)) {
      restartRow = distanceMatrixFilter.getRestartRow();
      restartColumn = distanceMatrixFilter.getRestartColumn();
    }

    int nRow = distMatrix.length;
    int nColumn = distMatrix[0].length;
    if (restartRow == nRow && restartColumn == nColumn) {
      logger.info(format(" Complete distance matrix found (%d x %d).", restartRow, restartColumn));
    } else {
      restartColumn = 0;
      logger.info(format(" Incomplete distance matrix found.\n Restarting at row %d, column %d.",
          restartRow + 1, restartColumn + 1));
    }
  }

  /**
   * This method calls <code>world.AllGather</code> to collect numProc PAC RMSD values.
   *
   * @param row Current row of the PAC RMSD matrix.
   * @param column Current column of the PAC RMSD matrix.
   */
  private void gatherRMSDs(int row, int column) {
    try {
      logger.finer(" Receiving results.");
      world.allGather(myBuffer, buffers);
      for (int i = 0; i < numProc; i++) {
        int c = (column + 1) - numProc + i;
        if (c < targetSize) {
          distMatrix[row][c] = distances[i][0];
          logger.finer(format(" %d %d %16.8f", row, c, distances[i][0]));
        }
      }
    } catch (Exception e) {
      logger.severe(" Exception collecting distance values." + e);
    }
  }

  /**
   * Perform first rotation to match center most molecules.
   *
   * @param baseXYZ Coordinates for system 1
   * @param targetXYZ Coordinates for system 2
   * @param mass Masses for atoms in a molecule
   * @param index1 Index of center most molecule in coords
   * @param index2 Index of center most molecule in coords
   */
  private static void firstRotation(double[] baseXYZ, double[] targetXYZ, double[] mass,
      int index1, int index2) {
    int nAtoms = mass.length;
    double[] baseMol = new double[nAtoms * 3];
    double[] targetMol = new double[nAtoms * 3];
    int baseIndex = index1 * nAtoms * 3;
    int targetIndex = index2 * nAtoms * 3;
    for (int i = 0; i < nAtoms; i++) {
      int atomIndex = i * 3;
      baseMol[atomIndex] = baseXYZ[baseIndex + atomIndex];
      baseMol[atomIndex + 1] = baseXYZ[baseIndex + atomIndex + 1];
      baseMol[atomIndex + 2] = baseXYZ[baseIndex + atomIndex + 2];
      targetMol[atomIndex] = targetXYZ[targetIndex + atomIndex];
      targetMol[atomIndex + 1] = targetXYZ[targetIndex + atomIndex + 1];
      targetMol[atomIndex + 2] = targetXYZ[targetIndex + atomIndex + 2];
    }
    double[] translation = calculateTranslation(baseMol, mass);
    applyTranslation(baseMol, translation);
    applyTranslation(baseXYZ, translation);
    translation = calculateTranslation(targetMol, mass);
    applyTranslation(targetMol, translation);
    applyTranslation(targetXYZ, translation);
    applyRotation(targetXYZ, calculateRotation(baseMol, targetMol, mass));
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
    double[] baseMol = new double[nAtoms * 3];
    int baseIndex = index * nAtoms * 3;
    for (int i = 0; i < nAtoms; i++) {
      int atomIndex = i * 3;
      baseMol[atomIndex] = baseXYZ[baseIndex + atomIndex];
      baseMol[atomIndex + 1] = baseXYZ[baseIndex + atomIndex + 1];
      baseMol[atomIndex + 2] = baseXYZ[baseIndex + atomIndex + 2];
    }
    double[] translateBase = calculateTranslation(baseMol, mass);
    applyTranslation(baseMol, translateBase);
    applyTranslation(baseXYZ, translateBase);
  }

  /**
   * Perform second rotation to better match crystal systems.
   *
   * @param baseXYZ Coordinates for system 1
   * @param targetXYZ Coordinates for system 2
   * @param mass Masses for atoms in one molecule
   * @param molDist12 Indices of molecules closest to center
   * @param matchMols Indices for system 2 that match molecules in system 1
   */
  private static void secondRotation(double[] baseXYZ, double[] targetXYZ, double[] mass,
      DoubleIndexPair[] molDist12, DoubleIndexPair[] matchMols) {
    int nAtoms = mass.length;
    double[] base3Mols = new double[nAtoms * 3 * 3];
    double[] target3Mols = new double[nAtoms * 3 * 3];
    for (int i = 0; i < 3; i++) {
      int index = i * nAtoms * 3;
      int baseIndex = molDist12[i].getIndex() * nAtoms * 3;
      int targetIndex = matchMols[i].getIndex() * nAtoms * 3;
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = j * 3;
        base3Mols[index + atomIndex] = baseXYZ[baseIndex + atomIndex];
        base3Mols[index + atomIndex + 1] = baseXYZ[baseIndex + atomIndex + 1];
        base3Mols[index + atomIndex + 2] = baseXYZ[baseIndex + atomIndex + 2];
        target3Mols[index + atomIndex] = targetXYZ[targetIndex + atomIndex];
        target3Mols[index + atomIndex + 1] = targetXYZ[targetIndex + atomIndex + 1];
        target3Mols[index + atomIndex + 2] = targetXYZ[targetIndex + atomIndex + 2];
      }
    }

    double[] mass3 = new double[nAtoms * 3];
    for (int i = 0; i < nAtoms; i++) {
      mass3[i] = mass[i];
      mass3[nAtoms + i] = mass[i];
      mass3[2 * nAtoms + i] = mass[i];
    }

    applyRotation(targetXYZ, calculateRotation(base3Mols, target3Mols, mass3));
  }

  /**
   * Match molecules between two crystal systems.
   *
   * @param matchMols Matched indices to be returned
   * @param baseXYZ Coordinates for system 1
   * @param targetXYZ Coordinates for system 2
   * @param mass Array of masses for atoms in molecule.
   * @param molDist12 Indices of system 1 closest to center.
   * @param molDist22 Indices of system 2 closest to center.
   */
  private static void matchMolecules(DoubleIndexPair[] matchMols, double[] baseXYZ,
      double[] targetXYZ,
      double[] mass, DoubleIndexPair[] molDist12, DoubleIndexPair[] molDist22) {
    int nAtoms = mass.length;
    double massSum = 0.0;
    for (double value : mass) {
      massSum += value;
    }
    // Determine the center of mass for each entity in system 1.
    int desiredMols = matchMols.length;
    double[][] baseMols = new double[desiredMols][3];
    centerOfMass(baseMols, baseXYZ, mass, molDist12);

    int nTargetMols = targetXYZ.length / (nAtoms * 3);
    // List of indexes for second system.
    List<Integer> targetIndex = new ArrayList<>(nTargetMols);
    for (DoubleIndexPair doubleIndexPair : molDist22) {
      targetIndex.add(doubleIndexPair.getIndex());
    }

    // Compare distances between center of masses from system 1 and 2.
    for (int i = 0; i < desiredMols; i++) {
      double minDist = Double.MAX_VALUE;
      Integer minIndex = -1;
      for (Integer target : targetIndex) {
        double[] targetMol = new double[3];
        int molIndex = target * nAtoms * 3;
        for (int k = 0; k < nAtoms; k++) {
          int atomIndex = k * 3;
          targetMol[0] += targetXYZ[molIndex + atomIndex] * mass[k];
          targetMol[1] += targetXYZ[molIndex + atomIndex + 1] * mass[k];
          targetMol[2] += targetXYZ[molIndex + atomIndex + 2] * mass[k];
        }
        for (int k = 0; k < 3; k++) {
          targetMol[k] /= massSum;
        }
        double dist = dist(baseMols[i], targetMol);
        if (dist < minDist) {
          minDist = dist;
          minIndex = target;
        }
        if (abs(minDist) < MATCH_TOLERANCE) {
          // If distance between center of masses us 0 that's as good as it gets assuming no internal overlaps.
          break;
        }
      }
      matchMols[i] = new DoubleIndexPair(minIndex, minDist);
      if (!targetIndex.remove(minIndex)) {
        logger.warning(format(" Index value of %d was not found (%4.4f).", minIndex, minDist));
      }
    }

    if (logger.isLoggable(Level.FINER)) {
      double[][] targetCenterMols = new double[nTargetMols][3];
      centerOfMass(targetCenterMols, targetXYZ, mass, matchMols);
      for (int i = 0; i < desiredMols; i++) {
        logger.finer(
            format(" Base position:   %d: %8.4f %8.4f %8.4f", i, baseMols[i][0], baseMols[i][1],
                baseMols[i][2]));
        logger.finer(format(" Match Distance:  %d: %8.4f", i, matchMols[i].getDoubleValue()));
        logger.finer(format(" Target position: %d: %8.4f %8.4f %8.4f", i, targetCenterMols[i][0],
            targetCenterMols[i][1], targetCenterMols[i][2]));
      }
    }
  }

  /**
   * Calculate the center of mass for a given set of masses for the asymmetric unit and coordinates
   * (xyz)
   *
   * @param centersOfMass Returned center of mass for each asymmetric unit
   * @param coords Coordinates of every atom in system.
   * @param mass Masses of each atom in asymetric unit.
   * @param order Order to calculate the center of masses.
   */
  private static void centerOfMass(double[][] centersOfMass, double[] coords, double[] mass,
      DoubleIndexPair[] order) {
    int size = centersOfMass.length;
    if (order.length < size) {
      size = order.length;
    }
    int nAtoms = mass.length;
    double massSum = 0.0;
    for (double value : mass) {
      massSum += value;
    }
    for (int i = 0; i < size; i++) {
      int index = order[i].getIndex();
      int molIndex = index * nAtoms * 3;
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = j * 3;
        centersOfMass[i][0] += coords[molIndex + atomIndex] * mass[j];
        centersOfMass[i][1] += coords[molIndex + atomIndex + 1] * mass[j];
        centersOfMass[i][2] += coords[molIndex + atomIndex + 2] * mass[j];
      }
      for (int j = 0; j < 3; j++) {
        centersOfMass[i][j] /= massSum;
      }
    }
  }

  /**
   * Generate and expanded sphere of asymetric unit with the intention of observing a crystals
   * distribution of replicates rather to facilitate comparisons that go beyond lattice parameters.
   *
   * @param molecularAssembly Asymmetric unit to be expanded.
   * @param inflatedAU Minimum number of replicates desired in resultant sphere.
   * @return MolecularAssembly containing the expanded sphere.
   */
  private static MolecularAssembly generateInflatedSphere(MolecularAssembly molecularAssembly,
      int inflatedAU) {
    // Move molecules into the unit cell.
    molecularAssembly.moveAllIntoUnitCell();

    Atom[] atoms = molecularAssembly.getAtomArray();
    int nAtoms = atoms.length;
    // Collect asymmetric unit atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];
    double[] xf = new double[nAtoms];
    double[] yf = new double[nAtoms];
    double[] zf = new double[nAtoms];
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      x[i] = atom.getX();
      y[i] = atom.getY();
      z[i] = atom.getZ();
    }

    // When the system was read in, a replicates crystal may have been created to satisfy the cutoff.
    // Retrieve a reference to the unit cell (not the replicates crystal).
    // Here we will use the unit cell, to create a new replicates crystal that may be
    // a different size (i.e. larger).
    Crystal unitCell = molecularAssembly.getCrystal().getUnitCell();

    unitCell.toFractionalCoordinates(nAtoms, x, y, z, xf, yf, zf);

    double asymmetricUnitVolume = unitCell.volume / unitCell.getNumSymOps();

    // Estimate a radius that will include desired number of asymmetric units (inflatedAU).
    double radius = cbrt(inflatedAU * asymmetricUnitVolume) / 2;
    // double radius = cbrt(((3.0 / (4.0 * PI)) * inflatedAU * asymmetricUnitVolume));

    logger.finer(format(" Desired copies in target sphere:     %3d", inflatedAU));
    logger.finer(format(" Asymmetric Unit Volume:  %4.2f", asymmetricUnitVolume));
    logger.finer(format(" Estimated spherical radius:  %4.2f", radius));

    Crystal replicatesCrystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, radius * 2.0);

    logger.finer(" Replicates crystal " + replicatesCrystal.toString());

    // Update the Molecular Assembly with the ReplicatesCrystal (used outside this method).
    molecularAssembly.setCrystal(replicatesCrystal);

    // Symmetry coordinates for each molecule in replicates crystal
    int nSymm = replicatesCrystal.getNumSymOps();
    logger.finer(format(" Number of replicates: %3d", nSymm));
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

      double[] xSf = new double[nAtoms];
      double[] ySf = new double[nAtoms];
      double[] zSf = new double[nAtoms];

      replicatesCrystal.toFractionalCoordinates(nAtoms, xS[iSym], yS[iSym], zS[iSym], xSf, ySf, zSf);

      // Compute center-of-mass (CoM) for Cartesian coordinates
      double[] centerOfMass = new double[3];
      int index = 0;
      double totalMass = 0.0;
      for (Atom atom : atoms) {
        double m = atom.getMass();
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
    logger.finer(format(" Expanded Crystal Center: %16.8f %16.8f %16.8f",
        cartCenter[0], cartCenter[1], cartCenter[2]));

    for (int i = 0; i < nSymm; i++) {
      // Then compute Euclidean distance from Cartesian center of the replicates cell
      molsDists[i] = new DoubleIndexPair(i, dist(cartCenter, centerMolsCart[i]));
    }

    // Sort the molecules by their distance from the center.
    // Note that the smallest distances are first in the array after the sort.
    Arrays.parallelSort(molsDists);

    // Create expanded assembly
    MolecularAssembly expandedAssembly = new MolecularAssembly(molecularAssembly.getName());
    List<Bond> bondList = molecularAssembly.getBondList();
    ArrayList<Atom> newAtomList = new ArrayList<>();
    int atomIndex = 0;
    logger.finest("\n Copy  SymOp        Distance");
    for (int n = 0; n < nSymm; n++) {
      // Current molecule
      int iSym = molsDists[n].getIndex();
      double distance = molsDists[n].getDoubleValue();
      logger.finest(format(" %4d  %5d  %16.8f", n, iSym, distance));

      ArrayList<Atom> atomList = new ArrayList<>();
      // Create a new set of Atoms for each SymOp molecule
      for (int i = 0; i < nAtoms; i++) {
        Atom a = atoms[i];
        double[] xyz = new double[3];
        xyz[0] = xS[iSym][i];
        xyz[1] = yS[iSym][i];
        xyz[2] = zS[iSym][i];
        Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
        atomList.add(atom);
      }
      // Create a new set of Bonds for each SymOp molecule
      for (Bond bond : bondList) {
        Atom a1 = bond.getAtom(0);
        Atom a2 = bond.getAtom(1);
        Atom newA1 = atomList.get(a1.getIndex() - 1);
        Atom newA2 = atomList.get(a2.getIndex() - 1);
        Bond b = new Bond(newA1, newA2);
        b.setBondType(bond.getBondType());
      }
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

    expandedAssembly.setForceField(forceField);

    // The biochemistry method is designed to load chemical entities into the
    // Polymer, Molecule, Water and Ion data structure.
    Utilities.biochemistry(expandedAssembly, newAtomList);

    expandedAssembly.setFile(molecularAssembly.getFile());
    expandedAssembly.setCenter(cartCenter);

    return expandedAssembly;
  }

  /**
   * Produce a translation vector necessary to move an object with the current center of mass (com)
   * into the provided crystal.
   *
   * @param cryst Replicates crystal within whom coordinates should be moved.
   * @param com Center of mass (x, y, z) for the object of concern
   * @return double[] translation vector to move the object within the provided crystal.
   */
  private static double[] moveIntoCrystal(Crystal cryst, double[] com) {

    double[] translate = new double[3];
    double[] currentCoM = new double[3];
    currentCoM[0] = com[0];
    currentCoM[1] = com[1];
    currentCoM[2] = com[2];

    // Move the COM to the Replicates Crystal.
    cryst.toFractionalCoordinates(com, translate);
    translate[0] = mod(translate[0], 1.0);
    translate[1] = mod(translate[1], 1.0);
    translate[2] = mod(translate[2], 1.0);
    cryst.toCartesianCoordinates(translate, translate);

    // Correct center of mass.
    com[0] = translate[0];
    com[1] = translate[1];
    com[2] = translate[2];

    // The translation vector is difference between the new location and the current COM.
    translate[0] -= currentCoM[0];
    translate[1] -= currentCoM[1];
    translate[2] -= currentCoM[2];

    logger.finest(format(" Center of Mass Prior: %16.8f %16.8f %16.8f", currentCoM[0], currentCoM[1],
        currentCoM[2]));
    logger.finest(format(" Center of Mass Post: %16.8f %16.8f %16.8f", com[0], com[1], com[2]));

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
      DoubleIndexPair[] molDists, DoubleIndexPair[] molDists2) {
    prioritizeReplicates(nAtoms, coordsXYZ, mass, massSum, nMols, molDists, 0, molDists2);
  }

  /**
   * Prioritize asymmetric units within the system based on distance to specified index.
   *
   * @param nAtoms Number of atoms being compared in asymmetric unit.
   * @param coordsXYZ Coordinates for expanded crystal.
   * @param mass Mass of atoms within asymmetric unit.
   * @param massSum Sum of mass values.
   * @param nMols Number of molecules in expanded crystal.
   * @param molDists Prioritization of molecules in expanded system on distance to center-most
   *     molecule.
   * @param index Index of molecules to be center.
   * @param molDists2 Prioritization of molecules based on two center-most atoms.
   */
  private static void prioritizeReplicates(int nAtoms, double[] coordsXYZ, double[] mass,
      double massSum, int nMols,
      DoubleIndexPair[] molDists, int index, DoubleIndexPair[] molDists2) {
    // Molecules added to system based on distance to center. 0th molecule should be closest.
    int centerMolIndex = index * nAtoms * 3;
    // Center most entitiy's center of mass
    double[] coordCenter = new double[3];
    for (int j = 0; j < nAtoms; j++) {
      int atomIndex = j * 3;
      coordCenter[0] += coordsXYZ[centerMolIndex + atomIndex] * mass[j];
      coordCenter[1] += coordsXYZ[centerMolIndex + atomIndex + 1] * mass[j];
      coordCenter[2] += coordsXYZ[centerMolIndex + atomIndex + 2] * mass[j];
    }
    for (int j = 0; j < 3; j++) {
      coordCenter[j] /= massSum;
    }

    for (int i = 0; i < nMols; i++) {
      double[] moleculeCenter = new double[3];
      int molIndex = i * nAtoms * 3;
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = j * 3;
        moleculeCenter[0] += coordsXYZ[molIndex + atomIndex] * mass[j];
        moleculeCenter[1] += coordsXYZ[molIndex + atomIndex + 1] * mass[j];
        moleculeCenter[2] += coordsXYZ[molIndex + atomIndex + 2] * mass[j];
      }
      for (int j = 0; j < 3; j++) {
        moleculeCenter[j] /= massSum;
      }
      molDists[i] = new DoubleIndexPair(i, dist(coordCenter, moleculeCenter));
    }
    Arrays.sort(molDists);

    if (logger.isLoggable(Level.FINER)) {
      int numCheck = 7;
      double[][] targetMol = new double[coordsXYZ.length / (nAtoms * 3)][3];
      centerOfMass(targetMol, coordsXYZ, mass, molDists);
      for (int i = 0; i < numCheck; i++) {
        logger.finer(format("Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
            molDists[i].getIndex(),
            molDists[i].getDoubleValue(), targetMol[i][0], targetMol[i][1], targetMol[i][2]));
      }
    }

    //Molecules in crystal sorted based on distance to center of mass of center most molecule
    //  Want the first two molecules chosen in this manner, but third molecule to be closest to both

    // Assigning distance from molDists ensures correct ordering of center most molecule.
    molDists2[0] = new DoubleIndexPair(molDists[0].getIndex(), molDists[0].getDoubleValue());
    double[] avgCenter = new double[3];
    for (int i = 0; i < 2; i++) {
      int molIndex = nAtoms * molDists[i].getIndex() * 3;
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = j * 3;
        avgCenter[0] += coordsXYZ[molIndex + atomIndex] * mass[j];
        avgCenter[1] += coordsXYZ[molIndex + atomIndex + 1] * mass[j];
        avgCenter[2] += coordsXYZ[molIndex + atomIndex + 2] * mass[j];
      }
    }
    for (int j = 0; j < 3; j++) {
      avgCenter[j] /= (2 * massSum);
    }

    for (int i = 1; i < nMols; i++) {
      double[] moleculeCenter = new double[3];
      int molIndex = nAtoms * molDists[i].getIndex() * 3;
      for (int j = 0; j < nAtoms; j++) {
        int atomIndex = j * 3;
        moleculeCenter[0] += coordsXYZ[molIndex + atomIndex] * mass[j];
        moleculeCenter[1] += coordsXYZ[molIndex + atomIndex + 1] * mass[j];
        moleculeCenter[2] += coordsXYZ[molIndex + atomIndex + 2] * mass[j];
      }
      for (int j = 0; j < 3; j++) {
        moleculeCenter[j] /= massSum;
      }
      molDists2[i] = new DoubleIndexPair(molDists[i].getIndex(), dist(avgCenter, moleculeCenter));
    }
    Arrays.sort(molDists2);
  }

  /**
   * Round a value to a given number of decimal places.
   *
   * @param value Value to be rounded.
   * @param places Desired number of decimal places.
   * @return Rounded value.
   */
  private static double round(double value, int places) {
    BigDecimal bd = new BigDecimal(Double.toString(value));
    bd = bd.setScale(places, RoundingMode.HALF_UP);
    return bd.doubleValue();
  }

  /**
   * Add a value to a list of long if it's difference to all listed values is greater than the
   * tolerance.
   *
   * @param values List of values already found.
   * @param value Potential new value if it is not already in list.
   */
  private static void addLooseUnequal(List<Double> values, double value) {
    boolean found = false;
    for (Double dbl : values) {
      if (abs(dbl - value) < 1.0E-6) {
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
      String description) {
    String fileName = FilenameUtils.removeExtension(molecularAssembly.getName());
    File saveLocationPDB = new File(fileName + description + ".pdb");
    int fileNumber = 0;
    while (saveLocationPDB.exists()) {
      saveLocationPDB = new File(fileName + description + "_" + fileNumber++ + ".pdb");
    }
    // Save aperiodic system of n_mol closest atoms for visualization
    MolecularAssembly currentAssembly = new MolecularAssembly(molecularAssembly.getName());
    List<Bond> bondList = molecularAssembly.getBondList();
    ArrayList<Atom> newAtomList = new ArrayList<>();
    Atom[] atoms = molecularAssembly.getAtomArray();
    int atomIndex = 0;
    int nAtoms = atoms.length;
    int numMols = coords.length / (3 * nAtoms);
    for (int n = 0; n < numMols; n++) {
      ArrayList<Atom> atomList = new ArrayList<>();
      // Create a new set of Atoms for each SymOp molecule
      for (int i = 0; i < nAtoms; i++) {
        Atom a = atoms[i];
        int atomValue = i * 3;
        double[] xyz = new double[3];
        xyz[0] = coords[n * nAtoms * 3 + atomValue];
        xyz[1] = coords[n * nAtoms * 3 + atomValue + 1];
        xyz[2] = coords[n * nAtoms * 3 + atomValue + 2];
        Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
        atomList.add(atom);
      }
      // Create a new set of Bonds for each SymOp molecule
      for (Bond bond : bondList) {
        Atom a1 = bond.getAtom(0);
        Atom a2 = bond.getAtom(1);
        Atom newA1 = atomList.get(a1.getIndex() - 1);
        Atom newA2 = atomList.get(a2.getIndex() - 1);
        Bond b = new Bond(newA1, newA2);
        b.setBondType(bond.getBondType());
      }
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
  }

  /**
   * Orient coordinates so that the second index is on the x axis, and the third index is on the X-Y
   * plane. First index should be at the origin (0, 0, 0).
   *
   * @param coordsXYZ An array of XYZ positions (e.g. [x0, y0, z0, x1, y1, z1, x2, y2, z2]
   * @param atomIndices Indices for three desired sets from the XYZ list (e.g. [0, 1, 2]). Index
   *     0 should be at origin.
   */
  public static void standardOrientation(double[] coordsXYZ, int[] atomIndices) {
    // Used in QEtoXYZ.groovy which is not ready for git which is why this method appears unused.
    int numCoords = coordsXYZ.length / 3;
    double[] atomCoords = new double[3 * 3];
    atomCoords[0] = coordsXYZ[atomIndices[0]];
    atomCoords[1] = coordsXYZ[atomIndices[0] + 1];
    atomCoords[2] = coordsXYZ[atomIndices[0] + 2];
    atomCoords[3] = coordsXYZ[atomIndices[1]];
    atomCoords[4] = coordsXYZ[atomIndices[1] + 1];
    atomCoords[5] = coordsXYZ[atomIndices[1] + 2];
    atomCoords[6] = coordsXYZ[atomIndices[2]];
    atomCoords[7] = coordsXYZ[atomIndices[2] + 1];
    atomCoords[8] = coordsXYZ[atomIndices[2] + 2];

    // TODO: Delete coordsXYZOrig?
    double[] coordsXYZOrig = new double[numCoords * 3];
    for (int i = 0; i < numCoords; i++) {
      int atomIndex = i * 3;
      coordsXYZOrig[atomIndex] = coordsXYZ[atomIndex];
      coordsXYZOrig[atomIndex + 1] = coordsXYZ[atomIndex + 1];
      coordsXYZOrig[atomIndex + 2] = coordsXYZ[atomIndex + 2];
    }

    // TODO: Delete atomsCoordsOrig?
    double[] atomsCoordsOrig = new double[3 * 3];
    arraycopy(atomCoords, 0, atomsCoordsOrig, 0, 9);

    logger.fine(
        format(" START: N1:\t%16.15f %16.15f %16.15f", atomCoords[0], atomCoords[1], atomCoords[2]));
    logger.fine(
        format(" START: N2:\t%16.15f %16.15f %16.15f", atomCoords[3], atomCoords[4], atomCoords[5]));
    logger.fine(
        format(" START: N3:\t%16.15f %16.15f %16.15f", atomCoords[6], atomCoords[7], atomCoords[8]));

    double p1n2 = coordsXYZ[atomIndices[1]];
    double q1n2 = coordsXYZ[atomIndices[1] + 1];
    double r1n2 = coordsXYZ[atomIndices[1] + 2];

    // Calculation of sigma, phai, and cita angles needed to get specified atoms to desired loci
    double cita0 = acos(p1n2 / sqrt(p1n2 * p1n2 + q1n2 * q1n2));
    double phai0 = acos(sqrt(p1n2 * p1n2 + q1n2 * q1n2) /
        sqrt(p1n2 * p1n2 + q1n2 * q1n2 + r1n2 * r1n2));
    if (q1n2 < 0.0) {
      cita0 = -cita0;
    }

    for (int i = 0; i < numCoords; i++) {
      int atomIndex = i * 3;
      double ptmp = coordsXYZ[atomIndex] * cos(cita0) + coordsXYZ[atomIndex + 1] * sin(cita0);
      double qtmp = -coordsXYZ[atomIndex] * sin(cita0) + coordsXYZ[atomIndex + 1] * cos(cita0);
      coordsXYZ[atomIndex] = ptmp;
      coordsXYZ[atomIndex + 1] = qtmp;
    }

    p1n2 = coordsXYZ[atomIndices[1]];
    q1n2 = coordsXYZ[atomIndices[1] + 1];

    if (r1n2 > 0.0) {
      phai0 = -phai0;
    }

    for (int i = 0; i < numCoords; i++) {
      int atomIndex = i * 3;
      double ptmp = coordsXYZ[atomIndex] * cos(phai0) - coordsXYZ[atomIndex + 2] * sin(phai0);
      double rtmp = coordsXYZ[atomIndex] * sin(phai0) + coordsXYZ[atomIndex + 2] * cos(phai0);
      coordsXYZ[atomIndex] = ptmp;
      coordsXYZ[atomIndex + 2] = rtmp;
    }

    p1n2 = coordsXYZ[atomIndices[1]];
    r1n2 = coordsXYZ[atomIndices[1] + 2];
    double p1n3 = coordsXYZ[atomIndices[2]];
    double q1n3 = coordsXYZ[atomIndices[2] + 1];
    double r1n3 = coordsXYZ[atomIndices[2] + 2];

    double sigma0 = acos(q1n3 / sqrt(q1n3 * q1n3 + r1n3 * r1n3));
    if (r1n3 < 0.0) {
      sigma0 = -sigma0;
    }

    for (int i = 0; i < numCoords; i++) {
      int atomIndex = i * 3;
      double qtmp = coordsXYZ[atomIndex + 1] * cos(sigma0) + coordsXYZ[atomIndex + 2] * sin(sigma0);
      double rtmp = -coordsXYZ[atomIndex + 1] * sin(sigma0) + coordsXYZ[atomIndex + 2] * cos(sigma0);
      coordsXYZ[atomIndex + 1] = qtmp;
      coordsXYZ[atomIndex + 2] = rtmp;
    }

    q1n2 = coordsXYZ[atomIndices[1] + 1];
    r1n2 = coordsXYZ[atomIndices[1] + 2];
    q1n3 = coordsXYZ[atomIndices[2] + 1];
    r1n3 = coordsXYZ[atomIndices[2] + 2];

    logger.finer(
        format(" DONE N1: %16.15f %16.15f %16.15f", atomCoords[0], atomCoords[1], atomCoords[2]));
    logger.finer(format(" DONE N2: %16.15f %16.15f %16.15f", p1n2, q1n2, r1n2));
    logger.finer(format(" DONE N3: %16.15f %16.15f %16.15f", p1n3, q1n3, r1n3));
  }
}
