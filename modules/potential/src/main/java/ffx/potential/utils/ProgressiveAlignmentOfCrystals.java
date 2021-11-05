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
import static ffx.potential.parsers.SystemFilter.setVersioning;
import static ffx.potential.parsers.SystemFilter.Versioning.PREFIX;
import static ffx.potential.parsers.SystemFilter.version;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static java.util.Arrays.sort;
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
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
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
   * Specify the number of digits to round intermediate comparisons.
   */
  private static final int ROUND_PLACES = 10;
  /**
   * If molecules between two crystals differ below this tolerance, it is assumed they are equivalent.
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
    setVersioning(PREFIX);
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

    // Initialize array as -1.0 as -1.0 is not a viable RMSD.
    fill(distRow, -1.0);

    // Each process will complete the following amount of work per row.
    distances = new double[numProc][numWorkItems];

    // Initialize each distance as -1.0.
    for (int i = 0; i < numProc; i++) {
      fill(distances[i], -1.0);
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
   * Compare the crystals within the SystemFilters that were inputted into the constructor of this
   * class.
   * @param nAU Number of asymmetric units to compare.
   * @param inflatedAU Minimum number of asymmetric units in inflated crystal
   * @param numSearch Number of loops to search for stereoisomers in first system.
   * @param numSearch2 Number of loops to search for stereoisomers in second system.
   * @param zPrime Number of asymmetric units in first crystal.
   * @param zPrime2 Number of asymmetric units in second crystal.
   * @param alphaCarbons Perform comparisons on only alpha carbons.
   * @param noHydrogen Perform comparisons without hydrogen atoms.
   * @param exhaustive Perform exhaustive number of comparisons (takes longer, but more information
   *     returned).
   * @param save Save out PDBs of the resulting superposition.
   * @param restart Try to restart from a previous job.
   * @param write Save out a PAC RMSD file.
   * @param machineLearning Save out CSV files for machine learning input (saves PDBs as well).
   * @param pacFileName The filename to use.
   * @return RunningStatistics Statistics for comparisons performed.
   */
  public RunningStatistics comparisons(int nAU, int inflatedAU, int numSearch, int numSearch2,
      int zPrime, int zPrime2, boolean alphaCarbons, boolean noHydrogen, boolean exhaustive, boolean save, boolean restart,
      boolean write, boolean machineLearning, String pacFileName) {

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

    // If Sohncke group must respect chirality (assuming mono-isomer AU) so numSearch = 1
    MolecularAssembly baseAssembly = baseFilter.getActiveMolecularSystem();
    Crystal baseCrystal = baseAssembly.getCrystal().getUnitCell();
    if (!exhaustive && baseCrystal.spaceGroup.respectsChirality() && zPrime < 2) {
      numSearch = 1;
    }
    Crystal targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
    if (!exhaustive && targetCrystal.spaceGroup.respectsChirality() && zPrime2 < 2) {
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

    // Atom arrays from the 1st assembly.
    Atom[] baseAtoms = baseAssembly.getAtomArray();
    int nAtoms = baseAtoms.length;

    // Collect selected atoms.
    ArrayList<Integer> atomList = new ArrayList<>();
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = baseAtoms[i];
      if (atom.isActive()) {
        String atomName = atom.getName();
        int atomAtNum = atom.getAtomicNumber();
        boolean proteinCheck = atomName.equals("CA") && atomAtNum == 6;
        boolean aminoCheck = (atomName.equals("N1") || atomName.equals("N9")) && atomAtNum == 7;
        if(alphaCarbons){
          if(proteinCheck||aminoCheck){
            atomList.add(i);
          }
        }else if (!noHydrogen || !atom.isHydrogen()) {
          atomList.add(i);
        }
      }
      // Reset all atoms to active once the selection is recorded.
      atom.setActive(true);
    }

    if (atomList.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD in first crystal.");
      return null;
    }

    // Atom arrays from the 2nd assembly.
    MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
    Atom[] targetAtoms = targetAssembly.getAtomArray();
    int nAtoms2 = targetAtoms.length;

    // Collect selected atoms.
    ArrayList<Integer> atomList2 = new ArrayList<>();
    for (int i = 0; i < nAtoms2; i++) {
      Atom atom = targetAtoms[i];
      if (atom.isActive()) {
        String atomName = atom.getName();
        int atomAtNum = atom.getAtomicNumber();
        boolean proteinCheck = atomName.equals("CA") && atomAtNum == 6;
        boolean aminoCheck = (atomName.equals("N1") || atomName.equals("N9")) && atomAtNum == 7;
        if(alphaCarbons){
          if(proteinCheck||aminoCheck){
            atomList2.add(i);
          }
        }else if (!noHydrogen || !atom.isHydrogen()) {
          atomList2.add(i);
        }
      }
      // Reset all atoms to active once the selection is recorded.
      atom.setActive(true);
    }
    if (atomList2.size() < 1) {
      logger.info("\n No atoms were selected for the PAC RMSD in second crystal.");
      return null;
    }
    int[] comparisonAtoms = atomList.stream().mapToInt(i->i).toArray();
    int[] comparisonAtoms2 = atomList2.stream().mapToInt(i->i).toArray();

    int compareAtomsSize = comparisonAtoms.length;
    int compareAtomsSize2 = comparisonAtoms2.length;

    //Determine number of species within asymmetric unit.
    int z1 = guessZPrime(zPrime, baseAssembly.getMolecules().size());
    int z2 = guessZPrime(zPrime2, targetAssembly.getMolecules().size());
    // Each ASU contains z * comparisonAtoms species so treat each species individually.
    if (z1 > 1) {
      compareAtomsSize /= z1;
    }
    if (z2 > 1) {
      compareAtomsSize2 /= z2;
    }

    if(compareAtomsSize != compareAtomsSize2){
      logger.warning("Selected atom sizes differ between crystals.");
    }

    // Number of atoms included in the PAC RMSD.
    logger.info(format("\n %d atoms will be used for the PAC RMSD out of %d in first crystal.", compareAtomsSize * z1, nAtoms));
    logger.info(format(" %d atoms will be used for the PAC RMSD out of %d in second crystal.\n", compareAtomsSize2 * z2, nAtoms2));

    // Label for logging.
    rmsdLabel = format("RMSD_%d", nAU);

    // Loop over conformations in the base assembly.
    for (int row = restartRow; row < baseSize; row++) {
      // Initialize the distance this rank is responsible for to zero.
      fill(myDistances, -1.0);
      int myIndex = 0;
      // Base unit cell for logging.
      baseCrystal = baseFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
      for (int column = restartColumn; column < targetSize; column++) {
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
            rmsd = compare(comparisonAtoms, comparisonAtoms2, compareAtomsSize, compareAtomsSize2, nAU, inflatedAU,
                numSearch, numSearch2, z1, z2, exhaustive, save, machineLearning);
            time += System.nanoTime();
            double timeSec = time * 1.0e-9;
            // Record the fastest comparison.
            if (timeSec < minTime) {
              minTime = timeSec;
            }
            // Log the final result.
            logger.info(format(" PAC %s: %12s %7.4f A (%5.3f sec)", rmsdLabel, "", rmsd, timeSec));
          }
          myDistances[myIndex] = rmsd;
          myIndex++;
        }
        targetFilter.readNext(false, false);
      }

      restartColumn = 0;
      targetFilter.readNext(true, false);
      baseFilter.readNext(false, false);

      // Gather RMSDs for this row.
      gatherRMSDs(row, runningStatistics);

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
   * @param comparisonAtomsStart List of indices for active atoms in first crystal.
   * @param comparisonAtoms2Start List of indices for active atoms in second crystal.
   * @param compareAtomsSize Number of active atoms in asymmetric unit of first crystal.
   * @param compareAtomsSize2 Number of active atoms in asymmetric unit of second crystal.
   * @param nAU Number of asymmetric units to compare.
   * @param inflatedAU Number of asymmetric units in expanded system.
   * @param numSearch Number of molecules to check for conformations in the first system.
   * @param numSearch2 Number of molecules to check for conformations in the second system.
   * @param exhaustive Perform exhaustive number of comparisons.
   * @param save Save out PDBs of compared crystals.
   * @param machineLearning Save out PDBs and CSVs of compared crystals.
   * @return the computed RMSD.
   */
  private double compare(int[] comparisonAtomsStart, int[] comparisonAtoms2Start, int compareAtomsSize,
      int compareAtomsSize2, int nAU, int inflatedAU, int numSearch, int numSearch2, int zPrime, int zPrime2,
      boolean exhaustive, boolean save, boolean machineLearning) {
    // TODO: Does PAC work for a combination of molecules and polymers?

    //Remove atoms not used in comparisons from the original molecular assembly (crystal 1).
    MolecularAssembly baseAssembly = baseFilter.getActiveMolecularSystem();
    baseAssembly.moveAllIntoUnitCell();
    double[] massStart = new double[comparisonAtomsStart.length];
    double[] reducedBaseCoords = reduceSystem(baseAssembly, comparisonAtomsStart, massStart);
    Crystal baseXtal = baseAssembly.getCrystal();
    double[] baseXYZOrig = generateInflatedSphere(baseXtal, reducedBaseCoords, massStart,
            inflatedAU);
    //Remove atoms not used in comparisons from the original molecular assembly (crystal 2).
    MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
    targetAssembly.moveAllIntoUnitCell();
    double[] massStart2 = new double[comparisonAtoms2Start.length];
    double[] reducedTargetCoords = reduceSystem(targetAssembly, comparisonAtoms2Start, massStart2);
    Crystal targetXtal = targetAssembly.getCrystal();
    double[] targetXYZOrig = generateInflatedSphere(targetXtal, reducedTargetCoords, massStart2,
            inflatedAU);

    double[] mass = new double[compareAtomsSize];
    double[] massTemp = new double[compareAtomsSize2];
    arraycopy(massStart, 0, mass, 0, compareAtomsSize);
    arraycopy(massStart2, 0, massTemp, 0, compareAtomsSize2);
    if (!Arrays.equals(mass, massTemp)) {
      if (logger.isLoggable(Level.FINER)) {
        for (int i = 0; i < compareAtomsSize; i++) {
          logger.finer(format("Mass: %4.4f MassTemp: %4.4f", mass[i], massTemp[i]));
        }
      }
      logger.warning("Atom masses are not equivalent between crystals.\n " +
              "Ensure atom ordering is same in both inputs.");
    }
    // Remove duplicated atoms from Z' > 1.
    int[] comparisonAtoms = new int[compareAtomsSize];
    arraycopy(comparisonAtomsStart, 0, comparisonAtoms, 0, compareAtomsSize);
    int[] comparisonAtoms2 = new int[compareAtomsSize2];
    arraycopy(comparisonAtoms2Start, 0, comparisonAtoms2, 0, compareAtomsSize2);

    if (!Arrays.equals(comparisonAtoms, comparisonAtoms2)) {
      if (logger.isLoggable(Level.FINER)) {
        for (int i = 0; i < compareAtomsSize; i++) {
          logger.finer(format("compAtoms: %d compAtoms: %d", comparisonAtoms[i], comparisonAtoms2[i]));
        }
      }
      logger.warning("Atoms to compare are not equivalent between crystals.\n " +
              "Ensure atom ordering is same in both inputs.");
    }

    // Number of used coordinates for atoms in one AU.
    int nCoords = compareAtomsSize * 3;
    //Number of species in expanded crystals.
    int nBaseMols = baseXYZOrig.length / nCoords;
    int nTargetMols = targetXYZOrig.length / nCoords;
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" Number of copies to compare:    %4d", nAU));
      logger.finer(format(" Number entities in base sphere: %4d", nBaseMols));
      logger.finer(format(" Number entities in target sphere: %d", nTargetMols));
    }
    //Mass of 3 species
    double[] mass3 = new double[nCoords];
    for (int i = 0; i < 3; i++) {
      arraycopy(mass, 0, mass3, i * compareAtomsSize, compareAtomsSize);
    }
    //Mass of N species
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
    // Sum of all masses for one species.
    double massSum = 0;
    for (double value : mass) {
      massSum += value;
    }
    // Center of Masses for crystal 1.
    double[][] baseCoMOrig = new double[nBaseMols][3];
    centerOfMass(baseCoMOrig, baseXYZOrig, mass, massSum);
    prioritizeReplicates(compareAtomsSize, baseXYZOrig, mass, massSum, nBaseMols, baseCoMOrig,
            molDist1);

    //Used for debugging. can be removed.
    if (logger.isLoggable(Level.FINEST)) {
      int printSize = 20;
      logger.finest(" System 1 distances to center of sphere:");
      for (int i = 0; i < printSize; i++) {
        logger.finest(format(" %d\t%16.8f", molDist1[i].getIndex(), molDist1[i].getDoubleValue()));
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
    prioritizeReplicates(compareAtomsSize2, targetXYZOrig, mass, massSum, nTargetMols, targetCoMOrig,
            molDist2);

    if (logger.isLoggable(Level.FINEST)) {
      int printSize = 20;
      logger.finest(" System 2 distances to center of sphere:");
      for (int i = 0; i < printSize; i++) {
        logger.finest(format(" %d\t%16.8f", molDist2[i].getIndex(), molDist2[i].getDoubleValue()));
      }
    }

    //  Frank-Kasper phases of metallic ions can reach a coordination number of 16...
    //Determine if AUs in first system are same hand as center most in first (stereoisomer handling).
    logger.finer(" Base Search Conformations:");
    ArrayList<double[]> base0XYZ = new ArrayList<>();
    double[] tempBase = new double[nCoords];
    arraycopy(baseXYZOrig, molDist1[0].getIndex() * nCoords, tempBase, 0, nCoords);
    List<Integer> baseindices = new ArrayList<>();
    List<Double> basediff = new ArrayList<>();
    baseindices.add(0);
    basediff.add(rmsd(tempBase, tempBase, mass));
    base0XYZ.add(tempBase);
    int baseSearchValue = (baseXtal.spaceGroup.respectsChirality()) ? zPrime : 2 * zPrime;
    // Start from 1 as zero is automatically added.
    int index = 1;

    //Determine number of conformations in first crystal
    // TODO: change below to while(baseDiff.size() < baseSearchValue) when convinced
    int numConfCheck = (exhaustive) ? nBaseMols : baseSearchValue;
    while (basediff.size() < numConfCheck) {
      double[] baseCheckMol = new double[nCoords];
      arraycopy(baseXYZOrig, molDist1[index].getIndex() * nCoords, baseCheckMol, 0,
              nCoords);
      double value = doMoleculesMatch(base0XYZ.get(0), baseCheckMol, mass);
      value = roundToPlaces(value);
      if (addLooseUnequal(basediff, value)) {
        baseindices.add(index);
        base0XYZ.add(baseCheckMol);
      } else if (Collections.frequency(basediff, value) < numSearch) {
        basediff.add(value);
        baseindices.add(index);
        base0XYZ.add(baseCheckMol);
      }
      index++;
      if (index >= molDist1.length) {
        break;
      }
    }
    double[] baseDiff = basediff.stream().mapToDouble(i -> i).toArray();
    double[][] baseXYZs = new double[base0XYZ.size()][nCoords];
    base0XYZ.toArray(baseXYZs);
    Integer[] baseIndices = new Integer[baseindices.size()];
    baseindices.toArray(baseIndices);
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" %d conformations detected out of %d in base crystal.", baseIndices.length, baseSearchValue * numSearch));
      logger.finer(" Target Search Conformations:");
    }
    int targetSearchValue = (targetXtal.spaceGroup.respectsChirality()) ? zPrime2 : 2 * zPrime2;

    ArrayList<double[]> target0XYZ = new ArrayList<>();
    List<Integer> targetindices = new ArrayList<>();
    List<Integer> targetBaseindices = new ArrayList<>();
    List<Double> targetdiff = new ArrayList<>();
    // Determine number of conformation in second.
    index = 0;
    // TODO: change below to while(targetDiff.size() < targetSearchValue) when convinced
    numConfCheck = (exhaustive)? nTargetMols : targetSearchValue;
    while (targetdiff.size() < numConfCheck) {
      double[] targetCheckMol = new double[nCoords];
      arraycopy(targetXYZOrig, molDist2[index].getIndex() * nCoords, targetCheckMol, 0, nCoords);
      int minIndex = -1;
      double minDiff = Double.MAX_VALUE;
      for (int i = 0; i < baseIndices.length; i++) {
        double value = doMoleculesMatch(baseXYZs[i], targetCheckMol, mass);
        value = roundToPlaces(value);
        if (value < minDiff) {
          minDiff = value;
          minIndex = i;
        }
      }
      if (!targetBaseindices.contains(baseIndices[minIndex])) {
        targetdiff.add(minDiff);
        targetBaseindices.add(baseIndices[minIndex]);
        targetindices.add(index);
        target0XYZ.add(targetCheckMol);
      } else if (Collections.frequency(targetBaseindices, baseIndices[minIndex]) < numSearch2) {
        targetBaseindices.add(baseIndices[minIndex]);
        targetdiff.add(minDiff);
        targetindices.add(index);
        target0XYZ.add(targetCheckMol);
      }
      index++;
      if (index >= molDist2.length) {
        break;
      }
    }
    double[] targetDiff = targetdiff.stream().mapToDouble(i -> i).toArray();
    double[][] targetXYZs = new double[target0XYZ.size()][nCoords];
    target0XYZ.toArray(targetXYZs);
    Integer[] targetIndices = new Integer[targetindices.size()];
    targetindices.toArray(targetIndices);
    Integer[] targetBaseIndices = new Integer[targetBaseindices.size()];
    targetBaseindices.toArray(targetBaseIndices);
    if (logger.isLoggable(Level.FINER)) {
      logger.finer(format(" %d conformations detected out of %d in target crystal.", targetIndices.length, targetSearchValue * numSearch2));
      logger.finer(" Base Matches:");
      for (int i = 0; i < baseIndices.length; i++) {
        logger.finer(format(" %d %4.4f %d", i, baseDiff[i], baseIndices[i]));
      }
      logger.finer(" Target Matches:");
      for (int i = 0; i < targetBaseIndices.length; i++) {
        logger.finer(format(" %d %4.4f %d", i, targetDiff[i], targetBaseIndices[i]));
      }
    }

    // Coordinate arrays to save out structures at the end.
    double[] n1TargetNMols = new double[nAU * nCoords];
    double[] n3TargetNMols = new double[nAU * nCoords];
    double[] bestBaseNMols = new double[nAU * nCoords];
    double[] bestTargetNMols = new double[nAU * nCoords];
    double bestRMSD = Double.MAX_VALUE;
    List<Double> crystalCheckRMSDs = new ArrayList<>();

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format("\n  Trial     RMSD_1 (%7s)  RMSD_3 (%7s)  %7s", rmsdLabel, rmsdLabel, rmsdLabel));
//      logger.fine(format("\n Trial      RMSD_1  RMSD_3 %7s", rmsdLabel));
    }
    // Begin comparison
    for (int l = 0; l < baseIndices.length; l++) {
      // Place rest of comparison code here.
      double[] baseXYZ = new double[nBaseMols * nCoords];
      double[][] baseCoM = new double[nBaseMols][3];
      arraycopy(baseXYZOrig, 0, baseXYZ, 0, baseXYZOrig.length);
      arraycopy(baseCoMOrig, 0, baseCoM, 0, baseCoMOrig.length);
      int center = molDist1[baseIndices[l]].getIndex();
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
      logger.finer(" Base 3 Conformations");
      double[] base3Mols = new double[nCoords * 3];
      for (int i = 0; i < 3; i++) {
        int baseIndex = molDist1_2[i].getIndex() * nCoords;
        arraycopy(baseXYZ, baseIndex, base3Mols, i * nCoords, nCoords);
      }
      int[] matchB3 = new int[3];
      if (logger.isLoggable(Level.FINER)) {
        for (int i = 0; i < 3; i++) {
          double[] base1of3Mols = new double[nCoords];
          arraycopy(base3Mols, i * nCoords, base1of3Mols, 0, nCoords);
          double minDiff = Double.MAX_VALUE;
          int minIndex = -1;
          for (int j = 0; j < baseIndices.length; j++) {
            double value = doMoleculesMatch(baseXYZs[j], base1of3Mols, mass);
            if (value < minDiff) {
              minDiff = value;
              minIndex = baseIndices[j];
            }
          }
          matchB3[i] = minIndex;
          logger.finer(format(" %d %4.4f %d", i, minDiff, matchB3[i]));
        }
      }

      // Acquire coordinates for final comparison
      double[] baseNMols = new double[nAU * nCoords];
      for (int i = 0; i < nAU; i++) {
        int molIndex = molDist1_2[i].getIndex() * nCoords;
        arraycopy(baseXYZ, molIndex, baseNMols, i * nCoords, nCoords);
      }
      int[] matchBN = new int[nAU];
      if (logger.isLoggable(Level.FINER)) {
        logger.finer(" Base N Conformations");
        for (int i = 0; i < nAU; i++) {
          double[] base1ofNMols = new double[nCoords];
          arraycopy(baseNMols, i * nCoords, base1ofNMols, 0, nCoords);
          double minDiff = Double.MAX_VALUE;
          int minIndex = -1;
          for (int j = 0; j < baseIndices.length; j++) {
            double value = doMoleculesMatch(baseXYZs[j], base1ofNMols, mass);
            if (value < minDiff) {
              minDiff = value;
              minIndex = baseIndices[j];
            }
          }
          matchBN[i] = minIndex;
          logger.finer(format(" %d %4.4f %d", i, minDiff, matchBN[i]));
        }
      }
      for (int m = 0; m < targetXYZs.length; m++) {
        if (exhaustive || Objects.equals(targetBaseIndices[m], baseIndices[l])) {
          double[] targetXYZ = new double[nTargetMols * nCoords];
          double[][] targetCoM = new double[nTargetMols][3];
          arraycopy(targetXYZOrig, 0, targetXYZ, 0, targetXYZOrig.length);
          arraycopy(targetCoMOrig, 0, targetCoM, 0, targetCoMOrig.length);
          // Switch m center most molecules (looking for stereoisomers)
          center = molDist2[targetIndices[m]].getIndex();
          DoubleIndexPair[] molDist2_2 = new DoubleIndexPair[nTargetMols];

          //Re-prioritize based on center most molecule.
          if (logger.isLoggable(Level.FINER)) {
            logger.finer(" Re-prioritize target system.");
          }
          prioritizeReplicates(compareAtomsSize2, targetXYZ, mass, massSum, nTargetMols, targetCoM,
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
            int numCheck = Math.min(5, nAU);
            double[][] baseCenters = new double[numCheck][3];
            double[][] targetCenters = new double[numCheck][3];
            for (int j = 0; j < numCheck; j++) {
              int baseCenterMols = molDist1_2[j].getIndex() * nCoords;
              int targetCenterMols = molDist2_2[j].getIndex() * nCoords;
              for (int i = 0; i < Math.min(compareAtomsSize, compareAtomsSize2); i++) {
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
          // TODO could have prioritization favor closest distance from molDist12 or molDist22 (was symmetric flag).
          if (logger.isLoggable(Level.FINER)) {
            logger.finer(" Match molecules between systems");
          }

          //Update center of masses with the first trans/rot
          centerOfMass(targetCoM, targetXYZ, mass, massSum);

          //TODO replace with commented out code. Increased for figure generation
//            DoubleIndexPair[] matchMols = new DoubleIndexPair[3];
          DoubleIndexPair[] matchMols = new DoubleIndexPair[Math.max(3, nAU)];
          matchMolecules(matchMols, baseCoM, targetCoM, molDist1_2, molDist2_2);

          //TODO remove following?
          for (int i = 0; i < nAU; i++) {
            int offset = i * nCoords;
            int molIndex = matchMols[i].getIndex() * nCoords;
            arraycopy(targetXYZ, molIndex, n1TargetNMols, offset, nCoords);
          }
          double n1RMSD = rmsd(baseNMols, n1TargetNMols, massN);
          if (logger.isLoggable(Level.FINEST)) {
            logger.finest("  Distance between pairs after rot 1:");
            for (DoubleIndexPair matchMol : matchMols) {
              logger.finest(
                      format(" %2d %16.8f", matchMol.getIndex(), matchMol.getDoubleValue()));
            }
            logger.finest(" Index  MolDist1    MolDist2    MatchMols");
            for (int i = 0; i < matchMols.length; i++) {
              logger.finest(format(" %2d base: %2d target: %2d Match: %2d", i,
                      molDist1_2[i].getIndex(), molDist2_2[i].getIndex(), matchMols[i].getIndex()));
            }
          }
          logger.finer(" Rotation 2:");

          secondRotation(base3Mols, targetXYZ, mass3, matchMols);
          double[] target3Mols = new double[nCoords * 3];
          for (int i = 0; i < 3; i++) {
            int molIndex = i * nCoords;
            targetCenterMol = matchMols[i].getIndex() * nCoords;
            arraycopy(targetXYZ, targetCenterMol, target3Mols, molIndex, nCoords);
          }
          double checkRMSD2 = rmsd(base3Mols, target3Mols, mass3);
          for (int i = 0; i < nAU; i++) {
            int offset = i * nCoords;
            int molIndex = matchMols[i].getIndex() * nCoords;
            arraycopy(targetXYZ, molIndex, n3TargetNMols, offset, nCoords);
          }
          double n3RMSD = rmsd(baseNMols, n3TargetNMols, massN);

          int[] matchT3 = new int[3];
          boolean hand3Match = true;
          if (logger.isLoggable(Level.FINER)) {
            logger.finer(" Target 3 Conformations");
            for (int i = 0; i < 3; i++) {
              double[] target1of3Mols = new double[nCoords];
              arraycopy(target3Mols, i * nCoords, target1of3Mols, 0, nCoords);
              double minDiff = Double.MAX_VALUE;
              int minIndex = -1;
              for (int j = 0; j < baseIndices.length; j++) {
                double value = doMoleculesMatch(baseXYZs[j], target1of3Mols, mass);
                if (value < minDiff) {
                  minDiff = value;
                  minIndex = baseIndices[j];
                }
              }
              matchT3[i] = minIndex;
              logger.finer(format(" %d %4.4f %d", i, minDiff, matchT3[i]));
              if (matchB3[i] != matchT3[i]) {
                hand3Match = false;
              }
            }
          }

          //Update center of masses with the second rot (only one crystal moves).
          centerOfMass(targetCoM, targetXYZ, mass, massSum);
          // Rotations 1 and 2 have been completed and both systems should be overlapped
          //  Isolate center most nAU from System 1 and matching molecules from System 2
          if (logger.isLoggable(Level.FINER)) {
            logger.finer(" Match Molecules:");
          }

          matchMols = new DoubleIndexPair[nAU];
          matchMolecules(matchMols, baseCoM, targetCoM, molDist1_2, molDist2_2);

          if (logger.isLoggable(Level.FINEST)) {
            int printSize = Math.min(nAU, 10);
            logger.finest(format("  Distance between %d pairs after rot 2:", printSize));
            for (int i = 0; i < printSize; i++) {
              logger.finest(
                      format(" %2d %16.8f", matchMols[i].getIndex(), matchMols[i].getDoubleValue()));
            }
          }

          double[] targetNMols = new double[nAU * nCoords];
          for (int i = 0; i < nAU; i++) {
            int offset = i * nCoords;
            int molIndex = matchMols[i].getIndex() * nCoords;
            arraycopy(targetXYZ, molIndex, targetNMols, offset, nCoords);
          }

          int[] matchTN = new int[nAU];
          boolean handNMatch = true;
          if (logger.isLoggable(Level.FINER)) {
            logger.finer(" Target N Conformations");

            for (int i = 0; i < nAU; i++) {
              double[] target1ofNMols = new double[nCoords];
              arraycopy(targetNMols, i * nCoords, target1ofNMols, 0, nCoords);
              double minDiff = Double.MAX_VALUE;
              int minIndex = -1;
              for (int j = 0; j < baseIndices.length; j++) {
                double value = doMoleculesMatch(baseXYZs[j], target1ofNMols, mass);
                if (value < minDiff) {
                  minDiff = value;
                  minIndex = baseIndices[j];
                }
              }
              matchTN[i] = minIndex;
              logger.finer(format(" %d %4.12f %d", i, minDiff, matchTN[i]));
              if (matchBN[i] != matchTN[i]) {
                handNMatch = false;
              }
            }
            logger.finer(" Final rotation:");
          }

          translate(baseNMols, massN, targetNMols, massN);
          rotate(baseNMols, targetNMols, massN);

          double rmsdSymOp = rmsd(baseNMols, targetNMols, massN);
          if (logger.isLoggable(Level.FINE)) {
            String output = format(" %2d of %2d: %7.4f (%7.4f) %7.4f (%7.4f) %7.4f",
                    l * targetSearchValue * numSearch2 + m + 1,
                    baseSearchValue * numSearch * targetSearchValue * numSearch2, checkRMSD1, n1RMSD, checkRMSD2, n3RMSD, rmsdSymOp);

            if (logger.isLoggable(Level.FINER)) {
              //TODO fix conformation comparisons (target system find first available match base has individual
              boolean conformationMatches =
                      Objects.equals(baseIndices[l], targetBaseIndices[m]) && hand3Match
                              && handNMatch;
              output += format(" %d=%d", baseIndices[l], targetIndices[m]);

              if (save) {
                int loop = l * targetSearchValue + m + 1;
                saveAssemblyPDB(baseAssembly, baseNMols, comparisonAtoms, "_" + loop + "_c1");
                saveAssemblyPDB(targetAssembly, targetNMols, comparisonAtoms, "_" + loop + "_c2");
              }
            }
            logger.fine(output);
          }

          if (rmsdSymOp < bestRMSD) {
            bestRMSD = rmsdSymOp;
            bestBaseNMols = baseNMols;
            bestTargetNMols = targetNMols;
          }
          addLooseUnequal(crystalCheckRMSDs, rmsdSymOp);
        }
      }
    }

    double finalRMSD = Double.NaN;
    if (bestRMSD < Double.MAX_VALUE) {
      finalRMSD = bestRMSD;
    } else {
      logger.warning(" This RMSD was filtered out! Try the --ex flag." +
              "\nAlternatively increase --ns and/or --ns2.");
      // TODO: Double.NaN causes an error in RunningStatistics... Set to -2.0 for now...
      finalRMSD = -2.0;
    }
    if (save) {
      if (machineLearning) {
        saveAssemblyPDB(baseAssembly, bestBaseNMols, comparisonAtoms, "_r" + rank + "_c1", 0.000);
        saveAssemblyPDB(targetAssembly, bestTargetNMols, comparisonAtoms, "_r" + rank + "_c2", finalRMSD);
      } else {
        saveAssemblyPDB(baseAssembly, bestBaseNMols, comparisonAtoms, "_r" + rank + "_c1");
        saveAssemblyPDB(targetAssembly, bestTargetNMols, comparisonAtoms, "_r" + rank + "_c2");
      }
    }

    // Logging to check number of RMSD values determined.
    StringBuilder dblOut = new StringBuilder();
    Double[] uniqueRMSDs = new Double[crystalCheckRMSDs.size()];
    crystalCheckRMSDs.toArray(uniqueRMSDs);
    sort(uniqueRMSDs);
    for (double dbl : uniqueRMSDs) {
      dblOut.append(" ").append(format("%4.4f", dbl));
    }
    String message = format(" Unique %s Values: %s", rmsdLabel, dblOut);

    int numUnique = crystalCheckRMSDs.size();
    if (!exhaustive && numUnique > (baseSearchValue * targetSearchValue) / 2 || exhaustive &&
            numUnique > baseSearchValue * targetSearchValue) {
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
   * @param runningStatistics Stats for the RMSDs.
   */
  private void gatherRMSDs(int row, RunningStatistics runningStatistics) {
    if (useMPI) {
      try {
        if (logger.isLoggable(Level.FINER)) {
          logger.finer(" Receiving results.");
        }
        world.allGather(myBuffer, buffers);
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
                logger.finer(format(" %d %d %16.8f", row, column, distances[proc][workItem]));
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
   * Try to automatically determine number of species in asymmetric unit (only works for molecules).
   * @param zPrime User input overrides detection method.
   * @param numSpecies Number of species detected.
   * @return Number of expected species in asymmetric unit.
   */
  private static int guessZPrime(int zPrime, int numSpecies){
    int z = (zPrime > 0) ? zPrime : Math.max(numSpecies, 0);
    if(z<1){
      logger.warning(
              " Number of species in asymmetric unit was not determined.\n" +
                      "Setting Z'=1. Use --zp/--zp2 flags to set manually.");
      z = 1;
    }
    logger.finer(format(" Number of species in asymmetric unit (Z'): %d", z));
    return z;
  }

  /**
   * Round a double value to ROUND_PLACES decimal points.
   * @param value Double to be rounded.
   * @return Rounded double.
   */
  private static double roundToPlaces(double value){
    BigDecimal bigDecimal = new BigDecimal(Double.toString(value));
    bigDecimal = bigDecimal.setScale(ROUND_PLACES, RoundingMode.HALF_UP);
    return bigDecimal.doubleValue();
  }

  /**
   * Perform first rotation to match center most molecules.
   *
   * @param targetXYZ Coordinates for second system.
   * @param mass Masses for atoms in a molecule.
   * @param baseXYZ Coordinates for first system.
   * @param index2 Index of
   */
  private static void firstRotation(double[] targetXYZ, double[] mass, double[] baseXYZ, int index2) {
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
    double[][] rotation = calculateRotation(baseXYZ, targetMol, mass);
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

  /**
   * Pair species between two crystals based on center of mass distances.
   * @param matchMols Mapping from base crystal to target.
   * @param baseCoM Center of Masses for base crystal.
   * @param targetCoM Center of Masses for target crystal.
   * @param molDist12 Prioritization of base crystal.
   * @param molDist22 Prioritization of target crystal.
   */
  private static void matchMolecules(DoubleIndexPair[] matchMols, double[][] baseCoM, double[][] targetCoM,
      DoubleIndexPair[] molDist12, DoubleIndexPair[] molDist22) {
    int desiredMols = matchMols.length;
    int nTargetMols = targetCoM.length;
    // List of indexes for second system.
    List<Integer> targetIndex = new ArrayList<>(nTargetMols);
    for (DoubleIndexPair doubleIndexPair : molDist22) {
      // Only search molecules within range of the desired number of molecules.
      // Must have enough molecules for matching (using exhaustive till better heuristic is determined)
      targetIndex.add(doubleIndexPair.getIndex());
    }

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
      if (logger.isLoggable(Level.FINEST)) {
        logger.finest(
            format(" Base position:   %d: %8.4f %8.4f %8.4f", i, baseCoM[molDist12[i].getIndex()][0],
                baseCoM[molDist12[i].getIndex()][1],
                baseCoM[molDist12[i].getIndex()][2]));
        logger.finest(format(" Match Distance:  %d: %8.4f", i, matchMols[i].getDoubleValue()));
        logger.finest(format(" Target position: %d: %8.4f %8.4f %8.4f", i,
            targetCoM[molDist22[i].getIndex()][0],
            targetCoM[molDist22[i].getIndex()][1], targetCoM[molDist22[i].getIndex()][2]));
      }
    }
  }

  /**
   * Determine if two species are of the same conformation (should be same species).
   * @param baseXYZ One species.
   * @param targetXYZ Second species.
   * @param mass Mass of atoms within species
   * @return RMSD between both species.
   */
  private static double doMoleculesMatch(double[] baseXYZ, double[] targetXYZ, double[] mass) {
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
   * @param mass Masses of each atom in asymmetric unit.
   * @param massSum Sum of masses within asymmetric unit.
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

  /**
   * Reduce asymmetric unit to atoms that are going to be used in final RMSD.
   * @param assembly Asymmetric unit we wish to reduce.
   * @param comparisonAtoms Atoms of interest within asymmetric unit.
   * @param mass Mass of atoms within asymmetric unit (filling values).
   * @return Linear coordinates for only atoms of interest.
   */
  private static double[] reduceSystem(MolecularAssembly assembly, int[] comparisonAtoms,
      double[] mass) {
    Atom[] atoms = assembly.getAtomArray();
    // Collect asymmetric unit atomic coordinates.
    double[] reducedCoords = new double[comparisonAtoms.length * 3];
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
   * @param reducedCoords Coordinates of asymmetric unit we wish to expand.
   * @param mass Masses for atoms within reduced asymmetric unit.
   * @param inflatedAU Number of asymmetric units after inflation.
   * @return double[] containing the coordinates for the expanded crystal.
   */
  private static double[] generateInflatedSphere(Crystal crystal, double[] reducedCoords,
      double[] mass, int inflatedAU) {
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
    sort(molsDists);

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
   * @param massSum Sum of atomic masses within asymmetric unit.
   * @param nMols Number of asymmetric units in expanded crystal.
   * @param centerOfMasses Center of masses for each replicate within inflated crystal.
   * @param molDists DoubleIndexPair of distances from center most asymmetric unit to asymmetric
   *     unit at index.
   */
  private static void prioritizeReplicates(int nAtoms, double[] coordsXYZ, double[] mass,
      double massSum, int nMols, double[][] centerOfMasses, DoubleIndexPair[] molDists) {
    // AUs added to system based on distance to center of all atoms. 0th AU should be closest to all atom center.
    prioritizeReplicates(nAtoms, coordsXYZ, mass, massSum, nMols, centerOfMasses, molDists, 0);
  }

  /**
   * Prioritize asymmetric units within the system based on distance to specified index.
   *
   * @param nAtoms Number of atoms being compared in asymmetric unit.
   * @param coordsXYZ Coordinates for expanded crystal.
   * @param mass Mass of atoms within asymmetric unit.
   * @param massSum Sum of atomic masses within asymmetric unit.
   * @param nMols Number of molecules in expanded crystal.
   * @param centerOfMasses Center of masses for each replicate within inflated crystal.
   * @param molDists Prioritization of molecules in expanded system on distance to center-most
   *     molecule.
   * @param index Index of molecules to be center.
   */
  private static void prioritizeReplicates(int nAtoms, double[] coordsXYZ, double[] mass,
      double massSum, int nMols, double[][] centerOfMasses, DoubleIndexPair[] molDists, int index) {
    // Find AU to be treated as the new center.
    double[] coordCenter = centerOfMasses[index];
    for (int i = 0; i < nMols; i++) {
      double[] moleculeCenter = centerOfMasses[i];
      molDists[i] = new DoubleIndexPair(i, dist(coordCenter, moleculeCenter));
    }
    // Reorder based on distance to AU closest to Index.
    sort(molDists);

    if (logger.isLoggable(Level.FINEST)) {
      int numCheck = Math.min(5, molDists.length);
      double[][] targetMol = new double[coordsXYZ.length / (nAtoms * 3)][3];
      centerOfMass(targetMol, coordsXYZ, mass, massSum);
      for (int i = 0; i < numCheck; i++) {
        logger.finest(format(" 1AU Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
            molDists[i].getIndex(),
            molDists[i].getDoubleValue(), targetMol[molDists[i].getIndex()][0],
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
        (centerOfMasses[molDists[0].getIndex()][0] + centerOfMasses[molDists[1].getIndex()][0]) / 2;
    avgCenter[1] =
        (centerOfMasses[molDists[0].getIndex()][1] + centerOfMasses[molDists[1].getIndex()][1]) / 2;
    avgCenter[2] =
        (centerOfMasses[molDists[0].getIndex()][2] + centerOfMasses[molDists[1].getIndex()][2]) / 2;

    for (int i = 1; i < nMols; i++) {
      double[] moleculeCenter = centerOfMasses[molDists[i].getIndex()];
      molDists2[i] = new DoubleIndexPair(molDists[i].getIndex(), dist(avgCenter, moleculeCenter));
    }
    sort(molDists2);
    if (logger.isLoggable(Level.FINEST)) {
      int numCheck = Math.min(5, molDists2.length);
      double[][] targetMol = new double[coordsXYZ.length / (nAtoms * 3)][3];
      centerOfMass(targetMol, coordsXYZ, mass, massSum);
      for (int i = 0; i < numCheck; i++) {
        logger.finest(format(" 2AU Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                molDists2[i].getIndex(),
                molDists2[i].getDoubleValue(), targetMol[molDists2[i].getIndex()][0],
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
      molDists3[i] = new DoubleIndexPair(molDists2[i].getIndex(), dist(avgCenter, moleculeCenter));
    }
    //Reorder based on center point between center-most AU to all atom center and closest AU to center-most AU.
    Arrays.sort(molDists3);
    arraycopy(molDists3, 0, molDists, 0, nMols);
    if (logger.isLoggable(Level.FINEST)) {
      int numCheck = Math.min(5, molDists3.length);
      double[][] targetMol = new double[coordsXYZ.length / (nAtoms * 3)][3];
      centerOfMass(targetMol, coordsXYZ, mass, massSum);
      for (int i = 0; i < numCheck; i++) {
        logger.finest(format(" 3AU Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                molDists3[i].getIndex(),
                molDists3[i].getDoubleValue(), targetMol[molDists3[i].getIndex()][0],
                targetMol[molDists3[i].getIndex()][1], targetMol[molDists3[i].getIndex()][2]));
      }
    }
  }

  /**
   * Add a value to a list of doubles if its difference to all listed values is greater than the
   * tolerance.
   *
   * @param values List of values already found.
   * @param value Potential new value if it is not already in list.
   */
  private static boolean addLooseUnequal(List<Double> values, double value) {
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
    // If value is not found it was added. Otherwise, not added.
    return !found;
  }

  /**
   * Save the provided coordinates as a PDB file.
   * @param molecularAssembly Asymmetric unit that forms the crystal of interest.
   * @param coords Coordinates to be saved within the PDB.
   * @param comparisonAtoms Atoms of interest within the initial asymmetric unit.
   * @param description Unique identifier that will be added to PDB file name.
   */
  private static void saveAssemblyPDB(MolecularAssembly molecularAssembly, double[] coords,
      int[] comparisonAtoms, String description) {
    String fileName = FilenameUtils.removeExtension(molecularAssembly.getName());
    File saveLocationPDB = version(new File(fileName + description + ".pdb"));
    // Save aperiodic system of n_mol closest atoms for visualization
    MolecularAssembly currentAssembly = new MolecularAssembly(molecularAssembly.getName());
    ArrayList<Atom> newAtomList = new ArrayList<>();
    Atom[] atoms = molecularAssembly.getAtomArray();
    int atomIndex = 0;
    int compareAtomsSize = comparisonAtoms.length;
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
   * Save the provided coordinates as a PDB file with accompanying CSV containing RMSD.
   * @param molecularAssembly Asymmetric unit that forms the crystal of interest.
   * @param coords Coordinates to be saved within the PDB.
   * @param comparisonAtoms Atoms of interest within the initial asymmetric unit.
   * @param description Unique identifier that will be added to PDB file name.
   * @param finalRMSD RMSD to be saved to CSV file.
   */
  private static void saveAssemblyPDB(MolecularAssembly molecularAssembly, double[] coords,
                                      int[] comparisonAtoms, String description, double finalRMSD) {
    saveAssemblyPDB(molecularAssembly, coords, comparisonAtoms, description);
    String fileName = FilenameUtils.removeExtension(molecularAssembly.getName());
    try {
      File csv = version(new File(fileName + description + ".csv"));
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