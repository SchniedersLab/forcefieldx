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

import static ffx.potential.parsers.DistanceMatrixFilter.writeDistanceMatrixRow;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.io.FilenameUtils.concat;
import static org.apache.commons.io.FilenameUtils.getBaseName;
import static org.apache.commons.io.FilenameUtils.getFullPath;
import static org.apache.commons.io.FilenameUtils.getPath;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.numerics.math.Double3;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parsers.DistanceMatrixFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import java.io.File;
import java.util.logging.Logger;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

public class Superpose {

  /**
   * Logger for the Superpose Class.
   */
  private static final Logger logger = Logger.getLogger(Superpose.class.getName());

  private final SystemFilter baseFilter;
  private final SystemFilter targetFilter;
  private final boolean isSymmetric;

  private final int baseSize;
  private final int targetSize;

  private int restartRow;
  private int restartColumn;
  private double[][] distMatrix;

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


  public Superpose(SystemFilter baseFilter, SystemFilter targetFilter, boolean isSymmetric) {
    this.baseFilter = baseFilter;
    this.targetFilter = targetFilter;
    this.isSymmetric = isSymmetric;

    // Number of models to be evaluated.
    baseSize = baseFilter.countNumModels();
    targetSize = targetFilter.countNumModels();

    // Initialize the distance matrix.
    distMatrix = new double[baseSize][targetSize];

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
   * This method calculates the all versus all RMSD of a multiple model pdb/arc file.
   *
   * @param usedIndices List of atom indices in use.
   * @param dRMSD Compute dRMSDs in addition to RMSD.
   * @param verbose Verbose logging.
   * @param restart Attempt to restart from an RMSD matrix.
   * @param write Write out an RMSD file.
   * @param saveSnapshots Save superposed snapshots.
   */
  public void calculateRMSDs(int[] usedIndices, boolean dRMSD, boolean verbose, boolean restart,
      boolean write, boolean saveSnapshots) {

    String filename = baseFilter.getFile().getAbsolutePath();

    // Prepare to write out superposed snapshots.
    File targetOutputFile = null;
    SystemFilter targetOutputFilter = null;
    if (saveSnapshots) {
      String targetOutputName = concat(getPath(filename), getBaseName(filename) + "_superposed.arc");
      targetOutputFile = SystemFilter.version(new File(targetOutputName));
      MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
      targetOutputFilter = new XYZFilter(targetOutputFile, targetAssembly,
          targetAssembly.getForceField(),
          targetAssembly.getProperties());
    }

    String matrixFilename = concat(getFullPath(filename), getBaseName(filename) + ".txt");
    if (restart) {
      // Define the filename to use for the RMSD values.
      readMatrix(matrixFilename);
    } else {
      File file = new File(matrixFilename);
      if (file.exists() && file.delete()) {
        logger.info(format(" RMSD file (%s) was deleted.", matrixFilename));
        logger.info(" To restart from a previous run, use the '-r' flag.");
      }
    }

    // Number of atoms used for the superposition
    int nUsed = usedIndices.length;
    int nUsedVars = nUsed * 3;

    // Load atomic coordinates for the base system.
    MolecularAssembly baseMolecularAssembly = baseFilter.getActiveMolecularSystem();
    ForceFieldEnergy baseForceFieldEnergy = baseMolecularAssembly.getPotentialEnergy();
    int nVars = baseForceFieldEnergy.getNumberOfVariables();
    double[] baseCoords = new double[nVars];
    baseForceFieldEnergy.getCoordinates(baseCoords);
    double[] baseUsedCoords = new double[nUsedVars];

    // Load an array for mass weighting.
    Atom[] atoms = baseMolecularAssembly.getAtomArray();
    double[] mass = new double[nUsed];
    for (int i = 0; i < nUsed; i++) {
      mass[i] = atoms[usedIndices[i]].getMass();
    }

    // Allocate space for the superposition
    MolecularAssembly targetMolecularAssembly = targetFilter.getActiveMolecularSystem();
    ForceFieldEnergy targetForceFieldEnergy = targetMolecularAssembly.getPotentialEnergy();
    double[] targetCoords = new double[nVars];
    double[] targetUsedCoords = new double[nUsedVars];

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

    // Loop over base structures.
    for (int row = restartRow; row < baseSize; row++) {

      // Initialize the distance this rank is responsible for to zero.
      myDistance[0] = -1.0;

      if (row == restartRow) {
        if (dRMSD) {
          logger.info(
              "\n Coordinate RMSD\n Snapshots       Original   After Translation   After Rotation     dRMSD");
        } else if (verbose) {
          logger.info(
              "\n Coordinate RMSD\n Snapshots       Original   After Translation   After Rotation");
        }
      }

      // Loop over target structures.
      for (int column = restartColumn; column < paddedTargetSize; column++) {

        // Make sure this is not a padded value of column.
        if (column < targetSize) {
          int targetRank = column % numProc;
          if (targetRank == rank) {
            if (isSymmetric && row == column) {
              // Fill the diagonal.
              myDistance[0] = 0.0;
              if (verbose) {
                logger.info(format(" %6d  %6d  %s                             %8.5f",
                    row + 1, column + 1, "Diagonal", myDistance[0]));
              }
            } else if (isSymmetric && row >= column) {
              // Fill the lower triangle from the upper triangle.
              myDistance[0] = distMatrix[column][row];
              if (verbose) {
                logger.info(format(" %6d  %6d  %s                  %8.5f",
                    row + 1, column + 1, "From Upper Triangle", myDistance[0]));
              }
            } else {
              // Calculate an upper triangle entry.

              // Load base coordinates.
              baseForceFieldEnergy.getCoordinates(baseCoords);

              // Save the current coordinates of the target
              AssemblyState origStateB = new AssemblyState(targetMolecularAssembly);

              // Load the target coordinates.
              targetForceFieldEnergy.getCoordinates(targetCoords);

              copyCoordinates(usedIndices, baseCoords, baseUsedCoords);
              copyCoordinates(usedIndices, targetCoords, targetUsedCoords);

              double origRMSD = rmsd(baseUsedCoords, targetUsedCoords, mass);

              // Calculate the translation on only the used subset, but apply it to the entire structure.
              applyTranslation(baseCoords, calculateTranslation(baseUsedCoords, mass));
              applyTranslation(targetCoords, calculateTranslation(targetUsedCoords, mass));
              // Copy the applied translation to baseUsedCoords and targetUsedCoords.
              copyCoordinates(usedIndices, baseCoords, baseUsedCoords);
              copyCoordinates(usedIndices, targetCoords, targetUsedCoords);
              double translatedRMSD = rmsd(baseUsedCoords, targetUsedCoords, mass);

              // Calculate the rotation on only the used subset, but apply it to the entire structure.
              applyRotation(targetCoords, calculateRotation(baseUsedCoords, targetUsedCoords, mass));
              // Copy the applied rotation to targetUsedCoords.
              copyCoordinates(usedIndices, targetCoords, targetUsedCoords);
              double rotatedRMSD = rmsd(baseUsedCoords, targetUsedCoords, mass);

              if (dRMSD) {
                double disRMSD = calcDRMSD(baseUsedCoords, targetUsedCoords, nUsed * 3);
                logger.info(format(" %6d  %6d  %8.5f            %8.5f         %8.5f  %8.5f",
                    row + 1, column + 1, origRMSD, translatedRMSD, rotatedRMSD, disRMSD));
              } else if (verbose) {
                logger.info(format(" %6d  %6d  %8.5f            %8.5f         %8.5f",
                    row + 1, column + 1, origRMSD, translatedRMSD, rotatedRMSD));
              }

              // Store the RMSD result.
              myDistance[0] = rotatedRMSD;

              // Save a PDB snapshot.
              if (saveSnapshots && numProc == 1) {
                MolecularAssembly molecularAssembly = targetOutputFilter.getActiveMolecularSystem();
                molecularAssembly.getPotentialEnergy().setCoordinates(targetCoords);
                targetOutputFilter.writeFile(targetOutputFile, true);
                origStateB.revertState();
              }
            }
          }

          // Read ahead to the next target structure.
          targetFilter.readNext(false, false);
        }

        // Every numProc iterations, send the results of each rank.
        if ((column + 1) % numProc == 0) {
          gatherRMSDs(row, column);
        }
      }

      // Reset the starting column.
      restartColumn = 0;

      // Reset the targetFilter and read the first structure.
      targetFilter.readNext(true, false);

      // Read the next base structure.
      baseFilter.readNext(false, false);

      // Only Rank 0 writes the distance matrix.
      if (rank == 0 && write) {
        writeDistanceMatrixRow(matrixFilename, distMatrix[row]);
      }
    }

    baseFilter.closeReader();
    targetFilter.closeReader();
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
   * Copy coordinates from the entire system to the used subset.
   *
   * @param usedIndices Mapping from the xUsed array to its source in x.
   * @param x All atomic coordinates.
   * @param xUsed The used subset of coordinates.
   */
  public static void copyCoordinates(int[] usedIndices, double[] x, double[] xUsed) {
    int nUsed = usedIndices.length;
    for (int u = 0; u < nUsed; u++) {
      int u3 = 3 * u;
      int i3 = 3 * usedIndices[u];
      arraycopy(x, i3, xUsed, u3, 3);
    }
  }

  /**
   * Calculates the dRMSD between to sets of coordinates.
   *
   * @param xUsed A double array containing the xyz coordinates for multiple atoms.
   * @param x2Used A double array containing the xyz coordinates for multiple atoms.
   * @param nUsed The number of atoms that dRMSD is calculated on.
   * @return A double containing the dRMSD value.
   */
  public static double calcDRMSD(double[] xUsed, double[] x2Used, int nUsed) {
    double disRMSD = 0.0;
    int counter = 0;
    for (int i = 0; i < nUsed; i = i + 3) {
      Double3 xi = new Double3(xUsed[i], xUsed[i + 1], xUsed[i + 2]);
      Double3 x2i = new Double3(x2Used[i], x2Used[i + 1], x2Used[i + 2]);
      for (int j = i + 3; j < nUsed; j = j + 3) {
        Double3 xj = new Double3(xUsed[j], xUsed[j + 1], xUsed[j + 2]);
        Double3 x2j = new Double3(x2Used[j], x2Used[j + 1], x2Used[j + 2]);
        double dis1 = xi.sub(xj).length();
        double dis2 = x2i.sub(x2j).length();
        double diff = dis1 - dis2;
        disRMSD += diff * diff;
        counter++;
      }
    }
    disRMSD = disRMSD / counter;
    return sqrt(disRMSD);
  }

  /**
   * Minimize the RMS distance between two sets of atoms using quaternions and a pre-calculated
   * rotation matrix; overlaps x2 onto x1.
   *
   * @param x2 Cartesian coordinates of the second system. Modified in-place.
   * @param rot A pre-calculated rotation matrix.
   */
  public static void applyRotation(double[] x2, double[][] rot) {
    int n = x2.length / 3;

    // Rotate second molecule to best fit with first molecule
    for (int i = 0; i < n; i++) {
      int k = i * 3;
      double xrot = x2[k] * rot[0][0] + x2[k + 1] * rot[0][1] + x2[k + 2] * rot[0][2];
      double yrot = x2[k] * rot[1][0] + x2[k + 1] * rot[1][1] + x2[k + 2] * rot[1][2];
      double zrot = x2[k] * rot[2][0] + x2[k + 1] * rot[2][1] + x2[k + 2] * rot[2][2];
      x2[k] = xrot;
      x2[k + 1] = yrot;
      x2[k + 2] = zrot;
    }
  }

  /**
   * Apply a translation matrix [dx,dy,dz] to a molecular system.
   *
   * @param x Coordinates to move; modified in-place.
   * @param translation Translation matrix.
   */
  public static void applyTranslation(double[] x, final double[] translation) {
    int n = x.length / 3;
    for (int i = 0; i < n; i++) {
      int i3 = 3 * i;
      for (int j = 0; j < 3; j++) {
        x[i3 + j] -= translation[j];
      }
    }
  }

  /**
   * Calculate a rotation to minimize RMS distance between two sets of atoms using quaternions,
   * overlapping x2 on x1.
   *
   * @param x1 Cartesian coordinates of the first system.
   * @param x2 Cartesian coordinates of the second system.
   * @param mass The mass of each particle in the system.
   * @return A rotation matrix.
   */
  public static double[][] calculateRotation(double[] x1, double[] x2, double[] mass) {

    // Build the upper triangle of the quadratic form matrix
    double xxyx = 0.0;
    double xxyy = 0.0;
    double xxyz = 0.0;
    double xyyx = 0.0;
    double xyyy = 0.0;
    double xyyz = 0.0;
    double xzyx = 0.0;
    double xzyy = 0.0;
    double xzyz = 0.0;

    int n = x1.length / 3;
    for (int i = 0; i < n; i++) {
      int k = i * 3;
      double weigh = mass[i];
      xxyx = xxyx + weigh * x1[k] * x2[k];
      xxyy = xxyy + weigh * x1[k + 1] * x2[k];
      xxyz = xxyz + weigh * x1[k + 2] * x2[k];
      xyyx = xyyx + weigh * x1[k] * x2[k + 1];
      xyyy = xyyy + weigh * x1[k + 1] * x2[k + 1];
      xyyz = xyyz + weigh * x1[k + 2] * x2[k + 1];
      xzyx = xzyx + weigh * x1[k] * x2[k + 2];
      xzyy = xzyy + weigh * x1[k + 1] * x2[k + 2];
      xzyz = xzyz + weigh * x1[k + 2] * x2[k + 2];
    }

    double[][] c = new double[4][4];
    c[0][0] = xxyx + xyyy + xzyz;
    c[0][1] = xzyy - xyyz;
    c[1][0] = c[0][1];
    c[1][1] = xxyx - xyyy - xzyz;
    c[0][2] = xxyz - xzyx;
    c[2][0] = c[0][2];
    c[1][2] = xxyy + xyyx;
    c[2][1] = c[1][2];
    c[2][2] = xyyy - xzyz - xxyx;
    c[0][3] = xyyx - xxyy;
    c[3][0] = c[0][3];
    c[1][3] = xzyx + xxyz;
    c[3][1] = c[1][3];
    c[2][3] = xyyz + xzyy;
    c[3][2] = c[2][3];
    c[3][3] = xzyz - xxyx - xyyy;

    // Diagonalize the quadratic form matrix
    Array2DRowRealMatrix cMatrix = new Array2DRowRealMatrix(c, false);
    EigenDecomposition eigenDecomposition = new EigenDecomposition(cMatrix);

    // Extract the quaternions.
    double[] q = eigenDecomposition.getEigenvector(0).toArray();

    // Assemble rotation matrix that superimposes the molecules
    double[][] rot = new double[3][3];
    double q02 = q[0] * q[0];
    double q12 = q[1] * q[1];
    double q22 = q[2] * q[2];
    double q32 = q[3] * q[3];
    rot[0][0] = q02 + q12 - q22 - q32;
    rot[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    rot[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
    rot[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
    rot[1][1] = q02 - q12 + q22 - q32;
    rot[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
    rot[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
    rot[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    rot[2][2] = q02 - q12 - q22 + q32;

    return rot;
  }

  /**
   * Calculate a translation matrix [dx,dy,dz] to return a molecular system's center of mass to the
   * origin.
   *
   * @param x Coordinates of the system.
   * @param mass Mass of each atom.
   * @return Translation to the origin.
   */
  public static double[] calculateTranslation(double[] x, final double[] mass) {
    double xmid = 0.0;
    double ymid = 0.0;
    double zmid = 0.0;
    double norm = 0.0;
    int n = x.length / 3;

    for (int i = 0; i < n; i++) {
      int k = 3 * i;
      double weigh = mass[i];
      xmid = xmid + x[k] * weigh;
      ymid = ymid + x[k + 1] * weigh;
      zmid = zmid + x[k + 2] * weigh;
      norm = norm + weigh;
    }
    xmid = xmid / norm;
    ymid = ymid / norm;
    zmid = zmid / norm;

    return new double[] {xmid, ymid, zmid};
  }

  /**
   * Compute the rms fit over superimposed atom pairs
   *
   * @param x1 Cartesian coordinates of the first system.
   * @param x2 Cartesian coordinates of the second system.
   * @param mass The mass of each particle in the system.
   * @return The RMSD.
   */
  public static double rmsd(double[] x1, double[] x2, double[] mass) {
    double rmsfit = 0.0;
    double norm = 0.0;
    int n = x1.length / 3;
    for (int i = 0; i < n; i++) {
      int k = 3 * i;
      double weigh = mass[i];
      double xr = x1[k] - x2[k];
      double yr = x1[k + 1] - x2[k + 1];
      double zr = x1[k + 2] - x2[k + 2];
      double dist2 = xr * xr + yr * yr + zr * zr;
      norm = norm + weigh;
      double rmsterm = dist2 * weigh;
      rmsfit = rmsfit + rmsterm;
    }
    return sqrt(rmsfit / norm);
  }

  /**
   * Minimize the RMS distance between two sets of atoms using quaternions.
   *
   * @param x1 Cartesian coordinates of the first system. Unmodified.
   * @param x2 Cartesian coordinates of the second system. Modified in-place.
   * @param mass The mass of each particle in the system.
   */
  public static void rotate(double[] x1, double[] x2, double[] mass) {
    double[][] rotation = calculateRotation(x1, x2, mass);
    applyRotation(x2, rotation);
  }

  /**
   * Move the center of mass for both sets of atoms to the origin.
   *
   * @param x1 Cartesian coordinates of the first system.
   * @param mass1 The mass of each particle in the first system.
   * @param x2 Cartesian coordinates of the second system.
   * @param mass2 The mass of each particle in the second system.
   */
  public static void translate(double[] x1, double[] mass1, double[] x2, double[] mass2) {
    // Move the first set of atoms.
    translate(x1, mass1);
    // Move the second set of atoms.
    translate(x2, mass2);
  }

  /**
   * Move the center of mass for a set of atoms to the origin.
   *
   * @param x Cartesian coordinates of the system; modified in-place.
   * @param mass The mass of each particle in the system.
   */
  private static void translate(double[] x, final double[] mass) {
    double[] translation = calculateTranslation(x, mass);
    applyTranslation(x, translation);
  }
}
