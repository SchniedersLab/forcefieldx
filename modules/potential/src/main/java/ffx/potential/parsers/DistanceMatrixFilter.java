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
package ffx.potential.parsers;

import static java.lang.Double.parseDouble;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import ffx.numerics.math.RunningStatistics;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * The DistanceMatrixFilter class parses a Distance Matrix (*.TXT) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DistanceMatrixFilter {

  private static final Logger logger = Logger.getLogger(DistanceMatrixFilter.class.getName());

  /**
   * No public constructor for DistanceMatrixFilter.
   */
  public DistanceMatrixFilter() {
  }

  private int nRows = 0;
  private int nColumns = 0;
  private boolean fillDistanceMatrix = false;
  private double[][] distanceMatrix = null;

  /**
   * Get the number of rows read in.
   *
   * @return The number of rows read in.
   */
  public int getRestartRow() {
    return nRows;
  }

  /**
   * Get the number of columns in the last row that was read in.
   *
   * @return The number of columns in the last row that was read in.
   */
  public int getRestartColumn() {
    return nColumns;
  }

  /**
   * Read in the distance matrix from a file. This method is used from the Clustering code, where the
   * size of the distance matrix is unknown, and the full N^2 matrix is needed.
   *
   * @param filename The filename to read from.
   * @param distanceList The distance matrix (if null, only the last row is returned).
   * @return Statistics for the parsed values.
   */
  public RunningStatistics readDistanceMatrix(String filename, List<double[]> distanceList) {
    fillDistanceMatrix = true;
    RunningStatistics runningStatistics = readDistanceMatrix(filename, -1, -1);
    if (distanceMatrix != null) {
      int size = distanceMatrix[0].length;
      boolean square = true;
      for (int i = 0; i < size; i++) {
        if (distanceMatrix[i] == null || distanceMatrix[i].length != size) {
          square = false;
          break;
        }
      }
      // The full NxN matrix has been read in.
      if (square) {
        // Add all rows of the distanceMatrix into the list
        Collections.addAll(distanceList, distanceMatrix);
      } else {
        // Assume only the upper triangle has been read in.
        for (int i = 0; i < size; i++) {
          double[] row = new double[size];
          // Fill the lower triangle from previous rows.
          for (int j = 0; j < i; j++) {
            double[] previousRow = distanceList.get(j);
            row[j] = previousRow[i];
          }
          // Fill the upper triangle using the current row of the distanceMatrix.
          arraycopy(distanceMatrix[i], 0, row, i, size - i);
          distanceList.add(row);
        }
      }
    }

    // Reset the DistanceMatrixFilter for reuse.
    distanceMatrix = null;
    fillDistanceMatrix = false;
    return runningStatistics;
  }

  /**
   * Read in the distance matrix from a file.
   *
   * @param filename The filename to read from.
   * @param expectedRows The number of rows to expect.
   * @param expectedColumns The number of columns to expect.
   * @return Statistics for the parsed values.
   */
  public RunningStatistics readDistanceMatrix(
      String filename, int expectedRows, int expectedColumns) {

    if (filename == null) {
      return null;
    }

    File distanceMatrixFile = new File(filename);
    if (!distanceMatrixFile.exists() || !distanceMatrixFile.canRead()) {
      return null;
    }

    // Read in the RMSD matrix.
    try (FileReader fr = new FileReader(distanceMatrixFile);
        BufferedReader br = new BufferedReader(fr)) {

      String data = br.readLine();

      // Check for blank lines at the top of the file
      while (data != null && data.trim().equals("")) {
        data = br.readLine();
      }

      if (data == null) {
        logger.info(format("\n No data in RMSD file %s.", distanceMatrixFile));
        return null;
      }

      String[] tokens = data.trim().split(" +");
      nColumns = tokens.length;

      // If the expectedRows is unknown (i.e. this is called from the Cluster script), set
      // the number of rows to the initial number of parsed tokens.
      if (expectedRows == -1) {
        expectedRows = nColumns;
      }
      if (expectedColumns == -1) {
        expectedColumns = nColumns;
      }

      if (nColumns != expectedColumns) {
        logger.info(
            format("\n Unexpected number of entries (%d) in the first row of the RMSD file %s.",
                nColumns, distanceMatrixFile));
        return null;
      }

      if (fillDistanceMatrix) {
        distanceMatrix = new double[expectedRows][];
      }

      RunningStatistics runningStatistics = new RunningStatistics();

      // Loop over the number of expected rows.
      for (int i = 0; i < expectedRows; i++) {
        double[] row = new double[nColumns];
        for (int j = 0; j < nColumns; j++) {
          row[j] = parseDouble(tokens[j]);
          runningStatistics.addValue(row[j]);
        }

        // Are we storing all rows?
        if (distanceMatrix != null && distanceMatrix.length > i) {
          distanceMatrix[i] = row;
        }

        nRows = i + 1;

        // Read the next line.
        data = br.readLine();
        if (data != null) {
          tokens = data.trim().split(" +");
        } else {
          break;
        }

        nColumns = tokens.length;
      }
      return runningStatistics;

    } catch (IOException e) {
      logger.info(format(" Exception reading %s:\n %s", distanceMatrixFile, e));
      return null;
    }
  }

  /**
   * Convert a distance matrix to a String.
   *
   * @param distanceMatrix The distance matrix.
   * @return Return a String representation (or null if the distanceMatrix is null).
   */
  public static String toDistanceMatrixString(List<double[]> distanceMatrix) {

    if (distanceMatrix == null) {
      return null;
    }

    StringBuilder sb = new StringBuilder("\n Distance Matrix:\n");
    for (double[] row : distanceMatrix) {
      sb.append("  ");
      for (int j = 0; j < row.length; j++) {
        if (row[j] == -2.0) {
          sb.append(format("%6.4f", Double.NaN));
        } else {
          sb.append(format("%6.4f", row[j]));
        }
        if (j == row.length - 1) {
          sb.append("\n");
        } else {
          sb.append(" ");
        }
      }
    }

    return sb.toString();
  }

  /**
   * Convert a distance matrix to a String.
   *
   * @param distanceMatrix The distance matrix.
   * @return Return a String representation (or null if the distanceMatrix is null).
   */
  public static String toDistanceMatrixString(double[][] distanceMatrix) {
    return toDistanceMatrixString(Arrays.asList(distanceMatrix));
  }

  /**
   * Write the distance matrix to a file.
   *
   * @param filename The filename to write to.
   * @param distanceMatrix a list containing rows of the distance matrix.
   * @return a boolean.
   */
  public static boolean writeDistanceMatrix(String filename, List<double[]> distanceMatrix) {
    // Check for null arguments.
    if (distanceMatrix == null) {
      return false;
    }

    for (double[] row : distanceMatrix) {
      boolean success = writeDistanceMatrixRow(filename, row, 0);
      if (!success) {
        return false;
      }
    }

    return true;
  }

  /**
   * Write the distance matrix to a file.
   *
   * @param filename The filename to write to.
   * @param distanceMatrixRow A row of the distance matrix.
   * @param firstColumn The first column of the distance matrix row to write.
   * @return a boolean.
   */
  public static boolean writeDistanceMatrixRow(String filename, double[] distanceMatrixRow,
      int firstColumn) {

    if (filename == null) {
      return false;
    }

    File distanceMatrixFile = new File(filename);

    // Check for null arguments.
    if (distanceMatrixFile == null || distanceMatrixRow == null) {
      return false;
    }

    // Write one row.
    try (FileWriter fw = new FileWriter(distanceMatrixFile, true);
        BufferedWriter bw = new BufferedWriter(fw)) {
      int nColumn = distanceMatrixRow.length;
      for (int column = firstColumn; column < nColumn; column++) {
        bw.write(format("%16.14f", distanceMatrixRow[column]));
        if (column < nColumn - 1) {
          bw.write(" ");
        }
      }
      bw.write("\n");
    } catch (Exception e) {
      logger.info(format(" Exception writing to %s:\n %s", distanceMatrixFile, e));
      return false;
    }

    return true;
  }
}
