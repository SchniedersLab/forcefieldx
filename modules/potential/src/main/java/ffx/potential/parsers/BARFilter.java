//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Logger;

import static ffx.potential.parsers.SystemFilter.version;
import static java.lang.Double.isNaN;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

/**
 * The BARFilter class parses TINKER bar(*.BAR) files.
 *
 * @author Rose A. Gogal
 * @since 1.0
 */

public class BARFilter {

  private static final Logger logger = Logger.getLogger(XYZFilter.class.getName());
  private final File barFile;

  private int snaps1;
  private double temperature1;
  private double[] e1l1;
  private double[] e1l2;
  private double[] volume1;

  private int snaps2;
  private double temperature2;
  private double[] e2l1;
  private double[] e2l2;
  private double[] volume2;

  private int startingSnap = 0;
  private int endingSnap = 0;
  private int count = 0;


  /**
   * BARFilter constructor
   *
   * @param barFile a {@link java.util.List} object.
   */
  public BARFilter(File barFile) {
    this.barFile = barFile;
  }

  /**
   * BARFilter constructor
   *
   * @param barFile      a {@link java.util.List} object.
   * @param startingSnap a {@link java.util.List} object.
   * @param endingSnap   a {@link java.util.List} object.
   */
  public BARFilter(File barFile, int startingSnap, int endingSnap) {
    this.barFile = barFile;
    this.startingSnap = startingSnap;
    this.endingSnap = endingSnap;
  }

  /**
   * BARFilter constructor
   *
   * @param xyzFile      a {@link java.util.List} object.
   * @param e1l1         energy in ensemble 1 at lambda 1
   * @param e1l2         energy in ensemble 1 at lambda 2
   * @param e2l1         energy in ensemble 2 at lambda 1
   * @param e2l2         energy in ensemble 2 at lambda 2
   * @param volume1      volume in ensemble 1
   * @param volume2      volume in ensemble 2
   * @param temperature1 temperature of ensemble 1
   * @param temperature2 temperature of ensemble 2
   */
  public BARFilter(File xyzFile, double[] e1l1, double[] e1l2, double[] e2l1, double[] e2l2,
                   double[] volume1, double[] volume2, double temperature1, double temperature2) {
    this.barFile = xyzFile;

    this.temperature1 = temperature1;
    this.e1l1 = e1l1;
    this.e1l2 = e1l2;
    this.volume1 = volume1;

    this.e2l1 = e2l1;
    this.e2l2 = e2l2;
    this.volume2 = volume2;
    this.temperature2 = temperature2;
  }


  /**
   * Read TINKER bar files and parse the snapshots into energy arrays
   *
   * @return True if the file was read successfully.
   */
  public boolean readFile() {
    ArrayList<Double> ens1lam1 = new ArrayList<>();
    ArrayList<Double> ens1lam2 = new ArrayList<>();
    ArrayList<Double> ens2lam1 = new ArrayList<>();
    ArrayList<Double> ens2lam2 = new ArrayList<>();
    ArrayList<Double> vol1 = new ArrayList<>();
    ArrayList<Double> vol2 = new ArrayList<>();
    // Processes all snapshots in a file.
    int snapshots = 0;
    int xyzCount = 0;
    try (BufferedReader br = new BufferedReader(new FileReader(barFile))) {
      String data;
      while ((data = br.readLine()) != null) {
        String[] tokens = data.trim().split(" +");
        int numTokens = tokens.length;
        if (data.contains(".xyz") || data.contains(".pdb") || numTokens < 3) {
          xyzCount++;
          if (xyzCount == 1) {
            snaps1 = parseInt(tokens[0]);
            temperature1 = parseDouble(tokens[1]);
          } else if (xyzCount == 2) {
            snaps2 = parseInt(tokens[0]);
            temperature2 = parseDouble(tokens[1]);
          }
        } else if (endingSnap != 0) {
          count++;
          snapshots = (endingSnap - startingSnap) + 1;
          if (count >= startingSnap + 1 && count <= endingSnap + 1) {
            if (count <= snaps1) {
              if (numTokens == 4) {
                vol1.add(parseDouble(tokens[3]));
              }
              ens1lam1.add(parseDouble(tokens[1]));
              ens1lam2.add(parseDouble(tokens[2]));
            } else {
              logger.warning(format(" BAR entry of (%3d) is larger than total entries (%3d).", count, snaps1));
            }
          } else if (count >= snaps1 + startingSnap + 1 && count <= snaps1 + endingSnap + 1) {
            if (count <= snaps1 + snaps2 + 1) {
              if (numTokens == 4) {
                vol2.add(parseDouble(tokens[3]));
              }
              ens2lam1.add(parseDouble(tokens[1]));
              ens2lam2.add(parseDouble(tokens[2]));
            } else {
              logger.warning(format(" BAR entry of (%3d) is larger than total entries (%3d).", count, snaps1 + snaps2));
            }
          }
        } else {
          count++;
          if (count <= snaps1) {
            if (numTokens == 4) {
              vol1.add(parseDouble(tokens[3]));
            }
            ens1lam1.add(parseDouble(tokens[1]));
            ens1lam2.add(parseDouble(tokens[2]));
          } else {
            if (numTokens == 4) {
              vol2.add(parseDouble(tokens[3]));
            }
            ens2lam1.add(parseDouble(tokens[1]));
            ens2lam2.add(parseDouble(tokens[2]));
          }
        }
      }
      if (snapshots != 0) {
        e1l1 = new double[snapshots];
        e1l2 = new double[snapshots];
        e2l1 = new double[snapshots];
        e2l2 = new double[snapshots];
        volume1 = new double[snapshots];
        volume2 = new double[snapshots];
        snaps1 = snapshots;
      } else {
        e1l1 = new double[snaps1];
        e1l2 = new double[snaps1];
        e2l1 = new double[snaps2];
        e2l2 = new double[snaps2];
        volume1 = new double[snaps1];
        volume2 = new double[snaps2];
      }
      for (int i = 0; i < ens1lam1.size(); i++) {
        e1l1[i] = ens1lam1.get(i);
        e1l2[i] = ens1lam2.get(i);
        if (!vol1.isEmpty()) {
          volume1[i] = vol1.get(i);
        }
      }
      for (int i = 0; i < ens2lam1.size(); i++) {
        e2l1[i] = ens2lam1.get(i);
        e2l2[i] = ens2lam2.get(i);
        if (!vol1.isEmpty()) {
          volume2[i] = vol2.get(i);
        }
      }
      // Read blank lines at the top of the file
      if (data == null) {
        return false;
      }
    } catch (IOException fileNotFoundException) {
      logger.warning(format(" Exception reading %s:\n %s", barFile, fileNotFoundException));
    }
    return true;
  }

  /**
   * Write TINKER bar files
   *
   * @param saveFile The file to write to.
   * @param isPBC    include volume in the output file.
   * @return True if successful.
   */
  public boolean writeFile(String saveFile, boolean isPBC) {
    return writeFile(saveFile, isPBC, true);
  }

  /**
   * Write TINKER bar files
   *
   * @param saveFile The file to write to.
   * @param isPBC    include volume in the output file.
   * @param append   If the append flag is true, "saveFile" will be appended to. Otherwise, the default versioning scheme will be applied.
   * @return True if successful.
   */
  public boolean writeFile(String saveFile, boolean isPBC, boolean append) {
    int snaps = e1l1.length;
    int snaps2 = e2l1.length;
    String name = barFile.getName();

    File newFile = new File(saveFile);
    if (!append) {
      newFile = version(newFile);
    }

    logger.info(format(" Writing BAR file: %s", newFile));
    try (FileWriter fw = new FileWriter(newFile,
        append && newFile.exists()); BufferedWriter bw = new BufferedWriter(fw)) {
      bw.write(format("%8d %9.3f %s\n", snaps, temperature1, name));
      for (int i = 0; i < snaps; i++) {
        if (isNaN(e1l1[i]) || isNaN(e1l2[i])) {
          continue;
        }
        if (isPBC) {
          bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e1l1[i], e1l2[i], volume2[i]));
        } else {
          bw.write(format("%8d %20.10f %20.10f\n", i + 1, e1l1[i], e1l2[i]));
        }
      }

      bw.write(format("%8d %9.3f  %s\n", snaps2, temperature1, name));
      for (int i = 0; i < snaps2; i++) {
        if (isNaN(e2l1[i]) || isNaN(e2l2[i])) {
          continue;
        }
        if (isPBC) {
          bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e2l1[i], e2l2[i], volume2[i]));
        } else {
          bw.write(format("%8d %20.10f %20.10f\n", i + 1, e2l1[i], e2l2[i]));
        }
      }
    } catch (IOException e) {
      logger.warning(format(" Exception writing %s", newFile));
      return false;
    }
    return true;
  }

  /**
   * Returns the temperature for ensemble 1.
   *
   * @return The temperature.
   */
  public double getTemperature1() {
    return temperature1;
  }

  /**
   * Returns the temperature for ensemble 2.
   *
   * @return The temperature.
   */
  public double getTemperature2() {
    return temperature2;
  }

  /**
   * Return the potential energy for each snapshot of ensemble 1 at lambda 1.
   *
   * @return The potential energy for each snapshot of ensemble 1 at lambda 1.
   */
  public double[] getE1l1() {
    return e1l1;
  }

  /**
   * Return the potential energy for each snapshot of ensemble 1 at lambda 2.
   *
   * @return The potential energy for each snapshot of ensemble 1 at lambda 2.
   */
  public double[] getE1l2() {
    return e1l2;
  }

  /**
   * Return the potential energy for each snapshot of ensemble 2 at lambda 1.
   *
   * @return The potential energy for each snapshot of ensemble 2 at lambda 1.
   */
  public double[] getE2l1() {
    return e2l1;
  }

  /**
   * Return the potential energy for each snapshot of ensemble 2 at lambda 2.
   *
   * @return The potential energy for each snapshot of ensemble 2 at lambda 2.
   */
  public double[] getE2l2() {
    return e2l2;
  }

  /**
   * Return the volume for each snapshot of ensemble 1.
   *
   * @return The volume for each snapshot of ensemble 1.
   */
  public double[] getVolume1() {
    return volume1;
  }

  /**
   * Return the volume for each snapshot of ensemble 2.
   *
   * @return The volume for each snapshot of ensemble 2.
   */
  public double[] getVolume2() {
    return volume2;
  }

  /**
   * Return the number of snapshots in the BAR file.
   *
   * @return The number of snapshots in the BAR file.
   */
  public int getNumberOfSnapshots() {
    return snaps1;
  }

  /**
   * Return the number of snapshots in the BAR file.
   *
   * @return The number of snapshots in the BAR file.
   */
  public File getFile() {
    return barFile;
  }
}
