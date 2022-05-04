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

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Double.parseDouble;
import static java.lang.String.format;

/**
 * The ESVFilter class parses Extended System Restart (*.ESV) files.
 *
 * @author Andrew Thiel
 * @since 1.0
 */
public class ESVFilter {

  private static final Logger logger = Logger.getLogger(ESVFilter.class.getName());
  private final String label;

  /**
   * Constructor for ESVFilter.
   *
   * @param label a Label for the this restart file.
   */
  public ESVFilter(String label) {
    this.label = label;
  }

  /**
   * readDYN
   *
   * @param esvFile a {@link File} object.
   * @param x an array of double.
   * @param v an array of double.
   * @param a an array of double.
   * @return a boolean.
   */
  public boolean readESV(
      File esvFile, double[] x, double[] v, double[] a) {
    if (!esvFile.exists() || !esvFile.canRead()) {
      return false;
    }
    try (BufferedReader br = new BufferedReader(new FileReader(esvFile))) {

      br.readLine();
      String data = br.readLine().trim();
      String[] tokens = data.split(" +");
      if (tokens.length == 0) {
        return false;
      }
      int numESVs = Integer.parseInt(tokens[0]);

      // Atomic coordinates
      br.readLine();
      for (int i = 0; i < numESVs; i++) {
        data = br.readLine().trim();
        tokens = data.split(" +");
        if (tokens.length != 1) {
          return false;
        }
        x[i] = parseDouble(tokens[0]);
      }

      // Velocities
      br.readLine();
      for (int i = 0; i < numESVs; i++) {
        data = br.readLine().trim();
        tokens = data.split(" +");
        if (tokens.length != 1) {
          return false;
        }
        v[i] = parseDouble(tokens[0]);
      }

      // Accelerations
      br.readLine();
      for (int i = 0; i < numESVs; i++) {
        data = br.readLine().trim();
        tokens = data.split(" +");
        if (tokens.length != 1) {
          return false;
        }
        a[i] = parseDouble(tokens[0]);
      }
    } catch (Exception e) {
      String message = "Exception reading ESV restart file: " + esvFile;
      logger.log(Level.WARNING, message, e);
    }
    return true;
  }

  /**
   * writeDYN
   *
   * @param dynFile a {@link File} object.
   * @param x an array of double.
   * @param v an array of double.
   * @param a an array of double.
   * @return a boolean.
   */
  public boolean writeESV(
      File dynFile, double[] x, double[] v, double[] a) {
    FileWriter fw = null;
    BufferedWriter bw = null;
    try {
      fw = new FileWriter(dynFile);
      bw = new BufferedWriter(fw);

      bw.write(" Number of ESVs and Title :\n");
      int numberOfAtoms = x.length;
      String output = format("%7d  %s\n", numberOfAtoms, label);
      bw.write(output);

      bw.write(" Current Theta Positions :\n");
      for (int i = 0; i < numberOfAtoms; i++) {
        bw.write(format("%26.16E\n", x[i]));
      }

      bw.write(" Current Atomic Velocities :\n");
      for (int i = 0; i < numberOfAtoms; i++) {
        bw.write(format("%26.16E\n", v[i]));
      }

      bw.write(" Current Atomic Accelerations :\n");
      for (int i = 0; i < numberOfAtoms; i++) {
        bw.write(format("%26.16E\n", a[i]));
      }
    } catch (IOException e) {
      String message = " Exception writing dynamic restart file " + dynFile;
      logger.log(Level.SEVERE, message, e);
      return false;
    } finally {
      try {
        bw.close();
        fw.close();
        return true;
      } catch (IOException e) {
        String message = " Exception closing dynamic restart file " + dynFile;
        logger.log(Level.WARNING, message, e);
        return false;
      }
    }
  }
}
