/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.potential.parsers;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import ffx.crystal.Crystal;

/**
 * The DYNFilter class parses TINKER Restart (*.DYN) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class DYNFilter {

    private static final Logger logger = Logger.getLogger(DYNFilter.class.getName());
    private String label;

    /**
     * <p>Constructor for DYNFilter.</p>
     *
     */
    public DYNFilter(String label) {
        this.label = label;
    }

    /**
     * <p>readDYN</p>
     *
     * @param dynFile a {@link java.io.File} object.
     * @param x an array of double.
     * @param v an array of double.
     * @param a an array of double.
     * @param ap an array of double.
     * @return a boolean.
     */
    public boolean readDYN(File dynFile, double x[], double v[], double a[], double ap[]) {
        if (!dynFile.exists() || !dynFile.canRead()) {
            return false;
        }
        FileReader fr = null;
        BufferedReader br = null;
        try {
            fr = new FileReader(dynFile);
            br = new BufferedReader(fr);
            br.readLine();
            String data = br.readLine().trim();
            String tokens[] = data.split(" +");
            if (tokens.length == 0) {
                return false;
            }
            int numatoms = Integer.parseInt(tokens[0]);

            // Box size and angles
            br.readLine();
            data = br.readLine().trim();
            tokens = data.split(" +");
            if (tokens.length != 3) {
                return false;
            }
            double d[] = new double[3];
            d[0] = Double.parseDouble(tokens[0]);
            d[1] = Double.parseDouble(tokens[1]);
            d[2] = Double.parseDouble(tokens[2]);

            data = br.readLine().trim();
            tokens = data.split(" +");
            if (tokens.length != 3) {
                return false;
            }
            d[0] = Double.parseDouble(tokens[0]);
            d[1] = Double.parseDouble(tokens[1]);
            d[2] = Double.parseDouble(tokens[2]);

            // Atomic coordinates
            br.readLine();
            for (int i = 0; i < numatoms; i++) {
                data = br.readLine().trim();
                tokens = data.split(" +");
                if (tokens.length != 3) {
                    return false;
                }
                int j = i * 3;
                x[j] = Double.parseDouble(tokens[0]);
                x[j + 1] = Double.parseDouble(tokens[1]);
                x[j + 2] = Double.parseDouble(tokens[2]);
            }

            // Velocities
            br.readLine();
            for (int i = 0; i < numatoms; i++) {
                data = br.readLine().trim();
                tokens = data.split(" +");
                if (tokens.length != 3) {
                    return false;
                }
                int j = i * 3;
                v[j] = Double.parseDouble(tokens[0]);
                v[j + 1] = Double.parseDouble(tokens[1]);
                v[j + 2] = Double.parseDouble(tokens[2]);
            }

            // Accelerations
            br.readLine();
            for (int i = 0; i < numatoms; i++) {
                data = br.readLine().trim();
                tokens = data.split(" +");
                if (tokens.length != 3) {
                    return false;
                }
                int j = i * 3;
                a[j] = Double.parseDouble(tokens[0]);
                a[j + 1] = Double.parseDouble(tokens[1]);
                a[j + 2] = Double.parseDouble(tokens[2]);
            }

            // Previous Accelerations
            br.readLine();
            for (int i = 0; i < numatoms; i++) {
                data = br.readLine().trim();
                tokens = data.split(" +");
                if (tokens.length != 3) {
                    return false;
                }
                int j = i * 3;
                ap[j] = Double.parseDouble(tokens[0]);
                ap[j + 1] = Double.parseDouble(tokens[1]);
                ap[j + 2] = Double.parseDouble(tokens[2]);
            }
        } catch (Exception e) {
            String message = "Exception reading dynamic restart file: " + dynFile;
            logger.log(Level.WARNING, message, e);
        } finally {
            try {
                br.close();
                fr.close();
            } catch (Exception e) {
                String message = "Exception closing restart file " + dynFile;
                logger.log(Level.WARNING, message, e);
                return false;
            }
        }
        return true;
    }

    /**
     * <p>writeDYN</p>
     *
     * @param dynFile a {@link java.io.File} object.
     * @param unitCell a {@link ffx.crystal.Crystal} object.
     * @param x an array of double.
     * @param v an array of double.
     * @param a an array of double.
     * @param ap an array of double.
     * @return a boolean.
     */
    public boolean writeDYN(File dynFile, Crystal unitCell, double x[], double v[],
            double[] a, double ap[]) {
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {
            fw = new FileWriter(dynFile);
            bw = new BufferedWriter(fw);

            bw.write(" Number of Atoms and Title :\n");
            assert (x.length % 3 == 0);
            int numberOfAtoms = x.length / 3;
            String output = format("%7d  %s\n", numberOfAtoms, label);
            bw.write(output);

            bw.write(" Periodic Box Dimensions :\n");
            bw.write(format("%26.16E%26.16E%26.16E\n", unitCell.a, unitCell.b, unitCell.c));
            bw.write(format("%26.16E%26.16E%26.16E\n", unitCell.alpha, unitCell.beta, unitCell.gamma));

            bw.write(" Current Atomic Positions :\n");
            for (int i = 0; i < numberOfAtoms; i++) {
                int k = i * 3;
                bw.write(format("%26.16E%26.16E%26.16E\n", x[k], x[k + 1], x[k + 2]));
            }

            bw.write(" Current Atomic Velocities :\n");
            for (int i = 0; i < numberOfAtoms; i++) {
                int k = i * 3;
                bw.write(format("%26.16E%26.16E%26.16E\n", v[k], v[k + 1], v[k + 2]));
            }

            bw.write(" Current Atomic Accelerations :\n");
            for (int i = 0; i < numberOfAtoms; i++) {
                int k = i * 3;
                bw.write(format("%26.16E%26.16E%26.16E\n", a[k], a[k + 1], a[k + 2]));
            }

            bw.write(" Previous Atomic Accelerations :\n");
            for (int i = 0; i < numberOfAtoms; i++) {
                int k = i * 3;
                bw.write(format("%26.16E%26.16E%26.16E\n", ap[k], ap[k + 1], ap[k + 2]));
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
