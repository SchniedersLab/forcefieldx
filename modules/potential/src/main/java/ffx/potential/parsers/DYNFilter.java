//******************************************************************************
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
//******************************************************************************
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import ffx.crystal.Crystal;

/**
 * The DYNFilter class parses TINKER Restart (*.DYN) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DYNFilter {

    private static final Logger logger = Logger.getLogger(DYNFilter.class.getName());
    private final String label;

    /**
     * <p>
     * Constructor for DYNFilter.</p>
     *
     * @param label a Label for the this restart file.
     */
    public DYNFilter(String label) {
        this.label = label;
    }

    /**
     * <p>
     * readDYN</p>
     *
     * @param dynFile a {@link java.io.File} object.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param x       an array of double.
     * @param v       an array of double.
     * @param a       an array of double.
     * @param ap      an array of double.
     * @return a boolean.
     */
    public boolean readDYN(File dynFile, Crystal crystal, double[] x, double[] v, double[] a, double[] ap) {
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
            String[] tokens = data.split(" +");
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
            double aaxis = Double.parseDouble(tokens[0]);
            double baxis = Double.parseDouble(tokens[1]);
            double caxis = Double.parseDouble(tokens[2]);

            data = br.readLine().trim();
            tokens = data.split(" +");
            if (tokens.length != 3) {
                return false;
            }
            double alpha = Double.parseDouble(tokens[0]);
            double beta = Double.parseDouble(tokens[1]);
            double gamma = Double.parseDouble(tokens[2]);

            crystal.changeUnitCellParameters(aaxis, baxis, caxis, alpha, beta, gamma);

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
     * <p>
     * writeDYN</p>
     *
     * @param dynFile a {@link java.io.File} object.
     * @param x       an array of double.
     * @param v       an array of double.
     * @param a       an array of double.
     * @param ap      an array of double.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @return a boolean.
     */
    public boolean writeDYN(File dynFile, Crystal crystal, double[] x, double[] v, double[] a, double[] ap) {
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
            Crystal unitCell = crystal.getUnitCell();
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
