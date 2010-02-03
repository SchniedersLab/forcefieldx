/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.crystal;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * assumptions:
 * - output map covers unit cell (NOT ASU!)
 * - map is set up for passage into FFTs (so X axis has +2 offset)
 *
 * @author fennt
 */
public class CNSMapWriter {

    private static final Logger logger = Logger.getLogger(CNSMapWriter.class.getName());
    private final String filename;
    private final Crystal crystal;
    private final int nx, ny, nz;

    public CNSMapWriter(int nx, int ny, int nz, Crystal crystal,
            String filename) {
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.crystal = crystal;
        this.filename = filename;
    }

    public void write(double data[]) {
        try {
            if (logger.isLoggable(Level.INFO)) {
                StringBuffer sb = new StringBuffer();
                sb.append(String.format("\nwriting CNS map file: \"%s\"\n", filename));
                logger.info(sb.toString());
            }

            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
            double mean = 0.0;
            double sd = 0.0;
            out.println();
            out.println("       1");
            out.println("map from ffx");
            out.printf("%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
                    nx, 0, nx - 1, ny, 0, ny - 1, nz, 0, nz - 1);
            out.printf("%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
                    crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma);
            out.println("ZYX");
            int n = 0;
            for (int k = 0; k < nz; k++) {
                out.printf("%8d\n", k);
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
                        out.printf("%12.5E", data[index]);
                        n++;
                        mean += (data[index] - mean) / n;
                        if ((n % 6) == 0) {
                            out.println();
                        }
                    }
                }
            }

            n = 0;
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
                        sd += pow(data[index] - mean, 2.0);
                        n++;
                    }
                }
            }
            sd = sqrt(sd / n);

            out.println("   -9999");
            out.printf("%12.4E%12.4E\n", mean, sd);
            out.close();
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }
}
