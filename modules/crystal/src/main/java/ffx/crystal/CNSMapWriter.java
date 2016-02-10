/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.crystal;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The CNSMapWriter class writes an output map that covers the unit cell (not
 * the asymmetric unit). The map is set up for passage into FFTs (the x-axis has
 * +2 offset).
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
 *
 * @see CCP4MapWriter
 */
public class CNSMapWriter {

    private static final Logger logger = Logger.getLogger(CNSMapWriter.class.getName());
    private final String filename;
    private final Crystal crystal;
    private final int nx, ny, nz;

    /**
     * <p>
     * Constructor for CNSMapWriter.</p>
     *
     * @param nx a int.
     * @param ny a int.
     * @param nz a int.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param filename a {@link java.lang.String} object.
     */
    public CNSMapWriter(int nx, int ny, int nz, Crystal crystal,
            String filename) {
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.crystal = crystal;
        this.filename = filename;
    }

    /**
     * <p>
     * write</p>
     *
     * @param data an array of double.
     */
    public void write(double data[]) {
        try {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\nwriting CNS map file: \"%s\"\n", filename));

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
                        int index = 2 * (i + nx * (j + ny * k));
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
                        int index = 2 * (i + nx * (j + ny * k));
                        sd += pow(data[index] - mean, 2.0);
                        n++;
                    }
                }
            }
            sd = sqrt(sd / n);

            out.println("   -9999");
            out.printf("%12.4E%12.4E\n", mean, sd);
            sb.append(String.format("map mean: %g standard dev.: %g",
                    mean, sd));
            out.close();
            if (logger.isLoggable(Level.INFO)) {
                logger.info(sb.toString());
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }
}
