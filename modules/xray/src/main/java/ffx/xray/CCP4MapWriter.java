/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.xray;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;

/**
 *
 * @author Tim Fenn
 *
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">
 */
public class CCP4MapWriter {

    private static final Logger logger = Logger.getLogger(CCP4MapWriter.class.getName());
    private final String filename;
    private final Crystal crystal;
    private final int nx, ny, nz;

    public CCP4MapWriter(int nx, int ny, int nz, Crystal crystal,
            String filename) {
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.crystal = crystal;
        this.filename = filename;
    }

    public void write(double data[]) {
        ByteOrder b = ByteOrder.nativeOrder();
        FileOutputStream fos;
        DataOutputStream dos;

        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double mean = 0.0;
        double sd = 0.0;

        int n = 0;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    int index = 2 * (i + nx * (j + ny * k));
                    // int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
                    n++;
                    if (data[index] < min) {
                        min = data[index];
                    }
                    if (data[index] > max) {
                        max = data[index];
                    }
                    mean += (data[index] - mean) / n;
                }
            }
        }

        n = 0;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    int index = 2 * (i + nx * (j + ny * k));
                    // int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
                    sd += pow(data[index] - mean, 2.0);
                    n++;
                }
            }
        }
        sd = sqrt(sd / n);

        try {
            if (logger.isLoggable(Level.INFO)) {
                StringBuffer sb = new StringBuffer();
                sb.append(String.format("\nwriting CCP4 map file: \"%s\"\n", filename));
                sb.append(String.format("map min: %g max: %g mean: %g standard dev.: %g",
                        min, max, mean, sd));
                logger.info(sb.toString());
            }

            fos = new FileOutputStream(filename);
            dos = new DataOutputStream(fos);

            byte bytes[] = new byte[2048];
            int offset = 0;

            int imapdata;
            float fmapdata;
            String mapstr;

            // header
            ByteBuffer bb = ByteBuffer.wrap(bytes);
            bb.order(b).putInt(nx);
            bb.order(b).putInt(ny);
            bb.order(b).putInt(nz);

            // mode (2 = reals, only one we accept)
            bb.order(b).putInt(2);

            for (int i = 0; i < 3; i++) {
                bb.order(b).putInt(0);
            }
            bb.order(b).putInt(nx);
            bb.order(b).putInt(ny);
            bb.order(b).putInt(nz);

            bb.order(b).putFloat((float) crystal.a);
            bb.order(b).putFloat((float) crystal.b);
            bb.order(b).putFloat((float) crystal.c);
            bb.order(b).putFloat((float) crystal.alpha);
            bb.order(b).putFloat((float) crystal.beta);
            bb.order(b).putFloat((float) crystal.gamma);

            bb.order(b).putInt(1);
            bb.order(b).putInt(2);
            bb.order(b).putInt(3);

            bb.order(b).putFloat((float) min);
            bb.order(b).putFloat((float) max);
            bb.order(b).putFloat((float) mean);

            bb.order(b).putInt(crystal.spaceGroup.number);

            // symmetry bytes - should set this up at some point
            // imapdata = swap ? ByteSwap.swap(320) : 320;
            bb.order(b).putInt(0);

            bb.order(b).putInt(0);

            for (int i = 0; i < 12; i++) {
                bb.order(b).putFloat(0.0f);
            }

            for (int i = 0; i < 15; i++) {
                bb.order(b).putInt(0);
            }
            dos.write(bytes, offset, 208);
            bb.rewind();

            mapstr = new String("MAP ");
            dos.writeBytes(mapstr);

            // machine code: double, float, int, uchar
            // 0x4441 for LE, 0x1111 for BE
            if (ByteOrder.nativeOrder().equals(ByteOrder.LITTLE_ENDIAN)) {
                imapdata = 0x4441;
            } else {
                imapdata = 0x1111;
            }
            bb.order(b).putInt(imapdata);

            bb.order(b).putFloat((float) sd);

            bb.order(b).putInt(1);
            dos.write(bytes, offset, 12);

            StringBuffer sb = new StringBuffer();
            sb.append("map data from ffx");
            sb.setLength(80);
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append(" ");
            sb.setLength(80);
            for (int i = 0; i < 9; i++) {
                dos.writeBytes(sb.toString());
            }

            /*
            sb = new StringBuffer();
            sb.append("x,y,z");
            sb.setLength(80);
            dos.writeBytes(sb.toString());
             */

            bb.rewind();
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        int index = 2 * (i + nx * (j + ny * k));
                        // int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
                        fmapdata = (float) data[index];
                        bb.order(b).putFloat(fmapdata);
                        if (!bb.hasRemaining()){
                            dos.write(bytes);
                            bb.rewind();
                        }
                    }
                }
            }

            dos.close();
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }
}
