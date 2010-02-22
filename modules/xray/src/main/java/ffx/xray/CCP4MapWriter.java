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
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.utilities.ByteSwap;

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
        Boolean swap = true;
        if (b.equals(ByteOrder.BIG_ENDIAN)) {
            swap = false;
        }
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
                    int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
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
                    int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
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
                logger.info(sb.toString());
            }

            fos = new FileOutputStream(filename);
            dos = new DataOutputStream(fos);

            byte bytes[] = new byte[80];
            int offset = 0;

            int imapdata;
            float fmapdata;
            String mapstr;

            // header
            imapdata = swap ? ByteSwap.swap(nx) : nx;
            dos.writeInt(imapdata);
            imapdata = swap ? ByteSwap.swap(ny) : ny;
            dos.writeInt(imapdata);
            imapdata = swap ? ByteSwap.swap(nz) : nz;
            dos.writeInt(imapdata);

            // mode (2 = reals, only one we accept)
            imapdata = swap ? ByteSwap.swap(2) : 2;
            dos.writeInt(imapdata);

            imapdata = swap ? ByteSwap.swap(0) : 0;
            for (int i = 0; i < 3; i++) {
                dos.writeInt(imapdata);
            }
            imapdata = swap ? ByteSwap.swap(nx) : nx;
            dos.writeInt(imapdata);
            imapdata = swap ? ByteSwap.swap(ny) : ny;
            dos.writeInt(imapdata);
            imapdata = swap ? ByteSwap.swap(nz) : nz;
            dos.writeInt(imapdata);

            fmapdata = (float) crystal.a;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) crystal.b;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) crystal.c;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) crystal.alpha;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) crystal.beta;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) crystal.gamma;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);

            imapdata = swap ? ByteSwap.swap(1) : 1;
            dos.writeInt(imapdata);
            imapdata = swap ? ByteSwap.swap(2) : 2;
            dos.writeInt(imapdata);
            imapdata = swap ? ByteSwap.swap(3) : 3;
            dos.writeInt(imapdata);

            fmapdata = (float) min;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) max;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);
            fmapdata = (float) mean;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);

            imapdata = swap ? ByteSwap.swap(crystal.spaceGroup.number)
                    : crystal.spaceGroup.number;
            dos.writeInt(imapdata);

            // symmetry bytes - should set this up at some point
            // imapdata = swap ? ByteSwap.swap(320) : 320;
            imapdata = swap ? ByteSwap.swap(0) : 0;
            dos.writeInt(imapdata);

            imapdata = swap ? ByteSwap.swap(0) : 0;
            dos.writeInt(imapdata);

            fmapdata = 0.0f;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            for (int i = 0; i < 12; i++) {
                dos.writeFloat(fmapdata);
            }

            for (int i = 0; i < 15; i++) {
                dos.writeInt(imapdata);
            }

            mapstr = new String("MAP ");
            dos.writeBytes(mapstr);

            // machine code: double, float, int, uchar
            // 0x4144 for LE, 0x1111 for BE
            imapdata = swap ? 0x4144 : 0x1111;
            imapdata = swap ? ByteSwap.swap(imapdata) : imapdata;
            dos.writeInt(imapdata);

            fmapdata = (float) sd;
            fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
            dos.writeFloat(fmapdata);

            imapdata = swap ? ByteSwap.swap(1) : 1;
            dos.writeInt(imapdata);

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

            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < ny; j++) {
                    for (int i = 0; i < nx; i++) {
                        int index = k * (ny * (nx + 2)) + j * (nx + 2) + i;
                        fmapdata = (float) data[index];
                        fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                        dos.writeFloat(fmapdata);
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
