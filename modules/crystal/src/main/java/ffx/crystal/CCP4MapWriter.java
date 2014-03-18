/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.crystal;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 * <p>CCP4MapWriter class.</p>
 *
 * @author Tim Fenn
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4
 * map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html"
 * target="_blank">CCP4 library documentation</a>
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4
 * map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html"
 * target="_blank">CCP4 library documentation</a>
 *
 */
public class CCP4MapWriter {

    private static final Logger logger = Logger.getLogger(CCP4MapWriter.class.getName());
    private final String filename;
    private final Crystal crystal;
    private final int extx, exty, extz;
    private final int orix, oriy, oriz, nx, ny, nz;
    private int stride;

    /**
     * construct mapwriter object
     *
     * @param extx slices in x
     * @param exty slices in y
     * @param extz slices in z
     * @param crystal {@link Crystal} object
     * @param filename output filename
     */
    public CCP4MapWriter(int extx, int exty, int extz, Crystal crystal,
            String filename) {
        this.orix = 0;
        this.oriy = 0;
        this.oriz = 0;
        this.extx = extx;
        this.exty = exty;
        this.extz = extz;
        this.nx = extx;
        this.ny = exty;
        this.nz = extz;
        this.crystal = crystal;
        this.filename = filename;
        this.stride = 2;
    }

    /**
     * <p>Constructor for CCP4MapWriter.</p>
     *
     * @param orix a int.
     * @param oriy a int.
     * @param oriz a int.
     * @param extx a int.
     * @param exty a int.
     * @param extz a int.
     * @param nx a int.
     * @param ny a int.
     * @param nz a int.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param filename a {@link java.lang.String} object.
     */
    public CCP4MapWriter(int orix, int oriy, int oriz, int extx, int exty, int extz,
            int nx, int ny, int nz, Crystal crystal, String filename) {
        this.orix = orix;
        this.oriy = oriy;
        this.oriz = oriz;
        this.extx = extx;
        this.exty = exty;
        this.extz = extz;
        this.nx = nx;
        this.ny = ny;
        this.nz = nz;
        this.crystal = crystal;
        this.filename = filename;
        this.stride = 2;
    }

    /**
     * set the stepping across the array (e.g. 2 if data is separated by 1
     * space)
     *
     * @param stride the step size desired
     */
    public void setStride(int stride) {
        this.stride = stride;
    }

    /**
     * write data to file, does not normalize
     *
     * @param data map data to write out
     */
    public void write(double data[]) {
        write(data, false);
    }

    /**
     * write data to file, does not normalize
     *
     * @param data map data to write out
     * @param norm should the data be normalized by mean/sd?
     */
    public void write(double data[], boolean norm) {
        ByteOrder b = ByteOrder.nativeOrder();
        FileOutputStream fos;
        DataOutputStream dos;

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double mean = 0.0;
        double sd = 0.0;

        int n = 0;
        for (int k = 0; k < extz; k++) {
            for (int j = 0; j < exty; j++) {
                for (int i = 0; i < extx; i++) {
                    int index = stride * (i + extx * (j + exty * k));
                    // int index = k * (exty * (extx + 2)) + j * (extx + 2) + i;
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
        for (int k = 0; k < extz; k++) {
            for (int j = 0; j < exty; j++) {
                for (int i = 0; i < extx; i++) {
                    int index = stride * (i + extx * (j + exty * k));
                    // int index = k * (exty * (extx + 2)) + j * (extx + 2) + i;
                    sd += pow(data[index] - mean, 2.0);
                    n++;
                }
            }
        }
        sd = sqrt(sd / n);

        if (norm) {
            for (int k = 0; k < extz; k++) {
                for (int j = 0; j < exty; j++) {
                    for (int i = 0; i < extx; i++) {
                        int index = stride * (i + extx * (j + exty * k));
                        data[index] = (data[index] - mean) / sd;
                    }
                }
            }
            // recurse
            write(data, false);
        }

        try {
            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
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
            bb.order(b).putInt(extx);
            bb.order(b).putInt(exty);
            bb.order(b).putInt(extz);

            // mode (2 = reals, only one we accept)
            bb.order(b).putInt(2);

            bb.order(b).putInt(orix);
            bb.order(b).putInt(oriy);
            bb.order(b).putInt(oriz);
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
            // bb.order(b).putInt(1);

            // symmetry bytes - should set this up at some point
            // imapdata = swap ? ByteSwap.swap(320) : 320;
            bb.order(b).putInt(80);

            bb.order(b).putInt(0);

            for (int i = 0; i < 12; i++) {
                bb.order(b).putFloat(0.0f);
            }

            for (int i = 0; i < 15; i++) {
                bb.order(b).putInt(0);
            }
            dos.write(bytes, offset, 208);
            bb.rewind();

            mapstr = "MAP ";
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

            StringBuilder sb = new StringBuilder();
            sb.append("map data from ffx");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            while (sb.length() < 80) {
                sb.append(" ");
            }
            for (int i = 0; i < 9; i++) {
                dos.writeBytes(sb.toString());
            }

            sb = new StringBuilder();
            sb.append("x,y,z");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            bb.rewind();
            for (int k = 0; k < extz; k++) {
                for (int j = 0; j < exty; j++) {
                    for (int i = 0; i < extx; i++) {
                        int index = stride * (i + extx * (j + exty * k));
                        // int index = k * (exty * (extx + 2)) + j * (extx + 2) + i;
                        fmapdata = (float) data[index];
                        bb.order(b).putFloat(fmapdata);
                        if (!bb.hasRemaining()) {
                            dos.write(bytes);
                            bb.rewind();
                        }
                    }
                }
            }
            if (bb.position() > 0) {
                dos.write(bytes);
                bb.rewind();
            }

            dos.close();
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }
}
