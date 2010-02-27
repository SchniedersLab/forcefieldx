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

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;
import ffx.utilities.ByteSwap;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;

/**
 *
 * @author Tim Fenn
 *
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">
 */
public class MTZWriter {

    private static final Logger logger = Logger.getLogger(MTZWriter.class.getName());
    private final String filename;
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final Crystal crystal;
    private final SpaceGroup sg;
    private final int n;
    private final int ncol;

    public MTZWriter(ReflectionList reflectionlist,
            RefinementData refinementdata, String filename) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.sg = crystal.spaceGroup;
        this.filename = filename;
        // ignore 0 0 0 reflection
        this.n = refinementdata.n - 1;
        // static for now
        this.ncol = 16;
    }

    public void write() {
        ByteOrder b = ByteOrder.nativeOrder();
        Boolean swap = true;
        FileOutputStream fos;
        DataOutputStream dos;

        try {
            if (logger.isLoggable(Level.INFO)) {
                StringBuffer sb = new StringBuffer();
                sb.append(String.format("\nwriting MTZ HKL file: \"%s\"\n",
                        filename));
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
            StringBuffer sb = new StringBuffer();
            sb.append("MTZ ");
            dos.writeBytes(sb.toString());

            // header offset
            int headeroffset = n * ncol + 21;
            imapdata = swap ? ByteSwap.swap(headeroffset) : headeroffset;
            dos.writeInt(imapdata);

            // machine code: double, float, int, uchar
            // 0x4144 for LE, 0x1111 for BE
            if (b.equals(ByteOrder.LITTLE_ENDIAN)) {
                imapdata = 0x4144;
            } else {
                imapdata = 0x1111;
            }
            imapdata = swap ? ByteSwap.swap(imapdata) : imapdata;
            dos.writeInt(imapdata);

            sb = new StringBuffer();
            sb.append(" ");
            sb.setLength(68);
            dos.writeBytes(sb.toString());

            // data
            Vector<String> colname = new Vector<String>(ncol);
            char coltype[] = new char[ncol];
            double res[] = new double[2];
            res[0] = Double.MAX_VALUE;
            res[1] = Double.MIN_VALUE;
            float colminmax[][] = new float[ncol][2];
            for (int i = 0; i < ncol; i++) {
                colminmax[i][0] = Float.MAX_VALUE;
                colminmax[i][1] = Float.MIN_VALUE;
            }
            for (HKL ih : reflectionlist.hkllist) {
                int i = ih.index();

                // skip the 0 0 0 reflection
                if (ih.h() == 0
                        && ih.k() == 0
                        && ih.l() == 0) {
                    continue;
                }

                double r = Crystal.invressq(crystal, ih);
                res[0] = Math.min(r, res[0]);
                res[1] = Math.max(r, res[1]);

                // HKL first (3)
                colname.add("H");
                coltype[0] = 'H';
                fmapdata = ih.h();
                colminmax[0][0] = Math.min(fmapdata, colminmax[0][0]);
                colminmax[0][1] = Math.max(fmapdata, colminmax[0][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                colname.add("K");
                coltype[1] = 'H';
                fmapdata = ih.k();
                colminmax[1][0] = Math.min(fmapdata, colminmax[1][0]);
                colminmax[1][1] = Math.max(fmapdata, colminmax[1][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                colname.add("L");
                coltype[2] = 'H';
                fmapdata = ih.l();
                colminmax[2][0] = Math.min(fmapdata, colminmax[2][0]);
                colminmax[2][1] = Math.max(fmapdata, colminmax[2][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                // F/sigF (2)
                // this should be user definable!
                colname.add("FO");
                coltype[3] = 'F';
                fmapdata = (float) refinementdata.f(i);
                if (!Double.isNaN(fmapdata)) {
                    colminmax[3][0] = Math.min(fmapdata, colminmax[3][0]);
                    colminmax[3][1] = Math.max(fmapdata, colminmax[3][1]);
                }
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("SIGFO");
                coltype[4] = 'Q';
                fmapdata = (float) refinementdata.sigf(i);
                if (!Double.isNaN(fmapdata)) {
                    colminmax[4][0] = Math.min(fmapdata, colminmax[4][0]);
                    colminmax[4][1] = Math.max(fmapdata, colminmax[4][1]);
                }
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                // free R (1)
                colname.add("FreeR");
                coltype[5] = 'I';
                fmapdata = (float) refinementdata.freer(i);
                if (!Double.isNaN(fmapdata)) {
                    colminmax[5][0] = Math.min(fmapdata, colminmax[5][0]);
                    colminmax[5][1] = Math.max(fmapdata, colminmax[5][1]);
                }
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                // Fs (2)
                colname.add("Fs");
                coltype[6] = 'F';
                fmapdata = (float) refinementdata.fs_f(i);
                colminmax[6][0] = Math.min(fmapdata, colminmax[6][0]);
                colminmax[6][1] = Math.max(fmapdata, colminmax[6][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("PHIFs");
                coltype[7] = 'P';
                fmapdata = (float) Math.toDegrees(refinementdata.fs_phi(i));
                colminmax[7][0] = Math.min(fmapdata, colminmax[7][0]);
                colminmax[7][1] = Math.max(fmapdata, colminmax[7][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                // Fctot (2)
                colname.add("Fc");
                coltype[8] = 'F';
                fmapdata = (float) refinementdata.fctot_f(i);
                colminmax[8][0] = Math.min(fmapdata, colminmax[8][0]);
                colminmax[8][1] = Math.max(fmapdata, colminmax[8][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("PHIFc");
                coltype[9] = 'P';
                fmapdata = (float) Math.toDegrees(refinementdata.fctot_phi(i));
                colminmax[9][0] = Math.min(fmapdata, colminmax[9][0]);
                colminmax[9][1] = Math.max(fmapdata, colminmax[9][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                // FOM/phase (2)
                colname.add("FOM");
                coltype[10] = 'W';
                fmapdata = (float) refinementdata.fomphi[i][0];
                if (!Double.isNaN(fmapdata)) {
                    colminmax[10][0] = Math.min(fmapdata, colminmax[10][0]);
                    colminmax[10][1] = Math.max(fmapdata, colminmax[10][1]);
                }
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("PHIW");
                coltype[11] = 'P';
                fmapdata = (float) Math.toDegrees(refinementdata.fomphi[i][1]);
                colminmax[11][0] = Math.min(fmapdata, colminmax[11][0]);
                colminmax[11][1] = Math.max(fmapdata, colminmax[11][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);

                // map coeffs (4)
                colname.add("FWT");
                coltype[12] = 'F';
                fmapdata = (float) refinementdata.fofc2_f(i);
                colminmax[12][0] = Math.min(fmapdata, colminmax[12][0]);
                colminmax[12][1] = Math.max(fmapdata, colminmax[12][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("PHWT");
                coltype[13] = 'P';
                fmapdata = (float) Math.toDegrees(refinementdata.fofc2_phi(i));
                colminmax[13][0] = Math.min(fmapdata, colminmax[13][0]);
                colminmax[13][1] = Math.max(fmapdata, colminmax[13][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("DELFWT");
                coltype[14] = 'F';
                fmapdata = (float) refinementdata.fofc1_f(i);
                colminmax[14][0] = Math.min(fmapdata, colminmax[14][0]);
                colminmax[14][1] = Math.max(fmapdata, colminmax[14][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
                colname.add("PHDELWT");
                coltype[15] = 'P';
                fmapdata = (float) Math.toDegrees(refinementdata.fofc1_phi(i));
                colminmax[15][0] = Math.min(fmapdata, colminmax[15][0]);
                colminmax[15][1] = Math.max(fmapdata, colminmax[15][1]);
                fmapdata = swap ? ByteSwap.swap(fmapdata) : fmapdata;
                dos.writeFloat(fmapdata);
            }

            // header
            sb = new StringBuffer();
            sb.append("VERS MTZ:V1.1 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            Date now = new Date();
            SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
            sb = new StringBuffer();
            sb.append("TITLE FFX output: " + sdf.format(now));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append(String.format("NCOL %8d %12d %8d", ncol, n, 0));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("SORT    0    0    0    0    0 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            char cdata = sg.shortName.charAt(0);
            if (cdata == 'H') {
                cdata = 'R';
            }
            sb.append(String.format("SYMINF %3d %2d %c %5d %22s %5s",
                    sg.getNumberOfSymOps(),
                    sg.numPrimitiveSymEquiv,
                    cdata,
                    sg.number,
                    "'" + sg.shortName + "'",
                    sg.pointGroupName));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            for (int i = 0; i < sg.symOps.size(); i++) {
                sb = new StringBuffer();
                sb.append("SYMM ");
                SymOp symop = sg.symOps.get(i);
                sb.append(symop.toXYZString());
                while (sb.length() < 80) {
                    sb.append(" ");
                }
                dos.writeBytes(sb.toString());
            }

            sb = new StringBuffer();
            sb.append(String.format("RESO %8.6f%13s%8.6f",
                    res[0], " ", res[1]));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("VALM NAN ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("NDIF        1 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("PROJECT       1 project ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("CRYSTAL       1 crystal ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("DATASET       1 dataset ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            for (int j = 0; j < ncol; j++) {
                sb = new StringBuffer();
                sb.append(String.format("COLUMN %-30s %c %17.4f %17.4f    1",
                        colname.get(j), coltype[j],
                        colminmax[j][0], colminmax[j][1]));
                dos.writeBytes(sb.toString());
            }

            sb = new StringBuffer();
            sb.append(String.format("CELL %10.4f %9.4f %9.4f %9.4f %9.4f %9.4f ",
                    crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append(String.format("DCELL %9d %10.4f %9.4f %9.4f %9.4f %9.4f %9.4f ",
                    1, crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("DWAVEL        1    1.00000 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("END ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuffer();
            sb.append("MTZENDOFHEADERS ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            dos.close();
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }

    private static double getMin(double minmax[], double val) {
        if (val > minmax[1]) {
            return val;
        } else {
            return minmax[1];
        }
    }
}
