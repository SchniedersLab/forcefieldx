/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.SimpleDateFormat;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Date;
import java.util.Vector;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;

/**
 *
 * @author Tim Fenn
 *
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4 map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">CCP4 library documentation</a>
 */
public class MTZWriter {

    private static final Logger logger = Logger.getLogger(MTZWriter.class.getName());
    private final String filename;
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final Crystal crystal;
    private final SpaceGroup sg;
    private final ReflectionSpline spline;
    private final int n;
    private final int ncol;
    private final boolean dataonly;

    public MTZWriter(ReflectionList reflectionlist,
            RefinementData refinementdata, String filename) {
        this(reflectionlist, refinementdata, filename, false);
    }

    public MTZWriter(ReflectionList reflectionlist,
            RefinementData refinementdata, String filename, boolean dataonly) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.sg = crystal.spaceGroup;
        this.spline = new ReflectionSpline(reflectionlist,
                refinementdata.spline.length);
        this.filename = filename;
        // ignore 0 0 0 reflection
        this.n = refinementdata.n - 1;
        this.dataonly = dataonly;
        if (dataonly) {
            this.ncol = 6;
        } else {
            this.ncol = 18;
        }
    }

    public void write() {
        ByteOrder b = ByteOrder.nativeOrder();
        FileOutputStream fos;
        DataOutputStream dos;

        try {
            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
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
            StringBuilder sb = new StringBuilder();
            sb.append("MTZ ");
            dos.writeBytes(sb.toString());

            // header offset
            int headeroffset = n * ncol + 21;
            ByteBuffer bb = ByteBuffer.wrap(bytes);
            bb.order(b).putInt(headeroffset);

            // machine code: double, float, int, uchar
            // 0x4441 for LE, 0x1111 for BE
            if (ByteOrder.nativeOrder().equals(ByteOrder.LITTLE_ENDIAN)) {
                imapdata = 0x4441;
            } else {
                imapdata = 0x1111;
            }
            bb.order(b).putInt(imapdata);
            dos.write(bytes, offset, 8);

            sb = new StringBuilder();
            sb.append(" ");
            sb.setLength(68);
            dos.writeBytes(sb.toString());

            // data
            Vector<String> colname = new Vector<String>(ncol);
            char coltype[] = new char[ncol];
            double res[] = new double[2];
            res[0] = Double.POSITIVE_INFINITY;
            res[1] = Double.NEGATIVE_INFINITY;
            float colminmax[][] = new float[ncol][2];
            for (int i = 0; i < ncol; i++) {
                colminmax[i][0] = Float.POSITIVE_INFINITY;
                colminmax[i][1] = Float.NEGATIVE_INFINITY;
            }
            ReflectionSpline sigmaaspline = new ReflectionSpline(reflectionlist,
                    refinementdata.sigmaa.length);
            colname.add("H");
            coltype[0] = 'H';
            colname.add("K");
            coltype[1] = 'H';
            colname.add("L");
            coltype[2] = 'H';
            colname.add("FO");
            coltype[3] = 'F';
            colname.add("SIGFO");
            coltype[4] = 'Q';
            colname.add("FreeR");
            coltype[5] = 'I';
            if (!dataonly) {
                colname.add("Fs");
                coltype[6] = 'F';
                colname.add("PHIFs");
                coltype[7] = 'P';
                colname.add("Fc");
                coltype[8] = 'F';
                colname.add("PHIFc");
                coltype[9] = 'P';
                colname.add("FOM");
                coltype[10] = 'W';
                colname.add("PHIW");
                coltype[11] = 'P';
                colname.add("SigmaAs");
                coltype[12] = 'F';
                colname.add("SigmaAw");
                coltype[13] = 'Q';
                colname.add("FWT");
                coltype[14] = 'F';
                colname.add("PHWT");
                coltype[15] = 'P';
                colname.add("DELFWT");
                coltype[16] = 'F';
                colname.add("PHDELWT");
                coltype[17] = 'P';
            }

            for (HKL ih : reflectionlist.hkllist) {
                int i = ih.index();

                // skip the 0 0 0 reflection
                if (ih.h() == 0
                        && ih.k() == 0
                        && ih.l() == 0) {
                    continue;
                }

                // spline setup
                double ss = Crystal.invressq(crystal, ih);
                double fh = spline.f(ss, refinementdata.spline);
                double sa = sigmaaspline.f(ss, refinementdata.sigmaa);
                double wa = sigmaaspline.f(ss, refinementdata.sigmaw);
                res[0] = Math.min(ss, res[0]);
                res[1] = Math.max(ss, res[1]);

                // HKL first (3)
                fmapdata = ih.h();
                colminmax[0][0] = Math.min(fmapdata, colminmax[0][0]);
                colminmax[0][1] = Math.max(fmapdata, colminmax[0][1]);
                bb.rewind();
                bb.order(b).putFloat(fmapdata);

                fmapdata = ih.k();
                colminmax[1][0] = Math.min(fmapdata, colminmax[1][0]);
                colminmax[1][1] = Math.max(fmapdata, colminmax[1][1]);
                bb.order(b).putFloat(fmapdata);

                fmapdata = ih.l();
                colminmax[2][0] = Math.min(fmapdata, colminmax[2][0]);
                colminmax[2][1] = Math.max(fmapdata, colminmax[2][1]);
                bb.order(b).putFloat(fmapdata);

                // F/sigF (2)
                // FIXME: this should be user definable!
                fmapdata = (float) refinementdata.get_f(i);
                if (!Double.isNaN(fmapdata)) {
                    colminmax[3][0] = Math.min(fmapdata, colminmax[3][0]);
                    colminmax[3][1] = Math.max(fmapdata, colminmax[3][1]);
                } else {
                    sa = Double.NaN;
                    wa = Double.NaN;
                }
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) refinementdata.get_sigf(i);
                if (!Double.isNaN(fmapdata)) {
                    colminmax[4][0] = Math.min(fmapdata, colminmax[4][0]);
                    colminmax[4][1] = Math.max(fmapdata, colminmax[4][1]);
                }
                bb.order(b).putFloat(fmapdata);

                // free R (1)
                fmapdata = (float) refinementdata.get_freer(i);
                if (!Double.isNaN(fmapdata)) {
                    colminmax[5][0] = Math.min(fmapdata, colminmax[5][0]);
                    colminmax[5][1] = Math.max(fmapdata, colminmax[5][1]);
                }
                bb.order(b).putFloat(fmapdata);

                // if we're only writing out the data, stop here
                if (dataonly) {
                    dos.write(bytes, offset, 24);
                    continue;
                }

                // Fs (2)
                fmapdata = (float) refinementdata.fs_f(i);
                colminmax[6][0] = Math.min(fmapdata, colminmax[6][0]);
                colminmax[6][1] = Math.max(fmapdata, colminmax[6][1]);
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) Math.toDegrees(refinementdata.fs_phi(i));
                colminmax[7][0] = Math.min(fmapdata, colminmax[7][0]);
                colminmax[7][1] = Math.max(fmapdata, colminmax[7][1]);
                bb.order(b).putFloat(fmapdata);

                // Fctot (2)
                fmapdata = (float) refinementdata.fctot_f(i);
                colminmax[8][0] = Math.min(fmapdata, colminmax[8][0]);
                colminmax[8][1] = Math.max(fmapdata, colminmax[8][1]);
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) Math.toDegrees(refinementdata.fctot_phi(i));
                colminmax[9][0] = Math.min(fmapdata, colminmax[9][0]);
                colminmax[9][1] = Math.max(fmapdata, colminmax[9][1]);
                bb.order(b).putFloat(fmapdata);

                // FOM/phase (2)
                fmapdata = (float) refinementdata.fomphi[i][0];
                if (!Double.isNaN(fmapdata)) {
                    colminmax[10][0] = Math.min(fmapdata, colminmax[10][0]);
                    colminmax[10][1] = Math.max(fmapdata, colminmax[10][1]);
                }
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) Math.toDegrees(refinementdata.fomphi[i][1]);
                colminmax[11][0] = Math.min(fmapdata, colminmax[11][0]);
                colminmax[11][1] = Math.max(fmapdata, colminmax[11][1]);
                bb.order(b).putFloat(fmapdata);

                // sigmaA/w (2)
                fmapdata = (float) sa;
                if (!Double.isNaN(fmapdata)) {
                    colminmax[12][0] = Math.min(fmapdata, colminmax[12][0]);
                    colminmax[12][1] = Math.max(fmapdata, colminmax[12][1]);
                }
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) wa;
                if (!Double.isNaN(fmapdata)) {
                    colminmax[13][0] = Math.min(fmapdata, colminmax[13][0]);
                    colminmax[13][1] = Math.max(fmapdata, colminmax[13][1]);
                }
                bb.order(b).putFloat(fmapdata);

                // map coeffs (4)
                fmapdata = (float) refinementdata.fofc2_f(i);
                colminmax[14][0] = Math.min(fmapdata, colminmax[14][0]);
                colminmax[14][1] = Math.max(fmapdata, colminmax[14][1]);
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) Math.toDegrees(refinementdata.fofc2_phi(i));
                colminmax[15][0] = Math.min(fmapdata, colminmax[15][0]);
                colminmax[15][1] = Math.max(fmapdata, colminmax[15][1]);
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) refinementdata.fofc1_f(i);
                colminmax[16][0] = Math.min(fmapdata, colminmax[16][0]);
                colminmax[16][1] = Math.max(fmapdata, colminmax[16][1]);
                bb.order(b).putFloat(fmapdata);
                fmapdata = (float) Math.toDegrees(refinementdata.fofc1_phi(i));
                colminmax[17][0] = Math.min(fmapdata, colminmax[17][0]);
                colminmax[17][1] = Math.max(fmapdata, colminmax[17][1]);
                bb.order(b).putFloat(fmapdata);

                dos.write(bytes, offset, 72);
            }

            // header
            sb = new StringBuilder();
            sb.append("VERS MTZ:V1.1 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            Date now = new Date();
            SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
            sb = new StringBuilder();
            sb.append("TITLE FFX output: " + sdf.format(now));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append(String.format("NCOL %8d %12d %8d", ncol, n, 0));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("SORT    0    0    0    0    0 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
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
                sb = new StringBuilder();
                sb.append("SYMM ");
                SymOp symop = sg.symOps.get(i);
                sb.append(symop.toXYZString());
                while (sb.length() < 80) {
                    sb.append(" ");
                }
                dos.writeBytes(sb.toString());
            }

            sb = new StringBuilder();
            sb.append(String.format("RESO %8.6f%13s%8.6f",
                    res[0], " ", res[1]));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("VALM NAN ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("NDIF        1 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("PROJECT       1 project ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("CRYSTAL       1 crystal ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("DATASET       1 dataset ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            for (int j = 0; j < ncol; j++) {
                sb = new StringBuilder();
                sb.append(String.format("COLUMN %-30s %c %17.4f %17.4f    1",
                        colname.get(j), coltype[j],
                        colminmax[j][0], colminmax[j][1]));
                dos.writeBytes(sb.toString());
            }

            sb = new StringBuilder();
            sb.append(String.format("CELL %10.4f %9.4f %9.4f %9.4f %9.4f %9.4f ",
                    crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append(String.format("DCELL %9d %10.4f %9.4f %9.4f %9.4f %9.4f %9.4f ",
                    1, crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("DWAVEL        1    1.00000 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("END ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dos.writeBytes(sb.toString());

            sb = new StringBuilder();
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
}
