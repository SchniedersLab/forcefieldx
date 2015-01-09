/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.xray;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;

/**
 * <p>
 * MTZWriter class.</p>
 *
 * @author Timothy D. Fenn
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
public class MTZWriter {

    /**
     * The possible output "styles"
     */
    public static interface MTZType {

        /**
         * output unscaled data only
         */
        public static final int DATAONLY = 1;
        /**
         * output unscaled Fcs only (still requires data to be read in)
         */
        public static final int FCONLY = 2;
        /**
         * everything, including fitted/scaled coefficients (e.g. sigmaA, map
         * coefficients)
         */
        public static final int ALL = 3;
    }

    private static final Logger logger = Logger.getLogger(MTZWriter.class.getName());
    private final String filename;
    private final ReflectionList reflectionlist;
    private final DiffractionRefinementData refinementdata;
    private final Crystal crystal;
    private final SpaceGroup sg;
    private final ReflectionSpline spline;
    private final int n;
    private final int ncol;
    private final int mtzType;

    /**
     * <p>
     * Constructor for MTZWriter.</p>
     *
     * @param reflectionlist a {@link ffx.crystal.ReflectionList} object.
     * @param refinementdata a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param filename a {@link java.lang.String} object.
     */
    public MTZWriter(ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata, String filename) {
        this(reflectionlist, refinementdata, filename, MTZType.ALL);
    }

    /**
     * <p>
     * Constructor for MTZWriter.</p>
     *
     * @param reflectionlist a {@link ffx.crystal.ReflectionList} object.
     * @param refinementdata a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param filename a {@link java.lang.String} object.
     * @param mtzType see MTZType interface types
     */
    public MTZWriter(ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata, String filename, int mtzType) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.sg = crystal.spaceGroup;
        this.spline = new ReflectionSpline(reflectionlist,
                refinementdata.spline.length);
        this.filename = filename;
        // ignore 0 0 0 reflection
        this.n = refinementdata.n - 1;
        this.mtzType = mtzType;
        switch (mtzType) {
            case MTZType.DATAONLY:
                this.ncol = 6;
                break;
            case MTZType.FCONLY:
                this.ncol = 7;
                break;
            case MTZType.ALL:
                this.ncol = 18;
                break;
            default:
                this.ncol = 18;
                break;
        }
    }

    /**
     * <p>
     * write</p>
     */
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
            int writelen = 0;

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
            int col = 0;

            colname.add("H");
            coltype[col++] = 'H';
            colname.add("K");
            coltype[col++] = 'H';
            colname.add("L");
            coltype[col++] = 'H';
            writelen += 12;

            if (mtzType != MTZType.FCONLY) {
                colname.add("FO");
                coltype[col++] = 'F';
                colname.add("SIGFO");
                coltype[col++] = 'Q';
                colname.add("FreeR");
                coltype[col++] = 'I';

                writelen += 12;
            }

            if (mtzType != MTZType.DATAONLY) {
                colname.add("Fs");
                coltype[col++] = 'F';
                colname.add("PHIFs");
                coltype[col++] = 'P';
                colname.add("Fc");
                coltype[col++] = 'F';
                colname.add("PHIFc");
                coltype[col++] = 'P';

                writelen += 16;
            }

            if (mtzType == MTZType.ALL) {
                colname.add("FOM");
                coltype[col++] = 'W';
                colname.add("PHIW");
                coltype[col++] = 'P';
                colname.add("SigmaAs");
                coltype[col++] = 'F';
                colname.add("SigmaAw");
                coltype[col++] = 'Q';
                colname.add("FWT");
                coltype[col++] = 'F';
                colname.add("PHWT");
                coltype[col++] = 'P';
                colname.add("DELFWT");
                coltype[col++] = 'F';
                colname.add("PHDELWT");
                coltype[col++] = 'P';

                writelen += 32;
            }

            for (HKL ih : reflectionlist.hkllist) {
                col = 0;
                int i = ih.index();

                // skip the 0 0 0 reflection
                if (ih.h() == 0
                        && ih.k() == 0
                        && ih.l() == 0) {
                    continue;
                }

                double ss = Crystal.invressq(crystal, ih);
                double sa = Double.NaN;
                double wa = Double.NaN;
                res[0] = Math.min(ss, res[0]);
                res[1] = Math.max(ss, res[1]);

                // HKL first (3)
                fmapdata = ih.h();
                colminmax[col][0] = Math.min(fmapdata, colminmax[0][0]);
                colminmax[col][1] = Math.max(fmapdata, colminmax[0][1]);
                bb.rewind();
                bb.order(b).putFloat(fmapdata);
                col++;

                fmapdata = ih.k();
                colminmax[col][0] = Math.min(fmapdata, colminmax[1][0]);
                colminmax[col][1] = Math.max(fmapdata, colminmax[1][1]);
                bb.order(b).putFloat(fmapdata);
                col++;

                fmapdata = ih.l();
                colminmax[col][0] = Math.min(fmapdata, colminmax[2][0]);
                colminmax[col][1] = Math.max(fmapdata, colminmax[2][1]);
                bb.order(b).putFloat(fmapdata);
                col++;

                if (mtzType != MTZType.FCONLY) {
                    // F/sigF (2)
                    fmapdata = (float) refinementdata.getF(i);
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    } else {
                        sa = Double.NaN;
                        wa = Double.NaN;
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) refinementdata.getSigF(i);
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;

                    // free R (1)
                    fmapdata = (float) refinementdata.getFreeR(i);
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                }

                if (mtzType == MTZType.FCONLY) {
                    // Fs (2)
                    fmapdata = (float) refinementdata.fsF(i);
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.fsPhi(i));
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;

                    // Fc (unscaled!) (2)
                    fmapdata = (float) refinementdata.fcF(i);
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.fcPhi(i));
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                }

                if (mtzType == MTZType.ALL) {
                    // Fs (2)
                    fmapdata = (float) refinementdata.fsF(i);
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.fsPhi(i));
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;

                    // Fctot (2)
                    fmapdata = (float) refinementdata.fcTotF(i);
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.fcTotPhi(i));
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;

                    // FOM/phase (2)
                    fmapdata = (float) refinementdata.fomphi[i][0];
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.fomphi[i][1]);
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;

                    // spline setup
                    double fh = spline.f(ss, refinementdata.spline);
                    sa = sigmaaspline.f(ss, refinementdata.sigmaa);
                    wa = sigmaaspline.f(ss, refinementdata.sigmaw);

                    // sigmaA/w (2)
                    fmapdata = (float) sa;
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) wa;
                    if (!Double.isNaN(fmapdata)) {
                        colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                        colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    }
                    bb.order(b).putFloat(fmapdata);
                    col++;

                    // map coeffs (4)
                    fmapdata = (float) refinementdata.FoFc2F(i);
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.FoFc2Phi(i));
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) refinementdata.foFc1F(i);
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;
                    fmapdata = (float) Math.toDegrees(refinementdata.foFc1Phi(i));
                    colminmax[col][0] = Math.min(fmapdata, colminmax[col][0]);
                    colminmax[col][1] = Math.max(fmapdata, colminmax[col][1]);
                    bb.order(b).putFloat(fmapdata);
                    col++;
                }

                dos.write(bytes, offset, writelen);
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
