/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.xray.parsers;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;
import ffx.xray.DiffractionRefinementData;

/**
 * <p>
 * MTZWriter class.</p>
 *
 * @author Timothy D. Fenn
 *
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4
 * map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html"
 * target="_blank">CCP4 library documentation</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4
 * map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html"
 * target="_blank">CCP4 library documentation</a>
 */
public class MTZWriter {

    /**
     * The possible output "styles".
     */
    public static interface MTZType {

        /**
         * Output unscaled data only.
         */
        public static final int DATAONLY = 1;
        /**
         * Output unscaled Fcs only (still requires data to be read in).
         */
        public static final int FCONLY = 2;
        /**
         * Everything, including fitted/scaled coefficients (e.g. sigmaA, map
         * coefficients).
         */
        public static final int ALL = 3;
    }

    private static final Logger logger = Logger.getLogger(MTZWriter.class.getName());
    private final String fileName;
    private final ReflectionList reflectionList;
    private final DiffractionRefinementData refinementData;
    private final Crystal crystal;
    private final SpaceGroup spaceGroup;
    private final ReflectionSpline spline;
    private final int n;
    private final int nCol;
    private final int mtzType;

    /**
     * <p>
     * Constructor for MTZWriter.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param filename a {@link java.lang.String} object.
     */
    public MTZWriter(ReflectionList reflectionList,
            DiffractionRefinementData refinementData, String filename) {
        this(reflectionList, refinementData, filename, MTZType.ALL);
    }

    /**
     * <p>
     * Constructor for MTZWriter.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param filename a {@link java.lang.String} object.
     * @param mtzType see MTZType interface types
     */
    public MTZWriter(ReflectionList reflectionList,
            DiffractionRefinementData refinementData, String filename, int mtzType) {
        this.reflectionList = reflectionList;
        this.refinementData = refinementData;
        this.crystal = reflectionList.crystal;
        this.spaceGroup = crystal.spaceGroup;
        this.spline = new ReflectionSpline(reflectionList,
                refinementData.spline.length);
        this.fileName = filename;

        // Ignore 0 0 0 reflection.
        this.n = refinementData.n - 1;
        this.mtzType = mtzType;
        switch (mtzType) {
            case MTZType.DATAONLY:
                this.nCol = 6;
                break;
            case MTZType.FCONLY:
                this.nCol = 7;
                break;
            case MTZType.ALL:
                this.nCol = 18;
                break;
            default:
                this.nCol = 18;
                break;
        }
    }

    /**
     * <p>
     * write</p>
     */
    public void write() {
        ByteOrder byteOrder = ByteOrder.nativeOrder();
        FileOutputStream fileOutputStream;
        DataOutputStream dataOutputStream;

        try {
            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
                sb.append(format("\n Writing MTZ HKL file: \"%s\"\n", fileName));
                logger.info(sb.toString());
            }

            fileOutputStream = new FileOutputStream(fileName);
            dataOutputStream = new DataOutputStream(fileOutputStream);

            byte bytes[] = new byte[80];
            int offset = 0;
            int writeLen = 0;

            int iMapData;
            float fMapData;

            // Header.
            StringBuilder sb = new StringBuilder();
            sb.append("MTZ ");
            dataOutputStream.writeBytes(sb.toString());

            // Header offset.
            int headerOffset = n * nCol + 21;
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);
            byteBuffer.order(byteOrder).putInt(headerOffset);

            // Machine code: double, float, int, uchar
            // 0x4441 for LE, 0x1111 for BE
            if (ByteOrder.nativeOrder().equals(ByteOrder.LITTLE_ENDIAN)) {
                iMapData = 0x4441;
            } else {
                iMapData = 0x1111;
            }
            byteBuffer.order(byteOrder).putInt(iMapData);
            dataOutputStream.write(bytes, offset, 8);

            sb = new StringBuilder();
            sb.append(" ");
            sb.setLength(68);
            dataOutputStream.writeBytes(sb.toString());

            // Data.
            Vector<String> colname = new Vector<>(nCol);
            char colType[] = new char[nCol];
            double res[] = new double[2];
            res[0] = Double.POSITIVE_INFINITY;
            res[1] = Double.NEGATIVE_INFINITY;
            float colMinMax[][] = new float[nCol][2];
            for (int i = 0; i < nCol; i++) {
                colMinMax[i][0] = Float.POSITIVE_INFINITY;
                colMinMax[i][1] = Float.NEGATIVE_INFINITY;
            }
            ReflectionSpline sigmaASpline = new ReflectionSpline(reflectionList,
                    refinementData.sigmaa.length);
            int col = 0;

            colname.add("H");
            colType[col++] = 'H';
            colname.add("K");
            colType[col++] = 'H';
            colname.add("L");
            colType[col++] = 'H';
            writeLen += 12;

            if (mtzType != MTZType.FCONLY) {
                colname.add("FO");
                colType[col++] = 'F';
                colname.add("SIGFO");
                colType[col++] = 'Q';
                colname.add("FreeR");
                colType[col++] = 'I';
                writeLen += 12;
            }

            if (mtzType != MTZType.DATAONLY) {
                colname.add("Fs");
                colType[col++] = 'F';
                colname.add("PHIFs");
                colType[col++] = 'P';
                colname.add("Fc");
                colType[col++] = 'F';
                colname.add("PHIFc");
                colType[col++] = 'P';
                writeLen += 16;
            }

            if (mtzType == MTZType.ALL) {
                colname.add("FOM");
                colType[col++] = 'W';
                colname.add("PHIW");
                colType[col++] = 'P';
                colname.add("SigmaAs");
                colType[col++] = 'F';
                colname.add("SigmaAw");
                colType[col++] = 'Q';
                colname.add("FWT");
                colType[col++] = 'F';
                colname.add("PHWT");
                colType[col++] = 'P';
                colname.add("DELFWT");
                colType[col++] = 'F';
                colname.add("PHDELWT");
                colType[col++] = 'P';
                writeLen += 32;
            }

            for (HKL ih : reflectionList.hkllist) {
                col = 0;
                int i = ih.index();

                // Skip the 0 0 0 reflection.
                if (ih.h() == 0 && ih.k() == 0 && ih.l() == 0) {
                    continue;
                }

                double ss = Crystal.invressq(crystal, ih);
                res[0] = min(ss, res[0]);
                res[1] = max(ss, res[1]);

                // HKL first (3)
                fMapData = ih.h();
                colMinMax[col][0] = min(fMapData, colMinMax[0][0]);
                colMinMax[col][1] = max(fMapData, colMinMax[0][1]);
                byteBuffer.rewind();
                byteBuffer.order(byteOrder).putFloat(fMapData);
                col++;

                fMapData = ih.k();
                colMinMax[col][0] = min(fMapData, colMinMax[1][0]);
                colMinMax[col][1] = max(fMapData, colMinMax[1][1]);
                byteBuffer.order(byteOrder).putFloat(fMapData);
                col++;

                fMapData = ih.l();
                colMinMax[col][0] = min(fMapData, colMinMax[2][0]);
                colMinMax[col][1] = max(fMapData, colMinMax[2][1]);
                byteBuffer.order(byteOrder).putFloat(fMapData);
                col++;

                if (mtzType != MTZType.FCONLY) {
                    // F/sigF (2)
                    fMapData = (float) refinementData.getF(i);
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) refinementData.getSigF(i);
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;

                    // Free R (1)
                    fMapData = (float) refinementData.getFreeR(i);
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                }

                if (mtzType == MTZType.FCONLY) {
                    // Fs (2)
                    fMapData = (float) refinementData.fsF(i);
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) toDegrees(refinementData.fsPhi(i));
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;

                    // Fc (unscaled!) (2)
                    fMapData = (float) refinementData.fcF(i);
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) Math.toDegrees(refinementData.fcPhi(i));
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                }

                if (mtzType == MTZType.ALL) {
                    // Fs (2)
                    fMapData = (float) refinementData.fsF(i);
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) toDegrees(refinementData.fsPhi(i));
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;

                    // Fctot (2)
                    fMapData = (float) refinementData.fcTotF(i);
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) toDegrees(refinementData.fcTotPhi(i));
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = Math.min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = Math.max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;

                    // FOM/phase (2)
                    fMapData = (float) refinementData.fomphi[i][0];
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = Math.min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = Math.max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) toDegrees(refinementData.fomphi[i][1]);
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;

                    // Spline setup.
                    double fh = spline.f(ss, refinementData.spline);
                    double sa = sigmaASpline.f(ss, refinementData.sigmaa);
                    double wa = sigmaASpline.f(ss, refinementData.sigmaw);

                    // sigmaA/w (2)
                    fMapData = (float) sa;
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) wa;
                    if (!isNaN(fMapData)) {
                        colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                        colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    }
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;

                    // Map coeffs (4).
                    fMapData = (float) refinementData.FoFc2F(i);
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) toDegrees(refinementData.FoFc2Phi(i));
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) refinementData.foFc1F(i);
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                    fMapData = (float) toDegrees(refinementData.foFc1Phi(i));
                    colMinMax[col][0] = min(fMapData, colMinMax[col][0]);
                    colMinMax[col][1] = max(fMapData, colMinMax[col][1]);
                    byteBuffer.order(byteOrder).putFloat(fMapData);
                    col++;
                }

                dataOutputStream.write(bytes, offset, writeLen);
            }

            // Header.
            sb = new StringBuilder();
            sb.append("VERS MTZ:V1.1 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            Date now = new Date();
            SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
            sb = new StringBuilder();
            sb.append("TITLE FFX output: " + sdf.format(now));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append(String.format("NCOL %8d %12d %8d", nCol, n, 0));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("SORT    0    0    0    0    0 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            char cdata = spaceGroup.shortName.charAt(0);
            if (cdata == 'H') {
                cdata = 'R';
            }
            sb.append(String.format("SYMINF %3d %2d %c %5d %22s %5s",
                    spaceGroup.getNumberOfSymOps(),
                    spaceGroup.numPrimitiveSymEquiv,
                    cdata,
                    spaceGroup.number,
                    "'" + spaceGroup.shortName + "'",
                    spaceGroup.pointGroupName));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            for (int i = 0; i < spaceGroup.symOps.size(); i++) {
                sb = new StringBuilder();
                sb.append("SYMM ");
                SymOp symop = spaceGroup.symOps.get(i);
                sb.append(symop.toXYZString());
                while (sb.length() < 80) {
                    sb.append(" ");
                }
                dataOutputStream.writeBytes(sb.toString());
            }

            sb = new StringBuilder();
            sb.append(String.format("RESO %8.6f%13s%8.6f",
                    res[0], " ", res[1]));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("VALM NAN ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("NDIF        1 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("PROJECT       1 project ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("CRYSTAL       1 crystal ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("DATASET       1 dataset ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            for (int j = 0; j < nCol; j++) {
                sb = new StringBuilder();
                sb.append(String.format("COLUMN %-30s %c %17.4f %17.4f    1",
                        colname.get(j), colType[j],
                        colMinMax[j][0], colMinMax[j][1]));
                dataOutputStream.writeBytes(sb.toString());
            }

            sb = new StringBuilder();
            sb.append(String.format("CELL %10.4f %9.4f %9.4f %9.4f %9.4f %9.4f ",
                    crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append(String.format("DCELL %9d %10.4f %9.4f %9.4f %9.4f %9.4f %9.4f ",
                    1, crystal.a, crystal.b, crystal.c,
                    crystal.alpha, crystal.beta, crystal.gamma));
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("DWAVEL        1    1.00000 ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("END ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            sb = new StringBuilder();
            sb.append("MTZENDOFHEADERS ");
            while (sb.length() < 80) {
                sb.append(" ");
            }
            dataOutputStream.writeBytes(sb.toString());

            dataOutputStream.close();
        } catch (Exception e) {
            String message = "Fatal exception evaluating structure factors.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }
}
