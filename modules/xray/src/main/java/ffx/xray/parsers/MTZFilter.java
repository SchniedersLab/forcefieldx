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
package ffx.xray.parsers;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Double.parseDouble;
import static java.lang.Float.parseFloat;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.lang3.StringUtils;

import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toRadians;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;
import ffx.numerics.ComplexNumber;
import ffx.xray.DiffractionRefinementData;
import ffx.xray.parsers.MTZWriter.MTZType;

/**
 * This class parses CCP4 MTZ files.<br>
 *
 * @author Timothy D. Fenn<br>
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
public class MTZFilter implements DiffractionFileFilter {

    private static final Logger logger = Logger.getLogger(MTZFilter.class.getName());

    private class Column {

        public String label;
        public char type;
        public int id;
        public double min, max;
    }

    private class Dataset {

        public String project;
        public String dataset;
        public double lambda;
        public double[] cell = new double[6];
    }

    private static enum Header {

        VERS, TITLE, NCOL, SORT, SYMINF, SYMM, RESO, VALM, COL, COLUMN, NDIF,
        PROJECT, CRYSTAL, DATASET, DCELL, DWAVEL, BATCH,
        END, NOVALUE;

        public static Header toHeader(String str) {
            try {
                return valueOf(str);
            } catch (Exception ex) {
                return NOVALUE;
            }
        }
    }

    final private ArrayList<Column> columns = new ArrayList<>();
    final private ArrayList<Dataset> dataSets = new ArrayList<>();
    private boolean headerParsed = false;
    private String title;
    private String foString, sigFoString, rFreeString;

    private int h, k, l, fo, sigFo, rFree;
    private int fPlus, sigFPlus, fMinus, sigFMinus, rFreePlus, rFreeMinus;
    private int fc, phiC, fs, phiS;
    private int dsetOffset = 1;

    private int nColumns;
    private int nReflections;
    private int nBatches;
    private int spaceGroupNum;
    private String spaceGroupName;
    private double resLow;
    private double resHigh;

    /**
     * <p>
     * Constructor for MTZFilter.</p>
     */
    public MTZFilter() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ReflectionList getReflectionList(File mtzFile) {
        return getReflectionList(mtzFile, null);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ReflectionList getReflectionList(File mtzFile, CompositeConfiguration properties) {
        ByteOrder byteOrder = ByteOrder.nativeOrder();
        FileInputStream fileInputStream;
        DataInputStream dataInputStream;
        try {
            fileInputStream = new FileInputStream(mtzFile);
            dataInputStream = new DataInputStream(fileInputStream);

            byte headerOffset[] = new byte[4];
            byte bytes[] = new byte[80];
            int offset = 0;

            // Eat "MTZ" title.
            dataInputStream.read(bytes, offset, 4);
            String mtzstr = new String(bytes);

            // Header offset.
            dataInputStream.read(headerOffset, offset, 4);

            // Machine stamp.
            dataInputStream.read(bytes, offset, 4);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);
            int stamp = byteBuffer.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampstr = Integer.toHexString(stamp);
            switch (stampstr.charAt(0)) {
                case '1':
                case '3':
                    if (byteOrder.equals(ByteOrder.LITTLE_ENDIAN)) {
                        byteOrder = ByteOrder.BIG_ENDIAN;
                    }
                    break;
                case '4':
                    if (byteOrder.equals(ByteOrder.BIG_ENDIAN)) {
                        byteOrder = ByteOrder.LITTLE_ENDIAN;
                    }
                    break;
            }

            byteBuffer = ByteBuffer.wrap(headerOffset);
            int headerOffsetI = byteBuffer.order(byteOrder).getInt();

            // skip to header and parse
            dataInputStream.skipBytes((headerOffsetI - 4) * 4);

            for (Boolean parsing = true; parsing; dataInputStream.read(bytes, offset, 80)) {
                mtzstr = new String(bytes);
                parsing = parseHeader(mtzstr);
            }
        } catch (EOFException e) {
            String message = " MTZ end of file reached.";
            logger.log(Level.WARNING, message, e);
            return null;
        } catch (IOException e) {
            String message = " MTZ IO exception.";
            logger.log(Level.WARNING, message, e);
            return null;
        }

        // column identifiers
        foString = sigFoString = rFreeString = null;
        if (properties != null) {
            foString = properties.getString("fostring", null);
            sigFoString = properties.getString("sigfostring", null);
            rFreeString = properties.getString("rfreestring", null);
        }
        h = k = l = fo = sigFo = rFree = -1;
        fPlus = sigFPlus = fMinus = sigFMinus = rFreePlus = rFreeMinus = -1;
        fc = phiC = -1;
        boolean print = false;
        parseColumns(print);
        parseFcColumns(print);

        if (fo < 0 && fPlus < 0 && sigFo < 0 && sigFPlus < 0 && fc < 0 && phiC < 0) {
            logger.info(" The MTZ header contains insufficient information to generate the reflection list.");
            logger.info(" For non-default column labels set fostring/sigfostring in the properties file.");
            return null;
        }

        Column column;
        if (fo > 0) {
            column = (Column) columns.get(fo);
        } else if (fPlus > 0) {
            column = (Column) columns.get(fPlus);
        } else {
            column = (Column) columns.get(fc);
        }
        Dataset dataSet = (Dataset) dataSets.get(column.id - dsetOffset);

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n Reading %s\n\n", mtzFile.getName()));
            sb.append(format(" Setting up reflection list based on MTZ file.\n"));
            sb.append(format("  Space group number: %d (name: %s)\n",
                    spaceGroupNum, SpaceGroup.spaceGroupNames[spaceGroupNum - 1]));
            sb.append(format("  Resolution:         %8.3f\n", 0.999999 * resHigh));
            sb.append(format("  Cell:               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    dataSet.cell[0], dataSet.cell[1], dataSet.cell[2],
                    dataSet.cell[3], dataSet.cell[4], dataSet.cell[5]));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(dataSet.cell[0], dataSet.cell[1], dataSet.cell[2],
                dataSet.cell[3], dataSet.cell[4], dataSet.cell[5], SpaceGroup.spaceGroupNames[spaceGroupNum - 1]);

        double sampling = 0.6;
        if (properties != null) {
            sampling = properties.getDouble("sampling", 0.6);
        }
        Resolution resolution = new Resolution(0.999999 * resHigh, sampling);

        return new ReflectionList(crystal, resolution, properties);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getResolution(File mtzFile, Crystal crystal) {
        ReflectionList reflectionList = getReflectionList(mtzFile, null);
        return reflectionList.maxres;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readFile(File mtzFile, ReflectionList reflectionList,
            DiffractionRefinementData refinementData, CompositeConfiguration properties) {
        int nRead, nIgnore, nRes, nFriedel, nCut;
        ByteOrder byteOrder = ByteOrder.nativeOrder();
        FileInputStream fileInputStream;
        DataInputStream dataInputStream;
        boolean transpose = false;

        StringBuilder sb = new StringBuilder();
        try {
            fileInputStream = new FileInputStream(mtzFile);
            dataInputStream = new DataInputStream(fileInputStream);

            byte headerOffset[] = new byte[4];
            byte bytes[] = new byte[80];
            int offset = 0;

            // Eat "MTZ" title.
            dataInputStream.read(bytes, offset, 4);
            String mtzstr = null;

            // Header offset.
            dataInputStream.read(headerOffset, offset, 4);

            // Machine stamp.
            dataInputStream.read(bytes, offset, 4);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);
            int stamp = byteBuffer.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampString = Integer.toHexString(stamp);
            switch (stampString.charAt(0)) {
                case '1':
                case '3':
                    if (byteOrder.equals(ByteOrder.LITTLE_ENDIAN)) {
                        byteOrder = ByteOrder.BIG_ENDIAN;
                    }
                    break;
                case '4':
                    if (byteOrder.equals(ByteOrder.BIG_ENDIAN)) {
                        byteOrder = ByteOrder.LITTLE_ENDIAN;
                    }
                    break;
            }

            byteBuffer = ByteBuffer.wrap(headerOffset);
            int headerOffsetI = byteBuffer.order(byteOrder).getInt();

            // skip to header and parse
            dataInputStream.skipBytes((headerOffsetI - 4) * 4);

            for (Boolean parsing = true; parsing; dataInputStream.read(bytes, offset, 80)) {
                mtzstr = new String(bytes);
                parsing = parseHeader(mtzstr);
            }

            // column identifiers
            foString = sigFoString = rFreeString = null;
            if (properties != null) {
                foString = properties.getString("fostring", null);
                sigFoString = properties.getString("sigfostring", null);
                rFreeString = properties.getString("rfreestring", null);
            }
            h = k = l = fo = sigFo = rFree = -1;
            fPlus = sigFPlus = fMinus = sigFMinus = rFreePlus = rFreeMinus = -1;
            boolean print = true;
            parseColumns(print);

            if (h < 0 || k < 0 || l < 0) {
                String message = "Fatal error in MTZ file - no H K L indexes?\n";
                logger.log(Level.SEVERE, message);
                return false;
            }

            // Reopen to start at beginning.
            fileInputStream = new FileInputStream(mtzFile);
            dataInputStream = new DataInputStream(fileInputStream);

            // Skip initial header.
            dataInputStream.skipBytes(80);

            // Check if HKLs need to be transposed or not.
            float data[] = new float[nColumns];
            HKL mate = new HKL();
            int nPosIgnore = 0;
            int nTransIgnore = 0;
            int nZero = 0;
            int none = 0;
            for (int i = 0; i < nReflections; i++) {
                for (int j = 0; j < nColumns; j++) {
                    dataInputStream.read(bytes, offset, 4);
                    byteBuffer = ByteBuffer.wrap(bytes);
                    data[j] = byteBuffer.order(byteOrder).getFloat();
                }
                int ih = (int) data[h];
                int ik = (int) data[k];
                int il = (int) data[l];
                boolean friedel = reflectionList.findSymHKL(ih, ik, il, mate, false);
                HKL hklpos = reflectionList.getHKL(mate);
                if (hklpos == null) {
                    nPosIgnore++;
                }

                friedel = reflectionList.findSymHKL(ih, ik, il, mate, true);
                HKL hkltrans = reflectionList.getHKL(mate);
                if (hkltrans == null) {
                    nTransIgnore++;
                }
                if (rFree > 0) {
                    if (((int) data[rFree]) == 0) {
                        nZero++;
                    } else if (((int) data[rFree]) == 1) {
                        none++;
                    }
                }
                if (rFreePlus > 0) {
                    if (((int) data[rFreePlus]) == 0) {
                        nZero++;
                    } else if (((int) data[rFreePlus]) == 1) {
                        none++;
                    }
                }
                if (rFreeMinus > 0) {
                    if (((int) data[rFreeMinus]) == 0) {
                        nZero++;
                    } else if (((int) data[rFreeMinus]) == 1) {
                        none++;
                    }
                }
            }
            if (nPosIgnore > nTransIgnore) {
                transpose = true;
            }

            if (none > (nZero * 2) && refinementData.rfreeflag < 0) {
                refinementData.setFreeRFlag(0);
                sb.append(format(" Setting R free flag to %d based on MTZ file data.\n", refinementData.rfreeflag));
            } else if (nZero > (none * 2) && refinementData.rfreeflag < 0) {
                refinementData.setFreeRFlag(1);
                sb.append(format(" Setting R free flag to %d based on MTZ file data.\n", refinementData.rfreeflag));
            } else if (refinementData.rfreeflag < 0) {
                refinementData.setFreeRFlag(0);
                sb.append(format(" Setting R free flag to MTZ default: %d\n", refinementData.rfreeflag));
            }

            // reopen to start at beginning
            fileInputStream = new FileInputStream(mtzFile);
            dataInputStream = new DataInputStream(fileInputStream);

            // skip initial header
            dataInputStream.skipBytes(80);

            // read in data
            double anofSigF[][] = new double[refinementData.n][4];
            for (int i = 0; i < refinementData.n; i++) {
                anofSigF[i][0] = anofSigF[i][1] = anofSigF[i][2] = anofSigF[i][3] = Double.NaN;
            }
            nRead = nIgnore = nRes = nFriedel = nCut = 0;
            for (int i = 0; i < nReflections; i++) {
                for (int j = 0; j < nColumns; j++) {
                    dataInputStream.read(bytes, offset, 4);
                    byteBuffer = ByteBuffer.wrap(bytes);
                    data[j] = byteBuffer.order(byteOrder).getFloat();
                }
                int ih = (int) data[h];
                int ik = (int) data[k];
                int il = (int) data[l];
                boolean friedel = reflectionList.findSymHKL(ih, ik, il, mate, transpose);
                HKL hkl = reflectionList.getHKL(mate);
                if (hkl != null) {
                    if (fo > 0 && sigFo > 0) {
                        if (refinementData.fsigfcutoff > 0.0) {
                            if ((data[fo] / data[sigFo]) < refinementData.fsigfcutoff) {
                                nCut++;
                                continue;
                            }
                        }
                        if (friedel) {
                            anofSigF[hkl.index()][2] = data[fo];
                            anofSigF[hkl.index()][3] = data[sigFo];
                            nFriedel++;
                        } else {
                            anofSigF[hkl.index()][0] = data[fo];
                            anofSigF[hkl.index()][1] = data[sigFo];
                        }
                    } else {
                        if (fPlus > 0 && sigFPlus > 0) {
                            if (refinementData.fsigfcutoff > 0.0) {
                                if ((data[fPlus] / data[sigFPlus]) < refinementData.fsigfcutoff) {
                                    nCut++;
                                    continue;
                                }
                            }
                            anofSigF[hkl.index()][0] = data[fPlus];
                            anofSigF[hkl.index()][1] = data[sigFPlus];
                        }
                        if (fMinus > 0 && sigFMinus > 0) {
                            if (refinementData.fsigfcutoff > 0.0) {
                                if ((data[fMinus] / data[sigFMinus]) < refinementData.fsigfcutoff) {
                                    nCut++;
                                    continue;
                                }
                            }
                            anofSigF[hkl.index()][2] = data[fMinus];
                            anofSigF[hkl.index()][3] = data[sigFMinus];
                        }
                    }
                    if (rFree > 0) {
                        refinementData.setFreeR(hkl.index(), (int) data[rFree]);
                    } else {
                        if (rFreePlus > 0 && rFreeMinus > 0) {
                            // not sure what the correct thing to do here is?
                            refinementData.setFreeR(hkl.index(), (int) data[rFreePlus]);
                        } else if (rFreePlus > 0) {
                            refinementData.setFreeR(hkl.index(), (int) data[rFreePlus]);
                        } else if (rFreeMinus > 0) {
                            refinementData.setFreeR(hkl.index(), (int) data[rFreeMinus]);
                        }
                    }
                    nRead++;
                } else {
                    HKL tmp = new HKL(ih, ik, il);
                    if (!reflectionList.resolution.inInverseResSqRange(
                            Crystal.invressq(reflectionList.crystal, tmp))) {
                        nRes++;
                    } else {
                        nIgnore++;
                    }
                }
            }

            // Set up fsigf from F+ and F-.
            refinementData.generate_fsigf_from_anofsigf(anofSigF);

            // Log results.
            if (logger.isLoggable(Level.INFO)) {
                sb.append(format(" MTZ file type (machine stamp): %s\n", stampString));
                sb.append(format(" HKL data is %s\n", transpose ? "transposed" : "not transposed"));
                sb.append(format(" HKL read in:                             %d\n", nRead));
                sb.append(format(" HKL read as friedel mates:               %d\n", nFriedel));
                sb.append(format(" HKL NOT read in (too high resolution):   %d\n", nRes));
                sb.append(format(" HKL NOT read in (not in internal list?): %d\n", nIgnore));
                sb.append(format(" HKL NOT read in (F/sigF cutoff):         %d\n", nCut));
                sb.append(format(" HKL in internal list:                    %d\n", reflectionList.hkllist.size()));
                logger.info(sb.toString());
            }
            if (rFree < 0 && rFreePlus < 0 && rFreeMinus < 0) {
                refinementData.generateRFree();
            }
        } catch (EOFException e) {
            String message = " MTZ end of file reached.";
            logger.log(Level.WARNING, message, e);
            return false;
        } catch (IOException e) {
            String message = " MTZ IO Exception.";
            logger.log(Level.WARNING, message, e);
            return false;
        }

        return true;
    }

    /**
     * Average the computed structure factors for two systems.
     *
     * @param mtzFile1 This file will be overwritten and become the new average.
     * @param mtzFile2 Second MTZ file.
     * @param reflectionlist List of HKLs.
     * @param iter The iteration in the running average.
     * @param properties The CompositeConfiguration defines the properties of
     * each system.
     */
    public void averageFcs(File mtzFile1, File mtzFile2, ReflectionList reflectionlist,
            int iter, CompositeConfiguration properties) {

        DiffractionRefinementData fcdata1 = new DiffractionRefinementData(properties, reflectionlist);
        DiffractionRefinementData fcdata2 = new DiffractionRefinementData(properties, reflectionlist);

        readFcs(mtzFile1, reflectionlist, fcdata1, properties);
        readFcs(mtzFile2, reflectionlist, fcdata2, properties);

        // compute running average using mtzFile1 as current average
        System.out.println("iteration for averaging: " + iter);
        for (int i = 0; i < reflectionlist.hkllist.size(); i++) {
            fcdata1.fc[i][0] += (fcdata2.fc[i][0] - fcdata1.fc[i][0]) / iter;
            fcdata1.fc[i][1] += (fcdata2.fc[i][1] - fcdata1.fc[i][1]) / iter;

            fcdata1.fs[i][0] += (fcdata2.fs[i][0] - fcdata1.fs[i][0]) / iter;
            fcdata1.fs[i][1] += (fcdata2.fs[i][1] - fcdata1.fs[i][1]) / iter;
        }

        // overwrite original MTZ
        MTZWriter mtzOut = new MTZWriter(reflectionlist, fcdata1,
                mtzFile1.getName(), MTZType.FCONLY);
        mtzOut.write();
    }

    /**
     * Read the structure factors.
     *
     * @param mtzFile
     * @param reflectionList
     * @param fcData
     * @param properties
     * @return
     */
    public boolean readFcs(File mtzFile, ReflectionList reflectionList,
            DiffractionRefinementData fcData, CompositeConfiguration properties) {

        int nRead, nIgnore, nRes, nFriedel, nCut;
        ByteOrder byteOrder = ByteOrder.nativeOrder();
        FileInputStream fileInputStream;
        DataInputStream dataInputStream;

        StringBuilder sb = new StringBuilder();
        try {
            fileInputStream = new FileInputStream(mtzFile);
            dataInputStream = new DataInputStream(fileInputStream);

            byte headerOffset[] = new byte[4];
            byte bytes[] = new byte[80];
            int offset = 0;

            // Eat "MTZ" title.
            dataInputStream.read(bytes, offset, 4);
            String mtzString = null;

            // Header offset.
            dataInputStream.read(headerOffset, offset, 4);

            // Machine stamp.
            dataInputStream.read(bytes, offset, 4);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);
            int stamp = byteBuffer.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampString = Integer.toHexString(stamp);
            switch (stampString.charAt(0)) {
                case '1':
                case '3':
                    if (byteOrder.equals(ByteOrder.LITTLE_ENDIAN)) {
                        byteOrder = ByteOrder.BIG_ENDIAN;
                    }
                    break;
                case '4':
                    if (byteOrder.equals(ByteOrder.BIG_ENDIAN)) {
                        byteOrder = ByteOrder.LITTLE_ENDIAN;
                    }
                    break;
            }

            byteBuffer = ByteBuffer.wrap(headerOffset);
            int headerOffsetI = byteBuffer.order(byteOrder).getInt();

            // Skip to header and parse.
            dataInputStream.skipBytes((headerOffsetI - 4) * 4);

            for (Boolean parsing = true; parsing; dataInputStream.read(bytes, offset, 80)) {
                mtzString = new String(bytes);
                parsing = parseHeader(mtzString);
            }

            // Column identifiers.
            fc = phiC = fs = phiS = -1;
            boolean print = true;
            parseFcColumns(print);

            if (h < 0 || k < 0 || l < 0) {
                String message = " Fatal error in MTZ file - no H K L indexes?\n";
                logger.log(Level.SEVERE, message);
                return false;
            }

            // Reopen to start at beginning.
            fileInputStream = new FileInputStream(mtzFile);
            dataInputStream = new DataInputStream(fileInputStream);

            // Skip initial header.
            dataInputStream.skipBytes(80);

            float data[] = new float[nColumns];
            HKL mate = new HKL();

            // Read in data.
            ComplexNumber complexNumber = new ComplexNumber();
            nRead = nIgnore = nRes = nFriedel = nCut = 0;
            for (int i = 0; i < nReflections; i++) {
                for (int j = 0; j < nColumns; j++) {
                    dataInputStream.read(bytes, offset, 4);
                    byteBuffer = ByteBuffer.wrap(bytes);
                    data[j] = byteBuffer.order(byteOrder).getFloat();
                }
                int ih = (int) data[h];
                int ik = (int) data[k];
                int il = (int) data[l];
                boolean friedel = reflectionList.findSymHKL(ih, ik, il, mate, false);
                HKL hkl = reflectionList.getHKL(mate);

                if (hkl != null) {
                    if (fc > 0 && phiC > 0) {
                        complexNumber.re(data[fc] * cos(toRadians(data[phiC])));
                        complexNumber.im(data[fc] * sin(toRadians(data[phiC])));
                        fcData.setFc(hkl.index(), complexNumber);
                    }
                    if (fs > 0 && phiS > 0) {
                        complexNumber.re(data[fs] * cos(toRadians(data[phiS])));
                        complexNumber.im(data[fs] * sin(toRadians(data[phiS])));
                        fcData.setFs(hkl.index(), complexNumber);
                    }
                    nRead++;
                } else {
                    HKL tmp = new HKL(ih, ik, il);
                    if (!reflectionList.resolution.inInverseResSqRange(
                            Crystal.invressq(reflectionList.crystal, tmp))) {
                        nRes++;
                    } else {
                        nIgnore++;
                    }
                }
            }

            if (logger.isLoggable(Level.INFO)) {
                sb.append(format(" MTZ file type (machine stamp): %s\n", stampString));
                sb.append(format(" Fc HKL read in:                             %d\n", nRead));
                sb.append(format(" Fc HKL read as friedel mates:               %d\n", nFriedel));
                sb.append(format(" Fc HKL NOT read in (too high resolution):   %d\n", nRes));
                sb.append(format(" Fc HKL NOT read in (not in internal list?): %d\n", nIgnore));
                sb.append(format(" Fc HKL NOT read in (F/sigF cutoff):         %d\n", nCut));
                sb.append(format(" HKL in internal list:                       %d\n",
                        reflectionList.hkllist.size()));
                logger.info(sb.toString());
            }
        } catch (EOFException e) {
            String message = " MTZ end of file reached.";
            logger.log(Level.WARNING, message, e);
            return false;
        } catch (IOException e) {
            String message = " MTZ IO Exception.";
            logger.log(Level.WARNING, message, e);
            return false;
        }

        return true;
    }

    /**
     * Parse the header.
     *
     * @param str
     * @return
     */
    private Boolean parseHeader(String str) {
        Boolean parsing = true;
        Column column;
        Dataset dataSet;

        int nDataSet;
        String[] strArray = str.split("\\s+");

        if (headerParsed) {
            return Header.toHeader(strArray[0]) != Header.END;
        }

        switch (Header.toHeader(strArray[0])) {
            case TITLE:
                title = str.substring(5);
                break;
            case NCOL:
                nColumns = parseInt(strArray[1]);
                nReflections = parseInt(strArray[2]);
                nBatches = parseInt(strArray[3]);
                break;
            case SORT:
                break;
            case SYMINF:
                String[] tmp = str.split("\'+");
                spaceGroupNum = parseInt(strArray[4]);
                if (tmp.length > 1) {
                    spaceGroupName = tmp[1];
                }
                break;
            case SYMM:
                break;
            case RESO:
                double r1 = sqrt(1.0 / parseFloat(strArray[1]));
                double r2 = sqrt(1.0 / parseFloat(strArray[2]));
                resLow = max(r1, r2);
                resHigh = min(r1, r2);
                break;
            case VALM:
                break;
            case NDIF:
                int ndif = parseInt(strArray[1]);
                break;
            case COL:
            case COLUMN:
                nDataSet = parseInt(strArray[5]);
                if (nDataSet == 0) {
                    dsetOffset = 0;
                }
                column = new Column();
                columns.add(column);
                column.label = strArray[1];
                column.type = strArray[2].charAt(0);
                column.id = nDataSet;
                column.min = parseDouble(strArray[3]);
                column.max = parseDouble(strArray[4]);
                break;
            case PROJECT:
                nDataSet = parseInt(strArray[1]);
                if (nDataSet == 0) {
                    dsetOffset = 0;
                }
                try {
                    dataSet = (Dataset) dataSets.get(nDataSet - dsetOffset);
                } catch (IndexOutOfBoundsException e) {
                    dataSet = new Dataset();
                    dataSets.add(dataSet);
                }
                dataSet.project = strArray[2];
                break;
            case CRYSTAL:
                break;
            case DATASET:
                nDataSet = parseInt(strArray[1]);
                if (nDataSet == 0) {
                    dsetOffset = 0;
                }
                try {
                    dataSet = (Dataset) dataSets.get(nDataSet - dsetOffset);
                } catch (IndexOutOfBoundsException e) {
                    dataSet = new Dataset();
                    dataSets.add(dataSet);
                }
                dataSet.dataset = strArray[2];
                break;
            case DCELL:
                nDataSet = Integer.parseInt(strArray[1]);
                if (nDataSet == 0) {
                    dsetOffset = 0;
                }
                try {
                    dataSet = (Dataset) dataSets.get(nDataSet - dsetOffset);
                } catch (IndexOutOfBoundsException e) {
                    dataSet = new Dataset();
                    dataSets.add(dataSet);
                }
                dataSet.cell[0] = parseDouble(strArray[2]);
                dataSet.cell[1] = parseDouble(strArray[3]);
                dataSet.cell[2] = parseDouble(strArray[4]);
                dataSet.cell[3] = parseDouble(strArray[5]);
                dataSet.cell[4] = parseDouble(strArray[6]);
                dataSet.cell[5] = parseDouble(strArray[7]);
                break;
            case DWAVEL:
                nDataSet = parseInt(strArray[1]);
                if (nDataSet == 0) {
                    dsetOffset = 0;
                }
                try {
                    dataSet = (Dataset) dataSets.get(nDataSet - dsetOffset);
                } catch (IndexOutOfBoundsException e) {
                    dataSet = new Dataset();
                    dataSets.add(dataSet);
                }
                dataSet.lambda = parseDouble(strArray[2]);
                break;
            case BATCH:
                break;
            case END:
                headerParsed = true;
                parsing = false;
                break;
            default:
                break;
        }

        return parsing;
    }

    /**
     * Parse columns.
     *
     * @param print
     */
    private void parseColumns(boolean print) {

        int nc = 0;
        StringBuilder sb = new StringBuilder();
        for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            Column column = (Column) i.next();
            String label = column.label.trim();
            if (label.equalsIgnoreCase("H") && column.type == 'H') {
                h = nc;
            } else if (label.equalsIgnoreCase("K") && column.type == 'H') {
                k = nc;
            } else if (label.equalsIgnoreCase("L") && column.type == 'H') {
                l = nc;
            } else if ((label.equalsIgnoreCase("free")
                    || label.equalsIgnoreCase("freer")
                    || label.equalsIgnoreCase("freerflag")
                    || label.equalsIgnoreCase("freer_flag")
                    || label.equalsIgnoreCase("rfree")
                    || label.equalsIgnoreCase("rfreeflag")
                    || label.equalsIgnoreCase("r-free-flags")
                    || label.equalsIgnoreCase("test")
                    || StringUtils.equalsIgnoreCase(label, rFreeString))
                    && column.type == 'I') {
                sb.append(format(" Reading R Free column: \"%s\"\n", column.label));
                rFree = nc;
            } else if ((label.equalsIgnoreCase("free(+)")
                    || label.equalsIgnoreCase("freer(+)")
                    || label.equalsIgnoreCase("freerflag(+)")
                    || label.equalsIgnoreCase("freer_flag(+)")
                    || label.equalsIgnoreCase("rfree(+)")
                    || label.equalsIgnoreCase("rfreeflag(+)")
                    || label.equalsIgnoreCase("r-free-flags(+)")
                    || label.equalsIgnoreCase("test(+)")
                    || StringUtils.equalsIgnoreCase(label + "(+)", rFreeString))
                    && column.type == 'I') {
                rFreePlus = nc;
            } else if ((label.equalsIgnoreCase("free(-)")
                    || label.equalsIgnoreCase("freer(-)")
                    || label.equalsIgnoreCase("freerflag(-)")
                    || label.equalsIgnoreCase("freer_flag(-)")
                    || label.equalsIgnoreCase("rfree(-)")
                    || label.equalsIgnoreCase("rfreeflag(-)")
                    || label.equalsIgnoreCase("r-free-flags(-)")
                    || label.equalsIgnoreCase("test(-)")
                    || StringUtils.equalsIgnoreCase(label + "(-)", rFreeString))
                    && column.type == 'I') {
                rFreeMinus = nc;
            } else if ((label.equalsIgnoreCase("f")
                    || label.equalsIgnoreCase("fp")
                    || label.equalsIgnoreCase("fo")
                    || label.equalsIgnoreCase("fobs")
                    || label.equalsIgnoreCase("f-obs")
                    || StringUtils.equalsIgnoreCase(label, foString))
                    && column.type == 'F') {
                sb.append(format(" Reading Fo column: \"%s\"\n", column.label));
                fo = nc;
            } else if ((label.equalsIgnoreCase("f(+)")
                    || label.equalsIgnoreCase("fp(+)")
                    || label.equalsIgnoreCase("fo(+)")
                    || label.equalsIgnoreCase("fobs(+)")
                    || label.equalsIgnoreCase("f-obs(+)")
                    || StringUtils.equalsIgnoreCase(label + "(+)", foString))
                    && column.type == 'G') {
                fPlus = nc;
            } else if ((label.equalsIgnoreCase("f(-)")
                    || label.equalsIgnoreCase("fp(-)")
                    || label.equalsIgnoreCase("fo(-)")
                    || label.equalsIgnoreCase("fobs(-)")
                    || label.equalsIgnoreCase("f-obs(-)")
                    || StringUtils.equalsIgnoreCase(label + "(-)", foString))
                    && column.type == 'G') {
                fMinus = nc;
            } else if ((label.equalsIgnoreCase("sigf")
                    || label.equalsIgnoreCase("sigfp")
                    || label.equalsIgnoreCase("sigfo")
                    || label.equalsIgnoreCase("sigfobs")
                    || label.equalsIgnoreCase("sigf-obs")
                    || StringUtils.equalsIgnoreCase(label, sigFoString))
                    && column.type == 'Q') {
                sb.append(format(" Reading sigFo column: \"%s\"\n", column.label));
                sigFo = nc;
            } else if ((label.equalsIgnoreCase("sigf(+)")
                    || label.equalsIgnoreCase("sigfp(+)")
                    || label.equalsIgnoreCase("sigfo(+)")
                    || label.equalsIgnoreCase("sigfobs(+)")
                    || label.equalsIgnoreCase("sigf-obs(+)")
                    || StringUtils.equalsIgnoreCase(label + "(+)", sigFoString))
                    && column.type == 'L') {
                sigFPlus = nc;
            } else if ((label.equalsIgnoreCase("sigf(-)")
                    || label.equalsIgnoreCase("sigfp(-)")
                    || label.equalsIgnoreCase("sigfo(-)")
                    || label.equalsIgnoreCase("sigfobs(-)")
                    || label.equalsIgnoreCase("sigf-obs(-)")
                    || StringUtils.equalsIgnoreCase(label + "(-)", sigFoString))
                    && column.type == 'L') {
                sigFMinus = nc;
            }
        }
        if (fo < 0 && sigFo < 0
                && fPlus > 0 && sigFPlus > 0
                && fMinus > 0 && sigFMinus > 0) {
            sb.append(format(" Reading Fplus/Fminus column to fill in Fo\n"));
        }
        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    /**
     * Parse Fc columns.
     *
     * @param print
     */
    private void parseFcColumns(boolean print) {

        int nc = 0;
        StringBuilder sb = new StringBuilder();
        for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            Column column = (Column) i.next();
            String label = column.label.trim();
            if (label.equalsIgnoreCase("H") && column.type == 'H') {
                h = nc;
            } else if (label.equalsIgnoreCase("K") && column.type == 'H') {
                k = nc;
            } else if (label.equalsIgnoreCase("L") && column.type == 'H') {
                l = nc;
            } else if ((label.equalsIgnoreCase("fc")
                    || label.equalsIgnoreCase("fcalc"))
                    && column.type == 'F') {
                sb.append(format(" Reading Fc column: \"%s\"\n", column.label));
                fc = nc;
            } else if ((label.equalsIgnoreCase("phic")
                    || label.equalsIgnoreCase("phifc")
                    || label.equalsIgnoreCase("phicalc")
                    || label.equalsIgnoreCase("phifcalc"))
                    && column.type == 'P') {
                sb.append(format(" Reading phiFc column: \"%s\"\n", column.label));
                phiC = nc;
            } else if ((label.equalsIgnoreCase("fs")
                    || label.equalsIgnoreCase("fscalc"))
                    && column.type == 'F') {
                sb.append(format(" Reading Fs column: \"%s\"\n", column.label));
                fs = nc;
            } else if ((label.equalsIgnoreCase("phis")
                    || label.equalsIgnoreCase("phifs")
                    || label.equalsIgnoreCase("phiscalc")
                    || label.equalsIgnoreCase("phifscalc"))
                    && column.type == 'P') {
                sb.append(format(" Reading phiFs column: \"%s\"\n", column.label));
                phiS = nc;
            }
        }
        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    /**
     * <p>
     * printHeader</p>
     */
    public void printHeader() {
        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(" MTZ title: ").append(title).append("\n");
            sb.append(" MTZ space group: ").append(spaceGroupName).
                    append(" space group number: ").append(spaceGroupNum).append(" (").
                    append(SpaceGroup.spaceGroupNames[spaceGroupNum - 1]).append(")\n");
            sb.append(" MTZ resolution: ").append(resLow).append(" - ").
                    append(resHigh).append("\n");
            sb.append(" Number of reflections: ").append(nReflections).append("\n");

            int ndset = 1;
            for (Iterator i = dataSets.iterator(); i.hasNext(); ndset++) {
                Dataset d = (Dataset) i.next();
                sb.append("  dataset ").append(ndset).append(": ").append(d.dataset).append("\n");
                sb.append("  project ").append(ndset).append(": ").append(d.project).append("\n");
                sb.append("  wavelength ").append(ndset).append(": ").
                        append(d.lambda).append("\n");
                sb.append("  cell ").append(ndset).append(": ").append(d.cell[0]).
                        append(" ").append(d.cell[1]).append(" ").append(d.cell[2]).
                        append(" ").append(d.cell[3]).append(" ").append(d.cell[4]).
                        append(" ").append(d.cell[5]).append("\n");
                sb.append("\n");
            }

            sb.append(" Number of columns: ").append(nColumns).append("\n");
            int nc = 0;
            for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
                Column c = (Column) i.next();
                sb.append(String.format(
                        "  column %d: dataset id: %d min: %9.2f max: %9.2f label: %s type: %c\n",
                        nc, c.id, c.min, c.max, c.label, c.type));
            }
            logger.info(sb.toString());
        }
    }

    /**
     * @return the nColumns
     */
    public int getnColumns() {
        return nColumns;
    }

    /**
     * @return the nReflections
     */
    public int getnReflections() {
        return nReflections;
    }

    /**
     * @return the spaceGroupNum
     */
    public int getSpaceGroupNum() {
        return spaceGroupNum;
    }

}
