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

import org.apache.commons.configuration.CompositeConfiguration;

import java.io.File;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;

/**
 * This class parses CCP4 MTZ files.<br>
 * 
 * @author Tim Fenn<br>
 *
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4 map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">CCP4 library documentation</a>
 */
public class MTZFilter implements DiffractionFileFilter {

    private static final Logger logger = Logger.getLogger(MTZFilter.class.getName());

    private class column {

        public String label;
        public char type;
        public int id;
        public double min, max;
    }

    private class dataset {

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
    final private ArrayList<column> columns = new ArrayList();
    final private ArrayList<dataset> datasets = new ArrayList();
    private boolean headerparsed = false;
    private String title;
    private int h, k, l, fo, sigfo, rfree, fplus, sigfplus, fminus, sigfminus, rfreeplus, rfreeminus;
    private int dsetoffset = 1;
    public int ncol;
    public int nrfl;
    public int nbatches;
    public int sgnum;
    public String sgname;
    public double reslow;
    public double reshigh;

    // null constructor
    public MTZFilter() {
    }

    /*
    public boolean readFile(){
    this(readFile(molecularAssembly.getFile()));
    }
     */
    @Override
    public ReflectionList getReflectionList(File mtzFile) {
        return getReflectionList(mtzFile, null);
    }

    @Override
    public ReflectionList getReflectionList(File mtzFile, CompositeConfiguration properties) {
        ByteOrder b = ByteOrder.nativeOrder();
        FileInputStream fis;
        DataInputStream dis;
        try {
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            byte headeroffset[] = new byte[4];
            byte bytes[] = new byte[80];
            int offset = 0;

            // eat "MTZ" title
            dis.read(bytes, offset, 4);
            String mtzstr = new String(bytes);

            // header offset
            dis.read(headeroffset, offset, 4);

            // machine stamp
            dis.read(bytes, offset, 4);
            ByteBuffer bb = ByteBuffer.wrap(bytes);
            int stamp = bb.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampstr = Integer.toHexString(stamp);
            switch (stampstr.charAt(0)) {
                case '1':
                case '3':
                    if (b.equals(ByteOrder.LITTLE_ENDIAN)) {
                        b = ByteOrder.BIG_ENDIAN;
                    }
                    break;
                case '4':
                    if (b.equals(ByteOrder.BIG_ENDIAN)) {
                        b = ByteOrder.LITTLE_ENDIAN;
                    }
                    break;
            }

            bb = ByteBuffer.wrap(headeroffset);
            int headeroffseti = bb.order(b).getInt();

            // skip to header and parse
            dis.skipBytes((headeroffseti - 4) * 4);

            for (Boolean parsing = true; parsing; dis.read(bytes, offset, 80)) {
                mtzstr = new String(bytes);
                parsing = parse_header(mtzstr);
            }
        } catch (EOFException eof) {
            System.out.println("EOF reached ");
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return null;
        }

        h = k = l = fo = sigfo = rfree = -1;
        fplus = sigfplus = fminus = sigfminus = rfreeplus = rfreeminus = -1;
        parse_columns();

        if (fo < 0 && fplus < 0) {
            logger.info("insufficient information in MTZ header to generate Reflection List");
            return null;
        }

        column c;
        if (fo > 0) {
            c = (column) columns.get(fo);
        } else {
            c = (column) columns.get(fplus);
        }
        dataset d = (dataset) datasets.get(c.id - dsetoffset);

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\nOpening %s\n", mtzFile.getName()));
            sb.append(String.format("setting up Reflection List based on MTZ:\n"));
            sb.append(String.format("  spacegroup #: %d (name: %s)\n",
                    sgnum, SpaceGroup.spaceGroupNames[sgnum - 1]));
            sb.append(String.format("  resolution: %8.3f\n", 0.9999 * reshigh));
            sb.append(String.format("  cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    d.cell[0], d.cell[1], d.cell[2], d.cell[3], d.cell[4], d.cell[5]));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(d.cell[0], d.cell[1], d.cell[2],
                d.cell[3], d.cell[4], d.cell[5], SpaceGroup.spaceGroupNames[sgnum - 1]);

        double sampling = 1.0 / 1.5;
        if (properties != null) {
            sampling = properties.getDouble("sampling", 1.0 / 1.5);
        }
        Resolution resolution = new Resolution(0.9999 * reshigh, sampling);

        return new ReflectionList(crystal, resolution, properties);
    }

    @Override
    public boolean readFile(File mtzFile, ReflectionList reflectionlist,
            RefinementData refinementdata) {
        int nread, nignore, nres, nfriedel, ncut;
        ByteOrder b = ByteOrder.nativeOrder();
        FileInputStream fis;
        DataInputStream dis;
        boolean transpose = false;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\nOpening %s\n", mtzFile.getName()));
        try {
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            byte headeroffset[] = new byte[4];
            byte bytes[] = new byte[80];
            int offset = 0;

            // eat "MTZ" title
            dis.read(bytes, offset, 4);
            String mtzstr = new String(bytes);

            // header offset
            dis.read(headeroffset, offset, 4);

            // machine stamp
            dis.read(bytes, offset, 4);
            ByteBuffer bb = ByteBuffer.wrap(bytes);
            int stamp = bb.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampstr = Integer.toHexString(stamp);
            switch (stampstr.charAt(0)) {
                case '1':
                case '3':
                    if (b.equals(ByteOrder.LITTLE_ENDIAN)) {
                        b = ByteOrder.BIG_ENDIAN;
                    }
                    break;
                case '4':
                    if (b.equals(ByteOrder.BIG_ENDIAN)) {
                        b = ByteOrder.LITTLE_ENDIAN;
                    }
                    break;
            }

            bb = ByteBuffer.wrap(headeroffset);
            int headeroffseti = bb.order(b).getInt();

            // skip to header and parse
            dis.skipBytes((headeroffseti - 4) * 4);

            for (Boolean parsing = true; parsing; dis.read(bytes, offset, 80)) {
                mtzstr = new String(bytes);
                parsing = parse_header(mtzstr);
            }

            // column identifiers
            h = k = l = fo = sigfo = rfree = -1;
            fplus = sigfplus = fminus = sigfminus = rfreeplus = rfreeminus = -1;
            parse_columns();

            if (h < 0 || k < 0 || l < 0) {
                String message = "Fatal error in MTZ file - no H K L indexes?\n";
                logger.log(Level.SEVERE, message);
                return false;
            }

            // reopen to start at beginning
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            // skip initial header
            dis.skipBytes(80);

            // check if HKLs need to be transposed or not
            float data[] = new float[ncol];
            HKL mate = new HKL();
            int nposignore = 0;
            int ntransignore = 0;
            int nzero = 0;
            int none = 0;
            for (int i = 0; i < nrfl; i++) {
                for (int j = 0; j < ncol; j++) {
                    dis.read(bytes, offset, 4);
                    bb = ByteBuffer.wrap(bytes);
                    data[j] = bb.order(b).getFloat();
                }
                int ih = (int) data[h];
                int ik = (int) data[k];
                int il = (int) data[l];
                boolean friedel = reflectionlist.findSymHKL(ih, ik, il, mate, false);
                HKL hklpos = reflectionlist.getHKL(mate);
                if (hklpos == null) {
                    nposignore++;
                }

                friedel = reflectionlist.findSymHKL(ih, ik, il, mate, true);
                HKL hkltrans = reflectionlist.getHKL(mate);
                if (hkltrans == null) {
                    ntransignore++;
                }
                if (rfree > 0) {
                    if (((int) data[rfree]) == 0) {
                        nzero++;
                    } else if (((int) data[rfree]) == 1) {
                        none++;
                    }
                }
                if (rfreeplus > 0) {
                    if (((int) data[rfreeplus]) == 0) {
                        nzero++;
                    } else if (((int) data[rfreeplus]) == 1) {
                        none++;
                    }
                }
                if (rfreeminus > 0) {
                    if (((int) data[rfreeminus]) == 0) {
                        nzero++;
                    } else if (((int) data[rfreeminus]) == 1) {
                        none++;
                    }
                }
            }
            if (nposignore > ntransignore) {
                transpose = true;
            }

            if (none > (nzero * 2)
                    && refinementdata.rfreeflag < 0) {
                refinementdata.set_freerflag(0);
                sb.append(String.format("Setting R free flag to %d based on MTZ file data\n", refinementdata.rfreeflag));
            } else if (nzero > (none * 2)
                    && refinementdata.rfreeflag < 0) {
                refinementdata.set_freerflag(1);
                sb.append(String.format("Setting R free flag to %d based on MTZ file data\n", refinementdata.rfreeflag));
            } else if (refinementdata.rfreeflag < 0) {
                refinementdata.set_freerflag(0);
                sb.append(String.format("Setting R free flag to MTZ default: %d\n", refinementdata.rfreeflag));
            }

            // reopen to start at beginning
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            // skip initial header
            dis.skipBytes(80);

            // read in data
            nread = nignore = nres = nfriedel = ncut = 0;
            for (int i = 0; i < nrfl; i++) {
                for (int j = 0; j < ncol; j++) {
                    dis.read(bytes, offset, 4);
                    bb = ByteBuffer.wrap(bytes);
                    data[j] = bb.order(b).getFloat();
                }
                int ih = (int) data[h];
                int ik = (int) data[k];
                int il = (int) data[l];
                boolean friedel = reflectionlist.findSymHKL(ih, ik, il, mate, transpose);
                HKL hkl = reflectionlist.getHKL(mate);

                if (hkl != null) {
                    if (fo > 0 && sigfo > 0) {
                        if (refinementdata.fsigfcutoff > 0.0) {
                            if ((data[fo] / data[sigfo]) < refinementdata.fsigfcutoff) {
                                ncut++;
                                continue;
                            }
                        }
                        if (friedel) {
                            refinementdata.set_ano_fsigfminus(hkl.index(), data[fo], data[sigfo]);
                            nfriedel++;
                        } else {
                            refinementdata.set_ano_fsigfplus(hkl.index(), data[fo], data[sigfo]);
                        }
                    } else {
                        if (fplus > 0 && sigfplus > 0) {
                            if (refinementdata.fsigfcutoff > 0.0) {
                                if ((data[fplus] / data[sigfplus]) < refinementdata.fsigfcutoff) {
                                    ncut++;
                                    continue;
                                }
                            }
                            refinementdata.set_ano_fsigfplus(hkl.index(), data[fplus], data[sigfplus]);
                        }
                        if (fminus > 0 && sigfminus > 0) {
                            if (refinementdata.fsigfcutoff > 0.0) {
                                if ((data[fminus] / data[sigfminus]) < refinementdata.fsigfcutoff) {
                                    ncut++;
                                    continue;
                                }
                            }
                            refinementdata.set_ano_fsigfminus(hkl.index(), data[fminus], data[sigfminus]);
                        }
                    }
                    if (rfree > 0) {
                        refinementdata.set_freer(hkl.index(), (int) data[rfree]);
                    } else {
                        if (rfreeplus > 0 && rfreeminus > 0) {
                            // not sure what the correct thing to do here is?
                            refinementdata.set_freer(hkl.index(), (int) data[rfreeplus]);
                        } else if (rfreeplus > 0) {
                            refinementdata.set_freer(hkl.index(), (int) data[rfreeplus]);
                        } else if (rfreeminus > 0) {
                            refinementdata.set_freer(hkl.index(), (int) data[rfreeminus]);
                        }
                    }
                    nread++;
                } else {
                    HKL tmp = new HKL(ih, ik, il);
                    if (Crystal.invressq(reflectionlist.crystal, tmp)
                            > reflectionlist.resolution.invressq_limit()) {
                        nres++;
                    } else {
                        nignore++;
                    }
                }
            }

            // set up fsigf from F+ and F-
            refinementdata.generate_fsigf_from_anofsigf();

            sb.append(String.format("MTZ file type (machine stamp): %s\n",
                    stampstr));
            sb.append(String.format("HKL data is %s\n",
                    transpose ? "transposed" : "not transposed"));
            sb.append(String.format("# HKL read in:                             %d\n",
                    nread));
            sb.append(String.format("# HKL read as friedel mates:               %d\n",
                    nfriedel));
            sb.append(String.format("# HKL NOT read in (too high resolution):   %d\n",
                    nres));
            sb.append(String.format("# HKL NOT read in (not in internal list?): %d\n",
                    nignore));
            sb.append(String.format("# HKL NOT read in (F/sigF cutoff):         %d\n",
                    ncut));
            sb.append(String.format("# HKL in internal list:                    %d\n",
                    reflectionlist.hkllist.size()));
            if (logger.isLoggable(Level.INFO)) {
                logger.info(sb.toString());
            }

            if (rfree < 0 && rfreeplus < 0 && rfreeminus < 0) {
                refinementdata.generateRFree();
            }
        } catch (EOFException eof) {
            System.out.println("EOF reached ");
        } catch (IOException ioe) {
            System.out.println("IO Exception: " + ioe.getMessage());
            return false;
        }

        return true;
    }

    private Boolean parse_header(String str) {
        Boolean parsing = true;
        column col;
        dataset dset;

        int ndset;
        String[] strarray = str.split("\\s+");

        if (headerparsed) {
            if (Header.toHeader(strarray[0]) == Header.END) {
                return false;
            } else {
                return true;
            }
        }

        switch (Header.toHeader(strarray[0])) {
            case TITLE:
                title = str.substring(5);
                break;
            case NCOL:
                ncol = Integer.parseInt(strarray[1]);
                nrfl = Integer.parseInt(strarray[2]);
                nbatches = Integer.parseInt(strarray[3]);
                break;
            case SORT:
                break;
            case SYMINF:
                String[] tmp = str.split("\'+");
                sgnum = Integer.parseInt(strarray[4]);
                if (tmp.length > 1) {
                    sgname = tmp[1];
                }
                break;
            case SYMM:
                break;
            case RESO:
                double r1 = Math.sqrt(1.0 / Float.parseFloat(strarray[1]));
                double r2 = Math.sqrt(1.0 / Float.parseFloat(strarray[2]));
                reslow = Math.max(r1, r2);
                reshigh = Math.min(r1, r2);
                break;
            case VALM:
                break;
            case NDIF:
                int ndif = Integer.parseInt(strarray[1]);
                break;
            case COL:
            case COLUMN:
                ndset = Integer.parseInt(strarray[5]);
                if (ndset == 0) {
                    dsetoffset = 0;
                }
                col = new column();
                columns.add(col);
                col.label = strarray[1];
                col.type = strarray[2].charAt(0);
                col.id = ndset;
                col.min = Double.parseDouble(strarray[3]);
                col.max = Double.parseDouble(strarray[4]);
                break;
            case PROJECT:
                ndset = Integer.parseInt(strarray[1]);
                if (ndset == 0) {
                    dsetoffset = 0;
                }
                try {
                    dset = (dataset) datasets.get(ndset - dsetoffset);
                } catch (IndexOutOfBoundsException e) {
                    dset = new dataset();
                    datasets.add(dset);
                }
                dset.project = strarray[2];
                break;
            case CRYSTAL:
                break;
            case DATASET:
                ndset = Integer.parseInt(strarray[1]);
                if (ndset == 0) {
                    dsetoffset = 0;
                }
                try {
                    dset = (dataset) datasets.get(ndset - dsetoffset);
                } catch (IndexOutOfBoundsException e) {
                    dset = new dataset();
                    datasets.add(dset);
                }
                dset.dataset = strarray[2];
                break;
            case DCELL:
                ndset = Integer.parseInt(strarray[1]);
                if (ndset == 0) {
                    dsetoffset = 0;
                }
                try {
                    dset = (dataset) datasets.get(ndset - dsetoffset);
                } catch (IndexOutOfBoundsException e) {
                    dset = new dataset();
                    datasets.add(dset);
                }
                dset.cell[0] = Double.parseDouble(strarray[2]);
                dset.cell[1] = Double.parseDouble(strarray[3]);
                dset.cell[2] = Double.parseDouble(strarray[4]);
                dset.cell[3] = Double.parseDouble(strarray[5]);
                dset.cell[4] = Double.parseDouble(strarray[6]);
                dset.cell[5] = Double.parseDouble(strarray[7]);
                break;
            case DWAVEL:
                ndset = Integer.parseInt(strarray[1]);
                if (ndset == 0) {
                    dsetoffset = 0;
                }
                try {
                    dset = (dataset) datasets.get(ndset - dsetoffset);
                } catch (IndexOutOfBoundsException e) {
                    dset = new dataset();
                    datasets.add(dset);
                }
                dset.lambda = Double.parseDouble(strarray[2]);
                break;
            case BATCH:
                break;
            case END:
                headerparsed = true;
                parsing = false;
                break;
            default:
                break;
        }

        return parsing;
    }

    private void parse_columns() {

        // TODO: allow user to set mtz strings to look for in properties
        int nc = 0;
        StringBuilder sb = new StringBuilder();
        for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            column c = (column) i.next();
            String label = c.label.trim();
            if (label.equalsIgnoreCase("H") && c.type == 'H') {
                h = nc;
            } else if (label.equalsIgnoreCase("K") && c.type == 'H') {
                k = nc;
            } else if (label.equalsIgnoreCase("L") && c.type == 'H') {
                l = nc;
            } else if ((label.equalsIgnoreCase("free")
                    || label.equalsIgnoreCase("freer")
                    || label.equalsIgnoreCase("freerflag")
                    || label.equalsIgnoreCase("freer_flag")
                    || label.equalsIgnoreCase("rfree")
                    || label.equalsIgnoreCase("rfreeflag")
                    || label.equalsIgnoreCase("r-free-flags")
                    || label.equalsIgnoreCase("test"))
                    && c.type == 'I') {
                sb.append(String.format("Reading R Free column: \"%s\"\n", c.label));
                rfree = nc;
            } else if ((label.equalsIgnoreCase("free(+)")
                    || label.equalsIgnoreCase("freer(+)")
                    || label.equalsIgnoreCase("freerflag(+)")
                    || label.equalsIgnoreCase("freer_flag(+)")
                    || label.equalsIgnoreCase("rfree(+)")
                    || label.equalsIgnoreCase("rfreeflag(+)")
                    || label.equalsIgnoreCase("r-free-flags(+)")
                    || label.equalsIgnoreCase("test(+)"))
                    && c.type == 'I') {
                rfreeplus = nc;
            } else if ((label.equalsIgnoreCase("free(-)")
                    || label.equalsIgnoreCase("freer(-)")
                    || label.equalsIgnoreCase("freerflag(-)")
                    || label.equalsIgnoreCase("freer_flag(-)")
                    || label.equalsIgnoreCase("rfree(-)")
                    || label.equalsIgnoreCase("rfreeflag(-)")
                    || label.equalsIgnoreCase("r-free-flags(-)")
                    || label.equalsIgnoreCase("test(-)"))
                    && c.type == 'I') {
                rfreeminus = nc;
            } else if ((label.equalsIgnoreCase("f")
                    || label.equalsIgnoreCase("fp")
                    || label.equalsIgnoreCase("fo")
                    || label.equalsIgnoreCase("fobs"))
                    && c.type == 'F') {
                sb.append(String.format("Reading Fo column: \"%s\"\n", c.label));
                fo = nc;
            } else if ((label.equalsIgnoreCase("f(+)")
                    || label.equalsIgnoreCase("fp(+)")
                    || label.equalsIgnoreCase("fo(+)")
                    || label.equalsIgnoreCase("fobs(+)"))
                    && c.type == 'G') {
                fplus = nc;
            } else if ((label.equalsIgnoreCase("f(-)")
                    || label.equalsIgnoreCase("fp(-)")
                    || label.equalsIgnoreCase("fo(-)")
                    || label.equalsIgnoreCase("fobs(-)"))
                    && c.type == 'G') {
                fminus = nc;
            } else if ((label.equalsIgnoreCase("sigf")
                    || label.equalsIgnoreCase("sigfp")
                    || label.equalsIgnoreCase("sigfo")
                    || label.equalsIgnoreCase("sigfobs"))
                    && c.type == 'Q') {
                sb.append(String.format("Reading sigFo column: \"%s\"\n", c.label));
                sigfo = nc;
            } else if ((label.equalsIgnoreCase("sigf(+)")
                    || label.equalsIgnoreCase("sigfp(+)")
                    || label.equalsIgnoreCase("sigfo(+)")
                    || label.equalsIgnoreCase("sigfobs(+)"))
                    && c.type == 'L') {
                sigfplus = nc;
            } else if ((label.equalsIgnoreCase("sigf(-)")
                    || label.equalsIgnoreCase("sigfp(-)")
                    || label.equalsIgnoreCase("sigfo(-)")
                    || label.equalsIgnoreCase("sigfobs(-)"))
                    && c.type == 'L') {
                sigfminus = nc;
            }
        }
        if (fo < 0 && sigfo < 0
                && fplus > 0 && sigfplus > 0
                && fminus > 0 && sigfminus > 0) {
            sb.append(String.format("Reading Fplus/Fminus column to fill in Fo\n"));
        }
        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }
    }

    public void print_header() {
        StringBuilder sb = new StringBuilder();
        sb.append("MTZ title: " + title + "\n");
        sb.append("MTZ space group: " + sgname + " space group number: " + sgnum
                + " (" + SpaceGroup.spaceGroupNames[sgnum - 1] + ")\n");
        sb.append("MTZ resolution: " + reslow + " - " + reshigh + "\n");
        sb.append("# reflections: " + nrfl + "\n");

        int ndset = 1;
        for (Iterator i = datasets.iterator(); i.hasNext(); ndset++) {
            dataset d = (dataset) i.next();

            sb.append("  dataset " + ndset + ": " + d.dataset + "\n");
            sb.append("  project " + ndset + ": " + d.project + "\n");
            sb.append("  wavelength " + ndset + ": " + d.lambda + "\n");
            sb.append("  cell " + ndset + ": "
                    + d.cell[0] + " " + d.cell[1] + " " + d.cell[2] + " "
                    + d.cell[3] + " " + d.cell[4] + " " + d.cell[5] + "\n");
            sb.append("\n");
        }

        sb.append("# columns: " + ncol + "\n");
        int nc = 0;
        for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            column c = (column) i.next();

            sb.append(String.format("  column %d: dataset id: %d min: %9.2f max: %9.2f label: %s type: %c\n",
                    nc, c.id, c.min, c.max, c.label, c.type));
        }

        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }
    }
}
