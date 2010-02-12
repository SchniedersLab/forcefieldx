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

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;
import ffx.utilities.ByteSwap;

import java.io.File;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author fennt
 */
public class MTZFilter {

    private static final Logger logger = Logger.getLogger(MTZFilter.class.getName());

    private class column {

        public String label;
        public char type;
        public int id;
    }

    private class dataset {

        public String project;
        public String dataset;
        public double lambda;
        public double[] cell = new double[6];
    }

    private enum Header {

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
    private String title;
    private int h, k, l, fo, sigfo, rfree;
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
    public ReflectionList getReflectionList(File mtzFile) {
        ByteOrder b = ByteOrder.nativeOrder();
        Boolean swap = true;
        if (b.equals(ByteOrder.BIG_ENDIAN)) {
            swap = false;
        }
        FileInputStream fis;
        DataInputStream dis;
        try {
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            byte bytes[] = new byte[80];
            int offset = 0;
            String mtzstr = new String(bytes);

            // eat "MTZ" title
            dis.read(bytes, offset, 4);

            // header offset
            int headeroffset = swap ? ByteSwap.swap(dis.readInt()) : dis.readInt();

            // ignore machine stamp
            dis.read(bytes, offset, 4);
            mtzstr = new String(bytes);

            // skip to header and parse
            dis.skipBytes((headeroffset - 4) * 4);

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
        parse_columns();

        if (fo < 0) {
            return null;
        }

        column c = (column) columns.get(fo);
        dataset d = (dataset) datasets.get(c.id - 1);

        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            sb.append(String.format("\nsetting up Reflection List using spacegroup #: %d (name: %s)\n",
                    sgnum, SpaceGroup.spaceGroupNames[sgnum - 1]));
            sb.append(String.format("and cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                    d.cell[0], d.cell[1], d.cell[2], d.cell[3], d.cell[4], d.cell[5]));
            logger.info(sb.toString());
        }

        Crystal crystal = new Crystal(d.cell[0], d.cell[1], d.cell[2],
                d.cell[3], d.cell[4], d.cell[5], SpaceGroup.spaceGroupNames[sgnum - 1]);
        Resolution resolution = new Resolution(reshigh);

        return new ReflectionList(crystal, resolution);
    }

    public boolean readFile(File mtzFile, ReflectionList reflectionlist,
            RefinementData refinementdata) {
        ByteOrder b = ByteOrder.nativeOrder();
        Boolean swap = true;
        if (b.equals(ByteOrder.BIG_ENDIAN)) {
            swap = false;
        }
        FileInputStream fis;
        DataInputStream dis;
        try {
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            byte bytes[] = new byte[80];
            int offset = 0;
            String mtzstr = new String(bytes);

            // eat "MTZ" title
            dis.read(bytes, offset, 4);

            // header offset
            int headeroffset = swap ? ByteSwap.swap(dis.readInt()) : dis.readInt();

            // ignore machine stamp
            dis.read(bytes, offset, 4);
            mtzstr = new String(bytes);

            // skip to header and parse
            dis.skipBytes((headeroffset - 4) * 4);

            for (Boolean parsing = true; parsing; dis.read(bytes, offset, 80)) {
                mtzstr = new String(bytes);
                parsing = parse_header(mtzstr);
            }

            // reopen to start at beginning
            fis = new FileInputStream(mtzFile);
            dis = new DataInputStream(fis);

            // skip initial header
            dis.skipBytes(80);

            // column identifiers
            h = k = l = fo = sigfo = rfree = -1;
            parse_columns();

            if (h < 0 || k < 0 || l < 0) {
                String message = "Fatal error in MTZ file - no H K L indexes?\n";
                logger.log(Level.SEVERE, message);
                return false;
            }

            // read in data
            int nread, nignore, nres;
            nread = nignore = nres = 0;
            float data[] = new float[ncol];
            for (int i = 0; i < nrfl; i++) {
                for (int j = 0; j < ncol; j++) {
                    data[j] = swap ? ByteSwap.swap(dis.readFloat()) : dis.readFloat();
                }
                int ih = (int) data[h];
                int ik = (int) data[k];
                int il = (int) data[l];
                HKL hkl = reflectionlist.getHKL(ih, ik, il);
                if (hkl != null) {
                    if (fo > 0 && sigfo > 0) {
                        refinementdata.fsigf(hkl.index(), data[fo], data[sigfo]);
                    }
                    if (rfree > 0) {
                        refinementdata.freer(hkl.index(), (int) data[rfree]);
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
            StringBuffer sb = new StringBuffer();
            sb.append(String.format("\n# HKL read in:                             %d\n", nread));
            sb.append(String.format("# HKL NOT read in (too high resolution):   %d\n", nres));
            sb.append(String.format("# HKL NOT read in (not in internal list?): %d\n", nignore));
            if (logger.isLoggable(Level.INFO)) {
                logger.info(sb.toString());
            }

            // print out some reflections
            /*
             * TODO: columns.id should match desired column to read in
             * then read in cell from datasets.dataset == columns.id
             * into datasets.cell
             * flag desired columns to read in, then read in data
             */
            /*
            int nc = 1;
            for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            column c = (column) i.next();

            System.out.format("%3d %10s ", nc, c.label);
            }
            System.out.println();
            nc = 1;
            for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            column c = (column) i.next();

            System.out.format("%3d %10s ", nc, c.type);
            }
            System.out.println();
            nc = 1;
            for (Iterator i = columns.iterator(); i.hasNext(); nc++) {
            column c = (column) i.next();

            System.out.format("%3d %10s ", nc, c.id);
            }
            System.out.println();

            for (int i = 0; i < nrfl; i++) {
            for (int j = 0; j < ncol; j++) {
            float fval = swap ? ByteSwap.swap(dis.readFloat()) : dis.readFloat();
            System.out.format("%14.3g ", fval);
            }
            System.out.println();
            }
             */
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
                sgnum = Integer.parseInt(strarray[4]);
                if (strarray[5].startsWith("'")) {
                    sgname = strarray[5].substring(1, strarray[5].length() - 1);
                } else {
                    sgname = strarray[5];
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
                col = new column();
                columns.add(col);
                col.label = strarray[1];
                col.type = strarray[2].charAt(0);
                col.id = Integer.parseInt(strarray[5]);
                break;
            case PROJECT:
                ndset = Integer.parseInt(strarray[1]);
                try {
                    dset = (dataset) datasets.get(ndset - 1);
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
                try {
                    dset = (dataset) datasets.get(ndset - 1);
                } catch (IndexOutOfBoundsException e) {
                    dset = new dataset();
                    datasets.add(dset);
                }
                dset.dataset = strarray[2];
                break;
            case DCELL:
                ndset = Integer.parseInt(strarray[1]);
                try {
                    dset = (dataset) datasets.get(ndset - 1);
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
                try {
                    dset = (dataset) datasets.get(ndset - 1);
                } catch (IndexOutOfBoundsException e) {
                    dset = new dataset();
                    datasets.add(dset);
                }
                dset.lambda = Double.parseDouble(strarray[2]);
                break;
            case BATCH:
                break;
            case END:
                parsing = false;
                break;
            default:
                break;
        }

        return parsing;
    }

    private void parse_columns() {

        int nc = 0;
        StringBuffer sb = new StringBuffer();
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
                    || label.equalsIgnoreCase("rfree")
                    || label.equalsIgnoreCase("rfreeflag"))
                    && c.type == 'I') {
                sb.append(String.format("Reading R Free column: \"%s\"\n", c.label));
                rfree = nc;
            } else if ((label.equalsIgnoreCase("f")
                    || label.equalsIgnoreCase("fp")
                    || label.equalsIgnoreCase("fo"))
                    && c.type == 'F') {
                sb.append(String.format("Reading Fo column: \"%s\"\n", c.label));
                fo = nc;
            } else if ((label.equalsIgnoreCase("sigf")
                    || label.equalsIgnoreCase("sigfp")
                    || label.equalsIgnoreCase("sigfo"))
                    && c.type == 'Q') {
                sb.append(String.format("Reading sigFo column: \"%s\"\n", c.label));
                sigfo = nc;
            }
        }
        if (logger.isLoggable(Level.INFO)) {
            logger.info(sb.toString());
        }
    }

    static void print_header(MTZFilter mfile) {
        System.out.println("title: " + mfile.title);
        System.out.println("sg: " + mfile.sgname + " sgnum: " + mfile.sgnum);
        System.out.println("res: " + mfile.reslow + " " + mfile.reshigh);
        System.out.println("nrfl: " + mfile.nrfl);
        System.out.println("ncol: " + mfile.ncol);

        int ndset = 1;
        for (Iterator i = mfile.datasets.iterator(); i.hasNext(); ndset++) {
            dataset d = (dataset) i.next();

            System.out.println("  dset " + ndset + ": " + d.dataset);
            System.out.println("  project " + ndset + ": " + d.project);
            System.out.println("  wavelength " + ndset + ": " + d.lambda);
            System.out.println("  cell " + ndset + ": "
                    + d.cell[0] + " " + d.cell[1] + " " + d.cell[2] + " "
                    + d.cell[3] + " " + d.cell[4] + " " + d.cell[5]);
            System.out.println();
        }

    }
}
