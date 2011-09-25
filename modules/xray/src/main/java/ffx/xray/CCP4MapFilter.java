/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.xray;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

/**
 * <p>CCP4MapFilter class.</p>
 *
 * @author Tim Fenn
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4 map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">CCP4 library documentation</a>
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4 map format</a>
 *
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">CCP4 library documentation</a>
 * @version $Id: $
 */
public class CCP4MapFilter implements RealSpaceFileFilter {

    private static final Logger logger = Logger.getLogger(CCP4MapFilter.class.getName());

    /** {@inheritDoc} */
    @Override
    public Crystal getCrystal(String filename, CompositeConfiguration properties) {
        int imapdata;
        int sg = -1;
        double cella = -1.0;
        double cellb = -1.0;
        double cellc = -1.0;
        double cellalpha = -1.0;
        double cellbeta = -1.0;
        double cellgamma = -1.0;

        ByteOrder b = ByteOrder.nativeOrder();

        FileInputStream fis;
        DataInputStream dis;

        // first determine byte order of file versus system
        try {
            fis = new FileInputStream(filename);
            dis = new DataInputStream(fis);

            dis.skipBytes(212);
            byte bytes[] = new byte[4];
            dis.read(bytes, 0, 4);
            ByteBuffer bb = ByteBuffer.wrap(bytes);

            imapdata = bb.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampstr = Integer.toHexString(imapdata);
            // System.out.println("stamp: " + stampstr);
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
            fis.close();
        } catch (Exception e) {
            String message = "Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }

        try {
            fis = new FileInputStream(filename);
            dis = new DataInputStream(fis);

            dis.skipBytes(40);
            byte bytes[] = new byte[80];
            dis.read(bytes, 0, 80);
            ByteBuffer bb = ByteBuffer.wrap(bytes);

            cella = bb.order(b).getFloat();
            cellb = bb.order(b).getFloat();
            cellc = bb.order(b).getFloat();
            cellalpha = bb.order(b).getFloat();
            cellbeta = bb.order(b).getFloat();
            cellgamma = bb.order(b).getFloat();

            for (int i = 0; i < 3; i++) {
                bb.order(b).getInt();
            }
            for (int i = 0; i < 3; i++) {
                bb.order(b).getFloat();
            }

            sg = bb.order(b).getInt();
            fis.close();
        } catch (Exception e) {
            String message = "Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }

        return new Crystal(cella, cellb, cellc,
                cellalpha, cellbeta, cellgamma,
                SpaceGroup.spaceGroupNames[sg - 1]);
    }

    /** {@inheritDoc} */
    @Override
    public boolean readFile(String filename, RealSpaceRefinementData refinementdata,
            CompositeConfiguration properties) {

        int imapdata;
        double cella, cellb, cellc, cellalpha, cellbeta, cellgamma;
        String stampstr;

        ByteOrder b = ByteOrder.nativeOrder();

        FileInputStream fis;
        DataInputStream dis;

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double mean = 0.0;
        double sd = 0.0;
        double rmsd = 0.0;

        // first determine byte order of file versus system
        try {
            fis = new FileInputStream(filename);
            dis = new DataInputStream(fis);

            dis.skipBytes(212);
            byte bytes[] = new byte[4];
            dis.read(bytes, 0, 4);
            ByteBuffer bb = ByteBuffer.wrap(bytes);

            imapdata = bb.order(ByteOrder.BIG_ENDIAN).getInt();
            stampstr = Integer.toHexString(imapdata);
            // System.out.println("stamp: " + stampstr);
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

            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("\nOpening CCP4 map: %s\n", filename));
                sb.append(String.format("file type (machine stamp): %s\n",
                        stampstr));
                logger.info(sb.toString());
            }

            fis.close();
        } catch (Exception e) {
            String message = "Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }

        try {
            fis = new FileInputStream(filename);
            dis = new DataInputStream(fis);

            byte bytes[] = new byte[2048];

            dis.read(bytes, 0, 1024);
            ByteBuffer bb = ByteBuffer.wrap(bytes);

            int ext[] = new int[3];
            ext[0] = bb.order(b).getInt();
            ext[1] = bb.order(b).getInt();
            ext[2] = bb.order(b).getInt();

            // mode (2 = reals, only one we accept)
            int mode = bb.order(b).getInt();

            int ori[] = new int[3];
            ori[0] = bb.order(b).getInt();
            ori[1] = bb.order(b).getInt();
            ori[2] = bb.order(b).getInt();

            int ni[] = new int[3];
            ni[0] = bb.order(b).getInt();
            ni[1] = bb.order(b).getInt();
            ni[2] = bb.order(b).getInt();

            cella = bb.order(b).getFloat();
            cellb = bb.order(b).getFloat();
            cellc = bb.order(b).getFloat();
            cellalpha = bb.order(b).getFloat();
            cellbeta = bb.order(b).getFloat();
            cellgamma = bb.order(b).getFloat();

            int axisi[] = new int[3];
            for (int i = 0; i < 3; i++) {
                int axis = bb.order(b).getInt();
                switch (axis) {
                    case 1:
                        axisi[0] = i;
                        break;
                    case 2:
                        axisi[1] = i;
                        break;
                    case 3:
                        axisi[2] = i;
                        break;
                }
            }

            min = bb.order(b).getFloat();
            max = bb.order(b).getFloat();
            mean = bb.order(b).getFloat();

            int sg = bb.order(b).getInt();

            int nsymb = bb.order(b).getInt();

            int skew = bb.order(b).getInt();

            for (int i = 0; i < 12; i++) {
                bb.order(b).getFloat();
            }

            for (int i = 0; i < 15; i++) {
                bb.order(b).getInt();
            }

            byte word[] = new byte[2048];
            bb.order(b).get(word, 0, 4);
            String mapstr = new String(word);
            // System.out.println("MAP?: " + mapstr);

            sd = bb.order(b).getFloat();
            rmsd = bb.order(b).getFloat();

            /*
            System.out.println("col: " + ori[0] + " " + ext[0] + " " + ni[0]);
            System.out.println("row: " + ori[1] + " " + ext[1] + " " + ni[1]);
            System.out.println("sec: " + ori[2] + " " + ext[2] + " " + ni[2]);
            System.out.println("order: " + axisi[0] + " " + axisi[1] + " " + axisi[2]);
            System.out.println("min: " + min + " max: " + max + " mean: " + mean);
            System.out.println("sd: " + sd + " rmsd: " + rmsd);
            System.out.println("sg: " + sg);
            System.out.println("a: " + cella + " b: " + cellb + " c: " + cellc
            + " alpha: " + cellalpha + " beta: " + cellbeta + " gamma: " + cellgamma);
             */

            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("  column origin: %d extent: %d\n",
                        ori[0], ext[0]));
                sb.append(String.format("  row origin: %d extent: %d\n",
                        ori[1], ext[1]));
                sb.append(String.format("  section origin: %d extent: %d\n",
                        ori[2], ext[2]));
                sb.append(String.format("  axis order: %d %d %d\n",
                        axisi[0], axisi[1], axisi[2]));
                sb.append(String.format("  number of X, Y, Z columns: %d %d %d\n",
                        ni[0], ni[1], ni[2]));
                sb.append(String.format("  spacegroup #: %d (name: %s)\n",
                        sg, SpaceGroup.spaceGroupNames[sg - 1]));
                sb.append(String.format("  cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                        cella, cellb, cellc, cellalpha, cellbeta, cellgamma));
                logger.info(sb.toString());
            }

            int nlabel = bb.order(b).getInt();

            // System.out.println("nsymb: " + nsymb + " nlabel: " + nlabel);
            for (int i = 0; i < 10; i++) {
                bb.order(b).get(word, 0, 80);
                mapstr = new String(word);
                // System.out.println("label " + i + " : " + mapstr);
            }

            if (nsymb > 0) {
                bb.rewind();
                dis.read(bytes, 0, nsymb);
                for (int i = 0; i < nsymb / 80; i += 80) {
                    bb.order(b).get(word, 0, 80);
                    mapstr = new String(word);
                    // System.out.println("symm: " + mapstr);
                }
            }

            bb.rewind();
            dis.read(bytes, 0, 2048);
            refinementdata.data = new double[ext[0] * ext[1] * ext[2]];
            int ijk[] = new int[3];
            int index, x, y, z;
            refinementdata.ori[0] = ori[axisi[0]];
            refinementdata.ori[1] = ori[axisi[1]];
            refinementdata.ori[2] = ori[axisi[2]];
            int nx = ext[axisi[0]];
            int ny = ext[axisi[1]];
            int nz = ext[axisi[2]];
            refinementdata.ext[0] = nx;
            refinementdata.ext[1] = ny;
            refinementdata.ext[2] = nz;
            refinementdata.ni[0] = ni[0];
            refinementdata.ni[1] = ni[1];
            refinementdata.ni[2] = ni[2];
            for (ijk[2] = 0; ijk[2] < ext[2]; ijk[2]++) {
                for (ijk[1] = 0; ijk[1] < ext[1]; ijk[1]++) {
                    for (ijk[0] = 0; ijk[0] < ext[0]; ijk[0]++) {
                        x = ijk[axisi[0]];
                        y = ijk[axisi[1]];
                        z = ijk[axisi[2]];
                        index = x + nx * (y + ny * z);
                        refinementdata.data[index] = bb.order(b).getFloat();
                        if (!bb.hasRemaining()) {
                            bb.rewind();
                            dis.read(bytes, 0, 2048);
                        }
                    }
                }
            }
            fis.close();
        } catch (Exception e) {
            String message = "Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }

        return true;
    }
}
