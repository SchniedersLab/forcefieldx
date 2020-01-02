//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.realspace.parsers;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.realspace.RealSpaceRefinementData;

/**
 * <p>
 * CCP4MapFilter class.</p>
 *
 * @author Timothy D. Fenn
 * @see <a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4
 * map format</a>
 * @see <a href="http://www.ccp4.ac.uk/dist/html/library.html"
 * target="_blank">CCP4 library documentation</a>
 * @since 1.0
 */
public class CCP4MapFilter implements RealSpaceFileFilter {

    private static final Logger logger = Logger.getLogger(CCP4MapFilter.class.getName());

    /**
     * {@inheritDoc}
     */
    @Override
    public Crystal getCrystal(String fileName, CompositeConfiguration properties) {
        int imapData;
        int spaceGroup = -1;
        double cellA = -1.0;
        double cellB = -1.0;
        double cellC = -1.0;
        double cellAlpha = -1.0;
        double cellBeta = -1.0;
        double cellGamma = -1.0;

        ByteOrder byteOrder = ByteOrder.nativeOrder();
        FileInputStream fileInputStream;
        DataInputStream dataInputStream;

        // First determine byte order of file versus system.
        try {
            fileInputStream = new FileInputStream(fileName);
            dataInputStream = new DataInputStream(fileInputStream);

            dataInputStream.skipBytes(212);
            byte[] bytes = new byte[4];
            dataInputStream.read(bytes, 0, 4);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);

            imapData = byteBuffer.order(ByteOrder.BIG_ENDIAN).getInt();
            String stampString = Integer.toHexString(imapData);

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
            fileInputStream.close();
        } catch (Exception e) {
            String message = " Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
        }

        try {
            fileInputStream = new FileInputStream(fileName);
            dataInputStream = new DataInputStream(fileInputStream);

            dataInputStream.skipBytes(40);
            byte[] bytes = new byte[80];
            dataInputStream.read(bytes, 0, 80);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);

            cellA = byteBuffer.order(byteOrder).getFloat();
            cellB = byteBuffer.order(byteOrder).getFloat();
            cellC = byteBuffer.order(byteOrder).getFloat();
            cellAlpha = byteBuffer.order(byteOrder).getFloat();
            cellBeta = byteBuffer.order(byteOrder).getFloat();
            cellGamma = byteBuffer.order(byteOrder).getFloat();

            for (int i = 0; i < 3; i++) {
                byteBuffer.order(byteOrder).getInt();
            }
            for (int i = 0; i < 3; i++) {
                byteBuffer.order(byteOrder).getFloat();
            }

            spaceGroup = byteBuffer.order(byteOrder).getInt();
            fileInputStream.close();
        } catch (Exception e) {
            String message = " Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
        }

        return new Crystal(cellA, cellB, cellC,
                cellAlpha, cellBeta, cellGamma,
                SpaceGroup.spaceGroupNames[spaceGroup - 1]);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readFile(String filename, RealSpaceRefinementData refinementdata,
                            CompositeConfiguration properties) {

        int imapData;
        double cellA, cellB, cellC, cellAlpha, cellBeta, cellGamma;
        String stampString;

        ByteOrder byteOrder = ByteOrder.nativeOrder();
        FileInputStream fileInputStream;
        DataInputStream dataInputStream;

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double mean = 0.0;
        double sd = 0.0;
        double rmsd = 0.0;

        // First determine byte order of file versus system
        try {
            fileInputStream = new FileInputStream(filename);
            dataInputStream = new DataInputStream(fileInputStream);

            dataInputStream.skipBytes(212);
            byte[] bytes = new byte[4];
            dataInputStream.read(bytes, 0, 4);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);

            imapData = byteBuffer.order(ByteOrder.BIG_ENDIAN).getInt();
            stampString = Integer.toHexString(imapData);
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
            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
                sb.append(format(" Opening CCP4 map: %s\n", filename));
                logger.info(sb.toString());
            }

            fileInputStream.close();
        } catch (Exception e) {
            String message = " Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
        }

        try {
            fileInputStream = new FileInputStream(filename);
            dataInputStream = new DataInputStream(fileInputStream);

            byte[] bytes = new byte[2048];

            dataInputStream.read(bytes, 0, 1024);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);

            int[] ext = new int[3];
            ext[0] = byteBuffer.order(byteOrder).getInt();
            ext[1] = byteBuffer.order(byteOrder).getInt();
            ext[2] = byteBuffer.order(byteOrder).getInt();

            // mode (2 = reals, only one we accept)
            int mode = byteBuffer.order(byteOrder).getInt();

            int[] ori = new int[3];
            ori[0] = byteBuffer.order(byteOrder).getInt();
            ori[1] = byteBuffer.order(byteOrder).getInt();
            ori[2] = byteBuffer.order(byteOrder).getInt();

            int[] ni = new int[3];
            ni[0] = byteBuffer.order(byteOrder).getInt();
            ni[1] = byteBuffer.order(byteOrder).getInt();
            ni[2] = byteBuffer.order(byteOrder).getInt();

            cellA = byteBuffer.order(byteOrder).getFloat();
            cellB = byteBuffer.order(byteOrder).getFloat();
            cellC = byteBuffer.order(byteOrder).getFloat();
            cellAlpha = byteBuffer.order(byteOrder).getFloat();
            cellBeta = byteBuffer.order(byteOrder).getFloat();
            cellGamma = byteBuffer.order(byteOrder).getFloat();

            int[] axisi = new int[3];
            for (int i = 0; i < 3; i++) {
                int axis = byteBuffer.order(byteOrder).getInt();
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

            min = byteBuffer.order(byteOrder).getFloat();
            max = byteBuffer.order(byteOrder).getFloat();
            mean = byteBuffer.order(byteOrder).getFloat();
            int sg = byteBuffer.order(byteOrder).getInt();
            int nsymb = byteBuffer.order(byteOrder).getInt();
            int skew = byteBuffer.order(byteOrder).getInt();

            for (int i = 0; i < 12; i++) {
                byteBuffer.order(byteOrder).getFloat();
            }

            for (int i = 0; i < 15; i++) {
                byteBuffer.order(byteOrder).getInt();
            }

            byte[] word = new byte[2048];
            byteBuffer.order(byteOrder).get(word, 0, 4);
            String mapString = new String(word);
            sd = byteBuffer.order(byteOrder).getFloat();
            rmsd = byteBuffer.order(byteOrder).getFloat();

            if (logger.isLoggable(Level.INFO)) {
                StringBuilder sb = new StringBuilder();
                sb.append(format("  Column origin:  %d\t Extent: %d\n",
                        ori[0], ext[0]));
                sb.append(format("  Row origin:     %d\t Extent: %d\n",
                        ori[1], ext[1]));
                sb.append(format("  Section origin: %d\t Extent: %d\n",
                        ori[2], ext[2]));
                sb.append(format("  Axis order:     %d %d %d\n",
                        axisi[0], axisi[1], axisi[2]));
                sb.append(format("  Number of X, Y, Z columns: %d %d %d\n",
                        ni[0], ni[1], ni[2]));
                sb.append(format("  Spacegroup:     %d (%s)\n",
                        sg, SpaceGroup.spaceGroupNames[sg - 1]));
                sb.append(format("  Cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                        cellA, cellB, cellC, cellAlpha, cellBeta, cellGamma));
                logger.info(sb.toString());
            }

            int nlabel = byteBuffer.order(byteOrder).getInt();
            for (int i = 0; i < 10; i++) {
                byteBuffer.order(byteOrder).get(word, 0, 80);
                mapString = new String(word);
            }

            if (nsymb > 0) {
                byteBuffer.rewind();
                dataInputStream.read(bytes, 0, nsymb);
                for (int i = 0; i < nsymb / 80; i += 80) {
                    byteBuffer.order(byteOrder).get(word, 0, 80);
                    mapString = new String(word);
                }
            }

            byteBuffer.rewind();
            dataInputStream.read(bytes, 0, 2048);
            refinementdata.setData(new double[ext[0] * ext[1] * ext[2]]);
            int[] ijk = new int[3];
            int index, x, y, z;
            refinementdata.setOrigin(ori[axisi[0]], ori[axisi[1]], ori[axisi[2]]);
            int nx = ext[axisi[0]];
            int ny = ext[axisi[1]];
            int nz = ext[axisi[2]];
            refinementdata.setExtent(nx, ny, nz);
            refinementdata.setNI(ni[0], ni[1], ni[2]);
            for (ijk[2] = 0; ijk[2] < ext[2]; ijk[2]++) {
                for (ijk[1] = 0; ijk[1] < ext[1]; ijk[1]++) {
                    for (ijk[0] = 0; ijk[0] < ext[0]; ijk[0]++) {
                        x = ijk[axisi[0]];
                        y = ijk[axisi[1]];
                        z = ijk[axisi[2]];
                        index = x + nx * (y + ny * z);
                        refinementdata.getData()[index] = byteBuffer.order(byteOrder).getFloat();
                        if (!byteBuffer.hasRemaining()) {
                            byteBuffer.rewind();
                            dataInputStream.read(bytes, 0, 2048);
                        }
                    }
                }
            }
            fileInputStream.close();
        } catch (Exception e) {
            String message = " Fatal exception reading CCP4 map.\n";
            logger.log(Level.SEVERE, message, e);
        }

        return true;
    }
}
