// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
// ******************************************************************************
package ffx.xray.parsers;

import ffx.crystal.Crystal;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.TinkerUtils.version;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * CCP4MapWriter class.
 *
 * @author Timothy D. Fenn
 * @see <ul>
 *     <li><a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4 map format</a>
 *     <li><a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">CCP4 library
 *     documentation </a>
 *     <li><a href="http://www.ccp4.ac.uk/html/maplib.html" target="_blank">CCP4 map format </a>
 *     <li><a href="http://www.ccp4.ac.uk/dist/html/library.html" target="_blank">CCP4 library
 *     documentation </a>
 *     </ul>
 * @since 1.0
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
  public CCP4MapWriter(int extx, int exty, int extz, Crystal crystal, String filename) {
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
   * Constructor for CCP4MapWriter.
   *
   * @param orix an int.
   * @param oriy an int.
   * @param oriz an int.
   * @param extx an int.
   * @param exty an int.
   * @param extz an int.
   * @param nx an int.
   * @param ny an int.
   * @param nz an int.
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @param filename a {@link java.lang.String} object.
   */
  public CCP4MapWriter(
      int orix, int oriy, int oriz,
      int extx, int exty, int extz,
      int nx, int ny, int nz,
      Crystal crystal, String filename) {
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
   * Set the stepping across the array (e.g. 2 if data is separated by 1 space)
   *
   * @param stride the step size desired
   */
  public void setStride(int stride) {
    this.stride = stride;
  }

  /**
   * Write data to file (does not normalize).
   *
   * @param data map data to write out
   */
  public void write(double[] data) {
    write(data, false);
  }

  /**
   * write data to file, does not normalize
   *
   * @param data map data to write out
   * @param norm Normalize the data by mean/sd.
   */
  public void write(double[] data, boolean norm) {
    double min = Double.POSITIVE_INFINITY;
    double max = Double.NEGATIVE_INFINITY;
    double mean = 0.0;
    double sd = 0.0;

    int n = 0;
    for (int k = 0; k < extz; k++) {
      for (int j = 0; j < exty; j++) {
        for (int i = 0; i < extx; i++) {
          int index = stride * (i + extx * (j + exty * k));
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
          sd += pow(data[index] - mean, 2.0);
          n++;
        }
      }
    }
    sd = sqrt(sd / n);

    if (norm) {
      if (logger.isLoggable(Level.INFO)) {
        StringBuilder sb = new StringBuilder();
        sb.append("\n Normalizing CCP4 Map");
        sb.append(format("\n  Min:           %8.4f", min));
        sb.append(format("\n  Max:           %8.4f", max));
        sb.append(format("\n  Mean:          %8.4f", mean));
        sb.append(format("\n  Standard Dev.: %8.4f", sd));
        logger.info(sb.toString());
      }

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
      return;
    }

    File file = version(new File(filename));
    String name = file.getName();

    try {
      if (logger.isLoggable(Level.INFO)) {
        StringBuilder sb = new StringBuilder();
        sb.append(format("\n Writing CCP4 map file: \"%s\"", name));
        sb.append(format("\n  Min:           %8.4f", min));
        sb.append(format("\n  Max:           %8.4f", max));
        sb.append(format("\n  Mean:          %8.4f", mean));
        sb.append(format("\n  Standard Dev.: %8.4f", sd));
        logger.info(sb.toString());
      }

      FileOutputStream fos = new FileOutputStream(file);
      DataOutputStream dos = new DataOutputStream(fos);

      byte[] bytes = new byte[2048];
      ByteBuffer bb = ByteBuffer.wrap(bytes);
      // Set the byte order to the native order of the operating system.
      bb.order(ByteOrder.nativeOrder());

      // Header
      bb.putInt(extx).putInt(exty).putInt(extz);
      // mode (2 = reals, only one we accept)
      bb.putInt(2);
      bb.putInt(orix).putInt(oriy).putInt(oriz);
      bb.putInt(nx).putInt(ny).putInt(nz);

      bb.putFloat((float) crystal.a);
      bb.putFloat((float) crystal.b);
      bb.putFloat((float) crystal.c);
      bb.putFloat((float) crystal.alpha);
      bb.putFloat((float) crystal.beta);
      bb.putFloat((float) crystal.gamma);

      bb.putInt(1).putInt(2).putInt(3);
      bb.putFloat((float) min);
      bb.putFloat((float) max);
      bb.putFloat((float) mean);

      bb.putInt(crystal.spaceGroup.number);
      // bb.putInt(1);

      // symmetry bytes - should set this up at some point
      // imapdata = swap ? ByteSwap.swap(320) : 320;
      bb.putInt(80);
      bb.putInt(0);

      for (int i = 0; i < 12; i++) {
        bb.putFloat(0.0f);
      }

      for (int i = 0; i < 15; i++) {
        bb.putInt(0);
      }

      int offset = 0;
      dos.write(bytes, offset, 208);
      bb.rewind();

      String mapstr = "MAP ";
      dos.writeBytes(mapstr);

      // machine code: double, float, int, uchar
      // 0x4441 for LE, 0x1111 for BE
      int imapdata;
      if (ByteOrder.nativeOrder().equals(ByteOrder.LITTLE_ENDIAN)) {
        imapdata = 0x4441;
      } else {
        imapdata = 0x1111;
      }
      bb.putInt(imapdata);
      bb.putFloat((float) sd);
      bb.putInt(1);
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
            float fmapdata = (float) data[index];
            bb.putFloat(fmapdata);
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
    }
  }
}
