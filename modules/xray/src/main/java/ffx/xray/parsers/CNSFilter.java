// ******************************************************************************
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
// ******************************************************************************
package ffx.xray.parsers;

import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.min;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SpaceGroup;
import ffx.xray.DiffractionRefinementData;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * CNSFilter class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class CNSFilter implements DiffractionFileFilter {

  private static final Logger logger = Logger.getLogger(CNSFilter.class.getName());

  private final double[] cell = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  private double resHigh = -1.0;
  private String spaceGroupName = null;
  private int spaceGroupNum = -1;

  /** Constructor for CNSFilter. */
  public CNSFilter() {}

  /** {@inheritDoc} */
  @Override
  public ReflectionList getReflectionList(File cnsFile) {
    return getReflectionList(cnsFile, null);
  }

  /** {@inheritDoc} */
  @Override
  public ReflectionList getReflectionList(File cnsFile, CompositeConfiguration properties) {
    try {
      BufferedReader br = new BufferedReader(new FileReader(cnsFile));
      String string;
      while ((string = br.readLine()) != null) {
        String[] strArray = string.split("\\s+");
        if (strArray[0].equalsIgnoreCase("{")) {
          if (strArray[1].toLowerCase().startsWith("sg=")) {
            spaceGroupName = strArray[1].substring(3);
            cell[0] = parseDouble(strArray[2].substring(2));
            cell[1] = parseDouble(strArray[3].substring(2));
            cell[2] = parseDouble(strArray[4].substring(2));
            cell[3] = parseDouble(strArray[5].substring(6));
            cell[4] = parseDouble(strArray[6].substring(5));
            cell[5] = parseDouble(strArray[7].substring(6));
          }
        } else if (strArray[0].equalsIgnoreCase("CRYST1")) {
          cell[0] = parseDouble(strArray[1]);
          cell[1] = parseDouble(strArray[2]);
          cell[2] = parseDouble(strArray[3]);
          cell[3] = parseDouble(strArray[4]);
          cell[4] = parseDouble(strArray[5]);
          cell[5] = parseDouble(strArray[6]);
          spaceGroupName = SpaceGroup.pdb2ShortName(string.substring(55, 65));
        } else if (strArray[0].toLowerCase().startsWith("inde")) {
          break;
        }
      }
      br.close();
    } catch (IOException e) {
      String message = " CNS IO Exception.";
      logger.log(Level.WARNING, message, e);
      return null;
    }

    Resolution resolution = null;
    if (properties != null) {
      resolution = Resolution.checkProperties(properties);
      resHigh = resolution.resolution;
    }

    if (spaceGroupName != null) {
      spaceGroupNum = SpaceGroup.spaceGroupNumber(spaceGroupName);
    }

    if (spaceGroupNum < 0 || cell[0] < 0 || resolution == null) {
      logger.info(
          " The CNS file contains insufficient information to generate the reflection list.");
      return null;
    }

    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder();
      sb.append(format("\n Opening %s\n", cnsFile.getName()));
      sb.append(" Setting up Reflection List based on CNS:\n");
      sb.append(
          format(
              "  Spacegroup #: %d (name: %s)\n",
              spaceGroupNum, SpaceGroup.spaceGroupNames[spaceGroupNum - 1]));
      sb.append(format("  Resolution:   %8.3f\n", resHigh));
      sb.append(
          format(
              "  Cell:         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
              cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]));
      logger.info(sb.toString());
    }

    Crystal crystal =
        new Crystal(
            cell[0],
            cell[1],
            cell[2],
            cell[3],
            cell[4],
            cell[5],
            SpaceGroup.spaceGroupNames[spaceGroupNum - 1]);

    return new ReflectionList(crystal, resolution, properties);
  }

  /** {@inheritDoc} */
  @Override
  public double getResolution(File cnsFile, Crystal crystal) {
    double res = Double.POSITIVE_INFINITY;

    try {
      BufferedReader br = new BufferedReader(new FileReader(cnsFile));
      HKL hkl = new HKL();
      String string;
      while ((string = br.readLine()) != null) {
        String[] strArray = string.split("\\s+");
        for (int i = 0; i < strArray.length; i++) {
          if (strArray[i].toLowerCase().startsWith("inde")) {
            if (i < strArray.length - 3) {
              int ih = parseInt(strArray[i + 1]);
              int ik = parseInt(strArray[i + 2]);
              int il = parseInt(strArray[i + 3]);
              hkl.h(ih);
              hkl.k(ik);
              hkl.l(il);
              res = min(res, Crystal.res(crystal, hkl));
            }
          }
        }
      }
    } catch (IOException e) {
      String message = " CNS IO Exception.";
      logger.log(Level.WARNING, message, e);
      return -1.0;
    }

    return res;
  }

  /** {@inheritDoc} */
  @Override
  public boolean readFile(
      File cnsFile,
      ReflectionList reflectionList,
      DiffractionRefinementData refinementData,
      CompositeConfiguration properties) {

    int nRead, nRes, nIgnore, nFriedel, nCut;
    boolean transpose = false;

    StringBuilder sb = new StringBuilder();
    sb.append(format("\n Opening %s\n", cnsFile.getName()));

    if (refinementData.rFreeFlag < 0) {
      refinementData.setFreeRFlag(1);
      sb.append(format(" Setting R free flag to CNS default: %d\n", refinementData.rFreeFlag));
    }

    try {
      BufferedReader br = new BufferedReader(new FileReader(cnsFile));

      boolean hasHKL, hasFo, hasSigFo, hasFree;
      int ih, ik, il, free;
      double fo, sigFo;

      hasHKL = hasFo = hasSigFo = hasFree = false;
      ih = ik = il = free = -1;
      fo = sigFo = -1.0;

      // Check if HKLs need to be transposed or not.
      HKL mate = new HKL();
      int nPosIgnore = 0;
      int nTransIgnore = 0;

      String string;
      while ((string = br.readLine()) != null) {
        String[] strArray = string.split("\\s+");

        for (int i = 0; i < strArray.length; i++) {
          if (strArray[i].toLowerCase().startsWith("inde")) {
            if (i < strArray.length - 3) {
              ih = Integer.parseInt(strArray[i + 1]);
              ik = Integer.parseInt(strArray[i + 2]);
              il = Integer.parseInt(strArray[i + 3]);
              reflectionList.findSymHKL(ih, ik, il, mate, false);
              HKL hklpos = reflectionList.getHKL(mate);
              if (hklpos == null) {
                nPosIgnore++;
              }

              reflectionList.findSymHKL(ih, ik, il, mate, true);
              HKL hkltrans = reflectionList.getHKL(mate);
              if (hkltrans == null) {
                nTransIgnore++;
              }
            }
          }
        }
      }
      if (nPosIgnore > nTransIgnore) {
        transpose = true;
      }

      // column identifiers
      String foString = null;
      String sigFoString = null;
      String rFreeString = null;
      if (properties != null) {
        foString = properties.getString("fostring", null);
        sigFoString = properties.getString("sigfostring", null);
        rFreeString = properties.getString("rfreestring", null);
      }

      // reopen to start at beginning
      br = new BufferedReader(new FileReader(cnsFile));

      // read in data
      double[][] anofSigF = new double[refinementData.n][4];
      for (int i = 0; i < refinementData.n; i++) {
        anofSigF[i][0] = anofSigF[i][1] = anofSigF[i][2] = anofSigF[i][3] = Double.NaN;
      }
      nRead = nRes = nIgnore = nFriedel = nCut = 0;

      while ((string = br.readLine()) != null) {
        String[] strArray = string.split("\\s+");
        for (int i = 0; i < strArray.length; i++) {
          if (strArray[i].toLowerCase().startsWith("inde")) {
            if (hasHKL && hasFo && hasSigFo && hasFree) {
              boolean friedel = reflectionList.findSymHKL(ih, ik, il, mate, transpose);
              HKL hkl = reflectionList.getHKL(mate);
              if (hkl != null) {
                if (refinementData.fSigFCutoff > 0.0 && (fo / sigFo) < refinementData.fSigFCutoff) {
                  nCut++;
                } else if (friedel) {
                  anofSigF[hkl.index()][2] = fo;
                  anofSigF[hkl.index()][3] = sigFo;
                  nFriedel++;
                } else {
                  anofSigF[hkl.index()][0] = fo;
                  anofSigF[hkl.index()][1] = sigFo;
                }
                refinementData.setFreeR(hkl.index(), free);
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
            hasHKL = false;
            hasFo = false;
            hasSigFo = false;
            hasFree = false;
            if (i < strArray.length - 3) {
              ih = parseInt(strArray[i + 1]);
              ik = parseInt(strArray[i + 2]);
              il = parseInt(strArray[i + 3]);
              hasHKL = true;
            }
          }
          if (strArray[i].toLowerCase().startsWith("fobs=")
              || strArray[i].equalsIgnoreCase(foString + "=")) {
            fo = parseDouble(strArray[i + 1]);
            hasFo = true;
          }
          if (strArray[i].toLowerCase().startsWith("sigma=")
              || strArray[i].equalsIgnoreCase(sigFoString + "=")) {
            sigFo = parseDouble(strArray[i + 1]);
            hasSigFo = true;
          }
          if (strArray[i].toLowerCase().startsWith("test=")
              || strArray[i].equalsIgnoreCase(rFreeString + "=")) {
            free = parseInt(strArray[i + 1]);
            hasFree = true;
          }
        }
      }
      br.close();

      // Set up fsigf from F+ and F-.
      refinementData.generateFsigFfromAnomalousFsigF(anofSigF);
    } catch (IOException e) {
      String message = "CNS IO Exception.";
      logger.log(Level.WARNING, message, e);
      return false;
    }

    if (logger.isLoggable(Level.INFO)) {
      sb.append(format(" HKL read in:                             %d\n", nRead));
      sb.append(format(" HKL read as friedel mates:               %d\n", nFriedel));
      sb.append(format(" HKL not read in (too high resolution):   %d\n", nRes));
      sb.append(format(" HKL not read in (not in internal list?): %d\n", nIgnore));
      sb.append(format(" HKL not read in (F/sigF cutoff):         %d\n", nCut));
      sb.append(
          format(" HKL in internal list:                    %d\n", reflectionList.hkllist.size()));
      logger.info(sb.toString());
    }

    return true;
  }
}
