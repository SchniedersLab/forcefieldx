//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.crystal;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import java.util.logging.Logger;

/**
 * This class holds the functionality to convert between equivalent space groups.
 */
public class SpaceGroupConversions {

  /**
   * The logger.
   */
  private static final Logger logger = Logger.getLogger(SpaceGroupConversions.class.getName());

  /**
   * Convert between hexagonal and rhombohedral space groups.
   *
   * @param crystal Crystal parameters to be converted.
   * @return Converted crystal.
   */
  public static Crystal hrConversion(Crystal crystal) {
    return hrConversion(crystal.a, crystal.b, crystal.c, crystal.alpha, crystal.beta, crystal.gamma,
        crystal.spaceGroup);
  }

  /**
   * Convert between hexagonal and rhombohedral space groups.
   *
   * @param a proposed axis length
   * @param b proposed axis length
   * @param c proposed axis length
   * @param alpha proposed angle
   * @param beta proposed angle
   * @param gamma proposed angle
   * @param currentSG Space group to be converted
   * @return Converted crystal satisfying other lattice system.
   */
  public static Crystal hrConversion(double a, double b, double c, double alpha, double beta,
      double gamma, SpaceGroup currentSG) {
    //Name for converted space group.
    String xtalName;
    // Going from hexagonal to rhombohedral (true), or visa versa (false).
    boolean hexStart = false;
    // Determine starting space group.
    switch (currentSG.shortName) {
      case ("H3") -> {
        logger.info(" Converting from H3 to R3:");
        xtalName = "R3";
        hexStart = true;
      }
      case ("H-3") -> {
        logger.info(" Converting from H-3 to R-3:");
        xtalName = "R-3";
        hexStart = true;
      }
      case ("H32") -> {
        logger.info(" Converting from H32 to R32:");
        xtalName = "R32";
        hexStart = true;
      }
      case ("H3m") -> {
        logger.info(" Converting from H3m to R3m:");
        xtalName = "R3m";
        hexStart = true;
      }
      case ("H3c") -> {
        logger.info(" Converting from H3c to R3c:");
        xtalName = "R3c";
        hexStart = true;
      }
      case ("H-3m") -> {
        logger.info(" Converting from H-3m to R-3m:");
        xtalName = "R-3m";
        hexStart = true;
      }
      case ("H-3c") -> {
        logger.info(" Converting from H-3c to R-3c:");
        xtalName = "R-3c";
        hexStart = true;
      }
      case ("R3") -> {
        logger.info(" Converting from R3 to H3:");
        xtalName = "H3";
      }
      case ("R-3") -> {
        logger.info(" Converting from R-3 to H-3:");
        xtalName = "H-3";
      }
      case ("R32") -> {
        logger.info(" Converting from R32 to H32:");
        xtalName = "H32";
      }
      case ("R3m") -> {
        logger.info(" Converting from R3m to H3m:");
        xtalName = "H3m";
      }
      case ("R3c") -> {
        logger.info(" Converting from R3c to H3c:");
        xtalName = "H3c";
      }
      case ("R-3m") -> {
        logger.info(" Converting from R-3m to H-3m:");
        xtalName = "H-3m";
      }
      case ("R-3c") -> {
        logger.info(" Converting from R-3c to H-3c:");
        xtalName = "H-3c";
      }
      default -> {
        logger.severe(
            format(" Unable to determine converted version for space group: %s", currentSG));
        return new Crystal(a, b, c, alpha, beta, gamma, currentSG.shortName);
      }
    }

    // Hexagonal and Rhombohedral space groups are frequently treated synonymously, therefore check if mislabeled.
    if (hexStart) {
      //Hexagonal aH = bH, alpha = beta = 90 gamma = 120
      if (LatticeSystem.RHOMBOHEDRAL_LATTICE.validParameters(a, b, c, alpha, beta, gamma)) {
        logger.info(" Crystal already has valid lattice parameters for new space group " + xtalName);
        return new Crystal(a, b, c, alpha, beta, gamma, xtalName);
      }

      double aR = sqrt(1.0 / 9.0 * (pow(c, 2) + 3 * pow(a, 2)));
      double aRAlpha = acos((2 * pow(c, 2) - 3 * pow(a, 2)) /
          (2 * pow(c, 2) + 6 * pow(a, 2))) / PI * 180;

      return new Crystal(aR, aR, aR, aRAlpha, aRAlpha, aRAlpha, xtalName);
    } else {
      if (LatticeSystem.HEXAGONAL_LATTICE.validParameters(a, b, c, alpha, beta, gamma)) {
        logger.info(" Crystal already has valid lattice parameters for new space group " + xtalName);
        return new Crystal(a, b, c, alpha, beta, gamma, xtalName);
      }

      double aH = 2 * pow(a, 2) * (1 - cos(alpha / 180 * PI));
      double cH = sqrt(3 * pow(a, 2) * (1 + 2 * cos(alpha / 180 * PI)));

      return new Crystal(aH, aH, cH, 90.00, 90.00, 120.00, xtalName);
    }
  }
}
