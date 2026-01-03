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
package ffx.xray.parallel;

import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * Enum representing the different methods of grid processing.
 */
public enum GridMethod {
  /**
   * Decompose into 3D special domains..
   */
  SPATIAL,
  /**
   * Decompose into 2D slices.
   */
  SLICE,
  /**
   * Decompose into 1D rows.
   */
  ROW;

  // Private logger for the parse method.
  private static final Logger logger = Logger.getLogger(GridMethod.class.getName());

  /**
   * Parses the provided method name and returns the corresponding GridMethod enum value.
   * If the provided method name is not recognized, the default value SLICE is returned.
   *
   * @param methodName the name of the grid method to parse (case-insensitive, trimming whitespace).
   * @return the GridMethod corresponding to the provided method name.
   */
  public static GridMethod parse(String methodName) {
    try {
      return valueOf(methodName.trim().toUpperCase());
    } catch (Exception e) {
      logger.info(format(" %s was not recognized; SLICE grid method selected.", methodName));
      return SLICE;
    }
  }
}
