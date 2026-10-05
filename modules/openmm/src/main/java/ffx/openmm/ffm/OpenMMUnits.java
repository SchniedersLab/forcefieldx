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
package ffx.openmm.ffm;

/**
 * Standard unit conversion constants used by OpenMM.
 *
 * <p>OpenMM consistently uses the following standard internal units:
 * <ul>
 *   <li>Length: nanometers (nm)</li>
 *   <li>Time: picoseconds (ps)</li>
 *   <li>Mass: atomic mass units / daltons (Da)</li>
 *   <li>Charge: proton charge (e)</li>
 *   <li>Temperature: Kelvin (K)</li>
 *   <li>Angle: radians</li>
 *   <li>Energy: kilojoules per mole (kJ/mol)</li>
 *   <li>Force: kilojoules per mole per nanometer (kJ/(mol*nm))</li>
 * </ul>
 *
 * <p>This utility class provides constants to convert between OpenMM internal units and other
 * common scientific units such as Angstroms, femtoseconds, degrees, and kilocalories per mole.</p>
 */
public final class OpenMMUnits {

  /**
   * The number of nanometers in an Angstrom (0.1).
   */
  public static final double NM_PER_ANGSTROM = 0.1;

  /**
   * The number of Angstroms in a nanometer (10.0).
   */
  public static final double ANGSTROMS_PER_NM = 10.0;

  /**
   * The number of picoseconds in a femtosecond (0.001).
   */
  public static final double PS_PER_FS = 0.001;

  /**
   * The number of femtoseconds in a picosecond (1000.0).
   */
  public static final double FS_PER_PS = 1000.0;

  /**
   * The number of kilojoules in a kilocalorie (4.184).
   */
  public static final double KJ_PER_KCAL = 4.184;

  /**
   * The number of kilocalories in a kilojoule (1.0 / 4.184).
   */
  public static final double KCAL_PER_KJ = 1.0 / 4.184;

  /**
   * The number of radians in a degree (pi / 180.0).
   */
  public static final double RADIANS_PER_DEGREE = 3.1415926535897932385 / 180.0;

  /**
   * The number of degrees in a radian (180.0 / pi).
   */
  public static final double DEGREES_PER_RADIAN = 180.0 / 3.1415926535897932385;

  /**
   * Lennard-Jones sigma per unit van der Waals radius: 2.0 / (2^(1/6)).
   */
  public static final double SIGMA_PER_VDW_RADIUS = 1.7817974362806786095;

  /**
   * Van der Waals radius per unit Lennard-Jones sigma: (2^(1/6)) / 2.0.
   */
  public static final double VDW_RADIUS_PER_SIGMA = 0.56123102415468649070;

  private OpenMMUnits() {
  }
}
