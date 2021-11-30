// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.utilities;

/**
 * Library class containing constants such as Avogadro's number.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Constants {

  // SI units: kg, m, s, C, K, mol, lm
  // Our typical units: g/mol, Angstrom, psec, elementary charges (+1 proton charge), K, mol, N/A
  // Below constants are the seven defining constants of SI as of May 20 2019 (BIPM).

  /**
   * Hyperfine transition frequency of cesium in Hertz, defining the second. <code>
   * DEL_V_Cs_SI=9192631770L</code>
   */
  public static final long DEL_V_Cs_SI = 9192631770L;
  /** Speed of light in m/s, defining the meter. <code>SPEED_OF_LIGHT_SI=299792458</code> */
  public static final int SPEED_OF_LIGHT_SI = 299792458;
  /**
   * Planck constant in J*s, defining the kilogram (by defining the derived Joule) <code>
   * PLANCK_CONSTANT_SI=6.62607015E-34d</code>
   */
  public static final double PLANCK_CONSTANT_SI = 6.62607015E-34d;
  /**
   * Elementary charge in Coulombs, defining the Coulomb. <code>
   *  ELEMENTARY_CHARGE_SI=1.602176634E-19d</code>
   */
  public static final double ELEMENTARY_CHARGE_SI = 1.602176634E-19d;
  /** Boltzmann's constant in J/K, defining the Kelvin. <code>BOLTZMANN_SI=1.380649E-23d</code> */
  public static final double BOLTZMANN_SI = 1.380649E-23d;
  /** Avogadro's number, defining the mol. <code>AVOGADRO=6.02214076E23d</code> */
  public static final double AVOGADRO = 6.02214076E23d;
  /** Luminous efficacy in lm/W, defining the lumen. <code>K_CD_SI=683</code> */
  public static final int K_CD_SI = 683;
  /** Constant <code>LITERS_PER_CUBIC_ANGSTROM=1E-30</code> */
  public static final double LITERS_PER_CUBIC_ANGSTROM = 1E-30;
  /** Constant <code>ATM_TO_BAR=1.01325</code> */
  public static final double ATM_TO_BAR = 1.01325;
  /** Constant <code>KCAL_TO_KJ=4.184</code> */
  public static final double KCAL_TO_KJ = 4.184;
  /** Constant <code>KJ_TO_KCAL=1.0 / KCAL_TO_KJ</code> */
  public static final double KJ_TO_KCAL = 1.0 / KCAL_TO_KJ;
  /** Constant <code>METERS_TO_ANG=1E10</code> */
  public static final double METERS_TO_ANG = 1E10;
  /** Constant <code>SEC_TO_PSEC=1E12</code> */
  public static final double SEC_TO_PSEC = 1E12;
  /** Constant <code>KG_TO_GRAMS=1000</code> */
  public static final double KG_TO_GRAMS = 1000;
  /** Constant <code>PSEC_TO_FSEC=1000</code> */
  public static final double PSEC_TO_FSEC = 1000;
  /** Constant <code>FSEC_TO_PSEC=0.001</code> */
  public static final double FSEC_TO_PSEC = 0.001;
  /**
   * Ideal gas constant in kcal/(mol*K) <code>R = BOLTZMANN_SI * AVOGADRO * 0.001 * KJ_TO_KCAL
   * </code>
   */
  public static final double R = BOLTZMANN_SI * AVOGADRO * 0.001 * KJ_TO_KCAL;
  /**
   * Boltzmann/ideal gas constant in units of g*Ang^2/(mol*psec^2*K). <code> kB = BOLTZMANN_SI *
   * KG_TO_GRAMS * METERS_TO_ANG * METERS_TO_ANG * AVOGADRO / (SEC_TO_PSEC * SEC_TO_PSEC)</code>
   */
  public static final double kB =
      BOLTZMANN_SI
          * KG_TO_GRAMS
          * METERS_TO_ANG
          * METERS_TO_ANG
          * AVOGADRO
          / (SEC_TO_PSEC * SEC_TO_PSEC);
  /** Conversion from kcal/mol/Ang^3 to Atm. <code>PRESCON=6.85684112e4</code> */
  public static final double PRESCON = 6.85684112e4;
  /** Permittivity of water at STP. <code>dWater=78.3</code> */
  public static final double dWater = 78.3;
  /** Convert nanoseconds to seconds. <code>NS2SEC=1e-9</code> */
  public static final double NS2SEC = 1e-9;
  /** Room temperature ~= 298.15 Kelvins. <code>ROOM_TEMPERATURE=298.15</code> */
  public static final double ROOM_TEMPERATURE = 298.15;
  /**
   * Coulomb constant in units of kcal*Ang/(mol*electron^2), as derived from CODATA 2018
   * permittivity of free space measured at 8.8541878128*10^-12 F/m <code>
   * ELECTRIC_CODATA_2018=332.063713299</code>
   */
  public static final double ELECTRIC_CODATA_2018 = 332.063713299;
  /**
   * Coulomb constant in units of kcal*Ang/(mol*electron^2)
   *
   * <p>Note -- this value varies slightly between force field definitions and can be set using the
   * ELECTRIC property. As such, SHOULD NOT ever be updated to SI/CODATA standards, but rather kept
   * up-to-date with the coulomb parameter in Tinker/source/units.f. At present, the Tinker value
   * appears to be a truncated version of the Coulomb constant derived from CODATA 2018. <code>
   * DEFAULT_ELECTRIC=332.063713</code>
   */
  public static final double DEFAULT_ELECTRIC = 332.063713;
  /** Conversion from electron-Angstroms to Debyes. <code>ELEC_ANG_TO_DEBYE=4.80321</code> */
  public static final double ELEC_ANG_TO_DEBYE = 4.80321;
  /**
   * Conversion from electron-Angstroms^2 to Buckinghams. <code> ELEC_ANG2_TO_BUCKINGHAMS =
   * ELEC_ANG_TO_DEBYE * ELEC_ANG_TO_DEBYE</code>
   */
  public static final double ELEC_ANG2_TO_BUCKINGHAMS = ELEC_ANG_TO_DEBYE * ELEC_ANG_TO_DEBYE;
  /**
   * Conversion from kcal/mole to g*Ang**2/ps**2. <code>KCAL_TO_GRAM_ANG2_PER_PS2=4.1840e2</code>
   */
  public static final double KCAL_TO_GRAM_ANG2_PER_PS2 = 4.1840e2;
  /** Conversion from Bohr to Angstroms. <code>BOHR=0.52917720859</code> */
  public static final double BOHR = 0.52917720859;
  /** Conversion from Bohr^2 to Angstroms^2. <code>BOHR2 = BOHR * BOHR</code> */
  public static final double BOHR2 = BOHR * BOHR;

  // Library class: make the default constructor private to ensure it's never constructed.
  private Constants() {}
}
