// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.parameters;

import static ffx.potential.parameters.ForceField.ForceFieldType.CHARGE;
import static ffx.utilities.KeywordGroup.PotentialFunctionParameter;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import ffx.utilities.FFXKeyword;
import java.util.Comparator;

/**
 * The ChargeType class defines a partial atomic charge type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXKeyword(name = "charge", clazz = String.class, keywordGroup = PotentialFunctionParameter, description =
    "[1 integer and 1 real] "
        + "Provides a value for a single atomic partial charge electrostatic parameter. "
        + "The integer modifier, if positive, gives the atom type number for which the charge parameter is to be defined. "
        + "Note that charge parameters are given for atom types, not atom classes. "
        + "If the integer modifier is negative, then the parameter value to follow applies only to the individual atom whose atom number is the negative of the modifier. "
        + "The real number modifier gives the values of the atomic partial charge in electrons.")
public final class ChargeType extends BaseType implements Comparator<String> {

  /** Partial atomic charge in units of electrons. */
  public final double charge;
  /** The atom type that uses this charge parameter. */
  public int atomType;

  /**
   * ChargeType constructor.
   *
   * @param atomType int
   * @param charge double
   */
  public ChargeType(int atomType, double charge) {
    super(CHARGE, "" + atomType);
    this.atomType = atomType;
    this.charge = charge;
  }

  /**
   * Average two ChargeType instances. The atom type that defines the new type must be supplied.
   *
   * @param chargeType1 a {@link ffx.potential.parameters.ChargeType} object.
   * @param chargeType2 a {@link ffx.potential.parameters.ChargeType} object.
   * @param atomType The atom type that defines the new type.
   * @return a {@link ffx.potential.parameters.ChargeType} object.
   */
  public static ChargeType average(ChargeType chargeType1, ChargeType chargeType2, int atomType) {
    if (chargeType1 == null || chargeType2 == null) {
      return null;
    }
    double charge = (chargeType1.charge + chargeType2.charge) / 2.0;
    return new ChargeType(atomType, charge);
  }

  /** {@inheritDoc} */
  @Override
  public int compare(String s1, String s2) {
    int t1 = parseInt(s1);
    int t2 = parseInt(s2);
    return Integer.compare(t1, t2);
  }

  /**
   * incrementType
   *
   * @param increment The amount to increment the atom type by.
   */
  public void incrementType(int increment) {
    this.atomType += increment;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted Charge type.
   */
  @Override
  public String toString() {
    return format("charge  %5d  % 7.5f", atomType, charge);
  }
}
