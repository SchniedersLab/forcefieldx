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
package ffx.potential.parameters;

import static ffx.potential.parameters.ForceField.ForceFieldType.SOLUTE;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import java.util.Comparator;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The SoluteType class defines one implicit solvent radius.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class SoluteType extends BaseType implements Comparator<String> {

  /** A Logger for the SoluteType class. */
  private static final Logger logger = Logger.getLogger(SoluteType.class.getName());
  /** Solute atomic diameter for PB. */
  public double pbDiameter;
  /** Solute atomic diameter for ddCOSMO. */
  public double cosDiameter;
  /** Solute atomic diameter for GK. */
  public double gkDiameter;

  /** Optional SMARTS description. */
  public String description;
  /** Atom type for this solute type. */
  private int atomType;

  /**
   * Constructor for SoluteType.
   *
   * @param atomType Atom type.
   * @param pbDiameter Diameter for PB continuum electrostatics.
   * @param cosDiameter Diameter for ddCOSMO continuum electrostatics.
   * @param gkDiameter Diameter for GK continuum electrostatics.
   */
  public SoluteType(int atomType, double pbDiameter, double cosDiameter, double gkDiameter) {
    super(SOLUTE, Integer.toString(atomType));
    this.atomType = atomType;
    this.pbDiameter = pbDiameter;
    this.cosDiameter = cosDiameter;
    this.gkDiameter = gkDiameter;
    this.description = null;
  }

  /**
   * Constructor for SoluteType.
   *
   * @param atomType Atom type.
   * @param description Smarts description.
   * @param pbDiameter Diameter for PB continuum electrostatics.
   * @param cosDiameter Diameter for ddCOSMO continuum electrostatics.
   * @param gkDiameter Diameter for GK continuum electrostatics.
   */
  public SoluteType(int atomType, String description, double pbDiameter, double cosDiameter, double gkDiameter) {
    this(atomType, pbDiameter, cosDiameter, gkDiameter);
    this.description = description;
  }

  /**
   * Construct a SoluteType from an input string.
   *
   * @param input The overall input String.
   * @param tokens The input String tokenized.
   * @return a SoluteType instance.
   */
  public static SoluteType parse(String input, String[] tokens) {
    if (tokens.length < 4) {
      logger.log(Level.WARNING, "Invalid SOLUTE type:\n{0}", input);
    } else {
      try {
        if (tokens.length == 5) {
          int atomType = parseInt(tokens[1].trim());
          double pbDiameter = parseDouble(tokens[2].trim());
          double cosDiameter = parseDouble(tokens[3].trim());
          double gkDiameter = parseDouble(tokens[4].trim());
          return new SoluteType(atomType, pbDiameter, cosDiameter, gkDiameter);
        } else {
          int atomType = parseInt(tokens[1].trim());
          String description = tokens[2].trim();
          double pbDiameter = parseDouble(tokens[3].trim());
          double cosDiameter = parseDouble(tokens[4].trim());
          double gkDiameter = parseDouble(tokens[5].trim());
          return new SoluteType(atomType, description, pbDiameter, cosDiameter, gkDiameter);
        }
      } catch (NumberFormatException e) {
        String message = "Exception parsing SOLUTE type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public int compare(String key1, String key2) {
    int type1 = parseInt(key1);
    int type2 = parseInt(key2);
    return Integer.compare(type1, type2);
  }

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    SoluteType soluteType = (SoluteType) o;
    return soluteType.atomType == this.atomType;
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(atomType);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted bond stretch string.
   */
  @Override
  public String toString() {
    if (description == null) {
      // TODO: updated solute lines to be "PB ddCOSMO GK" radii (paper order)
      return format("solute  %4d  %7.4f  %7.4f", atomType, pbDiameter, cosDiameter, gkDiameter);
    } else {
      return format("solute  %4d  %30s  %7.4f  %7.4f", atomType, description, pbDiameter, cosDiameter, gkDiameter);
    }
  }

  /**
   * incrementType
   *
   * @param increment a int.
   */
  void incrementType(int increment) {
    atomType += increment;
    setKey(Integer.toString(atomType));
  }
}
