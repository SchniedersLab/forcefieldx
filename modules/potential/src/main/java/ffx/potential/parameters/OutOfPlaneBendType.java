// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import ffx.utilities.FFXProperty;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.ForceFieldType.OPBEND;
import static ffx.utilities.Constants.DEGREES_PER_RADIAN;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.PropertyGroup.EnergyUnitConversion;
import static ffx.utilities.PropertyGroup.LocalGeometryFunctionalForm;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.String.valueOf;
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The OutOfPlaneBendType class defines one Allinger style out-of-plane angle bending energy type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "opbend", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [4 integers and 1 real]
    Provides the values for a single out-of-plane bending potential parameter.
    The first integer modifier is the atom class of the out-of-plane atom and the second integer is the atom class of the central trigonal atom.
    The third and fourth integers give the atom classes of the two remaining atoms attached to the trigonal atom.
    Values of zero for the third and fourth integers are treated as wildcards, and can represent any atom type.
    The real number modifier gives the force constant value for the out-of-plane angle.
    The default units for the force constant are kcal/mole/radian^2, but this can be controlled via the opbendunit keyword.
    """)
public final class OutOfPlaneBendType extends BaseType implements Comparator<String> {

  /**
   * Default cubic coefficient in out-of-plane angle bending potential.
   */
  public static final double DEFAULT_OPBEND_CUBIC = 0.0;
  /**
   * Default quartic coefficient in out-of-plane angle bending potential.
   */
  public static final double DEFAULT_OPBEND_QUARTIC = 0.0;
  /**
   * Default pentic coefficient in out-of-plane angle bending potential.
   */
  public static final double DEFAULT_OPBEND_PENTIC = 0.0;
  /**
   * Default quintic coefficient in out-of-plane angle bending potential.
   */
  public static final double DEFAULT_OPBEND_SEXTIC = 0.0;

  /**
   * Cubic coefficient in out-of-plane angle bending potential.
   */
  @FFXProperty(name = "opbend-cubic", propertyGroup = LocalGeometryFunctionalForm, defaultValue = "0.0", description = """ 
      Sets the value of the cubic term in the Taylor series expansion form of the out-of-plane bending potential energy.
      The real number modifier gives the value of the coefficient as a multiple of the quadratic coefficient.
      This term multiplied by the out-of-plane bending energy unit conversion factor, the force constant,
      and the cube of the deviation of the out-of-plane angle from zero gives the cubic contribution to the out-of-plane bending energy.
      The default value in the absence of the opbend-cubic keyword is zero; i.e., the cubic out-of-plane bending term is omitted.
      """)
  public double cubic = DEFAULT_OPBEND_CUBIC;

  /**
   * Quartic coefficient in out-of-plane angle bending potential.
   */
  @FFXProperty(name = "opbend-quartic", propertyGroup = LocalGeometryFunctionalForm, defaultValue = "0.0", description = """
      Sets the value of the quartic term in the Taylor series expansion form of the out-of-plane bending potential energy.
      The real number modifier gives the value of the coefficient as a multiple of the quadratic coefficient.
      This term multiplied by the out-of-plane bending energy unit conversion factor, the force constant,
      and the forth power of the deviation of the out-of-plane angle from zero gives the quartic contribution to the out-of-plane bending energy.
      The default value in the absence of the opbend-quartic keyword is zero; i.e., the quartic out-of-plane bending term is omitted.
      """)
  public double quartic = DEFAULT_OPBEND_QUARTIC;

  /**
   * Quintic coefficient in out-of-plane angle bending potential.
   */
  @FFXProperty(name = "opbend-pentic", propertyGroup = LocalGeometryFunctionalForm, defaultValue = "0.0", description = """
      Sets the value of the fifth power term in the Taylor series expansion form of the out-of-plane bending potential energy.
      The real number modifier gives the value of the coefficient as a multiple of the quadratic coefficient.
      This term multiplied by the out-of-plane bending energy unit conversion factor, the force constant,
      and the fifth power of the deviation of the out-of-plane angle from zero gives the pentic contribution to the out-of-plane bending energy.
      The default value in the absence of the opbend-pentic keyword is zero; i.e., the pentic out-of-plane bending term is omitted.
      """)
  public double pentic = DEFAULT_OPBEND_PENTIC;

  /**
   * Sextic coefficient in out-of-plane angle bending potential.
   */
  @FFXProperty(name = "opbend-sextic", propertyGroup = LocalGeometryFunctionalForm, defaultValue = "0.0", description = """
      Sets the value of the sixth power term in the Taylor series expansion form of the out-of-plane bending potential energy.
      The real number modifier gives the value of the coefficient as a multiple of the quadratic coefficient.
      This term multiplied by the out-of-plane bending energy unit conversion factor, the force constant,
      and the sixth power of the deviation of the out-of-plane angle from zero gives the sextic contribution to the out-of-plane bending energy.
      The default value in the absence of the opbend-sextic keyword is zero; i.e., the sextic out-of-plane bending term is omitted.
      """)
  public double sextic = DEFAULT_OPBEND_SEXTIC;

  /**
   * Convert Out-of-Plane bending energy to kcal/mole.
   * <p>
   * TINKER v.5 and v.6 Units: 1.0 / (180.0/PI)^2 = 0.00030461741979
   * <p>
   * TINKER v.4 Units: 0.02191418
   * <p>
   * Ratio of v.4 to v.5/6 = 0.02191418 / 1.0 / (180.0/PI)^2 = 71.94
   */
  @FFXProperty(name = "opbendunit", propertyGroup = EnergyUnitConversion, defaultValue = "(Pi/180)^2", description = """
      Sets the scale factor needed to convert the energy value computed by the out-of-plane bending potential into units of kcal/mole. "
      The correct value is force field dependent and typically provided in the header of the master force field parameter file.
      """)
  public double opBendUnit = DEFAULT_OPBEND_UNIT;

  public static final double DEFAULT_OPBEND_UNIT = pow(PI / 180.0, 2);
  /**
   * A Logger for the OutOfPlaneBendType class.
   */
  private static final Logger logger = Logger.getLogger(OutOfPlaneBendType.class.getName());
  /**
   * Atom classes for this out-of-plane angle bending type.
   */
  public final int[] atomClasses;
  /**
   * Force constant (Kcal/mol/Angstrom).
   */
  public final double forceConstant;

  /**
   * OutOfPlaneBendType Constructor.
   *
   * @param atomClasses   int[]
   * @param forceConstant double
   */
  public OutOfPlaneBendType(int[] atomClasses, double forceConstant) {
    super(OPBEND, sortKey(atomClasses));
    this.atomClasses = atomClasses;
    this.forceConstant = forceConstant;
  }

  /**
   * Average two OutOfPlaneBendType instances. The atom classes that define the new type must be
   * supplied.
   *
   * @param outOfPlaneBendType1 a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
   * @param outOfPlaneBendType2 a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
   * @param atomClasses         an array of {@link int} objects.
   * @return a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
   */
  public static OutOfPlaneBendType average(OutOfPlaneBendType outOfPlaneBendType1,
                                           OutOfPlaneBendType outOfPlaneBendType2, int[] atomClasses) {
    if (outOfPlaneBendType1 == null || outOfPlaneBendType2 == null || atomClasses == null) {
      return null;
    }

    double forceConstant =
        (outOfPlaneBendType1.forceConstant + outOfPlaneBendType2.forceConstant) / 2.0;

    return new OutOfPlaneBendType(atomClasses, forceConstant);
  }

  /**
   * Construct an OutOfPlaneBendType from an input string.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return an OutOfPlaneBendType instance.
   */
  public static OutOfPlaneBendType parse(String input, String[] tokens) {
    if (tokens.length < 6) {
      logger.log(Level.WARNING, "Invalid OPBEND type:\n{0}", input);
    } else {
      try {
        int[] atomClasses = new int[4];
        atomClasses[0] = parseInt(tokens[1]);
        atomClasses[1] = parseInt(tokens[2]);
        atomClasses[2] = parseInt(tokens[3]);
        atomClasses[3] = parseInt(tokens[4]);
        double forceConstant = parseDouble(tokens[5]);
        return new OutOfPlaneBendType(atomClasses, forceConstant);
      } catch (NumberFormatException e) {
        String message = "Exception parsing OPBEND type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * This method sorts the atom classes for the out-of-plane angle bending type.
   *
   * @param c atomClasses
   * @return lookup key
   */
  public static String sortKey(int[] c) {
    if (c == null || c.length != 4) {
      return null;
    }
    return c[0] + " " + c[1] + " " + c[2] + " " + c[3];
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compare(String s1, String s2) {
    String[] keys1 = s1.split(" ");
    String[] keys2 = s2.split(" ");

    for (int i = 0; i < 4; i++) {
      int c1 = parseInt(keys1[i]);
      int c2 = parseInt(keys2[i]);
      if (c1 < c2) {
        return -1;
      } else if (c1 > c2) {
        return 1;
      }
    }
    return 0;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    OutOfPlaneBendType outOfPlaneBendType = (OutOfPlaneBendType) o;
    return Arrays.equals(atomClasses, outOfPlaneBendType.atomClasses);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int hashCode() {
    return Arrays.hashCode(atomClasses);
  }

  /**
   * incrementClasses
   *
   * @param increment The increment to apply to the atom classes.
   */
  public void incrementClasses(int increment) {
    for (int i = 0; i < atomClasses.length; i++) {
      if (atomClasses[i] > 0) {
        atomClasses[i] += increment;
      }
    }
    setKey(sortKey(atomClasses));
  }

  /**
   * Remap new atom classes to known internal ones.
   *
   * @param typeMap a lookup between new atom types and known atom types.
   * @return a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
   */
  public OutOfPlaneBendType patchClasses(HashMap<AtomType, AtomType> typeMap) {
    int count = 0;
    int len = atomClasses.length;

    // Look for new OutOfPlaneBends that contain 1 mapped atom classes.
    for (AtomType newType : typeMap.keySet()) {

      for (int atomClass : atomClasses) {
        if (atomClass == newType.atomClass) {
          count++;
        }
      }
    }

    // If found, create a new OutOfPlaneBend that bridges to known classes.
    if (count == 1) {
      int[] newClasses = copyOf(atomClasses, len);
      for (AtomType newType : typeMap.keySet()) {
        for (int i = 0; i < len; i++) {
          if (atomClasses[i] == newType.atomClass) {
            AtomType knownType = typeMap.get(newType);
            newClasses[i] = knownType.atomClass;
          }
        }
      }
      return new OutOfPlaneBendType(newClasses, forceConstant);
    }
    return null;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted out-of-plane angle bending string.
   */
  @Override
  public String toString() {
    return String.format("opbend  %5d  %5d  %5d  %5d  %6.2f", atomClasses[0], atomClasses[1],
        atomClasses[2], atomClasses[3], forceConstant);
  }

  /**
   * Create an AmoebaOutOfPlaneBendForce instance.
   *
   * @param doc        the Document instance.
   * @param forceField the ForceField to grab constants from.
   */
  public static Element getXMLForce(Document doc, ForceField forceField) {
    Map<String, OutOfPlaneBendType> types = forceField.getOutOfPlaneBendTypes();
    if (!types.values().isEmpty()) {
      Element node = doc.createElement("AmoebaOutOfPlaneBendForce");
      node.setAttribute("type", forceField.getString("opbendtype", "ALLINGER"));
      node.setAttribute("opbend-cubic", valueOf(forceField.getDouble("opbend-cubic", DEFAULT_OPBEND_CUBIC)));
      node.setAttribute("opbend-quartic", valueOf(forceField.getDouble("opbend-quartic", DEFAULT_OPBEND_QUARTIC)));
      node.setAttribute("opbend-pentic", valueOf(forceField.getDouble("opbend-pentic", DEFAULT_OPBEND_PENTIC)));
      node.setAttribute("opbend-sextic", valueOf(forceField.getDouble("opbend-sextic", DEFAULT_OPBEND_SEXTIC)));
      for (OutOfPlaneBendType outOfPlaneBendType : types.values()) {
        node.appendChild(outOfPlaneBendType.toXML(doc));
      }
      return node;
    }
    return null;
  }

  /**
   * Write OutOfPlaneBendType to OpenMM XML format.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("Angle");
    int i = 1;
    for (int ac : atomClasses) {
      if (ac == 0) {
        node.setAttribute(format("class%d", i), "");
      } else {
        node.setAttribute(format("class%d", i), format("%d", ac));
      }
      i++;
    }
    // Convert Kcal/mol/radian^2 to KJ/mol/deg^2
    node.setAttribute("k", format("%f", forceConstant * KCAL_TO_KJ / (DEGREES_PER_RADIAN * DEGREES_PER_RADIAN)));
    return node;
  }
}
