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

import static ffx.potential.parameters.ForceField.ForceFieldType.BOND;
import static ffx.utilities.Constants.ANG_TO_NM;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.Constants.NM_TO_ANG;
import static ffx.utilities.PropertyGroup.EnergyUnitConversion;
import static ffx.utilities.PropertyGroup.LocalGeometryFunctionalForm;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

/**
 * The BondType class defines one harmonic bond stretch energy term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "bond", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [2 integers and 2 reals]
    Provides the values for a single bond stretching parameter.
    The integer modifiers give the atom class numbers for the two kinds of atoms involved in the bond which is to be defined.
    The real number modifiers give the force constant value for the bond and the ideal bond length in Angstroms.
    The default value of 1.0 is used, if the bondunit keyword is not given in the force field parameter file or the keyfile.
    """)
public final class BondType extends BaseType implements Comparator<String> {

  public static final double DEFAULT_BOND_UNIT = 1.0;

  /**
   * Default cubic coefficient in bond stretch potential.
   */
  public static final double DEFAULT_BOND_CUBIC = 0.0;
  /**
   * Default quartic coefficient in bond stretch potential.
   */
  public static final double DEFAULT_BOND_QUARTIC = 0.0;

  /**
   * Convert bond stretch energy to kcal/mole.
   */
  @FFXProperty(name = "bondunit", propertyGroup = EnergyUnitConversion, defaultValue = "1.0", description = """
      Sets the scale factor needed to convert the energy value computed by the bond stretching potential into units of kcal/mole.
      The correct value is force field dependent and typically provided in the header of the master force field parameter file.
      """)
  public double bondUnit = DEFAULT_BOND_UNIT;

  /**
   * Cubic coefficient in bond stretch potential.
   */
  @FFXProperty(name = "bond-cubic", propertyGroup = LocalGeometryFunctionalForm, defaultValue = "0.0", description = """
      Sets the value of the cubic term in the Taylor series expansion form of the bond stretching potential energy.
      The real number modifier gives the value of the coefficient as a multiple of the quadratic coefficient.
      This term multiplied by the bond stretching energy unit conversion factor, the force constant,
      and the cube of the deviation of the bond length from its ideal value gives the cubic contribution to the bond stretching energy.
      The default value in the absence of the bond-cubic keyword is zero; i.e., the cubic bond stretching term is omitted.
      """)
  public double cubic = DEFAULT_BOND_CUBIC;

  /**
   * Quartic coefficient in bond stretch potential.
   */
  @FFXProperty(name = "bond-quartic", propertyGroup = LocalGeometryFunctionalForm, defaultValue = "0.0", description = """
      Sets the value of the quartic term in the Taylor series expansion form of the bond stretching potential energy.
      The real number modifier gives the value of the coefficient as a multiple of the quadratic coefficient.
      This term multiplied by the bond stretching energy unit conversion factor, the force constant,
      and the forth power of the deviation of the bond length from its ideal value gives the quartic contribution to the bond stretching energy.
      The default value in the absence of the bond-quartic keyword is zero; i.e., the quartic bond stretching term is omitted.
      """)
  public double quartic = DEFAULT_BOND_QUARTIC;

  /**
   * A Logger for the BondType class.
   */
  private static final Logger logger = Logger.getLogger(BondType.class.getName());
  /**
   * Atom classes that form this bond stretch.
   */
  public final int[] atomClasses;
  /**
   * Force constant (Kcal/mol/A^2).
   */
  public final double forceConstant;
  /**
   * Equilibrium separation (Angstroms).
   */
  public final double distance;
  /**
   * Radius of a flat bottom where energy and force is 0; typically used for restraints. Will almost
   * always be 0.
   */
  public final double flatBottomRadius;
  /**
   * The function used by the bond: harmonic or quartic with flat-bottom variants.
   */
  public BondFunction bondFunction;

  /**
   * The default BondType constructor defines use of the Quartic BondFunction.
   *
   * @param atomClasses   int[]
   * @param forceConstant double
   * @param distance      double
   */
  public BondType(int[] atomClasses, double forceConstant, double distance) {
    this(atomClasses, forceConstant, distance, BondFunction.QUARTIC, 0.0);
  }

  /**
   * BondType constructor.
   *
   * @param atomClasses   int[]
   * @param forceConstant double
   * @param distance      double
   * @param bondFunction  the BondFunction type to apply.
   */
  public BondType(int[] atomClasses, double forceConstant, double distance,
                  BondFunction bondFunction) {
    this(atomClasses, forceConstant, distance, bondFunction, 0.0);
  }

  /**
   * BondType constructor.
   *
   * @param atomClasses      int[]
   * @param forceConstant    double
   * @param distance         double
   * @param bondFunction     the BondFunction type to apply.
   * @param flatBottomRadius a double.
   */
  public BondType(int[] atomClasses, double forceConstant, double distance, BondFunction bondFunction,
                  double flatBottomRadius) {
    super(BOND, sortKey(atomClasses));
    this.atomClasses = atomClasses;
    this.forceConstant = forceConstant;
    this.distance = distance;
    this.bondFunction = bondFunction;
    this.flatBottomRadius = flatBottomRadius;
    assert (flatBottomRadius == 0 || bondFunction == BondFunction.FLAT_BOTTOM_HARMONIC);
  }

  /**
   * Average two BondType instances. The atom classes that define the new type must be supplied.
   *
   * @param bondType1   a {@link ffx.potential.parameters.BondType} object.
   * @param bondType2   a {@link ffx.potential.parameters.BondType} object.
   * @param atomClasses an array of {@link int} objects.
   * @return a {@link ffx.potential.parameters.BondType} object.
   */
  public static BondType average(BondType bondType1, BondType bondType2, int[] atomClasses) {
    if (bondType1 == null || bondType2 == null || atomClasses == null) {
      return null;
    }
    BondType.BondFunction bondFunction = bondType1.bondFunction;
    if (bondFunction != bondType2.bondFunction) {
      return null;
    }

    double forceConstant = (bondType1.forceConstant + bondType2.forceConstant) / 2.0;
    double distance = (bondType1.distance + bondType2.distance) / 2.0;

    return new BondType(atomClasses, forceConstant, distance, bondFunction);
  }

  /**
   * Construct a BondType from an input string.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a BondType instance.
   */
  public static BondType parse(String input, String[] tokens) {
    if (tokens.length < 5) {
      logger.log(Level.WARNING, "Invalid BOND type:\n{0}", input);
    } else {
      try {
        int[] atomClasses = new int[2];
        atomClasses[0] = parseInt(tokens[1]);
        atomClasses[1] = parseInt(tokens[2]);
        double forceConstant = parseDouble(tokens[3]);
        double distance = parseDouble(tokens[4]);
        return new BondType(atomClasses, forceConstant, distance);
      } catch (NumberFormatException e) {
        String message = "Exception parsing BOND type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * This method sorts the atom classes as: min, max
   *
   * @param c atomClasses
   * @return lookup key
   */
  public static String sortKey(int[] c) {
    if (c == null || c.length != 2) {
      return null;
    }

    int temp;
    if (c[1] <= c[0]) {
      temp = c[1];
      c[1] = c[0];
      c[0] = temp;
    }

    return c[0] + " " + c[1];
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compare(String key1, String key2) {
    String[] keys1 = key1.split(" ");
    String[] keys2 = key2.split(" ");
    int[] c1 = new int[2];
    int[] c2 = new int[2];
    for (int i = 0; i < 2; i++) {
      c1[i] = parseInt(keys1[i]);
      c2[i] = parseInt(keys2[i]);
    }

    if (c1[0] < c2[0]) {
      return -1;
    } else if (c1[0] > c2[0]) {
      return 1;
    } else if (c1[1] < c2[1]) {
      return -1;
    } else if (c1[1] > c2[1]) {
      return 1;
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
    BondType bondType = (BondType) o;
    return Arrays.equals(atomClasses, bondType.atomClasses);
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
      atomClasses[i] += increment;
    }
    setKey(sortKey(atomClasses));
  }

  /**
   * Remap new atom classes to known internal ones.
   *
   * @param typeMap a lookup between new atom types and known atom types.
   * @return a {@link ffx.potential.parameters.BondType} object.
   */
  public BondType patchClasses(HashMap<AtomType, AtomType> typeMap) {

    int count = 0;
    int len = atomClasses.length;

    // Look for new BondTypes that contain one mapped atom class.
    for (AtomType newType : typeMap.keySet()) {

      for (int atomClass : atomClasses) {
        if (atomClass == newType.atomClass) {
          count++;
        }
      }
    }

    // If found, create a new BondType that bridges to known classes.
    if (count == 1) {
      int[] newClasses = Arrays.copyOf(atomClasses, len);
      for (AtomType newType : typeMap.keySet()) {
        for (int i = 0; i < len; i++) {
          if (atomClasses[i] == newType.atomClass) {
            AtomType knownType = typeMap.get(newType);
            newClasses[i] = knownType.atomClass;
          }
        }
      }

      return new BondType(newClasses, forceConstant, distance, bondFunction);
    }
    return null;
  }

  public void setBondFunction(BondFunction bondFunction) {
    this.bondFunction = bondFunction;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted bond stretch string.
   */
  @Override
  public String toString() {
    return format("bond  %5d  %5d  %7.2f  %7.4f", atomClasses[0], atomClasses[1], forceConstant,
        distance);
  }

  public static Element getXMLForce(Document doc, ForceField forceField) {
    Map<String, BondType> btMap = forceField.getBondTypes();
    if (!btMap.values().isEmpty()) {
      Element node = doc.createElement("AmoebaBondForce");
      node.setAttribute("bond-cubic", format("%f", forceField.getDouble("bond-cubic", DEFAULT_BOND_CUBIC) * NM_TO_ANG));
      node.setAttribute("bond-quartic", format("%f", forceField.getDouble("bond-quartic", DEFAULT_BOND_QUARTIC) * NM_TO_ANG * NM_TO_ANG));
      for (BondType bondType : btMap.values()) {
        node.appendChild(bondType.toXML(doc));
      }
      return node;
    }
    return null;
  }

  /**
   * Write BondType to OpenMM XML format.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("Bond");
    node.setAttribute("class1", format("%d", atomClasses[0]));
    node.setAttribute("class2", format("%d", atomClasses[1]));
    node.setAttribute("length", format("%f", distance * ANG_TO_NM));
    node.setAttribute("k", format("%f", forceConstant * NM_TO_ANG * NM_TO_ANG * KCAL_TO_KJ));
    return node;
  }

  /**
   * Describes the function used by the bond. Currently supported: harmonic and quartic with
   * flat-bottom variants.
   */
  public enum BondFunction {

    // Harmonic bond function.
    HARMONIC("0.5*k*(r-r0)^2", false),

    // Quartic bond function.
    QUARTIC("0.5*k*dv^2*((1+cubic)*dv+(1+quartic)*dv^2);dv=r-r0", false),

    // Flat bottom harmonic bond function for restraints.
    FLAT_BOTTOM_HARMONIC(
        "0.5*k*dv^2;dv=step(dv)*step(dv-fb)*(dv-fb)" + "+step(-dv)*step(-dv-fb)*(-dv-fb);dv=r-r0",
        true),

    // Flat bottom Quartic bond function for restraints.
    FLAT_BOTTOM_QUARTIC("0.5*k*dv^2*((1+cubic)*dv+(1+quartic)*dv^2);"
        + "dv=step(dv)*step(dv-fb)*(dv-fb)+step(-dv)*step(-dv-fb)*(-dv-fb);dv=r-r0", true);

    /**
     * String representation of the mathematical form.
     */
    private final String mathematicalForm;

    /**
     * Flag to indicate if the bond function has a flat bottom.
     */
    private final boolean hasFlatBottom;

    /**
     * BondFunction constructor.
     *
     * @param mathForm   String representation of the mathematical form.
     * @param flatBottom Flag to indicate if the bond function has a flat bottom.
     */
    BondFunction(String mathForm, boolean flatBottom) {
      this.mathematicalForm = mathForm;
      this.hasFlatBottom = flatBottom;
    }

    /**
     * Returns whether this BondFunction has a flat bottom.
     *
     * @return Flat bottom.
     */
    public boolean hasFlatBottom() {
      return hasFlatBottom;
    }

    /**
     * Returns the form of this bond as a mathematical expression parsable by OpenMM.
     *
     * @return A string describing mathematical form.
     */
    public String toMathematicalForm() {
      return mathematicalForm;
    }
  }
}
