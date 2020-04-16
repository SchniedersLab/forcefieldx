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

import static ffx.potential.parameters.ForceField.ForceFieldType.BOND;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The BondType class defines one harmonic bond stretch energy term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class BondType extends BaseType implements Comparator<String> {

  /** Convert bond stretch energy to kcal/mole. */
  public static final double units = 1.0;
  /** Cubic coefficient in bond stretch potential. */
  public static final double cubic = -2.55;
  /** Quartic coefficient in bond stretch potential. */
  public static final double quartic = 3.793125;
  /** A Logger for the BondType class. */
  private static final Logger logger = Logger.getLogger(BondType.class.getName());
  /** Atom classes that form this bond stretch. */
  public final int[] atomClasses;
  /** Force constant (Kcal/mol). */
  public final double forceConstant;
  /** Equilibrium separation (Angstroms). */
  public final double distance;
  /**
   * Radius of a flat bottom where energy and force is 0; typically used for restraints. Will almost
   * always be 0.
   */
  public final double flatBottomRadius;
  /** The function used by the bond: harmonic or quartic with flat-bottom variants. */
  public BondFunction bondFunction;

  /**
   * The default BondType constructor defines use of the Quartic BondFunction.
   *
   * @param atomClasses int[]
   * @param forceConstant double
   * @param distance double
   */
  public BondType(int[] atomClasses, double forceConstant, double distance) {
    this(atomClasses, forceConstant, distance, BondFunction.QUARTIC, 0.0);
  }

  /**
   * BondType constructor.
   *
   * @param atomClasses int[]
   * @param forceConstant double
   * @param distance double
   * @param bondFunction the BondFunction type to apply.
   */
  public BondType(
      int[] atomClasses, double forceConstant, double distance, BondFunction bondFunction) {
    this(atomClasses, forceConstant, distance, bondFunction, 0.0);
  }

  /**
   * BondType constructor.
   *
   * @param atomClasses int[]
   * @param forceConstant double
   * @param distance double
   * @param bondFunction the BondFunction type to apply.
   * @param flatBottomRadius a double.
   */
  public BondType(
      int[] atomClasses,
      double forceConstant,
      double distance,
      BondFunction bondFunction,
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
   * @param bondType1 a {@link ffx.potential.parameters.BondType} object.
   * @param bondType2 a {@link ffx.potential.parameters.BondType} object.
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
   * @param input The overall input String.
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

  /** {@inheritDoc} */
  @Override
  public int compare(String key1, String key2) {
    String[] keys1 = key1.split(" ");
    String[] keys2 = key2.split(" ");
    int[] c1 = new int[2];
    int[] c2 = new int[2];
    for (int i = 0; i < 2; i++) {
      c1[i] = Integer.parseInt(keys1[i]);
      c2[i] = Integer.parseInt(keys2[i]);
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

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    BondType bondType = (BondType) o;
    return Arrays.equals(atomClasses, bondType.atomClasses);
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Arrays.hashCode(atomClasses);
  }

  /**
   * incrementClasses
   *
   * @param increment a int.
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
    return String.format(
        "bond  %5d  %5d  %7.2f  %7.4f", atomClasses[0], atomClasses[1], forceConstant, distance);
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
    FLAT_BOTTOM_QUARTIC(
        "0.5*k*dv^2*((1+cubic)*dv+(1+quartic)*dv^2);"
            + "dv=step(dv)*step(dv-fb)*(dv-fb)+step(-dv)*step(-dv-fb)*(-dv-fb);dv=r-r0",
        true);

    /** String representation of the mathematical form. */
    private final String mathematicalForm;

    /** Flag to indicate if the bond function has a flat bottom. */
    private final boolean hasFlatBottom;

    /**
     * BondFunction constructor.
     *
     * @param mathForm String representation of the mathematical form.
     */
    BondFunction(String mathForm) {
      this(mathForm, false);
    }

    /**
     * BondFunction constructor.
     *
     * @param mathForm String representation of the mathematical form.
     * @param flatBottom Flag to indicate if the bond function has a flat bottom.
     */
    BondFunction(String mathForm, boolean flatBottom) {
      this.mathematicalForm = mathForm;
      this.hasFlatBottom = flatBottom;
    }

    /**
     * Returns whether or not this BondFunction has a flat bottom.
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
