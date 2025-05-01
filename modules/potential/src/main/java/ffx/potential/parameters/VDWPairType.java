// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import ffx.potential.parameters.ForceField.ForceFieldType;
import ffx.utilities.FFXProperty;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.Arrays;
import java.util.Comparator;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.Constants.ANG_TO_NM;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

/**
 * The VDWPairType class defines a van der Waals Pair type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "vdwpr", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [2 integers and 2 reals]
    Provides the values for the vdw parameters for a single special heteroatomic pair of atoms.
    The integer modifiers give the pair of atom class numbers for which special vdw parameters are to be defined.
    The two real number modifiers give the values of the minimum energy contact distance in Angstroms and the well depth at the minimum distance in kcal/mole.
    """)
public final class VDWPairType extends BaseType implements Comparator<String> {

  private static final Logger logger = Logger.getLogger(VDWPairType.class.getName());
  /**
   * The radius of the minimum well depth energy (angstroms).
   */
  public final double radius;
  /**
   * The minimum energy of the vdw function (kcal/mol).
   */
  public final double wellDepth;
  /**
   * Atom classes that form this bond stretch.
   */
  public final int[] atomClasses;

  /**
   * van der Waals constructor. If the reduction factor is .LE. 0.0, no reduction is used for this
   * atom type.
   *
   * @param atomClasses The atom class that uses this van der Waals Pair.
   * @param radius      The radius of the minimum well depth energy (angstroms).
   * @param wellDepth   The minimum energy of the vdw function (kcal/mol).
   */
  public VDWPairType(int[] atomClasses, double radius, double wellDepth) {
    super(ForceFieldType.VDWPR, sortKey(atomClasses));
    this.atomClasses = atomClasses;
    this.radius = radius;
    this.wellDepth = wellDepth;
  }

  /**
   * Average.
   *
   * @param vdwType1    a {@link VDWPairType} object.
   * @param vdwType2    a {@link VDWPairType} object.
   * @param atomClasses The atom classes that uses this van der Waals Pair.
   * @return a {@link VDWPairType} object.
   */
  public static VDWPairType average(VDWPairType vdwType1, VDWPairType vdwType2, int[] atomClasses) {
    if (vdwType1 == null || vdwType2 == null) {
      return null;
    }

    double radius = (vdwType1.radius + vdwType2.radius) / 2.0;
    double wellDepth = (vdwType1.wellDepth + vdwType2.wellDepth) / 2.0;

    return new VDWPairType(atomClasses, radius, wellDepth);
  }

  /**
   * Construct a VDWPairType from multiple input lines.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a VDWType instance.
   */
  public static VDWPairType parse(String input, String[] tokens) {
    if (tokens.length < 5) {
      logger.log(Level.WARNING, "Invalid VDWPR type:\n{0}", input);
    } else {
      try {
        int atomClass1 = parseInt(tokens[1]);
        int atomClass2 = parseInt(tokens[2]);
        double radius = parseDouble(tokens[3]);
        double wellDepth = parseDouble(tokens[4]);
        return new VDWPairType(new int[]{atomClass1, atomClass2}, radius, wellDepth);
      } catch (NumberFormatException e) {
        String message = "Exception parsing VDWPR type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
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
    VDWPairType vdwPairType = (VDWPairType) o;
    return Arrays.equals(atomClasses, vdwPairType.atomClasses);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int hashCode() {
    return Arrays.hashCode(atomClasses);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted van der Waals type.
   */
  @Override
  public String toString() {
    return format("vdwpr  %5d  %5d  %11.9f  %11.9f", atomClasses[0], atomClasses[1], radius, wellDepth);
  }

  /**
   * Write VDWPairType to OpenMM XML format.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("Pair");
    node.setAttribute("class1", format("%d", atomClasses[0]));
    node.setAttribute("class2", format("%d", atomClasses[1]));
    // Convert Angstroms to nm.
    node.setAttribute("sigma", format("%.17f", radius * ANG_TO_NM));
    // Convert Kcal/mol to KJ/mol
    node.setAttribute("epsilon", format("%.17f", wellDepth * KCAL_TO_KJ));
    return node;
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
   * Increment the atom classes by a specified amount.
   *
   * @param increment The increment to add to the atom classes.
   */
  public void incrementClasses(int increment) {
    for (int i = 0; i < atomClasses.length; i++) {
      atomClasses[i] += increment;
    }
    setKey(sortKey(atomClasses));
  }

}
