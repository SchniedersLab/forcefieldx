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

import static ffx.potential.parameters.ForceField.ForceFieldType.PITORS;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.PropertyGroup.EnergyUnitConversion;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.String.valueOf;

/**
 * The PiOrbitalTorsionType class defines a Pi-Orbital Torsion energy term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "pitors", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [2 integers and 1 real]
    Provides the values for a single pi-orbital torsional angle potential parameter.
    The two integer modifiers give the atom class numbers for the atoms involved in the central bond of the torsional angle to be parameterized.
    The real modifier gives the value of the 2-fold Fourier amplitude for the torsional angle between p-orbitals centered on the defined bond atom classes.
    The default units for the stretch-torsion force constant can be controlled via the pitorsunit keyword.
    """)
public final class PiOrbitalTorsionType extends BaseType implements Comparator<String> {

  /**
   * A Logger for the PiTorsionType class.
   */
  private static final Logger logger = Logger.getLogger(PiOrbitalTorsionType.class.getName());

  public static final double DEFAULT_PITORS_UNIT = 1.0;

  /**
   * Convert Pi-Torsion energy to kcal/mole.
   */
  @FFXProperty(name = "pitorsunit", propertyGroup = EnergyUnitConversion, defaultValue = "1.0", description = """
      Sets the scale factor needed to convert the energy value computed by the pi-orbital torsional angle potential into units of kcal/mole.
      The correct value is force field dependent and typically provided in the header of the master force field parameter file.
      """)
  public double piTorsUnit = DEFAULT_PITORS_UNIT;

  /**
   * Atom classes that form this Pi-Torsion.
   */
  public final int[] atomClasses;

  /**
   * Force constant.
   */
  public double forceConstant;

  /**
   * PiTorsionType Constructor.
   *
   * @param atomClasses   int[]
   * @param forceConstant double
   */
  public PiOrbitalTorsionType(int[] atomClasses, double forceConstant) {
    super(PITORS, sortKey(atomClasses));
    this.atomClasses = atomClasses;
    this.forceConstant = forceConstant;
  }

  /**
   * Average two PiTorsionType instances. The atom classes that define the new type must be
   * supplied.
   *
   * @param piOrbitalTorsionType1 a {@link PiOrbitalTorsionType} object.
   * @param piOrbitalTorsionType2 a {@link PiOrbitalTorsionType} object.
   * @param atomClasses           an array of {@link int} objects.
   * @return a {@link PiOrbitalTorsionType} object.
   */
  public static PiOrbitalTorsionType average(PiOrbitalTorsionType piOrbitalTorsionType1,
                                             PiOrbitalTorsionType piOrbitalTorsionType2, int[] atomClasses) {
    if (piOrbitalTorsionType1 == null || piOrbitalTorsionType2 == null || atomClasses == null) {
      return null;
    }

    double forceConstant =
        (piOrbitalTorsionType1.forceConstant + piOrbitalTorsionType2.forceConstant) / 2.0;

    return new PiOrbitalTorsionType(atomClasses, forceConstant);
  }

  /**
   * Construct a PiTorsionType from an input string.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a PiTorsionType instance.
   */
  public static PiOrbitalTorsionType parse(String input, String[] tokens) {
    if (tokens.length < 4) {
      logger.log(Level.WARNING, "Invalid PITORS type:\n{0}", input);
    } else {
      try {
        int[] atomClasses = new int[2];
        atomClasses[0] = parseInt(tokens[1]);
        atomClasses[1] = parseInt(tokens[2]);
        double forceConstant = parseDouble(tokens[3]);
        return new PiOrbitalTorsionType(atomClasses, forceConstant);
      } catch (NumberFormatException e) {
        String message = "Exception parsing PITORS type:\n" + input + "\n";
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
  public int compare(String s1, String s2) {
    String[] keys1 = s1.split(" ");
    String[] keys2 = s2.split(" ");

    for (int i = 0; i < 2; i++) {
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
    PiOrbitalTorsionType piOrbitalTorsionType = (PiOrbitalTorsionType) o;
    return Arrays.equals(atomClasses, piOrbitalTorsionType.atomClasses);
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
   * @param increment The increment to add to the atom classes.
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
   * @return a {@link PiOrbitalTorsionType} object.
   */
  public PiOrbitalTorsionType patchClasses(HashMap<AtomType, AtomType> typeMap) {
    int count = 0;
    int len = atomClasses.length;

    // Look for new PiTorsions that contain 1 mapped atom classes.
    for (AtomType newType : typeMap.keySet()) {

      for (int atomClass : atomClasses) {
        if (atomClass == newType.atomClass) {
          count++;
        }
      }
    }

    // If found, create a new PiTorsion that bridges to known classes.
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
      return new PiOrbitalTorsionType(newClasses, forceConstant);
    }
    return null;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted Pi-Torsion type.
   */
  @Override
  public String toString() {
    return format("pitors  %5d  %5d  %4.2f", atomClasses[0], atomClasses[1], forceConstant);
  }

  /**
   * Create an AmoebaPiTorsionForce Element.
   *
   * @param doc        the Document instance.
   * @param forceField the ForceField instance to grab constants from.
   * @return the AmoebaPiTorsionForce Element.
   */
  public static Element getXMLElement(Document doc, ForceField forceField) {
    Map<String, PiOrbitalTorsionType> types = forceField.getPiOrbitalTorsionTypes();
    if (!types.values().isEmpty()) {
      Element node = doc.createElement("AmoebaPiTorsionForce");
      node.setAttribute("piTorsionUnit", valueOf(forceField.getDouble("pitorsunit", DEFAULT_PITORS_UNIT)));
      for (PiOrbitalTorsionType piOrbitalTorsionType : types.values()) {
        node.appendChild(piOrbitalTorsionType.toXML(doc));
      }
      return node;
    }
    return null;
  }

  /**
   * Write PiOrbitalTorsionType to OpenMM XML format.
   *
   * @param doc the Document instance.
   * @return the PiTorsion Element.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("PiTorsion");
    node.setAttribute("class1", format("%d", atomClasses[0]));
    node.setAttribute("class2", format("%d", atomClasses[1]));
    // OpenMM has the piTorsUnit hard coded to 1.0 (the default).
    node.setAttribute("k", format("%f", forceConstant * KCAL_TO_KJ));
    return node;
  }

}
