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

import static ffx.potential.parameters.ForceField.ForceFieldType.IMPROPER;
import static ffx.potential.parameters.ForceField.ForceFieldType.TORSION;
import static ffx.utilities.Constants.DEGREES_PER_RADIAN;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.PropertyGroup.EnergyUnitConversion;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * The TorsionType class defines a torsional angle.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "improper", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [4 integers and 2 reals]"
    Provides the values for a single CHARMM-style improper dihedral angle parameter.
    The integer modifiers give the atom class numbers for the four kinds of atoms involved in the torsion which is to be defined.
    The real number modifiers give the force constant value for the deviation from the target improper torsional angle, and the target value for the torsional angle, respectively.
    The default units for the improper force constant are kcal/mole/radian^2, but this can be controlled via the impropunit keyword.
    """)
@FFXProperty(name = "torsion", clazz = String.class, propertyGroup = PotentialFunctionParameter, description = """
    [4 integers and up to 6 real/real/integer triples]
    Provides the values for a single torsional angle parameter.
    The first four integer modifiers give the atom class numbers for the atoms involved in the torsional angle to be defined.
    Each of the remaining triples of real/real/integer modifiers give the amplitude,
    phase offset in degrees and periodicity of a particular torsional function term, respectively.
    Periodicities through 6-fold are allowed for torsional parameters.
    """)
public final class TorsionType extends BaseType implements Comparator<String> {

  private static final Logger logger = Logger.getLogger(TorsionType.class.getName());
  /**
   * Atom classes that for this Torsion angle.
   */
  public final int[] atomClasses;
  /**
   * Number of terms in the Fourier series.
   */
  public final int terms;
  /**
   * Amplitudes of the Fourier series.
   */
  public final double[] amplitude;
  /**
   * Phases of the Fourier series in degrees.
   */
  public final double[] phase;
  /**
   * Cosine of the phase angle.
   */
  public final double[] cosine;
  /**
   * Sine of the phase angle.
   */
  public final double[] sine;
  /**
   * Periodicity of the Fourier series.
   */
  public final int[] periodicity;
  /**
   * The torsion mode in use.
   */
  private final TorsionMode torsionMode;

  public static final double DEFAULT_TORSION_UNIT = 1.0;
  /**
   * Unit conversion.
   */
  @FFXProperty(name = "torsionunit", propertyGroup = EnergyUnitConversion, defaultValue = "1.0", description = """
      Sets the scale factor needed to convert the energy value computed by the torsional angle potential into units of kcal/mole.
      The correct value is force field dependent and typically provided in the header of the master force field parameter file.
      """)
  public double torsionUnit = DEFAULT_TORSION_UNIT;

  /**
   * TorsionType Constructor.
   *
   * @param atomClasses Atom classes.
   * @param amplitude   Amplitudes of the Fourier series.
   * @param phase       Phases of the Fourier series in degrees.
   * @param periodicity Periodicity of the Fourier series.
   */
  public TorsionType(int[] atomClasses, double[] amplitude, double[] phase, int[] periodicity) {
    this(atomClasses, amplitude, phase, periodicity, TorsionMode.NORMAL);
  }

  /**
   * TorsionType Constructor.
   *
   * @param atomClasses Atom classes.
   * @param amplitude   Amplitudes of the Fourier series.
   * @param phase       Phases of the Fourier series in degrees.
   * @param periodicity Periodicity of the Fourier series.
   * @param torsionMode Define the TorsionMode for this TorsionType.
   */
  public TorsionType(int[] atomClasses, double[] amplitude, double[] phase, int[] periodicity,
                     TorsionMode torsionMode) {
    super(TORSION, sortKey(atomClasses));
    this.atomClasses = atomClasses;
    int max = 1;
    for (int i1 : periodicity) {
      if (i1 > max) {
        max = i1;
      }
    }
    terms = max;
    this.amplitude = new double[terms];
    this.phase = new double[terms];
    this.periodicity = new int[terms];
    for (int i = 0; i < terms; i++) {
      this.periodicity[i] = i + 1;
    }
    for (int i = 0; i < amplitude.length; i++) {
      int j = periodicity[i] - 1;
      if (j >= 0 && j < terms) {
        this.amplitude[j] = amplitude[i];
        this.phase[j] = phase[i];
      }
    }

    cosine = new double[terms];
    sine = new double[terms];
    for (int i = 0; i < terms; i++) {
      double angle = toRadians(this.phase[i]);
      cosine[i] = cos(angle);
      sine[i] = sin(angle);
    }

    this.torsionMode = torsionMode;
    if (torsionMode == TorsionMode.IMPROPER) {
      forceFieldType = IMPROPER;
    }
  }

  /**
   * average.
   *
   * @param torsionType1 a {@link ffx.potential.parameters.TorsionType} object.
   * @param torsionType2 a {@link ffx.potential.parameters.TorsionType} object.
   * @param atomClasses  an array of {@link int} objects.
   * @return a {@link ffx.potential.parameters.TorsionType} object.
   */
  public static TorsionType average(TorsionType torsionType1, TorsionType torsionType2,
                                    int[] atomClasses) {
    if (torsionType1 == null || torsionType2 == null || atomClasses == null) {
      return null;
    }
    int len = torsionType1.amplitude.length;
    if (len != torsionType2.amplitude.length) {
      return null;
    }
    double[] amplitude = new double[len];
    double[] phase = new double[len];
    int[] periodicity = new int[len];
    for (int i = 0; i < len; i++) {
      amplitude[i] = (torsionType1.amplitude[i] + torsionType2.amplitude[i]) / 2.0;
      phase[i] = (torsionType1.phase[i] + torsionType2.phase[i]) / 2.0;
      periodicity[i] = (torsionType1.periodicity[i] + torsionType2.periodicity[i]) / 2;
    }

    return new TorsionType(atomClasses, amplitude, phase, periodicity);
  }

  /**
   * Construct a TorsionType from an input string.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a TorsionType instance.
   */
  public static TorsionType parse(String input, String[] tokens) {
    if (tokens.length < 5) {
      logger.log(Level.WARNING, "Invalid TORSION type:\n{0}", input);
    } else {
      try {
        int[] atomClasses = new int[4];
        atomClasses[0] = parseInt(tokens[1]);
        atomClasses[1] = parseInt(tokens[2]);
        atomClasses[2] = parseInt(tokens[3]);
        atomClasses[3] = parseInt(tokens[4]);
        int terms = (tokens.length - 5) / 3;
        double[] amplitude = new double[terms];
        double[] phase = new double[terms];
        int[] periodicity = new int[terms];
        int index = 5;
        for (int i = 0; i < terms; i++) {
          amplitude[i] = parseDouble(tokens[index++]);
          phase[i] = parseDouble(tokens[index++]);
          periodicity[i] = parseInt(tokens[index++]);
        }
        return new TorsionType(atomClasses, amplitude, phase, periodicity);
      } catch (NumberFormatException e) {
        String message = "NumberFormatException parsing TORSION type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      } catch (Exception e) {
        String message = "Exception parsing TORSION type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * Construct a TorsionType with <code>TorsionMode.IMPROPER</code> from an input string.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a TorsionType instance.
   */
  public static TorsionType parseImproper(String input, String[] tokens) {
    if (tokens.length < 5) {
      logger.log(Level.WARNING, "Invalid IMPROPER type:\n{0}", input);
    } else {
      try {
        int[] atomClasses = new int[4];
        atomClasses[0] = parseInt(tokens[1]);
        atomClasses[1] = parseInt(tokens[2]);
        atomClasses[2] = parseInt(tokens[3]);
        atomClasses[3] = parseInt(tokens[4]);
        double[] amplitude = new double[1];
        double[] phase = new double[1];
        int[] periodicity = new int[1];
        int index = 5;
        amplitude[0] = parseDouble(tokens[index++]);
        phase[0] = parseDouble(tokens[index]);
        periodicity[0] = 1;
        return new TorsionType(atomClasses, amplitude, phase, periodicity, TorsionMode.IMPROPER);
      } catch (NumberFormatException e) {
        String message = "Exception parsing IMPROPER type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * This method sorts the atom classes for the torsion.
   *
   * @param c atomClasses
   * @return lookup key
   * @since 1.0
   */
  public static String sortKey(int[] c) {
    if (c == null || c.length != 4) {
      return null;
    }
    if (c[1] < c[2]) {
      // Do nothing.
    } else if (c[2] < c[1]) {
      // Reverse the order.
      int temp = c[0];
      c[0] = c[3];
      c[3] = temp;
      temp = c[1];
      c[1] = c[2];
      c[2] = temp;
    } else if (c[1] == c[2]) {
      if (c[0] > c[3]) {
        // Reverse the order.
        int temp = c[0];
        c[0] = c[3];
        c[3] = temp;
      }
    } else if (c[0] <= c[3]) {
      // Do nothing.
    } else {
      // Reverse the order.
      int temp = c[0];
      c[0] = c[3];
      c[3] = temp;
      temp = c[1];
      c[1] = c[2];
      c[2] = temp;
    }

    return c[0] + " " + c[1] + " " + c[2] + " " + c[3];
  }

  /**
   * {@inheritDoc}
   *
   * @since 1.0
   */
  @Override
  public int compare(String s1, String s2) {
    String[] keys1 = s1.split(" ");
    String[] keys2 = s2.split(" ");
    int[] c1 = new int[4];
    int[] c2 = new int[4];

    for (int i = 0; i < 4; i++) {
      c1[i] = parseInt(keys1[i]);
      c2[i] = parseInt(keys2[i]);
    }

    if (c1[1] < c2[1]) {
      return -1;
    } else if (c1[1] > c2[1]) {
      return 1;
    } else if (c1[2] < c2[2]) {
      return -1;
    } else if (c1[2] > c2[2]) {
      return 1;
    } else if (c1[0] < c2[0]) {
      return -1;
    } else if (c1[0] > c2[0]) {
      return 1;
    } else if (c1[3] < c2[3]) {
      return -1;
    } else if (c1[3] > c2[3]) {
      return 1;
    }

    return 0;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Override the default <code>equals</code> method.
   *
   * @since 1.0
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TorsionType torsionType = (TorsionType) o;
    return Arrays.equals(atomClasses, torsionType.atomClasses);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Implementation of the <code>hashCode</code> method.
   *
   * @since 1.0
   */
  @Override
  public int hashCode() {
    return Arrays.hashCode(atomClasses);
  }

  /**
   * Increment the atom classes by a specified amount.
   *
   * @param increment The increment to add to the atom classes.
   */
  public void incrementClasses(int increment) {
    for (int i = 0; i < atomClasses.length; i++) {
      if (atomClasses[i] != 0) {
        atomClasses[i] += increment;
      }
    }
    setKey(sortKey(atomClasses));
  }

  /**
   * Remap new atom classes to known internal ones.
   *
   * @param typeMap a lookup between new atom types and known atom types.
   * @return a {@link ffx.potential.parameters.TorsionType} object.
   */
  public TorsionType patchClasses(HashMap<AtomType, AtomType> typeMap) {
    int count = 0;
    int len = atomClasses.length;

    // Look for new TorsionTypes that contain 1 to 3 mapped atom classes.
    for (AtomType newType : typeMap.keySet()) {
      for (int atomClass : atomClasses) {
        if (atomClass == newType.atomClass) {
          count++;
        }
      }
    }

    // If found, create a new TorsionType that bridges to known classes.
    if (count == 1 || count == 2 || count == 3) {
      int[] newClasses = copyOf(atomClasses, len);
      for (AtomType newType : typeMap.keySet()) {
        for (int i = 0; i < len; i++) {
          if (atomClasses[i] == newType.atomClass) {
            AtomType knownType = typeMap.get(newType);
            newClasses[i] = knownType.atomClass;
          }
        }
      }
      return new TorsionType(newClasses, amplitude, phase, periodicity);
    }
    return null;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted Torsion angle.
   *
   * @since 1.0
   */
  @Override
  public String toString() {
    StringBuilder torsionBuffer;
    if (torsionMode == TorsionMode.IMPROPER) {
      torsionBuffer = new StringBuilder("improper");
    } else {
      torsionBuffer = new StringBuilder("torsion");
    }
    for (int i : atomClasses) {
      torsionBuffer.append(format(" %5d", i));
    }

    boolean nonZero = false;
    for (double v : amplitude) {
      if (v != 0.0) {
        nonZero = true;
        break;
      }
    }

    for (int i = 0; i < amplitude.length; i++) {
      if (amplitude[i] == 0.0 && (i > 0 || nonZero)) {
        continue;
      }
      if (torsionMode == TorsionMode.NORMAL) {
        torsionBuffer.append(format(" %7.3f %7.3f %1d", amplitude[i], phase[i], periodicity[i]));
      } else {
        torsionBuffer.append(format(" %7.3f %7.3f", amplitude[i], phase[i]));
      }
    }
    return torsionBuffer.toString();
  }


  /**
   * Create a PeriodicTorsionForce Element.
   *
   * @param doc        the Document instance.
   * @param forceField the ForceField that contains the Torsion types.
   * @return the PeriodicTorsionForce instance.
   */
  public static Element getXMLForce(Document doc, ForceField forceField) {
    Map<String, TorsionType> types = forceField.getTorsionTypes();
    if (!types.values().isEmpty()) {
      Element node = doc.createElement("PeriodicTorsionForce");
      double torsionUnit = forceField.getDouble("torsionunit", DEFAULT_TORSION_UNIT);
      for (TorsionType torsionType : types.values()) {
        node.appendChild(torsionType.toXML(doc, torsionUnit));
      }
      return node;
    }
    return null;
  }

  /**
   * Write TorsionType (Proper) to OpenMM XML format.
   *
   * @param doc      the Document instance.
   * @param torsUnit scale torsion force constants by this factor.
   * @return the Proper torsion Element.
   */
  public Element toXML(Document doc, Double torsUnit) {
    Element node = doc.createElement("Proper");
    // TODO specify proper?
    node.setAttribute("class1", format("%d", atomClasses[0]));
    node.setAttribute("class2", format("%d", atomClasses[1]));
    node.setAttribute("class3", format("%d", atomClasses[2]));
    node.setAttribute("class4", format("%d", atomClasses[3]));
    for (int i = 0; i < amplitude.length; i++) {
      // Convert from Kcal/mol to KJ/mol
      node.setAttribute(format("k%d", i + 1), format("%f", amplitude[i] * KCAL_TO_KJ * torsUnit));
      // Convert from degrees to radians.
      node.setAttribute(format("phase%d", i + 1), format("%f", phase[i] / DEGREES_PER_RADIAN));
      node.setAttribute(format("periodicity%d", i + 1), format("%d", periodicity[i]));
    }
    return node;
  }

  /**
   * Torsion modes include Normal or In-Plane
   */
  public enum TorsionMode {
    NORMAL, IMPROPER
  }
}
