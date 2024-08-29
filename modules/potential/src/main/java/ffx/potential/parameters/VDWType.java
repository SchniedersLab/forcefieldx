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

import java.util.Comparator;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.ForceFieldType.VDW;
import static ffx.potential.parameters.ForceField.ForceFieldType.VDW14;
import static ffx.utilities.Constants.ANG_TO_NM;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.PropertyGroup.PotentialFunctionParameter;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.StrictMath.abs;
import static java.lang.String.format;
import static java.lang.String.valueOf;

/**
 * The VDWType class defines van der Waals type for a normal interaction or a special 1-4
 * interaction.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@FFXProperty(name = "vdw", clazz = String.class, propertyGroup = PotentialFunctionParameter,
    description = """
        [1 integer and 3 reals]
        Provides values for a single van der Waals parameter. The integer modifier, if positive,
        gives the atom class number for which vdw parameters are to be defined. Note that vdw parameters are given for atom classes, not atom types.
        The three real number modifiers give the values of the atom size in Angstroms, homoatomic well depth in kcal/mole,
        and an optional reduction factor for univalent atoms.
        """)
@FFXProperty(name = "vdw14", clazz = String.class, propertyGroup = PotentialFunctionParameter,
    description = """
        [1 integer and 2 reals]
        Provides values for a single van der Waals parameter for use in 1-4 nonbonded interactions.
        The integer modifier, if positive, gives the atom class number for which vdw parameters are to be defined.
        Note that vdw parameters are given for atom classes, not atom types.
        The two real number modifiers give the values of the atom size in Angstroms and the homoatomic well depth in kcal/mole.
        Reduction factors, if used, are carried over from the vdw keyword for the same atom class.
        """)
public final class VDWType extends BaseType implements Comparator<String> {

  private static final Logger logger = Logger.getLogger(VDWType.class.getName());
  /**
   * The default gamma parameter in Halgren’s buffered 14-7 vdw potential energy functional form.
   */
  public static final double DEFAULT_GAMMA = 0.12;
  /**
   * The default delta parameter in Halgren’s buffered 14-7 vdw potential energy functional form.
   */
  public static final double DEFAULT_DELTA = 0.07;
  /**
   * The default epsilon combining rule.
   */
  public static final EPSILON_RULE DEFAULT_EPSILON_RULE = EPSILON_RULE.GEOMETRIC;
  /**
   * The default radius combining rule.
   */
  public static final RADIUS_RULE DEFAULT_RADIUS_RULE = RADIUS_RULE.ARITHMETIC;
  /**
   * The default radius size.
   */
  public static final RADIUS_SIZE DEFAULT_RADIUS_SIZE = RADIUS_SIZE.RADIUS;
  /**
   * The default radius type.
   */
  public static final RADIUS_TYPE DEFAULT_RADIUS_TYPE = RADIUS_TYPE.R_MIN;
  /**
   * The default van der Waals functional form type.
   */
  public static final VDW_TYPE DEFAULT_VDW_TYPE = VDW_TYPE.LENNARD_JONES;
  /**
   * The default van der Waals scale factor for 1-2 (bonded) interactions.
   */
  public static final double DEFAULT_VDW_12_SCALE = 0.0;
  /**
   * The default van der Waals scale factor for 1-3 (angle) interactions.
   */
  public static final double DEFAULT_VDW_13_SCALE = 0.0;
  /**
   * The default van der Waals scale factor for 1-4 (torisonal) interactions.
   */
  public static final double DEFAULT_VDW_14_SCALE = 1.0;

  /**
   * The radius of the minimum well depth energy (angstroms).
   */
  public final double radius;
  /**
   * The minimum energy of the vdw function (kcal/mol).
   */
  public final double wellDepth;
  /**
   * Reduction factor for evaluating van der Waals pairs. Valid range: 0.0 .GT. reduction .LE. 1.0
   * Usually only hydrogen atoms have a reduction factor. Setting the reduction to .LT. 0.0 indicates
   * it is not being used.
   */
  public final double reductionFactor;
  /**
   * The atom class that uses this van der Waals parameter.
   */
  public int atomClass;
  /**
   * Is this a normal vdW parameter or is it for 1-4 interactions.
   */
  private final VDWMode vdwMode;

  /**
   * van der Waals constructor. If the reduction factor is .LE. 0.0, no reduction is used for this
   * atom type.
   *
   * @param atomClass       The atom class that uses this van der Waals parameter.
   * @param radius          The radius of the minimum well depth energy (angstroms).
   * @param wellDepth       The minimum energy of the vdw function (kcal/mol).
   * @param reductionFactor Reduction factor for evaluating van der Waals pairs.
   */
  public VDWType(int atomClass, double radius, double wellDepth, double reductionFactor) {
    this(atomClass, radius, wellDepth, reductionFactor, VDWMode.NORMAL);
  }

  /**
   * van der Waals constructor. If the reduction factor is .LE. 0.0, no reduction is used for this
   * atom type.
   *
   * @param atomClass       The atom class that uses this van der Waals parameter.
   * @param radius          The radius of the minimum well depth energy (angstroms).
   * @param wellDepth       The minimum energy of the vdw function (kcal/mol).
   * @param reductionFactor Reduction factor for evaluating van der Waals pairs.
   * @param vdwMode         The VDWMode to use.
   */
  public VDWType(int atomClass, double radius, double wellDepth, double reductionFactor,
                 VDWMode vdwMode) {
    super(VDW, Integer.toString(atomClass));
    this.atomClass = atomClass;
    this.radius = radius;
    this.wellDepth = abs(wellDepth);
    this.reductionFactor = reductionFactor;
    this.vdwMode = vdwMode;
    if (vdwMode == VDWMode.VDW14) {
      forceFieldType = VDW14;
    }
  }

  /**
   * Average two VDWType objects.
   *
   * @param vdwType1  The first VDWType.
   * @param vdwType2  The second VDWType.
   * @param atomClass The new atom class.
   * @return The new averaged VDWType.
   */
  public static VDWType average(VDWType vdwType1, VDWType vdwType2, int atomClass) {
    if (vdwType1 == null || vdwType2 == null) {
      return null;
    }
    double radius = (vdwType1.radius + vdwType2.radius) / 2.0;
    double wellDepth = (vdwType1.wellDepth + vdwType2.wellDepth) / 2.0;
    double reductionFactor = (vdwType1.reductionFactor + vdwType2.reductionFactor) / 2.0;
    return new VDWType(atomClass, radius, wellDepth, reductionFactor);
  }

  /**
   * Construct a VDWType from multiple input lines.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a VDWType instance.
   */
  public static VDWType parse(String input, String[] tokens) {
    if (tokens.length < 4) {
      logger.log(Level.WARNING, "Invalid VDW type:\n{0}", input);
    } else {
      try {
        int atomType = parseInt(tokens[1]);
        double radius = parseDouble(tokens[2]);
        double wellDepth = parseDouble(tokens[3]);
        double reductionFactor = -1.0;
        if (tokens.length == 5) {
          reductionFactor = parseDouble(tokens[4]);
        }
        return new VDWType(atomType, radius, wellDepth, reductionFactor);
      } catch (NumberFormatException e) {
        String message = "Exception parsing VDW type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * Construct a 1-4 VDWType from multiple input lines.
   *
   * @param input  The overall input String.
   * @param tokens The input String tokenized.
   * @return a VDWType instance.
   */
  public static VDWType parseVDW14(String input, String[] tokens) {
    if (tokens.length < 4) {
      logger.log(Level.WARNING, "Invalid VDW type:\n{0}", input);
    } else {
      try {
        int atomType = parseInt(tokens[1]);
        double radius = parseDouble(tokens[2]);
        double wellDepth = parseDouble(tokens[3]);
        double reductionFactor = -1.0;
        if (tokens.length == 5) {
          reductionFactor = parseDouble(tokens[4]);
        }
        return new VDWType(atomType, radius, wellDepth, reductionFactor, VDWMode.VDW14);
      } catch (NumberFormatException e) {
        String message = "Exception parsing VDW14 type:\n" + input + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
    return null;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compare(String s1, String s2) {
    int t1 = parseInt(s1);
    int t2 = parseInt(s2);
    return Integer.compare(t1, t2);
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
    VDWType vdwType = (VDWType) o;
    return (vdwType.atomClass == this.atomClass);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int hashCode() {
    return Objects.hash(atomClass);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Nicely formatted van der Waals type.
   */
  @Override
  public String toString() {
    StringBuilder vdwString = new StringBuilder("vdw");
    if (vdwMode == VDWMode.VDW14) {
      vdwString.append("14");
    }

    // No reduction factor.
    if (reductionFactor <= 0e0) {
      vdwString.append(format("  %5d  %11.9f  %11.9f", atomClass, radius, wellDepth));
    } else {
      vdwString.append(
          format("  %5d  %11.9f  %11.9f  %5.3f", atomClass, radius, wellDepth, reductionFactor));
    }

    return vdwString.toString();
  }

  /**
   * Create an AmoebaVdwForce force node.
   *
   * @param doc        the Document instance.
   * @param forceField the ForceField to collect constants from.
   * @return The AmoebaVdwForce force node.
   */
  public static Element getXMLForce(Document doc, ForceField forceField) {
    Map<String, VDWType> vdwTypes = forceField.getVDWTypes();
    Map<String, VDWPairType> vdwPairTypes = forceField.getVDWPairTypes();
    if (!vdwTypes.values().isEmpty() || !vdwPairTypes.values().isEmpty()) {
      Element node = doc.createElement("AmoebaVdwForce");
      node.setAttribute("type", forceField.getString("vdwtype", DEFAULT_VDW_TYPE.toString()));
      node.setAttribute("radiusrule", forceField.getString("radiusrule", DEFAULT_RADIUS_RULE.toString()));
      node.setAttribute("radiustype", forceField.getString("radiustype", DEFAULT_RADIUS_TYPE.toString()));
      node.setAttribute("radiussize", forceField.getString("radiussize", DEFAULT_RADIUS_SIZE.toString()));
      node.setAttribute("epsilonrule", forceField.getString("epsilonrule", DEFAULT_EPSILON_RULE.toString()));
      // There is not currently support for a 1-2 scale in OpenMM (this is always zero).
      node.setAttribute("vdw-13-scale", valueOf(forceField.getDouble("vdw-13-scale", DEFAULT_VDW_13_SCALE)));
      node.setAttribute("vdw-14-scale", valueOf(forceField.getDouble("vdw-14-scale", DEFAULT_VDW_14_SCALE)));
      node.setAttribute("vdw-15-scale", valueOf(forceField.getDouble("vdw-15-scale", 1.0)));
      for (VDWType vdwType : vdwTypes.values()) {
        node.appendChild(vdwType.toXML(doc));
      }
      for (VDWPairType vdwPairType : vdwPairTypes.values()) {
        node.appendChild(vdwPairType.toXML(doc));
      }
      return node;
    }
    return null;
  }

  /**
   * Write VDWType to OpenMM XML format.
   */
  public Element toXML(Document doc) {
    Element node = doc.createElement("Vdw");
    node.setAttribute("class", format("%d", atomClass));
    // Convert Angstroms to nm
    node.setAttribute("sigma", format("%f", radius * ANG_TO_NM));
    // Convert Kcal/mol to KJ/mol
    node.setAttribute("epsilon", format("%f", wellDepth * KCAL_TO_KJ));
    if (reductionFactor <= 0e0) {
      node.setAttribute("reduction", "1.0");
    } else {
      node.setAttribute("reduction", format("%f", reductionFactor));
    }
    return node;
  }

  /**
   * Increment the atom class by a specified amount.
   *
   * @param increment The increment to add to the atom class.
   */
  void incrementClass(int increment) {
    atomClass += increment;
    setKey(Integer.toString(atomClass));
  }

  /**
   * Torsion modes include Normal or In-Plane
   */
  public enum VDWMode {
    NORMAL, VDW14
  }

  /**
   * VDW Type.
   */
  public enum VDW_TYPE {
    BUFFERED_14_7,
    LENNARD_JONES;

    public String toString() {
      return name().replace("_", "-");
    }
  }

  /**
   * Radius combining rule.
   */
  public enum RADIUS_RULE {
    ARITHMETIC,
    CUBIC_MEAN,
    GEOMETRIC;

    public String toString() {
      return name().replace("_", "-");
    }
  }

  /**
   * Radius size in the parameter file (radius or diameter).
   */
  public enum RADIUS_SIZE {
    DIAMETER,
    RADIUS
  }

  /**
   * Radius type in the parameter file (R-Min or Sigma).
   */
  public enum RADIUS_TYPE {
    R_MIN,
    SIGMA;

    public String toString() {
      return name().replace("_", "-");
    }
  }

  /**
   * Epsilon combining rule.
   */
  public enum EPSILON_RULE {
    GEOMETRIC,
    HHG,
    W_H;

    public String toString() {
      return name().replace("_", "-");
    }
  }
}
