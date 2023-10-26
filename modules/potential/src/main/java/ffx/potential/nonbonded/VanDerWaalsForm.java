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
package ffx.potential.nonbonded;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWPairType;
import ffx.potential.parameters.VDWType;
import ffx.utilities.FFXKeyword;

import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.toEnumForm;
import static ffx.utilities.KeywordGroup.VanDerWaalsFunctionalForm;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * This class contains fields and methods for maintaining details of the van der Waals functional
 * form.
 *
 * @author Michael J. Schnieders
 */
public class VanDerWaalsForm {

  /**
   * The logger.
   */
  private static final Logger logger = Logger.getLogger(VanDerWaalsForm.class.getName());

  /**
   * Constant <code>RADMIN=0</code>
   */
  static final byte RADMIN = 0;
  /**
   * Constant <code>EPS=1</code>
   */
  static final byte EPS = 1;

  /**
   * The default gamma parameter in Halgren’s buffered 14-7 vdw potential energy functional form.
   */
  private static final double DEFAULT_GAMMA = 0.12;

  /**
   * The default delta parameter in Halgren’s buffered 14-7 vdw potential energy functional form.
   */
  private static final double DEFAULT_DELTA = 0.07;

  /**
   * The default epsilon combining rule.
   */
  private static final EPSILON_RULE DEFAULT_EPSILON_RULE = EPSILON_RULE.GEOMETRIC;

  /**
   * The default radius combining rule.
   */
  private static final RADIUS_RULE DEFAULT_RADIUS_RULE = RADIUS_RULE.ARITHMETIC;

  /**
   * The default radius size.
   */
  private static final RADIUS_SIZE DEFAULT_RADIUS_SIZE = RADIUS_SIZE.RADIUS;

  /**
   * The default radius type.
   */
  private static final RADIUS_TYPE DEFAULT_RADIUS_TYPE = RADIUS_TYPE.R_MIN;

  /**
   * The default van der Waals functional form type.
   */
  private static final VDW_TYPE DEFAULT_VDW_TYPE = VDW_TYPE.LENNARD_JONES;

  /**
   * The default van der Waals scale factor for 1-2 (bonded) interactions.
   */
  private static final double DEFAULT_VDW_12_SCALE = 0.0;

  /**
   * The default van der Waals scale factor for 1-3 (angle) interactions.
   */
  private static final double DEFAULT_VDW_13_SCALE = 0.0;

  /**
   * The default van der Waals scale factor for 1-4 (torisonal) interactions.
   */
  private static final double DEFAULT_VDW_14_SCALE = 1.0;

  /**
   * The default van der Waals taper location is at 90% of the cut-off distance.
   */
  public static final double DEFAULT_VDW_TAPER = 0.9;

  /**
   * The default van der Waals cut-off radius is 12.0 Angstroms.
   */
  public static final double DEFAULT_VDW_CUTOFF = 12.0;

  /**
   * First constant suggested by Halgren for the Buffered-14-7 potential.
   */
  @FFXKeyword(name = "halgren-gamma", keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "0.12",
      description = """
          Sets the value of the gamma parameter in Halgren’s buffered 14-7 vdw potential energy functional form.
          In the absence of the gamma-halgren property, a default value of 0.12 is used."
          """)

  public final double gamma;

  /**
   * Second constant suggested by Halgren for the Buffered-14-7 potential.
   */
  @FFXKeyword(name = "halgren-delta", keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "0.07",
      description = """
          Sets the value of the delta parameter in Halgren’s buffered 14-7 vdw potential energy functional form.
          In the absence of the delta-halgren property, a default value of 0.07 is used.
          """)
  public final double delta;

  /**
   * van der Waals functional form.
   */
  @FFXKeyword(name = "vdwtype", clazz = String.class, keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "LENNARD-JONES",
      description = """
          [LENNARD-JONES / BUFFERED-14-7]
          Sets the functional form for the van der Waals potential energy term.
          The text modifier gives the name of the functional form to be used.
          The default in the absence of the vdwtype keyword is to use the standard two parameter Lennard-Jones function.
          """)
  public VDW_TYPE vdwType;

  /**
   * Epsilon combining rule.
   */
  @FFXKeyword(name = "epsilonrule", clazz = String.class, keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "GEOMETRIC",
      description = """
          [GEOMETRIC / HHG / W-H]
          Selects the combining rule used to derive the epsilon value for van der Waals interactions.
          The default in the absence of the epsilonrule keyword is to use the geometric mean of the
          individual epsilon values of the two atoms involved in the van der Waals interaction.
          """)
  public EPSILON_RULE epsilonRule;

  /**
   * Radius combining rule.
   */
  @FFXKeyword(name = "radiusrule", clazz = String.class, keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "ARITHMETIC",
      description = """
          [ARITHMETIC / GEOMETRIC / CUBIC-MEAN]
          Sets the functional form of the radius combining rule for heteroatomic van der Waals potential energy interactions.
          The default in the absence of the radiusrule keyword is to use the arithmetic mean combining rule to get radii for heteroatomic interactions.
          """)
  public RADIUS_RULE radiusRule;

  /**
   * Radius size in the parameter file (radius or diameter).
   */
  @FFXKeyword(name = "radiussize", clazz = String.class, keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "RADIUS",
      description = """
          [RADIUS / DIAMETER]
          Determines whether the atom size values given in van der Waals parameters read from
          VDW keyword statements are interpreted as atomic radius or diameter values.
          The default in the absence of the radiussize keyword is to assume that vdw size parameters are given as radius values.
          """)
  public RADIUS_SIZE radiusSize;

  /**
   * Radius type in the parameter file (R-Min or Sigma).
   */
  @FFXKeyword(name = "radiustype", clazz = String.class, keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "R-MIN",
      description = """
          [R-MIN / SIGMA]
          Determines whether atom size values given in van der Waals parameters read from VDW keyword
          statements are interpreted as potential minimum (Rmin) or LJ-style sigma values.
          The default in the absence of the radiustype keyword is to assume that vdw size parameters are given as Rmin values.
          """)
  public RADIUS_TYPE radiusType;

  /**
   * Define scale factors between 1-2 atoms.
   */
  @FFXKeyword(name = "vdw-12-scale", keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to van der Waals potential
          interactions between 1-2 connected atoms, i.e., atoms that are directly bonded.
          The default value of 0.0 is used to omit 1-2 interactions,
          if the vdw-12-scale property is not given in either the parameter file or the property file.
          """)
  protected double scale12;

  /**
   * Define scale factors between 1-3 atoms.
   */
  @FFXKeyword(name = "vdw-13-scale", keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to van der Waals potential
          interactions between 1-3 connected atoms, i.e., atoms separated by two covalent bonds.
          The default value of 0.0 is used to omit 1-3 interactions, if the vdw-13-scale property
          is not given in either the parameter file or the property file.
          """)
  protected double scale13;

  /**
   * Define scale factors between 1-4 atoms.
   */
  @FFXKeyword(name = "vdw-14-scale", keywordGroup = VanDerWaalsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to van der Waals potential
          interactions between 1-4 connected atoms, i.e., atoms separated by three covalent bonds.
          The default value of 1.0 is used, if the vdw-14-scale keyword is not given in either
          the parameter file or the property file.
          """)
  protected double scale14;

  /**
   * Maximum number of classes in the force field.
   */
  int maxClass;
  /**
   * Store the Rmin for each class.
   */
  double[] rMin;
  /**
   * Store the Eps for each class
   */
  double[] eps;
  /**
   * Store combined radius and epsilon values.
   */
  private final double[][] radEps;
  /**
   * Store combined radius and epsilon values for 1-4 interactions.
   *
   * <p>Currently this is specific to the CHARMM force fields.
   */
  private final double[][] radEps14;
  /**
   * Softcore vdW partial derivatives that are different between LJ and B14-7.
   */
  private final VDWPowers vdwPowers;
  /**
   * vdW Dispersive Power (e.g. 6).
   */
  final int dispersivePower;
  /**
   * vdW Dispersive Power minus 1 (e.g. 6-1 = 5).
   */
  private final int dispersivePower1;
  /**
   * Define some handy constants.
   */
  final double t1n;
  /**
   * Repulsive power minus dispersive power (e.g., 12-6 = 6).
   */
  final int repDispPower;
  /**
   * Repulsive power minus dispersive power minus 1 (e.g., 12-6-1 = 5).
   */
  private final int repDispPower1;
  /**
   * First constant suggested by Halgren for the Buffered-14-7 potential plus 1
   */
  final double gamma1;

  /**
   * Constructor for VanDerWaalsForm.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public VanDerWaalsForm(ForceField forceField) {

    // Define functional form.
    vdwType = DEFAULT_VDW_TYPE;
    String value = forceField.getString("VDWTYPE", DEFAULT_VDW_TYPE.name());
    try {
      vdwType = VDW_TYPE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized VDWTYPE %s; defaulting to %s.", value, vdwType));
    }

    switch (vdwType) {
      case BUFFERED_14_7 -> vdwPowers = new Buffered_14_7();
      case LENNARD_JONES -> vdwPowers = new LJ_6_12();
      default -> vdwPowers = new VDWPowers();
    }

    // Define epsilon combining rule.
    epsilonRule = DEFAULT_EPSILON_RULE;
    value = forceField.getString("EPSILONRULE", DEFAULT_EPSILON_RULE.name());
    try {
      epsilonRule = EPSILON_RULE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized EPSILONRULE %s; defaulting to %s.", value, epsilonRule));
    }

    // Define radius combining rule.
    radiusRule = DEFAULT_RADIUS_RULE;
    value = forceField.getString("RADIUSRULE", DEFAULT_RADIUS_RULE.name());
    try {
      radiusRule = RADIUS_RULE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized RADIUSRULE %s; defaulting to %s.", value, radiusRule));
    }

    // Define radius size.
    radiusSize = DEFAULT_RADIUS_SIZE;
    value = forceField.getString("RADIUSSIZE", DEFAULT_RADIUS_SIZE.name());
    try {
      radiusSize = RADIUS_SIZE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized RADIUSSIZE %s; defaulting to %s.", value, radiusSize));
    }

    // Define radius type.
    radiusType = DEFAULT_RADIUS_TYPE;
    value = forceField.getString("RADIUSTYPE", DEFAULT_RADIUS_TYPE.name());
    try {
      radiusType = RADIUS_TYPE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized RADIUSTYPE %s; defaulting to %s", value, radiusType));
    }

    // Configure van der Waals well shape parameters.
    int repulsivePower;
    switch (vdwType) {
      case LENNARD_JONES -> {
        repulsivePower = 12;
        dispersivePower = 6;
        delta = 0.0;
        gamma = 0.0;
      }
      default -> {
        repulsivePower = 14;
        dispersivePower = 7;
        delta = forceField.getDouble("DELTA-HALGREN", DEFAULT_DELTA);
        gamma = forceField.getDouble("GAMMA-HALGREN", DEFAULT_GAMMA);
      }
    }

    repDispPower = repulsivePower - dispersivePower;
    dispersivePower1 = dispersivePower - 1;
    repDispPower1 = repDispPower - 1;
    double delta1 = 1.0 + delta;
    t1n = pow(delta1, dispersivePower);
    gamma1 = 1.0 + gamma;

    scale12 = forceField.getDouble("VDW_12_SCALE", DEFAULT_VDW_12_SCALE);
    scale13 = forceField.getDouble("VDW_13_SCALE", DEFAULT_VDW_13_SCALE);
    scale14 = forceField.getDouble("VDW_14_SCALE", DEFAULT_VDW_14_SCALE);
    double scale15 = forceField.getDouble("VDW_15_SCALE", 1.0);
    if (scale15 != 1.0) {
      logger.severe(" Van Der Waals 1-5 masking rules are not supported.");
    }

    // The convention in TINKER is a vdw-14-scale factor of 2.0 means to scale by 0.5.
    if (scale12 > 1.0) {
      scale12 = 1.0 / scale12;
    }
    if (scale13 > 1.0) {
      scale13 = 1.0 / scale13;
    }
    if (scale14 > 1.0) {
      scale14 = 1.0 / scale14;
    }

    Map<String, VDWType> map = forceField.getVDWTypes();
    TreeMap<String, VDWType> vdwTypes = new TreeMap<>(map);
    maxClass = 0;
    for (VDWType currentType : vdwTypes.values()) {
      if (currentType.atomClass > maxClass) {
        maxClass = currentType.atomClass;
      }
    }
    radEps = new double[maxClass + 1][2 * (maxClass + 1)];
    radEps14 = new double[maxClass + 1][2 * (maxClass + 1)];

    // Scale factor to convert to vdW size to Rmin.
    double radScale = switch (radiusSize) {
      case DIAMETER -> 0.5;
      default -> 1.0;
    };
    switch (radiusType) {
      case SIGMA -> radScale *= 1.122462048309372981;
      default -> {
      }
    }

    rMin = new double[maxClass + 1];
    eps = new double[maxClass + 1];

    // Atom Class numbering starts at 1.
    for (VDWType vdwi : vdwTypes.values()) {
      int i = vdwi.atomClass;
      double ri = radScale * vdwi.radius;
      double e1 = vdwi.wellDepth;
      rMin[i] = ri;
      eps[i] = e1;
      for (VDWType vdwj : vdwTypes.tailMap(vdwi.getKey()).values()) {
        int j = vdwj.atomClass;
        double rj = radScale * vdwj.radius;
        double e2 = vdwj.wellDepth;
        double radmin = getCombinedRadius(ri, rj, radiusRule);
        double eps = getCombinedEps(e1, e2, ri, rj, epsilonRule);
        if (radmin > 0) {
          radEps[i][j * 2 + RADMIN] = 1.0 / radmin;
          radEps[j][i * 2 + RADMIN] = 1.0 / radmin;
          radEps14[i][j * 2 + RADMIN] = 1.0 / radmin;
          radEps14[j][i * 2 + RADMIN] = 1.0 / radmin;
        } else {
          radEps[i][j * 2 + RADMIN] = 0.0;
          radEps[j][i * 2 + RADMIN] = 0.0;
          radEps14[i][j * 2 + RADMIN] = 0.0;
          radEps14[j][i * 2 + RADMIN] = 0.0;
        }
        radEps[i][j * 2 + EPS] = eps;
        radEps[j][i * 2 + EPS] = eps;
        radEps14[i][j * 2 + EPS] = eps;
        radEps14[j][i * 2 + EPS] = eps;
      }
    }

    // Handle vdw14 types -- loop over VDW types.
    Map<String, VDWType> vdw14Types = forceField.getVDW14Types();
    for (VDWType vdwi : vdwTypes.values()) {
      // Replace a normal VDW type with a VDW14 type if available.
      VDWType vdw14 = forceField.getVDW14Type(vdwi.getKey());
      if (vdw14 != null) {
        vdwi = vdw14;
      }
      int i = vdwi.atomClass;
      double ri = radScale * vdwi.radius;
      double e1 = vdwi.wellDepth;
      // Loop over VDW14 types.
      for (VDWType vdwj : vdw14Types.values()) {
        int j = vdwj.atomClass;
        double rj = radScale * vdwj.radius;
        double e2 = vdwj.wellDepth;
        double radmin = getCombinedRadius(ri, rj, radiusRule);
        double eps = getCombinedEps(e1, e2, ri, rj, epsilonRule);
        if (radmin > 0) {
          radEps14[i][j * 2 + RADMIN] = 1.0 / radmin;
          radEps14[j][i * 2 + RADMIN] = 1.0 / radmin;
        } else {
          radEps14[i][j * 2 + RADMIN] = 0.0;
          radEps14[j][i * 2 + RADMIN] = 0.0;
        }
        radEps14[i][j * 2 + EPS] = eps;
        radEps14[j][i * 2 + EPS] = eps;
      }
    }

    // Replace combined VDW and VDW14 parameters with VDW Pair values.
    Map<String, VDWPairType> vdwPairTypes = forceField.getVDWPairTypes();
    for (VDWPairType vdwPairType : vdwPairTypes.values()) {
      int i = vdwPairType.atomClasses[0];
      int j = vdwPairType.atomClasses[1];
      double radmin = vdwPairType.radius;
      double eps = vdwPairType.wellDepth;
      if (radmin > 0) {
        radEps[i][j * 2 + RADMIN] = 1.0 / radmin;
        radEps[j][i * 2 + RADMIN] = 1.0 / radmin;
        radEps14[i][j * 2 + RADMIN] = 1.0 / radmin;
        radEps14[j][i * 2 + RADMIN] = 1.0 / radmin;
      } else {
        radEps[i][j * 2 + RADMIN] = 0.0;
        radEps[j][i * 2 + RADMIN] = 0.0;
        radEps14[i][j * 2 + RADMIN] = 0.0;
        radEps14[j][i * 2 + RADMIN] = 0.0;
      }
      radEps[i][j * 2 + EPS] = eps;
      radEps[j][i * 2 + EPS] = eps;
      radEps14[i][j * 2 + EPS] = eps;
      radEps14[j][i * 2 + EPS] = eps;
    }
  }

  /**
   * Get the combined EPS value.
   *
   * @param ei          The eps value of the first atom.
   * @param ej          The eps value of the second atom.
   * @param ri          The radius of the first atom.
   * @param rj          The radius of the second atom.
   * @param epsilonRule The epsilon rule to use.
   * @return The combined eps value.
   */
  public static double getCombinedEps(double ei, double ej, double ri, double rj,
                                      EPSILON_RULE epsilonRule) {
    double sei = sqrt(ei);
    double sej = sqrt(ej);
    switch (epsilonRule) {
      case GEOMETRIC -> {
        return sei * sej;
      }
      case W_H -> {
        // ep = 2.0d0 * (seps(i)*seps(k)) * (rad(i)*rad(k))**3 / (rad(i)**6+rad(k)**6)
        double rr = ri * rj;
        double rr3 = rr * rr * rr;
        double ri6 = pow(ri, 6);
        double rj6 = pow(rj, 6);
        return 2.0 * sei * sej * rr3 / (ri6 + rj6);
      }
      default -> {
        return 4.0 * (ei * ej) / ((sei + sej) * (sei + sej));
      }
    }
  }

  /**
   * Get the combined radius value.
   *
   * @param ri         The radius of the first atom.
   * @param rj         The radius of the second atom.
   * @param radiusRule The radius combining rule to use.
   * @return The combined radius value.
   */
  public static double getCombinedRadius(double ri, double rj, RADIUS_RULE radiusRule) {
    switch (radiusRule) {
      case ARITHMETIC -> {
        return ri + rj;
      }
      case GEOMETRIC -> {
        return 2.0 * sqrt(ri) * sqrt(rj);
      }
      default -> {
        double ri2 = ri * ri;
        double ri3 = ri * ri2;
        double rj2 = rj * rj;
        double rj3 = rj * rj2;
        return 2.0 * (ri3 + rj3) / (ri2 + rj2);
      }
    }
  }

  /**
   * Return the combined well depth (kcal/mol)
   *
   * @param class1 Class for atom 1.
   * @param class2 Class for atom 2.
   * @return Combined Eps.
   */
  public double getCombinedEps(int class1, int class2) {
    return radEps[class1][class2 * 2 + EPS];
  }

  /**
   * Return the combined well depth (kcal/mol) for special 1-4 interactions
   *
   * @param class1 Class for atom 1.
   * @param class2 Class for atom 2.
   * @return Combined Eps.
   */
  public double getCombinedEps14(int class1, int class2) {
    return radEps14[class1][class2 * 2 + EPS];
  }

  /**
   * Return the combined inverse Rmin value (1/Rmin).
   *
   * @param class1 Class for atom 1.
   * @param class2 Class for atom 2.
   * @return Combined inverse Rmin.
   */
  public double getCombinedInverseRmin(int class1, int class2) {
    return radEps[class1][class2 * 2 + RADMIN];
  }

  /**
   * Return the combined inverse Rmin value (1/Rmin) for special 1-4 interactions.
   *
   * @param class1 Class for atom 1.
   * @param class2 Class for atom 2.
   * @return Combined inverse Rmin.
   */
  public double getCombinedInverseRmin14(int class1, int class2) {
    return radEps14[class1][class2 * 2 + RADMIN];
  }

  /**
   * Return the eps value for each class.
   *
   * @return Returns the eps array.
   */
  public double[] getEps() {
    return eps;
  }

  /**
   * Return the Rmin value for each class.
   *
   * @return Returns the Rmin array.
   */
  public double[] getRmin() {
    return rMin;
  }

  /**
   * Getter for the field <code>scale14</code>.
   *
   * @return a double.
   */
  public double getScale14() {
    return scale14;
  }

  /**
   * rhoDisp1.
   *
   * @param rho a double.
   * @return a double.
   */
  double rhoDisp1(double rho) {
    return vdwPowers.rhoDisp1(rho);
  }

  /**
   * rhoDelta1.
   *
   * @param rhoDelta a double.
   * @return a double.
   */
  double rhoDelta1(double rhoDelta) {
    return vdwPowers.rhoDisp1(rhoDelta);
  }

  /**
   * VDW Type.
   */
  public enum VDW_TYPE {
    BUFFERED_14_7,
    LENNARD_JONES
  }

  /**
   * Radius combining rule.
   */
  public enum RADIUS_RULE {
    ARITHMETIC,
    CUBIC_MEAN,
    GEOMETRIC
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
    SIGMA
  }

  /**
   * Epsilon combining rule.
   */
  public enum EPSILON_RULE {
    GEOMETRIC,
    HHG,
    W_H
  }

  /**
   * Softcore vdW partial derivatives that are different between LJ and B14-7.
   */
  private class VDWPowers {

    public double rhoDelta1(double rhoDelta) {
      return pow(rhoDelta, repDispPower1);
    }

    public double rhoDisp1(double rho) {
      return pow(rho, dispersivePower1);
    }
  }

  /**
   * LJ softcore vdW partial derivatives.
   */
  private class LJ_6_12 extends VDWPowers {

    @Override
    public double rhoDelta1(double rhoDelta) {
      double rhoDelta2 = rhoDelta * rhoDelta;
      return rhoDelta2 * rhoDelta2 * rhoDelta;
    }

    @Override
    public double rhoDisp1(double rho) {
      double rho2 = rho * rho;
      return rho2 * rho2 * rho;
    }
  }

  /**
   * Buffered 14-7 softcore vdW partial derivatives.
   */
  private class Buffered_14_7 extends VDWPowers {

    @Override
    public double rhoDelta1(double rhoDelta) {
      double rhoDelta2 = rhoDelta * rhoDelta;
      return rhoDelta2 * rhoDelta2 * rhoDelta2;
    }

    @Override
    public double rhoDisp1(double rho) {
      double rho2 = rho * rho;
      return rho2 * rho2 * rho2;
    }
  }
}
