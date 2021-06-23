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
package ffx.potential.nonbonded;

import static ffx.potential.parameters.ForceField.toEnumForm;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWType;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 * This class contains fields and methods for maintaining details of the van der Waals functional
 * form.
 *
 * @author Michael J. Schnieders
 */
public class VanDerWaalsForm {

  /** Constant <code>RADMIN=0</code> */
  static final byte RADMIN = 0;
  /** Constant <code>EPS=1</code> */
  static final byte EPS = 1;
  /** The logger. */
  private static final Logger logger = Logger.getLogger(VanDerWaalsForm.class.getName());
  /** First constant suggested by Halgren for the Buffered-14-7 potential. */
  public final double gamma;
  /** Second constant suggested by Halgren for the Buffered-14-7 potential. */
  public final double delta;
  /** vdW Dispersive Power (e.g. 6). */
  final int dispersivePower;

  final double gamma1;
  /** Define some handy constants. */
  final double t1n;

  final int repDispPower;
  private final VDWPowers vdwPowers;
  private final int dispersivePower1;
  private final int repDispPower1;
  /** van der Waals functional form. */
  public VDW_TYPE vdwType;

  public EPSILON_RULE epsilonRule;
  public RADIUS_RULE radiusRule;
  public RADIUS_SIZE radiusSize;
  public RADIUS_TYPE radiusType;
  /** Define scale factors between 1-2, 1-3, etc. atoms. */
  protected double scale12;

  protected double scale13;
  double scale14;
  /** Maximum number of classes in the force field. */
  int maxClass;
  /** Store the Rmin for each class. */
  double[] rMin;
  /** Store the Eps for each class */
  double[] eps;
  /** Store combined radius and epsilon values. */
  private final double[][] radEps;
  /**
   * Store combined radius and epsilon values for 1-4 interactions.
   *
   * <p>Currently this is specific to the CHARMM force fields.
   */
  private final double[][] radEps14;
  /**
   * Constructor for VanDerWaalsForm.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public VanDerWaalsForm(ForceField forceField) {

    // Set-up default rules.
    vdwType = VDW_TYPE.BUFFERED_14_7;
    epsilonRule = EPSILON_RULE.HHG;
    radiusRule = RADIUS_RULE.CUBIC_MEAN;
    radiusSize = RADIUS_SIZE.DIAMETER;
    radiusType = RADIUS_TYPE.R_MIN;

    // Define functional form.
    String value = forceField.getString("VDWTYPE", vdwType.toString());
    try {
      vdwType = VDW_TYPE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized VDWTYPE %s; defaulting to %s.", value, vdwType));
    }

    switch (vdwType) {
      case BUFFERED_14_7:
        vdwPowers = new Buffered_14_7();
        break;
      case LENNARD_JONES:
        vdwPowers = new LJ_6_12();
        break;
      default:
        vdwPowers = new VDWPowers();
        break;
    }

    // Define epsilon combining rule.
    value = forceField.getString("EPSILONRULE", epsilonRule.toString());
    try {
      epsilonRule = EPSILON_RULE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized EPSILONRULE %s; defaulting to %s.", value, epsilonRule));
    }

    // Define radius combining rule.
    value = forceField.getString("RADIUSRULE", radiusRule.toString());
    try {
      radiusRule = RADIUS_RULE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized RADIUSRULE %s; defaulting to %s.", value, radiusRule));
    }

    // Define radius size.
    value = forceField.getString("RADIUSSIZE", radiusSize.toString());
    try {
      radiusSize = RADIUS_SIZE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized RADIUSSIZE %s; defaulting to %s.", value, radiusSize));
    }

    // Define radius type.
    value = forceField.getString("RADIUSTYPE", radiusType.toString());
    try {
      radiusType = RADIUS_TYPE.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized RADIUSTYPE %s; defaulting to %s", value, radiusType));
    }

    // Configure van der Waals well shape parameters.
    int repulsivePower;
    switch (vdwType) {
      case LENNARD_JONES:
        repulsivePower = 12;
        dispersivePower = 6;
        delta = 0.0;
        gamma = 0.0;
        break;
      case BUFFERED_14_7:
      default:
        repulsivePower = 14;
        dispersivePower = 7;
        delta = 0.07;
        gamma = 0.12;
        break;
    }

    repDispPower = repulsivePower - dispersivePower;
    dispersivePower1 = dispersivePower - 1;
    repDispPower1 = repDispPower - 1;
    double delta1 = 1.0 + delta;
    t1n = pow(delta1, dispersivePower);
    gamma1 = 1.0 + gamma;

    scale12 = forceField.getDouble("VDW_12_SCALE", 0.0);
    scale13 = forceField.getDouble("VDW_13_SCALE", 0.0);
    scale14 = forceField.getDouble("VDW_14_SCALE", 1.0);
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

    Map<String, VDWType> map14 = forceField.getVDW14Types();
    TreeMap<String, VDWType> vdw14Types = new TreeMap<>(map14);

    // Scale factor to convert to vdW size to Rmin.
    double radScale;
    switch (radiusSize) {
      case DIAMETER:
        radScale = 0.5;
        break;
      case RADIUS:
      default:
        radScale = 1.0;
        break;
    }
    switch (radiusType) {
      case SIGMA:
        radScale *= 1.122462048309372981;
        break;
      case R_MIN:
      default:
        break;
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
        double eps = getCombinedEps(e1, e2, epsilonRule);
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
        double eps = getCombinedEps(e1, e2, epsilonRule);
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
  }

  /**
   * Get the combined EPS value.
   *
   * @param ei The eps value of the first atom.
   * @param ej The eps value of the second atom.
   * @param epsilonRule The epsilon rule to use.
   * @return The combined eps value.
   */
  public static double getCombinedEps(double ei, double ej, EPSILON_RULE epsilonRule) {
    double sei = sqrt(ei);
    double sej = sqrt(ej);
    switch (epsilonRule) {
      case GEOMETRIC:
        return sei * sej;
      default:
      case HHG:
        return 4.0 * (ei * ej) / ((sei + sej) * (sei + sej));
    }
  }

  /**
   * Get the combined radius value.
   *
   * @param ri The radius of the first atom.
   * @param rj The radius of the second atom.
   * @param radiusRule The radius combining rule to use.
   * @return The combined radius value.
   */
  public static double getCombinedRadius(double ri, double rj, RADIUS_RULE radiusRule) {
    switch (radiusRule) {
      case ARITHMETIC:
        return ri + rj;
      case GEOMETRIC:
        return 2.0 * sqrt(ri) * sqrt(rj);
      default:
      case CUBIC_MEAN:
        double ri2 = ri * ri;
        double ri3 = ri * ri2;
        double rj2 = rj * rj;
        double rj3 = rj * rj2;
        return 2.0 * (ri3 + rj3) / (ri2 + rj2);
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

  public enum VDW_TYPE {
    BUFFERED_14_7,
    LENNARD_JONES
  }

  public enum RADIUS_RULE {
    ARITHMETIC,
    CUBIC_MEAN,
    GEOMETRIC
  }

  public enum RADIUS_SIZE {
    DIAMETER,
    RADIUS
  }

  public enum RADIUS_TYPE {
    R_MIN,
    SIGMA
  }

  public enum EPSILON_RULE {
    GEOMETRIC,
    HHG
  }

  private class VDWPowers {

    public double rhoDelta1(double rhoDelta) {
      return pow(rhoDelta, repDispPower1);
    }

    public double rhoDisp1(double rho) {
      return pow(rho, dispersivePower1);
    }
  }

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
