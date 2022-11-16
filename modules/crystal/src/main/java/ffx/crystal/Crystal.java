// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
package ffx.crystal;

import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.DoubleMath.length;
import static ffx.numerics.math.MatrixMath.mat3Mat3;
import static ffx.numerics.math.MatrixMath.mat3SymVec6;
import static ffx.numerics.math.MatrixMath.transpose3;
import static ffx.numerics.math.ScalarMath.mod;
import static ffx.utilities.Constants.AVOGADRO;
import static ffx.utilities.KeywordGroup.UnitCellAndSpaceGroup;
import static ffx.utilities.StringUtils.padRight;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.signum;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;
import static org.apache.commons.math3.util.FastMath.toRadians;

import ffx.utilities.FFXKeyword;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * The Crystal class encapsulates the lattice parameters and space group that describe the geometry
 * and symmetry of a crystal. Methods are available to apply the minimum image convention and space
 * group symmetry operators.
 *
 * @author Michael J. Schnieders
 * @see ReplicatesCrystal
 * @since 1.0
 */
public class Crystal {

  private static final Logger logger = Logger.getLogger(Crystal.class.getName());

  /** The space group of the crystal. */
  @FFXKeyword(name = "SpaceGroup", keywordGroup = UnitCellAndSpaceGroup, clazz = SpaceGroup.class, defaultValue = "P1",
      description = "This keyword selects the space group to be used in manipulation of crystal unit cells and asymmetric units.")
  public final SpaceGroup spaceGroup;

  /** Length of the cell edge in the direction of the <b>a</b> basis vector. */
  @FFXKeyword(name = "a-axis", keywordGroup = UnitCellAndSpaceGroup, defaultValue = "None",
      description = "Sets the value of the a-axis length for a crystal unit cell, or, equivalently, the X-axis length for a periodic box (Angstroms).")
  public double a;

  /** Length of the cell edge in the direction of the <b>b</b> basis vector. */
  @FFXKeyword(name = "b-axis", keywordGroup = UnitCellAndSpaceGroup, defaultValue = "A-Axis",
      description = "Sets the value of the b-axis length for a crystal unit cell, or, equivalently, the Y-axis length for a periodic box (Angstroms).")
  public double b;

  /** Length of the cell edge in the direction of the <b>c</b> basis vector. */
  @FFXKeyword(name = "c-axis", keywordGroup = UnitCellAndSpaceGroup, defaultValue = "A-Axis",
      description = "Sets the value of the c-axis length for a crystal unit cell, or, equivalently, the Z-axis length for a periodic box (Angstroms).")
  public double c;

  /** The interaxial lattice angle between <b>b</b> and <b>c</b>. */
  @FFXKeyword(name = "alpha", keywordGroup = UnitCellAndSpaceGroup, defaultValue = "90.0",
      description = "Sets the value of the α-angle of a crystal unit cell, i.e., the angle between the b-axis and c-axis of a unit cell, or, equivalently, the angle between the Y-axis and Z-axis of a periodic box.")
  public double alpha;

  /** The interaxial lattice angle between <b>a</b> and <b>c</b>. */
  @FFXKeyword(name = "beta", keywordGroup = UnitCellAndSpaceGroup, defaultValue = "Alpha",
      description = "Sets the value of the β-angle of a crystal unit cell, i.e., the angle between the a-axis and c-axis of a unit cell, or, equivalently, the angle between the X-axis and Z-axis of a periodic box.")
  public double beta;

  /** The interaxial lattice angle between <b>a</b> and <b>b</b>. */
  @FFXKeyword(name = "gamma", keywordGroup = UnitCellAndSpaceGroup, defaultValue = "Alpha",
      description = "Sets the value of the γ-angle of a crystal unit cell, i.e., the angle between the a-axis and b-axis of a unit cell, or, equivalently, the angle between the X-axis and Y-axis of a periodic box.")
  public double gamma;

  /** A mask equal to 0 for X-coordinates. */
  private static final int XX = 0;
  /** A mask equal to 1 for Y-coordinates. */
  private static final int YY = 1;
  /** A mask equal to 2 for Z-coordinates. */
  private static final int ZZ = 2;

  /**
   * Matrix to convert from fractional to Cartesian coordinates.
   * <br>a-axis vector is the first row of A^(-1).
   * <br>b-axis vector is the second row of A^(-1).
   * <br>c-axis vector is the third row of A^(-1).
   */
  public final double[][] Ai = new double[3][3];
  /** The direct space metric matrix. */
  public final double[][] G = new double[3][3];
  /** Reference to the space group crystal system. */
  private final CrystalSystem crystalSystem;
  /** Reference to the space group lattice system. */
  private final LatticeSystem latticeSystem;
  /** The crystal unit cell volume. */
  public double volume;
  /** Matrix to convert from Cartesian to fractional coordinates. */
  public double[][] A;
  /** Entry in the A matrix. */
  public double A00, A01, A02, A10, A11, A12, A20, A21, A22;
  /** Interfacial radius in the direction of the A-axis. */
  public double interfacialRadiusA;
  /** Interfacial radius in the direction of the B-axis. */
  public double interfacialRadiusB;
  /** Interfacial radius in the direction of the C-axis. */
  public double interfacialRadiusC;
  /**
   * Anisotropic bulk solvent B-factor scaling (0 or 1 for each component).
   */
  public int[] scaleB = new int[6];
  /**
   * Number of bulk solvent B-factor components.
   */
  public int scaleN;
  /** Entry in the Ai matrix. */
  public double Ai00, Ai01, Ai02, Ai10, Ai11, Ai12, Ai20, Ai21, Ai22;
  /** Change in the volume with respect to a. */
  public double dVdA;
  /** Change in the volume with respect to b. */
  public double dVdB;
  /** Change in the volume with respect to c. */
  public double dVdC;
  /**
   * Change in the volume with respect to alpha (in Radians). This is set to zero if alpha is fixed.
   */
  public double dVdAlpha;
  /**
   * Change in the volume with respect to beta (in Radians). This is set to zero if beta is fixed.
   */
  public double dVdBeta;
  /**
   * Change in the volume with respect to gamma (in Radians). This is set to zero if gamma is fixed.
   */
  public double dVdGamma;
  /**
   * For some finite-difference calculations, it's currently necessary to remove lattice system
   * restrictions.
   */
  boolean checkRestrictions = true;
  /**
   * An atom and one of its symmetry copies within the specialPositionCutoff should be flagged to be
   * at a special position.
   */
  private double specialPositionCutoff = 0.3;
  /** Copy of symmetry operators in Cartesian coordinates. */
  private List<SymOp> symOpsCartesian;
  /** The reciprocal space metric matrix. */
  private double[][] Gstar;
  /**
   * SpecialPositionCutoff squared.
   */
  private double specialPositionCutoff2 = specialPositionCutoff * specialPositionCutoff;
  /**
   * Flag to indicate an aperiodic system.
   */
  private boolean aperiodic;

  /**
   * The Crystal class encapsulates the lattice parameters and space group. Methods are available to
   * apply the minimum image convention and to apply space group operators.
   *
   * @param a The a-axis length.
   * @param b The b-axis length.
   * @param c The c-axis length.
   * @param alpha The alpha angle.
   * @param beta The beta angle.
   * @param gamma The gamma angle.
   * @param sgNumber The space group number.
   */
  public Crystal(double a, double b, double c, double alpha, double beta, double gamma,
      int sgNumber) {
    this(a, b, c, alpha, beta, gamma, SpaceGroupDefinitions.spaceGroupFactory(sgNumber).pdbName);
  }

  /**
   * The Crystal class encapsulates the lattice parameters and space group. Methods are available to
   * apply the minimum image convention and to apply space group operators.
   *
   * @param a The a-axis length.
   * @param b The b-axis length.
   * @param c The c-axis length.
   * @param alpha The alpha angle.
   * @param beta The beta angle.
   * @param gamma The gamma angle.
   * @param sg The space group symbol.
   */
  public Crystal(double a, double b, double c, double alpha, double beta, double gamma, String sg) {
    // Crystal SpaceGroup and LatticeSystem are final variables. Temp variable to delay assigning.
    SpaceGroup tempSG = SpaceGroupDefinitions.spaceGroupFactory(sg);
    LatticeSystem tempLS = tempSG.latticeSystem;
    this.a = a;
    this.b = b;
    this.c = c;
    this.alpha = alpha;
    this.beta = beta;
    this.gamma = gamma;
    aperiodic = false;
    if (!tempLS.validParameters(a, b, c, alpha, beta, gamma)) {
      // Invalid parameters... Start error/warning log and try to fix.
      StringBuilder sb = new StringBuilder(format(
          " The %s lattice parameters do not satisfy the %s lattice system restrictions.\n",
          tempSG.pdbName, tempLS));
      sb.append(format("  A-axis:                              %18.15e\n", a));
      sb.append(format("  B-axis:                              %18.15e\n", b));
      sb.append(format("  C-axis:                              %18.15e\n", c));
      sb.append(format("  Alpha:                               %18.15e\n", alpha));
      sb.append(format("  Beta:                                %18.15e\n", beta));
      sb.append(format("  Gamma:                               %18.15e\n", gamma));
      logger.info(sb.toString());
      sb = new StringBuilder();
      if (tempLS == LatticeSystem.HEXAGONAL_LATTICE
          || tempLS == LatticeSystem.RHOMBOHEDRAL_LATTICE) {
        // Try to convert between hexagonal and rhombohedral lattices to fix crystal.
        Crystal convertedCrystal = SpaceGroupConversions.hrConversion(a, b, c, alpha, beta, gamma,
            tempSG);
        this.a = convertedCrystal.a;
        this.b = convertedCrystal.b;
        this.c = convertedCrystal.c;
        this.alpha = convertedCrystal.alpha;
        this.beta = convertedCrystal.beta;
        this.gamma = convertedCrystal.gamma;
        spaceGroup = convertedCrystal.spaceGroup;
        crystalSystem = spaceGroup.crystalSystem;
        latticeSystem = spaceGroup.latticeSystem;
        sb.append(" Converted ").append(tempSG.pdbName).append(" to ").append(spaceGroup.pdbName);
        if (!latticeSystem.validParameters(this.a, this.b, this.c, this.alpha, this.beta,
            this.gamma)) {
          sb.append(format(
              " The %s lattice parameters do not satisfy the %s lattice system restrictions.\n",
              spaceGroup.pdbName, latticeSystem));
          sb.append(format("  A-axis:                              %18.15e\n", this.a));
          sb.append(format("  B-axis:                              %18.15e\n", this.b));
          sb.append(format("  C-axis:                              %18.15e\n", this.c));
          sb.append(format("  Alpha:                               %18.15e\n", this.alpha));
          sb.append(format("  Beta:                                %18.15e\n", this.beta));
          sb.append(format("  Gamma:                               %18.15e\n", this.gamma));
          logger.severe(sb.toString());
        } else {
          // Successfully converted space group between hexagonal and rhombohedral.
          logger.info(sb.toString());
        }
      } else {
        // Invalid lattice parameters. Update Crystal as much as possible, then print error message.
        this.a = a;
        this.b = b;
        this.c = c;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        spaceGroup = tempSG;
        crystalSystem = spaceGroup.crystalSystem;
        latticeSystem = spaceGroup.latticeSystem;
        logger.severe(sb.toString());
      }
    } else {
      // Valid parameters, update crystal and continue.
      this.a = a;
      this.b = b;
      this.c = c;
      this.alpha = alpha;
      this.beta = beta;
      this.gamma = gamma;
      spaceGroup = tempSG;
      crystalSystem = spaceGroup.crystalSystem;
      latticeSystem = spaceGroup.latticeSystem;
    }

    for (int i = 0; i < 6; i++) {
      scaleB[i] = -1;
    }

    SymOp symop;
    double[][] rot;
    int index = 0;
    switch (crystalSystem) {
      case TRICLINIC:
        for (int i = 0; i < 6; i++) {
          scaleB[i] = index++;
        }
        break;
      case MONOCLINIC:
        scaleB[0] = index++;
        scaleB[1] = index++;
        scaleB[2] = index++;
        // determine unique axis
        symop = spaceGroup.symOps.get(1);
        rot = symop.rot;
        if (rot[0][0] > 0) {
          scaleB[5] = index++;
        } else if (rot[1][1] > 0) {
          scaleB[4] = index++;
        } else {
          scaleB[3] = index++;
        }
        break;
      case ORTHORHOMBIC:
        scaleB[0] = index++;
        scaleB[1] = index++;
        scaleB[2] = index++;
        break;
      case TETRAGONAL:
        scaleB[0] = index++;
        scaleB[1] = scaleB[0];
        scaleB[2] = index++;
        break;
      case TRIGONAL:
      case HEXAGONAL:
        boolean hexagonal = false;
        for (int i = 1; i < spaceGroup.symOps.size(); i++) {
          symop = spaceGroup.symOps.get(i);
          rot = symop.rot;
          index = 0;
          if ((rot[1][1] * rot[1][2] == -1)
              || (rot[2][1] * rot[2][2] == -1)
              || (rot[1][1] * rot[1][2] == 1)
              || (rot[2][1] * rot[2][2] == 1)) {
            scaleB[0] = index++;
            scaleB[1] = index++;
            scaleB[2] = scaleB[1];
            hexagonal = true;
          } else if ((rot[0][0] * rot[0][2] == -1)
              || (rot[2][0] * rot[2][2] == -1)
              || (rot[0][0] * rot[0][2] == 1)
              || (rot[2][0] * rot[2][2] == 1)) {
            scaleB[0] = index++;
            scaleB[1] = index++;
            scaleB[2] = scaleB[0];
            hexagonal = true;
          } else if ((rot[0][0] * rot[0][1] == -1)
              || (rot[1][0] * rot[1][1] == -1)
              || (rot[0][0] * rot[0][1] == 1)
              || (rot[1][0] * rot[1][1] == 1)) {
            scaleB[0] = index++;
            scaleB[1] = scaleB[0];
            scaleB[2] = index++;
            hexagonal = true;
          }
          if (hexagonal) {
            break;
          }
        }
        if (!hexagonal) {
          // rhombohedral
          scaleB[3] = index++;
          scaleB[4] = scaleB[3];
          scaleB[5] = scaleB[3];
        }
        break;
      case CUBIC:
        break;
    }
    scaleN = index;

    updateCrystal();
  }

  /**
   * checkProperties
   *
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @return a {@link ffx.crystal.Crystal} object.
   */
  public static Crystal checkProperties(CompositeConfiguration properties) {
    double a = properties.getDouble("a-axis", -1.0);
    double b = properties.getDouble("b-axis", -1.0);
    double c = properties.getDouble("c-axis", -1.0);
    double alpha = properties.getDouble("alpha", -1.0);
    double beta = properties.getDouble("beta", -1.0);
    double gamma = properties.getDouble("gamma", -1.0);
    String sg = properties.getString("spacegroup", null);

    sg = SpaceGroupInfo.pdb2ShortName(sg);

    if (a < 0.0 || b < 0.0 || c < 0.0 || alpha < 0.0 || beta < 0.0 || gamma < 0.0 || sg == null) {
      return null;
    }

    // check the space group name is valid
    SpaceGroup spaceGroup = SpaceGroupDefinitions.spaceGroupFactory(sg);
    if (spaceGroup == null) {
      sg = sg.replaceAll(" ", "");
      spaceGroup = SpaceGroupDefinitions.spaceGroupFactory(sg);
      if (spaceGroup == null) {
        return null;
      }
    }

    return new Crystal(a, b, c, alpha, beta, gamma, sg);
  }

  /**
   * invressq
   *
   * @param hkl a {@link ffx.crystal.HKL} object.
   * @return a double.
   */
  public double invressq(HKL hkl) {
    return hkl.quadForm(Gstar);
  }

  /**
   * res
   *
   * @param hkl a {@link ffx.crystal.HKL} object.
   * @return a double.
   */
  public double res(HKL hkl) {
    return 1.0 / sqrt(hkl.quadForm(Gstar));
  }

  /**
   * aperiodic
   *
   * @return a boolean.
   */
  public boolean aperiodic() {
    return aperiodic;
  }

  /**
   * Apply a fractional symmetry operator to an array of Cartesian coordinates. If the arrays x, y or
   * z are null or not of length n, the method returns immediately. If mateX, mateY or mateZ are null
   * or not of length n, new arrays are allocated.
   *
   * @param n Number of atoms.
   * @param x Input fractional x-coordinates.
   * @param y Input fractional y-coordinates.
   * @param z Input fractional z-coordinates.
   * @param mateX Output fractional x-coordinates.
   * @param mateY Output fractional y-coordinates.
   * @param mateZ Output fractional z-coordinates.
   * @param symOp The fractional symmetry operator.
   */
  public void applySymOp(
      int n,
      double[] x,
      double[] y,
      double[] z,
      double[] mateX,
      double[] mateY,
      double[] mateZ,
      SymOp symOp) {
    if (x == null || y == null || z == null) {
      return;
    }
    if (x.length < n || y.length < n || z.length < n) {
      return;
    }
    if (mateX == null || mateX.length < n) {
      mateX = new double[n];
    }
    if (mateY == null || mateY.length < n) {
      mateY = new double[n];
    }
    if (mateZ == null || mateZ.length < n) {
      mateZ = new double[n];
    }

    final double[][] rot = symOp.rot;
    final double[] trans = symOp.tr;

    final double rot00 = rot[0][0];
    final double rot10 = rot[1][0];
    final double rot20 = rot[2][0];
    final double rot01 = rot[0][1];
    final double rot11 = rot[1][1];
    final double rot21 = rot[2][1];
    final double rot02 = rot[0][2];
    final double rot12 = rot[1][2];
    final double rot22 = rot[2][2];
    final double t0 = trans[0];
    final double t1 = trans[1];
    final double t2 = trans[2];
    for (int i = 0; i < n; i++) {
      double xc = x[i];
      double yc = y[i];
      double zc = z[i];
      // Convert to fractional coordinates.
      double xi = xc * A00 + yc * A10 + zc * A20;
      double yi = xc * A01 + yc * A11 + zc * A21;
      double zi = xc * A02 + yc * A12 + zc * A22;
      // Apply Symmetry Operator.
      double fx = rot00 * xi + rot01 * yi + rot02 * zi + t0;
      double fy = rot10 * xi + rot11 * yi + rot12 * zi + t1;
      double fz = rot20 * xi + rot21 * yi + rot22 * zi + t2;
      // Convert back to Cartesian coordinates.
      mateX[i] = fx * Ai00 + fy * Ai10 + fz * Ai20;
      mateY[i] = fx * Ai01 + fy * Ai11 + fz * Ai21;
      mateZ[i] = fx * Ai02 + fy * Ai12 + fz * Ai22;
    }
  }

  /**
   * Apply a fractional symmetry operator to one set of cartesian coordinates.
   *
   * @param xyz Input cartesian coordinates.
   * @param mate Symmetry mate cartesian coordinates.
   * @param symOp The fractional symmetry operator.
   */
  public void applySymOp(double[] xyz, double[] mate, SymOp symOp) {
    assert (xyz.length % 3 == 0);
    assert (xyz.length == mate.length);

    var rot = symOp.rot;
    var r00 = rot[0][0];
    var r01 = rot[0][1];
    var r02 = rot[0][2];
    var r10 = rot[1][0];
    var r11 = rot[1][1];
    var r12 = rot[1][2];
    var r20 = rot[2][0];
    var r21 = rot[2][1];
    var r22 = rot[2][2];

    var trans = symOp.tr;
    var t0 = trans[0];
    var t1 = trans[1];
    var t2 = trans[2];

    int len = xyz.length / 3;
    for (int i = 0; i < len; i++) {
      int index = i * 3;
      var xc = xyz[index + XX];
      var yc = xyz[index + YY];
      var zc = xyz[index + ZZ];
      // Convert to fractional coordinates.
      var xi = xc * A00 + yc * A10 + zc * A20;
      var yi = xc * A01 + yc * A11 + zc * A21;
      var zi = xc * A02 + yc * A12 + zc * A22;
      // Apply Symmetry Operator.
      var fx = r00 * xi + r01 * yi + r02 * zi + t0;
      var fy = r10 * xi + r11 * yi + r12 * zi + t1;
      var fz = r20 * xi + r21 * yi + r22 * zi + t2;
      // Convert back to Cartesian coordinates.
      mate[index + XX] = fx * Ai00 + fy * Ai10 + fz * Ai20;
      mate[index + YY] = fx * Ai01 + fy * Ai11 + fz * Ai21;
      mate[index + ZZ] = fx * Ai02 + fy * Ai12 + fz * Ai22;
    }
  }

  /**
   * Apply a fractional symmetry operator to one set of cartesian coordinates.
   *
   * @param xyz Input cartesian coordinates.
   * @param mate Symmetry mate cartesian coordinates.
   * @param symOp The fractional symmetry operator.
   */
  public void applySymRot(double[] xyz, double[] mate, SymOp symOp) {
    double[][] rot = symOp.rot;
    // Convert to fractional coordinates.
    double xc = xyz[0];
    double yc = xyz[1];
    double zc = xyz[2];
    // Convert to fractional coordinates.
    double xi = xc * A00 + yc * A10 + zc * A20;
    double yi = xc * A01 + yc * A11 + zc * A21;
    double zi = xc * A02 + yc * A12 + zc * A22;

    // Apply Symmetry Operator.
    double fx = rot[0][0] * xi + rot[0][1] * yi + rot[0][2] * zi;
    double fy = rot[1][0] * xi + rot[1][1] * yi + rot[1][2] * zi;
    double fz = rot[2][0] * xi + rot[2][1] * yi + rot[2][2] * zi;

    // Convert back to Cartesian coordinates.
    mate[0] = fx * Ai00 + fy * Ai10 + fz * Ai20;
    mate[1] = fx * Ai01 + fy * Ai11 + fz * Ai21;
    mate[2] = fx * Ai02 + fy * Ai12 + fz * Ai22;
  }

  /**
   * Apply the transpose of a symmetry rotation to an array of Cartesian coordinates. If the arrays
   * x, y or z are null or not of length n, the method returns immediately. If mateX, mateY or mateZ
   * are null or not of length n, new arrays are allocated.
   *
   * @param n Number of atoms.
   * @param x Input x-coordinates.
   * @param y Input y-coordinates.
   * @param z Input z-coordinates.
   * @param mateX Output x-coordinates.
   * @param mateY Output y-coordinates.
   * @param mateZ Output z-coordinates.
   * @param symOp The symmetry operator.
   * @param rotmat an array of double.
   */
  public void applyTransSymRot(
      int n, double[] x, double[] y, double[] z,
      double[] mateX, double[] mateY, double[] mateZ, SymOp symOp, double[][] rotmat) {

    if (x == null || y == null || z == null) {
      return;
    }
    if (x.length < n || y.length < n || z.length < n) {
      return;
    }
    if (mateX == null || mateX.length < n) {
      mateX = new double[n];
    }
    if (mateY == null || mateY.length < n) {
      mateY = new double[n];
    }
    if (mateZ == null || mateZ.length < n) {
      mateZ = new double[n];
    }

    // The transformation operator R = ToCart * Rot * ToFrac
    getTransformationOperator(symOp, rotmat);

    for (int i = 0; i < n; i++) {
      // Apply R^T (its transpose).
      double xc = x[i];
      double yc = y[i];
      double zc = z[i];
      mateX[i] = xc * rotmat[0][0] + yc * rotmat[1][0] + zc * rotmat[2][0];
      mateY[i] = xc * rotmat[0][1] + yc * rotmat[1][1] + zc * rotmat[2][1];
      mateZ[i] = xc * rotmat[0][2] + yc * rotmat[1][2] + zc * rotmat[2][2];
    }
  }

  /**
   * averageTensor
   *
   * @param m an array of double.
   * @param r an array of double.
   */
  public void averageTensor(double[][] m, double[][] r) {
    int n = spaceGroup.symOps.size();
    for (int i = 0; i < n; i++) {
      SymOp symop = spaceGroup.symOps.get(i);
      double[][] rot = symop.rot;
      double[][] rt = transpose3(rot);
      double[][] rmrt = mat3Mat3(mat3Mat3(rot, m), rt);
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          r[j][k] += rmrt[j][k];
        }
      }
    }
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        r[j][k] /= 6.0;
      }
    }
  }

  /**
   * averageTensor
   *
   * @param v an array of double.
   * @param r an array of double.
   */
  public void averageTensor(double[] v, double[][] r) {
    int n = spaceGroup.symOps.size();
    for (int i = 0; i < n; i++) {
      SymOp symop = spaceGroup.symOps.get(i);
      double[][] rot = symop.rot;
      double[][] rt = transpose3(rot);
      double[][] rmrt = mat3Mat3(mat3SymVec6(rot, v), rt);
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          r[j][k] += rmrt[j][k];
        }
      }
    }
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        r[j][k] /= 6.0;
      }
    }
  }

  /**
   * This method should be called to update the unit cell parameters of a crystal. The proposed
   * parameters will only be accepted if symmetry restrictions are satisfied. If so, all Crystal
   * variables that depend on the unit cell parameters will be updated.
   *
   * @param a length of the a-axis.
   * @param b length of the b-axis.
   * @param c length of the c-axis.
   * @param alpha Angle between b-axis and c-axis.
   * @param beta Angle between a-axis and c-axis.
   * @param gamma Angle between a-axis and b-axis.
   * @return The method return true if the parameters are accepted, false otherwise.
   */
  public boolean changeUnitCellParameters(
      double a, double b, double c, double alpha, double beta, double gamma) {
    if (checkRestrictions) {
      if (!latticeSystem.validParameters(a, b, c, alpha, beta, gamma)) {
        if (logger.isLoggable(Level.FINE)) {
          StringBuilder sb = new StringBuilder(
              " The proposed lattice parameters for " + spaceGroup.pdbName
                  + " do not satisfy the " + latticeSystem +
                  " lattice system restrictions and were ignored.\n");
          sb.append(format("  A-axis:                              %18.15e\n", a));
          sb.append(format("  B-axis:                              %18.15e\n", b));
          sb.append(format("  C-axis:                              %18.15e\n", c));
          sb.append(format("  Alpha:                               %18.15e\n", alpha));
          sb.append(format("  Beta:                                %18.15e\n", beta));
          sb.append(format("  Gamma:                               %18.15e\n", gamma));
          logger.fine(sb.toString());
        }
        return false;
      }
    }

    this.a = a;
    this.b = b;
    this.c = c;
    this.alpha = alpha;
    this.beta = beta;
    this.gamma = gamma;

    updateCrystal();

    return true;
  }

  /**
   * This method should be called to update the unit cell parameters of a crystal. The proposed
   * parameters will only be accepted if symmetry restrictions are satisfied. If so, all Crystal
   * variables that depend on the unit cell parameters will be updated.
   * <p>
   * If the new parameters are accepted, the target asymmetric unit volume is achieved by uniformly
   * scaling all lattice lengths.
   *
   * @param a length of the a-axis.
   * @param b length of the b-axis.
   * @param c length of the c-axis.
   * @param alpha Angle between b-axis and c-axis.
   * @param beta Angle between a-axis and c-axis.
   * @param gamma Angle between a-axis and b-axis.
   * @param targetAUVolume Target asymmetric unit volume.
   * @return The method return true if the parameters are accepted, false otherwise.
   */
  public boolean changeUnitCellParametersAndVolume(
      double a, double b, double c, double alpha, double beta, double gamma, double targetAUVolume) {
    if (changeUnitCellParameters(a, b, c, alpha, beta, gamma)) {
      double currentAUVolume = volume / getNumSymOps();
      double scale = cbrt(targetAUVolume / currentAUVolume);
      return changeUnitCellParameters(scale * a, scale * b, scale * c, alpha, beta, gamma);
    }
    return false;
  }

  /**
   * Two crystals are equal only if all unit cell parameters are exactly the same.
   *
   * @param o the Crystal to compare to.
   * @return true if all unit cell parameters are exactly the same.
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    Crystal crystal = (Crystal) o;
    return (a == crystal.a
        && b == crystal.b
        && c == crystal.c
        && alpha == crystal.alpha
        && beta == crystal.beta
        && gamma == crystal.gamma
        && spaceGroup.number == crystal.spaceGroup.number);
  }

  public boolean getCheckRestrictions() {
    return checkRestrictions;
  }

  public void setCheckRestrictions(boolean checkRestrictions) {
    this.checkRestrictions = checkRestrictions;
  }

  /**
   * Compute the density of the system.
   *
   * @param mass The total mass of the asymmetric unit.
   * @return The density (g/cc)
   */
  public double getDensity(double mass) {
    int nSymm = spaceGroup.symOps.size();
    return (mass * nSymm / AVOGADRO) * (1.0e24 / volume);
  }

  /**
   * Return the number of symmetry operators for this crystal.
   *
   * @return The number of symmetry operators.
   */
  public int getNumSymOps() {
    return spaceGroup.getNumberOfSymOps();
  }

  /**
   * Create a random Cartesian translation vector.
   *
   * <p>First, a random fractional translation is created. Second, the random fractional operator is
   * converted to Cartesian coordinates.
   *
   * @return A random Cartesian translation vector.
   */
  public double[] getRandomCartTranslation() {
    double[] coords = {random(), random(), random()};
    toCartesianCoordinates(coords, coords);
    return coords;
  }

  /**
   * Getter for the field <code>specialPositionCutoff</code>.
   *
   * @return a double.
   */
  public double getSpecialPositionCutoff() {
    return specialPositionCutoff;
  }

  /**
   * Getter for the field <code>specialPositionCutoff2</code>.
   *
   * @return a double.
   */
  public double getSpecialPositionCutoff2() {
    return specialPositionCutoff2;
  }

  /**
   * Setter for the field <code>specialPositionCutoff</code>.
   *
   * @param cutoff a double.
   */
  public void setSpecialPositionCutoff(double cutoff) {
    specialPositionCutoff = cutoff;
    specialPositionCutoff2 = cutoff * cutoff;
  }

  /**
   * Compute the total transformation operator R = ToCart * Rot * ToFrac.
   *
   * @param symOp Symmetry operator to apply.
   * @param rotmat Resulting transformation operator R.
   */
  public void getTransformationOperator(SymOp symOp, double[][] rotmat) {
    double[][] rot = symOp.rot;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rotmat[i][j] = 0.0;
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            rotmat[i][j] += Ai[k][i] * rot[k][l] * A[j][l];
          }
        }
      }
    }
  }

  /**
   * The ReplicatesCrystal over-rides this method to return the unit cell rather than the
   * ReplicateCell.
   *
   * @return The unit cell Crystal instance.
   */
  public Crystal getUnitCell() {
    return this;
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(a, b, c, alpha, beta, gamma, spaceGroup.number);
  }

  /**
   * Apply the minimum image convention.
   *
   * @param xyz input distances that are over-written.
   * @return the output distance squared.
   */
  public double image(final double[] xyz) {
    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];
    if (aperiodic) {
      return x * x + y * y + z * z;
    }
    double xf = x * A00 + y * A10 + z * A20;
    double yf = x * A01 + y * A11 + z * A21;
    double zf = x * A02 + y * A12 + z * A22;
    xf = floor(abs(xf) + 0.5) * signum(-xf) + xf;
    yf = floor(abs(yf) + 0.5) * signum(-yf) + yf;
    zf = floor(abs(zf) + 0.5) * signum(-zf) + zf;
    x = xf * Ai00 + yf * Ai10 + zf * Ai20;
    y = xf * Ai01 + yf * Ai11 + zf * Ai21;
    z = xf * Ai02 + yf * Ai12 + zf * Ai22;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    return x * x + y * y + z * z;
  }

  /**
   * Apply the minimum image convention.
   *
   * @param dx x-distance
   * @param dy y-distance
   * @param dz z-distance
   * @return the output distance squared.
   */
  public double image(double dx, double dy, double dz) {
    if (aperiodic) {
      return dx * dx + dy * dy + dz * dz;
    }
    double xf = dx * A00 + dy * A10 + dz * A20;
    double yf = dx * A01 + dy * A11 + dz * A21;
    double zf = dx * A02 + dy * A12 + dz * A22;
    xf = floor(abs(xf) + 0.5) * signum(-xf) + xf;
    yf = floor(abs(yf) + 0.5) * signum(-yf) + yf;
    zf = floor(abs(zf) + 0.5) * signum(-zf) + zf;
    dx = xf * Ai00 + yf * Ai10 + zf * Ai20;
    dy = xf * Ai01 + yf * Ai11 + zf * Ai21;
    dz = xf * Ai02 + yf * Ai12 + zf * Ai22;
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * Minimum distance between two coordinates over all symmetry operators.
   *
   * @param xyzA Coordinate A
   * @param xyzB Coordinate B
   * @return Minimum distance in crystal
   */
  public double minDistOverSymOps(double[] xyzA, double[] xyzB) {
    double dist = 0;
    for (int i = 0; i < 3; i++) {
      double dx = xyzA[i] - xyzB[i];
      dist += (dx * dx);
    }
    var symB = new double[3];
    for (SymOp symOp : spaceGroup.symOps) {
      applySymOp(xyzB, symB, symOp);
      for (int i = 0; i < 3; i++) {
        symB[i] -= xyzA[i];
      }
      double d = image(symB);
      dist = min(d, dist);
    }
    return sqrt(dist);
  }

  /**
   * Strain the unit cell vectors.
   *
   * @param dStrain a 3x3 matrix of unitless Strain percentages.
   * @return True if the perturbation of cell vectors succeeds.
   */
  public boolean perturbCellVectors(double[][] dStrain) {

    double[][] AA = new double[3][3];
    double[][] newAi = new double[3][3];

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        AA[i][j] = getUnitCell().Ai[i][j];
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          double Kronecker = 0.0;
          if (j == k) {
            Kronecker = 1.0;
          }
          // Aij' = Sum over K ( Kron_jk + dStrain_jk ) * Aij
          newAi[i][j] += (Kronecker + dStrain[j][k]) * AA[i][k];
        }
      }
    }

    return setCellVectors(newAi);
  }

  public boolean randomParameters(double dens, double mass) {
    double[] params = latticeSystem.resetUnitCellParams();
    boolean succeed =
        changeUnitCellParameters(params[0], params[1], params[2], params[3], params[4], params[5]);
    if (succeed) {
      setDensity(dens, mass);
    }
    return succeed;
  }

  /**
   * Is this a finite system - ie. one unit cell in isolation?
   *
   * @param aperiodic a boolean.
   */
  public void setAperiodic(boolean aperiodic) {
    this.aperiodic = aperiodic;
  }

  public double[] getCellParametersFromVectors(double[][] cellVectors) {
    // Update a-, b-, and c-axis lengths.
    double aa = length(cellVectors[0]);
    double bb = length(cellVectors[1]);
    double cc = length(cellVectors[2]);

    // Update alpha, beta and gamma angles.
    double aalpha = toDegrees(acos(dot(cellVectors[1], cellVectors[2]) / (bb * cc)));
    double bbeta = toDegrees(acos(dot(cellVectors[0], cellVectors[2]) / (aa * cc)));
    double ggamma = toDegrees(acos(dot(cellVectors[0], cellVectors[1]) / (aa * bb)));

    // Load and return new cell parameters.
    double[] params = new double[6];
    params[0] = aa;
    params[1] = bb;
    params[2] = cc;
    params[3] = aalpha;
    params[4] = bbeta;
    params[5] = ggamma;
    return params;
  }

  /**
   * Set the unit cell vectors.
   *
   * @param cellVectors 3x3 matrix of cell vectors.
   * @return True if the perturbation of cell vectors succeeds.
   */
  public boolean setCellVectors(double[][] cellVectors) {
    double[] p = getCellParametersFromVectors(cellVectors);
    return changeUnitCellParameters(p[0], p[1], p[2], p[3], p[4], p[5]);
  }

  /**
   * Set the unit cell vectors. Scale lattice lengths if necessary to hit the target volume.
   *
   * @param cellVectors 3x3 matrix of cell vectors.
   * @param targetAUVolume the target volume for the new cell Vectors.
   * @return True if the perturbation of cell vectors succeeds.
   */
  public boolean setCellVectorsAndVolume(double[][] cellVectors, double targetAUVolume) {
    if (setCellVectors(cellVectors)) {
      double currentAUVolume = volume / getNumSymOps();
      double scale = cbrt(targetAUVolume / currentAUVolume);
      return changeUnitCellParameters(scale * a, scale * b, scale * c, alpha, beta, gamma);
    } else {
      return false;
    }
  }

  public void setDensity(double dens, double mass) {
    double currentDensity = getDensity(mass);
    double scale = cbrt(currentDensity / dens);
    changeUnitCellParameters(a * scale, b * scale, c * scale, alpha, beta, gamma);
    currentDensity = getDensity(mass);
    logger.info(
        format(
            " Updated density %6.3f (g/cc) with unit cell %s.", currentDensity, toShortString()));
  }

  /**
   * Return a CRYST1 record useful for writing a PDB file.
   *
   * @return The CRYST1 record.
   */
  public String toCRYST1() {
    return format(
        "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s\n",
        a, b, c, alpha, beta, gamma, padRight(spaceGroup.pdbName, 10));
  }

  /**
   * toCartesianCoordinates
   *
   * @param n a int.
   * @param xf an array of double for fractional x-coordinates.
   * @param yf an array of double for fractional y-coordinates.
   * @param zf an array of double for fractional z-coordinates.
   * @param x an array of double for cartesian x-coordinates.
   * @param y an array of double for cartesian y-coordinates.
   * @param z an array of double for cartesian z-coordinates.
   */
  public void toCartesianCoordinates(
      int n, double[] xf, double[] yf, double[] zf, double[] x, double[] y, double[] z) {
    for (int i = 0; i < n; i++) {
      double xi = xf[i];
      double yi = yf[i];
      double zi = zf[i];
      x[i] = xi * Ai00 + yi * Ai10 + zi * Ai20;
      y[i] = xi * Ai01 + yi * Ai11 + zi * Ai21;
      z[i] = xi * Ai02 + yi * Ai12 + zi * Ai22;
    }
  }

  /**
   * toCartesianCoordinates
   *
   * @param n a int.
   * @param frac an array of double for fractional coordinates.
   * @param cart an array of double for cartesian coordinates.
   */
  public void toCartesianCoordinates(int n, double[] frac, double[] cart) {
    int i3 = 0;
    for (int i = 0; i < n; i++) {
      // Convert to cartesian coordinates.
      int iX = i3 + XX;
      int iY = i3 + YY;
      int iZ = i3 + ZZ;
      i3 += 3;
      double xf = frac[iX];
      double yf = frac[iY];
      double zf = frac[iZ];
      cart[iX] = xf * Ai00 + yf * Ai10 + zf * Ai20;
      cart[iY] = xf * Ai01 + yf * Ai11 + zf * Ai21;
      cart[iZ] = xf * Ai02 + yf * Ai12 + zf * Ai22;
    }
  }

  /**
   * toCartesianCoordinates
   *
   * @param xf an array of double for fractional coordinate.
   * @param x an array of double for cartesian coordinate.
   */
  public void toCartesianCoordinates(double[] xf, double[] x) {
    double fx = xf[0];
    double fy = xf[1];
    double fz = xf[2];
    x[0] = fx * Ai00 + fy * Ai10 + fz * Ai20;
    x[1] = fx * Ai01 + fy * Ai11 + fz * Ai21;
    x[2] = fx * Ai02 + fy * Ai12 + fz * Ai22;
  }

  /**
   * toFractionalCoordinates
   *
   * @param n a int.
   * @param x an array of double for cartesian x-coordinates.
   * @param y an array of double for cartesian y-coordinates.
   * @param z an array of double for cartesian z-coordinates.
   * @param xf an array of double for fractional x-coordinates.
   * @param yf an array of double for fractional y-coordinates.
   * @param zf an array of double for fractional z-coordinates.
   */
  public void toFractionalCoordinates(
      int n, double[] x, double[] y, double[] z, double[] xf, double[] yf, double[] zf) {
    for (int i = 0; i < n; i++) {
      double xc = x[i];
      double yc = y[i];
      double zc = z[i];
      xf[i] = xc * A00 + yc * A10 + zc * A20;
      yf[i] = xc * A01 + yc * A11 + zc * A21;
      zf[i] = xc * A02 + yc * A12 + zc * A22;
    }
  }

  /**
   * toFractionalCoordinates
   *
   * @param n a int.
   * @param cart an array of double for cartesian coordinates.
   * @param frac an array of double for fractional coordinates.
   */
  public void toFractionalCoordinates(int n, double[] cart, double[] frac) {
    int i3 = 0;
    for (int i = 0; i < n; i++) {
      // Convert to fractional coordinates.
      int iX = i3 + XX;
      int iY = i3 + YY;
      int iZ = i3 + ZZ;
      i3 += 3;
      double xc = cart[iX];
      double yc = cart[iY];
      double zc = cart[iZ];
      frac[iX] = xc * A00 + yc * A10 + zc * A20;
      frac[iY] = xc * A01 + yc * A11 + zc * A21;
      frac[iZ] = xc * A02 + yc * A12 + zc * A22;
    }
  }

  /**
   * toFractionalCoordinates
   *
   * @param x an array of double for cartesian coordinate.
   * @param xf an array of double for fractional coordinate.
   */
  public void toFractionalCoordinates(double[] x, double[] xf) {
    double xc = x[0];
    double yc = x[1];
    double zc = x[2];
    xf[0] = xc * A00 + yc * A10 + zc * A20;
    xf[1] = xc * A01 + yc * A11 + zc * A21;
    xf[2] = xc * A02 + yc * A12 + zc * A22;
  }

  /**
   * toPrimaryCell
   *
   * @param in an array of double.
   * @param out an array of double.
   */
  public void toPrimaryCell(double[] in, double[] out) {
    toFractionalCoordinates(in, out);
    out[0] = mod(out[0], 1.0);
    out[1] = mod(out[1], 1.0);
    out[2] = mod(out[2], 1.0);
    toCartesianCoordinates(out, out);
  }

  /**
   * A String containing the unit cell parameters.
   *
   * @return A string with the unit cell parameters.
   */
  public String toShortString() {
    return format("%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f", a, b, c, alpha, beta, gamma);
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder("\n Unit Cell\n");
    sb.append(
        format("  A-axis:                              %8.3f (%8.3f, %8.3f, %8.3f)\n", a, Ai00, Ai01,
            Ai02));
    sb.append(
        format("  B-axis:                              %8.3f (%8.3f, %8.3f, %8.3f)\n", b, Ai10, Ai11,
            Ai12));
    sb.append(
        format("  C-axis:                              %8.3f (%8.3f, %8.3f, %8.3f)\n", c, Ai20, Ai21,
            Ai22));
    sb.append(format("  Alpha:                               %8.3f\n", alpha));
    sb.append(format("  Beta:                                %8.3f\n", beta));
    sb.append(format("  Gamma:                               %8.3f\n", gamma));
    sb.append("  Space group\n");
    sb.append(format("   Number:                                  %3d\n", spaceGroup.number));
    sb.append(format("   Symbol:                             %8s\n", spaceGroup.shortName));
    sb.append(
        format("   Number of Symmetry Operators:            %3d", spaceGroup.getNumberOfSymOps()));
    return sb.toString();
  }

  /** Update all Crystal variables that are a function of unit cell parameters. */
  public void updateCrystal() {

    double cos_alpha;
    double sin_beta;
    double cos_beta;
    double sin_gamma;
    double cos_gamma;
    double beta_term;
    double gamma_term;

    switch (crystalSystem) {
      case CUBIC:
      case ORTHORHOMBIC:
      case TETRAGONAL:
        cos_alpha = 0.0;
        cos_beta = 0.0;
        sin_gamma = 1.0;
        cos_gamma = 0.0;
        beta_term = 0.0;
        gamma_term = 1.0;
        volume = a * b * c;
        dVdA = b * c;
        dVdB = a * c;
        dVdC = a * b;
        dVdAlpha = 0.0;
        dVdBeta = 0.0;
        dVdGamma = 0.0;
        break;
      case MONOCLINIC:
        cos_alpha = 0.0;
        sin_beta = sin(toRadians(beta));
        cos_beta = cos(toRadians(beta));
        sin_gamma = 1.0;
        cos_gamma = 0.0;
        beta_term = 0.0;
        gamma_term = sin_beta;
        volume = sin_beta * a * b * c;
        dVdA = sin_beta * b * c;
        dVdB = sin_beta * a * c;
        dVdC = sin_beta * a * b;
        dVdAlpha = 0.0;
        dVdBeta = cos_beta * a * b * c;
        dVdGamma = 0.0;
        break;
      case HEXAGONAL:
      case TRICLINIC:
      case TRIGONAL:
      default:
        double sin_alpha = sin(toRadians(alpha));
        cos_alpha = cos(toRadians(alpha));
        sin_beta = sin(toRadians(beta));
        cos_beta = cos(toRadians(beta));
        sin_gamma = sin(toRadians(gamma));
        cos_gamma = cos(toRadians(gamma));
        beta_term = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        gamma_term = sqrt(sin_beta * sin_beta - beta_term * beta_term);
        volume = sin_gamma * gamma_term * a * b * c;

        dVdA = sin_gamma * gamma_term * b * c;
        dVdB = sin_gamma * gamma_term * a * c;
        dVdC = sin_gamma * gamma_term * a * b;

        double dbeta =
            2.0 * sin_beta * cos_beta
                - (2.0 * (cos_alpha - cos_beta * cos_gamma) * sin_beta * cos_gamma)
                / (sin_gamma * sin_gamma);
        double dgamma1 = -2.0 * (cos_alpha - cos_beta * cos_gamma) * cos_beta / sin_gamma;
        double dgamma2 = cos_alpha - cos_beta * cos_gamma;
        dgamma2 *= dgamma2 * 2.0 * cos_gamma / (sin_gamma * sin_gamma * sin_gamma);

        dVdAlpha =
            (cos_alpha - cos_beta * cos_gamma) * sin_alpha / (sin_gamma * gamma_term) * a * b * c;
        dVdBeta = 0.5 * sin_gamma * dbeta / gamma_term * a * b * c;
        dVdGamma =
            (cos_gamma * gamma_term + 0.5 * sin_gamma * (dgamma1 + dgamma2) / gamma_term)
                * a
                * b
                * c;

        break;
    }

    G[0][0] = a * a;
    G[0][1] = a * b * cos_gamma;
    G[0][2] = a * c * cos_beta;
    G[1][0] = G[0][1];
    G[1][1] = b * b;
    G[1][2] = b * c * cos_alpha;
    G[2][0] = G[0][2];
    G[2][1] = G[1][2];
    G[2][2] = c * c;

    // invert G to yield Gstar
    RealMatrix m = new Array2DRowRealMatrix(G, true);
    m = new LUDecomposition(m).getSolver().getInverse();
    Gstar = m.getData();

    // a is the first row of A^(-1).
    Ai[0][0] = a;
    Ai[0][1] = 0.0;
    Ai[0][2] = 0.0;
    // b is the second row of A^(-1).
    Ai[1][0] = b * cos_gamma;
    Ai[1][1] = b * sin_gamma;
    Ai[1][2] = 0.0;
    // c is the third row of A^(-1).
    Ai[2][0] = c * cos_beta;
    Ai[2][1] = c * beta_term;
    Ai[2][2] = c * gamma_term;

    Ai00 = Ai[0][0];
    Ai01 = Ai[0][1];
    Ai02 = Ai[0][2];
    Ai10 = Ai[1][0];
    Ai11 = Ai[1][1];
    Ai12 = Ai[1][2];
    Ai20 = Ai[2][0];
    Ai21 = Ai[2][1];
    Ai22 = Ai[2][2];

    // Invert A^-1 to get A
    m = new Array2DRowRealMatrix(Ai, true);
    m = new LUDecomposition(m).getSolver().getInverse();
    A = m.getData();

    // The columns of A are the reciprocal basis vectors
    A00 = A[0][0];
    A10 = A[1][0];
    A20 = A[2][0];
    A01 = A[0][1];
    A11 = A[1][1];
    A21 = A[2][1];
    A02 = A[0][2];
    A12 = A[1][2];
    A22 = A[2][2];

    // Reciprocal basis vector lengths
    double aStar = 1.0 / sqrt(A00 * A00 + A10 * A10 + A20 * A20);
    double bStar = 1.0 / sqrt(A01 * A01 + A11 * A11 + A21 * A21);
    double cStar = 1.0 / sqrt(A02 * A02 + A12 * A12 + A22 * A22);
    if (logger.isLoggable(Level.FINEST)) {
      logger.finest(
          format(" Reciprocal Lattice Lengths: (%8.3f, %8.3f, %8.3f)", aStar, bStar, cStar));
    }

    // Interfacial diameters from the dot product of the real and reciprocal vectors
    interfacialRadiusA = (Ai00 * A00 + Ai01 * A10 + Ai02 * A20) * aStar;
    interfacialRadiusB = (Ai10 * A01 + Ai11 * A11 + Ai12 * A21) * bStar;
    interfacialRadiusC = (Ai20 * A02 + Ai21 * A12 + Ai22 * A22) * cStar;

    // Divide by 2 to get radii.
    interfacialRadiusA /= 2.0;
    interfacialRadiusB /= 2.0;
    interfacialRadiusC /= 2.0;

    if (logger.isLoggable(Level.FINEST)) {
      logger.finest(
          format(
              " Interfacial radii: (%8.3f, %8.3f, %8.3f)",
              interfacialRadiusA, interfacialRadiusB, interfacialRadiusC));
    }

    List<SymOp> symOps = spaceGroup.symOps;
    int nSymm = symOps.size();
    if (symOpsCartesian == null) {
      symOpsCartesian = new ArrayList<>(nSymm);
    } else {
      symOpsCartesian.clear();
    }

    RealMatrix toFrac = new Array2DRowRealMatrix(A);
    RealMatrix toCart = new Array2DRowRealMatrix(Ai);
    for (SymOp symOp : symOps) {
      m = new Array2DRowRealMatrix(symOp.rot);
      // rot_c = A^(-1).rot_f.A
      RealMatrix rotMat = m.preMultiply(toCart.transpose()).multiply(toFrac.transpose());
      // tr_c = tr_f.A^(-1)
      double[] tr = toCart.preMultiply(symOp.tr);
      symOpsCartesian.add(new SymOp(rotMat.getData(), tr));
    }
  }
}
