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
package ffx.potential.bonded;

import static ffx.numerics.math.DoubleMath.X;
import static ffx.numerics.math.DoubleMath.angle;
import static ffx.numerics.math.DoubleMath.length;
import static ffx.numerics.math.DoubleMath.normalize;
import static ffx.numerics.math.DoubleMath.scale;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.potential.parameters.BondType.cubic;
import static ffx.potential.parameters.BondType.quartic;
import static ffx.potential.parameters.BondType.units;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.math.Double3;
import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.RendererCache.ViewModel;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Geometry;
import org.jogamp.java3d.LineArray;
import org.jogamp.java3d.Shape3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.vecmath.AxisAngle4d;
import org.jogamp.vecmath.Vector3d;

/**
 * The Bond class represents a covalent bond formed between two atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("CloneableImplementsClone")
public class Bond extends BondedTerm {

  /**
   * Length in Angstroms that is added to Atomic Radii when determining if two Atoms are within
   * bonding distance
   */
  static final float BUFF = 0.7f;

  private static final Logger logger = Logger.getLogger(Bond.class.getName());
  // Some static variables used for computing cylinder orientations
  private static final float[] a0col = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
  };
  private static final float[] f4a = {0.0f, 0.0f, 0.0f, 0.9f};
  private static final float[] f4b = {0.0f, 0.0f, 0.0f, 0.9f};
  private static final float[] f16 = {
    0.0f, 0.0f, 0.0f, 0.9f, 0.0f, 0.0f, 0.0f, 0.9f, 0.0f, 0.0f, 0.0f, 0.9f, 0.0f, 0.0f, 0.0f, 0.9f
  };
  private static final double[] a13d = new double[3];
  private static final double[] a23d = new double[3];
  private static final double[] mid = new double[3];
  private static final double[] diff3d = new double[3];
  private static final double[] sum3d = new double[3];
  private static final double[] coord = new double[12];
  private static final double[] y = {0.0d, 1.0d, 0.0d};
  private static final AxisAngle4d axisAngle = new AxisAngle4d();
  private static final double[] bcross = new double[4];
  private static final double[] cstart = new double[3];
  private static final Vector3d pos3d = new Vector3d();
  /** List of Bonds that this Bond forms angles with */
  private final ArrayList<Bond> formsAngleWith = new ArrayList<>();
  /** The force field BondType for this bond. */
  public BondType bondType = null;
  /** Rigid Scale factor. */
  private double rigidScale = 1.0;
  // Java3D methods and variables for visualization of this Bond.
  private RendererCache.ViewModel viewModel = RendererCache.ViewModel.INVISIBLE;
  private BranchGroup branchGroup;
  private TransformGroup cy1tg, cy2tg;
  private Transform3D cy1t3d, cy2t3d;
  private Shape3D cy1, cy2;
  private Vector3d scale;
  private int detail = 3;
  private LineArray la;
  private int lineIndex;
  private boolean wireVisible = true;

  /**
   * Simple Bond constructor that is intended to be used with the equals method.
   *
   * @param n Bond id
   */
  public Bond(String n) {
    super(n);
  }

  /**
   * Bond constructor.
   *
   * @param a1 Atom number 1.
   * @param a2 Atom number 2.
   */
  public Bond(Atom a1, Atom a2) {
    atoms = new Atom[2];
    int i1 = a1.getIndex();
    int i2 = a2.getIndex();
    if (i1 < i2) {
      atoms[0] = a1;
      atoms[1] = a2;
    } else {
      atoms[0] = a2;
      atoms[1] = a1;
    }
    setID_Key(false);
    viewModel = RendererCache.ViewModel.WIREFRAME;
    a1.setBond(this);
    a2.setBond(this);
  }

  /**
   * Log that no BondType exists.
   *
   * @param a1 Atom 1.
   * @param a2 Atom 2.
   * @param key The class key.
   * @param forceField The force field in use.
   */
  public static void logNoBondType(Atom a1, Atom a2, String key, ForceField forceField) {
    AtomType atomType1 = a1.getAtomType();
    AtomType atomType2 = a2.getAtomType();
    StringBuilder sb = new StringBuilder(
        format(" No BondType for key: %s\n %s -> %s\n %s -> %s", key,
        a1, atomType1, a2, atomType2));
    int c1 = atomType1.atomClass;
    int c2 = atomType2.atomClass;
    List<AtomType> types1 = forceField.getSimilarAtomTypes(atomType1);
    List<AtomType> types2 = forceField.getSimilarAtomTypes(atomType2);
    List<BondType> bondTypes = new ArrayList<>();
    boolean match = false;
    for (AtomType type1 : types1) {
      for (AtomType type2 : types2) {
        // Similar bond type must match at least one class.
        if ((type1.atomClass != c1) && (type1.atomClass != c2) &&
            (type2.atomClass != c1) && (type2.atomClass != c2)) {
          continue;
        }
        int[] c = new int[2];
        c[0] = type1.atomClass;
        c[1] = type2.atomClass;
        String closeKey = BondType.sortKey(c);
        BondType bondType = forceField.getBondType(closeKey);
        if (bondType != null && !bondTypes.contains(bondType)) {
          if (!match) {
            match = true;
            sb.append("\n Similar Bond Types:");
          }
          bondTypes.add(bondType);
          sb.append(format("\n  %s", bondType));
        }
      }
    }

    logger.severe(sb.toString());
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(BondedTerm b) {
    if (b == null) {
      throw new NullPointerException();
    }
    if (b == this) {
      return 0;
    }
    if (!b.getClass().isInstance(this)) {
      return super.compareTo(b);
    }
    int this0 = atoms[0].getIndex();
    int a0 = b.atoms[0].getIndex();
    if (this0 < a0) {
      return -1;
    }
    if (this0 > a0) {
      return 1;
    }
    int this1 = atoms[1].getIndex();
    int a1 = b.atoms[1].getIndex();

    return Integer.compare(this1, a1);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Evaluate this Bond energy.
   */
  @Override
  public double energy(
      boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
    var atomA = atoms[0];
    var atomB = atoms[1];
    var va = atomA.getXYZ();
    var vb = atomB.getXYZ();
    var vab = va.sub(vb);
    value = vab.length();
    var prefactor = units * rigidScale * bondType.forceConstant;
    // dv is deviation from ideal
    var dv = value - bondType.distance;
    if (bondType.bondFunction.hasFlatBottom()) {
      if (dv > 0) {
        dv = max(0, dv - bondType.flatBottomRadius);
      } else if (dv < 0) {
        dv = min(0, dv + bondType.flatBottomRadius);
      } // Else, no adjustments needed.
    }
    var dv2 = dv * dv;
    switch (bondType.bondFunction) {
      case QUARTIC:
      case FLAT_BOTTOM_QUARTIC:
        {
          energy = prefactor * dv2 * (1.0 + cubic * dv + quartic * dv2);
          if (esvTerm) {
            setEsvDeriv(energy * dedesvChain);
            energy = energy * esvLambda;
          }
          if (gradient) {
            // Compute the magnitude of the gradient.
            var dedr = 2.0 * prefactor * esvLambda * dv * (1.0 + 1.5 * cubic * dv + 2.0 * quartic * dv2);
            computeGradient(threadID, grad, atomA, atomB, vab, dedr);
          }
          break;
        }
      case HARMONIC:
      case FLAT_BOTTOM_HARMONIC:
      default:
        {
          energy = prefactor * dv2;
          if (esvTerm) {
            setEsvDeriv(energy * dedesvChain);
            energy = energy * esvLambda;
          }
          if (gradient) {
            // Compute the magnitude of the gradient.
            var dedr = 2.0 * prefactor * esvLambda * dv;
            computeGradient(threadID, grad, atomA, atomB, vab, dedr);
          }
          break;
        }
    }
    value = dv;
    return energy;
  }

  /**
   * Find the other Atom in <b>this</b> Bond. These two atoms are said to be 1-2.
   *
   * @param a The known Atom.
   * @return The other Atom that makes up <b>this</b> Bond, or Null if Atom a is not part of
   *     <b>this</b> Bond.
   */
  public Atom get1_2(Atom a) {
    if (a == atoms[0]) {
      return atoms[1];
    }
    if (a == atoms[1]) {
      return atoms[0];
    }
    return null; // Atom not found in bond
  }

  /**
   * Gets the current distance between the two Atoms in this Bond.
   *
   * @return Current distance.
   */
  public double getCurrentDistance() {
    double[] x1 = new double[3];
    x1 = atoms[0].getXYZ(x1);
    double[] x2 = new double[3];
    x2 = atoms[1].getXYZ(x2);
    return DoubleMath.dist(x1, x2);
  }

  /** Log details for this Bond energy term. */
  public void log() {
    logger.info(
        format(
            " %-8s %6d-%s %6d-%s %6.4f  %6.4f  %10.4f",
            "Bond",
            atoms[0].getIndex(),
            atoms[0].getAtomType().name,
            atoms[1].getIndex(),
            atoms[1].getAtomType().name,
            bondType.distance,
            value,
            energy));
  }

  /** {@inheritDoc} */
  @Override
  public void removeFromParent() {
    super.removeFromParent();
    cy1 = null;
    cy2 = null;
    cy1tg = null;
    cy2tg = null;
    if (cy1t3d != null) {
      RendererCache.poolTransform3D(cy1t3d);
      RendererCache.poolTransform3D(cy2t3d);
      cy1t3d = null;
      cy2t3d = null;
    }
    if (branchGroup != null) {
      branchGroup.detach();
      branchGroup.setUserData(null);
      RendererCache.poolDoubleCylinder(branchGroup);
      branchGroup = null;
    }
  }

  /**
   * Set a reference to the force field parameters.
   *
   * @param bondType a {@link ffx.potential.parameters.BondType} object.
   */
  public void setBondType(BondType bondType) {
    this.bondType = bondType;
  }

  /**
   * Return the BondType for this Bond.
   * @return Returns the BondType.
   */
  public BondType getBondType() {
    return bondType;
  }

  /**
   * Set the color of this Bond's Java3D shapes based on the passed Atom.
   *
   * @param a Atom
   */
  public void setColor(Atom a) {
    if (viewModel != ViewModel.INVISIBLE
        && viewModel != ViewModel.WIREFRAME
        && branchGroup != null) {
      if (a == atoms[0]) {
        cy1.setAppearance(a.getAtomAppearance());
      } else if (a == atoms[1]) {
        cy2.setAppearance(a.getAtomAppearance());
      }
    }
    setWireVisible(wireVisible);
  }

  /**
   * Setter for the field <code>rigidScale</code>.
   *
   * @param rigidScale a double.
   */
  public void setRigidScale(double rigidScale) {
    this.rigidScale = rigidScale;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Polymorphic setView method.
   */
  @Override
  public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
    switch (newViewModel) {
      case WIREFRAME:
        viewModel = ViewModel.WIREFRAME;
        setWireVisible(true);
        setCylinderVisible(false, newShapes);
        break;
      case SPACEFILL:
      case INVISIBLE:
      case RMIN:
        viewModel = ViewModel.INVISIBLE;
        setWireVisible(false);
        setCylinderVisible(false, newShapes);
        break;
      case RESTRICT:
        if (!atoms[0].isSelected() || !atoms[1].isSelected()) {
          viewModel = ViewModel.INVISIBLE;
          setWireVisible(false);
          setCylinderVisible(false, newShapes);
        }
        break;
      case BALLANDSTICK:
      case TUBE:
        viewModel = newViewModel;
        // Get the radius to use
        double rad;
        double len = getValue() / 2.0d;
        if (viewModel == RendererCache.ViewModel.BALLANDSTICK) {
          rad = 0.1d * RendererCache.radius;
        } else {
          rad = 0.2d * RendererCache.radius;
        }
        if (scale == null) {
          scale = new Vector3d();
        }
        scale.set(rad, len, rad);
        setWireVisible(false);
        setCylinderVisible(true, newShapes);
        break;
      case DETAIL:
        int res = RendererCache.detail;
        if (res != detail) {
          detail = res;
          if (branchGroup != null) {
            Geometry geom1 = RendererCache.getCylinderGeom(0, detail);
            Geometry geom2 = RendererCache.getCylinderGeom(1, detail);
            Geometry geom3 = RendererCache.getCylinderGeom(2, detail);
            cy1.removeAllGeometries();
            cy2.removeAllGeometries();
            cy1.addGeometry(geom1);
            cy1.addGeometry(geom2);
            cy1.addGeometry(geom3);
            cy2.addGeometry(geom1);
            cy2.addGeometry(geom2);
            cy2.addGeometry(geom3);
          }
        }
        if (scale == null) {
          scale = new Vector3d();
        }
        double newRadius;
        if (viewModel == RendererCache.ViewModel.BALLANDSTICK) {
          newRadius = 0.1d * RendererCache.radius;
        } else if (viewModel == RendererCache.ViewModel.TUBE) {
          newRadius = 0.2d * RendererCache.radius;
        } else {
          break;
        }
        if (newRadius != scale.x) {
          scale.x = newRadius;
          scale.y = newRadius;
          if (branchGroup != null) {
            setView(viewModel, newShapes);
          }
        }
        break;
      case SHOWHYDROGENS:
        if (atoms[0].getAtomicNumber() == 1 || atoms[1].getAtomicNumber() == 1) {
          setView(viewModel, newShapes);
        }
        break;
      case HIDEHYDROGENS:
        if (atoms[0].getAtomicNumber() == 1 || atoms[1].getAtomicNumber() == 1) {
          viewModel = ViewModel.INVISIBLE;
          setWireVisible(false);
          setCylinderVisible(false, newShapes);
        }
        break;
      case FILL:
      case POINTS:
      case LINES:
        if (branchGroup != null && viewModel != ViewModel.INVISIBLE) {
          cy1.setAppearance(atoms[0].getAtomAppearance());
          cy2.setAppearance(atoms[1].getAtomAppearance());
        }
        break;
    }
  }

  /**
   * setWire
   *
   * @param l a {@link org.jogamp.java3d.LineArray} object.
   * @param i a int.
   */
  public void setWire(LineArray l, int i) {
    la = l;
    lineIndex = i;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Update recomputes the bonds length, Wireframe vertices, and Cylinder Transforms
   */
  @Override
  public void update() {
    // Update the Bond Length
    atoms[0].getXYZ(a13d);
    atoms[1].getXYZ(a23d);
    sub(a13d, a23d, diff3d);
    double d = length(diff3d);
    setValue(d);
    DoubleMath.add(a13d, a23d, sum3d);
    scale(sum3d, 0.5d, mid);
    // Update the Wireframe Model.
    if (la != null) {
      for (int i = 0; i < 3; i++) {
        coord[i] = a13d[i];
        coord[3 + i] = mid[i];
        coord[6 + i] = mid[i];
        coord[9 + i] = a23d[i];
      }
      la.setCoordinates(lineIndex, coord);
    }
    // Update the Bond cylinder transforms.
    if (branchGroup != null) {
      normalize(diff3d, diff3d);
      scale.y = d / 2.0d;
      setBondTransform3d(cy1t3d, mid, diff3d, d, true);
      scale(diff3d, -1.0d, diff3d);
      setBondTransform3d(cy2t3d, mid, diff3d, d, false);
      cy1tg.setTransform(cy1t3d);
      cy2tg.setTransform(cy2t3d);
    }
  }

  private void computeGradient(
      int threadID, AtomicDoubleArray3D grad, Atom atomA, Atom atomB, Double3 vab, double dedr) {
    var de = 0.0;
    if (value > 0.0) {
      de = dedr / value;
    }
    var ga = vab.scale(de);
    var ia = atomA.getIndex() - 1;
    var ib = atomB.getIndex() - 1;
    grad.add(threadID, ia, ga);
    grad.sub(threadID, ib, ga);
  }

  /**
   * Check to see if <b>this</b> Bond and another combine to form an angle
   *
   * @param b a {@link ffx.potential.bonded.Bond} object.
   * @return True if Bond b helps form an angle with <b>this</b> Bond
   */
  boolean formsAngleWith(Bond b) {
    for (Bond bond : formsAngleWith) {
      if (b == bond) {
        return true;
      }
    }
    return false;
  }

  /**
   * Finds the common Atom between <b>this</b> Bond and Bond b.
   *
   * @param b Bond to compare with.
   * @return The Atom the Bonds have in common or Null if they are the same Bond or have no atom in
   *     common
   */
  Atom getCommonAtom(Bond b) {
    if (b == this || b == null) {
      return null;
    }
    if (b.atoms[0] == atoms[0]) {
      return atoms[0];
    }
    if (b.atoms[0] == atoms[1]) {
      return atoms[1];
    }
    if (b.atoms[1] == atoms[0]) {
      return atoms[0];
    }
    if (b.atoms[1] == atoms[1]) {
      return atoms[1];
    }
    return null; // Common atom not found
  }

  /**
   * Find the Atom that <b>this</b> Bond and Bond b do not have in common.
   *
   * @param b Bond to compare with
   * @return The Atom that Bond b and <b>this</b> Bond do not have in common, or Null if they have
   *     no Atom in common
   */
  Atom getOtherAtom(Bond b) {
    if (b == this || b == null) {
      return null;
    }
    if (b.atoms[0] == atoms[0]) {
      return atoms[1];
    }
    if (b.atoms[0] == atoms[1]) {
      return atoms[0];
    }
    if (b.atoms[1] == atoms[0]) {
      return atoms[1];
    }
    if (b.atoms[1] == atoms[1]) {
      return atoms[0];
    }
    return null;
  }

  /**
   * sameGroup
   *
   * @return a boolean.
   */
  boolean sameGroup() {
    return atoms[0].getParent() == atoms[1].getParent();
  }

  /**
   * Specifies <b>this</b> Bond helps form an angle with the given Bond
   *
   * @param b Bond that forms an angle with <b>this</b> Bond
   */
  void setAngleWith(Bond b) {
    formsAngleWith.add(b);
  }

  /**
   * Create the Bond Scenegraph Objects.
   *
   * @param newShapes List
   */
  private void initJ3D(List<BranchGroup> newShapes) {
    detail = RendererCache.detail;
    branchGroup = RendererCache.doubleCylinderFactory(atoms[0], atoms[1], detail);
    cy1tg = (TransformGroup) branchGroup.getChild(0);
    cy2tg = (TransformGroup) branchGroup.getChild(1);
    cy1 = (Shape3D) cy1tg.getChild(0);
    cy2 = (Shape3D) cy2tg.getChild(0);
    newShapes.add(branchGroup);
    cy1t3d = RendererCache.transform3DFactory();
    cy2t3d = RendererCache.transform3DFactory();
    update();
  }

  /**
   * setBondTransform3d
   *
   * @param t3d a {@link org.jogamp.java3d.Transform3D} object.
   * @param pos an array of double.
   * @param orient an array of double.
   * @param len a double.
   * @param newRot a boolean.
   */
  private void setBondTransform3d(
      Transform3D t3d, double[] pos, double[] orient, double len, boolean newRot) {
    // Bond Orientation
    if (newRot) {
      double angle = angle(orient, y);
      X(y, orient, bcross);
      bcross[3] = angle - Math.PI;
      axisAngle.set(bcross);
    }
    // Scale the orientation vector to be a fourth the bond length
    // and add it to the position vector of the of the first atom
    scale(orient, len / 4.0d, cstart);
    DoubleMath.add(cstart, pos, cstart);
    pos3d.set(cstart);
    t3d.setTranslation(pos3d);
    t3d.setRotation(axisAngle);
    t3d.setScale(scale);
  }

  /**
   * Manage cylinder visibility.
   *
   * @param visible boolean
   * @param newShapes List
   */
  private void setCylinderVisible(boolean visible, List<BranchGroup> newShapes) {
    if (!visible) {
      // Make this Bond invisible.
      if (branchGroup != null) {
        cy1.setPickable(false);
        cy1.setAppearance(RendererCache.nullAp);
        cy2.setPickable(false);
        cy2.setAppearance(RendererCache.nullAp);
        // branchGroup = null;
      }
    } else if (branchGroup == null) {
      // Get Java3D primitives from the RendererCache
      initJ3D(newShapes);
    } else {
      // Scale the cylinders to match the current ViewModel
      cy1t3d.setScale(scale);
      cy1tg.setTransform(cy1t3d);
      cy2t3d.setScale(scale);
      cy2tg.setTransform(cy2t3d);
      cy1.setAppearance(atoms[0].getAtomAppearance());
      cy2.setAppearance(atoms[1].getAtomAppearance());
    }
  }

  /**
   * Manage wireframe visibility.
   *
   * @param visible a boolean.
   */
  private void setWireVisible(boolean visible) {
    if (!visible) {
      wireVisible = false;
      la.setColors(lineIndex, a0col);
    } else {
      wireVisible = true;
      float[] cols = f16;
      float[] col1 = f4a;
      float[] col2 = f4b;
      atoms[0].getAtomColor().get(col1);
      atoms[1].getAtomColor().get(col2);
      for (int i = 0; i < 3; i++) {
        cols[i] = col1[i];
        cols[4 + i] = col1[i];
        cols[8 + i] = col2[i];
        cols[12 + i] = col2[i];
      }
      la.setColors(lineIndex, cols);
    }
  }
}
