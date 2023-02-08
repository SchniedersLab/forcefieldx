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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.Objects;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Canvas3D;
import org.jogamp.java3d.J3DGraphics2D;
import org.jogamp.java3d.Material;
import org.jogamp.java3d.Node;
import org.jogamp.vecmath.Color3f;

/**
 * The MSNode class forms the basic unit that all data classes extend.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("CloneableImplementsClone")
public class MSNode extends DefaultMutableTreeNode implements ROLS {

  /** The multiscale level of this node. */
  private final int MultiScaleLevel;
  /** True if this node is selected. */
  protected boolean selected = false;
  /** The name of this node. */
  private String name;
  /** Total mass of this node and its children. */
  private double totalMass;

  /** Default MSNode Constructor */
  public MSNode() {
    name = "";
    MultiScaleLevel = ROLS.MaxLengthScale;
  }

  /**
   * Constructs a MSNode with the name n.
   *
   * @param n a {@link java.lang.String} object.
   */
  public MSNode(String n) {
    name = n;
    MultiScaleLevel = ROLS.MaxLengthScale;
  }

  /**
   * Constructor for MSNode.
   *
   * @param n a {@link java.lang.String} object.
   * @param multiScaleLevel The multiscale level of this node.
   */
  public MSNode(String n, int multiScaleLevel) {
    this.name = n;
    this.MultiScaleLevel = multiScaleLevel;
  }

  /**
   * If <code>this</code> MSNode or any MSNode below it <code>equals</code> the argument, that MSNode
   * is returned.
   *
   * @param msNode a {@link ffx.potential.bonded.MSNode} object.
   * @return a {@link ffx.potential.bonded.MSNode} object.
   */
  public MSNode contains(MSNode msNode) {
    @SuppressWarnings("unchecked")
    Enumeration<TreeNode> e = depthFirstEnumeration();
    List<TreeNode> list = Collections.list(e);
    for (TreeNode node : list) {
      if (node.equals(msNode)) {
        return (MSNode) node;
      }
    }
    return null;
  }

  /**
   * destroy
   *
   * @return a boolean.
   */
  public boolean destroy() {
    if (getParent() != null) {
      removeFromParent();
    }
    name = null;
    selected = false;
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public void drawLabel(Canvas3D graphics, J3DGraphics2D g2d, Node node) {
    if (!isSelected()) {
      return;
    }
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode dataNode = (MSNode) e.nextElement();
      dataNode.drawLabel(graphics, g2d, node);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>MSNode equality := same class and same name. Consider replacing with a
   * Comparator&lt;MSNode&gt; for cases where non-reference equality is desired.
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    MSNode msNode = (MSNode) o;
    return Objects.equals(name, msNode.getName());
  }

  /**
   * Returns a List of all Angles below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<Angle> getAngleList() {
    return getList(Angle.class);
  }

  /**
   * Returns a List of all AngleTorsions below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<AngleTorsion> getAngleTorsionList() {
    return getList(AngleTorsion.class);
  }

  /**
   * Returns a List of all Atoms below the present MSNode.
   *
   * @return a new {@link java.util.List} object.
   */
  public List<Atom> getAtomList() {
    List<Atom> atomList = getList(Atom.class);
    Collections.sort(atomList);
    return atomList;
  }

  /**
   * getAtomList.
   *
   * @param originalOrder a boolean.
   * @return a {@link java.util.List} object.
   */
  public List<Atom> getAtomList(boolean originalOrder) {
    // As of now, for generic MSNode objects, atoms remain in their original
    // order. It is presently only a concern for MultiResidue.
    return getAtomList();
  }

  /**
   * Returns a List of all Bonds below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<Bond> getBondList() {
    return getList(Bond.class);
  }

  /** {@inheritDoc} */
  @Override
  public double[] getCenter(boolean w) {
    double[] Rc = {0, 0, 0};
    double sum = 0, mass = 1;
    List<Atom> atomList = getAtomList();
    for (Atom a : atomList) {
      if (w) {
        mass = a.getMass();
        sum += mass;
      }
      Rc[0] += mass * a.getX();
      Rc[1] += mass * a.getY();
      Rc[2] += mass * a.getZ();
    }
    if (!w) {
      sum = atomList.size();
    }
    for (int i = 0; i < 3; i++) {
      Rc[i] /= sum;
    }
    return Rc;
  }

  /**
   * Returns a List of the MSNode's Children (instead of using an Enumeration).
   *
   * @return a {@link java.util.List} object.
   */
  public List<MSNode> getChildList() {
    List<MSNode> l = new ArrayList<>();
    Enumeration<TreeNode> e = children();
    while (e.hasMoreElements()) {
      l.add((MSNode) e.nextElement());
    }
    return l;
  }

  /**
   * getExtent
   *
   * @return a double.
   */
  public double getExtent() {
    double extent = 0.0;
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      double temp = node.getExtent();
      if (temp > extent) {
        extent = temp;
      }
    }
    return extent;
  }

  /**
   * Returns a List of all ImproperTorsions below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<ImproperTorsion> getImproperTorsionList() {
    return getList(ImproperTorsion.class);
  }

  /** {@inheritDoc} */
  public <T extends TreeNode> List<T> getList(Class<T> c) {
    return getList(c, new ArrayList<>());
  }

  /** {@inheritDoc} */
  @Override
  public <T extends TreeNode> List<T> getList(Class<T> c, List<T> nodes) {
    if (c.isInstance(this)) {
      nodes.add(c.cast(this));
    }
    if (isLeaf() || !canBeChild(c)) {
      return nodes;
    }
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      node.getList(c, nodes);
    }
    return nodes;
  }

  /** {@inheritDoc} */
  @Override
  public <T extends TreeNode> long getMSCount(Class<T> c, long count) {
    if (c.isInstance(this)) {
      count++;
    }
    if (!canBeChild(c)) {
      return count;
    }
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      count += node.getMSCount(c, count);
    }
    return count;
  }

  /** {@inheritDoc} */
  @Override
  public <T extends TreeNode> T getMSNode(Class<T> c) {
    TreeNode[] nodes = getPath();
    for (TreeNode n : nodes) {
      if (c.isInstance(n)) {
        return c.cast(n);
      }
    }
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public double getMW() {
    double weight = 0.0;
    for (Atom atom : getAtomList()) {
      weight += atom.getMass();
    }
    return weight;
  }

  /**
   * Returns the name of this MSNode.
   *
   * @return a {@link java.lang.String} object.
   */
  public String getName() {
    return name;
  }

  /**
   * Sets the name of this NodeObject to n.
   *
   * @param n a {@link java.lang.String} object.
   */
  public void setName(String n) {
    name = n;
  }

  /**
   * Returns a List of all Out-of-Plane Bends below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<OutOfPlaneBend> getOutOfPlaneBendList() {
    return getList(OutOfPlaneBend.class);
  }

  /**
   * Returns a List of all Pi-Orbital Torsions below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<PiOrbitalTorsion> getPiOrbitalTorsionList() {
    return getList(PiOrbitalTorsion.class);
  }

  /**
   * Returns a List of all Stretch-Bends below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<StretchBend> getStretchBendList() {
    return getList(StretchBend.class);
  }

  /**
   * Returns a List of all StretchTorsions below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<StretchTorsion> getStretchTorsionList() {
    return getList(StretchTorsion.class);
  }

  /**
   * Returns a List of all Torsions below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<Torsion> getTorsionList() {
    return getList(Torsion.class);
  }

  /**
   * Returns a List of all Torsion-Torsions below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<TorsionTorsion> getTorsionTorsionList() {
    return getList(TorsionTorsion.class);
  }

  /**
   * Returns the total mass of all atoms in the MolecularAssembly, calculating the mass if it has not
   * already been done, defaulting to simple addition.
   *
   * @return Total mass of atoms in system.
   */
  public double getTotalMass() {
    if (totalMass == 0.0) {
      return getTotalMass(true, false);
    }
    return totalMass;
  }

  /**
   * Returns a List of all Urey-Bradleys below the present MSNode.
   *
   * @return a {@link java.util.List} object.
   */
  public List<UreyBradley> getUreyBradleyList() {
    return getList(UreyBradley.class);
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(name);
  }

  /**
   * isSelected
   *
   * @return a boolean.
   */
  public boolean isSelected() {
    return selected;
  }

  /**
   * Setter for the field <code>selected</code>.
   *
   * @param b a boolean.
   */
  public void setSelected(boolean b) {
    selected = b;
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      node.setSelected(b);
    }
  }

  /** Prints the MSNode's name */
  public void print() {
    System.out.println(name);
  }

  /**
   * removeChild.
   *
   * @param child a {@link ffx.potential.bonded.MSNode} object.
   */
  public void removeChild(MSNode child) {
    if (child != null && child.getParent() == this) {
      remove(child);
    }
  }

  /** {@inheritDoc} */
  @Override
  public void setColor(RendererCache.ColorModel colorModel, Color3f color, Material mat) {
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      node.setColor(colorModel, color, mat);
    }
  }

  /** {@inheritDoc} */
  @Override
  public void setView(RendererCache.ViewModel viewModel, List<BranchGroup> newShapes) {
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      node.setView(viewModel, newShapes);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Overridden toString method returns the MSNode's name
   */
  @Override
  public String toString() {
    return name;
  }

  /** {@inheritDoc} */
  @Override
  public void update() {
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode node = (MSNode) e.nextElement();
      node.update();
    }
  }

  /**
   * Calculates the total mass of all atoms in the MolecularAssembly, using either a simple addition
   * or the Kahan summation algorithm. The simple algorithm is a standard loop to add up the masses.
   *
   * @param recalculate Force recalculation
   * @param useKahan Use Kahan or simple addition algorithms
   * @return Total mass of all atoms in system.
   */
  private double getTotalMass(boolean recalculate, boolean useKahan) {
    if (recalculate) {
      List<Atom> atoms = getAtomList();
      if (atoms.isEmpty()) {
        totalMass = 0.0;
      } else if (useKahan) {
        totalMass = kahanSumMasses(atoms);
      } else {
        totalMass = sumMasses(atoms);
      }
    }
    return totalMass;
  }

  /**
   * Iterative summation of atomic masses.
   *
   * @param atoms List of atoms.
   * @return Mass of atoms.
   */
  private double sumMasses(List<Atom> atoms) {
    double sumMasses = 0.0;
    for (Atom atom : atoms) {
      sumMasses += atom.getMass();
    }
    return sumMasses;
  }

  /**
   * Implements the Kahan algorithm to very accurately sum the masses of all the atoms provided,
   * minimizing rounding error.
   *
   * @param atoms Atoms to sum the mass of.
   * @return Total mass.
   */
  private double kahanSumMasses(List<Atom> atoms) {
    double sum = 0.0;
    double comp = 0.0; // Running compensation
    for (Atom atom : atoms) {
      double atomMass = atom.getMass() - comp;
      double temp = sum + atomMass;
      comp = (temp - sum) - atomMass;
      sum = temp;
    }
    return sum;
  }

  /**
   * Returns true if Class c can be below this Class in the Hierarchy
   *
   * @param c Class
   * @return boolean
   */
  private <T extends TreeNode> boolean canBeChild(Class<T> c) {
    try {
      int multiScaleLevel = c.getDeclaredField("MultiScaleLevel").getInt(null);
      if (multiScaleLevel >= this.MultiScaleLevel) {
        return false;
      }
    } catch (NoSuchFieldException
             | SecurityException
             | IllegalArgumentException
             | IllegalAccessException e) {
      return true;
    }
    return true;
  }
}
