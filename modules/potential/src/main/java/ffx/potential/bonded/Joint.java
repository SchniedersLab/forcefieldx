// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import org.jogamp.java3d.BranchGroup;

import javax.swing.tree.TreeNode;
import java.io.Serial;
import java.util.Enumeration;
import java.util.List;

/**
 * The Joint class contains the geometry produced by the FGroup Joiner method.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Joint extends MSNode {

  @Serial
  private static final long serialVersionUID = 1L;

  /** First group forming this Joint */
  private final MSGroup group1;
  /** Second group forming this Joint */
  private final MSGroup group2;

  /** Default Constructor */
  public Joint() {
    super("Joint");
    group1 = group2 = null;
  }

  /**
   * Constructor for Joint.
   *
   * @param name a {@link java.lang.String} object.
   */
  public Joint(String name) {
    super(name);
    group1 = group2 = null;
  }

  /**
   * Constructs a Joint between Group 1 and Group 2.
   *
   * @param group1 a {@link ffx.potential.bonded.MSGroup} object.
   * @param group2 a {@link ffx.potential.bonded.MSGroup} object.
   * @param bondNode a {@link ffx.potential.bonded.MSNode} object.
   * @param angleNode a {@link ffx.potential.bonded.MSNode} object.
   * @param stretchBendNode a {@link ffx.potential.bonded.MSNode} object.
   * @param ureyBradleyNode a {@link ffx.potential.bonded.MSNode} object.
   * @param outOfPlaneNode a {@link ffx.potential.bonded.MSNode} object.
   * @param torsionNode a {@link ffx.potential.bonded.MSNode} object.
   * @param improperTorsionNode a {@link ffx.potential.bonded.MSNode} object.
   * @param stretchTorsionNode a {@link ffx.potential.bonded.MSNode} object.
   * @param angleTorsionNode a {@link ffx.potential.bonded.MSNode} object.
   * @param piOrbitalTorsionNode a {@link ffx.potential.bonded.MSNode} object.
   * @param torsionTorsionNode a {@link ffx.potential.bonded.MSNode} object.
   */
  public Joint(MSGroup group1, MSGroup group2, MSNode bondNode, MSNode angleNode,
      MSNode stretchBendNode, MSNode ureyBradleyNode, MSNode outOfPlaneNode, MSNode torsionNode,
      MSNode improperTorsionNode, MSNode stretchTorsionNode, MSNode angleTorsionNode,
      MSNode piOrbitalTorsionNode, MSNode torsionTorsionNode) {
    super(group1 + "  " + group2);
    this.group1 = group1;
    this.group2 = group2;
    if (bondNode != null && bondNode.getChildCount() != 0) {
      add(bondNode);
    }
    if (angleNode != null && angleNode.getChildCount() != 0) {
      add(angleNode);
    }
    if (stretchBendNode != null && stretchBendNode.getChildCount() != 0) {
      add(stretchBendNode);
    }
    if (ureyBradleyNode != null && ureyBradleyNode.getChildCount() != 0) {
      add(ureyBradleyNode);
    }
    if (outOfPlaneNode != null && outOfPlaneNode.getChildCount() != 0) {
      add(outOfPlaneNode);
    }
    if (torsionNode != null && torsionNode.getChildCount() != 0) {
      add(torsionNode);
    }
    if (improperTorsionNode != null && improperTorsionNode.getChildCount() != 0) {
      add(improperTorsionNode);
    }
    if (stretchTorsionNode != null && stretchTorsionNode.getChildCount() != 0) {
      add(stretchTorsionNode);
    }
    if (angleTorsionNode != null && angleTorsionNode.getChildCount() != 0) {
      add(angleTorsionNode);
    }
    if (piOrbitalTorsionNode != null && piOrbitalTorsionNode.getChildCount() != 0) {
      add(piOrbitalTorsionNode);
    }
    if (torsionTorsionNode != null && torsionTorsionNode.getChildCount() != 0) {
      add(torsionTorsionNode);
    }
    refresh(null, null, null, null, null,
        null, null, null, null, null, null);
  }

  /**
   * merge
   *
   * @param j a {@link ffx.potential.bonded.Joint} object.
   */
  public void merge(Joint j) {
    if (!((group1 == j.group1 && group2 == j.group2)
        || (group2 == j.group1 && group1 == j.group2))) {
      return;
    }
    refresh(j.getNode(Bond.class), j.getNode(Angle.class), j.getNode(StretchBend.class),
        j.getNode(UreyBradley.class), j.getNode(OutOfPlaneBend.class), j.getNode(Torsion.class),
        j.getNode(ImproperTorsion.class), j.getNode(StretchTorsion.class),
        j.getNode(AngleTorsion.class), j.getNode(PiOrbitalTorsion.class),
        j.getNode(TorsionTorsion.class));
  }

  /**
   * setColor
   *
   * @param newColorModel a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
   */
  public void setColor(RendererCache.ColorModel newColorModel) {
    MSNode bonds = getNode(Bond.class);
    if (bonds == null) {
      return;
    }
    for (Enumeration<TreeNode> e = bonds.children(); e.hasMoreElements(); ) {
      Bond b = (Bond) e.nextElement();
      b.setColor(b.getAtom(0));
      b.setColor(b.getAtom(1));
    }
  }

  /** {@inheritDoc} */
  @Override
  public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
    MSNode bonds = getNode(Bond.class);
    if (bonds == null) {
      return;
    }
    for (Enumeration<TreeNode> e = bonds.children(); e.hasMoreElements(); ) {
      Bond b = (Bond) e.nextElement();
      b.setView(newViewModel, newShapes);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Overridden toString method returns: "Joint: m1 Name - m2 Name"
   */
  @Override
  public String toString() {
    return getName();
  }

  /**
   * Refresh the Joint.
   *
   * @param bonds Bond node.
   * @param angles Angle node.
   * @param stretchBends StretchBend node.
   * @param ureyBradleys UreyBradley node.
   * @param outOfPlaneBends OutOfPlaneBend node.
   * @param torsions Torsion node.
   * @param improperTorsions ImproperTorsion node.
   * @param stretchTorsions StretchTorsion node.
   * @param angleTorsions AngleTorsion node.
   * @param piOrbitalTorsions PiOrbitalTorsion node.
   * @param torsionTorsions TorsionTorsion node.
   */
  private void refresh(MSNode bonds, MSNode angles, MSNode stretchBends, MSNode ureyBradleys,
      MSNode outOfPlaneBends, MSNode torsions, MSNode improperTorsions, MSNode stretchTorsions,
      MSNode angleTorsions, MSNode piOrbitalTorsions, MSNode torsionTorsions) {

    // Loop over all children of the Joint.
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode jointChild = (MSNode) e.nextElement();
      if (jointChild.getChildCount() == 0) {
        jointChild.removeFromParent();
        continue;
      }
      MSNode node = (MSNode) jointChild.getChildAt(0);

      if (node instanceof Bond) {
        addAll("Bonds", jointChild, bonds);
        bonds = null;
      } else if (node instanceof Angle) {
        addAll("Angles", jointChild, angles);
        angles = null;
      } else if (node instanceof StretchBend) {
        addAll("Stretch-Bends", jointChild, stretchBends);
        stretchBends = null;
      } else if (node instanceof UreyBradley) {
        addAll("Urey-Bradleys", jointChild, ureyBradleys);
        ureyBradleys = null;
      } else if (node instanceof OutOfPlaneBend) {
        addAll("Out-of-Plane Bends", jointChild, outOfPlaneBends);
        outOfPlaneBends = null;
      } else if (node instanceof Torsion) {
        addAll("Torsional Angles", jointChild, torsions);
        torsions = null;
      } else if (node instanceof ImproperTorsion) {
        addAll("Improper Torsions", jointChild, improperTorsions);
        improperTorsions = null;
      } else if (node instanceof StretchTorsion) {
        addAll("Stretch-Torsions", jointChild, stretchTorsions);
        stretchTorsions = null;
      } else if (node instanceof AngleTorsion) {
        addAll("Angle-Torsions", jointChild, angleTorsions);
        angleTorsions = null;
      } else if (node instanceof PiOrbitalTorsion) {
        addAll("Pi-Orbital Torsions", jointChild, piOrbitalTorsions);
        piOrbitalTorsions = null;
      } else if (node instanceof TorsionTorsion) {
        addAll("Torsion-Torsions", jointChild, torsionTorsions);
        torsionTorsions = null;
      }
    }

    if (bonds != null) {
      add(bonds);
    }
    if (angles != null) {
      add(angles);
    }
    if (stretchBends != null) {
      add(stretchBends);
    }
    if (ureyBradleys != null) {
      add(ureyBradleys);
    }
    if (outOfPlaneBends != null) {
      add(outOfPlaneBends);
    }
    if (torsions != null) {
      add(torsions);
    }
    if (improperTorsions != null) {
      add(improperTorsions);
    }
    if (stretchTorsions != null) {
      add(stretchTorsions);
    }
    if (angleTorsions != null) {
      add(angleTorsions);
    }
    if (piOrbitalTorsions != null) {
      add(piOrbitalTorsions);
    }
    if (torsionTorsions != null) {
      add(torsionTorsions);
    }
  }

  /**
   * Add all the children of the node to the parent.
   *
   * @param name Name of the parent.
   * @param parent Parent node.
   * @param node Node to add.
   */
  private void addAll(String name, MSNode parent, MSNode node) {
    if (node == null) {
      return;
    }
    for (MSNode bond : node.getChildList()) {
      parent.add(bond);
    }
    parent.setName(name + " (" + parent.getChildCount() + ")");
  }

  /**
   * Get the Joint node for the specified class.
   *
   * @param clazz The class of the bonded term to search for.
   * @return The Joint node for the specified class.
   */
  private MSNode getNode(Class<? extends BondedTerm> clazz) {
    for (Enumeration<TreeNode> e = children(); e.hasMoreElements(); ) {
      MSNode m = (MSNode) e.nextElement();
      TreeNode node = m.getChildAt(0);
      if (node.getClass().equals(clazz)) {
        return m;
      }
    }
    return null;
  }

}
