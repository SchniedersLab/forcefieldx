/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.potential.bonded;

import java.util.Enumeration;
import java.util.List;

import javax.media.j3d.BranchGroup;
import javax.swing.tree.TreeNode;

/**
 * The Joint class contains the geometry produced by the FGroup Joiner method.
 *
 * @author Michael J. Schnieders
 *
 */
public class Joint extends MSNode {

    private static final long serialVersionUID = 1L;
    /**
     * One of two Moieties forming this Joint
     */
    protected MSGroup group1 = null;
    protected MSGroup group2 = null;

    /**
     * Default Constructor
     */
    public Joint() {
        super("Joint");
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
     * @param piOrbitalTorsionNode a {@link ffx.potential.bonded.MSNode} object.
     * @param torsionTorsionNode a {@link ffx.potential.bonded.MSNode} object.
     */
    public Joint(MSGroup group1, MSGroup group2, MSNode bondNode,
            MSNode angleNode, MSNode stretchBendNode, MSNode ureyBradleyNode,
            MSNode outOfPlaneNode, MSNode torsionNode,
            MSNode piOrbitalTorsionNode, MSNode torsionTorsionNode) {
        super(group1.toString() + "  " + group2.toString());
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
        if (piOrbitalTorsionNode != null
                && piOrbitalTorsionNode.getChildCount() != 0) {
            add(piOrbitalTorsionNode);
        }
        if (torsionTorsionNode != null
                && torsionTorsionNode.getChildCount() != 0) {
            add(torsionTorsionNode);
        }
        refresh(null, null, null, null, null, null, null, null);
    }

    /**
     * <p>getAngles</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAngles() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof Angle) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getBonds</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getBonds() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof Bond) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getStretchBends</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getStretchBends() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof StretchBend) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getUreyBradleys</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getUreyBradleys() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof UreyBradley) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getOutOfPlaneBends</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getOutOfPlaneBends() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof OutOfPlaneBend) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getTorsions() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof Torsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getPiOrbitalTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getPiOrbitalTorsions() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof PiOrbitalTorsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>getTorsionTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getTorsionTorsions() {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof TorsionTorsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>merge</p>
     *
     * @param j a {@link ffx.potential.bonded.Joint} object.
     */
    public void merge(Joint j) {
        if (!((group1 == j.group1 && group2 == j.group2) || (group2 == j.group1 && group1 == j.group2))) {
            return;
        }
        refresh(j.getBonds(), j.getAngles(), j.getStretchBends(), j.getUreyBradleys(), j.getOutOfPlaneBends(), j.getTorsions(), j.getPiOrbitalTorsions(), j.getTorsionTorsions());
    }

    private void refresh(MSNode bonds, MSNode angles, MSNode stretchBends,
            MSNode ureyBradleys, MSNode outOfPlaneBends, MSNode torsions,
            MSNode piOrbitalTorsions, MSNode torsionTorsions) {
        for (Enumeration e = children(); e.hasMoreElements();) {
            MSNode jointChild = (MSNode) e.nextElement();
            if (jointChild.getChildCount() == 0) {
                jointChild.removeFromParent();
                continue;
            }
            MSNode node = (MSNode) jointChild.getChildAt(0);
            if (node instanceof Bond) {
                jointChild.setName("Bonds (" + jointChild.getChildCount() + ")");
                if (bonds != null) {
                    for (MSNode bond : bonds.getChildList()) {
                        jointChild.add(bond);
                    }
                }
                bonds = null;
            } else if (node instanceof Angle) {
                jointChild.setName("Angles (" + jointChild.getChildCount()
                        + ")");
                if (angles != null) {
                    for (MSNode angle : angles.getChildList()) {
                        jointChild.add(angle);
                    }
                }
                angles = null;
            } else if (node instanceof StretchBend) {
                jointChild.setName("Stretch-Bends ("
                        + jointChild.getChildCount() + ")");
                if (stretchBends != null) {
                    for (MSNode sb : stretchBends.getChildList()) {
                        jointChild.add(sb);
                    }
                }
                stretchBends = null;
            } else if (node instanceof UreyBradley) {
                jointChild.setName("Urey-Bradleys ("
                        + jointChild.getChildCount() + ")");
                if (ureyBradleys != null) {
                    for (MSNode ureyBradley : ureyBradleys.getChildList()) {
                        jointChild.add(ureyBradley);
                    }
                }
                ureyBradleys = null;
            } else if (node instanceof OutOfPlaneBend) {
                jointChild.setName("Out-of-Plane Bends ("
                        + jointChild.getChildCount() + ")");
                if (outOfPlaneBends != null) {
                    for (MSNode outOfPlaneBend : outOfPlaneBends.getChildList()) {
                        jointChild.add(outOfPlaneBend);
                    }
                }
                outOfPlaneBends = null;
            } else if (node instanceof Torsion) {
                jointChild.setName("Torsional Angles ("
                        + jointChild.getChildCount() + ")");
                if (torsions != null) {
                    for (MSNode torsion : torsions.getChildList()) {
                        jointChild.add(torsion);
                    }
                }
                torsions = null;
            } else if (node instanceof PiOrbitalTorsion) {
                jointChild.setName("Pi-Orbital Torsions ("
                        + jointChild.getChildCount() + ")");
                if (piOrbitalTorsions != null) {
                    for (MSNode piOrbitalTorsion : piOrbitalTorsions.getChildList()) {
                        jointChild.add(piOrbitalTorsion);
                    }
                }
                torsions = null;
            } else if (node instanceof TorsionTorsion) {
                jointChild.setName("Torsion-Torsions ("
                        + jointChild.getChildCount() + ")");
                if (torsionTorsions != null) {
                    for (MSNode torsionTorsion : torsionTorsions.getChildList()) {
                        jointChild.add(torsionTorsion);
                    }
                }
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
        if (piOrbitalTorsions != null) {
            add(piOrbitalTorsions);
        }
        if (torsionTorsions != null) {
            add(torsionTorsions);
        }
    }

    /**
     * <p>setColor</p>
     *
     * @param newColorModel a
     * {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void setColor(RendererCache.ColorModel newColorModel) {
        MSNode bonds = getBonds();
        if (bonds == null) {
            return;
        }
        for (Enumeration e = bonds.children(); e.hasMoreElements();) {
            Bond b = (Bond) e.nextElement();
            b.setColor(b.getAtom(0));
            b.setColor(b.getAtom(1));
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel newViewModel,
            List<BranchGroup> newShapes) {
        MSNode bonds = getBonds();
        if (bonds == null) {
            return;
        }
        for (Enumeration e = bonds.children(); e.hasMoreElements();) {
            Bond b = (Bond) e.nextElement();
            b.setView(newViewModel, newShapes);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Overidden toString method returns: "Joint: m1 Name - m2 Name"
     */
    @Override
    public String toString() {
        return getName();
    }
}
