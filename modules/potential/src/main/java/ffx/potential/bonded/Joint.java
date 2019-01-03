/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.bonded;

import javax.media.j3d.BranchGroup;
import javax.swing.tree.TreeNode;

import java.util.Enumeration;
import java.util.List;

/**
 * The Joint class contains the geometry produced by the FGroup Joiner method.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Joint extends MSNode {

    private static final long serialVersionUID = 1L;

    /**
     * First group forming this Joint
     */
    protected MSGroup group1 = null;
    /**
     * Second group forming this Joint
     */
    protected MSGroup group2 = null;

    /**
     * Default Constructor
     */
    public Joint() {
        super("Joint");
        group1 = group2 = null;
    }

    /**
     * <p>Constructor for Joint.</p>
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
     * @param group1               a {@link ffx.potential.bonded.MSGroup} object.
     * @param group2               a {@link ffx.potential.bonded.MSGroup} object.
     * @param bondNode             a {@link ffx.potential.bonded.MSNode} object.
     * @param angleNode            a {@link ffx.potential.bonded.MSNode} object.
     * @param stretchBendNode      a {@link ffx.potential.bonded.MSNode} object.
     * @param ureyBradleyNode      a {@link ffx.potential.bonded.MSNode} object.
     * @param outOfPlaneNode       a {@link ffx.potential.bonded.MSNode} object.
     * @param torsionNode          a {@link ffx.potential.bonded.MSNode} object.
     * @param stretchTorsionNode   a {@link ffx.potential.bonded.MSNode} object.
     * @param angleTorsionNode     a {@link ffx.potential.bonded.MSNode} object.
     * @param piOrbitalTorsionNode a {@link ffx.potential.bonded.MSNode} object.
     * @param torsionTorsionNode   a {@link ffx.potential.bonded.MSNode} object.
     */
    public Joint(MSGroup group1, MSGroup group2, MSNode bondNode,
                 MSNode angleNode, MSNode stretchBendNode, MSNode ureyBradleyNode,
                 MSNode outOfPlaneNode, MSNode torsionNode,
                 MSNode stretchTorsionNode, MSNode angleTorsionNode,
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
        if (stretchTorsionNode != null && stretchTorsionNode.getChildCount() != 0) {
            add(stretchTorsionNode);
        }
        if (angleTorsionNode != null && angleTorsionNode.getChildCount() != 0) {
            add(angleTorsionNode);
        }
        if (piOrbitalTorsionNode != null
                && piOrbitalTorsionNode.getChildCount() != 0) {
            add(piOrbitalTorsionNode);
        }
        if (torsionTorsionNode != null
                && torsionTorsionNode.getChildCount() != 0) {
            add(torsionTorsionNode);
        }
        refresh(null, null, null, null, null,
                null, null, null, null, null);
    }

    /**
     * <p>
     * getAngles</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAngles() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof Angle) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getBonds</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getBonds() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof Bond) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getStretchBends</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getStretchBends() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof StretchBend) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getUreyBradleys</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getUreyBradleys() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof UreyBradley) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getOutOfPlaneBends</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getOutOfPlaneBends() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof OutOfPlaneBend) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getTorsions() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof Torsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getStretchTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getStretchTorsions() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof StretchTorsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getAngleTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAngleTorsions() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof AngleTorsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getPiOrbitalTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getPiOrbitalTorsions() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof PiOrbitalTorsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * getTorsionTorsions</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getTorsionTorsions() {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
            MSNode m = (MSNode) e.nextElement();
            TreeNode node = m.getChildAt(0);
            if (node instanceof TorsionTorsion) {
                return m;
            }
        }
        return null;
    }

    /**
     * <p>
     * merge</p>
     *
     * @param j a {@link ffx.potential.bonded.Joint} object.
     */
    public void merge(Joint j) {
        if (!((group1 == j.group1 && group2 == j.group2) || (group2 == j.group1 && group1 == j.group2))) {
            return;
        }
        refresh(j.getBonds(), j.getAngles(), j.getStretchBends(), j.getUreyBradleys(), j.getOutOfPlaneBends(),
                j.getTorsions(), j.getStretchTorsions(), j.getAngleTorsions(), j.getPiOrbitalTorsions(), j.getTorsionTorsions());
    }

    private void refresh(MSNode bonds, MSNode angles, MSNode stretchBends,
                         MSNode ureyBradleys, MSNode outOfPlaneBends,
                         MSNode torsions, MSNode stretchTorsions, MSNode angleTorsions,
                         MSNode piOrbitalTorsions, MSNode torsionTorsions) {
        for (Enumeration e = children(); e.hasMoreElements(); ) {
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
            } else if (node instanceof StretchTorsion) {
                jointChild.setName("Stretch-Torsions ("
                        + jointChild.getChildCount() + ")");
                if (torsions != null) {
                    for (MSNode stretchTorsion : stretchTorsions.getChildList()) {
                        jointChild.add(stretchTorsion);
                    }
                }
                stretchTorsions = null;
            } else if (node instanceof AngleTorsion) {
                jointChild.setName("Angle-Torsions ("
                        + jointChild.getChildCount() + ")");
                if (torsions != null) {
                    for (MSNode angleTorsion : angleTorsions.getChildList()) {
                        jointChild.add(angleTorsion);
                    }
                }
                angleTorsions = null;
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
     * <p>assignBonds.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     */
    public void assignBonds(Atom atom) {
        for (ROLS bond : getBondList()) {
            Bond b = (Bond) bond;
            if (b.containsAtom(atom)) {
                atom.setBond(b);
            }
        }
    }

    /**
     * <p>assignAngles.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     */
    public void assignAngles(Atom atom) {
        for (ROLS angle : getAngleList()) {
            Angle a = (Angle) angle;
            if (a.containsAtom(atom)) {
                atom.setAngle(a);
            }
        }
    }

    /**
     * <p>assignTorsions.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     */
    public void assignTorsions(Atom atom) {
        for (ROLS torsion : getTorsionList()) {
            Torsion t = (Torsion) torsion;
            if (t.containsAtom(atom)) {
                atom.setTorsion(t);
            }
        }
    }

    /**
     * <p>assignReferences.</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     */
    public void assignReferences(Atom atom) {
        assignBonds(atom);
        assignAngles(atom);
        assignTorsions(atom);
    }

    /**
     * <p>
     * setColor</p>
     *
     * @param newColorModel a
     *                      {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void setColor(RendererCache.ColorModel newColorModel) {
        MSNode bonds = getBonds();
        if (bonds == null) {
            return;
        }
        for (Enumeration e = bonds.children(); e.hasMoreElements(); ) {
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
        for (Enumeration e = bonds.children(); e.hasMoreElements(); ) {
            Bond b = (Bond) e.nextElement();
            b.setView(newViewModel, newShapes);
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overidden toString method returns: "Joint: m1 Name - m2 Name"
     */
    @Override
    public String toString() {
        return getName();
    }
}
