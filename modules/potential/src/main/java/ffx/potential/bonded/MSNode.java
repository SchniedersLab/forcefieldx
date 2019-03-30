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
import javax.media.j3d.Canvas3D;
import javax.media.j3d.J3DGraphics2D;
import javax.media.j3d.Material;
import javax.media.j3d.Node;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;
import javax.vecmath.Color3f;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.ListIterator;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

/**
 * The MSNode class forms the basic unit that all data classes extend.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("CloneableImplementsClone")
public class MSNode extends DefaultMutableTreeNode implements ROLS {

    private static final long serialVersionUID = 1L;
    /**
     * The multiscale level of this node.
     */
    private final int MultiScaleLevel;
    /**
     * The name of this node.
     */
    private String name;
    /**
     * True if this node is selected.
     */
    protected boolean selected = false;
    /**
     * Total mass of this node and its children.
     */
    private double totalMass;

    /**
     * Default MSNode Constructor
     */
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
     * <p>
     * Constructor for MSNode.</p>
     *
     * @param n               a {@link java.lang.String} object.
     * @param multiScaleLevel a int.
     */
    public MSNode(String n, int multiScaleLevel) {
        this.name = n;
        this.MultiScaleLevel = multiScaleLevel;
    }

    /**
     * <p>
     * destroy</p>
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

    /**
     * Returns an ArrayList of all Atoms below the present MSNode.
     *
     * @return a new {@link java.util.ArrayList} object.
     */
    public ArrayList<Atom> getAtomList() {
        @SuppressWarnings("unchecked")
        Enumeration<TreeNode> e = depthFirstEnumeration();
        List<TreeNode> list = Collections.list(e);
        ArrayList<Atom> arrayList = new ArrayList<>();
        for (TreeNode node : list) {
            if (node instanceof Atom) {
                arrayList.add((Atom) node);
            }
        }

        Collections.sort(arrayList);
        return arrayList;
    }

    /**
     * <p>getAtomList.</p>
     *
     * @param originalOrder a boolean.
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Atom> getAtomList(boolean originalOrder) {
        // As of now, for generic MSNode objects, atoms remain in their original
        // order. It is presently only a concern for MultiResidue.
        return getAtomList();
    }

    /**
     * If <code>this</code> MSNode or any MSNode below it <code>equals</code>
     * the argument, that MSNode is returned.
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
     * Returns an ArrayList of all Bonds below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getBondList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(Bond.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Angles below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getAngleList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(Angle.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Stretch-Bends below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getStretchBendList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(StretchBend.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Urey-Bradleys below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getUreyBradleyList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(UreyBradley.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Out-of-Plane Bends below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getOutOfPlaneBendList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(OutOfPlaneBend.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Torsions below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getTorsionList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(Torsion.class, arrayList);
    }

    /**
     * Returns an ArrayList of all AngleTorsions below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getAngleTorsionList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(AngleTorsion.class, arrayList);
    }

    /**
     * Returns an ArrayList of all StretchTorsions below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getStretchTorsionList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(StretchTorsion.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Pi-Orbital Torsions below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getPiOrbitalTorsionList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(PiOrbitalTorsion.class, arrayList);
    }

    /**
     * Returns an ArrayList of all Torsion-Torsions below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getTorsionTorsionList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(TorsionTorsion.class, arrayList);
    }

    /**
     * Returns an ArrayList of all ImproperTorsions below the present MSNode.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<ROLS> getImproperTorsionList() {
        ArrayList<ROLS> arrayList = new ArrayList<>();
        return getList(ImproperTorsion.class, arrayList);
    }

    /**
     * Returns an ArrayList of the MSNode's Children (instead of using an
     * Enumeration).
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<MSNode> getChildList() {
        ArrayList<MSNode> l = new ArrayList<>();
        Enumeration<?> e = children();
        while (e.hasMoreElements()) {
            l.add((MSNode) e.nextElement());
        }
        return l;
    }

    /**
     * Returns a ListIterator containing this MSNode's children.
     *
     * @return a {@link java.util.ListIterator} object.
     */
    ListIterator<MSNode> getChildListIterator() {
        return getChildList().listIterator();
    }

    /**
     * <p>
     * getExtent</p>
     *
     * @return a double.
     */
    public double getExtent() {
        double extent = 0.0;
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            double temp = node.getExtent();
            if (temp > extent) {
                extent = temp;
            }
        }
        return extent;
    }

    /**
     * Request all descendants of the given type; returned list will automatically
     * conform to any superclass thereof.
     *
     * @param c   a {@link java.lang.Class} object.
     * @param <U> a U object.
     * @param <T> a T object.
     * @return a {@link java.util.List} object.
     */
    public <U extends MSNode, T extends U> List<U> getDescendants(Class<T> c) {
        List<U> nodes = new ArrayList<>();
        castDescendants(c, nodes);
        return nodes;
    }

    private <U extends MSNode, T extends U> void castDescendants(Class<T> clas, List<U> nodes) {
        if (clas.isInstance(this)) {
            nodes.add(clas.cast(this));
        }
        if (isLeaf()) {
            return;
        }
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            node.castDescendants(clas, nodes);
        }
    }

    /**
     * <p>
     * getMultiScaleLevel</p>
     *
     * @return a int.
     */
    public int getMultiScaleLevel() {
        return MultiScaleLevel;
    }

    /**
     * Returns the total mass of all atoms in the MolecularAssembly, calculating
     * the mass if it has not already been done, defaulting to simple addition.
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
     * Returns the name of this MSNode.
     *
     * @return a {@link java.lang.String} object.
     */
    public String getName() {
        return name;
    }

    /**
     * <p>
     * isSelected</p>
     *
     * @return a boolean.
     */
    public boolean isSelected() {
        return selected;
    }

    /**
     * Prints the MSNode's name
     */
    public void print() {
        System.out.println(name);
    }

    /**
     * <p>removeChild.</p>
     *
     * @param child a {@link ffx.potential.bonded.MSNode} object.
     */
    public void removeChild(MSNode child) {
        if (child != null && child.getParent() == this) {
            remove(child);
        }
    }

    /**
     * Sets the name of this NodeObect to n.
     *
     * @param n a {@link java.lang.String} object.
     */
    public void setName(String n) {
        name = n;
    }

    /**
     * <p>
     * Setter for the field <code>selected</code>.</p>
     *
     * @param b a boolean.
     */
    public void setSelected(boolean b) {
        selected = b;
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            node.setSelected(b);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public long getMSCount(Class<?> c, long count) {
        if (c.isInstance(this)) {
            count++;
        }
        if (!canBeChild(c)) {
            return count;
        }
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            count += node.getMSCount(c, count);
        }
        return count;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ROLS getMSNode(Class<?> c) {
        TreeNode[] nodes = getPath();
        for (TreeNode n : nodes) {
            if (c.isInstance(n)) {
                return (ROLS) n;
            }
        }
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ArrayList<ROLS> getList(Class<?> c, ArrayList<ROLS> nodes) {
        if (c.isInstance(this)) {
            nodes.add(this);
        }
        if (isLeaf() || !canBeChild(c)) {
            return nodes;
        }
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            ROLS node = (ROLS) e.nextElement();
            node.getList(c, nodes);
        }
        return nodes;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCenter(boolean w) {
        double[] Rc = {0, 0, 0};
        double sum = 0, mass = 1;
        ArrayList<Atom> atomList = getAtomList();
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
     * {@inheritDoc}
     */
    @Override
    public double getMW() {
        double weight = 0.0;
        for (Atom atom : getAtomList()) {
            weight += atom.getMass();
        }
        return weight;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel colorModel, Color3f color,
                         Material mat) {
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            node.setColor(colorModel, color, mat);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel viewModel,
                        List<BranchGroup> newShapes) {
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            node.setView(viewModel, newShapes);
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overidden toString method returns the MSNode's name
     */
    @Override
    public String toString() {
        return name;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update() {
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            MSNode node = (MSNode) e.nextElement();
            node.update();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void drawLabel(Canvas3D graphics, J3DGraphics2D g2d, Node node) {
        if (!isSelected()) {
            return;
        }
        MSNode dataNode;
        for (Enumeration<?> e = children(); e.hasMoreElements(); ) {
            dataNode = (MSNode) e.nextElement();
            dataNode.drawLabel(graphics, g2d, node);
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * MSNode equality := same class and same name.
     * Consider replacing with a Comparator&lt;MSNode&gt; for cases where
     * non-reference equality is desired.
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        MSNode other = (MSNode) object;
        if (name == null && other.getName() == null) {
            return true;
        } else if (name != null && other.getName() != null) {
            return name.equals(other.getName());
        } else {
            return false;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        if (name != null) {
            return hash(SEED, name.hashCode());
        } else {
            return SEED;
        }
    }

    /**
     * Calculates the total mass of all atoms in the MolecularAssembly, using
     * either a simple addition or the Kahan summation algorithm. The simple
     * algorithm is a standard loop to add up the masses.
     *
     * @param recalculate Force recalculation
     * @param useKahan    Use Kahan or simple addition algorithms
     * @return Total mass of all atoms in system.
     */
    private double getTotalMass(boolean recalculate, boolean useKahan) {
        if (recalculate) {
            ArrayList<Atom> atoms = getAtomList();
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
     * Implements the Kahan algorithm to very accurately sum the masses of all
     * the atoms provided, minimizing rounding error.
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
    private boolean canBeChild(Class<?> c) {
        try {
            int multiScaleLevel = c.getDeclaredField("MultiScaleLevel").getInt(null);
            if (multiScaleLevel >= this.MultiScaleLevel) {
                return false;
            }
        } catch (NoSuchFieldException | SecurityException | IllegalArgumentException | IllegalAccessException e) {
            return true;
        }
        return true;
    }

}
