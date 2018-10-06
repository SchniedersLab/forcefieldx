/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import ffx.numerics.math.VectorMath;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.parameters.ForceField;
import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

/**
 * The Polymer class encapsulates a peptide or nucleotide chain.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Polymer extends MSGroup {

    private static final long serialVersionUID = 1L;
    /**
     * Constant <code>MultiScaleLevel=3</code>
     */
    public static final int MultiScaleLevel = 3;
    private static int count = 0;
    private static double[] da = new double[3];
    private static double[] db = new double[3];
    /**
     * Constant <code>polymerColor</code>
     */
    public final static Map<Integer, Color3f> polymerColor = new HashMap<>();

    static {
        polymerColor.put(0, RendererCache.RED);
        polymerColor.put(1, RendererCache.ORANGE);
        polymerColor.put(2, RendererCache.YELLOW);
        polymerColor.put(3, RendererCache.GREEN);
        polymerColor.put(4, RendererCache.BLUE);
        polymerColor.put(5, RendererCache.MAGENTA);
        polymerColor.put(6, RendererCache.CYAN);
        polymerColor.put(7, RendererCache.WHITE);
        polymerColor.put(8, RendererCache.GRAY);
        polymerColor.put(9, RendererCache.PINK);
    }

    private boolean link = false;
    private int polymerNumber;
    private Character chainID;

    /**
     * Polymer constructor.
     *
     * @param chainID Possibly redundant PDB chainID.
     * @param segID   Unique identifier from A-Z,0-9, then 1A-1Z,10-19, etc.
     */
    public Polymer(Character chainID, String segID) {
        super(segID);
        this.chainID = chainID;
        this.polymerNumber = ++count;
    }

    /**
     * Polymer constructor.
     *
     * @param chainID Possibly redundant PDB chainID.
     * @param segID   Unique identifier from A-Z,0-9, then 1A-1Z,10-19, etc.
     * @param link    a boolean.
     */
    public Polymer(Character chainID, String segID, boolean link) {
        this(chainID, segID);
        this.link = link;
    }

    /**
     * Polymer Constructor.
     *
     * @param segID    A unique identifier from A-Z,0-9, then 1A-1Z,10-19, etc.
     * @param residues Represents a MSNode where the Polymer's residues have
     *                 been attached.
     * @param chainID  a {@link java.lang.Character} object.
     */
    public Polymer(Character chainID, String segID, MSNode residues) {
        super(segID, residues);
        this.chainID = chainID;
        polymerNumber = ++count;
    }

    /**
     * {@inheritDoc}
     * <p>
     * A generic method for adding a MSNode to the Polymer.
     */
    @Override
    public MSNode addMSNode(MSNode msNode) {
        assert (msNode instanceof Residue);

        Residue residue = (Residue) msNode;
        int resNumber = residue.getResidueNumber();

        MSNode residueNode = getAtomNode();
        int n = residueNode.getChildCount();
        int childIndex = n;

        for (int i = 0; i < n; i++) {
            Residue current = (Residue) residueNode.getChildAt(i);
            if (current.getResidueNumber() > resNumber) {
                childIndex = i;
                break;
            }
        }

        getAtomNode().insert(residue, childIndex);
        return msNode;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Joiner joins Moieties m1 and m2 and returns the Geometry objects formed
     * in a Joint.
     */
    public Joint createJoint(Residue residue1, Residue residue2, ForceField forceField) {
        Joint joint = null;
        for (Enumeration e = residue1.getAtomNode().children(); e.hasMoreElements(); ) {
            Atom a1 = (Atom) e.nextElement();
            a1.getXYZ(da);
            for (Enumeration e2 = residue2.getAtomNode().children(); e2.hasMoreElements(); ) {
                Atom a2 = (Atom) e2.nextElement();
                a2.getXYZ(db);
                double d1 = VectorMath.dist(da, db);
                double d2 = Bond.BUFF + a1.getVDWR() / 2 + a2.getVDWR() / 2;
                if (d1 < d2) {
                    Bond b = new Bond(a1, a2);
                    Joint newJoint = createJoint(b, residue1, residue2, forceField);
                    if (joint != null) {
                        joint.merge(newJoint);
                    } else {
                        joint = newJoint;
                    }
                }
            }
        }
        return joint;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overidden equals method.
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        Polymer other = (Polymer) object;
        return getName().equals(other.getName());
    }

    /**
     * {@inheritDoc}
     * <p>
     * Finalize should be called after all the Residues have been added to the
     * Polymer. This method in turn calls the Finalize method of each Residue,
     * then forms Joints between adjacent Residues in the Polymer
     */
    @Override
    public void finalize(boolean finalizeGroups, ForceField forceField) {
        ListIterator li;
        ArrayList residues = getAtomNodeList();
        setFinalized(false);

        // Finalize the residues in the Polymer
        if (finalizeGroups) {
            for (li = residues.listIterator(); li.hasNext(); ) {
                ((Residue) li.next()).finalize(true, forceField);
            }
        }
        // Join the residues in the Polymer
        if (link) {
            Residue residue = getFirstResidue();
            if (residue.residueType == ResidueType.AA) {
                getAtomNode().setName("Amino Acids " + "(" + residues.size() + ")");
            } else if (residue.residueType == ResidueType.NA) {
                getAtomNode().setName("Nucleic Acids " + "(" + residues.size() + ")");
            } else {
                getAtomNode().setName("Residues " + "(" + residues.size() + ")");
            }
            Joint j;
            MSNode joints = getTermNode();
            joints.removeAllChildren();
            List<Atom> atoms = getAtomList();

            for (Atom a : atoms) {
                if (a.getNumBonds() > 0) {
                    for (Bond b : a.getBonds()) {
                        if (!b.sameGroup() && b.getParent() == null) {
                            Residue r1 = (Residue) a.getMSNode(Residue.class);
                            Residue r2 = (Residue) b.get1_2(a).getMSNode(
                                    Residue.class);
                            j = createJoint(b, r1, r2, forceField);
                            joints.add(j);
                        }
                    }
                }
            }

            if (residue != null) {
                if (residue.residueType == ResidueType.AA) {
                    getTermNode().setName(
                            "Peptide Bonds " + "(" + joints.getChildCount() + ")");
                } else {
                    getTermNode().setName(
                            "Linkages " + "(" + joints.getChildCount() + ")");
                }
            } else {
                getTermNode().setName(
                        "Linkages " + "(" + joints.getChildCount() + ")");
            }
        } else {
            getAtomNode().setName("Sub-Groups " + "(" + residues.size() + ")");
            if (getTermNode().getParent() != null) {
                removeChild(getTermNode());
            }
        }
        removeLeaves();
        setFinalized(true);
    }

    /**
     * <p>
     * Getter for the field <code>link</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getLink() {
        return link;
    }

    /**
     * <p>
     * Getter for the field <code>chainID</code>.</p>
     *
     * @return a {@link java.lang.Character} object.
     */
    public Character getChainID() {
        return chainID;
    }

    /**
     * <p>addMultiResidue.</p>
     *
     * @param multiResidue a {@link ffx.potential.bonded.MultiResidue} object.
     */
    public void addMultiResidue(MultiResidue multiResidue) {
        Residue residue = multiResidue.getActive();
        MSNode residueNode = getAtomNode();
        int index = residueNode.getIndex(residue);
        if (index < 0) {
            System.err.println("WARNING!  Polymer::addMultiResidue did not find a corresponding Residue on Polymer.");
            residueNode.add(multiResidue);
        } else {
            residueNode.remove(index);
            residueNode.insert(multiResidue, index);
        }
        multiResidue.add(residue);
    }

    /**
     * <p>addMultiTerminus.</p>
     *
     * @param residue       a {@link ffx.potential.bonded.Residue} object.
     * @param multiTerminus a {@link ffx.potential.bonded.MultiTerminus} object.
     */
    public void addMultiTerminus(Residue residue, MultiTerminus multiTerminus) {
        ArrayList<MSNode> children = residue.getChildList();
        for (MSNode child : children) {
            multiTerminus.add(child);
        }
        MSNode residueNode = getAtomNode();
        int index = residueNode.getIndex(residue);
        residueNode.remove(index);
        residueNode.insert(multiTerminus, index);
    }

    /**
     * TODO: Was the sole hook on BondedTerm equality definition via getID();
     * will rewrite with a simple Comparator soon.
     *
     * @return An ArrayList of Dihedral objects representing the Phi/Psi angles
     * of the Polymer, useful for creating Ramachandran plots
     */
    public List<ArrayList<Torsion>> getPhiPsiList() {
        MSNode dihedrals;
        ListIterator<MSNode> li, lj;
        List<ArrayList<Torsion>> phipsi = new ArrayList<>();
        ArrayList<Torsion> phi = new ArrayList<>();
        ArrayList<Torsion> psi = new ArrayList<>();
        phipsi.add(phi);
        phipsi.add(psi);
        MSNode joints = getTermNode();
        for (li = joints.getChildListIterator(); li.hasNext(); ) {
            dihedrals = ((Joint) li.next()).getTorsions();
            for (lj = dihedrals.getChildListIterator(); lj.hasNext(); ) {
                Torsion d = (Torsion) lj.next();
                String s = d.getID();
                // Phi
                if (s.equals("C-N-CA-C") || s.equals("C-CA-N-C")) {
                    phi.add(d);
                } // Psi
                else if (s.equals("N-C-CA-N") || s.equals("N-CA-C-N")) {
                    psi.add(d);
                }
            }
        }
        return phipsi;
    }

    /**
     * <p>
     * getFirstResidue</p>
     *
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public Residue getFirstResidue() {
        MSNode atomNode = getAtomNode();
        if (atomNode == null) {
            return null;
        }
        return (Residue) atomNode.getChildAt(0);
    }

    /**
     * <p>
     * getResidues</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Residue> getResidues() {
        ArrayList<Residue> residues = new ArrayList<Residue>();
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements(); ) {
            Residue r = (Residue) e.nextElement();
            residues.add(r);
        }
        return residues;
    }

    /**
     * <p>
     * getResidue</p>
     *
     * @param resNum a int.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public Residue getResidue(int resNum) {
        if (resNum > 0 && getAtomNode().getChildCount() >= resNum) {
            Residue r = (Residue) getAtomNode().getChildAt(resNum - 1);
            if (r.getResidueNumber() == resNum) {
                return r;
            }
        }
        // Fall back for non-ordered children
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements(); ) {
            Residue r = (Residue) e.nextElement();
            if (r.getResidueNumber() == resNum) {
                return r;
            }
        }
        return null;
    }

    /**
     * <p>
     * getResidue</p>
     *
     * @param resName a {@link java.lang.String} object.
     * @param resNum  a int.
     * @param create  a boolean.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public Residue getResidue(String resName, int resNum, boolean create) {
        return getResidue(resName, resNum, create, Residue.ResidueType.UNK);
    }

    /**
     * <p>
     * getResidue</p>
     *
     * @param resName   a {@link java.lang.String} object.
     * @param resNum    a int.
     * @param create    a boolean.
     * @param defaultRT Default ResidueType if it cannot be assigned.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public Residue getResidue(String resName, int resNum, boolean create, Residue.ResidueType defaultRT) {
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements(); ) {
            Residue r = (Residue) e.nextElement();
            if (r.getResidueNumber() == resNum && r.getName().equalsIgnoreCase(resName)) {
                return r;
            }
        }
        if (!create) {
            return null;
        }
        Residue residue = null;
        resName = resName.toUpperCase();
        if (resName.length() == 1) {
            try {
                ResidueEnumerations.NucleicAcid1 na = ResidueEnumerations.NucleicAcid1.valueOf(resName);
                residue = new Residue(resName, resNum, Residue.ResidueType.NA,
                        chainID, getName());
            } catch (Exception e) {
                try {
                    ResidueEnumerations.AminoAcid1 aa = ResidueEnumerations.AminoAcid1.valueOf(resName);
                    residue = new Residue(resName, resNum, Residue.ResidueType.AA,
                            chainID, getName());
                } catch (Exception ex) {
                }
            }
        } else if (resName.length() >= 2) {
            try {
                ResidueEnumerations.NucleicAcid3 na = ResidueEnumerations.NucleicAcid3.valueOf(resName);
                residue = new Residue(resName, resNum, Residue.ResidueType.NA,
                        chainID, getName());
            } catch (Exception e) {
                try {
                    ResidueEnumerations.AminoAcid3 aa = ResidueEnumerations.AminoAcid3.valueOf(resName);
                    residue = new Residue(resName, resNum, Residue.ResidueType.AA,
                            chainID, getName());
                } catch (Exception ex) {
                }
            }
        }
        if (residue == null) {
            residue = new Residue(resName, resNum, defaultRT,
                    chainID, getName());
        }
        addMSNode(residue);
        return residue;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return hash(SEED, polymerNumber);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
                         Material mat) {
        // If coloring by Polymer, pass this Polymer's color
        if (newColorModel == RendererCache.ColorModel.POLYMER) {
            int index = polymerNumber % 10;
            color = polymerColor.get(index);
            mat = RendererCache.materialFactory(color);
        }
        for (ListIterator li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSGroup atomGroup = (MSGroup) li.next();
            atomGroup.setColor(newColorModel, color, mat);
        }
        for (Enumeration e = getTermNode().children(); e.hasMoreElements(); ) {
            Joint joint = (Joint) e.nextElement();
            joint.setColor(newColorModel);
        }
    }

    /**
     * <p>
     * Setter for the field <code>link</code>.</p>
     *
     * @param t a boolean.
     */
    public void setLink(boolean t) {
        link = t;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel newViewModel,
                        List<BranchGroup> newShapes) {
        for (ListIterator li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSGroup atomGroup = (MSGroup) li.next();
            atomGroup.setView(newViewModel, newShapes);
        }
        for (Enumeration e = getTermNode().children(); e.hasMoreElements(); ) {
            Joint joint = (Joint) e.nextElement();
            joint.setView(newViewModel, newShapes);
        }
    }
}
