/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.ListIterator;
import java.util.Vector;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

import ffx.potential.bonded.Residue.ResidueType;
import ffx.numerics.VectorMath;

/**
 * The Polymer class encapsulates a polypeptide or polynucleotide chain.
 */
public class Polymer extends MSGroup {

    private static final long serialVersionUID = 1L;
    public static final int MultiScaleLevel = 3;
    private static int count = 0;
    private static double[] da = new double[3];
    private static double[] db = new double[3];
    public static Hashtable<Integer, Color3f> polymerColor = new Hashtable<Integer, Color3f>();

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

    /**
     * Default Polymer Construtor
     */
    public Polymer() {
        super();
        polymerNumber = ++count;
    }

    public Polymer(String n) {
        super(n);
        polymerNumber = ++count;
    }

    /**
     * Polymer Constructor
     *
     * @param n
     *            Polymer indentifier, generally a letter in PDB files (ie
     *            A,B,C,etc)
     */
    public Polymer(String n, boolean l) {
        this(n);
        link = l;
    }

    /**
     * Polymer Constructor
     *
     * @param n
     *            Polymer indentifier, generally a letter in PDB files (ie
     *            A,B,C,etc)
     * @param residues
     *            Represents a FNode where the Polymer's residues have been
     *            attached
     */
    public Polymer(String n, MSNode residues) {
        super(n, residues);
        polymerNumber = ++count;
    }

    /**
     * A generic method for adding a MSNode to the Polymer.
     *
     * @param o
     *            If the MSNode is a Residue, it will be added to the Polymer,
     *            so long as its sequence number is not already in use.
     */
    @Override
    public MSNode addMSNode(MSNode o) {
        assert (o instanceof Residue);
        getAtomNode().add(o);

        return o;
    }

    /**
     * Joiner joins Moieties m1 and m2 and returns the Geometry objects formed
     * in a Joint.
     */
    public Joint createJoint(Residue residue1, Residue residue2) {
        Joint joint = null;
        for (Enumeration e = residue1.getAtomNode().children(); e.hasMoreElements();) {
            Atom a1 = (Atom) e.nextElement();
            a1.getXYZ(da);
            for (Enumeration e2 = residue2.getAtomNode().children(); e2.hasMoreElements();) {
                Atom a2 = (Atom) e2.nextElement();
                a2.getXYZ(db);
                double d1 = VectorMath.dist(da, db);
                double d2 = Bond.BUFF + a1.getVDWR() / 2 + a2.getVDWR() / 2;
                if (d1 < d2) {
                    Bond b = new Bond(a1, a2);
                    Joint newJoint = createJoint(b, residue1, residue2);
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
     * Overidden equals method.
     *
     * @param object
     *            Object to compare
     * @return True if object is not <b>this</b> Polymer, is of Class Polymer,
     *         and both object and this Polymer have identical names
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
     * Finalize should be called after all the Residues have been added to the
     * Polymer. This method in turn calls the Finalize method of each Residue,
     * then forms Joints between adjacent Residues in the Polymer
     */
    @Override
    public void finalize(boolean finalizeGroups) {
        ListIterator li;
        ArrayList res = getAtomNodeList();
        setFinalized(false);

        // Finalize the residues in the Polymer
        if (finalizeGroups) {
            for (li = res.listIterator(); li.hasNext();) {
                ((Residue) li.next()).finalize(true);
            }
        }
        // Join the residues in the Polymer
        if (link) {
            Residue residue = getFirstResidue();
            if (residue.residueType == ResidueType.AA) {
                getAtomNode().setName("Amino Acids " + "(" + res.size() + ")");
            } else if (residue.residueType == ResidueType.NA) {
                getAtomNode().setName("Nucleic Acids " + "(" + res.size() + ")");
            } else {
                getAtomNode().setName("Residues " + "(" + res.size() + ")");
            }
            Joint j;
            MSNode joints = getTerms();
            joints.removeAllChildren();
            List<Atom> atoms = getAtomList();

            for (Atom a : atoms) {
                if (a.getNumBonds() > 0) {
                    for (Bond b : a.getBonds()) {
                        if (!b.sameGroup() && b.getParent() == null) {
                            Residue r1 = (Residue) a.getMSNode(Residue.class);
                            Residue r2 = (Residue) b.get1_2(a).getMSNode(
                                    Residue.class);
                            j = createJoint(b, r1, r2);
                            joints.add(j);
                        }
                    }
                }
            }

            if (residue != null) {
                if (residue.residueType == ResidueType.AA) {
                    getTerms().setName(
                            "Peptide Bonds " + "(" + joints.getChildCount() + ")");
                } else {
                    getTerms().setName(
                            "Linkages " + "(" + joints.getChildCount() + ")");
                }
            } else {
                getTerms().setName(
                        "Linkages " + "(" + joints.getChildCount() + ")");
            }
        } else {
            getAtomNode().setName("Sub-Groups " + "(" + res.size() + ")");
            if (getTerms().getParent() != null) {
                remove(getTerms());
            }
        }
        removeLeaves();
        setFinalized(true);
    }

    public boolean getLink() {
        return link;
    }

    /**
     * Get the Phi Psi List for the Polymer
     *
     * @return An ArrayList of Dihedral objects representing the Phi/Psi angles
     *         of the Polymer, useful for creating Ramachandran plots
     */
    public Vector<ArrayList<Torsion>> getPhiPsiList() {
        MSNode dihedrals;
        ListIterator li, lj;
        Vector<ArrayList<Torsion>> phipsi = new Vector<ArrayList<Torsion>>();
        ArrayList<Torsion> phi = new ArrayList<Torsion>();
        ArrayList<Torsion> psi = new ArrayList<Torsion>();
        phipsi.add(phi);
        phipsi.add(psi);
        MSNode joints = getTerms();
        for (li = joints.getChildListIterator(); li.hasNext();) {
            dihedrals = ((Joint) li.next()).getTorsions();
            for (lj = dihedrals.getChildListIterator(); lj.hasNext();) {
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

    public Residue getFirstResidue() {
        MSNode atomNode = getAtomNode();
        if (atomNode == null) {
            return null;
        }
        return (Residue) atomNode.getChildAt(0);
    }

    public ArrayList<Residue> getResidues() {
        ArrayList<Residue> residues = new ArrayList<Residue>();
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements();) {
            Residue r = (Residue) e.nextElement();
            residues.add(r);
        }
        return residues;
    }

    public Residue getResidue(int resNum) {
        if (resNum > 0 && getAtomNode().getChildCount() >= resNum) {
            Residue r = (Residue) getAtomNode().getChildAt(resNum - 1);
            if (r.getResidueNumber() == resNum) {
                return r;
            }
        }
        // Fall back for non-ordered children
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements();) {
            Residue r = (Residue) e.nextElement();
            if (r.getResidueNumber() == resNum) {
                return r;
            }
        }
        return null;
    }

    public Residue getResidue(String resName, int resNum, boolean create) {
        if (resNum > 0 && getAtomNode().getChildCount() >= resNum) {
            Residue r = (Residue) getAtomNode().getChildAt(resNum - 1);
            if (r.getResidueNumber() == resNum && r.getName().equalsIgnoreCase(resName)) {
                return r;
            }
        }
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements();) {
            Residue r = (Residue) e.nextElement();
            if (r.getResidueNumber() == resNum && r.getName().equalsIgnoreCase(resName)) {
                return r;
            }
        }
        //System.out.println(resName + ": " + resNum);
        if (!create) {
            return null;
        }
        Residue residue = null;
        resName = resName.toUpperCase();
        if (resName.length() == 1) {
            try {
                Residue.NA1.valueOf(resName);
                residue = new Residue(resName, resNum, Residue.ResidueType.NA);
            } catch (Exception e) {
                try {
                    Residue.AA1.valueOf(resName);
                    residue = new Residue(resName, resNum, Residue.ResidueType.AA);
                } catch (Exception ex) {
                }
            }
        } else if (resName.length() == 2 || resName.length() == 3) {
            try {
                Residue.NA3.valueOf(resName);
                residue = new Residue(resName, resNum, Residue.ResidueType.NA);
            } catch (Exception e) {
                try {
                    Residue.AA3.valueOf(resName);
                    residue = new Residue(resName, resNum, Residue.ResidueType.AA);
                } catch (Exception ex) {
                }
            }
        }
        if (residue == null) {
            residue = new Residue(resName, resNum, Residue.ResidueType.UNK);
        }
        addMSNode(residue);
        return residue;
    }

    @Override
    public int hashCode() {
        return HashCodeUtil.hash(HashCodeUtil.POLYMERSEED, polymerNumber);
    }

    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
            Material mat) {
        // If coloring by Polymer, pass this Polymer's color
        if (newColorModel == RendererCache.ColorModel.POLYMER) {
            int index = polymerNumber % 10;
            color = polymerColor.get(index);
            mat = RendererCache.materialFactory(color);
        }
        for (ListIterator li = getAtomNodeList().listIterator(); li.hasNext();) {
            MSGroup atomGroup = (MSGroup) li.next();
            atomGroup.setColor(newColorModel, color, mat);
        }
        for (Enumeration e = getTerms().children(); e.hasMoreElements();) {
            Joint joint = (Joint) e.nextElement();
            joint.setColor(newColorModel);
        }
    }

    public void setLink(boolean t) {
        link = t;
    }

    @Override
    public void setView(RendererCache.ViewModel newViewModel,
            List<BranchGroup> newShapes) {
        for (ListIterator li = getAtomNodeList().listIterator(); li.hasNext();) {
            MSGroup atomGroup = (MSGroup) li.next();
            atomGroup.setView(newViewModel, newShapes);
        }
        for (Enumeration e = getTerms().children(); e.hasMoreElements();) {
            Joint joint = (Joint) e.nextElement();
            joint.setView(newViewModel, newShapes);
        }
    }
}
