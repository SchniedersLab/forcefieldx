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
package ffx.potential.bonded;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.parameters.ForceField;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Logger;
import javax.swing.tree.TreeNode;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Material;
import org.jogamp.vecmath.Color3f;

/**
 * The Polymer class encapsulates a peptide or nucleotide chain.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Polymer extends MSGroup {

  private static final Logger logger = Logger.getLogger(Polymer.class.getName());

  /** Constant <code>polymerColor</code> */
  private static final Map<Integer, Color3f> polymerColor = new HashMap<>();
  /** Polymer count. */
  private static int count = 0;

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

  /** Flag to indicate the residues in the polymer should be joined. */
  private boolean link = false;
  /** The number of this Polymer. */
  private final int polymerNumber;
  /** The ChainID of this Polymer. */
  private final Character chainID;

  /**
   * Polymer constructor.
   *
   * @param chainID Possibly redundant PDB chainID.
   * @param segID Unique identifier from A-Z,0-9, then 1A-1Z,10-19, etc.
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
   * @param segID Unique identifier from A-Z,0-9, then 1A-1Z,10-19, etc.
   * @param link a boolean.
   */
  public Polymer(Character chainID, String segID, boolean link) {
    this(chainID, segID);
    this.link = link;
  }

  /**
   * Polymer Constructor.
   *
   * @param segID A unique identifier from A-Z,0-9, then 1A-1Z,10-19, etc.
   * @param residues Represents a MSNode where the Polymer's residues have been attached.
   * @param chainID a {@link java.lang.Character} object.
   */
  public Polymer(Character chainID, String segID, MSNode residues) {
    super(segID, residues);
    this.chainID = chainID;
    polymerNumber = ++count;
  }

  /**
   * {@inheritDoc}
   *
   * <p>A generic method for adding a MSNode to the Polymer.
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

    residueNode.insert(residue, childIndex);

    residue.setChainID(chainID);

    return msNode;
  }

  /**
   * addMultiResidue.
   *
   * @param multiResidue a {@link ffx.potential.bonded.MultiResidue} object.
   */
  public void addMultiResidue(MultiResidue multiResidue) {
    Residue residue = multiResidue.getActive();
    MSNode residueNode = getAtomNode();
    int index = residueNode.getIndex(residue);

    if (index < 0) {
      System.err.println(
          "WARNING!  Polymer::addMultiResidue did not find a corresponding Residue on Polymer.");
      residueNode.add(multiResidue);
    } else {
      residue.removeFromParent();
      residueNode.insert(multiResidue, index);
      multiResidue.add(residue);
    }
  }

  /**
   * addMultiTerminus.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param multiTerminus a {@link ffx.potential.bonded.MultiTerminus} object.
   */
  public void addMultiTerminus(Residue residue, MultiTerminus multiTerminus) {
    List<MSNode> children = residue.getChildList();
    for (MSNode child : children) {
      multiTerminus.add(child);
    }
    MSNode residueNode = getAtomNode();
    int index = residueNode.getIndex(residue);
    residueNode.remove(index);
    residueNode.insert(multiTerminus, index);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Joiner joins Moieties m1 and m2 and returns the Geometry objects formed in a Joint.
   */
  public Joint createJoint(Residue residue1, Residue residue2, ForceField forceField) {
    Joint joint = null;
    double[] da = new double[3];
    double[] db = new double[3];
    for (Enumeration<TreeNode> e = residue1.getAtomNode().children(); e.hasMoreElements(); ) {
      Atom a1 = (Atom) e.nextElement();
      a1.getXYZ(da);
      for (Enumeration<TreeNode> e2 = residue2.getAtomNode().children(); e2.hasMoreElements(); ) {
        Atom a2 = (Atom) e2.nextElement();
        a2.getXYZ(db);
        double d1 = DoubleMath.dist(da, db);
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
   *
   * <p>Overidden equals method.
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    Polymer polymer = (Polymer) o;
    return polymerNumber == polymer.polymerNumber && Objects.equals(getName(), polymer.getName());
  }

  /**
   * {@inheritDoc}
   *
   * <p>Finalize should be called after all the Residues have been added to the Polymer. This method
   * in turn calls the Finalize method of each Residue, then forms Joints between adjacent Residues
   * in the Polymer
   */
  @Override
  public void finalize(boolean finalizeGroups, ForceField forceField) {
    List<MSNode> residues = getAtomNodeList();
    setFinalized(false);

    // Finalize the residues in the Polymer
    if (finalizeGroups) {
      for (MSNode node : residues) {
        Residue residue = (Residue) node;
        residue.finalize(true, forceField);
      }
    }
    // Join the residues in the Polymer
    if (link) {
      Residue residue = getFirstResidue();
      if (residue.residueType == Residue.ResidueType.AA) {
        getAtomNode().setName("Amino Acids " + "(" + residues.size() + ")");
      } else if (residue.residueType == Residue.ResidueType.NA) {
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
              Residue r1 = a.getMSNode(Residue.class);
              Residue r2 = b.get1_2(a).getMSNode(Residue.class);
              j = createJoint(b, r1, r2, forceField);
              joints.add(j);
            }
          }
        }
      }

      if (residue.residueType == Residue.ResidueType.AA) {
        getTermNode().setName("Peptide Bonds " + "(" + joints.getChildCount() + ")");
      } else {
        getTermNode().setName("Linkages " + "(" + joints.getChildCount() + ")");
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
   * Getter for the field <code>chainID</code>.
   *
   * @return a {@link java.lang.Character} object.
   */
  public Character getChainID() {
    return chainID;
  }

  /**
   * getFirstResidue
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
   * Getter for the field <code>link</code>.
   *
   * @return a boolean.
   */
  public boolean getLink() {
    return link;
  }

  /**
   * Setter for the field <code>link</code>.
   *
   * @param link a boolean.
   */
  public void setLink(boolean link) {
    this.link = link;
  }

  /**
   * TODO: Was the sole hook on BondedTerm equality definition via getID(); will rewrite with a
   * simple Comparator soon.
   *
   * @return An List of Dihedral objects representing the Phi/Psi angles of the Polymer, useful for
   *     creating Ramachandran plots
   */
  public List<List<Torsion>> getPhiPsiList() {
    List<List<Torsion>> phipsi = new ArrayList<>();
    List<Torsion> phi = new ArrayList<>();
    List<Torsion> psi = new ArrayList<>();
    phipsi.add(phi);
    phipsi.add(psi);
    for (Residue residue : this.getResidues()) {
      for (Torsion torsion : residue.getTorsionList()) {
        Atom[] atoms = torsion.atoms;
        StringBuilder s = new StringBuilder(atoms[0].getName());
        for (int i=1; i<4; i++) {
          s.append("-").append(atoms[i].getName());
        }
        // Phi
        if (s.toString().equals("C-N-CA-C") || s.toString().equals("C-CA-N-C")) {
          phi.add(torsion);
        } // Psi
        else if (s.toString().equals("N-C-CA-N") || s.toString().equals("N-CA-C-N")) {
          psi.add(torsion);
        }
      }
    }
    return phipsi;
  }

  /**
   * getResidue
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
    for (Enumeration<TreeNode> e = getAtomNode().children(); e.hasMoreElements(); ) {
      Residue r = (Residue) e.nextElement();
      if (r.getResidueNumber() == resNum) {
        return r;
      }
    }
    return null;
  }

  /**
   * getResidue
   *
   * @param resName a {@link java.lang.String} object.
   * @param resNum a int.
   * @param create a boolean.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public Residue getResidue(String resName, int resNum, boolean create) {
    return getResidue(resName, resNum, create, Residue.ResidueType.UNK);
  }

  /**
   * getResidue
   *
   * @param resName a {@link java.lang.String} object.
   * @param resNum a int.
   * @param create a boolean.
   * @param defaultRT Default ResidueType if it cannot be assigned.
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public Residue getResidue(
      String resName, int resNum, boolean create, ResidueType defaultRT) {
    for (Enumeration<TreeNode> e = getAtomNode().children(); e.hasMoreElements(); ) {
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
        NucleicAcidUtils.NucleicAcid1.valueOf(resName);
        residue = new Residue(resName, resNum, Residue.ResidueType.NA, chainID, getName());
      } catch (Exception e) {
        try {
          AminoAcidUtils.AminoAcid1.valueOf(resName);
          residue = new Residue(resName, resNum, Residue.ResidueType.AA, chainID, getName());
        } catch (Exception ex) {
          //
        }
      }
    } else if (resName.length() >= 2) {
      try {
        NucleicAcidUtils.NucleicAcid3.valueOf(resName);
        residue = new Residue(resName, resNum, Residue.ResidueType.NA, chainID, getName());
      } catch (Exception e) {
        try {
          AminoAcidUtils.AminoAcid3.valueOf(resName);
          residue = new Residue(resName, resNum, Residue.ResidueType.AA, chainID, getName());
        } catch (Exception ex) {
          //
        }
      }
    }
    if (residue == null) {
      residue = new Residue(resName, resNum, defaultRT, chainID, getName());
    }
    addMSNode(residue);
    return residue;
  }

  /**
   * getResidues
   *
   * @return a {@link java.util.List} object.
   */
  public List<Residue> getResidues() {
    List<Residue> residues = new ArrayList<>();
    for (Enumeration<TreeNode> e = getAtomNode().children(); e.hasMoreElements(); ) {
      Residue r = (Residue) e.nextElement();
      residues.add(r);
    }
    return residues;
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(polymerNumber, getName());
  }

  /** {@inheritDoc} */
  @Override
  public void setColor(RendererCache.ColorModel newColorModel, Color3f color, Material mat) {
    // If coloring by Polymer, pass this Polymer's color
    if (newColorModel == RendererCache.ColorModel.POLYMER) {
      int index = polymerNumber % 10;
      color = polymerColor.get(index);
      mat = RendererCache.materialFactory(color);
    }
    for (MSNode node : getAtomNodeList()) {
      MSGroup atomGroup = (MSGroup) node;
      atomGroup.setColor(newColorModel, color, mat);
    }
    for (Enumeration<TreeNode> e = getTermNode().children(); e.hasMoreElements(); ) {
      Joint joint = (Joint) e.nextElement();
      joint.setColor(newColorModel);
    }
  }

  /** {@inheritDoc} */
  @Override
  public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
    for (MSNode node : getAtomNodeList()) {
      MSGroup atomGroup = (MSGroup) node;
      atomGroup.setView(newViewModel, newShapes);
    }
    for (Enumeration<TreeNode> e = getTermNode().children(); e.hasMoreElements(); ) {
      Joint joint = (Joint) e.nextElement();
      joint.setView(newViewModel, newShapes);
    }
  }
}
