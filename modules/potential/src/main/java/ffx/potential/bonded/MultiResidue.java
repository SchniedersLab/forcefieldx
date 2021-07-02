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

import static ffx.potential.bonded.AminoAcidUtils.assignAminoAcidAtomTypes;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.parameters.ForceField;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.logging.Logger;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Material;
import org.jogamp.vecmath.Color3f;

/**
 * The MultiResidue class allows switching between residues for uses such as sequence optimization.
 *
 * @author Will Tollefson
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MultiResidue extends Residue {

  private static final Logger logger = Logger.getLogger(MultiResidue.class.getName());

  /** The force field in use. */
  ForceField forceField;
  /** The active residue. */
  private Residue activeResidue;
  /** List of residues under consideration. */
  private final List<Residue> consideredResidues;
  /** Current rotamers. */
  private Rotamer[] rotamers;
  /** The original rotamer. */
  private Rotamer originalRotamer;

  /**
   * Constructor for MultiResidue.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public MultiResidue(Residue residue, ForceField forceField) {
    super(
        "" + residue.getResidueNumber() + "-" + "MultiResidue",
        residue.getResidueNumber(),
        residue.residueType,
        residue.getChainID(),
        residue.getChainID().toString());
    this.forceField = forceField;
    activeResidue = residue;
    // Initialize consideredResidue list.
    consideredResidues = new ArrayList<>();
    consideredResidues.add(residue);
    removeLeaves();
  }

  /** {@inheritDoc} */
  @Override
  public MSNode addMSNode(MSNode o) {
    if (o instanceof Residue) {
      add(o);
      return o;
    } else {
      return null;
    }
  }

  /**
   * addResidue.
   *
   * @param newResidue a {@link ffx.potential.bonded.Residue} object.
   */
  public void addResidue(Residue newResidue) {
    // Add the new residue to list.
    consideredResidues.add(newResidue);

    // Get references to nearby residues.
    Residue prevResidue = activeResidue.getPreviousResidue();
    Residue nextResidue = activeResidue.getNextResidue();
    Residue prev2Residue = null;
    if (prevResidue != null) {
      prev2Residue = prevResidue.getPreviousResidue();
    }
    Residue next2Residue = null;
    if (nextResidue != null) {
      next2Residue = nextResidue.getNextResidue();
    }

    moveBackBoneAtoms(activeResidue, newResidue);

    // Pass references of the active Residues' joints to the new Residue.
    List<Joint> joints = activeResidue.getJoints();
    for (Joint joint : joints) {
      newResidue.addJoint(joint);
    }

    // Make the new Residue active.
    activeResidue.removeFromParent();
    add(newResidue);
    activeResidue = newResidue;

    // Build side-chain atoms and assign atom types for the new Residue.
    try {
      assignAminoAcidAtomTypes(newResidue, prevResidue, nextResidue, forceField, null);
      if (nextResidue != null) {
        Atom C = (Atom) newResidue.getAtomNode("C");
        Atom nextN = (Atom) nextResidue.getAtomNode("N");
        for (Joint joint : joints) {
          Bond bond = joint.getBondList().get(0);
          if (bond.containsAtom(C) && bond.containsAtom(nextN)) {
            C.setBond(bond);
          }
        }
      }
    } catch (MissingHeavyAtomException | MissingAtomTypeException exception) {
      logger.severe(exception.toString());
    }
    newResidue.finalize(true, forceField);

    updateGeometry(newResidue, prevResidue, nextResidue, prev2Residue, next2Residue);
  }

  /** {@inheritDoc} */
  @Override
  public void assignBondedTerms(ForceField forceField) {
    activeResidue.assignBondedTerms(forceField);
  }

  /** {@inheritDoc} */
  @Override
  public void constructValenceTerms() {
    activeResidue.constructValenceTerms();
  }

  /** {@inheritDoc} */
  @Override
  public Joint createJoint(Bond bond, MSGroup group1, MSGroup group2, ForceField forceField) {
    return activeResidue.createJoint(bond, group1, group2, forceField);
  }

  /** {@inheritDoc} */
  @Override
  public Joint createJoint(MSGroup group1, MSGroup group2, ForceField forceField) {
    return activeResidue.createJoint(group1, group2, forceField);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Overidden equals method that return true if object is not equals to this, is of the same
   * class, has the same parent Polymer, the same sequence number, the same ResidueType, and the same
   * AA3/NA3.
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    MultiResidue multiResidue = (MultiResidue) o;
    return Objects.equals(getSegID(), multiResidue.getSegID())
        && Objects.equals(getResidueNumber(), multiResidue.getResidueNumber())
        && Objects.equals(getName(), multiResidue.getName());
  }

  /** {@inheritDoc} */
  @Override
  public void finalize(boolean finalizeGeometry, ForceField forceField) {
    activeResidue.finalize(finalizeGeometry, forceField);
  }

  /** {@inheritDoc} */
  @Override
  public void findDangelingAtoms() {
    activeResidue.findDangelingAtoms();
  }

  /**
   * getActive.
   *
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public Residue getActive() {
    return activeResidue;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns the AminoAcid3 of the active residue.
   */
  @Override
  public AminoAcid3 getAminoAcid3() {
    return activeResidue.getAminoAcid3();
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getAngles() {
    return activeResidue.getAngles();
  }

  /** {@inheritDoc} */
  @Override
  public void setAngles(MSNode t) {
    activeResidue.setAngles(t);
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getAtomNode() {
    return activeResidue.getAtomNode();
  }

  /** {@inheritDoc} */
  @Override
  public void setAtomNode(MSNode t) {
    activeResidue.setAtomNode(t);
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getAtomNode(int index) {
    return activeResidue.getAtomNode(index);
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getAtomNode(String n) {
    return activeResidue.getAtomNode(n);
  }

  /** {@inheritDoc} */
  @Override
  public List<MSNode> getAtomNodeList() {
    return activeResidue.getAtomNodeList();
  }

  /** {@inheritDoc} */
  @Override
  public Bond getBond(String id) {
    return activeResidue.getBond(id);
  }

  /** {@inheritDoc} */
  @Override
  public Bond getBond(int index) {
    return activeResidue.getBond(index);
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getBonds() {
    return activeResidue.getBonds();
  }

  /** {@inheritDoc} */
  @Override
  public void setBonds(MSNode t) {
    activeResidue.setBonds(t);
  }

  /** {@inheritDoc} */
  @Override
  public double[] getCenter() {
    return activeResidue.getCenter();
  }

  /** {@inheritDoc} */
  @Override
  public void setCenter(double[] d) {
    activeResidue.setCenter(d);
  }

  /**
   * Returns a copy of this MultiResidue's consideredResidues array.
   *
   * @return a new List of the considered residues.
   */
  public List<Residue> getConsideredResidues() {
    return new ArrayList<>(consideredResidues);
  }

  /** {@inheritDoc} */
  @Override
  public List<Atom> getDangelingAtoms() {
    return activeResidue.getDangelingAtoms();
  }

  /** {@inheritDoc} */
  @Override
  public void setDangelingAtoms(List<Atom> a) {
    activeResidue.setDangelingAtoms(a);
  }

  /**
   * Returns a list of this MultiResidue's inactive residues. Adding/removing from the returned list
   * does nothing.
   *
   * @return a new List of inactive residues.
   */
  public List<Residue> getInactive() {
    List<Residue> ret = new ArrayList<>();
    for (Residue res : consideredResidues) {
      if (res != activeResidue) {
        ret.add(res);
      }
    }
    return ret;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getMultiScaleCenter(boolean w) {
    return activeResidue.getMultiScaleCenter(w);
  }

  /** {@inheritDoc} */
  @Override
  public String getName() {
    if (activeResidue != null) {
      return activeResidue.getName();
    }
    return super.getName();
  }

  /**
   * getResidueCount.
   *
   * @return a int.
   */
  public int getResidueCount() {
    if (consideredResidues == null) {
      return 0;
    }
    return consideredResidues.size();
  }

  /** {@inheritDoc} */
  @Override
  public Rotamer[] getRotamers(RotamerLibrary library) {
    if (rotamers != null) {
      return rotamers;
    }

    List<Rotamer[]> usual = new ArrayList<>();
    int nRots = 0;

    for (Residue residue : consideredResidues) {
      Rotamer[] rotamers = library.getRotamers(residue);
      if (rotamers != null && rotamers.length > 0) {
        usual.add(rotamers);
        nRots += rotamers.length;
      }
    }

    if (library.getUsingOrigCoordsRotamer()) {
      if (originalRotamer == null
          && (residueType == ResidueType.AA || residueType == ResidueType.NA)) {
        ResidueState origState = storeState();
        double[] chi = RotamerLibrary.measureRotamer(activeResidue, false);
        if (residueType == ResidueType.AA) {
          AminoAcid3 aa3 = this.getAminoAcid3();
          originalRotamer = new Rotamer(aa3, origState, chi);
        } else if (residueType == ResidueType.NA) {
          NucleicAcid3 na3 = this.getNucleicAcid3();
          originalRotamer = new Rotamer(na3, origState, chi);
        }
      }
      Rotamer[] allRotamers;
      if (originalRotamer != null) {
        allRotamers = new Rotamer[nRots + 1];
        int index;

        if (origAtEnd) {
          index = 0;
          allRotamers[allRotamers.length - 1] = originalRotamer;
        } else {
          index = 1;
          allRotamers[0] = originalRotamer;
        }

        for (Rotamer[] rotamersI : usual) {
          int nrotamers = rotamersI.length;
          arraycopy(rotamersI, 0, allRotamers, index, nrotamers);
          index += nrotamers;
        }
      } else {
        allRotamers = addAllDefaultRotamers(usual, nRots);
      }
      rotamers = allRotamers;
      return rotamers;
    } else {
      rotamers = addAllDefaultRotamers(usual, nRots);
      return rotamers;
    }
  }

  /** {@inheritDoc} */
  @Override
  public List<Atom> getSideChainAtoms() {
    return activeResidue.getSideChainAtoms();
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getTermNode() {
    return activeResidue.getTermNode();
  }

  /** {@inheritDoc} */
  @Override
  public MSNode getTorsions() {
    return activeResidue.getTorsions();
  }

  /** {@inheritDoc} */
  @Override
  public void setTorsions(MSNode t) {
    activeResidue.setTorsions(t);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns all atoms (all atoms are variable during DEE).
   */
  @Override
  public List<Atom> getVariableAtoms() {
    return activeResidue.getAtomList();
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(getSegID(), getResidueNumber(), getName());
  }

  /** {@inheritDoc} */
  @Override
  public boolean isFinalized() {
    return activeResidue.isFinalized();
  }

  /** {@inheritDoc} */
  @Override
  public void setFinalized(boolean t) {
    activeResidue.setFinalized(t);
  }

  /** {@inheritDoc} */
  @Override
  public void reOrderAtoms() {
    activeResidue.reOrderAtoms();
  }

  /**
   * Request the passed amino acid be set as the active residue.
   *
   * @param aa a {@link ffx.potential.bonded.ResidueEnumerations.AminoAcid3} object.
   * @return true if the request is satisfied.
   */
  public boolean requestSetActiveResidue(AminoAcid3 aa) {
    if (aa == AminoAcid3.valueOf(activeResidue.getName())) {
      return true;
    }
    for (Residue res : consideredResidues) {
      if (aa == AminoAcid3.valueOf(res.getName())) {
        return setActiveResidue(res);
      }
    }
    logger.warning(
        format(" Couldn't assign residue %s to MultiResidue %s.", aa, this));
    return false;
  }

  /**
   * Request the ith residue be set active.
   *
   * @param i a int.
   * @return true if the ith residue was set active, false otherwise.
   */
  public boolean setActiveResidue(int i) {
    if (consideredResidues == null) {
      return false;
    }
    if (i >= consideredResidues.size()) {
      return false;
    }
    return setActiveResidue(consideredResidues.get(i));
  }

  /**
   * Setter for the field <code>activeResidue</code>.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @return a boolean.
   */
  public boolean setActiveResidue(Residue residue) {
    if (!consideredResidues.contains(residue)) {
      return false;
    }
    if (residue == activeResidue) {
      return true;
    }
    Residue prevResidue = activeResidue.getPreviousResidue();
    Residue nextResidue = activeResidue.getNextResidue();
    Residue prev2Residue = null;
    if (prevResidue != null) {
      prev2Residue = prevResidue.getPreviousResidue();
    }
    Residue next2Residue = null;
    if (nextResidue != null) {
      next2Residue = nextResidue.getNextResidue();
    }

    activeResidue.removeFromParent();

    // Move backbone atoms to the new active residue.
    moveBackBoneAtoms(activeResidue, residue);
    updateGeometry(residue, prevResidue, nextResidue, prev2Residue, next2Residue);
    activeResidue = residue;

    setName(toString());
    add(activeResidue);

    return true;
  }

  /**
   * Method may be redundant with requestSetActiveResidue. Will not function correctly if there is
   * more than one residue of type UNK (unknown).
   *
   * @param aa a {@link ffx.potential.bonded.ResidueEnumerations.AminoAcid3} object.
   * @return True if successful
   */
  public boolean setActiveResidue(AminoAcid3 aa) {
    Residue residue = null;
    for (Residue res : consideredResidues) {
      if (res.getAminoAcid3() == aa) {
        residue = res;
        break;
      }
    }
    if (residue == null) {
      return false;
    }
    return setActiveResidue(residue);
  }

  /** {@inheritDoc} */
  @Override
  public void setColor(RendererCache.ColorModel newColorModel, Color3f color, Material mat) {
    activeResidue.setColor(newColorModel, color, mat);
  }

  /** {@inheritDoc} */
  @Override
  public void setOutOfPlaneBends(MSNode t) {
    activeResidue.setOutOfPlaneBends(t);
  }

  /** {@inheritDoc} */
  @Override
  public void setPiOrbitalTorsions(MSNode t) {
    activeResidue.setPiOrbitalTorsions(t);
  }

  /** {@inheritDoc} */
  @Override
  public void setStretchBends(MSNode t) {
    activeResidue.setStretchBends(t);
  }

  /** {@inheritDoc} */
  @Override
  public void setTerms(MSNode t) {
    activeResidue.setTerms(t);
  }

  /** {@inheritDoc} */
  @Override
  public void setTorsionTorsions(MSNode t) {
    activeResidue.setTorsionTorsions(t);
  }

  /** {@inheritDoc} */
  @Override
  public void setUreyBradleys(MSNode t) {
    activeResidue.setUreyBradleys(t);
  }

  /** {@inheritDoc} */
  @Override
  public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
    activeResidue.setView(newViewModel, newShapes);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Publicly accessible method for storing a MultiResidue state.
   */
  @Override
  public ResidueState storeState() {
    return storeMultiResState();
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    int resNum = consideredResidues.get(0).getResidueNumber();
    StringBuilder sb = new StringBuilder();
    sb.append(resNum).append("-");
    for (Residue res : consideredResidues) {
      int num = ResidueEnumerations.getAminoAcidNumber(res.getName());
      String aa1 = ResidueEnumerations.AminoAcid1.values()[num].toString();
      if (res == activeResidue) {
        sb.append("[").append(aa1).append("]");
      } else {
        sb.append(aa1);
      }
    }
    return sb.toString();
  }

  /** {@inheritDoc} */
  @Override
  public void update() {
    activeResidue.update();
  }

  /** {@inheritDoc} */
  @Override
  public void updateAtoms() {
    activeResidue.updateAtoms();
  }

  /** {@inheritDoc} */
  @Override
  public void updateBonds() {
    activeResidue.updateBonds();
  }

  /**
   * Non-overrideable implementation method for storeState. Probably unnecessary, as I dropped the
   * idea of initializing the original-coordinates rotamer in the constructor.
   *
   * @return A ResidueState.
   */
  private ResidueState storeMultiResState() {
    return new ResidueState(this, activeResidue);
  }

  /**
   * Returns an array of all standard torsion-based Rotamers for this Multi-Residue.
   *
   * @param rotamerList List of Rotamer arrays to flatten
   * @param nRots Number of rotamers.
   * @return An array of all standard torsion-based Rotamers for this Multi-Residue.
   */
  private Rotamer[] addAllDefaultRotamers(List<Rotamer[]> rotamerList, int nRots) {
    Rotamer[] allRotamers = new Rotamer[nRots];
    int index = 0;
    for (Rotamer[] rotamers : rotamerList) {
      int nrotamers = rotamers.length;
      arraycopy(rotamers, 0, allRotamers, index, nrotamers);
      index += nrotamers;
    }
    return allRotamers;
  }

  @Override
  public void revertState(ResidueState state) {
    Residue res = state.getStateResidue();
    if (!setActiveResidue(res)) {
      throw new IllegalArgumentException(
          format(
              " Could not revert " + "multi-residue %s to residue identity %s",
              this, state.getStateResidue().toString()));
    }
    for (Atom atom : getAtomList()) {
      atom.moveTo(state.getAtomCoords(atom));
    }
  }

  private void moveBackBoneAtoms(Residue fromResidue, Residue toResidue) {
    Residue prevRes = this.getPreviousResidue();
    Residue nextRes = this.getNextResidue();

    // Begin with atoms common to all residues

    // Get references to the backbone atoms.
    Atom CA = (Atom) fromResidue.getAtomNode("CA");
    Atom C = (Atom) fromResidue.getAtomNode("C");
    Atom HA = (Atom) fromResidue.getAtomNode("HA");
    Atom N = (Atom) fromResidue.getAtomNode("N");
    Atom O = (Atom) fromResidue.getAtomNode("O");

    CA.removeFromParent();
    HA.removeFromParent();
    C.removeFromParent();
    O.removeFromParent();
    N.removeFromParent();

    // Clear their references to bonded geometry.
    CA.clearGeometry();
    HA.clearGeometry();
    C.clearGeometry();
    O.clearGeometry();
    N.clearGeometry();

    // Change their residue name.
    String resName = toResidue.getName();
    CA.setResName(resName);
    HA.setResName(resName);
    C.setResName(resName);
    O.setResName(resName);
    N.setResName(resName);

    // Add the backbone atoms to the new Residue.
    toResidue.addMSNode(CA);
    toResidue.addMSNode(HA);
    toResidue.addMSNode(C);
    toResidue.addMSNode(O);
    toResidue.addMSNode(N);

    if (prevRes == null) {
      Atom H1 = (Atom) fromResidue.getAtomNode("H1");
      Atom H2 = (Atom) fromResidue.getAtomNode("H2");
      Atom H3 = (Atom) fromResidue.getAtomNode("H3");

      H1.removeFromParent();
      H2.removeFromParent();

      H1.clearGeometry();
      H2.clearGeometry();
      H1.setResName(resName);
      H2.setResName(resName);
      toResidue.addMSNode(H1);
      toResidue.addMSNode(H2);

      if (H3 != null) {
        H3.removeFromParent();
        H3.clearGeometry();
        H3.setResName(resName);
        toResidue.addMSNode(H3);
      }
    } else {
      Atom H = (Atom) fromResidue.getAtomNode("H");
      H.removeFromParent();
      H.clearGeometry();
      H.setResName(resName);
      toResidue.addMSNode(H);
    }

    if (nextRes == null) {
      Atom OXT = (Atom) fromResidue.getAtomNode("OXT");
      if (OXT != null) {
        OXT.removeFromParent();
        OXT.clearGeometry();
        OXT.setResName(resName);
        toResidue.addMSNode(OXT);
      } else {
        Atom OH = (Atom) fromResidue.getAtomNode("OH");
        Atom HO = (Atom) fromResidue.getAtomNode("HO");
        OH.removeFromParent();
        HO.removeFromParent();
        OH.clearGeometry();
        HO.clearGeometry();
        OH.setResName(resName);
        HO.setResName(resName);
        toResidue.addMSNode(OH);
        toResidue.addMSNode(HO);
      }
    }
  }

  /**
   * Update Atom references to local geometry.
   *
   * @param residue Current residue.
   * @param prev Previous residue.
   * @param next Next residue.
   * @param prev2 Previous previous residue.
   * @param next2 Next next residue.
   */
  private void updateGeometry(
      Residue residue, Residue prev, Residue next, Residue prev2, Residue next2) {
    if (residue == null) {
      return;
    }

    // Update atom references to local geometry.
    List<Atom> atoms = residue.getAtomList();
    List<Bond> bonds = residue.getBondList();
    List<Angle> angles = residue.getAngleList();
    List<Torsion> torsions = residue.getTorsionList();
    if (prev != null) {
      atoms.addAll(prev.getAtomList());
      bonds.addAll(prev.getBondList());
      angles.addAll(prev.getAngleList());
      torsions.addAll(prev.getTorsionList());
      List<Joint> joints = prev.getJoints();
      for (Joint joint : joints) {
        bonds.addAll(joint.getBondList());
        angles.addAll(joint.getAngleList());
        torsions.addAll(joint.getTorsionList());
      }
    }
    if (prev2 != null) {
      bonds.addAll(prev2.getBondList());
      angles.addAll(prev2.getAngleList());
      torsions.addAll(prev2.getTorsionList());
    }
    if (next != null) {
      atoms.addAll(next.getAtomList());
      bonds.addAll(next.getBondList());
      angles.addAll(next.getAngleList());
      torsions.addAll(next.getTorsionList());
      List<Joint> joints = next.getJoints();
      for (Joint joint : joints) {
        bonds.addAll(joint.getBondList());
        angles.addAll(joint.getAngleList());
        torsions.addAll(joint.getTorsionList());
      }
    }
    if (next2 != null) {
      bonds.addAll(next2.getBondList());
      angles.addAll(next2.getAngleList());
      torsions.addAll(next2.getTorsionList());
    }

    for (Atom atom : atoms) {
      atom.clearGeometry();
    }
    for (Atom atom : atoms) {
      for (Bond b : bonds) {
        if (b.containsAtom(atom)) {
          atom.setBond(b);
        }
      }
    }
    for (Atom atom : atoms) {
      for (Angle a : angles) {
        if (a.containsAtom(atom)) {
          atom.setAngle(a);
        }
      }
    }
    for (Atom atom : atoms) {
      for (Torsion t : torsions) {
        if (t.containsAtom(atom)) {
          atom.setTorsion(t);
        }
      }
    }
  }

}
