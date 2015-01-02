/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

import ffx.potential.bonded.BondedUtils.MissingAtomTypeException;
import ffx.potential.bonded.BondedUtils.MissingHeavyAtomException;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parameters.ForceField;

import static ffx.potential.bonded.AminoAcidUtils.assignAminoAcidAtomTypes;
import static ffx.potential.bonded.BondedUtils.buildBond;

/**
 * @author Will Tollefson and Michael J. Schnieders
 */
public class MultiResidue extends Residue {

    private static final Logger logger = Logger.getLogger(MultiResidue.class.getName());

    /**
     * The active residue.
     */
    private Residue activeResidue = null;

    /**
     * List of residues under consideration.
     */
    ArrayList<Residue> consideredResidues;

    /**
     * Force field in use.
     */
    ForceField forceField;

    public MultiResidue(Residue residue, ForceField forceField) {
        super("MultiResidue", residue.getResidueNumber(), residue.residueType);
        this.forceField = forceField;
        activeResidue = residue;
        // Initialize consideredResidue list.
        consideredResidues = new ArrayList<>();
        consideredResidues.add(residue);
        removeLeaves();
    }

    @Override
    public MSNode addMSNode(MSNode o) {
        if (o instanceof Residue) {
            add(o);
            return o;
        } else {
            return null;
        }
    }

    @Override
    public void assignBondedTerms(ForceField forceField) {
        activeResidue.assignBondedTerms(forceField);
    }

    @Override
    public void constructValenceTerms() {
        activeResidue.constructValenceTerms();
    }

    @Override
    public Joint createJoint(Bond bond, MSGroup group1, MSGroup group2, ForceField forceField) {
        return activeResidue.createJoint(bond, group1, group2, forceField);
    }

    @Override
    public Joint createJoint(MSGroup group1, MSGroup group2, ForceField forceField) {
        return activeResidue.createJoint(group1, group2, forceField);
    }

    @Override
    public void finalize(boolean finalizeGeometry, ForceField forceField) {
        activeResidue.finalize(finalizeGeometry, forceField);
    }

    @Override
    public void findDangelingAtoms() {
        activeResidue.findDangelingAtoms();
    }

    public Residue getActive() {
        return activeResidue;
    }

    @Override
    public MSNode getAngles() {
        return activeResidue.getAngles();
    }

    @Override
    public MSNode getAtomNode() {
        return activeResidue.getAtomNode();
    }

    @Override
    public MSNode getAtomNode(int index) {
        return activeResidue.getAtomNode(index);
    }

    @Override
    public MSNode getAtomNode(String n) {
        return activeResidue.getAtomNode(n);
    }

    @Override
    public ArrayList<MSNode> getAtomNodeList() {
        return activeResidue.getAtomNodeList();
    }

    @Override
    public Bond getBond(String id) {
        return activeResidue.getBond(id);
    }

    @Override
    public Bond getBond(int index) {
        return activeResidue.getBond(index);
    }

    @Override
    public MSNode getBonds() {
        return activeResidue.getBonds();
    }

    @Override
    public double[] getCenter() {
        return activeResidue.getCenter();
    }

    @Override
    public ArrayList getDangelingAtoms() {
        return activeResidue.getDangelingAtoms();
    }

    @Override
    public double[] getMultiScaleCenter(boolean w) {
        return activeResidue.getMultiScaleCenter(w);
    }

    @Override
    public MSNode getTerms() {
        return activeResidue.getTerms();
    }

    @Override
    public MSNode getTorsions() {
        return activeResidue.getTorsions();
    }

    @Override
    public boolean isFinalized() {
        return activeResidue.isFinalized();
    }

    @Override
    public void reOrderAtoms() {
        activeResidue.reOrderAtoms();
    }

    @Override
    public void setAngles(MSNode t) {
        activeResidue.setAngles(t);
    }

    @Override
    public void setAtomNode(MSNode t) {
        activeResidue.setAtomNode(t);
    }

    @Override
    public void setBonds(MSNode t) {
        activeResidue.setBonds(t);
    }

    @Override
    public void setCenter(double[] d) {
        activeResidue.setCenter(d);
    }

    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
            Material mat) {
        activeResidue.setColor(newColorModel, color, mat);
    }

    @Override
    public void setDangelingAtoms(ArrayList<Atom> a) {
        activeResidue.setDangelingAtoms(a);
    }

    @Override
    public void setFinalized(boolean t) {
        activeResidue.setFinalized(t);
    }

    @Override
    public void setOutOfPlaneBends(MSNode t) {
        activeResidue.setOutOfPlaneBends(t);
    }

    @Override
    public void setPiOrbitalTorsions(MSNode t) {
        activeResidue.setPiOrbitalTorsions(t);
    }

    @Override
    public void setStretchBends(MSNode t) {
        activeResidue.setStretchBends(t);
    }

    @Override
    public void setTerms(MSNode t) {
        activeResidue.setTerms(t);
    }

    @Override
    public void setTorsions(MSNode t) {
        activeResidue.setTorsions(t);
    }

    @Override
    public void setTorsionTorsions(MSNode t) {
        activeResidue.setTorsionTorsions(t);
    }

    @Override
    public void setUreyBradleys(MSNode t) {
        activeResidue.setUreyBradleys(t);
    }

    @Override
    public void setView(RendererCache.ViewModel newViewModel,
            List<BranchGroup> newShapes) {
        activeResidue.setView(newViewModel, newShapes);
    }

    @Override
    public void update() {
        activeResidue.update();
    }

    @Override
    public void updateAtoms() {
        activeResidue.updateAtoms();
    }

    @Override
    public void updateBonds() {
        activeResidue.updateBonds();
    }

    @Override
    public Rotamer[] getRotamers(Residue residue) {
        if (residue == null) {
            return null;
        }
        Rotamer allRotamers[];
        Residue residueOptions[] = consideredResidues.toArray(new Residue[consideredResidues.size()]);
        int nResidues = residueOptions.length;
        int rotamerTotal = 0;
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residueOptions[i];
            Rotamer rotamersi[] = RotamerLibrary.getRotamers(residuei);
            if (rotamersi == null) {
                continue;
            }
            rotamerTotal += rotamersi.length;
        }
        allRotamers = new Rotamer[rotamerTotal];
        for (int i = 0; i < nResidues; i++) {
            Residue residuei = residueOptions[i];
            Rotamer rotamersi[] = RotamerLibrary.getRotamers(residuei);
            if (rotamersi == null) {
                continue;
            }
            int shift = 0;
            for (int j = 0; j < rotamersi.length; j++) {
                allRotamers[j + shift] = rotamersi[j];
            }
            shift += rotamersi.length;
        }
        logger.info(consideredResidues.size() + " residue options with " + rotamerTotal + " rotamers.");
        return allRotamers;
    }

    private void moveBackBoneAtoms(Residue fromResidue, Residue toResidue) {
        /**
         * Get references to the backbone atoms.
         */
        Atom CA = (Atom) fromResidue.getAtomNode("CA");
        Atom HA = (Atom) fromResidue.getAtomNode("HA");
        Atom C = (Atom) fromResidue.getAtomNode("C");
        Atom O = (Atom) fromResidue.getAtomNode("O");
        Atom N = (Atom) fromResidue.getAtomNode("N");
        Atom H = (Atom) fromResidue.getAtomNode("H");
        /**
         * Detach them from their parent Residue.
         */
        CA.removeFromParent();
        HA.removeFromParent();
        C.removeFromParent();
        O.removeFromParent();
        N.removeFromParent();
        H.removeFromParent();
        /**
         * Clear their references to bonded geometry.
         */
        CA.clearGeometry();
        HA.clearGeometry();
        C.clearGeometry();
        O.clearGeometry();
        N.clearGeometry();
        H.clearGeometry();
        /**
         * Change their residue name.
         */
        String resName = toResidue.getName();
        CA.setResName(resName);
        HA.setResName(resName);
        C.setResName(resName);
        O.setResName(resName);
        N.setResName(resName);
        H.setResName(resName);
        /**
         * Add the backbone atoms to the new Residue.
         */
        toResidue.addMSNode(CA);
        toResidue.addMSNode(HA);
        toResidue.addMSNode(C);
        toResidue.addMSNode(O);
        toResidue.addMSNode(N);
        toResidue.addMSNode(H);
    }

    /**
     * Update Atom references to local geometry.
     *
     * @param residue
     */
    private void updateGeometry(Residue residue, Residue prev, Residue next,
            Residue prev2, Residue next2) {
        if (residue == null) {
            return;
        }
        /**
         * Update atom references to local geometry.
         */
        ArrayList<Atom> atoms = residue.getAtomList();
        ArrayList<ROLS> bonds = residue.getBondList();
        ArrayList<ROLS> angles = residue.getAngleList();
        ArrayList<ROLS> torsions = residue.getTorsionList();
        if (prev != null) {
            atoms.addAll(prev.getAtomList());
            bonds.addAll(prev.getBondList());
            angles.addAll(prev.getAngleList());
            torsions.addAll(prev.getTorsionList());
            ArrayList<Joint> joints = prev.getJoints();
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
            ArrayList<Joint> joints = next.getJoints();
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
            for (ROLS bond : bonds) {
                Bond b = (Bond) bond;
                if (b.containsAtom(atom)) {
                    atom.setBond(b);
                }
            }
        }
        for (Atom atom : atoms) {
            for (ROLS angle : angles) {
                Angle a = (Angle) angle;
                if (a.containsAtom(atom)) {
                    atom.setAngle(a);
                }
            }
        }
        for (Atom atom : atoms) {
            for (ROLS torsion : torsions) {
                Torsion t = (Torsion) torsion;
                if (t.containsAtom(atom)) {
                    atom.setTorsion(t);
                }
            }
        }
    }

    public void addResidue(Residue newResidue) {
        /**
         * Add the new residue to list.
         */
        consideredResidues.add(newResidue);
        /**
         * Get references to nearby residues.
         */
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

        /**
         * Move atoms from the active Residue to the new Residue.
         */
        moveBackBoneAtoms(activeResidue, newResidue);
        /**
         * Pass references of the active Residues' joints to the new Residue.
         */
        ArrayList<Joint> joints = activeResidue.getJoints();
        for (Joint joint : joints) {
            newResidue.addJoint(joint);
        }
        /**
         * Make the new Residue active.
         */
        activeResidue.removeFromParent();
        activeResidue = newResidue;
        add(activeResidue);
        /**
         * Build side-chain atoms and assign atom types for the new Residue.
         */
        try {
            assignAminoAcidAtomTypes(newResidue, prevResidue, nextResidue, forceField, null);
            if (nextResidue != null) {
                Atom C = (Atom) newResidue.getAtomNode("C");
                Atom nextN = (Atom) nextResidue.getAtomNode("N");
                for (Joint joint : joints) {
                    Bond bond = (Bond) joint.getBondList().get(0);
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

    public boolean setActiveResidue(int i) {
        if (consideredResidues == null) {
            return false;
        }
        if (i >= consideredResidues.size()) {
            return false;
        }
        return setActiveResidue(consideredResidues.get(i));
    }

    public boolean setActiveResidue(Residue residue) {
        if (!consideredResidues.contains(residue)) {
            return false;
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

        /**
         * Move backbone atoms to the new active residue.
         */
        moveBackBoneAtoms(activeResidue, residue);
        updateGeometry(residue, prevResidue, nextResidue, prev2Residue, next2Residue);
        activeResidue = residue;
        add(activeResidue);

        return true;
    }

    public int getResidueCount() {
        if (consideredResidues == null) {
            return 0;
        }
        return consideredResidues.size();
    }

    @Override
    public String toString() {
        if (activeResidue == null) {
            return null;
        } else {
            return activeResidue.toString();
        }
    }

}
