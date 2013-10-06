/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
import java.util.logging.Logger;

import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.RotamerLibrary;
import ffx.potential.Rotamer;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.PDBFilter.ResiduePosition;
import ffx.potential.parsers.PDBFilter.MissingHeavyAtomException;
import java.util.List;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

/**
 *
 * @author wtollefson
 */
public class MultiResidue extends Residue {

    private static final Logger logger = Logger.getLogger(MultiResidue.class.getName());
    /**
     * Which residue is active.
     */
    private Residue activeResidue = null;
    /**
     * List of residues under consideration.
     */
    ArrayList<Residue> consideredResidues;
    PDBFilter pdbFilter;

    public MultiResidue(Residue residue, MolecularAssembly molecularAssembly) {
        super("MultiResidue", residue.getResidueNumber(), residue.residueType);
        pdbFilter = new PDBFilter(molecularAssembly.getFile(), molecularAssembly, molecularAssembly.getForceField(), null);
        activeResidue = residue;
        // Initialize consideredResidue list.
        consideredResidues = new ArrayList<Residue>();
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
    public void collectValenceTerms() {
        activeResidue.collectValenceTerms();
    }

    @Override
    public void constructValenceTerms() {
        activeResidue.constructValenceTerms();
    }

    @Override
    public Joint createJoint(Bond bond, MSGroup group1, MSGroup group2) {
        return activeResidue.createJoint(bond, group1, group2);
    }

    @Override
    public Joint createJoint(MSGroup group1, MSGroup group2) {
        return activeResidue.createJoint(group1, group2);
    }

    @Override
    public void finalize(boolean finalizeGeometry) {
        activeResidue.finalize(finalizeGeometry);
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

    public void addResidue(Residue residue) {
        int number = residue.getResidueNumber();
        ResiduePosition position = pdbFilter.getResiduePosition(number);
        AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
        // Create copies of CA, N, C

        Atom CA = (Atom) activeResidue.getAtomNode("CA");
        Atom HA = (Atom) activeResidue.getAtomNode("HA");
        Atom C = (Atom) activeResidue.getAtomNode("C");
        Atom O = (Atom) activeResidue.getAtomNode("O");
        Atom N = (Atom) activeResidue.getAtomNode("N");
        Atom H = (Atom) activeResidue.getAtomNode("H");
        Atom newN = N.copy();
        newN.setResName(residue.getName());
        Atom newH = H.copy();
        newH.setResName(residue.getName());
        Atom newCA = CA.copy();
        newCA.setResName(residue.getName());
        Atom newHA = HA.copy();
        newHA.setResName(residue.getName());
        Atom newC = C.copy();
        newC.setResName(residue.getName());
        Atom newO = O.copy();
        newO.setResName(residue.getName());
        pdbFilter.buildBond(newN, newH);
        pdbFilter.buildBond(newN, newCA);
        pdbFilter.buildBond(newCA, newHA);
        pdbFilter.buildBond(newCA, newC);
        pdbFilter.buildBond(newC, newO);
        // Add them to residue
        residue.addMSNode(newN);
        residue.addMSNode(newH);
        residue.addMSNode(newCA);
        residue.addMSNode(newHA);
        residue.addMSNode(newC);
        residue.addMSNode(newO);
        try {
            pdbFilter.assignAminoAcidSideChain(position, name, residue, CA, N, C);
            add(residue);
            residue.finalize(true);
        } catch (MissingHeavyAtomException missingHeavyAtomException) {
            logger.severe(missingHeavyAtomException.toString());
        }

    }
}
