/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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

import java.util.logging.Logger;

import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.BiojavaFilter;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.swing.tree.TreeNode;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;

/**
 * The Molecule class is a general container used for simple compounds or in
 * cases where more specialized classes have not been implemented.
 *
 * @author Michael J. Schnieders
 *
 */
public class Molecule extends MSGroup implements Group {

    private Logger logger = Logger.getLogger(Molecule.class.getName());
    private static final long serialVersionUID = 1L;
    /**
     * Constant <code>MultiScaleLevel=2</code>
     */
    public static final int MultiScaleLevel = 2;
    /**
     * Residue number assigned in PDB files.
     */
    private int residueNum = -1;
    /**
     * Residue name assigned in PDB files.
     */
    private String residueName = null;
    /**
     * Possibly redundant chainID assigned in PDB files.
     */
    private Character chainID = null;
    /**
     * Unique segID.
     */
    private String segID = null;
    /**
     * The Chain to which this Molecule belongs.
     */
    private Chain parentChain;
    /**
     * String-mapped Biojava-related properties.
     */
    private Map<String, Object> properties;
    /**
     * List of Molecules matching alternative locations.
     */
    private Map<Character, Group> altLocGroups;
    /**
     * Biojava residue identifier.
     */
    private ResidueNumber resNum;
    /**
     * Chemical component definition.
     */
    private ChemComp chemComp;

    /**
     * <p>
     * Constructor for Molecule.</p>
     */
    public Molecule() {
    }

    /**
     * <p>
     * Constructor for Molecule.</p>
     *
     * @param name a {@link java.lang.String} object.
     */
    public Molecule(String name) {
        super(name);
        residueName = name;
        chainID = 'A';
    }

    /**
     * <p>
     * Constructor for Molecule.</p>
     *
     * @param name a {@link java.lang.String} object.
     * @param residueNum a int.
     * @param chainID a {@link java.lang.Character} object.
     * @param segID a {@link java.lang.String} object.
     */
    public Molecule(String name, int residueNum,
            Character chainID, String segID) {
        super(name + "-" + residueNum + " " + segID);
        this.residueName = name;
        this.residueNum = residueNum;
        this.chainID = chainID;
        this.segID = segID;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setName(String name) {
        if (segID != null) {
            super.setName(name + "-" + residueNum + " " + segID);
        } else {
            super.setName(name);
        }
        this.residueName = name;
    }

    /**
     * <p>
     * Getter for the field <code>residueName</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getResidueName() {
        return residueName;
    }

    /**
     * <p>
     * getResidueNumber</p>
     *
     * @return a int.
     */
    public int getResidueIndex() {
        return residueNum;
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
     * <p>
     * Getter for the field <code>segID</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getSegID() {
        return segID;
    }

    /**
     * <p>
     * getAtom</p>
     *
     * @param name a {@link java.lang.String} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getMoleculeAtom(String name) {
        for (Atom a : getAtomList()) {
            if (a.getName().equalsIgnoreCase(name)) {
                return a;
            }
        }
        return null;
    }

    /**
     * {@inheritDoc}
     *
     * Allows adding Atom MSNodes to the Molecule.
     */
    @Override
    public MSNode addMSNode(MSNode o) {
        Atom currentAtom = null;
        if (o instanceof Atom) {
            Atom newAtom = (Atom) o;
            Character newAlt = newAtom.getAltLoc();
            MSNode atoms = getAtomNode();
            currentAtom = (Atom) atoms.contains(newAtom);
            if (currentAtom == null) {
                newAtom.setGroup(this, true);
                currentAtom = newAtom;
                atoms.add(newAtom);
                setFinalized(false);
            } else {
                /**
                 * Allow overwriting of the root alternate conformer (' ' or
                 * 'A').
                 */
                Character currentAlt = currentAtom.getAltLoc();
                if (currentAlt.equals(' ') || currentAlt.equals('A')) {
                    if (!newAlt.equals(' ') && !newAlt.equals('A')) {
                        newAtom.xyzIndex = currentAtom.xyzIndex;
                        atoms.remove(currentAtom);
                        currentAtom = newAtom;
                        atoms.add(currentAtom);
                        setFinalized(false);
                    }
                }
            }
        } else {
            logger.warning("Only an Atom can be added to a Residue.");
        }
        return currentAtom;
    }

    /**
     * {@inheritDoc}
     */
    public void finalize(boolean finalizeGeometry, ForceField forceField) {
        setFinalized(false);
        getAtomNode().setName("Atoms (" + getAtomList().size() + ")");
        if (finalizeGeometry) {
            //constructValenceTerms();
            assignBondedTerms(forceField);
            removeLeaves();
        }
        // findDangelingAtoms();
        setCenter(getMultiScaleCenter(false));
        setFinalized(true);
    }
    
    @Override
    public void setChain(Chain chain) {
        if (parentChain instanceof Polymer) {
            removeFromParent();
        }
        if (chain instanceof Polymer) {
            ((Polymer) chain).addMSNode(this);
        }
        this.parentChain = chain;
    }
    
    /**
     * Sets the parent Chain; if onlySetRef is true, does not detach or attach to
     * FFX data structure.
     * @param chain Chain to set parentChain reference to
     * @param onlySetRef If true, only shallowly sets reference
     */
    public void setChain(Chain chain, boolean onlySetRef) {
        if (onlySetRef) {
            this.parentChain = chain;
        } else {
            setChain(chain);
        }
    }
    
    /**
     * Returns a copy with basic information (force field information generally
     * not copied).
     * @return 
     */
    public Molecule clone() {
        Molecule ret = new Molecule(getName(), residueNum, chainID, segID);
        for (Atom atom : getAtomList()) {
            Atom newAtom = atom.clone();
            ret.addMSNode(newAtom);
        }
        ret.setChemComp(chemComp);
        ret.setResidueNumber(copyResNum());
        return ret;
    }
    
    public void findParentPolymer() {
        TreeNode parentNode = getParent();
        while (parentNode != null) {
            if (parentNode instanceof Polymer) {
                this.parentChain = (Chain) parentChain;
                break;
            } else {
                parentNode = parentNode.getParent();
            }
        }
    }
    
    public Chain getParentChain() {
        if (parentChain == null) {
            findParentPolymer();
        }
        return parentChain;
    }
    
    @Override
    public int size() {
        return getAtomList().size();
    }

    /**
     * FFX data structures have 3D coordinates by default.
     * @return Always true.
     */
    @Override
    public boolean has3D() {
        return true;
    }

    /**
     * Does nothing (FFX data structures have 3D coordinates by default).
     * @param bln Discarded.
     */
    @Override
    public void setPDBFlag(boolean bln) {
        logger.fine(" FFX atoms always have coordinates; setPDBFlag is meaningless.");
        // throw new UnsupportedOperationException("Force Field X atoms always have coordinates");
    }

    @Override
    public GroupType getType() {
        return org.biojava.nbio.structure.GroupType.HETATM;
    }

    @Override
    public void addAtom(org.biojava.nbio.structure.Atom atom) {
        if (atom instanceof Atom) {
            addMSNode((Atom) atom);
        } else {
            Atom at = BiojavaFilter.readAtom(atom, segID);
            addMSNode(at);
            /*Polymer parentPolymer = (Polymer) parentChain;
            if (parentPolymer.hasFFXParents()) {
                parentPolymer.addExteriorAtom(atom);
            }*/
        }
    }

    @Override
    public List<org.biojava.nbio.structure.Atom> getAtoms() {
        List<org.biojava.nbio.structure.Atom> retList = new ArrayList<>();
        retList.addAll(getAtomList());
        return retList;
    }

    @Override
    public void setAtoms(List<org.biojava.nbio.structure.Atom> list) {
        clearAtoms();
        for (org.biojava.nbio.structure.Atom atom : list) {
            try {
                this.addAtom(atom);
            } catch (IllegalArgumentException ex) {
                logger.fine(String.format(" Failure to add atom %s", atom.toString()));
            }
        }
    }

    @Override
    public void clearAtoms() {
        List<Atom> atoms = this.getAtomList();
        for (Atom atom : atoms) {
            atom.setGroup(null, true);
            this.remove(atom);
        }
    }

    @Override
    public org.biojava.nbio.structure.Atom getAtom(String string) {
        Atom atom = (Atom) this.getAtomNode(string);
        if (atom != null) {
            return atom;
        }
        return null;
    }

    @Override
    public org.biojava.nbio.structure.Atom getAtom(int i) {
        return (org.biojava.nbio.structure.Atom) this.getAtomNode(i);
    }

    @Override
    public boolean hasAtom(String string) {
        return (getMoleculeAtom(string) != null);
    }

    @Override
    public String getPDBName() {
        return getName();
    }

    @Override
    public void setPDBName(String string) {
        setName(string);
    }

    @Override
    public boolean hasAminoAtoms() {
		// if this method call is performed too often, it should become a
        // private method and provide a flag for Group object ...

        return hasAtom(StructureTools.CA_ATOM_NAME)
                && hasAtom(StructureTools.C_ATOM_NAME)
                && hasAtom(StructureTools.N_ATOM_NAME)
                && hasAtom(StructureTools.O_ATOM_NAME);

    }

    @Override
    public void setProperties(Map<String, Object> map) {
        this.properties = map;
    }

    @Override
    public Map<String, Object> getProperties() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setProperty(String string, Object o) {
        properties.put(string, o);
    }

    @Override
    public Object getProperty(String string) {
        return properties.get(string);
    }

    @Override
    public Iterator<org.biojava.nbio.structure.Atom> iterator() {
        return new AtomIterator(this);
    }
    
    @Override
    public Chain getChain() {
        if (parentChain == null) {
            findParentPolymer();
        }
        return parentChain;
    }

    @Override
    public ResidueNumber getResidueNumber() {
        if (resNum == null) {
            generateResNum();
        }
        return resNum;
    }
    
    private void generateResNum() {
        char insCode = ' ';
        for (Atom atom : getAtomList()) {
            if (atom.getInsertionCode() != ' ') {
                insCode = atom.getInsertionCode();
                break;
            }
        }
        resNum = new ResidueNumber("" + chainID, residueNum, insCode);
    }
    
    private ResidueNumber copyResNum() {
        if (resNum == null) {
            generateResNum();
        }
        return new ResidueNumber(resNum.getChainId(), resNum.getSeqNum(), resNum.getInsCode());
    }

    @Override
    public void setResidueNumber(ResidueNumber rn) {
        this.resNum = rn;
    }

    @Override
    public void setResidueNumber(String chnID, Integer rNum, Character iCode) {
        resNum = new ResidueNumber(chnID, rNum, iCode);
    }

    @Override
    public String getChainId() {
        return "" + chainID;
    }

    @Override
    public void setChemComp(ChemComp cc) {
        this.chemComp = cc;
    }

    @Override
    public ChemComp getChemComp() {
        return chemComp;
    }

    @Override
    public boolean hasAltLoc() {
        /*if (altLocGroups == null) {
            findAltLocs();
        }
        return !(altLocGroups.isEmpty());*/
        return false;
    }
    
    private void findAltLocs() {
        // TO BE IMPLEMENTAZORLALIZATIONED
    }

    @Override
    public List<Group> getAltLocs() {
        /*List<Group> altLocs = new ArrayList<>();
        if (altLocGroups == null) {
            findAltLocs();
        }
        altLocs.addAll(altLocGroups.values());
        return altLocs;*/
        throw new UnsupportedOperationException("FFX does not yet support references"
                + " to alternate locations.");
    }

    @Override
    public void addAltLoc(Group group) {
        /*if (altLocGroups == null) {
            findAltLocs();
        }
        boolean altLocFound = false;
        for (org.biojava.nbio.structure.Atom atom : group.getAtoms()) {
            char aLoc = atom.getAltLoc();
            if (atom.getAltLoc() != ' ' && !altLocGroups.containsKey(aLoc)) {
                altLocGroups.put(aLoc, group);
                altLocFound = true;
                break;
            }
        }
        if (!altLocFound) {
            for (char alpha = 'A'; alpha <= 'Z'; alpha++) {
                if (!altLocGroups.containsKey(alpha)) {
                    altLocGroups.put(alpha, group);
                    for (org.biojava.nbio.structure.Atom atom : group.getAtoms()) {
                        atom.setAltLoc(alpha);
                    }
                    logger.info(String.format(" Alternate location group %s does "
                            + "not have a unique alternate location code; setting"
                            + "the group to altloc %c", group.toString(), alpha));
                }
            }
        }*/
        throw new UnsupportedOperationException("FFX does not yet support references"
                + " to alternate locations.");
    }

    @Override
    public boolean isWater() {
        return GroupType.WATERNAMES.contains(getName());
    }

    @Override
    public Group getAltLocGroup(Character aLoc) {
        /*if (altLocGroups == null) {
            findAltLocs();
        }
        return altLocGroups.get(aLoc);*/
        throw new UnsupportedOperationException("FFX does not yet support references"
                + " to alternate locations.");
    }

    @Override
    public void trimToSize() {
        logger.fine(" Operation trimToSize() not yet supported.");
    }
}
