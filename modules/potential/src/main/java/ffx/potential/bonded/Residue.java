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

import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Iterator;
import java.util.Map;

import javax.media.j3d.Canvas3D;
import javax.media.j3d.J3DGraphics2D;
import javax.media.j3d.Material;
import javax.media.j3d.Node;
import javax.vecmath.Color3f;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import javax.swing.tree.TreeNode;

import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.BiojavaFilter;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;

/**
 * The Residue class represents individual amino acids or nucleic acid bases.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public class Residue extends MSGroup implements Group {

    private static final Logger logger = Logger.getLogger(Residue.class.getName());

    /**
     * The residue number of this residue in a chain.
     */
    private int resIndex;
    /**
     * Possibly redundant PDB chain ID.
     */
    private Character chainID;
    /**
     * The Polymer to which this Residue belongs.
     */
    private Chain parentChain;
    /**
     * String-mapped Biojava-related properties.
     */
    private Map<String, Object> properties;
    /**
     * List of Residues matching alternative locations.
     */
    //private Map<Character, Group> altLocGroups;
    /**
     * Biojava residue identifier.
     */
    private ResidueNumber resNum;
    /**
     * Chemical component definition
     */
    private ChemComp chemComp;
    /**
     * Unique segID.
     */
    private String segID;
    /**
     * Secondary structure type.
     */
    private SSType ssType = SSType.NONE;
    /**
     * Residue type.
     */
    protected ResidueType residueType = ResidueType.UNK;
    private AA3 aa;
    private NA3 na;
    /**
     * These arrays store default coordinates for certain atoms in nucleic acid
     * Residues. C1', O4', and C4' are the critical sugar atoms off which every
     * other atom is drawn when applyRotamer is called. The backbone
     * corrections, however, move these atoms, so they must be reverted to their
     * original coordinates each time applyRotamer is called.
     *
     * O3' North and South coordinates are technically non-essential, as they
     * could be derived from C1', O4', C4', and a given sugar pucker, however,
     * it is much less computationally expensive to calculate them once and then
     * store them.
     *
     * TODO: Add O3' coordinates for the DNA C3'-exo configuration.
     */
    private double[] O3sNorthCoords = null;
    private double[] O3sSouthCoords = null;
    private double[] C1sCoords = null;
    private double[] O4sCoords = null;
    private double[] C4sCoords = null;

    private Rotamer currentRotamer = null;
    private Rotamer originalRotamer = null;
    protected static final boolean origAtEnd;
    
    static {
        String origAtEndStr = System.getProperty("ro-origAtEnd");
        if (origAtEndStr != null) {
            origAtEnd = Boolean.parseBoolean(origAtEndStr);
        } else {
            origAtEnd = false;
        }
    }

    /**
     * Default Constructor where num is this Residue's position in the Polymer.
     *
     * @param num a int.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(int num, ResidueType rt) {
        super();
        resIndex = num;
        residueType = rt;
        assignResidueType();
    }

    /**
     * <p>
     * Constructor for Residue.</p>
     *
     * @param name a {@link java.lang.String} object.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(String name, ResidueType rt) {
        super(name);
        residueType = rt;
        assignResidueType();
    }

    /**
     * Name is the residue's 3 letter abbreviation and num is its position in
     * the Polymer.
     *
     * @param name a {@link java.lang.String} object.
     * @param num a int.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(String name, int num, ResidueType rt) {
        this(name, rt);
        resIndex = num;
    }

    /**
     * Name is the residue's 3 letter abbreviation and num is its position in
     * the Polymer.
     *
     * @param name a {@link java.lang.String} object.
     * @param resNumber a int.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     * @param chainID a {@link java.lang.Character} object.
     * @param segID a {@link java.lang.String} object.
     */
    public Residue(String name, int resNumber, ResidueType rt, Character chainID,
            String segID) {
        this(name, rt);
        this.resIndex = resNumber;
        this.chainID = chainID;
        this.segID = segID;
    }

    /**
     * As above, with atoms being a FNode with this Residue's atoms as child
     * nodes
     *
     * @param name a {@link java.lang.String} object.
     * @param num a int.
     * @param atoms a {@link ffx.potential.bonded.MSNode} object.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     * @param forceField the ForceField to use when created bonded terms.
     */
    public Residue(String name, int num, MSNode atoms, ResidueType rt, ForceField forceField) {
        super(name, atoms);
        resIndex = num;
        residueType = rt;
        assignResidueType();
        finalize(true, forceField);
    }

    /*public Rotamer[] getRotamers(Residue residue) {
        return RotamerLibrary.getRotamers(residue);
    }*/
    
    public Rotamer[] getRotamers() {
        if (RotamerLibrary.getUsingOrigCoordsRotamer()) {
            Rotamer[] libRotamers = RotamerLibrary.getRotamers(this);
            if (libRotamers == null) {
                return null;
            }
            int nRots = libRotamers.length;
            Rotamer[] rotamers = new Rotamer[nRots + 1];
            if (originalRotamer == null) {
                double[][] origCoordinates = storeCoordinateArray();
                double[] chi = RotamerLibrary.measureRotamer(this, false);
                switch (residueType) {
                    case AA:
                        AminoAcid3 aa = AminoAcid3.valueOf(getName());
                        originalRotamer = new Rotamer(aa, origCoordinates, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0);
                        break;
                    case NA:
                        NucleicAcid3 na = NucleicAcid3.valueOf(getName());
                        originalRotamer = new Rotamer(na, origCoordinates, chi[0], 0, chi[1], 0, chi[2], 0, chi[3], 0, chi[4], 0, chi[5], 0);
                        break;
                    default:
                        originalRotamer = null;
                        return libRotamers;
                }
            }
            if (origAtEnd) {
                System.arraycopy(libRotamers, 0, rotamers, 0, nRots);
                rotamers[rotamers.length - 1] = originalRotamer;
            } else {
                System.arraycopy(libRotamers, 0, rotamers, 1, nRots);
                rotamers[0] = originalRotamer;
            }
            return rotamers;
        } else {
            return RotamerLibrary.getRotamers(this);
        }
    }

    public ResidueType getResidueType() {
        return residueType;
    }
    
    public boolean isDeoxy() {
        if (getResidueType() == ResidueType.NA) {
            Atom HOs = (Atom) getAtomNode("HO\'");
            if (HOs == null) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the Residue bonded to this Residue at this Residue's 3' or
     * C-terminal end. Any use of this method to add Residues to a sliding
     * window or similar MUST not add that residue if that residue has no
     * Rotamers, as several algorithms (such as the distance matrix) assume that
     * all Residues being optimized have Rotamers.
     *
     * @return The next Residue.
     */
    public Residue getNextResidue() {
        switch (residueType) {
            case AA: {
                Atom carbon = (Atom) getAtomNode("C");
                if (carbon == null) {
                    return null;
                }
                ArrayList<Bond> bonds = carbon.getFFXBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(carbon);
                    if (other.getName().equalsIgnoreCase("N")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            }
            case NA:
                Atom oxygen = (Atom) getAtomNode("O3\'");
                if (oxygen == null) {
                    return null;
                }
                ArrayList<Bond> bonds = oxygen.getFFXBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(oxygen);
                    if (other.getName().equalsIgnoreCase("P")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            default:
                return null;
        }
        return null;
        // Will generally indicate that you passed in a chain-terminal residue.
    }

    /**
     * Returns the Residue bonded to this Residue at this Residue's 5' or
     * N-terminal end. Any use of this method to add Residues to a sliding
     * window or similar MUST not add that residue if that residue has no
     * Rotamers, as several algorithms (such as the distance matrix) assume that
     * all Residues being optimized have Rotamers.
     *
     * @return The previous Residue.
     */
    public Residue getPreviousResidue() {
        switch (residueType) {
            case AA: {
                Atom nitrogen = (Atom) getAtomNode("N");
                if (nitrogen == null) {
                    return null;
                }
                ArrayList<Bond> bonds = nitrogen.getFFXBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(nitrogen);
                    if (other.getName().equalsIgnoreCase("C")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            }
            case NA:
                Atom phosphate = (Atom) getAtomNode("P");
                if (phosphate == null) {
                    return null;
                }
                ArrayList<Bond> bonds = phosphate.getFFXBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(phosphate);
                    if (other.getName().equalsIgnoreCase("O3\'")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            default:
                return null;
        }
        return null;
        // Will generally indicate that you passed in a chain-starting residue.
    }

    /**
     * Returns a reference Atom for a Residue, primarily intended for rough
     * distance calculations. This atom should be roughly centrally located
     * within the residue, and be invariant.
     *
     * @return A reference Atom.
     */
    public Atom getReferenceAtom() {
        Atom atom = null;
        switch (this.getResidueType()) {
            case AA:
                atom = (Atom) this.getAtomNode("CA");
                break;
            case NA:
                // If pyrimidine, atom will be N1.  Else, if purine,
                // N1 will return null, so grab N9.
                atom = (Atom) this.getAtomNode("N1");
                if (atom == null) {
                    atom = (Atom) this.getAtomNode("N9");
                }
                break;
            default:
                break;
        }
        if (atom == null) {
            atom = (Atom) this.getAtomNode(0);
        }
        return atom;
    }

    public ResidueState storeCoordinates() {
        return new ResidueState(this, this);
    }

    public void revertCoordinates(ResidueState state) {
        List<Atom> atomList = getAtomList();
        for (Atom atom : atomList) {
            atom.moveTo(state.getAtomCoords(atom));
        }
    }

    public double[][] storeCoordinateArray() {
        List<Atom> atomList = getAtomList();
        int nAtoms = atomList.size();
        double[][] coords = new double[nAtoms][3];
        int i = 0;
        for (Atom atom : atomList) {
            atom.getXYZ(coords[i++]);
        }
        return coords;
    }

    public void setRotamer(Rotamer rotamer) {
        this.currentRotamer = rotamer;
    }

    public Rotamer getRotamer() {
        return currentRotamer;
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
    public Residue clone() {
        Residue ret = new Residue(getName(), resIndex, residueType, chainID, segID);
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
        switch (residueType) {
            case AA:
                return org.biojava.nbio.structure.GroupType.AMINOACID;
            case NA:
                return org.biojava.nbio.structure.GroupType.NUCLEOTIDE;
            default:
                return org.biojava.nbio.structure.GroupType.HETATM;
        }
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

    /**
     * May not be reliable, as FFX data structure may reorder atoms.
     * @param i
     * @return 
     */
    @Override
    public org.biojava.nbio.structure.Atom getAtom(int i) {
        return this.getAtomList().get(i);
        //return (org.biojava.nbio.structure.Atom) this.getAtomNode(i);
    }

    @Override
    public boolean hasAtom(String string) {
        return (getAtom(string) != null);
    }

    @Override
    public String getPDBName() {
        String name = getName();
        AminoAcid3 aa3 = AminoAcid3.valueOf(name);
        switch (aa3) {
            case GLH:
                return "GLU";
            case HIE:
            case HID:
                return "HIS";
            case CYD:
                return "CYS";
            case TYD:
                return "TYR";
            case ASH:
                return "ASP";
            case LYD:
                return "LYS";
            default:
                return name;
        }
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
        return properties;
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
    
    /**
     * Generates a ResidueNumber object to act as unique identifier.
     */
    private void generateResNum() {
        char insCode = ' ';
        for (Atom atom : getAtomList()) {
            if (atom.getInsertionCode() != ' ') {
                insCode = atom.getInsertionCode();
                break;
            }
        }
        resNum = new ResidueNumber("" + chainID, resIndex, insCode);
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
        /*throw new UnsupportedOperationException("FFX does not yet support references"
                + " to alternate locations.");*/
        return false;
    }
    
    private void findAltLocs() {
        // TO BE IMPLEMENTAZORLALIZATIONED
    }
    
    /**
     * As FFX treats alternate locations on a per-atom basis, not a per-residue
     * basis, this method was developed to create generic Groups in a similar
     * fashion.
     */
    /*@Override
    public void updateAltLocs() {
        List<Atom> allAtoms = this.getAtomList();
        altLocGroups = new HashMap<>();

        List<Atom> nonAltAtoms = new ArrayList<>();
        for (Atom atom : allAtoms) {
            char altLoc = atom.getAltLoc();
            if (altLoc != ' ') {
                if (altLocGroups.containsKey(altLoc)) {
                    altLocGroups.get(altLoc).addAtom(atom);
                } else {
                    Group altGroup;
                    switch (this.residueType) {
                        case NA:
                            altGroup = new NucleotideImpl();
                            break;
                        case AA:
                            altGroup = new AminoAcidImpl();
                            break;
                        default:
                            altGroup = new HetatomImpl();
                            break;
                    }
                    altGroup.addAtom(atom);
                    altGroup.setChain(this.getChain());
                    altGroup.setPDBFlag(true);
                    altGroup.setPDBName(getName());
                    altGroup.setProperties(properties);
                    altGroup.setResidueNumber(getResidueIndex());
                    altLocGroups.put(altLoc, altGroup);
                }
            } else {
                nonAltAtoms.add(atom);
            }
        }
        if (!altLocGroups.isEmpty()) {
            for (Group altGroup : altLocGroups.values()) {
                for (Atom atom : nonAltAtoms) {
                    altGroup.addAtom(atom);
                }
            }
        }
    }*/

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

    /**
     * Always returns false (a Residue is never a water; only a Molecule can be 
     * a water).
     * @return False.
     */
    @Override
    public boolean isWater() {
        return false;
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

    /**
     * {@inheritDoc}
     *
     * Allows adding Atoms to the Residue.
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
     * <p>
     * deleteAtom</p>
     *
     * @param atomToDelete a {@link ffx.potential.bonded.Atom} object.
     */
    public void deleteAtom(Atom atomToDelete) {
        MSNode atoms = getAtomNode();
        if (atoms.contains(atomToDelete) != null) {
            logger.info(" The following atom is being deleted from the model:\n"
                    + atomToDelete.toString());
            atoms.remove(atomToDelete);
        }
    }

    /**
     * Returns a list of side chain atoms; for our purposes, nucleic acid side
     * chain atoms are the sugar and the phosphate.
     *
     * @return ArrayList of side chain (or nucleic backbone) atoms.
     */
    public ArrayList<Atom> getSideChainAtoms() {
        ArrayList<Atom> atoms = getAtomList();
        ArrayList<Atom> ret;
        switch (residueType) {
            case NA:
                ret = new ArrayList<>();
                for (Atom atom : atoms) {
                    String name = atom.getName().toUpperCase();
                    /*
                     * Very conveniently for our purposes, the entire sugar is
                     * denoted with ' at the end. Note that we add the side chain
                     * atoms for NAs, instead of removing the backbone.
                     */
                    if (name.contains("\'") || name.equals("P") || name.startsWith("OP") || name.equals("H5T")
                            || name.equals("H3T")) {
                        ret.add(atom);
                    }
                }
                return ret;
            case AA:
                ret = new ArrayList<>(atoms);
                for (Atom atom : atoms) {
                    String name = atom.getName().toUpperCase();
                    if (name.equals("N") || name.equals("H") || name.equals("H1") || name.equals("H2") || name.equals("H3")
                            || name.equals("CA") || name.startsWith("HA")
                            || name.equals("C") || name.equals("O") || name.equals("OXT") || name.equals("OT2")) {
                        ret.remove(atom);
                    }
                }
                return ret;
            default:
                return null;
        }
    }
    
    /**
     * Returns a list of atoms liable to change during dead-end elimination repacking.
     * For ordinary amino acids: side chain atoms. For ordinary nucleic acids:
     * sugar/phosphate backbone atoms. MultiResidue over-rides this to return all
     * atoms (as backbone atom types are nonconstant).
     * 
     * @return Atoms changeable during DEE.
     */
    public List<Atom> getVariableAtoms() {
        return getSideChainAtoms();
    }

    /**
     * Returns a list of backbone atoms; for our purposes, nucleic acid backbone
     * atoms are those of the nucleobase. Protein backbone atoms will be
     * ordered:
     *
     *
     * @return ArrayList of backbone (or nucleobase) atoms.
     */
    public ArrayList<Atom> getBackboneAtoms() {
        ArrayList<Atom> atoms = getAtomList();
        ArrayList<Atom> ret;
        switch (residueType) {
            case NA:
                ret = new ArrayList<>(atoms);
                for (Atom atom : atoms) {
                    String name = atom.getName().toUpperCase();
                    if (name.contains("\'") || name.equals("P") || name.startsWith("OP") || name.equals("H5T")
                            || name.equals("H3T")) {
                        ret.remove(atom);
                    }
                }
                return ret;
            case AA:
                ret = new ArrayList<>();
                tryAddAtom(ret, "N");
                tryAddAtom(ret, "CA");
                tryAddAtom(ret, "C");
                tryAddAtom(ret, "O");
                tryAddAtom(ret, "OXT"); // C-terminal residues
                tryAddAtom(ret, "OT2"); // Probably alternate name for OXT.
                tryAddAtom(ret, "H1"); // N-terminal residues
                tryAddAtom(ret, "H2");
                tryAddAtom(ret, "H3");
                tryAddAtom(ret, "H");
                tryAddAtom(ret, "HA");
                tryAddAtom(ret, "HA2"); // Glycines
                tryAddAtom(ret, "HA3");
                return ret;
            default:
                return null;
        }
    }

    /**
     * Uses a name to add an Atom to a List<Atom> if the Atom exists for this
     * residue.
     *
     * @param atList List to add to.
     * @param name Atom to add.
     * @return If atom exists.
     */
    private boolean tryAddAtom(List<Atom> atList, String name) {
        try {
            Atom at = (Atom) getAtomNode(name);
            if (at != null) {
                atList.add(at);
                return true;
            } else {
                return false;
            }
        } catch (Exception ex) {
            return false;
        }
    }

    private void assignResidueType() {
        String name = getName().toUpperCase();
        switch (residueType) {
            case AA:
                aa = null;
                try {
                    if (name.length() >= 2) {
                        aa = AA3.valueOf(name);
                    } else if (name.length() == 1) {
                        AA1 aa1 = AA1.valueOf(name);
                        aa = AA1toAA3.get(aa1);
                    }
                } catch (Exception e) {
                    aa = AA3.UNK;
                }
                break;
            case NA:
                na = null;
                try {
                    if (name.length() >= 2) {
                        na = NA3.valueOf(name);
                    } else if (name.length() == 1) {
                        NA1 na1 = NA1.valueOf(name);
                        na = NA1toNA3.get(na1);
                    }
                } catch (Exception e) {
                    na = NA3.UNK;
                }
                break;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void drawLabel(Canvas3D canvas, J3DGraphics2D g2d, Node node) {
        if (RendererCache.labelResidues) {
            double d[] = getCenter();
            point3d.x = d[0];
            point3d.y = d[1];
            point3d.z = d[2];
            RendererCache.getScreenCoordinate(canvas, node, point3d, point2d);
            g2d.drawString(getName(), (float) point2d.x, (float) point2d.y);
        }
        if (RendererCache.labelAtoms) {
            super.drawLabel(canvas, g2d, node);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Overidden equals method that return true if object is not equals to this,
     * is of the same class, has the same parent Polymer, the same sequence
     * number, the same ResidueType, and the same AA3/NA3.
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        Residue other = (Residue) object;
        ResidueType otherType = other.getResidueType();
        if (residueType != otherType) {
            return false;
        }
        switch (residueType) {
            case AA:
                if (aa != other.aa) {
                    return false;
                }
                break;
            case NA:
                if (na != other.na) {
                    return false;
                }
                break;
            default:
                break;
        }
        if (!getName().equals(other.getName())) {
            return false;
        }
        if (getParent() == null || other.getParent() == null) {
            return (getResidueIndex() == other.getResidueIndex());
        } else if (getParent() == other.getParent()) {
            return (getResidueIndex() == other.getResidueIndex());
        } else {
            return false;
        }
    }

    /**
     * {@inheritDoc}
     *
     * The Finalize method should be called once all atoms have been added to
     * the Residue. Geometry objects (Bonds, Angles, etc) are then formed,
     * followed by a determination of under-constrained (Dangeling) atoms.
     */
    @Override
    public void finalize(boolean finalizeGeometry, ForceField forceField) {
        setFinalized(false);
        getAtomNode().setName("Atoms (" + getAtomList().size() + ")");
        if (finalizeGeometry) {
            assignBondedTerms(forceField);
            removeLeaves();
        }
        setCenter(getMultiScaleCenter(false));
        setFinalized(true);
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's default
     * O3' coordinates given a North pucker.
     *
     * @return a new double[] with default XYZ coordinates for O3' in a North
     * pucker.
     */
    public double[] getO3sNorth() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = O3sNorthCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's default
     * O3' coordinates given a South pucker.
     *
     * @return a new double[] with default XYZ coordinates for O3' in a South
     * pucker.
     */
    public double[] getO3sSouth() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = O3sSouthCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's original
     * C1' coordinates.
     *
     * @return a new double[] with original XYZ coordinates for C1'.
     */
    public double[] getC1sCoords() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = C1sCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's original
     * O4' coordinates.
     *
     * @return a new double[] with original XYZ coordinates for O4'.
     */
    public double[] getO4sCoords() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = O4sCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's original
     * C4' coordinates.
     *
     * @return a new double[] with original XYZ coordinates for C4'.
     */
    public double[] getC4sCoords() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = C4sCoords[i];
        }
        return ret;
    }

    /**
     * Initializes this (presumably nucleic acid) Residue's C1s, O4s, C4s,
     * O3sNorth, and O3sSouth default coordinates based on default PDB atom
     * locations; to preserve rotamer independence, this must be called before
     * any NA rotamers are applied.
     */
    public void initializeDefaultAtomicCoordinates() {
        if (residueType != ResidueType.NA) {
            return;
        }
        boolean isDeoxy;
        try {
            switch (NucleicAcid3.valueOf(this.getName())) {
                case DAD:
                case DCY:
                case DGU:
                case DTY:
                    isDeoxy = true;
                    break;
                case CYT:
                case ADE:
                case THY:
                case URI:
                case GUA:
                default:
                    isDeoxy = false;
                    break;
            }
            C1sCoords = new double[3];
            ((Atom) getAtomNode("C1\'")).getXYZ(C1sCoords);
            O4sCoords = new double[3];
            ((Atom) getAtomNode("O4\'")).getXYZ(O4sCoords);
            C4sCoords = new double[3];
            ((Atom) getAtomNode("C4\'")).getXYZ(C4sCoords);

            /**
             * With the place flag set false, applySugarPucker returns
             * hypothetical O3' coordinates based on default atom positions and
             * the supplied sugar pucker.
             */
            O3sNorthCoords = RotamerLibrary.applySugarPucker(this, 1, isDeoxy, false);
            O3sSouthCoords = RotamerLibrary.applySugarPucker(this, 2, isDeoxy, false);
        } catch (Exception e) {
            logger.log(Level.WARNING, toString(), e);
        }
    }

    /**
     * Returns the position of this Residue's Alpha Carbon (if it is an amino
     * acid).
     *
     * @return a {@link javax.vecmath.Vector3d} object.
     */
    public Vector3d getAlpha3d() {
        Atom a = (Atom) getAtomNode("CA");
        if (a == null) {
            return null;
        }
        Vector3d temp = new Vector3d();
        a.getV3D(temp);
        return temp;
    }

    /**
     * Returns this Residues Parent Polymer name.
     *
     * @return a {@link java.lang.Character} object.
     */
    public Character getChainID() {
        return chainID;
    }

    // Public data access methods
    /**
     * Returns this Residue's sequence number.
     *
     * @return a int.
     */
    public int getResidueIndex() {
        return resIndex;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = hash(SEED, segID);
        hash = hash(hash, getResidueIndex());
        hash = hash(hash, residueType);
        if (residueType == ResidueType.AA) {
            hash = hash(hash, aa);
        } else if (residueType == ResidueType.NA) {
            hash = hash(hash, na);
        }
        return hash(hash, getName());
    }

    /**
     * {@inheritDoc}
     *
     * Prints "Residue Number: x" to stdout.
     */
    @Override
    public void print() {
        logger.info(" " + toString());
        for (Atom a : getAtomNode().getAtomList()) {
            a.print();
        }

    }

    public double[] getSideChainCOM() {
        Vector3d v = new Vector3d();
        Vector3d v2 = new Vector3d();
        int count = 0;
        for (ListIterator li = getSideChainAtoms().listIterator(); li.hasNext();) {
            Atom a = (Atom) li.next();
            String id = a.getName();
            if (!id.equals("CA") && !id.equals("N") && !id.equals("C")
                    && !id.equals("O")) {
                a.getV3D(v2);
                v.add(v2);
                count++;
            } else if (id.equals("CA")) {
                a.print();
            }
        }
        v.scale(1.0 / count);
        double ret[] = new double[3];
        v.get(ret);
        return ret;
    }

    /**
     * <p>
     * printSideChainCOM</p>
     */
    public void logSideChainCOM() {
        double com[] = this.getSideChainCOM();
        logger.info(String.format(" %s %8.3f %8.3f %8.3f ", getName(), com[0], com[1], com[2]));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
            Material mat) {
// If Color by Residue, pass this Residue's Color
        if (newColorModel == RendererCache.ColorModel.RESIDUE) {
            switch (residueType) {
                case AA:
                    color = AA3Color.get(aa);
                    break;
                case NA:
                    color = NA3Color.get(na);
                    break;
                default:
                    color = null;
            }
            if (color == null) {
                return;
            }
            mat = RendererCache.materialFactory(color);
        } else if (newColorModel == RendererCache.ColorModel.STRUCTURE) {
            color = SSTypeColor.get(ssType);
            mat = RendererCache.materialFactory(color);
        }
        super.setColor(newColorModel, color, mat);
    }

    /**
     * <p>
     * setNumber</p>
     *
     * @param n a int.
     */
    public void setNumber(int n) {
        resIndex = n;
    }

    /**
     * <p>
     * Setter for the field <code>chainID</code>.</p>
     *
     * @param c a {@link java.lang.Character} object.
     */
    public void setChainID(Character c) {
        chainID = c;
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
     * setSSType</p>
     *
     * @param ss a {@link ffx.potential.bonded.Residue.SSType} object.
     */
    public void setSSType(SSType ss) {
        ssType = ss;
    }
    private String shortString = null;

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (shortString == null) {
            shortString = new String("" + resIndex + "-" + getName());
        }
        return shortString;
    }

    private static Point3d point3d = new Point3d();
    private static Point2d point2d = new Point2d();
    /**
     * Constant <code>NA1Set</code>
     */
    public static final EnumSet NA1Set = EnumSet.allOf(NA1.class);
    /**
     * Constant <code>NA3Set</code>
     */
    public static final EnumSet NA3Set = EnumSet.allOf(NA3.class);
    /**
     * Constant <code>NASet</code>
     */
    public static final EnumSet NASet = EnumSet.allOf(NA.class);
    /**
     * Constant <code>NA1toNA3</code>
     */
    public static final HashMap<NA1, NA3> NA1toNA3 = new HashMap<>();
    /**
     * Constant <code>NA3Color</code>
     */
    public static final HashMap<NA3, Color3f> NA3Color = new HashMap<>();
    /**
     * Constant <code>AA1Set</code>
     */
    public static final EnumSet AA1Set = EnumSet.allOf(AA1.class);
    /**
     * Constant <code>AA3Set</code>
     */
    public static final EnumSet AA3Set = EnumSet.allOf(AA3.class);
    /**
     * Constant <code>AASet</code>
     */
    public static final EnumSet AASet = EnumSet.allOf(AA.class);
    /**
     * Constant <code>AA1toAA3</code>
     */
    public static final HashMap<AA1, AA3> AA1toAA3 = new HashMap<>();
    /**
     * Constant <code>AA3toAA1</code>
     */
    public static final HashMap<AA3, AA1> AA3toAA1 = new HashMap<>();
    /**
     * Constant <code>AA3Color</code>
     */
    public static final HashMap<AA3, Color3f> AA3Color = new HashMap<>();
    /**
     * Constant <code>SSTypeColor</code>
     */
    public static final HashMap<SSType, Color3f> SSTypeColor = new HashMap<>();
    /**
     * Constant <code>Ramachandran="new String[17]"</code>
     */
    public static String Ramachandran[] = new String[17];
    
    /**
     * Converts an NA3 enum to an equivalent NA1; if simpleCodes is true, ignores
     * the differences between DNA and RNA (deoxy-cytosine and cytosine are both
     * returned as C, for example).
     * @param na3 To convert
     * @param simpleCodes Whether to use the same codes for DNA and RNA
     * @return NA1 code
     */
    public static NA1 NucleicAcid3toNA1(NucleicAcid3 na3, boolean simpleCodes) {
        if (simpleCodes) {
            switch (na3) {
                case ADE:
                case DAD:
                    return NA1.A;
                case CYT:
                case DCY:
                    return NA1.C;
                case GUA:
                case DGU:
                    return NA1.G;
                case URI:
                    return NA1.U;
                case THY:
                case DTY:
                    return NA1.T;
                default:
                    return NA1.X;
            }
        } else {
            switch (na3) {
                case ADE:
                    return NA1.A;
                case DAD:
                    return NA1.D;
                case CYT:
                    return NA1.C;
                case DCY:
                    return NA1.I;
                case GUA:
                    return NA1.G;
                case DGU:
                    return NA1.B;
                case URI:
                    return NA1.U;
                case THY:
                    return NA1.T;
                case DTY:
                    return NA1.T;
                default:
                    return NA1.X;
            }
        }
    }
    /*
    /**
     * Since enumeration values must start with a letter, an 'M' is added to
     * modified bases whose IUPAC name starts with an integer.
     *
    public enum NucleicAcid3 {

        ADE, GUA, CYT, URI, DAD, DGU, DCY, DTY, THY, MP1, DP2, TP3, UNK, M2MG,
        H2U, M2G, OMC, OMG, PSU, M5MC, M7MG, M5MU, M1MA, YYG
    };
    public enum NA1 {

        A, C, G, U, D, I, B, T, P, Q, R, X;
    }

    public enum NA3 {

        A, C, G, U, DA, DC, DG, DT, MPO, DPO, TPO, UNK;
    }
    */

    static {
        NA1 na1[] = NA1.values();
        NA3 na3[] = NA3.values();
        for (int i = 0; i < NA1.values().length; i++) {
            NA1toNA3.put(na1[i], na3[i]);
        }
    }

    static {
        NA3Color.put(NA3.A, RendererCache.RED);
        NA3Color.put(NA3.C, RendererCache.MAGENTA);
        NA3Color.put(NA3.G, RendererCache.BLUE);
        NA3Color.put(NA3.U, RendererCache.YELLOW);
        NA3Color.put(NA3.DA, RendererCache.RED);
        NA3Color.put(NA3.DC, RendererCache.MAGENTA);
        NA3Color.put(NA3.DG, RendererCache.BLUE);
        NA3Color.put(NA3.DT, RendererCache.ORANGE);
        NA3Color.put(NA3.MPO, RendererCache.GREEN);
        NA3Color.put(NA3.DPO, RendererCache.GREEN);
        NA3Color.put(NA3.TPO, RendererCache.GREEN);
        NA3Color.put(NA3.UNK, RendererCache.CYAN);
    }

    static {
        AA1 aa1[] = AA1.values();
        AA3 aa3[] = AA3.values();
        for (int i = 0; i < AA1.values().length; i++) {
            AA1toAA3.put(aa1[i], aa3[i]);
            AA3toAA1.put(aa3[i], aa1[i]);
        }
    }

    static {
        AA3Color.put(AA3.ALA, RendererCache.GRAY);
        AA3Color.put(AA3.ARG, RendererCache.BLUE);
        AA3Color.put(AA3.ASN, RendererCache.BLUE);
        AA3Color.put(AA3.ASP, RendererCache.RED);
        AA3Color.put(AA3.CYS, RendererCache.YELLOW);
        AA3Color.put(AA3.GLN, RendererCache.BLUE);
        AA3Color.put(AA3.GLU, RendererCache.RED);
        AA3Color.put(AA3.GLY, RendererCache.GRAY);
        AA3Color.put(AA3.ILE, RendererCache.GRAY);
        AA3Color.put(AA3.LEU, RendererCache.GRAY);
        AA3Color.put(AA3.LYS, RendererCache.BLUE);
        AA3Color.put(AA3.MET, RendererCache.YELLOW);
        AA3Color.put(AA3.PHE, RendererCache.GREEN);
        AA3Color.put(AA3.PRO, RendererCache.ORANGE);
        AA3Color.put(AA3.SER, RendererCache.BLUE);
        AA3Color.put(AA3.THR, RendererCache.BLUE);
        AA3Color.put(AA3.TRP, RendererCache.GREEN);
        AA3Color.put(AA3.TYR, RendererCache.GREEN);
        AA3Color.put(AA3.VAL, RendererCache.GRAY);
        AA3Color.put(AA3.HIS, RendererCache.BLUE);
        AA3Color.put(AA3.HIE, RendererCache.BLUE);
        AA3Color.put(AA3.HID, RendererCache.BLUE);
        AA3Color.put(AA3.ORN, RendererCache.ORANGE);
        AA3Color.put(AA3.AIB, RendererCache.ORANGE);
        AA3Color.put(AA3.PCA, RendererCache.ORANGE);
        AA3Color.put(AA3.FOR, RendererCache.RED);
        AA3Color.put(AA3.ACE, RendererCache.RED);
        AA3Color.put(AA3.NH2, RendererCache.BLUE);
        AA3Color.put(AA3.NME, RendererCache.BLUE);
        AA3Color.put(AA3.UNK, RendererCache.MAGENTA);
    }

    static {
        SSTypeColor.put(SSType.NONE, RendererCache.WHITE);
        SSTypeColor.put(SSType.SHEET, RendererCache.PINK);
        SSTypeColor.put(SSType.HELIX, RendererCache.BLUE);
        SSTypeColor.put(SSType.TURN, RendererCache.YELLOW);
    }

    static {
        Ramachandran[0] = "Default (Extended)       [-135.0  135.0]";
        Ramachandran[1] = "Alpha Helix (R)          [ -57.0  -47.0]";
        Ramachandran[2] = "Alpha Helix (L)          [  57.0   47.0]";
        Ramachandran[3] = "3-10 Helix               [ -49.0  -26.0]";
        Ramachandran[4] = "Pi Helix                 [ -57.0  -70.0]";
        Ramachandran[5] = "Polyproline II Helix     [ -79.0  149.0]";
        Ramachandran[6] = "Parallel Beta Strand     [-119.0  113.0]";
        Ramachandran[7] = "Antiparallel Beta Strand [-139.0  135.0]";
        Ramachandran[8] = "Beta-Hairpin 2' (i+1)    [  90.0 -170.0]";
        Ramachandran[9] = "Beta-Hairpin 2' (i+2)    [ -80.0  -10.0]";
        Ramachandran[10] = "Beta-Hairpin 1' (i+1)    [  57.0   47.0]";
        Ramachandran[11] = "Beta-Hairpin 1' (i+2)    [  57.0   47.0]";
        Ramachandran[12] = "Beta-Hairpin 1  (i+1)    [ -57.0  -47.0]";
        Ramachandran[13] = "Beta-Hairpin 1  (i+2)    [ -57.0  -47.0]";
        Ramachandran[14] = "Beta-Hairpin 1  (i+3)    [  90.0 -170.0]";
        Ramachandran[15] = "Beta-Hairpin 3' (i+1)    [  57.0   47.0]";
        Ramachandran[16] = "Beta-Hairpin 3' (i+2)    [ -80.0  -10.0]";
    }

    /**
     * The location of a residue within a chain.
     */
    public enum ResiduePosition {

        FIRST_RESIDUE, MIDDLE_RESIDUE, LAST_RESIDUE
    };

    public enum AA {

        GLYCINE, ALANINE, VALINE, LEUCINE, ILLUECINE, SERINE, THREONINE,
        CYSTIENE, PROLINE, PHENYLALANINE, TYROSINE, TYPTOPHAN, ASPARTATE,
        ASPARTAMINE, GLUTAMATE, GLUTAMINE, METHIONINE, LYSINE, ARGININE,
        HISTIDINE;
    }

    public enum AA1 {

        G, A, V, L, I, S, T, C, P, F, Y, W, D, N, E, Q, M, K, R, H, U, Z, O, B,
        J, f, a, n, m, X;
    }

    public enum AA3 {

        GLY, ALA, VAL, LEU, ILE, SER, THR, CYS, PRO, PHE, TYR, TRP, ASP, ASN,
        GLU, GLN, MET, LYS, ARG, HIS, HID, HIE, ORN, AIB, PCA, FOR, ACE, NH2,
        NME, UNK;
    }

    public enum NA {

        ADENINE, CYTOSINE, GUANINE, URACIL, DEOXYADENINE, DEOXYCYTOSINE,
        DEOXYGUANINE, THYMINE, MONOPHOSPHATE, DIPHOSPHATE, TRIPHOSPHATE;
    }

    public enum NA1 {

        A, C, G, U, D, I, B, T, P, Q, R, X;
    }

    public enum NA3 {

        A, C, G, U, DA, DC, DG, DT, MPO, DPO, TPO, UNK;
    }

    public enum ResidueType {

        NA, AA, UNK;
    }

    public enum SSType {

        NONE, HELIX, SHEET, TURN;
    }
}
