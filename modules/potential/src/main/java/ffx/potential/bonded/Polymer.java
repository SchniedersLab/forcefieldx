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

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.logging.Logger;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

import ffx.numerics.VectorMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.BiojavaFilter;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Compound;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.FileConvert;

/**
 * The Polymer class encapsulates a peptide or nucleotide chain.
 *
 * @author Michael J. Schnieders
 *
 */
public class Polymer extends MSGroup implements Chain {

    private static final Logger logger = Logger.getLogger(Polymer.class.getName());
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
    private Structure parentStructure;
    private long hibID; // Hibernate ID
    private String asym_id; // mmCIF chain ID.
    private String swissprotID;
    private List<Molecule> associatedMolecules;
    private Compound compound;

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
     * @param residues Represents a MSNode where the Polymer's residues have
     * been attached.
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
     * A generic method for adding a MSNode to the Polymer.
     */
    @Override
    public MSNode addMSNode(MSNode msNode) {
        /*assert (msNode instanceof Residue || msNode instanceof Molecule);
        getAtomNode().add(msNode);*/
        assert (msNode instanceof Residue);

        Residue residue = (Residue) msNode;
        int resNumber = residue.getResidueIndex();

        MSNode residueNode = getAtomNode();
        int n = residueNode.getChildCount();
        int childIndex = n;

        for (int i=0; i<n; i++) {
            Residue current = (Residue) residueNode.getChildAt(i);
            if (current.getResidueIndex() > resNumber) {
                childIndex = i;
                break;
            }
        }

        getAtomNode().insert(residue, childIndex);
        return msNode;
    }

    /**
     * Joiner joins Moieties m1 and m2 and returns the Geometry objects formed
     * in a Joint.
     *
     * @param residue1 a {@link ffx.potential.bonded.Residue} object.
     * @param residue2 a {@link ffx.potential.bonded.Residue} object.
     * @param forceField the ForceField to use when creating joint bonded terms.
     * @return a {@link ffx.potential.bonded.Joint} object.
     */
    public Joint createJoint(Residue residue1, Residue residue2, ForceField forceField) {
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
     *
     * Finalize should be called after all the Residues have been added to the
     * Polymer. This method in turn calls the Finalize method of each Residue,
     * then forms Joints between adjacent Residues in the Polymer
     */
    @Override
    public void finalize(boolean finalizeGroups, ForceField forceField) {
        ListIterator li;
        List<MSNode> residues = getAtomNodeList();
        setFinalized(false);

        // Finalize the residues in the Polymer
        if (finalizeGroups) {
            for (li = residues.listIterator(); li.hasNext();) {
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
            MSNode joints = getTerms();
            joints.removeAllChildren();
            List<Atom> atoms = getAtomList();

            for (Atom a : atoms) {
                if (a.getNumBonds() > 0) {
                    for (Bond b : a.getFFXBonds()) {
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
            getAtomNode().setName("Sub-Groups " + "(" + residues.size() + ")");
            if (getTerms().getParent() != null) {
                remove(getTerms());
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
    public Character getChainIDChar() {
        return chainID;
    }
    
    @Override
    public void addGroup(Group group) {
        if (!group.has3D()) {
            throw new IllegalArgumentException(" FFX cannot accept any atoms which lack coordinates.");
        }
        if (group instanceof Residue) {
            addMSNode((Residue) group);
        } else {
            List<org.biojava.nbio.structure.Atom> atoms = group.getAtoms();
            List<Residue> newResidues = new ArrayList<>(atoms.size());
            for (org.biojava.nbio.structure.Atom atom : atoms) {
                Atom newAtom = BiojavaFilter.readAtom(atom, this.getChainID());
                String resName = newAtom.getResidueName();
                int resNum = newAtom.getResidueNumber();
                Residue res = this.getResidue(resName, resNum, true);
                res.setChain(this, true);
                res.addAtom(newAtom);
            }
        }
    }
    
    public Polymer clone() {
        Polymer polymer = new Polymer(chainID, getName(), link);
        polymer.setCompound(compound);
        polymer.setId(hibID);
        polymer.setInternalChainID(asym_id);
        polymer.setSwissprotId(swissprotID);
        for (Group group : getAtomGroups()) {
            polymer.addGroup((Group) group.clone());
        }
        return polymer;
    }

    @Override
    public Long getId() {
        return hibID;
    }

    @Override
    public void setId(Long id) {
        hibID = id;
    }

    @Override
    public Group getAtomGroup(int position) {
        return (Group) getResidue(position);
    }

    /**
     * Presently, acts the same as getAtomGroup (all FFX data structures have
     * physical coordinates).
     * @param position Residue to get
     * @return Residue at position
     */
    @Override
    public Group getSeqResGroup(int position) {
        return (Group) getResidue(position);
    }

    @Override
    public List<Group> getAtomGroups() {
        List<Residue> residues = getResidues();
        List<Group> molecules = getAtomLigands();
        List<Group> groupList = new ArrayList<>(residues.size() + molecules.size());
        groupList.addAll(residues);
        groupList.addAll(getAtomLigands());
        return groupList;
    }

    @Override
    public void setAtomGroups(List<Group> groups) {
        for (Group group : getAtomGroups()) {
            if (group instanceof Residue) {
                ((Residue) group).setChain(null, true);
            } else if (group instanceof Molecule) {
                ((Molecule) group).setChain(null, true);
            } else {
                group.setChain(null);
            }
        }
        this.removeAllChildren();
        for (Group group : groups) {
            addGroup(group);
        }
        /*for (Group group : groups) {
            if (group instanceof MSGroup) {
                this.addMSNode((MSGroup) group);
            } else {
                List<org.biojava.nbio.structure.Atom> atomList = group.getAtoms();
                List<Residue> residues = new ArrayList<>();
                for (org.biojava.nbio.structure.Atom atom : atomList) {
                    Atom newAtom = BiojavaFilter.readAtom(atom, this.getName());
                    String resName = newAtom.getResidueName();
                    int resNum = newAtom.getResidueNumber();
                    Residue res = this.getResidue(resName, resNum, true);
                    if (!residues.contains(res)) {
                        residues.add(res);
                    }
                    res.addMSNode(newAtom);
                }
                if (residues.size() != 1) {
                    logger.fine(String.format(" Group %s created nonzero number of Residues %d", group.toString(), residues.size()));
                }
                for (Residue residue : residues) {
                    this.addMSNode(residue);
                }
            }
        }*/
    }

    @Override
    public List<Group> getAtomGroups(GroupType type) {
        List<Group> ret = new ArrayList<>();
        List<Group> groups = getAtomGroups();
        for (Group group : groups) {
            if (group.getType() == type) {
                ret.add(group);
            }
        }
        return ret;
    }

    @Override
    public Group getGroupByPDB(ResidueNumber resNum) throws StructureException {
        int seqNum = resNum.getSeqNum();
        Residue res = getResidue(seqNum);
        if (res == null) {
            throw new StructureException("unknown PDB residue number " + seqNum + " in chain " + chainID);
        }
        return res;
    }

    @Override
    public Group[] getGroupsByPDB(ResidueNumber pdbresnumStart, ResidueNumber pdbresnumEnd) throws StructureException {
        return getGroupsByPDB(pdbresnumStart, pdbresnumEnd, false);
    }

    @Override
    public Group[] getGroupsByPDB(ResidueNumber pdbresnumStart, ResidueNumber pdbresnumEnd, boolean ignoreMissing) throws StructureException {
        int start = pdbresnumStart.getSeqNum();
        int end = pdbresnumEnd.getSeqNum();
        List<Group> groups = new ArrayList<>();
        if (ignoreMissing) {
            for (int i = start; i <= end; i++) {
                Residue res = getResidue(i);
                if (res != null) {
                    groups.add(res);
                }
            }
        } else {
            Residue res = getResidue(start);
            if (res == null) {
                throw new StructureException("did not find start PDB residue number " + pdbresnumStart + " in chain " + chainID);
            }
            groups.add(res);
            res = getResidue(end);
            if (res == null) {
                throw new StructureException("did not find end PDB residue number " + pdbresnumEnd + " in chain " + chainID);
            }
            for (int i = start + 1; i <= end; i++) {
                res = getResidue(i);
                if (res != null) {
                    groups.add(res);
                }
            }
        }
        Group[] ret = new Group[groups.size()];
        groups.toArray(ret);
        return ret;
    }

    @Override
    public int getAtomLength() {
        return getResidues().size();
    }

    /**
     * Identical to getAtomLength(), as FFX currently only supports Groups with
     * physical coordinates (seqres records are ignored).
     * @return Number of residues.
     */
    @Override
    public int getSeqResLength() {
        logger.fine(" Seqres length is equal to atom length for FFX Polymer, as all Groups must have coordinates");
        return getResidues().size();
    }

    @Override
    public void setCompound(Compound compound) {
        this.compound = compound;
    }

    @Override
    public Compound getCompound() {
        return compound;
    }

    /**
     * Only takes the first character (FFX chain IDs are currently single characters).
     * @param name 
     */
    @Override
    public void setChainID(String name) {
        this.chainID = name.charAt(0);
    }

    @Override
    public String getChainID() {
        return "" + chainID;
    }

    @Override
    public String getInternalChainID() {
        return asym_id;
    }

    @Override
    public void setInternalChainID(String internalChainID) {
        asym_id = internalChainID;
    }

    @Override
    public Sequence<?> getBJSequence() {
        String seq = getSeqResSequence();
        ResidueType rtype;
        Residue firstRes = getResidues().get(0);
        try {
            rtype = firstRes.getResidueType();
        } catch (Exception ex) {
            rtype = ResidueType.AA;
        }
        switch (rtype) {
            case AA:
                Sequence<AminoAcidCompound> protSeq = null;

                try {
                    protSeq = new ProteinSequence(seq);
                } catch (CompoundNotFoundException e) {
                    logger.warning(String.format("Could not create sequence object "
                            + "from sequence. Some unknown compound: %s", e.toString()));
                }
                return protSeq;
            case NA:
                Sequence<NucleotideCompound> naSeq = null;
                boolean deoxy = firstRes.isDeoxy();
                if (deoxy) {
                    try {
                        naSeq = new DNASequence(seq);
                    } catch (CompoundNotFoundException e) {
                        logger.warning(String.format("Could not create sequence object "
                                + "from sequence. Some unknown compound: %s", e.toString()));
                    }
                } else {
                    try {
                        naSeq = new RNASequence(seq);
                    } catch (CompoundNotFoundException e) {
                        logger.warning(String.format("Could not create sequence object "
                                + "from sequence. Some unknown compound: %s", e.toString()));
                    }
                }
                return naSeq;
            case UNK:
            default:
                return null;
        }
        

        //TODO: return a DNA sequence if the content is DNA...
        /*Sequence seq;
        List<Residue> residues = getResidues();
        switch (residues.get(0).residueType) {
            case NA:
                StringBuilder sb = new StringBuilder();
                for (Residue res : residues) {
                    try {
                        ResidueEnumerations.NucleicAcid3 rescode = ResidueEnumerations.NucleicAcid3.valueOf(res.getName());
                        sb.append(Residue.NucleicAcid3toNA1(rescode, true).toString());
                    } catch (Exception ex) {
                        logger.fine(String.format(" Exception in getting NA1 value of residue %s: %s", res.toString(), ex.toString()));
                    }
                }
                try {
                    String seqString = sb.toString();
                    if (seqString.contains("U")) {
                        seq = new RNASequence(seqString);
                    } else {
                        seq = new DNASequence(seqString);
                    }
                    return seq;
                } catch (CompoundNotFoundException ex) {
                    logger.warning(String.format("Could not create sequence object "
                            + "from sequence. Some unknown compound: %s", ex.toString()));
                }
                return null;
            case AA:
                sb = new StringBuilder();
                for (Residue res : residues) {
                    try {
                        Residue.AA3 rescode = Residue.AA3.valueOf(res.getName());
                        sb.append(Residue.AA3toAA1.get(rescode));
                    } catch (Exception ex) {
                        logger.fine(String.format(" Exception in getting AA1 value of residue %s: %s", res.toString(), ex.toString()));
                    }
                }
                try {
                    seq = new ProteinSequence(sb.toString());
                    return seq;
                } catch (CompoundNotFoundException ex) {
                    logger.warning(String.format("Could not create sequence object "
                            + "from sequence. Some unknown compound: %s", ex.toString()));
                }
                return null;
            case UNK:
            default:
                return null;
        }*/
    }

    @Override
    public String getAtomSequence() {
        return StructureTools.getAtomSequence(this);
    }

    @Override
    public String getSeqResSequence() {
        return StructureTools.getSeqResSequence(this);
    }

    @Override
    public void setSwissprotId(String sp_id) {
        swissprotID = sp_id;
    }

    @Override
    public String getSwissprotId() {
        return swissprotID;
    }

    /**
     * Acts as a wrapper for getGroups(type); FFX does not support groups without 
     * physical coordinates.
     * @param type
     * @return 
     */
    @Override
    public List<Group> getSeqResGroups(GroupType type) {
        return getAtomGroups(type);
    }

    /**
     * Acts as a wrapper for getGroups(); FFX does not support groups without 
     * physical coordinates.
     * @return 
     */
    @Override
    public List<Group> getSeqResGroups() {
        return getAtomGroups();
    }

    /**
     * Acts as a proxy for setGroups (FFX does not support groups without physical
     * coordinates.
     * @param seqResGroups Groups to set.
     */
    @Override
    public void setSeqResGroups(List<Group> seqResGroups) {
        logger.fine(" FFX is not compatible with seqres Groups; all Groups are assumed to have coordinates.");
        setAtomGroups(seqResGroups);
    }

    @Override
    public void setStructure(Structure parent) {
        if (parentStructure instanceof MolecularAssembly) {
            removeFromParent();
        }
        parentStructure = parent;
        if (parent instanceof MolecularAssembly) {
            ((MolecularAssembly) parent).addMSNode(this);
        }
    }
    
    public void setStructure(Structure parent, boolean onlySetRef) {
        if (onlySetRef) {
            parentStructure = parent;
        } else {
            setStructure(parent);
        }
    }

    @Override
    public Structure getStructure() {
        return parentStructure;
    }

    @Override
    public List<Group> getAtomLigands() {
        if (associatedMolecules == null) {
            if (!(parentStructure instanceof MolecularAssembly)) {
                return null;
            } else {
                associatedMolecules = new ArrayList<>();
                List<Molecule> molecules = ((MolecularAssembly) parentStructure).getMolecules();
                for (Molecule molecule : molecules) {
                    if (molecule.getChainID().equals(chainID)) {
                        associatedMolecules.add(molecule);
                    }
                }
            }
        }
        List<Group> retList = new ArrayList<>();
        retList.addAll(associatedMolecules);
        return retList;
    }

    @Override
    public String toPDB() {
        return FileConvert.toPDB(this);
    }

    @Override
    public String toMMCIF() {
        return FileConvert.toMMCIF(this, true);
    }

    public void addMultiResidue(MultiResidue multiResidue) {
        Residue residue = multiResidue.getActive();
        MSNode residueNode = getAtomNode();
        int index = residueNode.getIndex(residue);
        residueNode.remove(index);
        residueNode.insert(multiResidue, index);
        multiResidue.add(residue);
    }

    /**
     * Get the Phi Psi List for the Polymer
     *
     * @return An ArrayList of Dihedral objects representing the Phi/Psi angles
     * of the Polymer, useful for creating Ramachandran plots
     */
    public List<ArrayList<Torsion>> getPhiPsiList() {
        MSNode dihedrals;
        ListIterator li, lj;
        List<ArrayList<Torsion>> phipsi = new ArrayList<ArrayList<Torsion>>();
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
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements();) {
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
            if (r.getResidueIndex() == resNum) {
                return r;
            }
        }
        // Fall back for non-ordered children
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements();) {
            Residue r = (Residue) e.nextElement();
            if (r.getResidueIndex() == resNum) {
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
     * @param resNum a int.
     * @param create a boolean.
     * @return a {@link ffx.potential.bonded.Residue} object.
     */
    public Residue getResidue(String resName, int resNum, boolean create) {
        if (resNum > 0 && getAtomNode().getChildCount() >= resNum) {
            Residue r = (Residue) getAtomNode().getChildAt(resNum - 1);
            if (r.getResidueIndex() == resNum && r.getName().equalsIgnoreCase(resName)) {
                return r;
            }
        }
        for (Enumeration e = getAtomNode().children(); e.hasMoreElements();) {
            Residue r = (Residue) e.nextElement();
            if (r.getResidueIndex() == resNum && r.getName().equalsIgnoreCase(resName)) {
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
            residue = new Residue(resName, resNum, Residue.ResidueType.UNK,
                    chainID, getName());
        }
        addMSNode(residue);
        residue.setChain(this, true);
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
        for (ListIterator li = getAtomNodeList().listIterator(); li.hasNext();) {
            MSGroup atomGroup = (MSGroup) li.next();
            atomGroup.setColor(newColorModel, color, mat);
        }
        for (Enumeration e = getTerms().children(); e.hasMoreElements();) {
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
