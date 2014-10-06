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

import java.util.logging.Logger;

import ffx.potential.parameters.ForceField;

/**
 * The Molecule class is a general container used for simple compounds or in
 * cases where more specialized classes have not been implemented.
 *
 * @author Michael J. Schnieders
 *
 */
public class Molecule extends MSGroup {

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
    public int getResidueNumber() {
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
    public Atom getAtom(String name) {
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
}
