/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import java.util.logging.Logger;

/**
 * The Molecule class is a general container used for simple compounds or in
 * cases where more specialized classes have not been implemented.
 */
public class Molecule extends MSGroup {

    private Logger logger = Logger.getLogger(Molecule.class.getName());
    private static final long serialVersionUID = 1L;
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

    public Molecule() {
    }

    public Molecule(String name) {
        super(name);
        residueName = name;
        chainID = 'A';
    }

    public Molecule(String name, int residueNum, 
            Character chainID, String segID) {
        super(name + "-" + residueNum + " " + segID);
        this.residueName = name;
        this.residueNum = residueNum;
        this.chainID = chainID;
        this.segID = segID;
    }

    @Override
    public void setName(String name) {
        super.setName(name + "-" + residueNum + " " + segID);
        this.residueName = name;
    }

    public String getResidueName() {
        return residueName;
    }

    public int getResidueNumber() {
        return residueNum;
    }

    public Character getChainID() {
        return chainID;
    }

    public String getSegID() {
        return segID;
    }

    public Atom getAtom(String name) {
        for (Atom a : getAtomList()) {
            if (a.getName().equalsIgnoreCase(name)) {
                return a;
            }
        }
        return null;
    }

    /**
     * Allows adding Atom FNodes to the Molecule.
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
                 * Allow overwriting of the root alternate
                 * conformer (' ' or 'A').
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

    @Override
    public void finalize(boolean finalizeGeometry) {
        setFinalized(false);
        getAtomNode().setName("Atoms (" + getAtomList().size() + ")");
        if (finalizeGeometry) {
            //constructValenceTerms();
            collectValenceTerms();
            removeLeaves();
        }
        // findDangelingAtoms();
        setCenter(getMultiScaleCenter(false));
        setFinalized(true);
    }
}
