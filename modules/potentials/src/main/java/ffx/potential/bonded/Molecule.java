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
     * Residue name assigned in PDB fles.
     */
    private String residueName = null;
    /**
     * Polymer name assigned in PDB files.
     */
    private String polymerName = null;

    public Molecule() {
    }

    public Molecule(String name) {
        super(name);
        residueName = name;
        polymerName = "A";
    }

    public Molecule(String name, int residueNum, String polymerName) {
        super(name + "-" + residueNum);
        this.residueName = name;
        this.residueNum = residueNum;
        this.polymerName = polymerName;
        if (!polymerName.equalsIgnoreCase(" ") && !polymerName.equalsIgnoreCase("Blank")) {
            this.setName(name + "-" + residueNum + " " + polymerName);
        }
    }

    public String getResidueName() {
        return residueName;
    }

    public int getResidueNumber() {
        return residueNum;
    }

    public String getPolymerName() {
        return polymerName;
    }

    /**
     * Allows adding Atom FNodes to the Molecule.
     */
    @Override
    public MSNode addMSNode(MSNode o) {
        if (o instanceof Atom) {
            MSNode node = getAtomNode().contains(o);
            if (node == null) {
                getAtomNode().add(o);
                setFinalized(false);
            } else {
                return node;
            }
        } else {
            logger.warning("Can not add " + o.getClass() + " to a Molecule, not of type Atom");
        }
        return o;
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
