/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.potential;

import java.util.ArrayList;
import java.util.List;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;

/**
 * The AssemblyState class stores the chemical and coordinate state of a Molecular
 * Assembly. Not robust to any chemical perturbation except for mutation of
 * MultiResidues.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class AssemblyState {
    private final MolecularAssembly mola;
    private final Residue[] residues;
    private final ResidueState[] resStates;
    private final Atom[] otherAtoms;
    private final double[][] otherCoords;

    /**
     * Construct a snapshot of a MolecularAssembly. Currently accounts for the
     * coordinates of all entities, and the chemical state of MultiResidues. Does
     * not include velocities, etc.
     *
     * @param assembly To store state of.
     */
    public AssemblyState(MolecularAssembly assembly) {
        mola = assembly;

        List<Residue> residueList = mola.getResidueList();
        List<Residue> copyResList = new ArrayList<>(residueList);
        // Remove all Residue nodes, including subnodes of a MultiResidue, from otherNodeList.
        List<MSNode> otherNodeList = mola.getNodeList();
        otherNodeList.removeAll(residueList);

        // Remove any Residue nodes which are beneath a MultiResidue.
        for (Residue res : copyResList) {
            if (res instanceof MultiResidue) {
                List<Residue> subResidues = ((MultiResidue) res).getConsideredResidues();
                residueList.removeAll(subResidues);
            }
        }

        residues = residueList.toArray(new Residue[residueList.size()]);
        resStates = ResidueState.storeAllCoordinates(residues);

        List<Atom> otherAtomList = new ArrayList<>();
        for (MSNode node : otherNodeList) {
            otherAtomList.addAll(node.getAtomList());
        }
        int nOtherAtoms = otherAtomList.size();
        otherAtoms = new Atom[nOtherAtoms];
        otherAtomList.toArray(otherAtoms);
        otherCoords = ResidueState.storeAtomicCoordinates(otherAtoms);
    }

    /**
     * Revert the state of the associated MolecularAssembly, assuming no chemical
     * changes were made except to MultiResidues.
     */
    public void revertState() {
        int nResidues = residues.length;
        for (int i = 0; i < nResidues; i++) {
            residues[i].revertState(resStates[i]);
        }

        int nOthers = otherAtoms.length;
        for (int i = 0; i < nOthers; i++) {
            otherAtoms[i].setXYZ(otherCoords[i]);
        }
    }
}
