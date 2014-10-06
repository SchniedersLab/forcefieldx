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
package ffx.xray;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;

/**
 * <p>
 * RefinementModel class.</p>
 *
 * @author Timothy D. Fenn
 */
public class RefinementModel {

    private static final Logger logger = Logger.getLogger(RefinementModel.class.getName());
    protected List<Atom> atomList;
    protected final Atom[] atomArray;
    protected List<Integer>[] xIndex;
    protected ArrayList<ArrayList<Residue>> altResidues;
    protected ArrayList<ArrayList<Molecule>> altMolecules;

    /**
     * <p>
     * Constructor for RefinementModel.</p>
     *
     * @param assembly an array of
     * {@link ffx.potential.MolecularAssembly} objects.
     */
    public RefinementModel(MolecularAssembly assembly[]) {
        this(assembly, false);
    }

    /**
     * <p>
     * Constructor for RefinementModel.</p>
     *
     * @param assembly an array of
     * {@link ffx.potential.MolecularAssembly} objects.
     * @param refinemolocc a boolean.
     */
    @SuppressWarnings("unchecked")
    public RefinementModel(MolecularAssembly assembly[], boolean refinemolocc) {
        List<Atom> alist;

        /**
         * Build alternate conformer list for occupancy refinement (if necessary).
         */
        altResidues = new ArrayList<>();
        altMolecules = new ArrayList<>();
        ArrayList<MSNode> nodeList0 = assembly[0].getNodeList();
        ArrayList<Residue> tempResidues = null;
        ArrayList<Molecule> tempMolecules = null;
        boolean altconf;

        /**
         * By residue/molecule.
         */
        for (int i = 0; i < nodeList0.size(); i++) {
            altconf = false;
            MSNode iNode = nodeList0.get(i);

            /**
             * first set up alternate residue restraint list.
             */
            for (Atom a : iNode.getAtomList()) {
                if (!a.getAltLoc().equals(' ')
                        || a.getOccupancy() < 1.0) {
                    if (iNode instanceof Residue) {
                        tempResidues = new ArrayList<>();
                        tempResidues.add((Residue) iNode);
                        altResidues.add(tempResidues);
                        altconf = true;
                        break;
                    } else if (iNode instanceof Molecule) {
                        if (refinemolocc) {
                            tempMolecules = new ArrayList<>();
                            tempMolecules.add((Molecule) iNode);
                            altMolecules.add(tempMolecules);
                        }
                        altconf = true;
                        break;
                    }
                }
            }
            if (altconf) {
                for (int j = 1; j < assembly.length; j++) {
                    ArrayList<MSNode> nlist = assembly[j].getNodeList();
                    MSNode node = nlist.get(i);

                    for (Atom a : node.getAtomList()) {
                        if (!a.getAltLoc().equals(' ')
                                && !a.getAltLoc().equals('A')) {
                            if (node instanceof Residue) {
                                if (tempResidues != null) {
                                    tempResidues.add((Residue) node);
                                }
                                break;
                            } else if (node instanceof Molecule) {
                                if (tempMolecules != null) {
                                    tempMolecules.add((Molecule) node);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }

        /**
         * For mapping between atoms between different molecular assemblies.
         */
        xIndex = new List[assembly.length];
        for (int i = 0; i < assembly.length; i++) {
            xIndex[i] = new ArrayList<>();
        }
        int index = 0;
        // also set up atomList that will be used for SF calc
        atomList = new ArrayList<>();

        // root list
        alist = assembly[0].getAtomList();
        for (Atom a : alist) {
            a.setFormFactorIndex(index);
            xIndex[0].add(index);
            atomList.add(a);
            index++;
        }
        // now add cross references to root and any alternate atoms not in root
        for (int i = 1; i < assembly.length; i++) {
            alist = assembly[i].getAtomList();
            for (Atom a : alist) {
                Atom root = assembly[0].findAtom(a);
                if (root != null
                        && root.getAltLoc().equals(a.getAltLoc())) {
                    xIndex[i].add(root.getFormFactorIndex());
                    a.setFormFactorIndex(root.getFormFactorIndex());
                } else {
                    xIndex[i].add(index);
                    atomList.add(a);
                    index++;
                }
            }
        }
        atomArray = atomList.toArray(new Atom[atomList.size()]);

        for (ArrayList<Residue> list : altResidues) {
            if (list.size() == 1) {
                Residue r = list.get(0);
                logger.log(Level.INFO, "residue: {0}: single conformer, non-unity occupancy: occupancy will be refined independently!", r.toString());
            }
        }

        for (ArrayList<Molecule> list : altMolecules) {
            if (list.size() == 1) {
                Molecule m = list.get(0);
                logger.log(Level.INFO, "molecule: {0}: single conformer, non-unity occupancy: occupancy will be refined independently!", m.toString());
            }
        }
    }
}
