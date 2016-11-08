/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.xray;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
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

    /**
     * Map between atom in different molecular assemblies.
     */
    private List<Integer>[] xIndex;
    private ArrayList<ArrayList<Residue>> altResidues;
    private ArrayList<ArrayList<Molecule>> altMolecules;

    /**
     * An atom list.
     */
    private final List<Atom> totalAtomList;
    /**
     * An atom array.
     */
    private final Atom[] totalAtomArray;

    /**
     * An atom list.
     */
    private final List<Atom> activeAtomList;
    /**
     * An atom array.
     */
    private final Atom[] activeAtomArray;

    /**
     * <p>
     * Constructor for RefinementModel.</p>
     *
     * @param assembly an array of {@link ffx.potential.MolecularAssembly}
     * objects.
     */
    public RefinementModel(MolecularAssembly assembly[]) {
        this(assembly, false);
    }

    /**
     * <p>
     * Constructor for RefinementModel.</p>
     *
     * @param assembly an array of {@link ffx.potential.MolecularAssembly}
     * objects.
     * @param refinemolocc a boolean.
     */
    @SuppressWarnings("unchecked")
    public RefinementModel(MolecularAssembly assembly[], boolean refinemolocc) {
        List<Atom> atomList;

        /**
         * Build alternate conformer list for occupancy refinement (if
         * necessary).
         */
        altResidues = new ArrayList<>();
        altMolecules = new ArrayList<>();
        ArrayList<MSNode> nodeList0 = assembly[0].getNodeList();
        ArrayList<Residue> tempResidues = null;
        ArrayList<Molecule> tempMolecules = null;
        boolean alternateConformer;

        /**
         * By residue/molecule.
         */
        for (int i = 0; i < nodeList0.size(); i++) {
            alternateConformer = false;
            MSNode iNode = nodeList0.get(i);

            /**
             * First set up alternate residue restraint list.
             */
            for (Atom a : iNode.getAtomList()) {
                if (!a.getUse()) {
                    continue;
                }
                if (!a.getAltLoc().equals(' ')
                        || a.getOccupancy() < 1.0) {
                    if (iNode instanceof Residue) {
                        tempResidues = new ArrayList<>();
                        tempResidues.add((Residue) iNode);
                        altResidues.add(tempResidues);
                        alternateConformer = true;
                        break;
                    } else if (iNode instanceof Molecule) {
                        if (refinemolocc) {
                            tempMolecules = new ArrayList<>();
                            tempMolecules.add((Molecule) iNode);
                            altMolecules.add(tempMolecules);
                        }
                        alternateConformer = true;
                        break;
                    }
                }
            }
            if (alternateConformer) {
                for (int j = 1; j < assembly.length; j++) {
                    ArrayList<MSNode> nlist = assembly[j].getNodeList();
                    MSNode node = nlist.get(i);
                    for (Atom a : node.getAtomList()) {
                        if (!a.getUse()) {
                            continue;
                        }
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

        /**
         * Create the atomList to be used for SF calculations.
         */
        totalAtomList = new ArrayList<>();

        /**
         * Create the activeAtomList (i.e. atoms that can move).
         */
        activeAtomList = new ArrayList<>();

        /**
         * Root list.
         */
        atomList = assembly[0].getAtomList();
        int index = 0;
        for (Atom a : atomList) {
            if (!a.getUse()) {
                continue;
            }
            a.setFormFactorIndex(index);
            xIndex[0].add(index);
            index++;
            totalAtomList.add(a);
            if (a.isActive()) {
                activeAtomList.add(a);
            }
        }

        // Now add cross references to root and any alternate atoms not in root
        for (int i = 1; i < assembly.length; i++) {
            atomList = assembly[i].getAtomList();
            for (Atom a : atomList) {
                if (!a.getUse()) {
                    continue;
                }
                Atom root = assembly[0].findAtom(a);
                if (root != null
                        && root.getAltLoc().equals(a.getAltLoc())) {
                    xIndex[i].add(root.getFormFactorIndex());
                    a.setFormFactorIndex(root.getFormFactorIndex());
                } else {
                    xIndex[i].add(index);
                    index++;
                    totalAtomList.add(a);
                    if (a.isActive()) {
                        activeAtomList.add(a);
                    }
                }
            }
        }

        totalAtomArray = totalAtomList.toArray(new Atom[totalAtomList.size()]);
        activeAtomArray = activeAtomList.toArray(new Atom[activeAtomList.size()]);

        for (ArrayList<Residue> list : altResidues) {
            if (list.size() == 1) {
                Residue r = list.get(0);
                logger.log(Level.INFO, " Residue {0} is a single conformer with non-unity occupancy.\n Occupancy will be refined independently!", r.toString());
            }
        }

        for (ArrayList<Molecule> list : altMolecules) {
            if (list.size() == 1) {
                Molecule m = list.get(0);
                logger.log(Level.INFO, " Molecule {0} is a single conformer with non-unity occupancy.\n Occupancy will be refined independently!", m.toString());
            }
        }
    }

    /**
     * @return the xIndex
     */
    public List<Integer>[] getxIndex() {
        return xIndex;
    }

    /**
     * @return the altResidues
     */
    public ArrayList<ArrayList<Residue>> getAltResidues() {
        return altResidues;
    }

    /**
     * @return the altMolecules
     */
    public ArrayList<ArrayList<Molecule>> getAltMolecules() {
        return altMolecules;
    }

    /**
     * @return the totalAtomList
     */
    public List<Atom> getTotalAtomList() {
        return totalAtomList;
    }

    /**
     * @return the totalAtomArray
     */
    public Atom[] getTotalAtomArray() {
        return totalAtomArray;
    }

    /**
     * @return the activeAtomList
     */
    public List<Atom> getActiveAtomList() {
        return activeAtomList;
    }

    /**
     * @return the activeAtomArray
     */
    public Atom[] getActiveAtomArray() {
        return activeAtomArray;
    }
}
