/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
package ffx.xray;

import java.util.ArrayList;
import java.util.List;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import java.util.logging.Logger;

/**
 *
 * @author Tim Fenn
 */
public class RefinementModel {

    private static final Logger logger = Logger.getLogger(RefinementModel.class.getName());
    protected List<Atom> atomlist;
    protected final Atom atomarray[];
    protected List<Integer> xindex[];
    protected ArrayList<ArrayList<Residue>> altresidues;
    protected ArrayList<ArrayList<Molecule>> altmolecules;

    public RefinementModel(MolecularAssembly assembly[]){
        this(assembly, false);
    }

    public RefinementModel(MolecularAssembly assembly[], boolean refinemolocc) {
        List<Atom> alist;

        // build alternate conformer list for occupancy refinement (if necessary)
        altresidues = new ArrayList<ArrayList<Residue>>();
        altmolecules = new ArrayList<ArrayList<Molecule>>();
        ArrayList<MSNode> nlist0 = assembly[0].getNodeList();
        ArrayList<Residue> rtmp = null;
        ArrayList<Molecule> mtmp = null;
        boolean altconf;

        // by residue/molecule
        for (int i = 0; i < nlist0.size(); i++) {
            altconf = false;
            MSNode node0 = nlist0.get(i);

            // first set up alternate residue restraint list
            for (Atom a : node0.getAtomList()) {
                if (!a.getAltLoc().equals(' ')
                        || a.getOccupancy() < 1.0) {
                    if (node0 instanceof Residue) {
                        rtmp = new ArrayList<Residue>();
                        rtmp.add((Residue) node0);
                        altresidues.add(rtmp);
                        altconf = true;
                        break;
                    } else if (node0 instanceof Molecule) {
                        if (refinemolocc) {
                            mtmp = new ArrayList<Molecule>();
                            mtmp.add((Molecule) node0);
                            altmolecules.add(mtmp);
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
                                if (rtmp != null) {
                                    rtmp.add((Residue) node);
                                }
                                break;
                            } else if (node instanceof Molecule) {
                                if (mtmp != null) {
                                    mtmp.add((Molecule) node);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }

        // for mapping between atoms between different molecular assemblies
        xindex = new List[assembly.length];
        for (int i = 0; i < assembly.length; i++) {
            xindex[i] = new ArrayList<Integer>();
        }
        int index = 0;
        // also set up atomlist that will be used for SF calc
        atomlist = new ArrayList<Atom>();

        // root list
        alist = assembly[0].getAtomList();
        for (Atom a : alist) {
            a.setFormFactorIndex(index);
            xindex[0].add(index);
            atomlist.add(a);
            index++;
        }
        // now add cross references to root and any alternate atoms not in root
        for (int i = 1; i < assembly.length; i++) {
            alist = assembly[i].getAtomList();
            for (Atom a : alist) {
                Atom root = assembly[0].findAtom(a);
                if (root != null
                        && root.getAltLoc().equals(a.getAltLoc())) {
                    xindex[i].add(root.getFormFactorIndex());
                    a.setFormFactorIndex(root.getFormFactorIndex());
                } else {
                    xindex[i].add(index);
                    atomlist.add(a);
                    index++;
                }
            }
        }
        atomarray = atomlist.toArray(new Atom[atomlist.size()]);

        for (ArrayList<Residue> list : altresidues) {
            if (list.size() == 1) {
                Residue r = list.get(0);
                logger.info("residue: " + r.toString() + ": single conformer, non-unity occupancy: occupancy will be refined independently!");
            }
        }

        for (ArrayList<Molecule> list : altmolecules) {
            if (list.size() == 1) {
                Molecule m = list.get(0);
                logger.info("molecule: " + m.toString() + ": single conformer, non-unity occupancy: occupancy will be refined independently!");
            }
        }
    }
}
