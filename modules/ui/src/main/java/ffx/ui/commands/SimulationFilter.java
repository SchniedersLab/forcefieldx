/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
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
package ffx.ui.commands;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import java.io.File;

import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities.FileType;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;

/**
 * The SimulationFilter class parses system data sent by TINKER to Force Field
 * Xplor.
 *
 * @author schnied
 * @version $Id: $
 */
public final class SimulationFilter extends SystemFilter {

    TinkerSystem system;
    Hashtable<Integer, AtomType> atomTypes = new Hashtable<Integer, AtomType>();

    /**
     * <p>Constructor for SimulationFilter.</p>
     *
     * @param sys a {@link ffx.ui.commands.TinkerSystem} object.
     * @param m a {@link ffx.potential.bonded.MolecularAssembly} object.
     */
    public SimulationFilter(TinkerSystem sys, MolecularAssembly m) {
        super(new File(""), m, null, null);
        system = sys;
        fileType = FileType.SIM;
        fileRead = false;
    }

    /** {@inheritDoc} */
    @Override
    public boolean readFile() {
        // Create Molecular Mechanics Data Objects from the TinkerSystem
        // information
        for (int i = 0; i < system.numatoms; i++) {
            AtomType atomType = atomTypes.get(system.types[i]);
            if (atomType == null) {
                atomType = new AtomType(system.types[i], -1, system.name[i],
                        system.story[i], system.atomic[i], system.mass[i], 0);
                atomTypes.put(system.types[i], atomType);
            }
        }
        atomList = new ArrayList<Atom>();
        Vector<Integer> bonds1 = new Vector<Integer>();
        Vector<Integer> bonds2 = new Vector<Integer>();
        double d[] = new double[3];
        int b[] = new int[4];
        for (int i = 0; i < system.numatoms; i++) {
            d[0] = system.coordinates[0][i];
            d[1] = system.coordinates[1][i];
            d[2] = system.coordinates[2][i];
            String s = new String("" + system.types[i]);
            AtomType atomType = atomTypes.get(s);
            Atom a = new Atom(i + 1, new String("" + atomType.type), atomType,
                    d);
            atomList.add(a);
            int b1 = i + 1;
            b[0] = system.connectivity[0][i];
            b[1] = system.connectivity[1][i];
            b[2] = system.connectivity[2][i];
            b[3] = system.connectivity[3][i];
            int j = 0;
            while (j < 4 && b[j] != 0) {
                int b2 = b[j];
                bonds1.add(b1);
                bonds2.add(b2);
                j++;
            }
        }
        bondList = new ArrayList<Bond>();
        for (int i = 0; i < bonds1.size(); i++) {
            int a1 = bonds1.get(i);
            int a2 = bonds2.get(i);
            if (a1 < a2) {
                Atom atom1 = atomList.get(a1 - 1);
                Atom atom2 = atomList.get(a2 - 1);
                bondList.add(new Bond(atom1, atom2));
            }
        }
        setFileRead(true);
        return true;
    }

    /** {@inheritDoc} */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        return false;
    }
}
