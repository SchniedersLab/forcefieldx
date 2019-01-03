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
package ffx.ui.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parsers.SystemFilter;

/**
 * The SimulationFilter class parses system data sent by FFXServer to FFXClient.
 *
 * @author Michael J. Schnieders
 */
public final class SimulationFilter extends SystemFilter {

    SimulationDefinition system;
    Hashtable<Integer, AtomType> atomTypes = new Hashtable<Integer, AtomType>();

    /**
     * <p>
     * Constructor for SimulationFilter.</p>
     *
     * @param sys a {@link ffx.ui.commands.SimulationDefinition} object.
     * @param m a {@link ffx.potential.MolecularAssembly} object.
     */
    public SimulationFilter(SimulationDefinition sys, MolecularAssembly m) {
        super(new File(""), m, null, null);
        system = sys;
        fileType = FileType.SIM;
        fileRead = false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean readFile() {
        // Create Molecular Mechanics Data Objects from the SimulationDefinition
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

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean writeFile(File saveFile, boolean append) {
        return false;
    }

    @Override
    public boolean readNext() {
        return readNext(false);
    }

    @Override
    public boolean readNext(boolean resetPosition) {
        return false;
    }

    @Override
    public void closeReader() {
        //logger.fine(" Reading trajectories not yet supported for MergeFilter");
        // No logger set for SimulationFilter.
    }
}
