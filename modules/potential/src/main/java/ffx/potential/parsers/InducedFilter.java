/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.logging.Logger;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;

/**
 * The InducedFilter class parses TINKER Induced Dipole (*.*U) files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class InducedFilter {

    private static final Logger logger = Logger.getLogger(InducedFilter.class.getName());
    MolecularAssembly fsystem;
    File file;

    /**
     * <p>
     * Constructor for InducedFilter.</p>
     *
     * @param s a {@link ffx.potential.MolecularAssembly} object.
     * @param f a {@link java.io.File} object.
     */
    public InducedFilter(MolecularAssembly s, File f) {
        fsystem = s;
        file = f;
    }

    /**
     * <p>
     * read</p>
     *
     * @return a boolean.
     */
    public boolean read() {
        if (!file.exists() || !file.canRead()) {
            return false;
        }
        try {
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String data = br.readLine();
            String tokens[] = data.trim().split(" +");
            if (tokens.length == 0) {
                return false;
            }
            int numatoms = Integer.parseInt(tokens[0]);
            if (numatoms != fsystem.getAtomList().size()) {
                return false;
            }
            // Read the Induced Dipoles
            double x[][] = new double[numatoms][3];
            for (int i = 0; i < numatoms; i++) {
                data = br.readLine().trim();
                tokens = data.split(" +");
                if (tokens.length != 5) {
                    return false;
                }
                x[i][0] = Double.parseDouble(tokens[2]);
                x[i][1] = Double.parseDouble(tokens[3]);
                x[i][2] = Double.parseDouble(tokens[4]);
            }
            List<Atom> atoms = fsystem.getAtomList();
            double max = 0.0d;
            for (Atom a : atoms) {
                int j = a.getXYZIndex() - 1;
                // a.setInducedDipole(-x[j][0], -x[j][1], -x[j][2]);
            }
            logger.warning("Max Induced: " + max);
            br.close();
            fr.close();
        } catch (Exception e) {
            return false;
        }
        return true;
    }
}
