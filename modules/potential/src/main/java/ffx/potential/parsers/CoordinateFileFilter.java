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

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * The CoordinateFileFilter class globally determines if a file is a valid
 * coordinate file (PDB, XYZ, INT, or ARC formats).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public final class CoordinateFileFilter extends FileFilter {

    /**
     * Public Constructor.
     */
    public CoordinateFileFilter() {
    }

    /**
     * {@inheritDoc}
     *
     * This method returns <code>true</code> if the file is a valid coordinate
     * file (xyz, pdb).
     */
    @Override
    public boolean accept(File file) {
        XYZFileFilter xyzFileFilter = new XYZFileFilter();
        if (xyzFileFilter.accept(file)) {
            return true;
        }
        PDBFileFilter pdbFileFilter = new PDBFileFilter();
        if (pdbFileFilter.accept(file)) {
            return true;
        }
        INTFileFilter intFileFilter = new INTFileFilter();
        if (intFileFilter.accept(file)) {
            return true;
        }
        ARCFileFilter arcFileFilter = new ARCFileFilter();
        return arcFileFilter.accept(file);
    }

    /**
     * <p>
     * acceptDeep</p>
     *
     * @param file a {@link java.io.File} object.
     * @return Whether a valid coordinate file.
     */
    public boolean acceptDeep(File file) {
        XYZFileFilter xyzFileFilter = new XYZFileFilter();
        if (xyzFileFilter.acceptDeep(file)) {
            return true;
            //ARC is concat of XYZ files, and so will also be handled.
        }
        PDBFileFilter pdbFileFilter = new PDBFileFilter();
        if (pdbFileFilter.acceptDeep(file)) {
            return true;
        }
        INTFileFilter intFileFilter = new INTFileFilter();
        return intFileFilter.acceptDeep(file);
    }

    /**
     * {@inheritDoc}
     *
     * Provides a description of the CoordinateFileFilter.
     */
    @Override
    public String getDescription() {
        return new String("Protein Databank (*.PDB), TINKER Cartesian Coordinates (*.XYZ),"
                + " TINKER Internal Coordinates (*.INT), or TINKER Archive (*.ARC)");
    }
}
