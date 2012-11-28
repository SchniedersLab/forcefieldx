/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;

import javax.swing.filechooser.FileFilter;

/**
 * The MTZFileFilter class is used to choose CCP4 MTZ files
 *
 * @author Michael J. Schnieders
 *
 */
public final class MTZFileFilter extends FileFilter {

    /**
     * Default Constructor
     */
    public MTZFileFilter() {
    }

    /**
     * {@inheritDoc}
     *
     * This method determines whether or not the parm File parameter is a Tinker
     * *.xyz or not, returning true if it is. (Also returns true for any
     * directory)
     */
    @Override
    public boolean accept(File parm) {
        if (parm.isDirectory()) {
            return true;
        }
        String filename = parm.getName().toLowerCase();
        return filename.endsWith(".mtz");
    }

    /**
     * <p>acceptDeep</p>
     *
     * @param file a {@link java.io.File} object.
     * @return a boolean.
     */
    public boolean acceptDeep(File file) {
        try {
            if (file == null || file.isDirectory() || !file.canRead()) {
                return false;
            }
            FileInputStream fis = new FileInputStream(file);
            DataInputStream dis = new DataInputStream(fis);

            byte bytes[] = new byte[80];
            int offset = 0;

            // is it an MTZ file?
            dis.read(bytes, offset, 4);
            String mtzstr = new String(bytes);
            if (!mtzstr.trim().equals("MTZ")) {
                return false;
            }

            return true;
        } catch (Exception e) {
            return true;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Provides a description of this FileFilter
     */
    @Override
    public String getDescription() {
        return new String("CCP4 MTZ Reflection Files: *.mtz");
    }
}
