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
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import javax.swing.filechooser.FileFilter;

import org.apache.commons.io.FilenameUtils;

/**
 * The INTFileFilter class is used to choose TINKER Internal Coordinate (*.INT)
 * files.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class INTFileFilter extends FileFilter {

    /**
     * Default Constructor.
     */
    public INTFileFilter() {
    }

    /**
     * {@inheritDoc}
     *
     * This method return
     * <code>true</code> if the file is a directory or matches the extension for
     * TINKER internal coordinates (*.INT).
     */
    @Override
    public boolean accept(File file) {
        if (file.isDirectory()) {
            return true;
        }
        String ext = FilenameUtils.getExtension(file.getName());
        return ext.toUpperCase().startsWith("INT");
    }

    /**
     * <p>acceptDeep</p>
     *
     * @param parm a {@link java.io.File} object.
     * @return a boolean.
     */
    public boolean acceptDeep(File parm) {
        try {
            if (parm == null || parm.isDirectory() || !parm.canRead()) {
                return false;
            }
            FileReader fr = new FileReader(parm);
            BufferedReader br = new BufferedReader(fr);
            if (!br.ready()) {
                // Empty File?
                return false;
            }
            // If the first token is not an integer this file is not
            // an Internal Coordinates File.
            String rawdata = br.readLine();
            String header[] = rawdata.trim().split(" +");
            if (header == null || header.length == 0) {
                return false;
            }
            try {
                Integer.parseInt(header[0]);
            } catch (Exception e) {
                return false;
            }
            // If the the first Atom line does not begin with an integer and
            // contain
            // three tokens, it is not an internal coordinate file.
            String firstAtom = br.readLine();
            if (firstAtom == null) {
                return false;
            }
            br.close();
            fr.close();
            String data[] = firstAtom.trim().split(" +");
            if (data == null || data.length != 3) {
                return false;
            }
            try {
                Integer.parseInt(data[0]);
            } catch (Exception e) {
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
     * Provides a description of the INTFileFilter.
     */
    @Override
    public String getDescription() {
        return new String("TINKER Internal Coordinates (*.INT)");
    }
}
