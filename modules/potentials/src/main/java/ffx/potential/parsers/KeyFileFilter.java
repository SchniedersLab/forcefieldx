/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
package ffx.potential.parsers;

import java.io.File;

import javax.swing.filechooser.FileFilter;

import org.apache.commons.io.FilenameUtils;

/**
 * The KeyFileFilter class is used to choose a Force Field X keyword (*.KEY) or
 * property (*.properties) file.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public final class KeyFileFilter extends FileFilter {

    /**
     * Default Constructor.
     */
    public KeyFileFilter() {
    }

    /**
     * {@inheritDoc}
     *
     * This method return <code>true</code> if the file is a directory or
     * Force Field X script (*.FFX).
     */
    @Override
    public boolean accept(File file) {
        if (file.isDirectory()) {
            return true;
        }
        String ext = FilenameUtils.getExtension(file.getName()).toUpperCase();

        return ext.startsWith("KEY") || ext.startsWith("PRM") || ext.startsWith("PROPERTIES");
    }

    /**
     * {@inheritDoc}
     *
     * Provides a description of the KeyFileFilter.
     */
    @Override
    public String getDescription() {
        return new String("Force Field X Properties (*.KEY, *.PRM, *.PROPERTIES)");
    }
}
