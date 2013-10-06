/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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

import java.io.File;

import javax.swing.filechooser.FileFilter;

import org.apache.commons.io.FilenameUtils;

/**
 * The PDBXFileFilter class is used to choose a PDBML File (*.XML).
 *
 * @author Michael J. Schnieders
 * @since 1.0.
 *
 */
public final class PDBMLFileFilter extends FileFilter {

    /**
     * Public Constructor.
     */
    public PDBMLFileFilter() {
    }

    /**
     * {@inheritDoc}
     *
     * This method return <code>true</code> if the file is a directory or
     * Protein Databank XML File (*.XML).
     */
    @Override
    public boolean accept(File file) {
        if (file.isDirectory()) {
            return true;
        }
        String ext = FilenameUtils.getExtension(file.getName());
        return ext.toUpperCase().startsWith("XML");
    }

    /**
     * {@inheritDoc}
     *
     * Provides a description of the PDBMLFileFilter.
     */
    @Override
    public String getDescription() {
        return new String("Protein Databank XML (*.XML)");
    }
}
