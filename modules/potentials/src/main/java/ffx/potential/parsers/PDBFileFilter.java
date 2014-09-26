/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
import java.io.IOException;

import javax.swing.filechooser.FileFilter;

import org.apache.commons.io.FilenameUtils;

/**
 * The PDBFileFilter class is used to choose a Protein Databank (*.PDB) file.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class PDBFileFilter extends FileFilter {

    /**
     * Public Constructor.
     */
    public PDBFileFilter() {
    }

    /**
     * {@inheritDoc}
     *
     * This method return <code>true</code> if the file is a directory or
     * Protein Databank File (*.PDB).
     */
    @Override
    public boolean accept(File file) {
        if (file.isDirectory()) {
            return true;
        }
        String ext = FilenameUtils.getExtension(file.getName());
        return ext.toUpperCase().startsWith("PDB");
    }

    /**
     * <p>
     * acceptDeep</p> Accepts a PDB file if it finds at least one parseable ATOM
     * line and one parseable TER line.
     *
     * @param file a {@link java.io.File} object.
     * @return Whether a valid PDB file.
     */
    public boolean acceptDeep(File file) {
        try {
            if (file == null || file.isDirectory() || !file.canRead()) {
                return false;
            }
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                if (!br.ready()) {
                    return false;
                }
                String line = br.readLine();
                if (line != null) {
                    line = line.trim();
                } else {
                    return false;
                }
                boolean validAtomLine = false;
                boolean validTerLine = false;
                while (line != null) {
                    line = line.trim();
                    if (!validAtomLine && line.startsWith("ATOM  ")) {
                        try {
                            Integer.parseInt(line.substring(6, 11).trim());
                            Integer.parseInt(line.substring(22, 26).trim());
                            String coordOccTempVals[] = line.substring(30, 66).trim().split(" +");
                            for (String value : coordOccTempVals) {
                                Double.parseDouble(value);
                            }
                            validAtomLine = true;
                        } catch (NumberFormatException | StringIndexOutOfBoundsException ex) {
                            // Do nothing.
                        }
                    } else if (line.startsWith("TER")) {
                        try {
                            /* In a perfect world, every PDB file which claims to be at the 3.3 standard
                             * will actually be at the 3.3 standard.
                             Integer.parseInt(line.substring(6, 11).trim());
                             Integer.parseInt(line.substring(22, 26).trim());*/
                            validTerLine = true;
                        } catch (NumberFormatException | StringIndexOutOfBoundsException ex) {
                            // Do nothing.
                        }
                    }
                    if (validAtomLine && validTerLine) {
                        return true;
                    }
                    line = br.readLine();
                }
            }
        } catch (IOException e) {
            return false;
        }
        return false;
    }

    /**
     * {@inheritDoc}
     *
     * Provides a description of the PDBFileFilter.
     */
    @Override
    public String getDescription() {
        return new String("Protein Databank (*.PDB)");
    }
}
