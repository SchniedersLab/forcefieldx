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
     * line.
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
                while (line != null) {
                    line = line.trim();
                    if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
                        try {
                            Integer.parseInt(line.substring(6, 11).trim());
                            Integer.parseInt(line.substring(22, 26).trim());
                            String coordOccTempVals[] = line.substring(30, 66).trim().split(" +");
                            for (String value : coordOccTempVals) {
                                Double.parseDouble(value);
                            }
                            br.close();
                            return true;
                        } catch (NumberFormatException | StringIndexOutOfBoundsException ex) {
                            // Do nothing.
                        }
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
