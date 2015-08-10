/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.utilities;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

/**
 * Exists primarily to download CIF files from the PDB; may be replaced if the 
 * Biojava API grows to contain a method for downloading CIF files.
 * 
 * @author Michael Schnieders
 * @author Jacob M Litman
 */
public class DownloadUtilities {
    
    private static final Logger logger = Logger.getLogger(DownloadUtilities.class.getName());
    
    /**
     * Copied almost exactly out of MainPanel.
     * @param fromString
     * @return 
     */
    public static File downloadURL(String fromString) {
        /**
         * Check for null input.
         */
        if (fromString == null) {
            return null;
        }

        /**
         * Convert the string to a URL instance.
         */
        URL fromURL = null;
        try {
            fromURL = new URL(fromString);
        } catch (MalformedURLException e) {
            String message = String.format(" URL incorrectly formatted %s.", fromString);
            logger.log(Level.INFO, message, e);
            return null;
        }

        /**
         * Download the URL to a local file.
         */
        logger.info(String.format(" Downloading %s", fromString));
        try {
            File toFile = new File(FilenameUtils.getName(fromURL.getPath()));
            FileUtils.copyURLToFile(fromURL, toFile, 1000, 1000);
            logger.info(String.format(" Saved to %s\n", toFile.getPath()));
            return toFile;
        } catch (IOException ex) {
            logger.log(Level.INFO, " Failed to read URL " + fromURL.getPath(), ex);
            return null;
        }
    }
}
