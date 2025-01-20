// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.utilities;

import java.io.File;
import java.net.URI;
import java.util.logging.Logger;

import static org.apache.commons.io.FileUtils.copyURLToFile;

/**
 * Download utilities.
 *
 * @author Michal J. Schnieders
 * @since 1.0
 */
public class DownloadUtils {

  private static final Logger logger = Logger.getLogger(DownloadUtils.class.getName());

  /**
   * Default constructor.
   */
  private DownloadUtils() {
    // Prevent instantiation.
  }

  /**
   * Download a PDB file.
   *
   * @param pdbID The PDB ID.
   * @return PDB The downloaded PDB filename (or null if there was an error).
   */
  public static String downloadPDB(String pdbID) {
    try {
      String url = "https://files.rcsb.org/download/" + pdbID + ".pdb";
      logger.fine(" Attempting to download PDB from:\n " + url);
      String fileName = pdbID + ".pdb";
      File file = new File(fileName);
      // Time-out after 5 seconds.
      int connectionTimeoutMilliseconds = 5000;
      int readTimeoutMilliseconds = 5000;
      copyURLToFile(new URI(url).toURL(), file, connectionTimeoutMilliseconds, readTimeoutMilliseconds);
      logger.fine(" PDB saved to:\n " + file.getAbsolutePath());
      return fileName;
    } catch (Exception e) {
      logger.warning(" Error downloading PDB ID: " + pdbID + "\n " + e);
      return null;
    }
  }
}
