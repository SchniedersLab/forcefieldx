// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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

import static java.io.File.createTempFile;
import static java.lang.String.format;
import static java.nio.file.Paths.get;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Path;

/**
 * FileUtils class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class FileUtils {

  /**
   * Returns the file name of a temporary copy of <code>input</code> content.
   *
   * @param input The input stream contents are copied to a temporary file.
   * @param prefix Temporary file prefix.
   * @param name Temporary file name.
   * @param suffix Temporary file suffix.
   * @return The temporary file.
   * @throws java.io.IOException An IOException is thrown if a temporary file could not be created.
   */
  public static String copyInputStreamToTmpFile(
      final InputStream input, String prefix, String name, final String suffix) throws IOException {
    File tmpFile = null;

    name = prefix + "." + name + ".";
    try {
      tmpFile = createTempFile(name, "." + suffix);
    } catch (IOException e) {
      System.out.println(format(" Could not extract %s.", name));
      System.err.println(e.toString());
      System.exit(-1);
    }
    tmpFile.deleteOnExit();

    try (input;
        OutputStream output = new BufferedOutputStream(new FileOutputStream(tmpFile))) {
      byte[] buffer = new byte[8192];
      int size;
      while ((size = input.read(buffer)) != -1) {
        output.write(buffer, 0, size);
      }
    }

    return tmpFile.toString();
  }

  /**
   * Constructs a relative path from the present working directory to a file.
   *
   * @param file Construct a relative path to File file.
   * @return Relative path to file.
   */
  public static Path relativePathTo(File file) {
    File pwd = new File(".");
    Path pwdPath = get(pwd.getAbsolutePath());
    Path otherPath = get(file.getAbsolutePath());
    return pwdPath.relativize(otherPath);
  }
}
