/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * @author Michael Schnieders
 */
public class StringUtils {

    /**
     * <p>
     * padRight</p>
     *
     * @param s a {@link java.lang.String} object.
     * @param n a int.
     * @return a {@link java.lang.String} object.
     */
    public static String padRight(String s, int n) {
        return String.format("%-" + n + "s", s);
    }

    /**
     * <p>
     * padLeft</p>
     *
     * @param s a {@link java.lang.String} object.
     * @param n a int.
     * @return a {@link java.lang.String} object.
     */
    public static String padLeft(String s, int n) {
        return String.format("%" + n + "s", s);
    }

    /**
     * <p>
     * pdbForID</p>
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public static String pdbForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".pdb";
    }

    /**
     * <p>
     * cifForID</p>
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public static String cifForID(String id) {
        if (id.length() != 4) {
            return null;
        }
        return "http://www.rcsb.org/pdb/files/" + id.toLowerCase() + ".cif";
    }

    /**
     * Returns the file name of a temporary copy of <code>input</code> content.
     *
     * @param input a {@link java.io.InputStream} object.
     * @param name a {@link java.lang.String} object.
     * @param suffix a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public static String copyInputStreamToTmpFile(final InputStream input,
            String name, final String suffix) throws IOException {
        File tmpFile = null;
        try {
            name = "ffx." + name + ".";
            tmpFile = File.createTempFile(name, suffix);
        } catch (Exception e) {
            System.out.println(" Could not extract a Force Field X library.");
            System.err.println(e.toString());
            System.exit(-1);
        }

        tmpFile.deleteOnExit();
        OutputStream output = null;
        try {
            output = new BufferedOutputStream(new FileOutputStream(tmpFile));
            byte[] buffer = new byte[8192];
            int size;
            while ((size = input.read(buffer)) != -1) {
                output.write(buffer, 0, size);
            }
        } finally {
            if (input != null) {
                input.close();
            }
            if (output != null) {
                output.close();
            }
        }

        return tmpFile.toString();
    }

}
