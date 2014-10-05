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
package ffx.utilities;

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

}
