/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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

package ffx.xray;

import java.io.File;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.ReflectionList;

/**
 *
 * @author Tim Fenn
 */
public interface DiffractionFileFilter {

    /*
     * get reflection information from a reflection file
     * 
     * @param file file to read in
     * @return the reflection list, or null if not enough
     *         information present in the reflection file
     */
    ReflectionList getReflectionList(File file);

    /*
     * get reflection information from a reflection file
     *
     * @param file file to read in
     * @param properties system properties
     * @return the reflection list, or null if not enough
     *         information present in the reflection file
     */
    ReflectionList getReflectionList(File file, CompositeConfiguration properties);

    /*
     * read in reflection file
     *
     * @param file file to read in
     * @param reflectionlist list of reflections to find data indices
     * @param refinementdata data to fill in
     */
    boolean readFile(File file, ReflectionList reflectionlist, RefinementData refinementdata);
}
