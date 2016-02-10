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
package ffx.xray;

import java.io.File;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;

/**
 * <p>
 * DiffractionFileFilter interface.</p>
 *
 * @author Timothy D. Fenn
 *
 */
public interface DiffractionFileFilter {

    /*
     * get reflection information from a reflection file
     *
     * @param file file to read in
     * @return the {@link ReflectionList}, or null if not enough
     *         information present in the reflection file
     */
    /**
     * <p>
     * getReflectionList</p>
     *
     * @param file a {@link java.io.File} object.
     * @return a {@link ffx.crystal.ReflectionList} object.
     */
    ReflectionList getReflectionList(File file);

    /*
     * get reflection information from a reflection file
     *
     * @param file file to read in
     * @param properties system properties
     * @return the {@link ReflectionList}, or null if not enough
     *         information present in the reflection file
     */
    /**
     * <p>
     * getReflectionList</p>
     *
     * @param file a {@link java.io.File} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @return a {@link ffx.crystal.ReflectionList} object.
     */
    ReflectionList getReflectionList(File file, CompositeConfiguration properties);

    /*
     * read in reflection file
     *
     * @param file file to read in
     * @param reflectionlist the {@link ReflectionList} to find data indices
     * @param refinementdata the {@link RefinementData} object to fill in
     * @param properties system properties
     * @return true if read in properly
     */
    /**
     * <p>
     * readFile</p>
     *
     * @param file a {@link java.io.File} object.
     * @param reflectionlist a {@link ffx.crystal.ReflectionList} object.
     * @param refinementdata a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @return a boolean.
     */
    boolean readFile(File file, ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata, CompositeConfiguration properties);

    /*
     * attempt to determine resolution of reflection file
     *
     * @param file file to read in
     * @param crystal crystal system to determine resolution information from
     */
    /**
     * <p>
     * getResolution</p>
     *
     * @param file a {@link java.io.File} object.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @return a double.
     */
    double getResolution(File file, Crystal crystal);
}
