//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.xray.parsers;

import java.io.File;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.xray.DiffractionRefinementData;

/**
 * <p>
 * DiffractionFileFilter interface.</p>
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
 */
public interface DiffractionFileFilter {

    /**
     * Get reflection information from a reflection file.
     *
     * @param file File to read in.
     * @return The {@link ffx.crystal.ReflectionList}, or null if not enough information
     * present in the reflection file.
     */
    ReflectionList getReflectionList(File file);

    /**
     * Get reflection information from a reflection file.
     *
     * @param file File to read in.
     * @param properties System properties.
     * @return The {@link ffx.crystal.ReflectionList}, or null if not enough information
     * present in the reflection file.
     */
    ReflectionList getReflectionList(File file, CompositeConfiguration properties);

    /**
     * Read in reflection file.
     *
     * @param file File to read in.
     * @param reflectionList The {@link ffx.crystal.ReflectionList} to find data indices.
     * @param refinementData The {@link DiffractionRefinementData} object to fill in.
     * @param properties System properties.
     * @return True if read in properly.
     */
    boolean readFile(File file, ReflectionList reflectionList,
                     DiffractionRefinementData refinementData, CompositeConfiguration properties);

    /**
     * Attempt to determine resolution of reflection file.
     *
     * @param file File to read in.
     * @param crystal Crystal system to determine resolution information from.
     * @return The resolution.
     */
    double getResolution(File file, Crystal crystal);
}
