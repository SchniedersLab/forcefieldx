/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import ffx.crystal.Crystal;
import org.apache.commons.configuration.CompositeConfiguration;

/**
 *
 * @author Tim Fenn
 */
public interface RealSpaceFileFilter {

    Crystal getCrystal(String filename, CompositeConfiguration properties);

    /*
     * read in Real Space file
     *
     * @param filename file to read in
     * @param refinementdata the {@link RefinementData} object to fill in
     * @param properties system properties
     * @return true if read in properly
     */
    boolean readFile(String filename, RealSpaceRefinementData refinementdata,
            CompositeConfiguration properties);
}
