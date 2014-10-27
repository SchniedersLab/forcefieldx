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
package ffx.xray;

import org.apache.commons.configuration.CompositeConfiguration;

/**
 * <p>
 * RealSpaceRefinementData class.</p>
 *
 * @author Timothy D. Fenn
 *
 */
public class RealSpaceRefinementData {

    protected final int ori[];
    protected final int ext[];
    protected final int ni[];
    protected boolean periodic;
    protected double data[];
    protected double densityscore;

    /**
     * <p>
     * Constructor for RealSpaceRefinementData.</p>
     *
     */
    public RealSpaceRefinementData() {
        ori = new int[3];
        ext = new int[3];
        ni = new int[3];
        periodic = false;
    }

    /**
     * <p>
     * getDataIndex</p>
     *
     * @param x a int.
     * @param y a int.
     * @param z a int.
     * @return a double.
     */
    public double getDataIndex(int x, int y, int z) {
        int index = x + ext[0] * (y + ext[1] * z);
        return data[index];
    }
}
