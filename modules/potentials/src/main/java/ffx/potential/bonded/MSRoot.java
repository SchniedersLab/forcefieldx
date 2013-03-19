/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
package ffx.potential.bonded;

/**
 * The MSRoot class is the root of the Force Field X data structure.
 *
 * @author Michael J. Schnieders
 *
 */
public class MSRoot extends MSNode {

    private static final long serialVersionUID = 1L;
    /**
     * Constant
     * <code>MultiScaleLevel=ROLS.MaxLengthScale</code>
     */
    public static final int MultiScaleLevel = ROLS.MaxLengthScale;

    /**
     * Default MSRoot Constructor
     */
    public MSRoot() {
        super("Structural Heirarchy");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return "Structural Heirarchy";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 7;
        return hash;
    }
}
