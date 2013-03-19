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
package ffx.potential.nonbonded;

/**
 * By implementing the MaskingInterface interface, interaction pairs can be
 * excluded during {@link NeighborList} construction.
 *
 * During nonbonded calculations, the interactions between bonded atoms (1-2,
 * 1-3 or 1-4, for example) may be masked to zero. The interaction energy
 * between excluded atoms is described by the bonded energy terms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public interface MaskingInterface {

    /**
     * Interactions with atom i that should not be included in the NeighborList
     * should be set to 0.
     *
     * @param mask The masking array.
     * @param i The atom whose masking rules should be applied.
     * @since 1.0
     */
    public void applyMask(final byte mask[], final int i);

    /**
     * After calling removeMask, all entries in the mask array should be 1.
     *
     * @param mask The masking array.
     * @param i The atom whose masking rules should be removed.
     * @since 1.0
     */
    public void removeMask(final byte mask[], final int i);
}
