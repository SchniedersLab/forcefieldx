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
package ffx.potential.nonbonded;

/**
 * By implementing the MaskingInterface interface, interaction pairs can be
 * excluded during {@link ffx.potential.nonbonded.NeighborList} construction.
 * <p>
 * During non-bonded calculations, the interactions between bonded atoms (1-2,
 * 1-3 or 1-4, for example) may be masked to zero. The interaction energy
 * between excluded atoms is described by the bonded energy terms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface MaskingInterface {

    /**
     * Interactions with atom i that should not be included in the NeighborList
     * should be set to 0.
     *
     * @param i     The atom whose masking rules should be applied.
     * @param is14  True if atom i and the current atom are 1-4 to each other.
     * @param masks One or more masking arrays.
     * @since 1.0
     */
    void applyMask(final int i, final boolean[] is14, final double[]... masks);

    /**
     * After calling removeMask, all entries in the mask array should be 1 and is14 array false.
     *
     * @param i     The atom whose masking rules should be removed.
     * @param is14  True if atom i and the current atom are 1-4 to each other.
     * @param masks One or more masking arrays.
     * @since 1.0
     */
    void removeMask(final int i, final boolean[] is14, final double[]... masks);
}
