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
package ffx.algorithms;

import ffx.potential.bonded.MolecularAssembly;

/**
 * The AlgorithmListener will be notified at regular intervals during an
 * algorithm. This interface is useful for updating the user interface or
 * terminating the algorithm.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface AlgorithmListener {

    /**
     * After a successful step or interval of an algorithm,
     * this method of the listener will be called.
     *
     * @param active The system the algorithm is operating on.
     * @return A return of <code>true</code> indicates the algorithm continues.
     */
    public abstract boolean algorithmUpdate(MolecularAssembly active);
}
