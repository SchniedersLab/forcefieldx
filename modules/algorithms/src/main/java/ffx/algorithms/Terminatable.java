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

/**
 * Algorithms that can be gracefully terminated early should implement this
 * interface.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Terminatable {

    /**
     * This should be implemented as a blocking interrupt; when the method
     * returns the <code>InterruptableAlgorithm</code> has reached a clean
     * termination point. For example, between minimize or molecular dynamics
     * steps.
     */
    public abstract void terminate();
}
