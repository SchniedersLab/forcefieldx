/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.algorithms;

/**
 * The MetropolisMC interface defines the basic methods of a Metropolis Monte 
 * Carlo application.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 *
 */
public interface MetropolisMC {
    
    /**
     * Returns true if the move from e1 to e2 is accepted. For molecular systems, 
     * this will mean accept if e2 <= e1, or accept with probability of 
     * exp(-dU/kbT) if e2 > e1.
     * 
     * @param e1 Initial energy
     * @param e2 Trial energy
     * @return If move accepted
     */
    public boolean tryMove (double e1, double e2);
    
    /**
     * Sets temperature of Monte Carlo criterion.
     * @param temp 
     */
    public void setTemperature (double temp);
    
    /**
     * Sets whether the implementation prints its own messages.
     * @param print Print energies, accept/reject, etc.
     */
    public void setPrint(boolean print);
    // Might want to set false for non-master nodes in certain algorithms.
    
    /**
     * Return starting energy from last attempted step.
     * @return e1
     */
    public double getE1();
    
    /**
     * Return trial energy from last attempted step.
     * @return e2
     */
    public double getE2();
    
    // Implementations should follow the example and style of DefaultMC
}
