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
package ffx.algorithms.mc;

import java.util.List;

/**
 * The MetropolisMC interface defines the basic methods of a Metropolis Monte 
 * Carlo application. Intended to be a framework to allow you to take and revert
 * a single step, and return basic information about the last step taken. Many 
 * of methods will only return a meaningful result after the first step.
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
    public boolean evaluateMove (double e1, double e2);
    
    /**
     * Sets temperature of Monte Carlo criterion.
     * @param temp 
     */
    public void setTemperature (double temp);
    
    /**
     * Returns temperature of the Monte Carlo criterion.
     * @return temperature
     */
    public double getTemperature();
    
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
    
    /**
     * Returns the energy as of the last step taken (not including any extra-
     * potential energy adjustments).
     * @return Last step's energy
     */
    public double lastEnergy();
    
    /**
     * If possible, reverts the last successful Monte Carlo step taken.
     */
    public void revertStep();
    
    // The easiest way to implement the next three methods is to follow the 
    // pattern of BoltzmannMC; calculate en1 via currentEnergy() if necessary,
    // and construct a single-member List<MCMove> if necessary, and pass it
    // through to mcStep(List<MCMove> moves, double en1).
    
    /**
     * Calculates the current system energy and performs an MCMove.
     * @param move MCMove to perform
     * @return If move accepted
     */
    public boolean mcStep(MCMove move);
    
    /**
     * Performs an MCMove.
     * @param move MCMove to perform
     * @param en1 Initial energy
     * @return If move accepted
     */
    public boolean mcStep(MCMove move, double en1);
    
    /**
     * Calculates the current system energy and performs a series of moves
     * sequentially as a single hybrid step.
     * @param moves MCMoves to perform
     * @return If move/moves accepted
     */
    public boolean mcStep(List<MCMove> moves);
    
    /**
     * Performs a series of moves sequentially, as a single hybrid step. Should
     * also work with single-member lists.
     *
     * @param moves MCMoves to perform
     * @param en1 Initial energy
     * @return If move/moves accepted.
     */
    public boolean mcStep(List<MCMove> moves, double en1);
    
    /**
     * If last step taken was a success.
     * @return Acceptance of last move
     */
    public boolean getAccept();

}
