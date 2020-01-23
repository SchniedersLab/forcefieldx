//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.mc;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import static ffx.utilities.Constants.R;

import static org.apache.commons.math3.util.FastMath.exp;

/**
 * The BoltzmannMC abstract class is a skeleton for Boltzmann-weighted
 * Metropolis Monte Carlo simulations.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public abstract class BoltzmannMC implements MetropolisMC {

    /**
     * Constant <code>logger</code>
     */
    private static final Logger logger = Logger.getLogger(BoltzmannMC.class.getName());

    private double temperature = 298.15; // Room temperature (also STP).
    private double kbTinv = -1.0 / (R * temperature); // Constant factor for Monte Carlo moves (-1/kbT)
    private boolean print = true;

    private double e1 = 0.0;
    private double e2 = 0.0;
    private double lastE = 0.0;
    protected Random random = new Random();

    private boolean lastAccept = false;

    /**
     * {@inheritDoc}
     * <p>
     * Criterion for accept/reject a move; intended to be used mostly
     * internally.
     */
    @Override
    public boolean evaluateMove(double e1, double e2) {
        if (e2 <= e1) {
            return true;
        } else {
            // p(X) = exp(-U(X)/kb*T)
            double prob = exp(kbTinv * (e2 - e1));

            assert (prob >= 0.0 && prob <= 1.0) : "Probability of a Monte Carlo move up in energy should be 0-1";

            double trial = random.nextDouble();

            return (trial <= prob);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setTemperature(double temp) {
        temperature = temp;
        kbTinv = -1.0 / (R * temperature);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setPrint(boolean print) {
        this.print = print;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getE1() {
        return e1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getE2() {
        return e2;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double lastEnergy() {
        return lastE;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean mcStep(MCMove move) {
        return mcStep(move, currentEnergy());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean mcStep(MCMove move, double en1) {
        List<MCMove> moveList = new ArrayList<>(1);
        moveList.add(move);
        return mcStep(moveList, en1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean mcStep(List<MCMove> moves) {
        return mcStep(moves, currentEnergy());
    }

    /**
     * {@inheritDoc}
     * <p>
     * Performs a Boltzmann-weighted Monte Carlo step with an arbitrary list of
     * moves and a defined starting energy. The list of MCMoves should be of a
     * type with O(1) element access, as the current implementation utilizes an
     * indexed for loop.
     */
    @Override
    public boolean mcStep(List<MCMove> moves, double en1) {
        storeState();
        e1 = en1;

        int nMoves = moves.size();
        for (MCMove move : moves) {
            move.move();
        }

        lastE = currentEnergy(); // Is reset to e1 if move rejected.
        e2 = lastE;
        if (evaluateMove(e1, e2)) {
            lastAccept = true;
            if (print) {
                logger.info(String.format(" Monte Carlo step accepted: e1 -> e2 %10.6f -> %10.6f", e1, e2));
            }
            return true;
        } else {
            lastAccept = false;
            for (int i = nMoves - 1; i >= 0; i--) {
                moves.get(i).revertMove();
            }
            lastE = e1;
            if (print) {
                logger.info(String.format(" Monte Carlo step rejected: e1 -> e2 %10.6f -> %10.6f", e1, e2));
            }
            return false;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTemperature() {
        return temperature;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean getAccept() {
        return lastAccept;
    }

    /**
     * Set the random seed.
     *
     * @param randomseed The seed.
     */
    protected void setRandomSeed(int randomseed) {
        random.setSeed(randomseed);
    }

    /**
     * Must return the current energy of the system.
     *
     * @return Current system energy
     */
    protected abstract double currentEnergy();

    /**
     * Store the state for reverting a move. Must be properly implemented for
     * revertStep() to function properly; otherwise, the implementation of
     * revertStep() should throw an OperationNotSupportedException.
     */
    protected abstract void storeState();
}
