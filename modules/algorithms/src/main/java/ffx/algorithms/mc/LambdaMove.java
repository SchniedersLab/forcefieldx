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

import java.util.Random;
import java.util.logging.Logger;
import static java.lang.Math.abs;

import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering;

/**
 * Define an MC move to update lambda.
 *
 * @author Mallory R. Tollefson
 */
public class LambdaMove implements MCMove {

    private static final Logger logger = Logger.getLogger(LambdaMove.class.getName());

    /**
     * Current value of lambda, which always refreshed from the OST instance.
     */
    private double currentLambda;
    /**
     * Apply the Lambda move to an OST instance.
     */
    private final OrthogonalSpaceTempering orthogonalSpaceTempering;
    /**
     * Random number generator.
     */
    private Random random;
    /**
     * Lambda move size:
     * 1) The standard deviation for continuous moves from a Gaussian distribution.
     * 2) The step size for discrete moves.
     */
    private double moveSize = 0.1;
    /**
     * If true, do continuous moves. Otherwise, use discrete moves.
     */
    private boolean isContinuous = true;

    /**
     * <p>Constructor for LambdaMove.</p>
     *
     * @param orthogonalSpaceTempering a {@link OrthogonalSpaceTempering} object.
     */
    public LambdaMove(OrthogonalSpaceTempering orthogonalSpaceTempering) {
        this.orthogonalSpaceTempering = orthogonalSpaceTempering;
        random = new Random();
    }

    /**
     * <p>Constructor for LambdaMove.</p>
     *
     * @param randomSeed               Random seed to use.
     * @param orthogonalSpaceTempering OrthogonalSpaceTempering instance.
     */
    public LambdaMove(int randomSeed, OrthogonalSpaceTempering orthogonalSpaceTempering) {
        this.orthogonalSpaceTempering = orthogonalSpaceTempering;
        random = new Random(randomSeed);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void move() {
        currentLambda = orthogonalSpaceTempering.getLambda();

        // Draw a trial move from the distribution.
        double dL = isContinuous ? continuousMove() : discreteMove();
        double newLambda = mirror(currentLambda, dL);

        // Update the OST instance.
        orthogonalSpaceTempering.setLambda(newLambda);
    }

    /**
     * Applies 0-1 mirroring conditions to lam + dL. Skips any moves where
     * dL is greater than 1 or less than -1, and skips 50% of moves from
     * 0 or 1 (exact).
     *
     * @param lam Initial lambda.
     * @param dL  Change in lambda.
     * @return    Correctly mirrored lam + dL
     */
    private double mirror(double lam, double dL) {
        // Telescope to public static method because a public static method
        // may be useful in the future.
        return mirror(random, lam, dL);
    }

    /**
     * Applies 0-1 mirroring conditions to lam + dL. Skips any moves where
     * dL is greater than 1 or less than -1, and skips 50% of moves from
     * 0 or 1 (exact).
     *
     * @param random Source of randomness.
     * @param lam    Initial lambda.
     * @param dL     Change in lambda.
     * @return       Correctly mirrored lam + dL
     */
    public static double mirror(Random random, double lam, double dL) {
        if (lam == 0 || lam == 1) {
            boolean skip = random.nextBoolean();
            if (skip) {
                return lam;
            }
        }
        // Eliminate really weird edge cases.
        if (Math.abs(dL) > 1.0) {
            logger.warning(String.format(" Skipping extremely large lambda move of %.3f not between -1 and +1", dL));
            return lam;
        }
        // Math.abs to mirror negative values.
        double newLam = Math.abs(lam + dL);
        // If greater than 1, mirror via 2.0 - val
        return newLam <= 1.0 ? newLam : 2.0 - newLam;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void revertMove() {
        orthogonalSpaceTempering.setLambda(currentLambda);
    }

    /**
     * Get the Lambda move size, which is a standard deviation for continuous moves or step size for discrete moves.
     *
     * @param moveSize a double.
     */
    public void setMoveSize(double moveSize) {
        this.moveSize = moveSize;
    }

    /**
     * Get the Lambda move size, which is a standard deviation for continuous moves or step size for discrete moves.
     *
     * @return The lambda move size.
     */
    public double getMoveSize() {
        return moveSize;
    }

    /**
     * If true, do continuous moves. Otherwise, use discrete moves.
     */
    public boolean isContinuous() {
        return isContinuous;
    }

    /**
     * If true, do continuous moves. Otherwise, use discrete moves.
     */
    public void setContinuous(boolean continuous) {
        isContinuous = continuous;
    }

    private double continuousMove() {
        // Draw a trial move from the distribution.
        return random.nextGaussian() * moveSize;
    }

    private double discreteMove() {
        // Make a discrete move.
        double dL = moveSize;
        if (random.nextBoolean()) {
            dL = -moveSize;
        }
        return dL;
    }
}
