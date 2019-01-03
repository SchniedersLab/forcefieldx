/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.potential.nonbonded;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.util.Arrays.fill;

import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.util.Pair;

import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.ParticleMeshEwald.LambdaMode;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import static ffx.numerics.math.VectorMath.binomial;

/**
 * Predict Mutual Induced Dipoles based on previous steps.
 *
 * @author Stephen LuCore
 * @since 1.0
 */
public class ScfPredictor {

    private static final Logger logger = Logger.getLogger(ScfPredictor.class.getName());

    /**
     * Number of atoms.
     */
    private int nAtoms;

    /**
     * Maps LambdaMode to array indices: OFF/CONDENSED=0, CONDENSED_NOLIGAND=1,
     * VAPOR=2
     */
    private int mode;
    private final PredictorMode predictorMode;

    public enum PredictorMode {
        NONE, LS, POLY, ASPC
    }

    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    protected double inducedDipole[][][];
    protected double inducedDipoleCR[][][];
    /**
     * Dimensions of [mode][predictorOrder][nAtoms][3]
     */
    private double predictorInducedDipole[][][][];
    private double predictorInducedDipoleCR[][][][];
    /**
     * Induced dipole predictor order.
     */
    private final int predictorOrder;
    /**
     * Induced dipole predictor index.
     */
    private int predictorStartIndex;
    /**
     * Induced dipole predictor count.
     */
    private int predictorCount;
    /**
     * Predicts induced dipoles by locally minimizing and testing eps against
     * squared change in parameters.
     */
    private LeastSquaresPredictor leastSquaresPredictor;
    /**
     * Convergence tolerance of the LeastSquares optimizer.
     */
    private static final double eps = 1.0e-4;

    /**
     * <p>Constructor for ScfPredictor.</p>
     *
     * @param mode  a {@link ffx.potential.nonbonded.ScfPredictor.PredictorMode} object.
     * @param order a int.
     * @param ff    a {@link ffx.potential.parameters.ForceField} object.
     */
    public ScfPredictor(PredictorMode mode, int order, ForceField ff) {
        predictorMode = mode;
        predictorOrder = order;
        predictorCount = 0;
        predictorStartIndex = 0;
        if (predictorMode != PredictorMode.NONE) {
            if (predictorMode == PredictorMode.LS) {
                leastSquaresPredictor = new LeastSquaresPredictor(eps);
            }
            if (ff.getBoolean(ForceFieldBoolean.LAMBDATERM, false) || ExtendedSystem.esvSystemActive) {
                predictorInducedDipole = new double[3][predictorOrder][nAtoms][3];
                predictorInducedDipoleCR = new double[3][predictorOrder][nAtoms][3];
            } else {
                predictorInducedDipole = new double[1][predictorOrder][nAtoms][3];
                predictorInducedDipoleCR = new double[1][predictorOrder][nAtoms][3];
            }
        }
    }

    /**
     * To be called upon initialization and update of inducedDipole arrays in
     * parent.
     *
     * @param inducedDipole   an array of induced dipoles.
     * @param inducedDipoleCR an array of induced dipoles chain rule terms.
     * @param lambdaTerm      a boolean.
     */
    public void setInducedDipoleReferences(double[][][] inducedDipole, double[][][] inducedDipoleCR,
                                           boolean lambdaTerm) {
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        if (lambdaTerm) {
            predictorInducedDipole = new double[3][predictorOrder][nAtoms][3];
            predictorInducedDipoleCR = new double[3][predictorOrder][nAtoms][3];
        } else {
            predictorInducedDipole = new double[1][predictorOrder][nAtoms][3];
            predictorInducedDipoleCR = new double[1][predictorOrder][nAtoms][3];
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return predictorMode.toString();
    }

    /**
     * <p>run.</p>
     *
     * @param lambdaMode a {@link ffx.potential.nonbonded.ParticleMeshEwald.LambdaMode} object.
     */
    public void run(LambdaMode lambdaMode) {
        if (predictorMode == PredictorMode.NONE) {
            return;
        }
        switch (lambdaMode) {
            case OFF:
            case CONDENSED:
                mode = 0;
                break;
            case CONDENSED_NO_LIGAND:
                mode = 1;
                break;
            case VAPOR:
                mode = 2;
                break;
            default:
                mode = 0;
        }
        if (predictorMode != PredictorMode.NONE) {
            switch (predictorMode) {
                case ASPC:
                    aspcPredictor();
                    break;
                case LS:
                    leastSquaresPredictor();
                    break;
                case POLY:
                    polynomialPredictor();
                    break;
                case NONE:
                default:
                    break;
            }
        }
    }

    /**
     * Save the current converged mutual induced dipoles.
     *
     * @param inducedDipole   an array of induced dipoles.
     * @param inducedDipoleCR an array of induced dipoles chain rule terms.
     * @param directDipole    an array of direct dipoles.
     * @param directDipoleCR  an array of direct dipoles chain rule terms.
     */
    public void saveMutualInducedDipoles(
            double[][][] inducedDipole, double[][][] inducedDipoleCR,
            double[][] directDipole, double[][] directDipoleCR) {

        // Current induced dipoles are saved before those from the previous step.
        predictorStartIndex--;
        if (predictorStartIndex < 0) {
            predictorStartIndex = predictorOrder - 1;
        }

        if (predictorCount < predictorOrder) {
            predictorCount++;
        }

        for (int i = 0; i < nAtoms; i++) {
            for (int j = 0; j < 3; j++) {
                predictorInducedDipole[mode][predictorStartIndex][i][j]
                        = inducedDipole[0][i][j] - directDipole[i][j];
                predictorInducedDipoleCR[mode][predictorStartIndex][i][j]
                        = inducedDipoleCR[0][i][j] - directDipoleCR[i][j];
            }
        }
    }

    /**
     * The least-squares predictor with induced dipole information from 8-10
     * previous steps reduces the number SCF iterations by ~50%.
     */
    private void leastSquaresPredictor() {
        if (predictorCount < 2) {
            return;
        }
        try {
            /**
             * The Jacobian and target do not change during the LS optimization,
             * so it's most efficient to update them once before the
             * Least-Squares optimizer starts.
             */
            leastSquaresPredictor.updateJacobianAndTarget();
            int maxEvals, maxIter;
            maxEvals = maxIter = 1000;
            LeastSquaresOptimizer.Optimum optimum = leastSquaresPredictor.predict(maxEvals, maxIter);
            double[] optimalValues = optimum.getPoint().toArray();
            if (logger.isLoggable(Level.FINEST)) {
                logger.finest(String.format("\n LS RMS:            %10.6f", optimum.getRMS()));
                logger.finest(String.format(" LS Iterations:     %10d", optimum.getIterations()));
                logger.finest(String.format(" Jacobian Evals:    %10d", optimum.getEvaluations()));
                logger.finest(String.format(" Root Mean Square:  %10.6f", optimum.getRMS()));
                logger.finest(String.format(" LS Coefficients"));
                for (int i = 0; i < predictorOrder - 1; i++) {
                    logger.finest(String.format(" %2d  %10.6f", i + 1, optimalValues[i]));
                }
            }

            /**
             * Initialize a pointer into predictor induced dipole array.
             */
            int index = predictorStartIndex;
            /**
             * Apply the LS coefficients in order to provide an initial guess at
             * the converged induced dipoles.
             */
            for (int k = 0; k < predictorOrder - 1; k++) {
                /**
                 * Set the current coefficient.
                 */
                double c = optimalValues[k];
                for (int i = 0; i < nAtoms; i++) {
                    for (int j = 0; j < 3; j++) {
                        inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
                        inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
                    }
                }
                index++;
                if (index >= predictorOrder) {
                    index = 0;
                }
            }
        } catch (Exception e) {
            logger.log(Level.WARNING, " Exception computing predictor coefficients", e);

        }
    }

    private class LeastSquaresPredictor {

        double weights[];
        double target[];
        double jacobian[][];
        double initialSolution[];
        double tolerance = 1.0;
        RealVector valuesVector;
        RealVector targetVector;
        RealMatrix jacobianMatrix;
        LeastSquaresOptimizer optimizer;
        ConvergenceChecker<LeastSquaresProblem.Evaluation> checker;

        public LeastSquaresPredictor(double eps) {
            tolerance = eps;
            weights = new double[2 * nAtoms * 3];
            target = new double[2 * nAtoms * 3];
            jacobian = new double[2 * nAtoms * 3][predictorOrder - 1];
            initialSolution = new double[predictorOrder - 1];
            fill(weights, 1.0);
            fill(initialSolution, 0.0);
            initialSolution[0] = 1.0;
            optimizer = new LevenbergMarquardtOptimizer().withParameterRelativeTolerance(eps);
            checker = new EvaluationRmsChecker(eps);
        }

        public void updateJacobianAndTarget() {

            // Update the target.
            int index = 0;
            for (int i = 0; i < nAtoms; i++) {
                target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][0];
                target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][1];
                target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][2];
                target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][0];
                target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][1];
                target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][2];
            }
            targetVector = new ArrayRealVector(target);

            // Update the Jacobian.
            index = predictorStartIndex + 1;
            if (index >= predictorOrder) {
                index = 0;
            }
            for (int j = 0; j < predictorOrder - 1; j++) {
                int ji = 0;
                for (int i = 0; i < nAtoms; i++) {
                    jacobian[ji++][j] = predictorInducedDipole[mode][index][i][0];
                    jacobian[ji++][j] = predictorInducedDipole[mode][index][i][1];
                    jacobian[ji++][j] = predictorInducedDipole[mode][index][i][2];
                    jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][0];
                    jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][1];
                    jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][2];
                }
                index++;
                if (index >= predictorOrder) {
                    index = 0;
                }
            }
            jacobianMatrix = new Array2DRowRealMatrix(jacobian);
        }

        private RealVector value(double[] variables) {

            double[] values = new double[2 * nAtoms * 3];
            for (int i = 0; i < nAtoms; i++) {
                int index = 6 * i;
                values[index] = 0;
                values[index + 1] = 0;
                values[index + 2] = 0;
                values[index + 3] = 0;
                values[index + 4] = 0;
                values[index + 5] = 0;
                int pi = predictorStartIndex + 1;
                if (pi >= predictorOrder) {
                    pi = 0;
                }
                for (int j = 0; j < predictorOrder - 1; j++) {
                    values[index] += variables[j] * predictorInducedDipole[mode][pi][i][0];
                    values[index + 1] += variables[j] * predictorInducedDipole[mode][pi][i][1];
                    values[index + 2] += variables[j] * predictorInducedDipole[mode][pi][i][2];
                    values[index + 3] += variables[j] * predictorInducedDipoleCR[mode][pi][i][0];
                    values[index + 4] += variables[j] * predictorInducedDipoleCR[mode][pi][i][1];
                    values[index + 5] += variables[j] * predictorInducedDipoleCR[mode][pi][i][2];
                    pi++;
                    if (pi >= predictorOrder) {
                        pi = 0;
                    }
                }
            }
            return new ArrayRealVector(values);
        }

        Pair<RealVector, RealMatrix> test = new Pair<>(targetVector, jacobianMatrix);

        public LeastSquaresOptimizer.Optimum predict(int maxEval, int maxIter) {
            RealVector start = new ArrayRealVector(initialSolution);
            LeastSquaresProblem lsp = LeastSquaresFactory.create(
                    function, targetVector, start, checker, maxEval, maxIter);

            LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
            if (true)
                logger.info(String.format(" LS Optimization parameters:"
                                + "  %s %s\n"
                                + "  %s %s\n"
                                + "  %d %d", function, targetVector.toString(),
                        start.toString(), checker.toString(), maxIter, maxEval));
            return optimum;
        }

        MultivariateJacobianFunction function = new MultivariateJacobianFunction() {
            @Override
            public Pair<RealVector, RealMatrix> value(RealVector point) {
                return new Pair<>(targetVector, jacobianMatrix);
            }
        };
    }

    /**
     * Always-stable predictor-corrector for the mutual induced dipoles.
     */
    private void aspcPredictor() {

        if (predictorCount < 6) {
            return;
        }

        final double aspc[] = {22.0 / 7.0, -55.0 / 14.0, 55.0 / 21.0, -22.0 / 21.0, 5.0 / 21.0, -1.0 / 42.0};
        /**
         * Initialize a pointer into predictor induced dipole array.
         */
        int index = predictorStartIndex;
        /**
         * Expansion loop.
         */
        for (int k = 0; k < 6; k++) {
            /**
             * Set the current predictor coefficient.
             */
            double c = aspc[k];
            for (int i = 0; i < nAtoms; i++) {
                for (int j = 0; j < 3; j++) {
                    inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
                    inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
                }
            }
            index++;
            if (index >= predictorOrder) {
                index = 0;

            }
        }
    }

    /**
     * Polynomial predictor for the mutual induced dipoles.
     */
    private void polynomialPredictor() {

        if (predictorCount == 0) {
            return;
        }

        /**
         * Check the number of previous induced dipole vectors available.
         */
        int n = predictorOrder;
        if (predictorCount < predictorOrder) {
            n = predictorCount;
        }
        /**
         * Initialize a pointer into predictor induced dipole array.
         */
        int index = predictorStartIndex;
        /**
         * Initialize the sign of the polynomial expansion.
         */
        double sign = -1.0;
        /**
         * Expansion loop.
         */
        for (int k = 0; k < n; k++) {
            /**
             * Set the current predictor sign and coefficient.
             */
            sign *= -1.0;
            double c = sign * binomial(n, k);
            for (int i = 0; i < nAtoms; i++) {
                for (int j = 0; j < 3; j++) {
                    inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
                    inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
                }
            }
            index++;
            if (index >= predictorOrder) {
                index = 0;
            }
        }
    }
}
