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
package ffx.potential.nonbonded.pme;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;

/**
 * Parallel conversion of torques into forces, and then reduce them.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ReduceRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(ReduceRegion.class.getName());

    /**
     * If set to false, multipoles are fixed in their local frame and torques
     * are zero, which is useful for narrowing down discrepancies between
     * analytic and finite-difference derivatives(default is true).
     */
    private final boolean rotateMultipoles;
    /**
     * If lambdaTerm is true, some ligand atom interactions with the environment
     * are being turned on/off.
     */
    private boolean lambdaTerm;
    /**
     * If true, compute coordinate gradient.
     */
    private boolean gradient;
    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    private double[][][] coordinates;
    /**
     * Multipole frame definition.
     */
    private MultipoleFrameDefinition[] frame;
    /**
     * Multipole frame defining atoms.
     */
    private int[][] axisAtom;
    /**
     * Atomic Gradient array.
     */
    private AtomicDoubleArray3D grad;
    /**
     * Atomic Torque array.
     */
    private AtomicDoubleArray3D torque;
    /**
     * Partial derivative of the gradient with respect to Lambda.
     */
    private AtomicDoubleArray3D lambdaGrad;
    /**
     * Partial derivative of the torque with respect to Lambda.
     */
    private AtomicDoubleArray3D lambdaTorque;

    private final TorqueLoop[] torqueLoop;
    private final ReduceLoop[] reduceLoop;

    public ReduceRegion(int threadCount, ForceField forceField) {
        torqueLoop = new TorqueLoop[threadCount];
        reduceLoop = new ReduceLoop[threadCount];
        rotateMultipoles = forceField.getBoolean("ROTATE_MULTIPOLES", true);
    }

    public void init(boolean lambdaTerm, boolean gradient,
                     Atom[] atoms, double[][][] coordinates,
                     MultipoleFrameDefinition[] frame, int[][] axisAtom,
                     AtomicDoubleArray3D grad, AtomicDoubleArray3D torque,
                     AtomicDoubleArray3D lambdaGrad, AtomicDoubleArray3D lambdaTorque) {
        this.lambdaTerm = lambdaTerm;
        this.gradient = gradient;
        this.atoms = atoms;
        this.coordinates = coordinates;
        this.frame = frame;
        this.axisAtom = axisAtom;
        this.grad = grad;
        this.torque = torque;
        this.lambdaGrad = lambdaGrad;
        this.lambdaTorque = lambdaTorque;
    }

    /**
     * Execute the ReduceRegion with the passed ParallelTeam.
     *
     * @param parallelTeam The ParallelTeam instance to execute with.
     */
    public void excuteWith(ParallelTeam parallelTeam) {
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = "Exception calculating torques.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    @Override
    public void run() {
        int nAtoms = atoms.length;
        try {
            int threadIndex = getThreadIndex();
            if (torqueLoop[threadIndex] == null) {
                torqueLoop[threadIndex] = new TorqueLoop();
                reduceLoop[threadIndex] = new ReduceLoop();
            }
            if (rotateMultipoles) {
                execute(0, nAtoms - 1, torqueLoop[threadIndex]);
            }
            execute(0, nAtoms - 1, reduceLoop[threadIndex]);
        } catch (Exception e) {
            String message = "Fatal exception computing torque in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private class TorqueLoop extends IntegerForLoop {

        private Torque torques;
        private int threadID;

        TorqueLoop() {
            torques = new Torque();
        }

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void start() {
            threadID = getThreadIndex();
            torques.init(axisAtom, frame, coordinates);
        }

        @Override
        public void run(int lb, int ub) {
            if (gradient) {
                torque.reduce(lb, ub);
                double[] trq = new double[3];
                for (int i = lb; i <= ub; i++) {
                    double[][] g = new double[4][3];
                    int[] frameIndex = {-1,-1,-1,-1};
                    trq[0] = torque.getX(i);
                    trq[1] = torque.getY(i);
                    trq[2] = torque.getZ(i);
                    torques.torque(i, 0, trq, frameIndex, g);
                    for (int j = 0; j<4; j++) {
                        int index = frameIndex[j];
                        if (index >= 0) {
                            double[] gj = g[j];
                            grad.add(threadID, index, gj[0], gj[1], gj[2]);
                        }
                    }
                }
            }
            if (lambdaTerm) {
                lambdaTorque.reduce(lb, ub);
                double[] trq = new double[3];
                for (int i = lb; i <= ub; i++) {
                    double[][] g = new double[4][3];
                    int[] frameIndex = {-1,-1,-1,-1};
                    trq[0] = lambdaTorque.getX(i);
                    trq[1] = lambdaTorque.getY(i);
                    trq[2] = lambdaTorque.getZ(i);
                    torques.torque(i, 0, trq, frameIndex, g);
                    for (int j = 0; j<4; j++) {
                        int index = frameIndex[j];
                        if (index >= 0) {
                            double[] gj = g[j];
                            lambdaGrad.add(threadID, index, gj[0], gj[1], gj[2]);
                        }
                    }
                }
            }
        }
    }

    private class ReduceLoop extends IntegerForLoop {

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            if (gradient) {
                grad.reduce(lb, ub);
                for (int i = lb; i <= ub; i++) {
                    Atom ai = atoms[i];
                    ai.addToXYZGradient(grad.getX(i), grad.getY(i), grad.getZ(i));
                }
            }
            if (lambdaTerm) {
                lambdaGrad.reduce(lb, ub);
            }
        }
    }
}
