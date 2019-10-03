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

import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import static ffx.numerics.math.VectorMath.cross;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;

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
        rotateMultipoles = forceField.getBoolean(ForceField.ForceFieldBoolean.ROTATE_MULTIPOLES, true);
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

        private final double[] trq = new double[3];
        private final double[] u = new double[3];
        private final double[] v = new double[3];
        private final double[] w = new double[3];
        private final double[] r = new double[3];
        private final double[] s = new double[3];
        private final double[] uv = new double[3];
        private final double[] uw = new double[3];
        private final double[] vw = new double[3];
        private final double[] ur = new double[3];
        private final double[] us = new double[3];
        private final double[] vs = new double[3];
        private final double[] ws = new double[3];
        private final double[] t1 = new double[3];
        private final double[] t2 = new double[3];
        private final double[] localOrigin = new double[3];
        private int threadID;

        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void start() {
            threadID = getThreadIndex();
        }

        @Override
        public void run(int lb, int ub) {
            if (gradient) {
                torque.reduce(lb, ub);
                for (int i = lb; i <= ub; i++) {
                    torque(i, torque, grad);
                }
            }
            if (lambdaTerm) {
                lambdaTorque.reduce(lb, ub);
                for (int i = lb; i <= ub; i++) {
                    torque(i, lambdaTorque, lambdaGrad);
                }
            }
        }

        public void torque(int i, AtomicDoubleArray3D tq, AtomicDoubleArray3D gd) {
            final int[] ax = axisAtom[i];
            // Ions, for example, have no torque.
            if (ax == null || ax.length < 2) {
                return;
            }
            final int ia = ax[0];
            final int ib = i;
            final int ic = ax[1];
            int id = 0;

            trq[0] = tq.getX(i);
            trq[1] = tq.getY(i);
            trq[2] = tq.getZ(i);

            double[] x = coordinates[0][0];
            double[] y = coordinates[0][1];
            double[] z = coordinates[0][2];
            localOrigin[0] = x[ib];
            localOrigin[1] = y[ib];
            localOrigin[2] = z[ib];
            u[0] = x[ia];
            u[1] = y[ia];
            u[2] = z[ia];
            v[0] = x[ic];
            v[1] = y[ic];
            v[2] = z[ic];

            // Construct the three rotation axes for the local frame
            diff(u, localOrigin, u);
            diff(v, localOrigin, v);
            switch (frame[i]) {
                default:
                case ZTHENX:
                case BISECTOR:
                    cross(u, v, w);
                    break;
                case THREEFOLD:
                case ZTHENBISECTOR:
                    id = ax[2];
                    w[0] = x[id];
                    w[1] = y[id];
                    w[2] = z[id];
                    diff(w, localOrigin, w);
            }

            double ru = r(u);
            double rv = r(v);
            double rw = r(w);
            scalar(u, 1.0 / ru, u);
            scalar(v, 1.0 / rv, v);
            scalar(w, 1.0 / rw, w);

            // Find the perpendicular and angle for each pair of axes.
            cross(v, u, uv);
            cross(w, u, uw);
            cross(w, v, vw);
            double ruv = r(uv);
            double ruw = r(uw);
            double rvw = r(vw);
            scalar(uv, 1.0 / ruv, uv);
            scalar(uw, 1.0 / ruw, uw);
            scalar(vw, 1.0 / rvw, vw);

            // Compute the sine of the angle between the rotation axes.
            double uvcos = dot(u, v);
            double uvsin = sqrt(1.0 - uvcos * uvcos);
            //double uwcos = dot(u, w);
            //double uwsin = sqrt(1.0 - uwcos * uwcos);
            //double vwcos = dot(v, w);
            //double vwsin = sqrt(1.0 - vwcos * vwcos);
            /*
             * Negative of dot product of torque with unit vectors gives
             * result of infinitesimal rotation along these vectors.
             */
            double dphidu = -(trq[0] * u[0] + trq[1] * u[1] + trq[2] * u[2]);
            double dphidv = -(trq[0] * v[0] + trq[1] * v[1] + trq[2] * v[2]);
            double dphidw = -(trq[0] * w[0] + trq[1] * w[1] + trq[2] * w[2]);
            double[][] g = new double[4][3];
            switch (frame[i]) {
                case ZTHENBISECTOR:
                    // Build some additional axes needed for the Z-then-Bisector method
                    sum(v, w, r);
                    cross(u, r, s);
                    double rr = r(r);
                    double rs = r(s);
                    scalar(r, 1.0 / rr, r);
                    scalar(s, 1.0 / rs, s);
                    // Find the perpendicular and angle for each pair of axes.
                    cross(r, u, ur);
                    cross(s, u, us);
                    cross(s, v, vs);
                    cross(s, w, ws);
                    double rur = r(ur);
                    double rus = r(us);
                    double rvs = r(vs);
                    double rws = r(ws);
                    scalar(ur, 1.0 / rur, ur);
                    scalar(us, 1.0 / rus, us);
                    scalar(vs, 1.0 / rvs, vs);
                    scalar(ws, 1.0 / rws, ws);
                    // Compute the sine of the angle between the rotation axes
                    double urcos = dot(u, r);
                    double ursin = sqrt(1.0 - urcos * urcos);
                    //double uscos = dot(u, s);
                    //double ussin = sqrt(1.0 - uscos * uscos);
                    double vscos = dot(v, s);
                    double vssin = sqrt(1.0 - vscos * vscos);
                    double wscos = dot(w, s);
                    double wssin = sqrt(1.0 - wscos * wscos);
                    // Compute the projection of v and w onto the ru-plane
                    scalar(s, -vscos, t1);
                    scalar(s, -wscos, t2);
                    sum(v, t1, t1);
                    sum(w, t2, t2);
                    double rt1 = r(t1);
                    double rt2 = r(t2);
                    scalar(t1, 1.0 / rt1, t1);
                    scalar(t2, 1.0 / rt2, t2);
                    double ut1cos = dot(u, t1);
                    double ut1sin = sqrt(1.0 - ut1cos * ut1cos);
                    double ut2cos = dot(u, t2);
                    double ut2sin = sqrt(1.0 - ut2cos * ut2cos);
                    double dphidr = -(trq[0] * r[0] + trq[1] * r[1] + trq[2] * r[2]);
                    double dphids = -(trq[0] * s[0] + trq[1] * s[1] + trq[2] * s[2]);
                    for (int j = 0; j < 3; j++) {
                        double du = ur[j] * dphidr / (ru * ursin) + us[j] * dphids / ru;
                        double dv = (vssin * s[j] - vscos * t1[j]) * dphidu / (rv * (ut1sin + ut2sin));
                        double dw = (wssin * s[j] - wscos * t2[j]) * dphidu / (rw * (ut1sin + ut2sin));
                        g[0][j] = du; // Atom A
                        g[1][j] = dv; // Atom C
                        g[2][j] = dw; // Atom D
                        g[3][j] = -(du + dv + dw); // Atom B
                    }
                    gd.add(threadID, ia, g[0][0], g[0][1], g[0][2]);
                    gd.add(threadID, ic, g[1][0], g[1][1], g[1][2]);
                    gd.add(threadID, id, g[2][0], g[2][1], g[2][2]);
                    gd.add(threadID, ib, g[3][0], g[3][1], g[3][2]);
                    break;
                case ZTHENX:
                    for (int j = 0; j < 3; j++) {
                        double du = uv[j] * dphidv / (ru * uvsin) + uw[j] * dphidw / ru;
                        double dv = -uv[j] * dphidu / (rv * uvsin);
                        g[0][j] = du; // Atom A
                        g[1][j] = dv; // Atom C
                        g[2][j] = -(du + dv); // Atom B
                    }
                    gd.add(threadID, ia, g[0][0], g[0][1], g[0][2]);
                    gd.add(threadID, ic, g[1][0], g[1][1], g[1][2]);
                    gd.add(threadID, ib, g[2][0], g[2][1], g[2][2]);
                    break;
                case BISECTOR:
                    for (int j = 0; j < 3; j++) {
                        double du = uv[j] * dphidv / (ru * uvsin) + 0.5 * uw[j] * dphidw / ru;
                        double dv = -uv[j] * dphidu / (rv * uvsin) + 0.5 * vw[j] * dphidw / rv;
                        g[0][j] = du; // Atom A
                        g[1][j] = dv; // Atom C
                        g[2][j] = -(du + dv); // Atom B
                    }
                    gd.add(threadID, ia, g[0][0], g[0][1], g[0][2]);
                    gd.add(threadID, ic, g[1][0], g[1][1], g[1][2]);
                    gd.add(threadID, ib, g[2][0], g[2][1], g[2][2]);
                    break;
                default:
                    String message = "Fatal exception: Unknown frame definition: " + frame[i] + "\n";
                    logger.log(Level.SEVERE, message);
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
