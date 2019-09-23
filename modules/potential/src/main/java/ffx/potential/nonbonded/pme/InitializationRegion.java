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

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.max;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.norm;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;
import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t200;

/**
 * Parallel initialization of accumulation arrays, expand atomic coordinates
 * and rotation of multipoles into the global frame.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class InitializationRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(InitializationRegion.class.getName());

    /**
     * If set to false, multipoles are fixed in their local frame and torques
     * are zero, which is useful for narrowing down discrepancies between
     * analytic and finite-difference derivatives(default is true).
     */
    private final boolean rotateMultipoles;
    /**
     * If set to false, multipole charges are set to zero (default is true).
     */
    private final boolean useCharges;
    /**
     * If set to false, multipole dipoles are set to zero (default is true).
     */
    private final boolean useDipoles;
    /**
     * If set to false, multipole quadrupoles are set to zero (default is true).
     */
    private final boolean useQuadrupoles;
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
     * Scale multipole moments by a lambda scale factor.
     */
    private double lambdaScaleMultipoles;
    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    private double[][][] coordinates;
    /**
     * Unit cell and spacegroup information.
     */
    private Crystal crystal;
    /**
     * Multipole frame definition.
     */
    private MultipoleFrameDefinition[] frame;
    /**
     * Multipole frame defining atoms.
     */
    private int[][] axisAtom;
    /**
     * Permanent multipoles in their local frame.
     */
    private double[][] localMultipole;
    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    private double[][][] globalMultipole;
    private double[] polarizability;
    /**
     * When computing the polarization energy at Lambda there are 3 pieces.
     * <p>
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand).
     * <p>
     * 2.) Uenv = The polarization energy of the system without the ligand.
     * <p>
     * 3.) Uligand = The polarization energy of the ligand by itself.
     * <p>
     * Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     * <p>
     * Set the "use" array to true for all atoms for part 1. Set the "use" array
     * to true for all atoms except the ligand for part 2. Set the "use" array
     * to true only for the ligand atoms for part 3.
     * <p>
     * The "use" array can also be employed to turn off atoms for computing the
     * electrostatic energy of sub-structures.
     */
    private boolean[] use;
    /**
     * Neighbor lists, including atoms beyond the real space cutoff.
     * [nsymm][nAtoms][nAllNeighbors]
     */
    private int[][][] neighborLists;
    /**
     * Neighbor lists, without atoms beyond the real space cutoff.
     * [nSymm][nAtoms][nIncludedNeighbors]
     */
    private int[][][] realSpaceLists;
    private int[][][] vaporLists;
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
    private final InitializationLoop[] initializationLoop;
    private final RotateMultipolesLoop[] rotateMultipolesLoop;

    public InitializationRegion(int maxThreads, ForceField forceField) {
        initializationLoop = new InitializationLoop[maxThreads];
        rotateMultipolesLoop = new RotateMultipolesLoop[maxThreads];
        useCharges = forceField.getBoolean(ForceField.ForceFieldBoolean.USE_CHARGES, true);
        useDipoles = forceField.getBoolean(ForceField.ForceFieldBoolean.USE_DIPOLES, true);
        useQuadrupoles = forceField.getBoolean(ForceField.ForceFieldBoolean.USE_QUADRUPOLES, true);
        rotateMultipoles = forceField.getBoolean(ForceField.ForceFieldBoolean.ROTATE_MULTIPOLES, true);
    }

    public void init(boolean lambdaTerm, boolean gradient, double lambdaScaleMultipoles,
                     Atom[] atoms, double[][][] coordinates, Crystal crystal,
                     MultipoleFrameDefinition[] frame, int[][] axisAtom,
                     double[][] localMultipole, double[][][] globalMultipole, double[] polarizability,
                     boolean[] use, int[][][] neighborLists, int[][][] realSpaceLists, int[][][] vaporLists,
                     AtomicDoubleArray3D grad, AtomicDoubleArray3D torque,
                     AtomicDoubleArray3D lambdaGrad, AtomicDoubleArray3D lambdaTorque) {
        this.lambdaTerm = lambdaTerm;
        this.gradient = gradient;
        this.lambdaScaleMultipoles = lambdaScaleMultipoles;
        this.atoms = atoms;
        this.coordinates = coordinates;
        this.crystal = crystal;
        this.frame = frame;
        this.axisAtom = axisAtom;
        this.localMultipole = localMultipole;
        this.globalMultipole = globalMultipole;
        this.polarizability = polarizability;
        this.use = use;
        this.neighborLists = neighborLists;
        this.realSpaceLists = realSpaceLists;
        this.vaporLists = vaporLists;
        this.grad = grad;
        this.torque = torque;
        this.lambdaGrad = lambdaGrad;
        this.lambdaTorque = lambdaTorque;
    }

    /**
     * Execute the InitializationRegion with the passed ParallelTeam.
     *
     * @param parallelTeam The ParallelTeam instance to execute with.
     */
    public void executeWith(ParallelTeam parallelTeam) {
        try {
            parallelTeam.execute(this);
        } catch (RuntimeException e) {
            String message = "RuntimeException expanding coordinates and rotating multipoles.\n";
            logger.log(Level.WARNING, message, e);
            throw e;
        } catch (Exception e) {
            String message = "Fatal exception expanding coordinates and rotating multipoles.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    @Override
    public void run() {

        int nAtoms = atoms.length;
        int threadIndex = getThreadIndex();
        if (initializationLoop[threadIndex] == null) {
            initializationLoop[threadIndex] = new InitializationLoop();
            rotateMultipolesLoop[threadIndex] = new RotateMultipolesLoop();
        }
        try {
            execute(0, nAtoms - 1, initializationLoop[threadIndex]);
            execute(0, nAtoms - 1, rotateMultipolesLoop[threadIndex]);
        } catch (Exception e) {
            String message = "Fatal exception initializing coordinates in thread: " + threadIndex + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private class InitializationLoop extends IntegerForLoop {

        private final double[] in = new double[3];
        private final double[] out = new double[3];
        private double[] x;
        private double[] y;
        private double[] z;
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
            x = coordinates[0][0];
            y = coordinates[0][1];
            z = coordinates[0][2];
            threadID = getThreadIndex();
        }

        @Override
        public void run(int lb, int ub) {
            grad.reset(threadID, lb, ub);
            torque.reset(threadID, lb, ub);
            if (lambdaTerm) {
                lambdaGrad.reset(threadID, lb, ub);
                lambdaTorque.reset(threadID, lb, ub);
            }

            // Initialize the local coordinate arrays.
            for (int i = lb; i <= ub; i++) {
                Atom atom = atoms[i];
                x[i] = atom.getX();
                y[i] = atom.getY();
                z[i] = atom.getZ();
                use[i] = atom.getUse();

                    /*
                      Real space Ewald is cutoff at ~7 A, compared to ~12 A for
                      vdW, so the number of neighbors is much more compact. A
                      specific list for real space Ewald is filled during
                      computation of the permanent real space field that
                      includes only evaluated interactions. Subsequent real
                      space loops, especially the SCF, then do not spend time
                      evaluating pairwise distances outside the cutoff.
                     */
                int size = neighborLists[0][i].length;
                if (vaporLists != null) {
                    size = max(size, vaporLists[0][i].length);
                }
                if (realSpaceLists[0][i] == null || realSpaceLists[0][i].length < size) {
                    realSpaceLists[0][i] = new int[size];
                }
            }

            // Expand coordinates.
            List<SymOp> symOps = crystal.spaceGroup.symOps;
            int nSymm = symOps.size();
            for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                SymOp symOp = symOps.get(iSymm);
                double[] xs = coordinates[iSymm][0];
                double[] ys = coordinates[iSymm][1];
                double[] zs = coordinates[iSymm][2];
                for (int i = lb; i <= ub; i++) {
                    in[0] = x[i];
                    in[1] = y[i];
                    in[2] = z[i];
                    crystal.applySymOp(in, out, symOp);
                    xs[i] = out[0];
                    ys[i] = out[1];
                    zs[i] = out[2];
                    int size = neighborLists[iSymm][i].length;
                    if (realSpaceLists[iSymm][i] == null || realSpaceLists[iSymm][i].length < size) {
                        realSpaceLists[iSymm][i] = new int[size];
                    }
                }
            }
        }
    }

    private class RotateMultipolesLoop extends IntegerForLoop {

        // Local variables
        private final double[] localOrigin = new double[3];
        private final double[] xAxis = new double[3];
        private final double[] yAxis = new double[3];
        private final double[] zAxis = new double[3];
        private final double[][] rotmat = new double[3][3];
        private final double[] tempDipole = new double[3];
        private final double[][] tempQuadrupole = new double[3][3];
        private final double[] dipole = new double[3];
        private final double[][] quadrupole = new double[3][3];
        private double chargeScale, dipoleScale, quadrupoleScale, polarizabilityScale;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) {
            List<SymOp> symOps = crystal.spaceGroup.symOps;
            int nSymm = symOps.size();
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                final double[] x = coordinates[iSymm][0];
                final double[] y = coordinates[iSymm][1];
                final double[] z = coordinates[iSymm][2];
                for (int ii = lb; ii <= ub; ii++) {
                    Atom atom = atoms[ii];
                    chargeScale = 1.0;
                    dipoleScale = 1.0;
                    quadrupoleScale = 1.0;
                    polarizabilityScale = 1.0;
                    if (atom.applyLambda()) {
                        chargeScale = lambdaScaleMultipoles;
                        dipoleScale = lambdaScaleMultipoles;
                        quadrupoleScale = lambdaScaleMultipoles;
                        polarizabilityScale = lambdaScaleMultipoles;
                    }
                    if (!useCharges) {
                        chargeScale = 0.0;
                    }
                    if (!useDipoles) {
                        dipoleScale = 0.0;
                    }
                    if (!useQuadrupoles) {
                        quadrupoleScale = 0.0;
                    }
                    final double[] in = localMultipole[ii];
                    final double[] out = globalMultipole[iSymm][ii];
                    double elecScale = 1.0;
                    if (!atom.getElectrostatics()) {
                        elecScale = 0.0;
                    }
                    if (rotateMultipoles) {
                        localOrigin[0] = x[ii];
                        localOrigin[1] = y[ii];
                        localOrigin[2] = z[ii];
                        int[] referenceSites = axisAtom[ii];
                        for (int i = 0; i < 3; i++) {
                            zAxis[i] = 0.0;
                            xAxis[i] = 0.0;
                            dipole[i] = 0.0;
                            for (int j = 0; j < 3; j++) {
                                quadrupole[i][j] = 0.0;
                            }
                        }
                        if (referenceSites == null || referenceSites.length < 2) {
                            out[t000] = in[0] * chargeScale * elecScale;
                            out[t100] = 0.0;
                            out[t010] = 0.0;
                            out[t001] = 0.0;
                            out[t200] = 0.0;
                            out[t020] = 0.0;
                            out[t002] = 0.0;
                            out[t110] = 0.0;
                            out[t101] = 0.0;
                            out[t011] = 0.0;
                            PolarizeType polarizeType = atoms[ii].getPolarizeType();
                            polarizability[ii] = polarizeType.polarizability * polarizabilityScale * elecScale;
                            continue;
                        }
                        switch (frame[ii]) {
                            case BISECTOR:
                                int index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                xAxis[0] = x[index];
                                xAxis[1] = y[index];
                                xAxis[2] = z[index];
                                diff(zAxis, localOrigin, zAxis);
                                norm(zAxis, zAxis);
                                diff(xAxis, localOrigin, xAxis);
                                norm(xAxis, xAxis);
                                sum(xAxis, zAxis, zAxis);
                                norm(zAxis, zAxis);
                                rotmat[0][2] = zAxis[0];
                                rotmat[1][2] = zAxis[1];
                                rotmat[2][2] = zAxis[2];
                                double dot = dot(xAxis, zAxis);
                                scalar(zAxis, dot, zAxis);
                                diff(xAxis, zAxis, xAxis);
                                norm(xAxis, xAxis);
                                rotmat[0][0] = xAxis[0];
                                rotmat[1][0] = xAxis[1];
                                rotmat[2][0] = xAxis[2];
                                break;
                            case ZTHENBISECTOR:
                                index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                xAxis[0] = x[index];
                                xAxis[1] = y[index];
                                xAxis[2] = z[index];
                                index = referenceSites[2];
                                yAxis[0] = x[index];
                                yAxis[1] = y[index];
                                yAxis[2] = z[index];
                                diff(zAxis, localOrigin, zAxis);
                                norm(zAxis, zAxis);
                                rotmat[0][2] = zAxis[0];
                                rotmat[1][2] = zAxis[1];
                                rotmat[2][2] = zAxis[2];
                                diff(xAxis, localOrigin, xAxis);
                                norm(xAxis, xAxis);
                                diff(yAxis, localOrigin, yAxis);
                                norm(yAxis, yAxis);
                                sum(xAxis, yAxis, xAxis);
                                norm(xAxis, xAxis);
                                dot = dot(xAxis, zAxis);
                                scalar(zAxis, dot, zAxis);
                                diff(xAxis, zAxis, xAxis);
                                norm(xAxis, xAxis);
                                rotmat[0][0] = xAxis[0];
                                rotmat[1][0] = xAxis[1];
                                rotmat[2][0] = xAxis[2];
                                break;
                            case ZTHENX:
                            default:
                                index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                xAxis[0] = x[index];
                                xAxis[1] = y[index];
                                xAxis[2] = z[index];
                                diff(zAxis, localOrigin, zAxis);
                                norm(zAxis, zAxis);
                                rotmat[0][2] = zAxis[0];
                                rotmat[1][2] = zAxis[1];
                                rotmat[2][2] = zAxis[2];
                                diff(xAxis, localOrigin, xAxis);
                                dot = dot(xAxis, zAxis);
                                scalar(zAxis, dot, zAxis);
                                diff(xAxis, zAxis, xAxis);
                                norm(xAxis, xAxis);
                                rotmat[0][0] = xAxis[0];
                                rotmat[1][0] = xAxis[1];
                                rotmat[2][0] = xAxis[2];
                        }
                        // Finally the Y elements.
                        rotmat[0][1] = rotmat[2][0] * rotmat[1][2] - rotmat[1][0] * rotmat[2][2];
                        rotmat[1][1] = rotmat[0][0] * rotmat[2][2] - rotmat[2][0] * rotmat[0][2];
                        rotmat[2][1] = rotmat[1][0] * rotmat[0][2] - rotmat[0][0] * rotmat[1][2];
                        // Do the rotation.
                        tempDipole[0] = in[t100];
                        tempDipole[1] = in[t010];
                        tempDipole[2] = in[t001];
                        tempQuadrupole[0][0] = in[t200];
                        tempQuadrupole[1][1] = in[t020];
                        tempQuadrupole[2][2] = in[t002];
                        tempQuadrupole[0][1] = in[t110];
                        tempQuadrupole[0][2] = in[t101];
                        tempQuadrupole[1][2] = in[t011];
                        tempQuadrupole[1][0] = in[t110];
                        tempQuadrupole[2][0] = in[t101];
                        tempQuadrupole[2][1] = in[t011];

                        // Check for chiral flipping.
                        if (frame[ii] == MultipoleFrameDefinition.ZTHENX
                                && referenceSites.length == 3) {
                            localOrigin[0] = x[ii];
                            localOrigin[1] = y[ii];
                            localOrigin[2] = z[ii];
                            int index = referenceSites[0];
                            zAxis[0] = x[index];
                            zAxis[1] = y[index];
                            zAxis[2] = z[index];
                            index = referenceSites[1];
                            xAxis[0] = x[index];
                            xAxis[1] = y[index];
                            xAxis[2] = z[index];
                            index = referenceSites[2];
                            yAxis[0] = x[index];
                            yAxis[1] = y[index];
                            yAxis[2] = z[index];
                            diff(localOrigin, yAxis, localOrigin);
                            diff(zAxis, yAxis, zAxis);
                            diff(xAxis, yAxis, xAxis);
                            double c1 = zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1];
                            double c2 = xAxis[1] * localOrigin[2] - xAxis[2] * localOrigin[1];
                            double c3 = localOrigin[1] * zAxis[2] - localOrigin[2] * zAxis[1];
                            double vol = localOrigin[0] * c1 + zAxis[0] * c2 + xAxis[0] * c3;
                            if (vol < 0.0) {
                                tempDipole[1] = -tempDipole[1];
                                tempQuadrupole[0][1] = -tempQuadrupole[0][1];
                                tempQuadrupole[1][0] = -tempQuadrupole[1][0];
                                tempQuadrupole[1][2] = -tempQuadrupole[1][2];
                                tempQuadrupole[2][1] = -tempQuadrupole[2][1];
                            }
                        }
                        for (int i = 0; i < 3; i++) {
                            double[] rotmati = rotmat[i];
                            double[] quadrupolei = quadrupole[i];
                            for (int j = 0; j < 3; j++) {
                                double[] rotmatj = rotmat[j];
                                dipole[i] += rotmati[j] * tempDipole[j];
                                if (j < i) {
                                    quadrupolei[j] = quadrupole[j][i];
                                } else {
                                    for (int k = 0; k < 3; k++) {
                                        double[] localQuadrupolek = tempQuadrupole[k];
                                        quadrupolei[j] += rotmati[k]
                                                * (rotmatj[0] * localQuadrupolek[0]
                                                + rotmatj[1] * localQuadrupolek[1]
                                                + rotmatj[2] * localQuadrupolek[2]);
                                    }
                                }
                            }
                        }
                        out[t000] = in[0] * chargeScale * elecScale;
                        out[t100] = dipole[0] * dipoleScale * elecScale;
                        out[t010] = dipole[1] * dipoleScale * elecScale;
                        out[t001] = dipole[2] * dipoleScale * elecScale;
                        out[t200] = quadrupole[0][0] * quadrupoleScale * elecScale;
                        out[t020] = quadrupole[1][1] * quadrupoleScale * elecScale;
                        out[t002] = quadrupole[2][2] * quadrupoleScale * elecScale;
                        out[t110] = quadrupole[0][1] * quadrupoleScale * elecScale;
                        out[t101] = quadrupole[0][2] * quadrupoleScale * elecScale;
                        out[t011] = quadrupole[1][2] * quadrupoleScale * elecScale;
                    } else {
                        // No multipole rotation for isolating torque vs. non-torque pieces of the multipole energy gradient.
                        out[t000] = in[t000] * chargeScale * elecScale;
                        out[t100] = in[t100] * dipoleScale * elecScale;
                        out[t010] = in[t010] * dipoleScale * elecScale;
                        out[t001] = in[t001] * dipoleScale * elecScale;
                        out[t200] = in[t200] * quadrupoleScale * elecScale;
                        out[t020] = in[t020] * quadrupoleScale * elecScale;
                        out[t002] = in[t002] * quadrupoleScale * elecScale;
                        out[t110] = in[t110] * quadrupoleScale * elecScale;
                        out[t101] = in[t101] * quadrupoleScale * elecScale;
                        out[t011] = in[t011] * quadrupoleScale * elecScale;
                    }
                    PolarizeType polarizeType = atoms[ii].getPolarizeType();
                    polarizability[ii] = polarizeType.polarizability * polarizabilityScale * elecScale;
                }
            }
        }
    }
}
