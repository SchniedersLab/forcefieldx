// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.potential.nonbonded.implicit;

import static ffx.potential.nonbonded.implicit.BornTanhRescaling.tanhRescalingChainRule;
import static ffx.potential.nonbonded.implicit.NeckIntegral.getNeckConstants;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.*;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.EnergyException;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel computation of Born radii chain rule terms via the Grycuk method.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class BornGradRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(BornRadiiRegion.class.getName());
    private static final double PI4_3 = 4.0 / 3.0 * PI;
    /**
     * Constant factor used with quadrupoles.
     */
    private static final double oneThird = 1.0 / 3.0;

    private final BornCRLoop[] bornCRLoop;
    /**
     * An ordered array of atoms in the system.
     */
    protected Atom[] atoms;
    /**
     * Periodic boundary conditions and symmetry.
     */
    private Crystal crystal;
    /**
     * Atomic coordinates for each symmetry operator.
     */
    private double[][][] sXYZ;
    /**
     * Neighbor lists for each atom and symmetry operator.
     */
    private int[][][] neighborLists;
    /**
     * Base radius of each atom.
     */
    private double[] baseRadius;
    /**
     * Descreen radius of each atom.
     */
    private double[] descreenRadius;
    /**
     * Overlap scale factor for each atom, when using the Hawkins, Cramer & Truhlar pairwise
     * descreening algorithm.
     *
     * <p>G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized Models of Aqueous Free Energies
     * of Solvation Based on Pairwise Descreening of Solute Atomic Charges from a Dielectric Medium",
     * J. Phys. Chem., 100, 19824-19839 (1996).
     */
    private double[] overlapScale;
    /**
     * Sneck scaling parameter for each atom. Set based on maximum Sneck scaling parameter and
     * number of bound non-hydrogen atoms
     */
    private double[] neckScale;
    private double descreenOffset;
    /**
     * If true, the descreening integral includes overlaps with the volume of the descreened atom
     */
    private final boolean perfectHCTScale;
    /**
     * Flag to indicate if an atom should be included.
     */
    private boolean[] use;
    /**
     * GK cut-off distance squared.
     */
    private double cut2;
    /**
     * Forces all atoms to be considered during Born radius updates.
     */
    private boolean nativeEnvironmentApproximation;
    /**
     * Born radius of each atom.
     */
    private double[] born;
    /**
     * Gradient array for each thread.
     */
    private AtomicDoubleArray3D grad;
    /**
     * Shared array for computation of Born radii gradient.
     */
    private AtomicDoubleArray sharedBornGrad;
    private final double factor = -pow(PI, oneThird) * pow(6.0, (2.0 * oneThird)) / 9.0;
    private double[] term;
    /**
     * Flag for using neck correction, default to false for now.
     */
    private boolean neckCorrection = false;
    /**
     * Flag for using tanh correction, default to false for now.
     */
    private boolean tanhCorrection = false;
    /**
     * This is the Born Integral prior to rescaling with a tanh function, which is a quantity
     * needed for the computing the derivative of the energy with respect to atomic coordinates.
     */
    private double[] unscaledBornIntegral;

    /**
     * Compute the gradient due to changes in Born radii.
     * @param nt Number of threads.
     * @param neckCorrection Include a neck correction.
     * @param tanhCorrection Include a tanh correction.
     * @param perfectHCTScale Use "perfect" HCT scale factors.
     */
    public BornGradRegion(int nt, boolean neckCorrection, boolean tanhCorrection, boolean perfectHCTScale) {
        bornCRLoop = new BornCRLoop[nt];
        for (int i = 0; i < nt; i++) {
            bornCRLoop[i] = new BornCRLoop();
        }
        this.neckCorrection = neckCorrection;
        this.tanhCorrection = tanhCorrection;
        this.perfectHCTScale = perfectHCTScale;
    }

    /**
     * Execute the InitializationRegion with the passed ParallelTeam.
     *
     * @param parallelTeam The ParallelTeam instance to execute with.
     */
    public void executeWith(ParallelTeam parallelTeam) {
        sharedBornGrad.reduce(parallelTeam, 0, atoms.length - 1);
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = " Exception evaluating Born radii chain rule gradient.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void init(
            Atom[] atoms,
            Crystal crystal,
            double[][][] sXYZ,
            int[][][] neighborLists,
            double[] baseRadius,
            double[] descreenRadius,
            double[] overlapScale,
            double[] neckScale,
            double descreenOffset,
            double[] unscaledBornIntegral,
            boolean[] use,
            double cut2,
            boolean nativeEnvironmentApproximation,
            double[] born,
            AtomicDoubleArray3D grad,
            AtomicDoubleArray sharedBornGrad) {
        this.atoms = atoms;
        this.crystal = crystal;
        this.sXYZ = sXYZ;
        this.neighborLists = neighborLists;
        this.baseRadius = baseRadius;
        this.descreenRadius = descreenRadius;
        this.overlapScale = overlapScale;
        this.neckScale = neckScale;
        this.descreenOffset = descreenOffset;
        this.unscaledBornIntegral = unscaledBornIntegral;
        this.use = use;
        this.cut2 = cut2;
        this.nativeEnvironmentApproximation = nativeEnvironmentApproximation;
        this.born = born;
        this.grad = grad;
        this.sharedBornGrad = sharedBornGrad;
    }

    @Override
    public void start() {
        int nAtoms = atoms.length;
        if (term == null || term.length < nAtoms) {
            term = new double[nAtoms];
        }

        // The "term" array, initialized here, provides the Grycuk method chain rule term for each atom in BornGradRegion
        // The base Grycuk correction sets the Born radius equal to (3/(4*PI) * integral(1/r^6))^(-1/3)
        for (int i = 0; i < nAtoms; i++) {
            double rbi = born[i];
            term[i] = PI4_3 / (rbi * rbi * rbi);
            term[i] = factor / pow(term[i], (4.0 * oneThird));
            if (tanhCorrection) {
                term[i] = term[i] * tanhRescalingChainRule(unscaledBornIntegral[i], baseRadius[i]);
            }
        }
    }

    @Override
    public void run() {
        try {
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, bornCRLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message =
                    "Fatal exception computing Born radii chain rule term in thread "
                            + getThreadIndex()
                            + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute Born radii chain rule terms for a range of atoms via the Grycuk method.
     *
     * @since 1.0
     */
    private class BornCRLoop extends IntegerForLoop {

        private final double[] dx_local;
        private int threadID;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        BornCRLoop() {
            dx_local = new double[3];
        }

        @Override
        public void run(int lb, int ub) {

            double[] x = sXYZ[0][0];
            double[] y = sXYZ[0][1];
            double[] z = sXYZ[0][2];

            int nSymm = crystal.spaceGroup.symOps.size();
            for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                SymOp symOp = crystal.spaceGroup.symOps.get(iSymOp);
                double[][] transOp = new double[3][3];
                double[][] xyz = sXYZ[iSymOp];
                crystal.getTransformationOperator(symOp, transOp);
                for (int i = lb; i <= ub; i++) {
                    if (!nativeEnvironmentApproximation && !use[i]) {
                        continue;
                    }

                    // Check the value of the Born radii chain rule term.
                    double bornGrad = sharedBornGrad.get(i);
                    if (isInfinite(bornGrad) || isNaN(bornGrad)) {
                        throw new EnergyException(
                                format(" %s\n Born radii CR %d %8.3f", atoms[i], i, bornGrad), true);
                    }
                    final double integralStartI = max(baseRadius[i], descreenRadius[i]) + descreenOffset;
                    final double descreenRi = descreenRadius[i];
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double rbi = born[i];
                    int[] list = neighborLists[iSymOp][i];
                    for (int k : list) {
                        if (!nativeEnvironmentApproximation && !use[k]) {
                            continue;
                        }
                        final double integralStartK = max(baseRadius[k], descreenRadius[k]) + descreenOffset;
                        final double descreenRk = descreenRadius[k];
                        double mixedNeckScale =  0.5 * (neckScale[i] + neckScale[k]);

                        if (k != i) {
                            dx_local[0] = xyz[0][k] - xi;
                            dx_local[1] = xyz[1][k] - yi;
                            dx_local[2] = xyz[2][k] - zi;
                            double r2 = crystal.image(dx_local);
                            if (r2 > cut2) {
                                continue;
                            }
                            final double xr = dx_local[0];
                            final double yr = dx_local[1];
                            final double zr = dx_local[2];
                            final double r = sqrt(r2);

                            // Atom i being descreeened by atom k.
                            double sk = overlapScale[k];
                            if (sk > 0.0 && rbi < 50.0 && descreenRk > 0.0) {
                                double de = descreenDerivative(r, r2, integralStartI, descreenRk, sk);
                                if (neckCorrection) {
                                    de += neckDescreenDerivative(r, integralStartI, descreenRk, mixedNeckScale);
                                }
                                if (isInfinite(de) || isNaN(de)) {
                                    logger.warning(format(" Born radii chain rule term is unstable %d %d %16.8f", i, k, de));
                                }
                                double dbr = term[i] * de / r;
                                de = dbr * sharedBornGrad.get(i);
                                incrementGradient(i, k, de, xr, yr, zr, transOp);
                            }

                            // Atom k being descreeened by atom i.
                            double rbk = born[k];
                            double si = overlapScale[i];
                            if (si > 0.0 && rbk < 50.0 && descreenRi > 0.0) {
                                double de = descreenDerivative(r, r2, integralStartK, descreenRi, si);
                                if (neckCorrection) {
                                    de += neckDescreenDerivative(r, integralStartK, descreenRi, mixedNeckScale);
                                }
                                if (isInfinite(de) || isNaN(de)) {
                                    logger.warning(format(" Born radii chain rule term is unstable %d %d %16.8f", k, i, de));
                                }
                                double dbr = term[k] * de / r;
                                de = dbr * sharedBornGrad.get(k);
                                incrementGradient(i, k, de, xr, yr, zr, transOp);
                            }
                        } else if (iSymOp > 0 && rbi < 50.0) {
                            dx_local[0] = xyz[0][k] - xi;
                            dx_local[1] = xyz[1][k] - yi;
                            dx_local[2] = xyz[2][k] - zi;
                            double r2 = crystal.image(dx_local);
                            double sk = overlapScale[k];
                            if (sk > 0.0 && r2 < cut2 && descreenRk > 0.0) {
                                final double xr = dx_local[0];
                                final double yr = dx_local[1];
                                final double zr = dx_local[2];
                                final double r = sqrt(r2);
                                // Atom i being descreeened by atom k.
                                double de = descreenDerivative(r, r2, integralStartI, descreenRk, sk);
                                if (neckCorrection) {
                                    de += neckDescreenDerivative(r, integralStartI, descreenRk, mixedNeckScale);
                                }
                                if (isInfinite(de) || isNaN(de)) {
                                    logger.warning(
                                            format(" Born radii chain rule term is unstable %d %d %d %16.8f", iSymOp, i,
                                                    k, de));
                                }
                                double dbr = term[i] * de / r;
                                de = dbr * sharedBornGrad.get(i);
                                incrementGradient(i, k, de, xr, yr, zr, transOp);
                                // For symmetry mates, atom k is not descreeened by atom i.
                            }
                        }
                    }
                }
            }
        }

        @Override
        public void start() {
            threadID = getThreadIndex();
        }

        private double neckDescreenDerivative(double r, double radius, double radiusK, double sneck) {
            double radiusWater = 1.4;

            if (r > radius + radiusK + 2 * radiusWater) {
                return 0.0;
            }

            // Get Aij and Bij
            double[] constants = getNeckConstants(radius, radiusK);
            //logger.info(format("Grad Inputs: Ri %2.4f Rk %2.4f\nGrad Outputs: Aij %2.4f Bij %2.4f",radius,radiusK,constants[0],constants[1]));

            double Aij = constants[0];
            double Bij = constants[1];

            // Use Aij and Bij to get neck value using derivative of Equation 13 from Aguilar/Onufriev 2010 paper
            double rMinusBij = r - Bij;
            double rMinusBij3 = rMinusBij * rMinusBij * rMinusBij;
            double rMinusBij4 = rMinusBij3 * rMinusBij;
            double radiiMinusr = radius + radiusK + 2.0 * radiusWater - r;
            double radiiMinusr3 = radiiMinusr * radiiMinusr * radiiMinusr;
            double radiiMinusr4 = radiiMinusr3 * radiiMinusr;

            return 4.0 * (sneck * Aij * rMinusBij3 * radiiMinusr4 - sneck * Aij * rMinusBij4 * radiiMinusr3);
        }

        private double descreenDerivative(double r, double r2, double radius, double radiusK,
                                          double hctScale) {
            if (perfectHCTScale) {
                return perfectHCTIntegralDerivative(r, r2, radius, radiusK, hctScale);
            } else {
                return integralDerivative(r, r2, radius, radiusK * hctScale);
            }
        }

        /**
         * Use pairwise descreening to compute derivative of the integral of 1/r^6 with respect to r.
         *
         * @param r            separation distance.
         * @param r2           separation distance squared.
         * @param radius       base radius of descreened atom.
         * @param scaledRadius scaled radius descreening atom.
         * @return the derivative.
         */
        private double integralDerivative(double r, double r2, double radius, double scaledRadius) {
            double de = 0.0;
            // Descreen only if the scaledRadius is greater than zero.
            // and atom I does not engulf atom K.
            if (scaledRadius > 0.0 && (radius < r + scaledRadius)) {
                // Atom i is engulfed by atom k.
                if (radius + r < scaledRadius) {
                    double uik = scaledRadius - r;
                    double uik2 = uik * uik;
                    double uik4 = uik2 * uik2;
                    de = -4.0 * PI / uik4;
                }

                // Lower integration bound depends on atoms sizes and separation.
                double sk2 = scaledRadius * scaledRadius;
                if (radius + r < scaledRadius) {
                    // Atom i is engulfed by atom k.
                    double lik = scaledRadius - r;
                    double lik2 = lik * lik;
                    double lik4 = lik2 * lik2;
                    de = de + 0.25 * PI * (sk2 - 4.0 * scaledRadius * r + 17.0 * r2) / (r2 * lik4);
                } else if (r < radius + scaledRadius) {
                    // Atoms are overlapped, begin integration from ri.
                    double lik = radius;
                    double lik2 = lik * lik;
                    double lik4 = lik2 * lik2;
                    de = de + 0.25 * PI * (2.0 * radius * radius - sk2 - r2) / (r2 * lik4);
                } else {
                    // No overlap between atoms.
                    double lik = r - scaledRadius;
                    double lik2 = lik * lik;
                    double lik4 = lik2 * lik2;
                    de = de + 0.25 * PI * (sk2 - 4.0 * scaledRadius * r + r2) / (r2 * lik4);
                }

                // Upper integration bound is always the same.
                double uik = r + scaledRadius;
                double uik2 = uik * uik;
                double uik4 = uik2 * uik2;
                de = de - 0.25 * PI * (sk2 + 4.0 * scaledRadius * r + r2) / (r2 * uik4);
            }
            return de;
        }

        /**
         * Use pairwise descreening to compute derivative of the integral of 1/r^6 with respect to r.
         *
         * @param r          separation distance.
         * @param r2         separation distance squared.
         * @param radius     base radius of descreened atom.
         * @param radiusK    base radius of descreening atom.
         * @param perfectHCT perfect HCT scale factor.
         * @return the derivative.
         */
        private double perfectHCTIntegralDerivative(double r, double r2, double radius, double radiusK,
                                                    double perfectHCT) {
            double de = 0.0;
            // Descreen only if the scaledRadius is greater than zero.
            // and atom I does not engulf atom K.
            if (radiusK > 0.0 && (radius < r + radiusK)) {
                // Atom i is engulfed by atom k.
                if (radius + r < radiusK) {
                    double uik = radiusK - r;
                    double uik2 = uik * uik;
                    double uik4 = uik2 * uik2;
                    de = -4.0 * PI / uik4;
                }

                // Lower integration bound depends on atoms sizes and separation.
                double sk2 = radiusK * radiusK;
                if (radius + r < radiusK) {
                    // Atom i is engulfed by atom k.
                    double lik = radiusK - r;
                    double lik2 = lik * lik;
                    double lik4 = lik2 * lik2;
                    de = de + 0.25 * PI * (sk2 - 4.0 * radiusK * r + 17.0 * r2) / (r2 * lik4);
                } else if (r < radius + radiusK) {
                    // Atoms are overlapped, begin integration from ri.
                    double lik = radius;
                    double lik2 = lik * lik;
                    double lik4 = lik2 * lik2;
                    de = de + 0.25 * PI * (2.0 * radius * radius - sk2 - r2) / (r2 * lik4);
                } else {
                    // No overlap between atoms.
                    double lik = r - radiusK;
                    double lik2 = lik * lik;
                    double lik4 = lik2 * lik2;
                    de = de + 0.25 * PI * (sk2 - 4.0 * radiusK * r + r2) / (r2 * lik4);
                }

                // Upper integration bound is always the same.
                double uik = r + radiusK;
                double uik2 = uik * uik;
                double uik4 = uik2 * uik2;
                de = de - 0.25 * PI * (sk2 + 4.0 * radiusK * r + r2) / (r2 * uik4);
            }
            return perfectHCT * de;
        }

        /**
         * Accumulate a contribution to the gradient and dU/dX/dL.
         *
         * @param i  index of atom i.
         * @param k  index of atom k.
         * @param dE partial derivative of the energy with respect to R.
         * @param xr x-component of the separation vector.
         * @param yr y-component of the separation vector.
         * @param zr z-component of the separation vector.
         */
        private void incrementGradient(
                int i, int k, double dE, double xr, double yr, double zr, double[][] transOp) {
            double dedx = dE * xr;
            double dedy = dE * yr;
            double dedz = dE * zr;
            grad.add(threadID, i, dedx, dedy, dedz);
            final double dedxk = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
            final double dedyk = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
            final double dedzk = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];
            grad.sub(threadID, k, dedxk, dedyk, dedzk);
        }
    }
}
