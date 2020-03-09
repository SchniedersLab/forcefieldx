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
package ffx.potential.nonbonded.implicit;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;

import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE;
import ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWType;
import static ffx.potential.nonbonded.VanDerWaalsForm.getCombinedEps;
import static ffx.potential.nonbonded.VanDerWaalsForm.getCombinedRadius;

/**
 * Parallel calculation of continuum dispersion energy via pairwise descreening.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DispersionRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(DispersionRegion.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    protected Atom[] atoms;
    /**
     * Periodic boundary conditions and symmetry.
     */
    private Crystal crystal;
    /**
     * Flag to indicate if an atom should be included.
     */
    private boolean[] use = null;
    /**
     * Neighbor lists for each atom and symmetry operator.
     */
    private int[][][] neighborLists;
    /**
     * Cartesian coordinates of each atom.
     */
    private double[] x, y, z;
    /**
     * GK cut-off distance squared.
     */
    private double cut2;
    /**
     * Boolean flag to indicate GK will be scaled by the lambda state variable.
     */
    private boolean gradient = false;
    /**
     * Gradient array for each thread.
     */
    private AtomicDoubleArray3D grad;
    /**
     * Radius of each atom for calculation of dispersion energy.
     */
    private double[] rDisp;
    private double[] cDisp;
    private final DispersionLoop[] dispersionLoop;
    private final SharedDouble sharedDispersion;
    /**
     * The dispersion integral HCT overlap scale factor.
     */
    private static final double DEFAULT_DISP_OVERLAP_FACTOR = 0.81;
    /**
     * The dispersion integral begin for each atom at:
     * Rmin + DISPERSION_OFFSET
     */
    public static final double DEFAULT_DISPERSION_OFFSET = 0.826;
    /**
     * Each solute atom blocks dispersion interactions with solvent:
     * Rmin + SOLUTE_OFFSET
     */
    public static final double DEFAULT_SOLUTE_OFFSET = 0.0;
    /**
     * Conversion between solute-solvent dispersion enthalpy and solute-solvent dispersion free energy.
     */
    private static final double SLEVY = 1.0;
    /**
     * Number density of water (water per A^3).
     */
    private static final double AWATER = 0.033428;
    /**
     * AMOEBA '03 Water oxygen epsilon.
     */
    private static final double EPSO = 0.1100;
    /**
     * AMOEBA '03 Water hydrogen epsilon.
     */
    private static final double EPSH = 0.0135;
    /**
     * AMOEBA '03 Water oxygen Rmin (A).
     */
    private static final double RMINO = 1.7025;
    /**
     * AMOEBA '03 Water hydrogen Rmin (A).
     */
    private static final double RMINH = 1.3275;
    /**
     * AMOEBA epsilon rule.
     */
    private static final EPSILON_RULE epsilonRule = EPSILON_RULE.HHG;
    /**
     * AMOEBA radius rule.
     */
    private static final RADIUS_RULE radiusRule = RADIUS_RULE.CUBIC_MEAN;

    /**
     * Where the dispersion integral begins for each atom (A):
     * Rmin + dispersionOffset
     */
    private double dispersionOffest;
    /**
     * Each solute atom blocks dispersion interactions with solvent with this radius (A):
     * Rmin + soluteOffset
     */
    private double soluteOffset;
    /**
     * The dispersion integral HCT overlap scale factor (unitless).
     */
    private double dispersionOverlapFactor;

    /**
     * DispersionRegion constructor.
     *
     * @param nt         Number of threads.
     * @param atoms      Atom array.
     * @param forceField ForceField in use.
     */
    public DispersionRegion(int nt, Atom[] atoms, ForceField forceField) {
        dispersionLoop = new DispersionLoop[nt];
        for (int i = 0; i < nt; i++) {
            dispersionLoop[i] = new DispersionLoop();
        }
        sharedDispersion = new SharedDouble();

        dispersionOffest = forceField.getDouble("DISPERSION_OFFSET", DEFAULT_DISPERSION_OFFSET);
        soluteOffset = forceField.getDouble("SOLUTE_OFFSET", DEFAULT_SOLUTE_OFFSET);
        dispersionOverlapFactor = forceField.getDouble("DISP_OVERLAP_FACTOR", DEFAULT_DISP_OVERLAP_FACTOR);

        allocate(atoms);
    }

    public double getDispersionOffest() {
        return dispersionOffest;
    }

    /**
     * The dispersion integral begins offset from the vdW radius.
     *
     * @param dispersionOffest The dispersion integral offset.
     */
    public void setDispersionOffest(double dispersionOffest) {
        this.dispersionOffest = dispersionOffest;
        // Update the maximum dispersion energy.
        maxDispersionEnergy();
    }

    public double getDispersionOverlapFactor() {
        return dispersionOverlapFactor;
    }

    /**
     * Set the dispersion overlap HCT scale factor.
     *
     * @param dispersionOverlapFactor The dispersion integral HCT scale factor.
     */
    public void setDispersionOverlapFactor(double dispersionOverlapFactor) {
        this.dispersionOverlapFactor = dispersionOverlapFactor;
    }

    public double getSoluteOffset() {
        return soluteOffset;
    }

    public void setSoluteOffset(double soluteOffset) {
        this.soluteOffset = soluteOffset;
    }

    /**
     * The dispersion integral begins offset from the vdW radius.
     *
     * @return the dispersion integral offset.
     */
    public double getDispersionOffset() {
        return dispersionOffest;
    }

    /**
     * Allocate storage given the Atom array.
     *
     * @param atoms Atom array in use.
     */
    public void allocate(Atom[] atoms) {
        this.atoms = atoms;
        int nAtoms = atoms.length;
        cDisp = new double[nAtoms];
        rDisp = new double[nAtoms];
        if (logger.isLoggable(Level.FINEST)) {
            logger.finest(" Dispersion radii:");
        }

        for (int i = 0; i < nAtoms; i++) {
            VDWType type = atoms[i].getVDWType();
            double rmini = type.radius;
            rDisp[i] = rmini / 2.0;
            if (logger.isLoggable(Level.FINEST)) {
                logger.finest(format(" %d %s %8.6f", i, atoms[i].toString(), rDisp[i]));
            }
        }
        maxDispersionEnergy();
    }

    /**
     * Initialize the DispersionRegion for energy calculation.
     *
     * @param atoms         Atom array.
     * @param crystal       Crystal for periodic boundary conditions.
     * @param use           Flag to indicate an atom is to be used.
     * @param neighborLists Neighbor-list for each atom.
     * @param x             X-coordinate array.
     * @param y             Y-coordinate array.
     * @param z             Z-coordinate array.
     * @param cut2          The cut-off distance squared.
     * @param gradient      If true, compute the gradient.
     * @param grad          Array to store the gradient.
     */
    public void init(Atom[] atoms, Crystal crystal, boolean[] use, int[][][] neighborLists,
                     double[] x, double[] y, double[] z, double cut2,
                     boolean gradient, AtomicDoubleArray3D grad) {
        this.atoms = atoms;
        this.crystal = crystal;
        this.use = use;
        this.neighborLists = neighborLists;
        this.x = x;
        this.y = y;
        this.z = z;
        this.cut2 = cut2;
        this.gradient = gradient;
        this.grad = grad;
    }

    public double getEnergy() {
        return sharedDispersion.get();
    }

    /**
     * Compute the maximum Dispersion energy for each atom in isolation. The
     * loss of dispersion energy due to descreening of other atoms is then
     * calculated in the DispersionLoop.
     */
    private void maxDispersionEnergy() {
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            VDWType type = atoms[i].getVDWType();
            double epsi = type.wellDepth;
            double rmini = type.radius / 2.0;
            if (rDisp[i] > 0.0 && epsi > 0.0) {
                double emixo = getCombinedEps(EPSO, epsi, epsilonRule);
                double emixh = getCombinedEps(EPSH, epsi, epsilonRule);
                double rmixo = getCombinedRadius(RMINO, rmini, radiusRule);
                double rmixh = getCombinedRadius(RMINH, rmini, radiusRule);
                // Apply the dispersion offset to start the integral beyond the atomic radius of atom i.
                double ri = rDisp[i] + dispersionOffest;
                // Integral with two water hydrogen atoms.
                cDisp[i] = 2.0 * tailCorrection(ri, emixh, rmixh);
                // Integral with a water oxygen atom.
                cDisp[i] += tailCorrection(ri, emixo, rmixo);
            }
            cDisp[i] = SLEVY * AWATER * cDisp[i];
        }
    }

    /**
     * Compute a Buffered-14-7 tail correction.
     *
     * @param ri   The separation distance where the integal begins.
     * @param eps  The mixed eps value.
     * @param rmin The mixed rmin values.
     * @return The tail correction.
     */
    private static double tailCorrection(double ri, double eps, double rmin) {
        if (ri < rmin) {
            double r3 = ri * ri * ri;
            double rmin3 = rmin * rmin * rmin;
            return -4.0 * PI * eps * (rmin3 - r3) / 3.0 - eps * 18.0 * PI * rmin3 / 11.0;
        } else {
            double ri2 = ri * ri;
            double ri4 = ri2 * ri2;
            double ri7 = ri * ri2 * ri4;
            double ri11 = ri7 * ri4;
            double rmin2 = rmin * rmin;
            double rmin4 = rmin2 * rmin2;
            double rmin7 = rmin * rmin2 * rmin4;
            return 2.0 * PI * eps * rmin7 * (2.0 * rmin7 - 11.0 * ri7) / (11.0 * ri11);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void start() {
        sharedDispersion.set(0.0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        try {
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, dispersionLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception computing Dispersion energy in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute Dispersion energy for a range of atoms via pairwise
     * descreening.
     *
     * @since 1.0
     */
    private class DispersionLoop extends IntegerForLoop {
        private double edisp;
        private final double[] dx_local;
        private double r, r2, r3;
        private double xr, yr, zr;
        private int threadID;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        DispersionLoop() {
            dx_local = new double[3];
        }

        @Override
        public void start() {
            threadID = getThreadIndex();
            edisp = 0;
        }

        @Override
        public void finish() {
            sharedDispersion.addAndGet(edisp);
        }

        @Override
        public void run(int lb, int ub) {
            for (int i = lb; i <= ub; i++) {
                if (!use[i]) {
                    continue;
                }

                // Begin with the limit of atom alone in solvent.
                edisp += cDisp[i];

                // Now descreen over neighbors.
                double sum = 0.0;
                final double xi = x[i];
                final double yi = y[i];
                final double zi = z[i];
                int[] list = neighborLists[0][i];
                for (int k : list) {
                    final double rk = rDisp[k];
                    if (i != k && rk > 0.0 && use[k]) {
                        dx_local[0] = xi - x[k];
                        dx_local[1] = yi - y[k];
                        dx_local[2] = zi - z[k];
                        r2 = crystal.image(dx_local);
                        if (r2 > cut2) {
                            continue;
                        }
                        xr = dx_local[0];
                        yr = dx_local[1];
                        zr = dx_local[2];
                        r = sqrt(r2);
                        r3 = r * r2;
                        // Atom i descreened by atom k.
                        sum += removeSoluteDispersion(i, k);
                        // Flip the sign on {xr, yr, zr};
                        xr = -xr;
                        yr = -yr;
                        zr = -zr;
                        // Atom k descreened by atom i.
                        sum += removeSoluteDispersion(k, i);
                    }
                }
                // Subtract descreening.
                edisp -= SLEVY * AWATER * sum;
            }
        }

        /**
         * Remove solute dispersion between atom i and solvent due to blocking of solveny by atom k.
         *
         * @param i Atom i index.
         * @param k Atom k index.
         * @return Reduction of dispersion.
         */
        private double removeSoluteDispersion(int i, int k) {
            double sum = 0.0;
            VDWType type = atoms[i].getVDWType();
            double epsi = type.wellDepth;
            double emixo = getCombinedEps(EPSO, epsi, epsilonRule);
            double emixh = getCombinedEps(EPSH, epsi, epsilonRule);
            double rmini = type.radius / 2.0;
            double rmixo = getCombinedRadius(RMINO, rmini, radiusRule);
            double rmixh = getCombinedRadius(RMINH, rmini, radiusRule);
            double rmixo7 = pow(rmixo, 7);
            double rmixh7 = pow(rmixh, 7);
            // Apply the offset to start the integral beyond the atomic radius of atom i.
            double ri = rDisp[i] + dispersionOffest;
            // Atom k blocks interaction of atom i with solvent.
            double rk = rDisp[k] + soluteOffset;
            double sk = rk * dispersionOverlapFactor;
            double sk2 = sk * sk;
            if (ri < r + sk) {
                double de = 0.0;
                double rmax = max(ri, r - sk);
                double lik = rmax;
                double lik2 = lik * lik;
                double lik3 = lik2 * lik;
                double lik4 = lik3 * lik;
                // Interaction with water oxygen from lik to Rmin.
                if (lik < rmixo) {
                    double uik = min(r + sk, rmixo);
                    double uik2 = uik * uik;
                    double uik3 = uik2 * uik;
                    double uik4 = uik3 * uik;
                    sum += integralBeforeRMin(emixo, r, r2, sk2, lik2, lik3, lik4, uik2, uik3, uik4);
                    if (gradient) {
                        de += integralBeforeRminDerivative(ri, emixo, rmixo, r, r2, r3, sk, sk2,
                                lik, lik2, lik3, uik, uik2, uik3);
                    }
                }
                // Interaction with water hydrogen from lik to Rmin.
                if (lik < rmixh) {
                    double uik = min(r + sk, rmixh);
                    double uik2 = uik * uik;
                    double uik3 = uik2 * uik;
                    double uik4 = uik3 * uik;
                    sum += 2.0 * integralBeforeRMin(emixh, r, r2, sk2, lik2, lik3, lik4, uik2, uik3, uik4);
                    if (gradient) {
                        de += 2.0 * integralBeforeRminDerivative(ri, emixh, rmixh, r, r2, r3, sk, sk2,
                                lik, lik2, lik3, uik, uik2, uik3);
                    }
                }
                double uik = r + sk;
                double uik2 = uik * uik;
                double uik3 = uik2 * uik;
                double uik4 = uik3 * uik;
                double uik5 = uik4 * uik;
                double uik6 = uik5 * uik;
                double uik10 = uik5 * uik5;
                double uik11 = uik10 * uik;
                double uik12 = uik11 * uik;
                double uik13 = uik12 * uik;
                // Interaction with water oxygen beyond Rmin, from lik to uik = r + sk.
                if (uik > rmixo) {
                    lik = max(rmax, rmixo);
                    lik2 = lik * lik;
                    lik3 = lik2 * lik;
                    lik4 = lik3 * lik;
                    double lik5 = lik4 * lik;
                    double lik6 = lik5 * lik;
                    double lik10 = lik5 * lik5;
                    double lik11 = lik10 * lik;
                    double lik12 = lik11 * lik;
                    sum += integratlAfterRmin(emixo, rmixo7, r, r2, sk2,
                            lik, lik2, lik3, lik4, lik5, lik10, lik11, lik12,
                            uik, uik2, uik3, uik4, uik5, uik10, uik11, uik12);
                    if (gradient) {
                        double lik13 = lik12 * lik;
                        de += integratlAfterRminDerivative(ri, emixo, rmixo, rmixo7, rmax, r, r2, r3, sk, sk2,
                                lik, lik2, lik3, lik5, lik6, lik12, lik13, uik, uik2, uik3, uik6, uik13);
                    }

                }
                // Interaction with water hydrogen beyond Rmin, from lik to uik = r + sk.
                if (uik > rmixh) {
                    lik = max(rmax, rmixh);
                    lik2 = lik * lik;
                    lik3 = lik2 * lik;
                    lik4 = lik3 * lik;
                    double lik5 = lik4 * lik;
                    double lik6 = lik5 * lik;
                    double lik10 = lik5 * lik5;
                    double lik11 = lik10 * lik;
                    double lik12 = lik11 * lik;
                    sum += 2.0 * integratlAfterRmin(emixh, rmixh7, r, r2, sk2,
                            lik, lik2, lik3, lik4, lik5, lik10, lik11, lik12,
                            uik, uik2, uik3, uik4, uik5, uik10, uik11, uik12);
                    if (gradient) {
                        double lik13 = lik12 * lik;
                        de += 2.0 * integratlAfterRminDerivative(ri, emixh, rmixh, rmixh7, rmax, r, r2, r3, sk, sk2,
                                lik, lik2, lik3, lik5, lik6, lik12, lik13, uik, uik2, uik3, uik6, uik13);
                    }
                }
                // Increment the individual dispersion gradient components.
                if (gradient) {
                    de = -de / r * SLEVY * AWATER;
                    double dedx = de * xr;
                    double dedy = de * yr;
                    double dedz = de * zr;
                    grad.add(threadID, i, dedx, dedy, dedz);
                    grad.sub(threadID, k, dedx, dedy, dedz);
                }
            }
            return sum;
        }

        /**
         * Integrate over the constant portion of the WCA disperion interaction Uwca(x) = eps; x < rmin.
         *
         * @param eps  The well depth.
         * @param r    The separation between the current atom and solvent blocking atom.
         * @param r2   The separation squared.
         * @param sk2  The scaled size of the solvent blocking atom squared.
         * @param lik2 The beginning of the integral squared.
         * @param lik3 The beginning of the integral to the third.
         * @param lik4 The beginning of the integral to the fourth.
         * @param uik2 The end of the integral squared.
         * @param uik3 The end of the integral to the third.
         * @param uik4 The end of the integral to the fourth.
         * @return The integral.
         */
        private double integralBeforeRMin(double eps, double r, double r2, double sk2,
                                          double lik2, double lik3, double lik4,
                                          double uik2, double uik3, double uik4) {
            return -eps * (4.0 * PI / (48.0 * r) * (3.0 * (lik4 - uik4)
                    - 8.0 * r * (lik3 - uik3) + 6.0 * (r2 - sk2) * (lik2 - uik2)));

        }

        /**
         * Derivative  of the integral over the constant portion of the WCA disperion interaction Uwca(x) = eps; x < rmin.
         *
         * @param ri   The beginning of the integral.
         * @param eps  The well depth.
         * @param rmin The Rmin value.
         * @param r    The separation between the current atom and solvent blocking atom.
         * @param r2   The separation squared.
         * @param r3   The separation to the third.
         * @param sk   The scaled size of the solvent blocking atom.
         * @param sk2  The scaled size of the solvent blocking atom squared.
         * @param lik  The beginning of the integral.
         * @param lik2 The beginning of the integral squared.
         * @param lik3 The beginning of the integral to the third.
         * @param uik  The end of the integral.
         * @param uik2 The end of the integral squared.
         * @param uik3 The end of the integral to the third.
         * @return The derivative of the integral.
         */
        private double integralBeforeRminDerivative(double ri, double eps, double rmin,
                                                    double r, double r2, double r3,
                                                    double sk, double sk2,
                                                    double lik, double lik2, double lik3,
                                                    double uik, double uik2, double uik3) {
            double dl;
            if (ri > r - sk) {
                dl = (-lik2 + 2.0 * r2 + 2.0 * sk2) * lik2;
            } else {
                dl = (-lik3 + 4.0 * lik2 * r - 6.0 * lik * r2 + 2.0 * lik * sk2 + 4.0 * r3 - 4.0 * r * sk2) * lik;
            }
            double du;
            if (r + sk > rmin) {
                du = -(-uik2 + 2.0 * r2 + 2.0 * sk2) * uik2;
            } else {
                du = -(-uik3 + 4.0 * uik2 * r - 6.0 * uik * r2 + 2.0 * uik * sk2 + 4.0 * r3 - 4.0 * r * sk2) * uik;
            }
            return -eps * PI * (dl + du) / (4.0 * r2);
        }

        /**
         * @param eps   The well depth.
         * @param rmin7 The rmin value to the seventh.
         * @param r     The separation between the current atom and solvent blocking atom.
         * @param r2    The separation squared.
         * @param sk2   The scaled size of the solvent blocking atom squared.
         * @param lik   The beginning of the integral.
         * @param lik2  The beginning of the integral squared.
         * @param lik3  The beginning of the integral to the third.
         * @param lik4  The beginning of the integral to the fourth.
         * @param lik5  The beginning of the integral to the fifth.
         * @param lik10 The beginning of the integral to the tenth.
         * @param lik11 The beginning of the integral to the eleventh.
         * @param lik12 The beginning of the integral to the twelfth.
         * @param uik   The end of the integral.
         * @param uik2  The end of the integral squared.
         * @param uik3  The end of the integral to the third.
         * @param uik4  The end of the integral to the fourth.
         * @param uik5  The end of the integral to the fifth.
         * @param uik10 The end of the integral to the tenth.
         * @param uik11 The end of the integral to the eleventh.
         * @param uik12 The end of the integral to the twelfth.
         * @return The value of the integral.
         */
        private double integratlAfterRmin(double eps, double rmin7, double r, double r2, double sk2,
                                          double lik, double lik2, double lik3, double lik4, double lik5, double lik10, double lik11, double lik12,
                                          double uik, double uik2, double uik3, double uik4, double uik5, double uik10, double uik11, double uik12) {
            double er7 = eps * rmin7;
            double term = 4.0 * PI / (120.0 * r * lik5 * uik5) * (15.0 * uik * lik * r * (uik4 - lik4)
                    - 10.0 * uik2 * lik2 * (uik3 - lik3) + 6.0 * (sk2 - r2) * (uik5 - lik5));
            double term2 = 4.0 * PI / (2640.0 * r * lik12 * uik12) * (120.0 * uik * lik * r * (uik11 - lik11)
                    - 66.0 * uik2 * lik2 * (uik10 - lik10) + 55.0 * (sk2 - r2) * (uik12 - lik12));
            double idisp = -2.0 * er7 * term;
            double irep = er7 * rmin7 * term2;
            return irep + idisp;
        }

        /**
         * @param ri    The beginning of the integral.
         * @param eps   The eps value.
         * @param rmin  The rmin value to the seventh.
         * @param rmin7 The rmin value to the seventh.
         * @param r     The separation between the current atom and solvent blocking atom.
         * @param r2    The separation squared.
         * @param r3    The separation cubed.
         * @param sk    The scaled size of the solvent blocking atom.
         * @param sk2   The scaled size of the solvent blocking atom squared.
         * @param lik   The beginning of the integral.
         * @param lik2  The beginning of the integral squared.
         * @param lik3  The beginning of the integral to the third.
         * @param lik5  The beginning of the integral to the fifth.
         * @param lik6  The beginning of the integral to the sixth.
         * @param lik12 The beginning of the integral to the twelfth.
         * @param lik13 The beginning of the integral to the thirteenth.
         * @param uik   The end of the integral.
         * @param uik2  The end of the integral squared.
         * @param uik3  The end of the integral to the third.
         * @param uik6  The end of the integral to the sixth.
         * @param uik13 The end of the integral to the thirteenth.
         * @return The value of the integral derivative.
         */
        private double integratlAfterRminDerivative(double ri, double eps, double rmin, double rmin7, double rmax,
                                                    double r, double r2, double r3, double sk, double sk2,
                                                    double lik, double lik2, double lik3, double lik5, double lik6, double lik12, double lik13,
                                                    double uik, double uik2, double uik3, double uik6, double uik13) {
            double er7 = eps * rmin7;
            double lowerTerm = lik2 * r + r3 - r * sk2;
            double upperTerm = uik2 * r + r3 - r * sk2;

            double dl;
            if (ri > r - sk || rmax < rmin) {
                dl = -(-5.0 * lik2 + 3.0 * r2 + 3.0 * sk2) / lik5;
            } else {
                dl = (5.0 * lik3 - 33.0 * lik * r2 - 3.0 * lik * sk2 + 15.0 * lowerTerm) / lik6;
            }
            double du = -(5.0 * uik3 - 33.0 * uik * r2 - 3.0 * uik * sk2 + 15.0 * upperTerm) / uik6;
            double de = -2.0 * PI * er7 * (dl + du) / (15.0 * r2);

            if (ri > r - sk || rmax < rmin) {
                dl = -(-6.0 * lik2 + 5.0 * r2 + 5.0 * sk2) / lik12;
            } else {
                dl = (6.0 * lik3 - 125.0 * lik * r2 - 5.0 * lik * sk2 + 60.0 * lowerTerm) / lik13;
            }
            du = -(6.0 * uik3 - 125.0 * uik * r2 - 5.0 * uik * sk2 + 60.0 * upperTerm) / uik13;
            de += PI * er7 * rmin7 * (dl + du) / (60.0 * r2);

            return de;
        }

    }
}
