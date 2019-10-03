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
package ffx.potential.nonbonded.implicit;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import static ffx.potential.nonbonded.GeneralizedKirkwood.dWater;

/**
 * Parallel calculation of the Generalized Kirkwood induced reaction field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class InducedGKFieldRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(InducedGKFieldRegion.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    protected Atom[] atoms;
    /**
     * Periodic boundary conditions and symmetry.
     */
    private Crystal crystal;
    /**
     * Induced dipoles for each symmetry operator.
     */
    private double[][][] inducedDipole;
    /**
     * Induced dipole chain rule terms for each symmetry operator.
     */
    private double[][][] inducedDipoleCR;
    /**
     * Atomic coordinates for each symmetry operator.
     */
    private double[][][] sXYZ;
    /**
     * Neighbor lists for each atom and symmetry operator.
     */
    private int[][][] neighborLists;
    /**
     * Flag to indicate if an atom should be included.
     */
    private boolean[] use = null;
    /**
     * GK cut-off distance squared.
     */
    private double cut2;
    /**
     * Born radius of each atom.
     */
    private double[] born;
    /**
     * Atomic GK field array.
     */
    private AtomicDoubleArray3D sharedGKField;
    /**
     * Atomic GK field chain-rule array.
     */
    private AtomicDoubleArray3D sharedGKFieldCR;

    /**
     * Empirical constant that controls the GK cross-term.
     */
    private static final double gkc = 2.455;
    /**
     * Kirkwood dipole reaction field constant.
     */
    private final double fd;
    /**
     * Kirkwood quadrupole reaction field constant.
     */
    private final double fq;
    private final InducedGKFieldLoop[] inducedGKFieldLoop;

    public InducedGKFieldRegion(int nt, ForceField forceField) {
        // Set the Kirkwood multipolar reaction field constants.
        double epsilon = forceField.getDouble("GK_EPSILON", dWater);
        fd = 2.0 * (1.0 - epsilon) / (1.0 + 2.0 * epsilon);
        fq = 3.0 * (1.0 - epsilon) / (2.0 + 3.0 * epsilon);

        inducedGKFieldLoop = new InducedGKFieldLoop[nt];
        for (int i = 0; i < nt; i++) {
            inducedGKFieldLoop[i] = new InducedGKFieldLoop();
        }
    }

    public void init(Atom[] atoms, double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     Crystal crystal, double[][][] sXYZ, int[][][] neighborLists,
                     boolean[] use, double cut2, double[] born,
                     AtomicDoubleArray3D sharedGKField, AtomicDoubleArray3D sharedGKFieldCR) {
        this.atoms = atoms;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.crystal = crystal;
        this.sXYZ = sXYZ;
        this.neighborLists = neighborLists;
        this.use = use;
        this.cut2 = cut2;
        this.born = born;
        this.sharedGKField = sharedGKField;
        this.sharedGKFieldCR = sharedGKFieldCR;
    }

    @Override
    public void run() {
        try {
            int nAtoms = atoms.length;
            int threadIndex = getThreadIndex();
            execute(0, nAtoms - 1, inducedGKFieldLoop[threadIndex]);
        } catch (Exception e) {
            String message = "Fatal exception computing GK field in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute the Generalized Kirkwood induced reaction field.
     *
     * @since 1.0
     */
    private class InducedGKFieldLoop extends IntegerForLoop {

        private int threadID;
        private final double[][] a;
        private final double[] gux;
        private final double[] guy;
        private final double[] guz;
        private final double[] dx_local;
        private double xi, yi, zi;
        private double uix, uiy, uiz;
        private double uixCR, uiyCR, uizCR;
        private double rbi;
        private int iSymm;
        private double[][] transOp;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        InducedGKFieldLoop() {
            a = new double[3][2];
            gux = new double[5];
            guy = new double[5];
            guz = new double[5];
            dx_local = new double[3];
            transOp = new double[3][3];
        }

        @Override
        public void run(int lb, int ub) {

            threadID = getThreadIndex();

            double[] x = sXYZ[0][0];
            double[] y = sXYZ[0][1];
            double[] z = sXYZ[0][2];

            int nSymm = crystal.spaceGroup.symOps.size();
            List<SymOp> symOps = crystal.spaceGroup.symOps;
            for (iSymm = 0; iSymm < nSymm; iSymm++) {
                SymOp symOp = symOps.get(iSymm);
                crystal.getTransformationOperator(symOp, transOp);
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    xi = x[i];
                    yi = y[i];
                    zi = z[i];
                    uix = inducedDipole[0][i][0];
                    uiy = inducedDipole[0][i][1];
                    uiz = inducedDipole[0][i][2];
                    uixCR = inducedDipoleCR[0][i][0];
                    uiyCR = inducedDipoleCR[0][i][1];
                    uizCR = inducedDipoleCR[0][i][2];
                    rbi = born[i];
                    int[] list = neighborLists[iSymm][i];
                    for (int k : list) {
                        if (!use[k]) {
                            continue;
                        }
                        inducedGKField(i, k);
                    }

                    // Include the self induced reaction field, which is not in the neighbor list.
                    if (iSymm == 0) {
                        inducedGKField(i, i);
                    }
                }
            }
        }

        private void inducedGKField(int i, int k) {
            dx_local[0] = sXYZ[iSymm][0][k] - xi;
            dx_local[1] = sXYZ[iSymm][1][k] - yi;
            dx_local[2] = sXYZ[iSymm][2][k] - zi;
            final double r2 = crystal.image(dx_local);
            if (r2 > cut2) {
                return;
            }
            double xr = dx_local[0];
            double yr = dx_local[1];
            double zr = dx_local[2];
            double xr2 = xr * xr;
            double yr2 = yr * yr;
            double zr2 = zr * zr;
            final double ukx = inducedDipole[iSymm][k][0];
            final double uky = inducedDipole[iSymm][k][1];
            final double ukz = inducedDipole[iSymm][k][2];
            final double ukxCR = inducedDipoleCR[iSymm][k][0];
            final double ukyCR = inducedDipoleCR[iSymm][k][1];
            final double ukzCR = inducedDipoleCR[iSymm][k][2];
            final double rbk = born[k];
            final double rb2 = rbi * rbk;
            final double expterm = exp(-r2 / (gkc * rb2));
            final double expc = expterm / gkc;
            final double expc1 = 1.0 - expc;
            final double gf2 = 1.0 / (r2 + rb2 * expterm);
            final double gf = sqrt(gf2);
            final double gf3 = gf2 * gf;
            final double gf5 = gf3 * gf2;

            // Reaction potential auxiliary terms.
            a[1][0] = -gf3;
            a[2][0] = 3.0 * gf5;

            // Reaction potential gradient auxiliary term.
            a[1][1] = expc1 * a[2][0];

            // Multiply the potential auxiliary terms by their dielectric functions.
            a[1][0] = fd * a[1][0];
            a[1][1] = fd * a[1][1];
            a[2][0] = fq * a[2][0];

            // Unweighted reaction potential gradient tensor.
            gux[2] = a[1][0] + xr2 * a[1][1];
            gux[3] = xr * yr * a[1][1];
            gux[4] = xr * zr * a[1][1];
            guy[2] = gux[3];
            guy[3] = a[1][0] + yr2 * a[1][1];
            guy[4] = yr * zr * a[1][1];
            guz[2] = gux[4];
            guz[3] = guy[4];
            guz[4] = a[1][0] + zr2 * a[1][1];

            // Compute the reaction field due to induced dipoles.
            double fix = ukx * gux[2] + uky * guy[2] + ukz * guz[2];
            double fiy = ukx * gux[3] + uky * guy[3] + ukz * guz[3];
            double fiz = ukx * gux[4] + uky * guy[4] + ukz * guz[4];
            double fkx = uix * gux[2] + uiy * guy[2] + uiz * guz[2];
            double fky = uix * gux[3] + uiy * guy[3] + uiz * guz[3];
            double fkz = uix * gux[4] + uiy * guy[4] + uiz * guz[4];
            double fixCR = ukxCR * gux[2] + ukyCR * guy[2] + ukzCR * guz[2];
            double fiyCR = ukxCR * gux[3] + ukyCR * guy[3] + ukzCR * guz[3];
            double fizCR = ukxCR * gux[4] + ukyCR * guy[4] + ukzCR * guz[4];
            double fkxCR = uixCR * gux[2] + uiyCR * guy[2] + uizCR * guz[2];
            double fkyCR = uixCR * gux[3] + uiyCR * guy[3] + uizCR * guz[3];
            double fkzCR = uixCR * gux[4] + uiyCR * guy[4] + uizCR * guz[4];

            // Scale the self-field by half, such that it sums to one below.
            if (i == k) {
                fix *= 0.5;
                fiy *= 0.5;
                fiz *= 0.5;
                fkx *= 0.5;
                fky *= 0.5;
                fkz *= 0.5;
                fixCR *= 0.5;
                fiyCR *= 0.5;
                fizCR *= 0.5;
                fkxCR *= 0.5;
                fkyCR *= 0.5;
                fkzCR *= 0.5;
            }

            double xc = fkx;
            double yc = fky;
            double zc = fkz;
            fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
            fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
            fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
            sharedGKField.add(threadID, i, fix, fiy, fiz);
            sharedGKField.add(threadID, k, fkx, fky, fkz);

            xc = fkxCR;
            yc = fkyCR;
            zc = fkzCR;
            fkxCR = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
            fkyCR = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
            fkzCR = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
            sharedGKFieldCR.add(threadID, i, fixCR, fiyCR, fizCR);
            sharedGKFieldCR.add(threadID, k, fkxCR, fkyCR, fkzCR);
        }
    }
}
