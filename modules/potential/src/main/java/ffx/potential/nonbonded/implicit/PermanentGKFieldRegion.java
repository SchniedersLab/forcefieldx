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
import static ffx.utilities.Constants.dWater;

/**
 * Parallel computation of the Generalized Kirkwood permanent reaction field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PermanentGKFieldRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PermanentGKFieldRegion.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    protected Atom[] atoms;
    /**
     * Multipole moments for each symmetry operator.
     */
    private double[][][] globalMultipole;
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
     * Empirical constant that controls the GK cross-term.
     */
    private static final double gkc = 2.455;
    /**
     * Kirkwood monopole reaction field constant.
     */
    private final double fc;
    /**
     * Kirkwood dipole reaction field constant.
     */
    private final double fd;
    /**
     * Kirkwood quadrupole reaction field constant.
     */
    private final double fq;
    private final PermanentGKFieldLoop[] permanentGKFieldLoop;

    public PermanentGKFieldRegion(int nt, ForceField forceField) {

        // Set the Kirkwood multipolar reaction field constants.
        double epsilon = forceField.getDouble("GK_EPSILON", dWater);
        fc = 1.0 * (1.0 - epsilon) / (0.0 + 1.0 * epsilon);
        fd = 2.0 * (1.0 - epsilon) / (1.0 + 2.0 * epsilon);
        fq = 3.0 * (1.0 - epsilon) / (2.0 + 3.0 * epsilon);

        permanentGKFieldLoop = new PermanentGKFieldLoop[nt];
        for (int i = 0; i < nt; i++) {
            permanentGKFieldLoop[i] = new PermanentGKFieldLoop();
        }
    }

    public void init(Atom[] atoms, double[][][] globalMultipole,
                     Crystal crystal, double[][][] sXYZ, int[][][] neighborLists,
                     boolean[] use, double cut2, double[] born, AtomicDoubleArray3D sharedGKField) {
        this.atoms = atoms;
        this.globalMultipole = globalMultipole;
        this.crystal = crystal;
        this.sXYZ = sXYZ;
        this.neighborLists = neighborLists;
        this.use = use;
        this.cut2 = cut2;
        this.born = born;
        this.sharedGKField = sharedGKField;
    }

    @Override
    public void run() {
        try {
            int threadIndex = getThreadIndex();
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, permanentGKFieldLoop[threadIndex]);
        } catch (Exception e) {
            String message = "Fatal exception computing GK Energy in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute the Generalized Kirkwood permanent reaction field.
     *
     * @since 1.0
     */
    private class PermanentGKFieldLoop extends IntegerForLoop {

        private int threadID;
        private final double[][] a;
        private final double[] gc;
        private final double[] gux, guy, guz;
        private final double[] gqxx, gqyy, gqzz;
        private final double[] gqxy, gqxz, gqyz;
        private final double[] dx_local;
        private double xi, yi, zi;
        private double ci, uxi, uyi, uzi, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi;
        private double rbi;
        private int iSymm;
        private double[][] transOp;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        PermanentGKFieldLoop() {
            a = new double[4][3];
            gc = new double[11];
            gux = new double[11];
            guy = new double[11];
            guz = new double[11];
            gqxx = new double[11];
            gqyy = new double[11];
            gqzz = new double[11];
            gqxy = new double[11];
            gqxz = new double[11];
            gqyz = new double[11];
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
                    double[] multipolei = globalMultipole[0][i];
                    ci = multipolei[t000];
                    uxi = multipolei[t100];
                    uyi = multipolei[t010];
                    uzi = multipolei[t001];
                    qxxi = multipolei[t200] * oneThird;
                    qxyi = multipolei[t110] * oneThird;
                    qxzi = multipolei[t101] * oneThird;
                    qyyi = multipolei[t020] * oneThird;
                    qyzi = multipolei[t011] * oneThird;
                    qzzi = multipolei[t002] * oneThird;
                    rbi = born[i];
                    int[] list = neighborLists[iSymm][i];
                    for (int k : list) {
                        if (!use[k]) {
                            continue;
                        }
                        permanentGKField(i, k);
                    }

                    // Include the self permanent reaction field, which is not in the neighbor list.
                    if (iSymm == 0) {
                        permanentGKField(i, i);
                    }
                }
            }
        }

        private void permanentGKField(int i, int k) {
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
            final double rbk = born[k];
            final double[] multipolek = globalMultipole[iSymm][k];
            final double ck = multipolek[t000];
            final double uxk = multipolek[t100];
            final double uyk = multipolek[t010];
            final double uzk = multipolek[t001];
            final double qxxk = multipolek[t200] * oneThird;
            final double qxyk = multipolek[t110] * oneThird;
            final double qxzk = multipolek[t101] * oneThird;
            final double qyyk = multipolek[t020] * oneThird;
            final double qyzk = multipolek[t011] * oneThird;
            final double qzzk = multipolek[t002] * oneThird;
            final double rb2 = rbi * rbk;
            final double expterm = exp(-r2 / (gkc * rb2));
            final double expc = expterm / gkc;
            final double expc1 = 1.0 - expc;
            final double dexpc = -2.0 / (gkc * rb2);
            final double expcdexpc = -expc * dexpc;
            final double gf2 = 1.0 / (r2 + rb2 * expterm);
            final double gf = sqrt(gf2);
            final double gf3 = gf2 * gf;
            final double gf5 = gf3 * gf2;
            final double gf7 = gf5 * gf2;

            // Reaction potential auxiliary terms.
            a[0][0] = gf;
            a[1][0] = -gf3;
            a[2][0] = 3.0 * gf5;
            a[3][0] = -15.0 * gf7;

            // Reaction potential gradient auxiliary terms.
            a[0][1] = expc1 * a[1][0];
            a[1][1] = expc1 * a[2][0];
            a[2][1] = expc1 * a[3][0];
            // 2nd reaction potential gradient auxiliary terms.
            a[1][2] = expc1 * a[2][1] + expcdexpc * a[2][0];

            // Multiply the potential auxiliary terms by their dielectric functions.
            a[0][1] = fc * a[0][1];
            a[1][0] = fd * a[1][0];
            a[1][1] = fd * a[1][1];
            a[1][2] = fd * a[1][2];
            a[2][0] = fq * a[2][0];
            a[2][1] = fq * a[2][1];

            // Unweighted reaction potential tensor.
            gux[1] = xr * a[1][0];
            guy[1] = yr * a[1][0];
            guz[1] = zr * a[1][0];

            // Unweighted reaction potential gradient tensor.
            gc[2] = xr * a[0][1];
            gc[3] = yr * a[0][1];
            gc[4] = zr * a[0][1];
            gux[2] = a[1][0] + xr2 * a[1][1];
            gux[3] = xr * yr * a[1][1];
            gux[4] = xr * zr * a[1][1];
            guy[2] = gux[3];
            guy[3] = a[1][0] + yr2 * a[1][1];
            guy[4] = yr * zr * a[1][1];
            guz[2] = gux[4];
            guz[3] = guy[4];
            guz[4] = a[1][0] + zr2 * a[1][1];
            gqxx[2] = xr * (2.0 * a[2][0] + xr2 * a[2][1]);
            gqxx[3] = yr * xr2 * a[2][1];
            gqxx[4] = zr * xr2 * a[2][1];
            gqyy[2] = xr * yr2 * a[2][1];
            gqyy[3] = yr * (2.0 * a[2][0] + yr2 * a[2][1]);
            gqyy[4] = zr * yr2 * a[2][1];
            gqzz[2] = xr * zr2 * a[2][1];
            gqzz[3] = yr * zr2 * a[2][1];
            gqzz[4] = zr * (2.0 * a[2][0] + zr2 * a[2][1]);
            gqxy[2] = yr * (a[2][0] + xr2 * a[2][1]);
            gqxy[3] = xr * (a[2][0] + yr2 * a[2][1]);
            gqxy[4] = zr * xr * yr * a[2][1];
            gqxz[2] = zr * (a[2][0] + xr2 * a[2][1]);
            gqxz[3] = gqxy[4];
            gqxz[4] = xr * (a[2][0] + zr2 * a[2][1]);
            gqyz[2] = gqxy[4];
            gqyz[3] = zr * (a[2][0] + yr2 * a[2][1]);
            gqyz[4] = yr * (a[2][0] + zr2 * a[2][1]);

            // Unweighted 2nd reaction potential gradient tensor.
            gux[5] = xr * (3.0 * a[1][1] + xr2 * a[1][2]);
            gux[6] = yr * (a[1][1] + xr2 * a[1][2]);
            gux[7] = zr * (a[1][1] + xr2 * a[1][2]);
            gux[8] = xr * (a[1][1] + yr2 * a[1][2]);
            gux[9] = zr * xr * yr * a[1][2];
            gux[10] = xr * (a[1][1] + zr2 * a[1][2]);
            guy[5] = yr * (a[1][1] + xr2 * a[1][2]);
            guy[6] = xr * (a[1][1] + yr2 * a[1][2]);
            guy[7] = gux[9];
            guy[8] = yr * (3.0 * a[1][1] + yr2 * a[1][2]);
            guy[9] = zr * (a[1][1] + yr2 * a[1][2]);
            guy[10] = yr * (a[1][1] + zr2 * a[1][2]);
            guz[5] = zr * (a[1][1] + xr2 * a[1][2]);
            guz[6] = gux[9];
            guz[7] = xr * (a[1][1] + zr2 * a[1][2]);
            guz[8] = zr * (a[1][1] + yr2 * a[1][2]);
            guz[9] = yr * (a[1][1] + zr2 * a[1][2]);
            guz[10] = zr * (3.0 * a[1][1] + zr2 * a[1][2]);

            // Generalized Kirkwood permanent reaction field.
            double fix = uxk * gux[2] + uyk * gux[3] + uzk * gux[4]
                    + 0.5 * (ck * gux[1] + qxxk * gux[5]
                    + qyyk * gux[8] + qzzk * gux[10]
                    + 2.0 * (qxyk * gux[6] + qxzk * gux[7]
                    + qyzk * gux[9]))
                    + 0.5 * (ck * gc[2] + qxxk * gqxx[2]
                    + qyyk * gqyy[2] + qzzk * gqzz[2]
                    + 2.0 * (qxyk * gqxy[2] + qxzk * gqxz[2]
                    + qyzk * gqyz[2]));
            double fiy = uxk * guy[2] + uyk * guy[3] + uzk * guy[4]
                    + 0.5 * (ck * guy[1] + qxxk * guy[5]
                    + qyyk * guy[8] + qzzk * guy[10]
                    + 2.0 * (qxyk * guy[6] + qxzk * guy[7]
                    + qyzk * guy[9]))
                    + 0.5 * (ck * gc[3] + qxxk * gqxx[3]
                    + qyyk * gqyy[3] + qzzk * gqzz[3]
                    + 2.0 * (qxyk * gqxy[3] + qxzk * gqxz[3]
                    + qyzk * gqyz[3]));
            double fiz = uxk * guz[2] + uyk * guz[3] + uzk * guz[4]
                    + 0.5 * (ck * guz[1] + qxxk * guz[5]
                    + qyyk * guz[8] + qzzk * guz[10]
                    + 2.0 * (qxyk * guz[6] + qxzk * guz[7]
                    + qyzk * guz[9]))
                    + 0.5 * (ck * gc[4] + qxxk * gqxx[4]
                    + qyyk * gqyy[4] + qzzk * gqzz[4]
                    + 2.0 * (qxyk * gqxy[4] + qxzk * gqxz[4]
                    + qyzk * gqyz[4]));
            double fkx = uxi * gux[2] + uyi * gux[3] + uzi * gux[4]
                    - 0.5 * (ci * gux[1] + qxxi * gux[5]
                    + qyyi * gux[8] + qzzi * gux[10]
                    + 2.0 * (qxyi * gux[6] + qxzi * gux[7]
                    + qyzi * gux[9]))
                    - 0.5 * (ci * gc[2] + qxxi * gqxx[2]
                    + qyyi * gqyy[2] + qzzi * gqzz[2]
                    + 2.0 * (qxyi * gqxy[2] + qxzi * gqxz[2]
                    + qyzi * gqyz[2]));
            double fky = uxi * guy[2] + uyi * guy[3] + uzi * guy[4]
                    - 0.5 * (ci * guy[1] + qxxi * guy[5]
                    + qyyi * guy[8] + qzzi * guy[10]
                    + 2.0 * (qxyi * guy[6] + qxzi * guy[7]
                    + qyzi * guy[9]))
                    - 0.5 * (ci * gc[3] + qxxi * gqxx[3]
                    + qyyi * gqyy[3] + qzzi * gqzz[3]
                    + 2.0 * (qxyi * gqxy[3] + qxzi * gqxz[3]
                    + qyzi * gqyz[3]));
            double fkz = uxi * guz[2] + uyi * guz[3] + uzi * guz[4]
                    - 0.5 * (ci * guz[1] + qxxi * guz[5]
                    + qyyi * guz[8] + qzzi * guz[10]
                    + 2.0 * (qxyi * guz[6] + qxzi * guz[7]
                    + qyzi * guz[9]))
                    - 0.5 * (ci * gc[4] + qxxi * gqxx[4]
                    + qyyi * gqyy[4] + qzzi * gqzz[4]
                    + 2.0 * (qxyi * gqxy[4] + qxzi * gqxz[4]
                    + qyzi * gqyz[4]));

            // Scale the self-field by half, such that it sums to one below.
            if (i == k) {
                fix *= 0.5;
                fiy *= 0.5;
                fiz *= 0.5;
                fkx *= 0.5;
                fky *= 0.5;
                fkz *= 0.5;
            }

            double xc = fkx;
            double yc = fky;
            double zc = fkz;
            fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
            fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
            fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];

            sharedGKField.add(threadID, i, fix, fiy, fiz);
            sharedGKField.add(threadID, k, fkx, fky, fkz);
        }
    }

    /**
     * Constant factor used with quadrupoles.
     */
    private static final double oneThird = 1.0 / 3.0;
}
