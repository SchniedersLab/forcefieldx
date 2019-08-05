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
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.PCGVariables;
import static ffx.numerics.special.Erf.erfc;

/**
 * Evaluate the real space field due to induced dipoles using a short cutoff
 * (~3-4 A).
 */
public class InducedDipolePreconditionerRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(InducedDipolePreconditionerRegion.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    private double[][][] coordinates;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public double[][][] inducedDipole;
    public double[][][] inducedDipoleCR;
    /**
     * Unit cell and spacegroup information.
     */
    private Crystal crystal;
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
    private double[] ipdamp;
    private double[] thole;

    /**
     * Neighbor lists, without atoms beyond the preconditioner cutoff.
     * [nSymm][nAtoms][nIncludedNeighbors]
     */
    private int[][][] preconditionerLists;
    /**
     * Number of neighboring atoms within the preconditioner cutoff.
     * [nSymm][nAtoms]
     */
    private int[][] preconditionerCounts;
    /**
     * Pairwise schedule for load balancing.
     */
    private IntegerSchedule realSpaceSchedule;
    /**
     * Field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] field;
    /**
     * Chain rule field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] fieldCR;
    private long[] realSpaceSCFTime;
    private double aewald;
    private double an0, an1;

    private final int maxThreads;
    private final InducedPreconditionerFieldLoop[] inducedPreconditionerFieldLoop;
    private final ReduceLoop[] reduceLoop;

    public InducedDipolePreconditionerRegion(int threadCount) {
        maxThreads = threadCount;
        inducedPreconditionerFieldLoop = new InducedPreconditionerFieldLoop[threadCount];
        reduceLoop = new ReduceLoop[threadCount];
    }

    public void init(Atom[] atoms, double[][][] coordinates,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     Crystal crystal, boolean[] use, double[] ipdamp, double[] thole,
                     PCGVariables pcgVariables, IntegerSchedule realSpaceSchedule,
                     double[][][] field, double[][][] fieldCR, long[] realSpaceSCFTime,
                     EwaldParameters ewaldParameters) {
        this.atoms = atoms;
        this.coordinates = coordinates;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.crystal = crystal;
        this.use = use;
        this.ipdamp = ipdamp;
        this.thole = thole;
        this.preconditionerLists = pcgVariables.preconditionerLists;
        this.preconditionerCounts = pcgVariables.preconditionerCounts;
        this.realSpaceSchedule = realSpaceSchedule;
        this.field = field;
        this.fieldCR = fieldCR;
        this.realSpaceSCFTime = realSpaceSCFTime;
        this.aewald = ewaldParameters.aewald;
        this.an0 = ewaldParameters.an0;
        this.an1 = ewaldParameters.an1;
    }

    @Override
    public void run() {
        int threadIndex = getThreadIndex();
        if (inducedPreconditionerFieldLoop[threadIndex] == null) {
            inducedPreconditionerFieldLoop[threadIndex] = new InducedPreconditionerFieldLoop();
            reduceLoop[threadIndex] = new ReduceLoop();
        }
        try {
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, inducedPreconditionerFieldLoop[threadIndex]);
            execute(0, nAtoms - 1, reduceLoop[threadIndex]);
        } catch (Exception e) {
            String message = "Fatal exception computing the induced real space field in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private class InducedPreconditionerFieldLoop extends IntegerForLoop {

        private double[] x, y, z;
        private double[][] ind, indCR;
        private double[] fX, fY, fZ;
        private double[] fXCR, fYCR, fZCR;

        InducedPreconditionerFieldLoop() {
        }

        @Override
        public IntegerSchedule schedule() {
            return realSpaceSchedule;
        }

        @Override
        public void start() {
            int threadIndex = getThreadIndex();
            realSpaceSCFTime[threadIndex] -= System.nanoTime();
            fX = field[threadIndex][0];
            fY = field[threadIndex][1];
            fZ = field[threadIndex][2];
            fXCR = fieldCR[threadIndex][0];
            fYCR = fieldCR[threadIndex][1];
            fZCR = fieldCR[threadIndex][2];
            fill(fX, 0.0);
            fill(fY, 0.0);
            fill(fZ, 0.0);
            fill(fXCR, 0.0);
            fill(fYCR, 0.0);
            fill(fZCR, 0.0);
            x = coordinates[0][0];
            y = coordinates[0][1];
            z = coordinates[0][2];
            ind = inducedDipole[0];
            indCR = inducedDipoleCR[0];
        }

        @Override
        public void finish() {
            int threadIndex = getThreadIndex();
            realSpaceSCFTime[threadIndex] += System.nanoTime();
        }

        @Override
        public void run(int lb, int ub) {
            final double[] dx = new double[3];
            final double[][] transOp = new double[3][3];

            // Loop over a chunk of atoms.
            int[][] lists = preconditionerLists[0];
            int[] counts = preconditionerCounts[0];
            for (int i = lb; i <= ub; i++) {
                if (!use[i]) {
                    continue;
                }
                double fx = 0.0;
                double fy = 0.0;
                double fz = 0.0;
                double px = 0.0;
                double py = 0.0;
                double pz = 0.0;
                final double xi = x[i];
                final double yi = y[i];
                final double zi = z[i];
                final double[] dipolei = ind[i];
                final double uix = dipolei[0];
                final double uiy = dipolei[1];
                final double uiz = dipolei[2];
                final double[] dipoleCRi = indCR[i];
                final double pix = dipoleCRi[0];
                final double piy = dipoleCRi[1];
                final double piz = dipoleCRi[2];
                final double pdi = ipdamp[i];
                final double pti = thole[i];

                // Loop over the neighbor list.
                final int[] list = lists[i];
                final int npair = counts[i];
                for (int j = 0; j < npair; j++) {
                    final int k = list[j];
                    if (!use[k]) {
                        continue;
                    }
                    final double pdk = ipdamp[k];
                    final double ptk = thole[k];
                    dx[0] = x[k] - xi;
                    dx[1] = y[k] - yi;
                    dx[2] = z[k] - zi;
                    final double r2 = crystal.image(dx);

                    // Calculate the error function damping terms.
                    final double r = sqrt(r2);
                    final double rr1 = 1.0 / r;
                    final double rr2 = rr1 * rr1;
                    final double ralpha = aewald * r;
                    final double exp2a = exp(-ralpha * ralpha);
                    final double bn0 = erfc(ralpha) * rr1;
                    // final double exp2a = 1.0;
                    // final double bn0 = rr1;
                    final double bn1 = (bn0 + an0 * exp2a) * rr2;
                    final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                    double scale3 = 1.0;
                    double scale5 = 1.0;
                    double damp = pdi * pdk;
                    final double pgamma = min(pti, ptk);
                    final double rdamp = r * damp;
                    damp = -pgamma * rdamp * rdamp * rdamp;
                    if (damp > -50.0) {
                        final double expdamp = exp(damp);
                        scale3 = 1.0 - expdamp;
                        scale5 = 1.0 - expdamp * (1.0 - damp);
                    }
                    double rr3 = rr1 * rr2;
                    double rr5 = 3.0 * rr3 * rr2;
                    rr3 *= (1.0 - scale3);
                    rr5 *= (1.0 - scale5);
                    final double xr = dx[0];
                    final double yr = dx[1];
                    final double zr = dx[2];
                    final double[] dipolek = ind[k];
                    final double ukx = dipolek[0];
                    final double uky = dipolek[1];
                    final double ukz = dipolek[2];
                    final double ukr = ukx * xr + uky * yr + ukz * zr;
                    final double bn2ukr = bn2 * ukr;
                    final double fimx = -bn1 * ukx + bn2ukr * xr;
                    final double fimy = -bn1 * uky + bn2ukr * yr;
                    final double fimz = -bn1 * ukz + bn2ukr * zr;
                    final double rr5ukr = rr5 * ukr;
                    final double fidx = -rr3 * ukx + rr5ukr * xr;
                    final double fidy = -rr3 * uky + rr5ukr * yr;
                    final double fidz = -rr3 * ukz + rr5ukr * zr;
                    fx += (fimx - fidx);
                    fy += (fimy - fidy);
                    fz += (fimz - fidz);
                    final double[] dipolepk = indCR[k];
                    final double pkx = dipolepk[0];
                    final double pky = dipolepk[1];
                    final double pkz = dipolepk[2];
                    final double pkr = pkx * xr + pky * yr + pkz * zr;
                    final double bn2pkr = bn2 * pkr;
                    final double pimx = -bn1 * pkx + bn2pkr * xr;
                    final double pimy = -bn1 * pky + bn2pkr * yr;
                    final double pimz = -bn1 * pkz + bn2pkr * zr;
                    final double rr5pkr = rr5 * pkr;
                    final double pidx = -rr3 * pkx + rr5pkr * xr;
                    final double pidy = -rr3 * pky + rr5pkr * yr;
                    final double pidz = -rr3 * pkz + rr5pkr * zr;
                    px += (pimx - pidx);
                    py += (pimy - pidy);
                    pz += (pimz - pidz);
                    final double uir = uix * xr + uiy * yr + uiz * zr;
                    final double bn2uir = bn2 * uir;
                    final double fkmx = -bn1 * uix + bn2uir * xr;
                    final double fkmy = -bn1 * uiy + bn2uir * yr;
                    final double fkmz = -bn1 * uiz + bn2uir * zr;
                    final double rr5uir = rr5 * uir;
                    final double fkdx = -rr3 * uix + rr5uir * xr;
                    final double fkdy = -rr3 * uiy + rr5uir * yr;
                    final double fkdz = -rr3 * uiz + rr5uir * zr;
                    fX[k] += (fkmx - fkdx);
                    fY[k] += (fkmy - fkdy);
                    fZ[k] += (fkmz - fkdz);
                    final double pir = pix * xr + piy * yr + piz * zr;
                    final double bn2pir = bn2 * pir;
                    final double pkmx = -bn1 * pix + bn2pir * xr;
                    final double pkmy = -bn1 * piy + bn2pir * yr;
                    final double pkmz = -bn1 * piz + bn2pir * zr;
                    final double rr5pir = rr5 * pir;
                    final double pkdx = -rr3 * pix + rr5pir * xr;
                    final double pkdy = -rr3 * piy + rr5pir * yr;
                    final double pkdz = -rr3 * piz + rr5pir * zr;
                    fXCR[k] += (pkmx - pkdx);
                    fYCR[k] += (pkmy - pkdy);
                    fZCR[k] += (pkmz - pkdz);
                }
                fX[i] += fx;
                fY[i] += fy;
                fZ[i] += fz;
                fXCR[i] += px;
                fYCR[i] += py;
                fZCR[i] += pz;
            }

            // Loop over symmetry mates.
            List<SymOp> symOps = crystal.spaceGroup.symOps;
            int nSymm = symOps.size();
            for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                crystal.getTransformationOperator(symOp, transOp);
                lists = preconditionerLists[iSymm];
                counts = preconditionerCounts[iSymm];
                final double[] xs = coordinates[iSymm][0];
                final double[] ys = coordinates[iSymm][1];
                final double[] zs = coordinates[iSymm][2];
                final double[][] inds = inducedDipole[iSymm];
                final double[][] indCRs = inducedDipoleCR[iSymm];

                // Loop over a chunk of atoms.
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double px = 0.0;
                    double py = 0.0;
                    double pz = 0.0;
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double[] dipolei = ind[i];
                    final double uix = dipolei[0];
                    final double uiy = dipolei[1];
                    final double uiz = dipolei[2];
                    final double[] dipoleCRi = indCR[i];
                    final double pix = dipoleCRi[0];
                    final double piy = dipoleCRi[1];
                    final double piz = dipoleCRi[2];
                    final double pdi = ipdamp[i];
                    final double pti = thole[i];

                    // Loop over the neighbor list.
                    final int[] list = lists[i];
                    final int npair = counts[i];
                    for (int j = 0; j < npair; j++) {
                        final int k = list[j];
                        if (!use[k]) {
                            continue;
                        }
                        double selfScale = 1.0;
                        if (i == k) {
                            selfScale = 0.5;
                        }
                        final double pdk = ipdamp[k];
                        final double ptk = thole[k];
                        dx[0] = xs[k] - xi;
                        dx[1] = ys[k] - yi;
                        dx[2] = zs[k] - zi;
                        final double r2 = crystal.image(dx);

                        // Calculate the error function damping terms.
                        final double r = sqrt(r2);
                        final double rr1 = 1.0 / r;
                        final double rr2 = rr1 * rr1;
                        final double ralpha = aewald * r;
                        final double exp2a = exp(-ralpha * ralpha);
                        final double bn0 = erfc(ralpha) * rr1;
                        //final double exp2a = 1.0;
                        //final double bn0 = rr1;
                        final double bn1 = (bn0 + an0 * exp2a) * rr2;
                        final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                        double scale3 = 1.0;
                        double scale5 = 1.0;
                        double damp = pdi * pdk;
                        final double pgamma = min(pti, ptk);
                        final double rdamp = r * damp;
                        damp = -pgamma * rdamp * rdamp * rdamp;
                        if (damp > -50.0) {
                            final double expdamp = exp(damp);
                            scale3 = 1.0 - expdamp;
                            scale5 = 1.0 - expdamp * (1.0 - damp);
                        }
                        double rr3 = rr1 * rr2;
                        double rr5 = 3.0 * rr3 * rr2;
                        rr3 *= (1.0 - scale3);
                        rr5 *= (1.0 - scale5);
                        final double xr = dx[0];
                        final double yr = dx[1];
                        final double zr = dx[2];
                        final double[] dipolek = inds[k];
                        final double ukx = dipolek[0];
                        final double uky = dipolek[1];
                        final double ukz = dipolek[2];
                        final double[] dipolepk = indCRs[k];
                        final double pkx = dipolepk[0];
                        final double pky = dipolepk[1];
                        final double pkz = dipolepk[2];
                        final double ukr = ukx * xr + uky * yr + ukz * zr;
                        final double bn2ukr = bn2 * ukr;
                        final double fimx = -bn1 * ukx + bn2ukr * xr;
                        final double fimy = -bn1 * uky + bn2ukr * yr;
                        final double fimz = -bn1 * ukz + bn2ukr * zr;
                        final double rr5ukr = rr5 * ukr;
                        final double fidx = -rr3 * ukx + rr5ukr * xr;
                        final double fidy = -rr3 * uky + rr5ukr * yr;
                        final double fidz = -rr3 * ukz + rr5ukr * zr;
                        fx += selfScale * (fimx - fidx);
                        fy += selfScale * (fimy - fidy);
                        fz += selfScale * (fimz - fidz);
                        final double pkr = pkx * xr + pky * yr + pkz * zr;
                        final double bn2pkr = bn2 * pkr;
                        final double pimx = -bn1 * pkx + bn2pkr * xr;
                        final double pimy = -bn1 * pky + bn2pkr * yr;
                        final double pimz = -bn1 * pkz + bn2pkr * zr;
                        final double rr5pkr = rr5 * pkr;
                        final double pidx = -rr3 * pkx + rr5pkr * xr;
                        final double pidy = -rr3 * pky + rr5pkr * yr;
                        final double pidz = -rr3 * pkz + rr5pkr * zr;
                        px += selfScale * (pimx - pidx);
                        py += selfScale * (pimy - pidy);
                        pz += selfScale * (pimz - pidz);
                        final double uir = uix * xr + uiy * yr + uiz * zr;
                        final double bn2uir = bn2 * uir;
                        final double fkmx = -bn1 * uix + bn2uir * xr;
                        final double fkmy = -bn1 * uiy + bn2uir * yr;
                        final double fkmz = -bn1 * uiz + bn2uir * zr;
                        final double rr5uir = rr5 * uir;
                        final double fkdx = -rr3 * uix + rr5uir * xr;
                        final double fkdy = -rr3 * uiy + rr5uir * yr;
                        final double fkdz = -rr3 * uiz + rr5uir * zr;
                        double xc = selfScale * (fkmx - fkdx);
                        double yc = selfScale * (fkmy - fkdy);
                        double zc = selfScale * (fkmz - fkdz);
                        fX[k] += (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                        fY[k] += (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                        fZ[k] += (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                        final double pir = pix * xr + piy * yr + piz * zr;
                        final double bn2pir = bn2 * pir;
                        final double pkmx = -bn1 * pix + bn2pir * xr;
                        final double pkmy = -bn1 * piy + bn2pir * yr;
                        final double pkmz = -bn1 * piz + bn2pir * zr;
                        final double rr5pir = rr5 * pir;
                        final double pkdx = -rr3 * pix + rr5pir * xr;
                        final double pkdy = -rr3 * piy + rr5pir * yr;
                        final double pkdz = -rr3 * piz + rr5pir * zr;
                        xc = selfScale * (pkmx - pkdx);
                        yc = selfScale * (pkmy - pkdy);
                        zc = selfScale * (pkmz - pkdz);
                        fXCR[k] += (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                        fYCR[k] += (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                        fZCR[k] += (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                    }
                    fX[i] += fx;
                    fY[i] += fy;
                    fZ[i] += fz;
                    fXCR[i] += px;
                    fYCR[i] += py;
                    fZCR[i] += pz;
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
            // Reduce the real space field.
            for (int i = lb; i <= ub; i++) {
                double fx = 0.0;
                double fy = 0.0;
                double fz = 0.0;
                double fxCR = 0.0;
                double fyCR = 0.0;
                double fzCR = 0.0;
                for (int j = 1; j < maxThreads; j++) {
                    fx += field[j][0][i];
                    fy += field[j][1][i];
                    fz += field[j][2][i];
                    fxCR += fieldCR[j][0][i];
                    fyCR += fieldCR[j][1][i];
                    fzCR += fieldCR[j][2][i];
                }
                field[0][0][i] += fx;
                field[0][1][i] += fy;
                field[0][2][i] += fz;
                fieldCR[0][0][i] += fxCR;
                fieldCR[0][1][i] += fyCR;
                fieldCR[0][2][i] += fzCR;
            }
        }
    }
}
