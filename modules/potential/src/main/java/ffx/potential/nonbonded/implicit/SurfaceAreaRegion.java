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
import static java.util.Arrays.fill;
import static java.util.Arrays.sort;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedBooleanArray;
import edu.rit.pj.reduction.SharedDouble;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.VDWType;
import ffx.potential.utils.EnergyException;

/**
 * SurfaceAreaRegion performs an analytical computation of the weighted
 * solvent accessible surface area of each atom and the first
 * derivatives of the area with respect to Cartesian coordinates
 * <p>
 * Literature references:
 * <p>
 * T. J. Richmond, "Solvent Accessible Surface Area and
 * Excluded Volume in Proteins", Journal of Molecular Biology,
 * 178, 63-89 (1984)
 * <p>
 * L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
 * Applied to Molecular Dynamics of Proteins in Solution",
 * Protein Science, 1, 227-235 (1992)
 * <p>
 * This was ported from the TINKER "surface.f" code.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class SurfaceAreaRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(SurfaceAreaRegion.class.getName());

    private final AtomOverlapLoop[] atomOverlapLoop;
    private final CavitationLoop[] cavitationLoop;
    private final InitLoop[] initLoop;
    private final SharedDouble sharedCavitation;
    private final int[][][] neighborLists;
    private final Atom[] atoms;
    private final double[] x, y, z;
    private final boolean[] use;

    private final int nAtoms;
    private final double probe;
    private final int maxarc = 150;
    private final static double delta = 1.0e-8;
    private final static double delta2 = delta * delta;

    private double[][] xc1;
    private double[][] yc1;
    private double[][] zc1;
    private double[][] dsq1;
    private double[][] bsq1;
    private double[][] b1;
    private IndexedDouble[][] gr;
    private int[][] intag1;
    private Integer[] count;
    private boolean[] buried;
    private SharedBooleanArray skip;
    private double[] area;
    private double[] r;
    private AtomicDoubleArray3D grad;
    private double surfaceTension;

    /**
     * This class is a port of the Cavitation code in TINKER.
     * <p>
     * GK implicit solvents are moving to use Gaussian based definitions of surface area for efficiency.
     *
     * @param atoms          Atom array.
     * @param x              X-coordinate array.
     * @param y              Y-coordinate array.
     * @param z              Z-coordinate array.
     * @param use            Specifies if the atom is used in the potential.
     * @param neighborLists  Neighbor-list array.
     * @param grad           Gradient array.
     * @param nt             Number of threads.
     * @param probe          Solvent probe radius.
     * @param surfaceTension Surface tension.
     */
    public SurfaceAreaRegion(Atom[] atoms, double[] x, double[] y, double[] z,
                             boolean[] use, int[][][] neighborLists,
                             AtomicDoubleArray3D grad,
                             int nt, double probe, double surfaceTension) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.x = x;
        this.y = y;
        this.z = z;
        this.use = use;
        this.neighborLists = neighborLists;
        this.probe = probe;
        this.grad = grad;
        this.surfaceTension = surfaceTension;

        atomOverlapLoop = new AtomOverlapLoop[nt];
        cavitationLoop = new CavitationLoop[nt];
        initLoop = new InitLoop[nt];
        for (int i = 0; i < nt; i++) {
            atomOverlapLoop[i] = new AtomOverlapLoop();
            cavitationLoop[i] = new CavitationLoop();
            initLoop[i] = new InitLoop();
        }
        sharedCavitation = new SharedDouble();
        init();
    }

    public double getEnergy() {
        return sharedCavitation.get();
    }

    public final void init() {
        if (count == null || count.length < nAtoms) {
            count = new Integer[nAtoms];
            xc1 = new double[nAtoms][maxarc];
            yc1 = new double[nAtoms][maxarc];
            zc1 = new double[nAtoms][maxarc];
            dsq1 = new double[nAtoms][maxarc];
            bsq1 = new double[nAtoms][maxarc];
            b1 = new double[nAtoms][maxarc];
            gr = new IndexedDouble[nAtoms][maxarc];
            intag1 = new int[nAtoms][maxarc];
            buried = new boolean[nAtoms];
            skip = new SharedBooleanArray(nAtoms);
            area = new double[nAtoms];
            r = new double[nAtoms];
        }

        // Set the sphere radii.
        for (int i = 0; i < nAtoms; i++) {
            VDWType type = atoms[i].getVDWType();
            double rmini = type.radius;
            r[i] = rmini / 2.0;
            if (r[i] != 0.0) {
                r[i] = r[i] + probe;
            }
            skip.set(i, true);
        }
    }

    @Override
    public void start() {
        sharedCavitation.set(0.0);
    }

    @Override
    public void finish() {
        if (logger.isLoggable(Level.FINE)) {
            int n = initLoop.length;
            long initTime = 0;
            long overlapTime = 0;
            long cavTime = 0;
            for (int i = 0; i < n; i++) {
                initTime = max(initLoop[i].time, initTime);
                overlapTime = max(atomOverlapLoop[i].time, overlapTime);
                cavTime = max(cavitationLoop[i].time, cavTime);
            }
            logger.fine(format(" Cavitation Init: %10.3f Overlap: %10.3f Cav:  %10.3f",
                    initTime * 1e-9, overlapTime * 1e-9, cavTime * 1e-9));
        }
    }

    @Override
    public void run() {
        try {
            execute(0, nAtoms - 1, initLoop[getThreadIndex()]);
            execute(0, nAtoms - 1, atomOverlapLoop[getThreadIndex()]);
            execute(0, nAtoms - 1, cavitationLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception computing Cavitation energy in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    static private class IndexedDouble implements Comparable {

        public double value;
        public int key;

        IndexedDouble(double value, int key) {
            this.value = value;
            this.key = key;
        }

        @Override
        public int compareTo(Object o) {
            if (!(o instanceof IndexedDouble)) {
                return 0;
            }
            IndexedDouble d = (IndexedDouble) o;
            return Double.compare(value, d.value);
        }
    }

    /**
     * Initialize arrays for Cavitation calculation.
     */
    private class InitLoop extends IntegerForLoop {

        public long time;

        @Override
        public void start() {
            time = -System.nanoTime();
        }

        @Override
        public void finish() {
            time += System.nanoTime();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            // Set the "skip" array to exclude all inactive atoms
            // that do not overlap any of the current active atoms
            for (int i = lb; i <= ub; i++) {
                buried[i] = false;
                area[i] = 0.0;
                count[i] = 0;
                double xr = x[i];
                double yr = y[i];
                double zr = z[i];
                double rri = r[i];
                final int[] list = neighborLists[0][i];
                for (int k : list) {
                    double rrik = rri + r[k];
                    double dx = x[k] - xr;
                    double dy = y[k] - yr;
                    double dz = z[k] - zr;
                    double ccsq = dx * dx + dy * dy + dz * dz;
                    if (ccsq <= rrik * rrik) {
                        if (use[i] || use[k]) {
                            skip.set(k, false);
                            skip.set(i, false);
                        }
                    }
                }
            }
        }
    }

    /**
     * Compute Cavitation energy for a range of atoms.
     *
     * @since 1.0
     */
    private class AtomOverlapLoop extends IntegerForLoop {

        public long time;

        @Override
        public void start() {
            time = -System.nanoTime();
        }

        @Override
        public void finish() {
            time += System.nanoTime();
        }

        @Override
        public void run(int lb, int ub) {

            // Find overlaps with the current sphere.
            for (int i = lb; i <= ub; i++) {
                if (skip.get(i) || !use[i]) {
                    continue;
                }
                int[] list = neighborLists[0][i];
                for (int k : list) {
                    if (k == i) {
                        continue;
                    }
                    pair(i, k);
                    if (!skip.get(k) || !use[k]) {
                        pair(k, i);
                    }
                }
            }
        }

        private void pair(int i, int k) {
            double xi = x[i];
            double yi = y[i];
            double zi = z[i];
            double rri = r[i];
            double rri2 = 2.0 * rri;
            double rplus = rri + r[k];
            double dx = x[k] - xi;
            double dy = y[k] - yi;
            double dz = z[k] - zi;
            if (abs(dx) >= rplus || abs(dy) >= rplus || abs(dz) >= rplus) {
                return;
            }

            // Check for overlap of spheres by testing center to center
            // distance against sum and difference of radii.
            double xysq = dx * dx + dy * dy;
            if (xysq < delta2) {
                dx = delta;
                dy = 0.0;
                xysq = delta2;
            }
            double r2 = xysq + dz * dz;
            double dr = sqrt(r2);
            if (rplus - dr <= delta) {
                return;
            }
            double rminus = rri - r[k];

            // Calculate overlap parameters between "i" and "ir" sphere.
            synchronized (count) {
                // Check for a completely buried "ir" sphere.
                if (dr - abs(rminus) <= delta) {
                    if (rminus <= 0.0) {
                        // SA for this atom is zero.
                        buried[i] = true;
                        return;
                    }
                    return;
                }
                int n = count[i];
                xc1[i][n] = dx;
                yc1[i][n] = dy;
                zc1[i][n] = dz;
                dsq1[i][n] = xysq;
                bsq1[i][n] = r2;
                b1[i][n] = dr;
                gr[i][n] = new IndexedDouble((r2 + rplus * rminus) / (rri2 * b1[i][n]), n);
                intag1[i][n] = k;
                count[i]++;
                if (count[i] >= maxarc) {
                    throw new EnergyException(format(" Increase the value of MAXARC to (%d).", count[i]), false);
                }
            }
        }
    }

    /**
     * Compute Cavitation energy for a range of atoms.
     *
     * @since 1.0
     */
    private class CavitationLoop extends IntegerForLoop {

        private double thec = 0;
        private double localSurfaceEnergy;
        private IndexedDouble[] arci;
        private boolean[] omit;
        private double[] xc;
        private double[] yc;
        private double[] zc;
        private double[] dsq;
        private double[] b;
        private double[] bsq;
        private double[] bg;
        private double[] risq;
        private double[] ri;
        private double[] ther;
        private double[] ider;
        private double[] sign_yder;
        private double[] arcf;
        private double[] ex;
        private double[] ux;
        private double[] uy;
        private double[] uz;
        private int[] kent;
        private int[] kout;
        private int[] intag;
        private int[] lt;
        private int i;
        private int j;
        private int ib;
        private int threadID;

        // Set pi multiples, overlap criterion and tolerances.
        private final static double pix2 = 2.0 * PI;
        private final static double pix4 = 4.0 * PI;
        private final static double pid2 = PI / 2.0;
        private final static double eps = 1.0e-8;
        public long time;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        CavitationLoop() {
            allocateMemory(maxarc);
        }

        private void allocateMemory(int maxarc) {
            arci = new IndexedDouble[maxarc];
            arcf = new double[maxarc];
            risq = new double[maxarc];
            ri = new double[maxarc];
            dsq = new double[maxarc];
            bsq = new double[maxarc];
            intag = new int[maxarc];
            lt = new int[maxarc];
            kent = new int[maxarc];
            kout = new int[maxarc];
            ider = new double[maxarc];
            sign_yder = new double[maxarc];
            xc = new double[maxarc];
            yc = new double[maxarc];
            zc = new double[maxarc];
            b = new double[maxarc];
            omit = new boolean[maxarc];
            bg = new double[maxarc];
            ther = new double[maxarc];
            ex = new double[maxarc];
            ux = new double[maxarc];
            uy = new double[maxarc];
            uz = new double[maxarc];
        }

        @Override
        public void start() {
            time = -System.nanoTime();
            threadID = getThreadIndex();
            localSurfaceEnergy = 0.0;
            fill(ider, 0);
            fill(sign_yder, 0);
        }

        @Override
        public void finish() {
            sharedCavitation.addAndGet(localSurfaceEnergy);
            time += System.nanoTime();
        }

        @Override
        public void run(int lb, int ub) {
            // Compute the area and derivatives of current "ir" sphere
            for (int ir = lb; ir <= ub; ir++) {
                if (skip.get(ir) || !use[ir]) {
                    continue;
                }
                double rri = r[ir];
                double rri2 = 2.0 * rri;
                double rrisq = rri * rri;
                double wght = surfaceTension;
                boolean moved = false;
                surface(rri, rri2, rrisq, wght, moved, ir);
                if (area[ir] < 0.0) {
                    logger.log(Level.WARNING, format(" Negative surface area set to 0 for atom %d.", ir));
                    area[ir] = 0.0;
                }
                area[ir] *= rrisq * wght;
                localSurfaceEnergy += area[ir];
            }
        }

        public void surface(double rri, double rri2, double rrisq, double wght, boolean moved, int ir) {

            ib = 0;
            int jb = 0;
            double arclen = 0.0;
            double exang = 0.0;

            // Case where no other spheres overlap the current sphere.
            if (count[ir] == 0) {
                area[ir] = pix4;
                return;
            }
            // Case where only one sphere overlaps the current sphere.
            if (count[ir] == 1) {
                int k = 0;
                double txk = xc1[ir][0];
                double tyk = yc1[ir][0];
                double tzk = zc1[ir][0];
                double bsqk = bsq1[ir][0];
                double bk = b1[ir][0];
                intag[0] = intag1[ir][0];
                double arcsum = pix2;
                ib = ib + 1;
                arclen += gr[ir][k].value * arcsum;
                if (!moved) {
                    int in = intag[k];
                    double t1 = arcsum * rrisq * (bsqk - rrisq + r[in] * r[in])
                            / (rri2 * bsqk * bk);
                    grad.sub(threadID, ir, txk * t1 * wght, tyk * t1 * wght, tzk * t1 * wght);
                    grad.add(threadID, in, txk * t1 * wght, tyk * t1 * wght, tzk * t1 * wght);
                }
                area[ir] = ib * pix2 + exang + arclen;
                area[ir] = area[ir] % pix4;
                return;
            }
            /*
              General case where more than one sphere intersects the
              current sphere; sort intersecting spheres by their degree of
              overlap with the current main sphere
             */
            sort(gr[ir], 0, count[ir]);
            for (int j = 0; j < count[ir]; j++) {
                int k = gr[ir][j].key;
                intag[j] = intag1[ir][k];
                xc[j] = xc1[ir][k];
                yc[j] = yc1[ir][k];
                zc[j] = zc1[ir][k];
                dsq[j] = dsq1[ir][k];
                b[j] = b1[ir][k];
                bsq[j] = bsq1[ir][k];
                omit[j] = false;
            }

            // Radius of the each circle on the surface of the "ir" sphere.
            for (int i = 0; i < count[ir]; i++) {
                double gi = gr[ir][i].value * rri;
                bg[i] = b[i] * gi;
                risq[i] = rrisq - gi * gi;
                ri[i] = sqrt(risq[i]);
                ther[i] = pid2 - asin(min(1.0, max(-1.0, gr[ir][i].value)));
            }

            // Find boundary of inaccessible area on "ir" sphere.
            for (int k = 0; k < count[ir] - 1; k++) {
                if (omit[k]) {
                    continue;
                }
                double txk = xc[k];
                double tyk = yc[k];
                double tzk = zc[k];
                double bk = b[k];
                double therk = ther[k];
                for (j = k + 1; j < count[ir]; j++) {
                    if (omit[j]) {
                        continue;
                    }
                    /*
                      Check to see if J circle is intersecting K circle;
                      get distance between circle centers and sum of radii.
                     */
                    double cc = (txk * xc[j] + tyk * yc[j] + tzk * zc[j]) / (bk * b[j]);

                    // Check acos FORTRAN vs. Java.
                    cc = acos(min(1.0, max(-1.0, cc)));
                    double td = therk + ther[j];
                    // Check to see if circles enclose separate regions
                    if (cc >= td) {
                        continue;
                    }
                    // Check for circle J completely inside circle K
                    if (cc + ther[j] < therk) {
                        omit[j] = true;
                        continue;
                    }
                    // Check for circles essentially parallel.
                    if (cc > delta) {
                        if (pix2 - cc <= td) {
                            area[ir] = 0.0;
                            return;
                        }
                    }
                }
            }

            // Find T value of circle intersections.
            for (int k = 0; k < count[ir]; k++) {
                if (omit[k]) {
                    continue; // goto 110
                }
                boolean komit = omit[k];
                omit[k] = true;
                int narc = 0;
                boolean top = false;
                double txk = xc[k];
                double tyk = yc[k];
                double tzk = zc[k];
                double dk = sqrt(dsq[k]);
                double bsqk = bsq[k];
                double bk = b[k];
                double gk = gr[ir][k].value * rri;
                double risqk = risq[k];
                double rik = ri[k];
                double therk = ther[k];

                // Rotation matrix elements.
                double t1 = tzk / (bk * dk);
                double axx = txk * t1;
                double axy = tyk * t1;
                double axz = dk / bk;
                double ayx = tyk / dk;
                double ayy = txk / dk;
                double azx = txk / bk;
                double azy = tyk / bk;
                double azz = tzk / bk;
                for (int l = 0; l < count[ir]; l++) {
                    if (omit[l]) {
                        continue;
                    }
                    double txl = xc[l];
                    double tyl = yc[l];
                    double tzl = zc[l];

                    // Rotate spheres so K vector collinear with z-axis.
                    double uxl = txl * axx + tyl * axy - tzl * axz;
                    double uyl = tyl * ayy - txl * ayx;
                    double uzl = txl * azx + tyl * azy + tzl * azz;
                    double cosine = min(1.0, max(-1.0, uzl / b[l]));
                    if (acos(cosine) < therk + ther[l]) {
                        double dsql = uxl * uxl + uyl * uyl;
                        double tb = uzl * gk - bg[l];
                        double txb = uxl * tb;
                        double tyb = uyl * tb;
                        double td = rik * dsql;
                        double tr2 = risqk * dsql - tb * tb;
                        tr2 = max(eps, tr2);
                        double tr = sqrt(tr2);
                        double txr = uxl * tr;
                        double tyr = uyl * tr;

                        // Get T values of intersection for K circle.
                        tb = (txb + tyr) / td;
                        tb = min(1.0, max(-1.0, tb));
                        double tk1 = acos(tb);
                        if (tyb - txr < 0.0) {
                            tk1 = pix2 - tk1;
                        }
                        tb = (txb - tyr) / td;
                        tb = min(1.0, max(-1.0, tb));
                        double tk2 = acos(tb);
                        if (tyb + txr < 0.0) {
                            tk2 = pix2 - tk2;
                        }
                        thec = (rrisq * uzl - gk * bg[l])
                                / (rik * ri[l] * b[l]);
                        double the = 0.0;
                        if (abs(thec) < 1.0) {
                            the = -acos(thec);
                        } else if (thec >= 1.0) {
                            the = 0.0;
                        } else if (thec <= -1.0) {
                            the = -PI;
                        }
                        /*
                          See if "tk1" is entry or exit point; check t=0
                          point; "ti" is exit point, "tf" is entry point.
                         */
                        cosine = min(1.0, max(-1.0, (uzl * gk - uxl * rik)
                                / (b[l] * rri)));
                        double ti, tf;
                        if ((acos(cosine) - ther[l]) * (tk2 - tk1) <= 0.0) {
                            ti = tk2;
                            tf = tk1;
                        } else {
                            ti = tk1;
                            tf = tk2;
                        }
                        narc += 1;
                        if (narc > maxarc) {
                            throw new EnergyException(format(" Increase value of MAXARC %d.", narc), false);
                        }
                        int narc1 = narc - 1;
                        if (tf <= ti) {
                            arcf[narc1] = tf;
                            arci[narc1] = new IndexedDouble(0.0, narc1);
                            tf = pix2;
                            lt[narc1] = l;
                            ex[narc1] = the;
                            top = true;
                            narc = narc + 1;
                            narc1 = narc - 1;
                        }
                        arcf[narc1] = tf;
                        arci[narc1] = new IndexedDouble(ti, narc1);
                        lt[narc1] = l;
                        ex[narc1] = the;
                        ux[l] = uxl;
                        uy[l] = uyl;
                        uz[l] = uzl;
                    }
                }
                omit[k] = komit;

                // Special case; K circle without intersections.
                if (narc <= 0) {
                    double arcsum = pix2;
                    ib += 1;
                    arclen += gr[ir][k].value * arcsum;
                    if (!moved) {
                        int in = intag[k];
                        t1 = arcsum * rrisq * (bsqk - rrisq + r[in] * r[in])
                                / (rri2 * bsqk * bk);
                        grad.sub(threadID, ir, txk * t1 * wght, tyk * t1 * wght, tzk * t1 * wght);
                        grad.add(threadID, in, txk * t1 * wght, tyk * t1 * wght, tzk * t1 * wght);
                    }
                    continue;
                }

                // General case; sum up arclength and set connectivity code.
                sort(arci, 0, narc);
                double arcsum = arci[0].value;
                int mi = arci[0].key;
                double t = arcf[mi];
                int ni = mi;
                for (j = 1; j < narc; j++) {
                    int m = arci[j].key;
                    if (t < arci[j].value) {
                        arcsum += (arci[j].value - t);
                        exang += ex[ni];
                        jb += 1;
                        if (jb >= maxarc) {
                            throw new EnergyException(format("Increase the value of MAXARC (%d).", jb), false);
                        }
                        int l = lt[ni];
                        ider[l] += 1;
                        sign_yder[l] += 1;
                        kent[jb] = maxarc * (l + 1) + (k + 1);
                        l = lt[m];
                        ider[l] += 1;
                        sign_yder[l] -= 1;
                        kout[jb] = maxarc * (k + 1) + (l + 1);
                    }
                    double tt = arcf[m];
                    if (tt >= t) {
                        t = tt;
                        ni = m;
                    }
                }
                arcsum += (pix2 - t);
                if (!top) {
                    exang += ex[ni];
                    jb = jb + 1;
                    int l = lt[ni];
                    ider[l] += 1;
                    sign_yder[l] += 1;
                    kent[jb] = maxarc * (l + 1) + (k + 1);
                    l = lt[mi];
                    ider[l] += 1;
                    sign_yder[l] -= 1;
                    kout[jb] = maxarc * (k + 1) + (l + 1);
                }

                // Calculate the surface area derivatives.
                for (int l = 0; l <= count[ir]; l++) {
                    if (ider[l] == 0) {
                        continue;
                    }
                    double rcn = ider[l] * rrisq;
                    ider[l] = 0;
                    double uzl = uz[l];
                    double gl = gr[ir][l].value * rri;
                    double bgl = bg[l];
                    double bsql = bsq[l];
                    double risql = risq[l];
                    double wxlsq = bsql - uzl * uzl;
                    double wxl = sqrt(wxlsq);
                    double p = bgl - gk * uzl;
                    double v = risqk * wxlsq - p * p;
                    v = max(eps, v);
                    v = sqrt(v);
                    t1 = rri * (gk * (bgl - bsql) + uzl * (bgl - rrisq))
                            / (v * risql * bsql);
                    double deal = -wxl * t1;
                    double decl = -uzl * t1 - rri / v;
                    double dtkal = (wxlsq - p) / (wxl * v);
                    double dtkcl = (uzl - gk) / v;
                    double s = gk * b[l] - gl * uzl;
                    t1 = 2.0 * gk - uzl;
                    double t2 = rrisq - bgl;
                    double dtlal = -(risql * wxlsq * b[l] * t1 - s * (wxlsq * t2 + risql * bsql))
                            / (risql * wxl * bsql * v);
                    double dtlcl = -(risql * b[l] * (uzl * t1 - bgl) - uzl * t2 * s)
                            / (risql * bsql * v);
                    double gaca = rcn * (deal - (gk * dtkal - gl * dtlal) / rri) / wxl;
                    double gacb = (gk - uzl * gl / b[l]) * sign_yder[l] * rri / wxlsq;
                    sign_yder[l] = 0;
                    if (!moved) {
                        double faca = ux[l] * gaca - uy[l] * gacb;
                        double facb = uy[l] * gaca + ux[l] * gacb;
                        double facc = rcn * (decl - (gk * dtkcl - gl * dtlcl) / rri);
                        double dax = axx * faca - ayx * facb + azx * facc;
                        double day = axy * faca + ayy * facb + azy * facc;
                        double daz = azz * facc - axz * faca;
                        int in = intag[l];
                        grad.add(threadID, ir, dax * wght, day * wght, daz * wght);
                        grad.sub(threadID, in, dax * wght, day * wght, daz * wght);
                    }

                }
                arclen += gr[ir][k].value * arcsum;
                if (!moved) {
                    int in = intag[k];
                    t1 = arcsum * rrisq * (bsqk - rrisq + r[in] * r[in]) / (rri2 * bsqk * bk);
                    grad.sub(threadID, ir, txk * t1 * wght, tyk * t1 * wght, tzk * t1 * wght);
                    grad.add(threadID, in, txk * t1 * wght, tyk * t1 * wght, tzk * t1 * wght);
                }
            }
            if (arclen == 0.0) {
                area[ir] = 0.0;
                return;
            }
            if (jb == 0) {
                area[ir] = ib * pix2 + exang + arclen;
                area[ir] = area[ir] % pix4;
                return;
            }

            // Find number of independent boundaries and check connectivity.
            j = 0;
            for (int k = 1; k <= jb; k++) {
                if (kout[k] == -1) {
                    continue;
                }
                i = k;
                boolean success = independentBoundaries(k, exang, jb, ir, arclen);
                if (success) {
                    return;
                }
            }
            ib = ib + 1;
            logger.log(Level.WARNING, format(" Connectivity error at atom %d", ir));
            area[ir] = 0.0;
        }

        /**
         * Find number of independent boundaries and check connectivity.
         *
         * @param k
         * @param exang
         * @param jb
         * @param ir
         * @param arclen
         */
        boolean independentBoundaries(int k, double exang, int jb, int ir, double arclen) {
            int m = kout[i];
            kout[i] = -1;
            j = j + 1;
            for (int ii = 1; ii <= jb; ii++) {
                if (m == kent[ii]) {
                    if (ii == k) {
                        ib++;
                        if (j == jb) {
                            area[ir] = ib * 2.0 * PI + exang + arclen;
                            area[ir] = area[ir] % (4.0 * PI);
                            return true;
                        }
                        return false;
                    }
                    i = ii;
                    return independentBoundaries(k, exang, jb, ir, arclen);
                }
            }
            return false;
        }
    }
}
