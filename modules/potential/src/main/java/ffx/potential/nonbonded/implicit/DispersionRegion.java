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
import ffx.potential.parameters.VDWType;

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
    private static final double DISP_OVERLAP_SCALE_FACTOR = 0.81;
    /**
     * This value was described as 0.36 in the original 2007 model (see Schnieders thesis)
     * and more recently the value was reduced to 0.26.
     */
    public static final double DEFAULT_DISP_OFFSET = 1.4;
    private static final double SLEVY = 1.0;
    private static final double AWATER = 0.033428;
    private static final double EPSO = 0.1100;
    private static final double EPSH = 0.0135;
    private static final double RMINO = 1.7025;
    private static final double RMINH = 1.3275;
    private double dispersionOffest = DEFAULT_DISP_OFFSET;

    public DispersionRegion(int nt, Atom[] atoms) {
        dispersionLoop = new DispersionLoop[nt];
        for (int i = 0; i < nt; i++) {
            dispersionLoop[i] = new DispersionLoop();
        }
        sharedDispersion = new SharedDouble();
        allocate(atoms);
    }

    /**
     * The dispersion integral begins offset from the vdW radius.
     *
     * @param dispersionOffest the dispersion integral offset.
     */
    public void setDispersionOffest(double dispersionOffest) {
        this.dispersionOffest = dispersionOffest;
        // Update the maximum dispersion energy.
        maxDispersionEnergy();
    }

    /**
     * The dispersion integral begins offset from the vdW radius.
     *
     * @return the dispersion integral offset.
     */
    public double getDispersionOffset() {
        return dispersionOffest;
    }

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
            // Do not apply DISP_OFFSET here -- we don't want to change the shape of the Buffered-14-7 curve.
            double rmini = type.radius / 2.0;
            if (rDisp[i] > 0.0 && epsi > 0.0) {
                double sqEpsoEpsi = sqrt(EPSO) + sqrt(epsi);
                double sqEpshEpsi = sqrt(EPSH) + sqrt(epsi);
                double emixo = 4.0 * EPSO * epsi / (pow(sqEpsoEpsi, 2));
                double rmixo = 2.0 * (pow(RMINO, 3) + pow(rmini, 3)) / (pow(RMINO, 2) + pow(rmini, 2));
                double rmixo3 = pow(rmixo, 3);
                double rmixo7 = pow(rmixo, 7);
                double ao = emixo * rmixo7;
                double emixh = 4.0 * EPSH * epsi / (pow(sqEpshEpsi, 2));
                double rmixh = 2.0 * (pow(RMINH, 3) + pow(rmini, 3)) / (pow(RMINH, 2) + pow(rmini, 2));
                double rmixh3 = pow(rmixh, 3);
                double rmixh7 = pow(rmixh, 7);
                double ah = emixh * rmixh7;
                // Apply the DISP_OFFSET here to start the integral beyond the atomic radius of atom i.
                double ri = rDisp[i] + dispersionOffest;
                double ri3 = pow(ri, 3);
                double ri7 = pow(ri, 7);
                double ri11 = pow(ri, 11);
                if (ri < rmixh) {
                    cDisp[i] = -4.0 * PI * emixh * (rmixh3 - ri3) / 3.0;
                    cDisp[i] = cDisp[i] - emixh * 18.0 / 11.0 * rmixh3 * PI;
                } else {
                    cDisp[i] = 2.0 * PI * (2.0 * rmixh7 - 11.0 * ri7) * ah;
                    cDisp[i] = cDisp[i] / (11.0 * ri11);
                }
                cDisp[i] = 2.0 * cDisp[i];
                if (ri < rmixo) {
                    cDisp[i] = cDisp[i] - 4.0 * PI * emixo * (rmixo3 - ri3) / 3.0;
                    cDisp[i] = cDisp[i] - emixo * 18.0 / 11.0 * rmixo3 * PI;
                } else {
                    cDisp[i] = cDisp[i] + 2.0 * PI * (2.0 * rmixo7 - 11.0 * ri7) * ao / (11.0 * ri11);
                }
            }
            cDisp[i] = SLEVY * AWATER * cDisp[i];
            // logger.info(format(" Max dispersion: %d %8.6f", i, cDisp[i]));
        }
    }

    @Override
    public void start() {
        sharedDispersion.set(0.0);
    }

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
                        sum += descreen(i, k);

                        // Flip the sign on {xr, yr, zr};
                        xr = -xr;
                        yr = -yr;
                        zr = -zr;

                        // Atom k descreened by atom i.
                        sum += descreen(k, i);
                    }
                }

                // Subtract descreening.
                edisp -= SLEVY * AWATER * sum;
            }
        }

        private double descreen(int i, int k) {
            double sum = 0.0;
            VDWType type = atoms[i].getVDWType();
            double epsi = type.wellDepth;
            // Do not apply DISP_OFFSET to rmini -- we don't want to change the shape of the Buffered-14-7 curve.
            double rmini = type.radius / 2.0;
            double emixo = (4.0 * EPSO * epsi) / (pow(sqrt(EPSO) + sqrt(epsi), 2));
            double rmixo = 2.0 * (pow(RMINO, 3) + pow(rmini, 3)) / (pow(RMINO, 2) + pow(rmini, 2));
            double rmixo7 = pow(rmixo, 7);
            double ao = emixo * rmixo7;
            double emixh = 4.0 * EPSH * epsi / (pow(sqrt(EPSH) + sqrt(epsi), 2));
            double rmixh = 2.0 * (pow(RMINH, 3) + pow(rmini, 3)) / (pow(RMINH, 2) + pow(rmini, 2));
            double rmixh7 = pow(rmixh, 7);
            double ah = emixh * rmixh7;
            // Apply the DISP_OFFSET here to start the integral beyond the atomic radius of atom i.
            double ri = rDisp[i] + dispersionOffest;
            // Atom k descreens with no offset applied.
            double rk = rDisp[k];
            double sk = rk * DISP_OVERLAP_SCALE_FACTOR;
            double sk2 = sk * sk;
            if (ri < r + sk) {
                double de = 0.0;
                double rmax = max(ri, r - sk);
                double lik = rmax;
                double lik2 = lik * lik;
                double lik3 = lik2 * lik;
                double lik4 = lik3 * lik;
                if (lik < rmixo) {
                    double uik = min(r + sk, rmixo);
                    double uik2 = uik * uik;
                    double uik3 = uik2 * uik;
                    double uik4 = uik3 * uik;
                    double term = 4.0 * PI / (48.0 * r) * (3.0 * (lik4 - uik4)
                            - 8.0 * r * (lik3 - uik3) + 6.0 * (r2 - sk2) * (lik2 - uik2));
                    double iwca = -emixo * term;
                    sum = sum + iwca;
                    if (gradient) {
                        double dl;
                        if (ri > r - sk) {
                            dl = -lik2 + 2.0 * r2 + 2.0 * sk2;
                            dl = dl * lik2;
                        } else {
                            dl = -lik3 + 4.0 * lik2 * r
                                    - 6.0 * lik * r2
                                    + 2.0 * lik * sk2 + 4.0 * r3
                                    - 4.0 * r * sk2;
                            dl = dl * lik;
                        }
                        double du;
                        if (r + sk > rmixo) {
                            du = -uik2 + 2.0 * r2 + 2.0 * sk2;
                            du = -du * uik2;
                        } else {
                            du = -uik3 + 4.0 * uik2 * r
                                    - 6.0 * uik * r2
                                    + 2.0 * uik * sk2 + 4.0 * r3
                                    - 4.0 * r * sk2;
                            du = -du * uik;
                        }
                        de = de - emixo * PI * (dl + du) / (4.0 * r2);
                    }
                }
                if (lik < rmixh) {
                    double uik = min(r + sk, rmixh);
                    double uik2 = uik * uik;
                    double uik3 = uik2 * uik;
                    double uik4 = uik3 * uik;
                    double term = 4.0 * PI / (48.0 * r) * (3.0 * (lik4 - uik4)
                            - 8.0 * r * (lik3 - uik3) + 6.0 * (r2 - sk2) * (lik2 - uik2));
                    double iwca = -2.0 * emixh * term;
                    sum = sum + iwca;
                    if (gradient) {
                        double dl;
                        if (ri > r - sk) {
                            dl = -lik2 + 2.0 * r2 + 2.0 * sk2;
                            dl = dl * lik2;
                        } else {
                            dl = -lik3 + 4.0 * lik2 * r - 6.0 * lik * r2
                                    + 2.0 * lik * sk2 + 4.0 * r3 - 4.0 * r * sk2;
                            dl = dl * lik;
                        }
                        double du;
                        if (r + sk > rmixh) {
                            du = -uik2 + 2.0 * r2 + 2.0 * sk2;
                            du = -du * uik2;
                        } else {
                            du = -uik3 + 4.0 * uik2 * r - 6.0 * uik * r2
                                    + 2.0 * uik * sk2 + 4.0 * r3 - 4.0 * r * sk2;
                            du = -du * uik;
                        }
                        de = de - 2.0 * emixh * PI * (dl + du) / (4.0 * r2);
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
                    double lik13 = lik12 * lik;
                    double term = 4.0 * PI / (120.0 * r * lik5 * uik5) * (15.0 * uik * lik * r * (uik4 - lik4)
                            - 10.0 * uik2 * lik2 * (uik3 - lik3) + 6.0 * (sk2 - r2) * (uik5 - lik5));
                    double term2 = 4.0 * PI / (2640.0 * r * lik12 * uik12) * (120.0 * uik * lik * r * (uik11 - lik11)
                            - 66.0 * uik2 * lik2 * (uik10 - lik10) + 55.0 * (sk2 - r2) * (uik12 - lik12));
                    double idisp = -2.0 * ao * term;
                    double irep = ao * rmixo7 * term2;
                    sum = sum + irep + idisp;
                    if (gradient) {
                        double dl;
                        if (ri > r - sk || rmax < rmixo) {
                            dl = -5.0 * lik2 + 3.0 * r2 + 3.0 * sk2;
                            dl = -dl / lik5;
                        } else {
                            dl = 5.0 * lik3 - 33.0 * lik * r2 - 3.0 * lik * sk2
                                    + 15.0 * (lik2 * r + r3 - r * sk2);
                            dl = dl / lik6;
                        }
                        double du;
                        du = 5.0 * uik3 - 33.0 * uik * r2 - 3.0 * uik * sk2
                                + 15.0 * (uik2 * r + r3 - r * sk2);
                        du = -du / uik6;
                        de = de - 2.0 * ao * PI * (dl + du) / (15.0 * r2);

                        if (ri > r - sk || rmax < rmixo) {
                            dl = -6.0 * lik2 + 5.0 * r2 + 5.0 * sk2;
                            dl = -dl / lik12;
                        } else {
                            dl = 6.0 * lik3 - 125.0 * lik * r2 - 5.0 * lik * sk2 + 60.0 * (lik2 * r + r3 - r * sk2);
                            dl = dl / lik13;
                        }
                        du = 6.0 * uik3 - 125.0 * uik * r2 - 5.0 * uik * sk2 + 60.0 * (uik2 * r + r3 - r * sk2);
                        du = -du / uik13;
                        de = de + ao * rmixo7 * PI * (dl + du) / (60.0 * r2);
                    }

                }
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
                    double lik13 = lik12 * lik;
                    double term = 4.0 * PI / (120.0 * r * lik5 * uik5) * (15.0 * uik * lik * r * (uik4 - lik4)
                            - 10.0 * uik2 * lik2 * (uik3 - lik3) + 6.0 * (sk2 - r2) * (uik5 - lik5));
                    double term2 = 4.0 * PI / (2640.0 * r * lik12 * uik12) * (120.0 * uik * lik * r * (uik11 - lik11)
                            - 66.0 * uik2 * lik2 * (uik10 - lik10) + 55.0 * (sk2 - r2) * (uik12 - lik12));
                    double idisp = -4.0 * ah * term;
                    double irep = 2.0 * ah * rmixh7 * term2;
                    sum = sum + irep + idisp;
                    if (gradient) {
                        double dl;
                        if (ri > r - sk || rmax < rmixh) {
                            dl = -5.0 * lik2 + 3.0 * r2 + 3.0 * sk2;
                            dl = -dl / lik5;
                        } else {
                            dl = 5.0 * lik3 - 33.0 * lik * r2
                                    - 3.0 * lik * sk2 + 15.0 * (lik2 * r + r3 - r * sk2);
                            dl = dl / lik6;
                        }
                        double du;
                        du = 5.0 * uik3 - 33.0 * uik * r2
                                - 3.0 * uik * sk2 + 15.0 * (uik2 * r + r3 - r * sk2);
                        du = -du / uik6;
                        de = de - 4.0 * ah * PI * (dl + du) / (15.0 * r2);
                        if (ri > r - sk || rmax < rmixh) {
                            dl = -6.0 * lik2 + 5.0 * r2 + 5.0 * sk2;
                            dl = -dl / lik12;
                        } else {
                            dl = 6.0 * lik3 - 125.0 * lik * r2 - 5.0 * lik * sk2 + 60.0 * (lik2 * r + r3 - r * sk2);
                            dl = dl / lik13;
                        }
                        du = 6.0 * uik3 - 125.0 * uik * r2 - 5.0 * uik * sk2 + 60.0 * (uik2 * r + r3 - r * sk2);
                        du = -du / uik13;
                        de = de + ah * rmixh7 * PI * (dl + du) / (30.0 * r2);
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
    }
}