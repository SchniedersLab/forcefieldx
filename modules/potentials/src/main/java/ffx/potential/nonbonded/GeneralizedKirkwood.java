/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.nonbonded;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;

import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import java.util.logging.Level;

/**
 * This Generalized Kirkwood class implements GK for the AMOEBA polarizable
 * mutlipole force field in parallel using a {@link NeighborList}.
 *
 * @author Michael J. Schnieders<br>
 *         derived from:<br>
 *         TINKER code by Michael J. Schnieders and Jay W. Ponder<br>
 *
 * @see <a href="http://dx.doi.org/10.1021/ct7001336" target="_blank">M. J.
 *      Schnieders and J. W. Ponder, Polarizable atomic multipole solutes in a
 *      generalized Kirkwood continuum, Journal of Chemical Theory and
 *      Computation 2007, 3, (6), 2083-2097.</a><br>
 */
public class GeneralizedKirkwood {

    private static final Logger logger = Logger.getLogger(GeneralizedKirkwood.class.getName());
    private final Atom atoms[];
    private final double x[], y[], z[];
    private final double baseRadius[];
    private final double overlapScale[];
    private final double born[];
    private final int nAtoms;
    private final ForceField forceField;
    private final ParticleMeshEwald particleMeshEwald;
    private final ParallelTeam parallelTeam;
    private final BornRadiiRegion bornRadiiRegion;

    public GeneralizedKirkwood(ForceField forceField, Atom[] atoms,
                               ParticleMeshEwald particleMeshEwald,
                               ParallelTeam parallelTeam) {

        this.forceField = forceField;
        this.atoms = atoms;
        this.particleMeshEwald = particleMeshEwald;
        this.parallelTeam = parallelTeam;
        nAtoms = atoms.length;

        x = new double[nAtoms];
        y = new double[nAtoms];
        z = new double[nAtoms];
        baseRadius = new double[nAtoms];
        overlapScale = new double[nAtoms];
        born = new double[nAtoms];

        bornRadiiRegion = new BornRadiiRegion(parallelTeam.getThreadCount());
    }

    public void computeBornRadii() {
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            double xyz[] = ai.getXYZ();
            x[i] = xyz[0];
            y[i] = xyz[1];
            z[i] = xyz[2];
        }

        try {
            parallelTeam.execute(bornRadiiRegion);
        } catch (Exception e) {
            String message = "Fatal exception computing Born radii.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private void bornRadiiGradients() {
    }
    private static final double third = 1.0 / 3.0;
    private static final double pi43 = 4.0 * third * PI;
    private static final double pi12 = PI / 12.0;

    /**
     * Compute Born radii in parallel via the Grycuk method.
     *
     * @since 1.0
     */
    private class BornRadiiRegion extends ParallelRegion {

        private final BornRadiiLoop bornRadiiLoop[];

        public BornRadiiRegion(int nt) {
            bornRadiiLoop = new BornRadiiLoop[nt];
            for (int i = 0; i < nt; i++) {
                bornRadiiLoop[i] = new BornRadiiLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, bornRadiiLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class BornRadiiLoop extends IntegerForLoop {

            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    born[i] = 0.0;
                    final double ri = baseRadius[i];
                    if (ri > 0.0) {
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        double sum = pi43 / (ri * ri * ri);
                        for (int k = 0; k < nAtoms; k++) {
                            final double rk = baseRadius[k];
                            if (i != k && rk > 0.0) {
                                final double xr = x[k] - xi;
                                final double yr = y[k] - yi;
                                final double zr = z[k] - zi;
                                final double r2 = xr * xr + yr * yr + zr * zr;
                                final double r = sqrt(r2);
                                final double sk = rk * overlapScale[k];
                                final double sk2 = sk * sk;
                                if (ri + r < sk) {
                                    final double lik = ri;
                                    final double uik = sk - r;
                                    sum = sum + pi43 * (1.0 / (uik * uik * uik) - 1.0 / (lik * lik * lik));
                                }
                                final double uik = r + sk;
                                double lik;
                                if (ri + r < sk) {
                                    lik = sk - r;
                                } else if (r < ri + sk) {
                                    lik = ri;
                                } else {
                                    lik = r - sk;
                                }
                                final double l2 = lik * lik;
                                final double l4 = l2 * l2;
                                final double lr = lik * r;
                                final double l4r = l4 * r;
                                final double u2 = uik * uik;
                                final double u4 = u2 * u2;
                                final double ur = uik * r;
                                final double u4r = u4 * r;
                                final double term = (3.0 * (r2 - sk2) + 6.0 * u2 - 8.0 * ur) / u4r
                                                    - (3.0 * (r2 - sk2) + 6.0 * l2 - 8.0 * lr) / l4r;
                                sum = sum - pi12 * term;
                            }
                            born[i] = pow(sum / pi43, third);
                            if (born[i] <= 0.0) {
                                born[i] = 0.0001;
                            }
                            born[i] = 1.0 / born[i];
                        }
                    }
                }
            }
        }
    }
}
