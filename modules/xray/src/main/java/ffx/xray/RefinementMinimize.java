/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
package ffx.xray;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.Terminatable;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;

/**
 *
 * @author fennt
 */
public class RefinementMinimize implements OptimizationListener, Terminatable {

    public enum RefinementMode {

        COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS
    }
    private static final Logger logger = Logger.getLogger(RefinementMinimize.class.getName());
    private static double toSeconds = 0.000000001;
    private final MolecularAssembly molecularAssembly[];
    private final XRayStructure xraystructure;
    private final Atom atomarray[];
    private final int nAtoms;
    private final RefinementEnergy refinementenergy;
    private RefinementMode refinementMode;
    private int n;
    private int nxyz;
    private int nb;
    private int nocc;
    private final double x[];
    private final double grad[];
    private final double scaling[];
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    public RefinementMinimize(MolecularAssembly molecularAssembly,
            XRayStructure xraystructure) {
        this(new MolecularAssembly[]{molecularAssembly}, xraystructure,
                RefinementMode.COORDINATES_AND_BFACTORS);
    }

    public RefinementMinimize(MolecularAssembly molecularAssembly,
            XRayStructure xraystructure, RefinementMode refinementmode) {
        this(new MolecularAssembly[]{molecularAssembly}, xraystructure,
                refinementmode);
    }

    public RefinementMinimize(MolecularAssembly molecularAssembly[],
            XRayStructure xraystructure, RefinementMode refinementmode) {
        this.molecularAssembly = molecularAssembly;
        this.xraystructure = xraystructure;
        this.refinementMode = refinementmode;
        this.atomarray = xraystructure.atomarray;
        this.nAtoms = atomarray.length;
        if (!xraystructure.scaled) {
            xraystructure.scalebulkfit();
            xraystructure.printstats();
        }

        // determine size of fit
        n = nxyz = nb = nocc = 0;
        switch (refinementmode) {
            case COORDINATES:
                nxyz = nAtoms * 3;
                break;
            case BFACTORS:
                for (Atom a : atomarray) {
                    // ignore hydrogens!!!
                    if (a.getAtomicNumber() == 1) {
                        continue;
                    }
                    if (a.getAnisou() == null) {
                        nb++;
                    } else {
                        nb += 6;
                    }
                }
                break;
            case COORDINATES_AND_BFACTORS:
                nxyz = nAtoms * 3;
                for (Atom a : atomarray) {
                    // ignore hydrogens!!!
                    if (a.getAtomicNumber() == 1) {
                        continue;
                    }
                    if (a.getAnisou() == null) {
                        nb++;
                    } else {
                        nb += 6;
                    }
                }
                break;
        }
        n = nxyz + nb + nocc;
        refinementenergy = new RefinementEnergy(molecularAssembly,
                xraystructure, nxyz, nb, nocc, refinementmode);

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];

        if (refinementmode == RefinementMode.COORDINATES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS) {
            refinementenergy.xrayEnergy.getCoordinates(x);
        }

        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS) {
            refinementenergy.xrayEnergy.getBFactors(x, nxyz);
        }

        for (int i = 0; i < n; i++) {
            scaling[i] = 1.0;
        }

        refinementenergy.setScaling(scaling);
    }

    public RefinementEnergy minimize() {
        return minimize(1.0);
    }

    public RefinementEnergy minimize(double eps) {
        return minimize(7, eps);
    }

    public RefinementEnergy minimize(int m, double eps) {

        switch (refinementMode) {
            case COORDINATES:
                logger.info("Beginning X-ray refinement - mode: coordinates nparams: " + n);
                break;
            case BFACTORS:
                logger.info("Beginning X-ray refinement - mode: bfactors nparams: " + n);
                break;
            case COORDINATES_AND_BFACTORS:
                logger.info("Beginning X-ray refinement - mode: coordinates and bfactors nparams: " + n);
                break;
        }

        double e = refinementenergy.energyAndGradient(x, grad);

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, refinementenergy, this);
        done = true;
        switch (status) {
            case 0:
                logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
                break;
            case 1:
                logger.info(String.format("\n Optimization terminated at step %d.\n", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.\n");
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            mtime += System.nanoTime();
            sb.append(String.format("minimizer time: %g\n", mtime * toSeconds));
            logger.info(sb.toString());
        }

        return refinementenergy;
    }

    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df, double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time      R  Rfree\n");
        }
        if (info == null) {
            logger.info(String.format("%6d %13.4g %11.4g\n",
                    iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8.3g %6.2f %6.2f\n",
                        iter, f, grms, df, xrms, angle, nfun, seconds,
                        xraystructure.crystalstats.get_r(),
                        xraystructure.crystalstats.get_rfree()));
            } else {
                logger.info(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8s\n",
                        iter, f, grms, df, xrms, angle, nfun, info.toString()));
            }
        }
        if (terminate) {
            logger.info(" The optimization recieved a termination request.");
            // Tell the L-BFGS optimizer to terminate.
            return false;
        }
        return true;
    }

    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Exception terminating minimization.\n", e);
                }
            }
        }
    }
}
