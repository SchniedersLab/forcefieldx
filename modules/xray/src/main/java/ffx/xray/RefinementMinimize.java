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
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;

/**
 *
 * @author fennt
 */
public class RefinementMinimize implements OptimizationListener, Terminatable {

    public enum RefinementMode {

        COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES,
        BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES,
        COORDINATES_AND_BFACTORS_AND_OCCUPANCIES
    }
    private static final Logger logger = Logger.getLogger(RefinementMinimize.class.getName());
    private static double toSeconds = 0.000000001;
    private final MolecularAssembly molecularAssembly[];
    private final XRayStructure xraystructure;
    private final RefinementData refinementdata;
    private final Atom atomarray[];
    private final int nAtoms;
    private final RefinementEnergy refinementenergy;
    private RefinementMode refinementMode;
    private final int nxyz;
    private final int nb;
    private final int nocc;
    private final int n;
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
            XRayStructure xraystructure) {
        this(molecularAssembly, xraystructure,
                RefinementMode.COORDINATES_AND_BFACTORS);
    }

    public RefinementMinimize(MolecularAssembly molecularAssembly[],
            XRayStructure xraystructure, RefinementMode refinementmode) {
        this.molecularAssembly = molecularAssembly;
        this.xraystructure = xraystructure;
        this.refinementdata = xraystructure.refinementdata;
        this.refinementMode = refinementmode;
        this.atomarray = xraystructure.atomarray;
        this.nAtoms = atomarray.length;
        if (!xraystructure.scaled) {
            xraystructure.scalebulkfit();
            xraystructure.printstats();
        }
        refinementenergy = new RefinementEnergy(molecularAssembly,
                xraystructure, refinementmode, null);

        this.nxyz = refinementenergy.nxyz;
        this.nb = refinementenergy.nb;
        this.nocc = refinementenergy.nocc;
        this.n = refinementenergy.getNumberOfVariables();

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];

        refinementenergy.xrayEnergy.getCoordinates(x);

        double xyzscale = 1.0;
        double anisouscale = 80.0;
        double bisoscale = 1.0;
        double occscale = 15.0;
        if (refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            bisoscale = 0.2;
        }

        if (refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            occscale = 15.0;
        }

        // set up scaling
        if (refinementmode == RefinementMode.COORDINATES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            for (int i = 0; i < nxyz; i++) {
                scaling[i] = xyzscale;
                x[i] *= xyzscale;
            }
        }

        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            int i = nxyz;
            int resnum = -1;
            int nres = refinementdata.nresiduebfactor + 1;
            for (Atom a : atomarray) {
                // ignore hydrogens!!!
                if (a.getAtomicNumber() == 1) {
                    continue;
                }
                if (a.getAnisou() != null) {
                    for (int j = 0; j < 6; j++) {
                        scaling[i + j] = anisouscale;
                        x[i + j] *= anisouscale;
                    }
                    i += 6;
                } else if (refinementdata.residuebfactor) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= refinementdata.nresiduebfactor) {
                            if (resnum > -1
                                    && i < nxyz + nb - 1) {
                                i++;
                            }
                            if (i < nxyz + nb) {
                                scaling[i] = bisoscale;
                                x[i] *= bisoscale;
                            }
                            nres = 1;
                        } else {
                            nres++;
                        }
                        resnum = a.getResidueNumber();
                    }
                } else {
                    scaling[i] = bisoscale;
                    x[i] *= bisoscale;
                    i++;
                }
            }
        }

        if (refinementmode == RefinementMode.OCCUPANCIES
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            int i = nxyz + nb;
            for (ArrayList<Residue> list : xraystructure.altresidues) {
                for (int j = 0; j < list.size(); j++) {
                    scaling[i] = occscale;
                    x[i] *= occscale;
                    i++;
                }
            }
            for (ArrayList<Molecule> list : xraystructure.altmolecules) {
                for (int j = 0; j < list.size(); j++) {
                    scaling[i] = occscale;
                    x[i] *= occscale;
                    i++;
                }
            }
        }

        refinementenergy.setScaling(scaling);
    }

    public RefinementEnergy minimize() {
        return minimize(1.0);
    }

    public RefinementEnergy minimize(double eps) {
        return minimize(7, eps, Integer.MAX_VALUE - 2);
    }

    public RefinementEnergy minimize(int maxiter) {
        return minimize(7, 1.0, maxiter);
    }

    public RefinementEnergy minimize(double eps, int maxiter) {
        return minimize(7, eps, maxiter);
    }

    public RefinementEnergy minimize(int m, double eps, int maxiter) {
        refinementenergy.xrayEnergy.setRefinementMode(refinementMode);
        xraystructure.setXRayEnergy(refinementenergy.xrayEnergy);

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
            case OCCUPANCIES:
                logger.info("Beginning X-ray refinement - mode: occupancies nparams: " + n);
                break;
            case COORDINATES_AND_OCCUPANCIES:
                logger.info("Beginning X-ray refinement - mode: coordinates and occupancies nparams: " + n);
                break;
            case BFACTORS_AND_OCCUPANCIES:
                logger.info("Beginning X-ray refinement - mode: bfactors and occupancies nparams: " + n);
                break;
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                logger.info("Beginning X-ray refinement - mode: coordinates and bfactors and occupancies nparams: " + n);
                break;
        }

        double e = refinementenergy.energyAndGradient(x, grad);

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = 0;
        status = LBFGS.minimize(n, m, x, e, grad, eps, maxiter, refinementenergy, this);
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
