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
 * Refinement minimization class using {@link OptimizationListener} interface,
 * constructs a {@link RefinementEnergy} object for this purpose
 *
 * @author Tim Fenn
 */
public class RefinementMinimize implements OptimizationListener, Terminatable {

    /**
     * Different refinement mode selection types
     */
    public enum RefinementMode {

        /**
         * refine coordinates only
         */
        COORDINATES,
        /**
         * refine B factors only (if anisotropic, refined as such)
         */
        BFACTORS,
        /**
         * refine coordinates and B factors (if anisotropic, refined as such)
         */
        COORDINATES_AND_BFACTORS,
        /**
         * refine occupancies only (alternate conformers are constrained)
         */
        OCCUPANCIES,
        /**
         * refine B factors and occupancies
         */
        BFACTORS_AND_OCCUPANCIES,
        /**
         * refine coordinates and occupancies
         */
        COORDINATES_AND_OCCUPANCIES,
        /**
         * refine all
         */
        COORDINATES_AND_BFACTORS_AND_OCCUPANCIES
    }
    private static final Logger logger = Logger.getLogger(RefinementMinimize.class.getName());
    private static double toSeconds = 0.000000001;
    private final DataContainer data;
    private final RefinementModel refinementmodel;
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
    // recommended eps - accessible to groovy
    private double eps = 0.1;

    /**
     * constructor for refinement, assumes coordinates and B factor optimization
     *
     * @param data input {@link DataContainer} that will be used as the model,
     * must contain a {@link RefinementModel}
     */
    public RefinementMinimize(DataContainer data) {
        this(data, RefinementMode.COORDINATES_AND_BFACTORS);
    }

    /**
     * constructor for refinement
     *
     * @param data input {@link DataContainer} that will be used as the model,
     * must contain a {@link RefinementModel} and either {@link DiffractionData}
     * or {@link RealSpaceData}
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for refinement
     */
    public RefinementMinimize(DataContainer data, RefinementMode refinementmode) {
        this.data = data;
        this.refinementmodel = data.getRefinementModel();
        this.refinementMode = refinementmode;
        this.atomarray = data.getAtomArray();
        this.nAtoms = atomarray.length;

        refinementenergy = new RefinementEnergy(data, refinementmode, null);

        this.nxyz = refinementenergy.nxyz;
        this.nb = refinementenergy.nb;
        this.nocc = refinementenergy.nocc;
        this.n = refinementenergy.getNumberOfVariables();

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];

        refinementenergy.getCoordinates(x);

        double xyzscale = 1.0;
        double bisoscale = 1.0;
        double anisouscale = 50.0;
        double occscale = 15.0;

        if (refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            bisoscale = 0.4;
            anisouscale = 40.0;
        }

        if (refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            occscale = 10.0;
        }

        // set up scaling
        if (refinementmode == RefinementMode.COORDINATES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            for (int i = 0; i < nxyz; i++) {
                scaling[i] = xyzscale;
            }
        }

        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            int i = nxyz;
            int resnum = -1;
            if (data instanceof DiffractionData) {
                DiffractionData diffractiondata = (DiffractionData) data;
                int nres = diffractiondata.nresiduebfactor + 1;
                for (Atom a : atomarray) {
                    // ignore hydrogens!!!
                    if (a.getAtomicNumber() == 1) {
                        continue;
                    }
                    if (a.getAnisou() != null) {
                        for (int j = 0; j < 6; j++) {
                            scaling[i + j] = anisouscale;
                        }
                        i += 6;
                    } else if (diffractiondata.residuebfactor) {
                        if (resnum != a.getResidueNumber()) {
                            if (nres >= diffractiondata.nresiduebfactor) {
                                if (resnum > -1
                                        && i < nxyz + nb - 1) {
                                    i++;
                                }
                                if (i < nxyz + nb) {
                                    scaling[i] = bisoscale;
                                }
                                nres = 1;
                            } else {
                                nres++;
                            }
                            resnum = a.getResidueNumber();
                        }
                    } else {
                        scaling[i] = bisoscale;
                        i++;
                    }
                }
            } else {
                logger.severe("B refinement not supported for this data type!");
            }
        }

        if (refinementmode == RefinementMode.OCCUPANCIES
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            if (data instanceof DiffractionData) {
                int i = nxyz + nb;
                for (ArrayList<Residue> list : refinementmodel.altresidues) {
                    for (int j = 0; j < list.size(); j++) {
                        scaling[i] = occscale;
                        i++;
                    }
                }
                for (ArrayList<Molecule> list : refinementmodel.altmolecules) {
                    for (int j = 0; j < list.size(); j++) {
                        scaling[i] = occscale;
                        i++;
                    }
                }
            } else {
                logger.severe("occupancy refinement not supported for this data type!");
            }
        }

        refinementenergy.setScaling(scaling);
    }

    public Double getEps() {
        boolean hasaniso = false;

        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() != null) {
                hasaniso = true;
                break;
            }
        }

        switch (refinementMode) {
            case COORDINATES:
                eps = 0.1;
                break;
            case BFACTORS:
                if (hasaniso) {
                    eps = 20.0;
                } else {
                    eps = 0.01;
                }
                break;
            case COORDINATES_AND_BFACTORS:
                if (hasaniso) {
                    eps = 20.0;
                } else {
                    eps = 0.1;
                }
                break;
            case OCCUPANCIES:
                eps = 0.1;
                break;
            case COORDINATES_AND_OCCUPANCIES:
                eps = 0.1;
                break;
            case BFACTORS_AND_OCCUPANCIES:
                if (hasaniso) {
                    eps = 20.0;
                } else {
                    eps = 0.01;
                }
                break;
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                if (hasaniso) {
                    eps = 20.0;
                } else {
                    eps = 0.1;
                }
                break;
        }

        return eps;
    }

    /**
     * get the number of xyz parameters being fit
     * @return the number of xyz parameters
     */
    public int getNXYZ() {
        return nxyz;
    }

    /**
     * get the number of B factor parameters being fit
     * @return the number of B factor parameters
     */
    public int getNB() {
        return nb;
    }

    /**
     * get the number of occupancy parameters being fit
     * @return the number of occupancy parameters
     */
    public int getNOcc() {
        return nocc;
    }

    /**
     * minimize assuming an eps of 1.0 and Integer.MAX_VALUE cycles
     *
     * @return {@link RefinementEnergy} result
     */
    public RefinementEnergy minimize() {
        return minimize(1.0);
    }

    /**
     * minimize assuming Integer.MAX_VALUE cycles
     *
     * @param eps input gradient rms desired
     * @return {@link RefinementEnergy} result
     */
    public RefinementEnergy minimize(double eps) {
        return minimize(7, eps, Integer.MAX_VALUE - 2);
    }

    /**
     * minimize assuming an eps of 1.0 and limited cycles
     *
     * @param maxiter maximum iterations allowed
     * @return {@link RefinementEnergy} result
     */
    public RefinementEnergy minimize(int maxiter) {
        return minimize(7, 1.0, maxiter);
    }

    /**
     * minimize with input eps and cycles
     *
     * @param eps input gradient rms desired
     * @param maxiter maximum iterations allowed
     * @return {@link RefinementEnergy} result
     */
    public RefinementEnergy minimize(double eps, int maxiter) {
        return minimize(7, eps, maxiter);
    }

    /**
     * minimize with input cycles for matrix conditioning, eps and cycles
     *
     * @param m number of cycles of matrix updates
     * @param eps input gradient rms desired
     * @param maxiter maximum iterations allowed
     * @return {@link RefinementEnergy} result
     */
    public RefinementEnergy minimize(int m, double eps, int maxiter) {
        String typestring;
        if (data instanceof DiffractionData) {
            typestring = "X-ray";
        } else if (data instanceof RealSpaceData) {
            typestring = "Real Space";
        } else {
            typestring = "null";
        }

        switch (refinementMode) {
            case COORDINATES:
                logger.info("Beginning " + typestring + " refinement - mode: coordinates nparams: " + n);
                break;
            case BFACTORS:
                logger.info("Beginning " + typestring + " refinement - mode: bfactors nparams: " + n);
                break;
            case COORDINATES_AND_BFACTORS:
                logger.info("Beginning " + typestring + " refinement - mode: coordinates and bfactors nparams: " + n);
                break;
            case OCCUPANCIES:
                logger.info("Beginning " + typestring + " refinement - mode: occupancies nparams: " + n);
                break;
            case COORDINATES_AND_OCCUPANCIES:
                logger.info("Beginning " + typestring + " refinement - mode: coordinates and occupancies nparams: " + n);
                break;
            case BFACTORS_AND_OCCUPANCIES:
                logger.info("Beginning " + typestring + " refinement - mode: bfactors and occupancies nparams: " + n);
                break;
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                logger.info("Beginning " + typestring + " refinement - mode: coordinates and bfactors and occupancies nparams: " + n);
                break;
        }

        refinementenergy.getCoordinates(x);
        /**
         * Scale coordinates.
         */
        for (int i = 0; i < n; i++) {
            x[i] *= scaling[i];
        }

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = 0;
        double e = refinementenergy.energyAndGradient(x, grad);
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
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time      "
                    + data.printOptimizationHeader() + "\n");
        }
        if (info == null) {
            logger.info(String.format("%6d %13.4g %11.4g\n",
                    iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8.3g ",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
                sb.append(data.printOptimizationUpdate());
                sb.append(String.format("\n"));
                logger.info(sb.toString());
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
