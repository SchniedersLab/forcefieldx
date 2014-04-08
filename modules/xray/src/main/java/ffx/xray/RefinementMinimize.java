/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.xray;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;

/**
 * Refinement minimization class using {@link OptimizationListener} interface,
 * constructs a {@link RefinementEnergy} object for this purpose
 *
 * @author Timothy D. Fenn
 *
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
    public final RefinementEnergy refinementEnergy;
    private final AlgorithmListener listener;
    private static final Logger logger = Logger.getLogger(RefinementMinimize.class.getName());
    private static double toSeconds = 1.0e-9;
    private final DataContainer dataContainer;
    private final RefinementModel refinementModel;
    private final Atom atomArray[];
    private final int nAtoms;
    private RefinementMode refinementMode;
    private final int nXYZ;
    private final int nB;
    private final int nOcc;
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
        this(data, RefinementMode.COORDINATES_AND_BFACTORS, null);
    }

    /**
     * constructor for refinement
     *
     * @param data input {@link DataContainer} that will be used as the model,
     * must contain a {@link RefinementModel} and either {@link DiffractionData}
     * or {@link RealSpaceData}
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for
     * refinement
     */
    public RefinementMinimize(DataContainer data, RefinementMode refinementmode) {
        this(data, refinementmode, null);
    }

    /**
     * constructor for refinement
     *
     * @param data input {@link DataContainer} that will be used as the model,
     * must contain a {@link RefinementModel} and either {@link DiffractionData}
     * or {@link RealSpaceData}
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for
     * refinement
     * @param listener {@link AlgorithmListener} a listener for updates
     */
    public RefinementMinimize(DataContainer data, RefinementMode refinementmode,
            AlgorithmListener listener) {
        this.listener = listener;
        this.dataContainer = data;
        this.refinementModel = data.getRefinementModel();
        this.refinementMode = refinementmode;
        this.atomArray = data.getAtomArray();
        this.nAtoms = atomArray.length;

        refinementEnergy = new RefinementEnergy(data, refinementmode, null);

        this.nXYZ = refinementEnergy.nxyz;
        this.nB = refinementEnergy.nb;
        this.nOcc = refinementEnergy.nocc;
        this.n = refinementEnergy.getNumberOfVariables();

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];

        refinementEnergy.getCoordinates(x);

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
            for (int i = 0; i < nXYZ; i++) {
                scaling[i] = xyzscale;
            }
        }

        if (refinementmode == RefinementMode.BFACTORS
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            int i = nXYZ;
            int resnum = -1;
            if (data instanceof DiffractionData) {
                DiffractionData diffractiondata = (DiffractionData) data;
                int nres = diffractiondata.nresiduebfactor + 1;
                for (Atom a : atomArray) {
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
                                        && i < nXYZ + nB - 1) {
                                    i++;
                                }
                                if (i < nXYZ + nB) {
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
                logger.severe(" B refinement not supported for this data type!");
            }
        }

        if (refinementmode == RefinementMode.OCCUPANCIES
                || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            if (data instanceof DiffractionData) {
                int i = nXYZ + nB;
                for (ArrayList<Residue> list : refinementModel.altresidues) {
                    for (int j = 0; j < list.size(); j++) {
                        scaling[i] = occscale;
                        i++;
                    }
                }
                for (ArrayList<Molecule> list : refinementModel.altmolecules) {
                    for (int j = 0; j < list.size(); j++) {
                        scaling[i] = occscale;
                        i++;
                    }
                }
            } else {
                logger.severe(" Occupancy refinement not supported for this data type!");
            }
        }

        refinementEnergy.setScaling(scaling);
    }

    /**
     * <p>Getter for the field
     * <code>eps</code>.</p>
     *
     * @return a {@link java.lang.Double} object.
     */
    public Double getEps() {
        boolean hasaniso = false;

        for (Atom a : atomArray) {
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
                eps = 0.4;
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
                    eps = 0.2;
                }
                break;
            case OCCUPANCIES:
                eps = 0.1;
                break;
            case COORDINATES_AND_OCCUPANCIES:
                eps = 0.2;
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
                    eps = 0.2;
                }
                break;
        }

        return eps;
    }

    /**
     * get the number of xyz parameters being fit
     *
     * @return the number of xyz parameters
     */
    public int getNXYZ() {
        return nXYZ;
    }

    /**
     * get the number of B factor parameters being fit
     *
     * @return the number of B factor parameters
     */
    public int getNB() {
        return nB;
    }

    /**
     * get the number of occupancy parameters being fit
     *
     * @return the number of occupancy parameters
     */
    public int getNOcc() {
        return nOcc;
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
        if (dataContainer instanceof DiffractionData) {
            typestring = "X-ray";
        } else if (dataContainer instanceof RealSpaceData) {
            typestring = "Real Space";
        } else {
            typestring = "null";
        }

        logger.info(" Beginning " + typestring + " Refinement");
        switch (refinementMode) {
            case COORDINATES:
                logger.info(" Mode: Coordinates");
                break;
            case BFACTORS:
                logger.info(" Mode: B-Factors");
                break;
            case COORDINATES_AND_BFACTORS:
                logger.info(" Mode: Coordinates and B-Factors");
                break;
            case OCCUPANCIES:
                logger.info(" Mode: Occupancies");
            case COORDINATES_AND_OCCUPANCIES:
                logger.info(" Mode: Coordinates and Occupancies");
                break;
            case BFACTORS_AND_OCCUPANCIES:
                logger.info(" Mode: B-Factors and Occupancies");
                break;
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                logger.info(" Mode: Coordinates, B-Factors and Occupancies");
                break;
        }
        logger.info(" Number of Parameters: " + n);

        refinementEnergy.getCoordinates(x);
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
        double e = refinementEnergy.energyAndGradient(x, grad);
        status = LBFGS.minimize(n, m, x, e, grad, eps, maxiter, refinementEnergy, this);
        done = true;
        switch (status) {
            case 0:
                logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f", grms));
                break;
            case 1:
                logger.info(String.format("\n Optimization terminated at step %d.", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.");
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            mtime += System.nanoTime();
            sb.append(String.format(" Optimization time: %g (sec)", mtime * toSeconds));
            logger.info(sb.toString());
        }

        return refinementEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms,
            double xrms, double f, double df, double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        // update display
        if (listener != null) {
            MolecularAssembly molecularAssembly[];
            molecularAssembly = dataContainer.getMolecularAssembly();
            for (MolecularAssembly ma : molecularAssembly) {
                listener.algorithmUpdate(ma);
            }
        }

        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time      "
                    + dataContainer.printOptimizationHeader());
        }
        if (info == null) {
            logger.info(String.format("%6d %12.3f %10.3f",
                    iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("%6d %12.3f %10.3f %10.3f %9.4f %8.2f %6d %8.3f ",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
                sb.append(dataContainer.printOptimizationUpdate());
                logger.info(sb.toString());
            } else {
                logger.info(String.format("%6d %12.3g %10.3f %10.3f %9.4f %8.2f %6d %8s",
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

    /**
     * {@inheritDoc}
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating minimization.\n", e);
                }
            }
        }
    }
}
