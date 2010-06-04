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

import ffx.numerics.Optimizable;
import ffx.potential.bonded.Atom;
import ffx.xray.RefinementEnergy.RefinementMode;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class XRayEnergy implements Optimizable {

    private static final Logger logger = Logger.getLogger(XRayEnergy.class.getName());
    private final XRayStructure xraystructure;
    private final SigmaAMinimize sigmaaminimize;
    private final CrystalReciprocalSpace crs_fc;
    private final CrystalReciprocalSpace crs_fs;
    private final RefinementData refinementdata;
    private final Atom atomarray[];
    private final int nAtoms;
    private RefinementMode refinementMode;
    protected double[] optimizationScaling = null;

    public XRayEnergy(XRayStructure xraystructure, RefinementMode refinementmode) {
        this.xraystructure = xraystructure;
        this.refinementdata = xraystructure.refinementdata;
        this.crs_fc = refinementdata.crs_fc;
        this.crs_fs = refinementdata.crs_fs;
        this.sigmaaminimize = xraystructure.sigmaaminimize;
        this.refinementMode = refinementmode;
        this.atomarray = xraystructure.atomarray;
        this.nAtoms = atomarray.length;
    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double e = 0.0;
        int natoms = atomarray.length;
        switch (refinementMode) {
            case COORDINATES:
                for (int i = 0; i < natoms; i++) {
                    atomarray[i].setXYZGradient(0.0, 0.0, 0.0);
                }
                // update coordinates
                crs_fc.setCoordinates(x);
                crs_fs.setCoordinates(x);

                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);

                // compute crystal likelihood
                e = sigmaaminimize.calculateLikelihood();

                // compute the crystal gradients (requires inverse FFT)
                crs_fc.computeAtomicGradients(refinementdata.dfc,
                        refinementdata.freer, refinementdata.rfreeflag);
                crs_fs.computeAtomicGradients(refinementdata.dfs,
                        refinementdata.freer, refinementdata.rfreeflag);
                getXYZGradients(g);

                break;
            case BFACTORS:
                break;
            case COORDINATES_AND_BFACTORS:
                break;
            default:
                String message = "refinement mode not implemented.";
                logger.log(Level.SEVERE, message);
                break;
        }
        return e;
    }

    public void getXYZGradients(double g[]) {
        assert (g != null && g.length == nAtoms * 3);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : atomarray) {
            a.getXYZGradient(grad);
            g[index++] = grad[0];
            g[index++] = grad[1];
            g[index++] = grad[2];
        }
    }

    @Override
    public void setOptimizationScaling(double[] scaling) {
        if (scaling != null && scaling.length == nAtoms * 3) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getOptimizationScaling() {
        return optimizationScaling;
    }
}
