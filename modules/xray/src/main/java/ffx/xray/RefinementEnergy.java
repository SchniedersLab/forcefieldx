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
import ffx.potential.PotentialEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class RefinementEnergy implements Optimizable {

    private static final Logger logger = Logger.getLogger(RefinementEnergy.class.getName());

    enum RefinementMode {

        COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS
    }
    private final MolecularAssembly molecularAssembly;
    private final XRayStructure xraystructure;
    private final SigmaAMinimize sigmaaminimize;
    private final CrystalReciprocalSpace crs_fc;
    private final CrystalReciprocalSpace crs_fs;
    private final RefinementData refinementdata;
    private final Atom atomArray[];
    private PotentialEnergy potentialEnergy;
    private RefinementMode refinementMode;
    private double weight = 1.0;

    public RefinementEnergy(MolecularAssembly molecularAssembly,
            XRayStructure xraystructure) {
        this.molecularAssembly = molecularAssembly;
        this.xraystructure = xraystructure;
        this.refinementdata = xraystructure.refinementdata;
        this.crs_fc = refinementdata.crs_fc;
        this.crs_fs = refinementdata.crs_fs;
        this.sigmaaminimize = xraystructure.sigmaaminimize;
        this.refinementMode = RefinementMode.COORDINATES;
        this.atomArray = molecularAssembly.getAtomArray();
        potentialEnergy = new PotentialEnergy(molecularAssembly);
    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double e = 0.0;
        switch (refinementMode) {
            case COORDINATES:
            case BFACTORS:
                // Compute the energy and gradient for the chemical term.
                e = potentialEnergy.energyAndGradient(x, g);
                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);
                // update sigmaA
                // sigmaaminimize.minimize(7, 1e-1);
                // e += refinementdata.llkr;
                // compute crystal likelihood
                e += sigmaaminimize.calculateLikelihood();
                // compute the crystal gradients (requires inverse FFT)
                crs_fc.computeAtomicGradients(refinementdata.dfc,
                        refinementdata.freer, refinementdata.rfreeflag);
                if (crs_fs.solventmodel == SolventModel.GAUSSIAN
                        || crs_fs.solventmodel == SolventModel.POLYNOMIAL) {
                    crs_fs.computeAtomicGradients(refinementdata.dfs,
                            refinementdata.freer, refinementdata.rfreeflag);
                }
                return e;
            default:
                String message = "Joint coordinate + bfactor refinement is not implemented.";
                logger.log(Level.SEVERE, message);
        }
        return e;
    }

    @Override
    public void setOptimizationScaling(double[] scaling) {
        switch (refinementMode) {
            case COORDINATES:
                potentialEnergy.setOptimizationScaling(scaling);
                break;
            default:
                String message = "Only coordinate refinement is implemented.";
                logger.log(Level.SEVERE, message);
        }
    }

    @Override
    public double[] getOptimizationScaling() {
        switch (refinementMode) {
            case COORDINATES:
                return potentialEnergy.getOptimizationScaling();
            default:
                String message = "Only coordinate refinement is implemented.";
                logger.log(Level.SEVERE, message);
        }
        return null;
    }
}
