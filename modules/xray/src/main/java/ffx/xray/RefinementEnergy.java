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
import ffx.potential.bonded.MolecularAssembly;
import ffx.xray.RefinementMinimize.RefinementMode;
import java.util.Arrays;
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
    private final MolecularAssembly molecularAssembly;
    private final XRayStructure xraystructure;
    protected PotentialEnergy potentialEnergy;
    protected XRayEnergy xrayEnergy;
    private RefinementMode refinementMode;
    private final int nxyz;
    private final int nb;
    private final int nocc;
    private final int n;
    private double weight = 1.0;
    private double gChemical[];
    private double gXray[];

    public RefinementEnergy(MolecularAssembly molecularAssembly,
            XRayStructure xraystructure, int nxyz, int nb, int nocc,
            RefinementMode refinementmode) {
        this.molecularAssembly = molecularAssembly;
        this.xraystructure = xraystructure;
        this.refinementMode = refinementmode;
        this.nxyz = nxyz;
        this.nb = nb;
        this.nocc = nocc;
        this.n = nxyz + nb + nocc;
        potentialEnergy = molecularAssembly.getPotentialEnergy();
        if (potentialEnergy == null) {
            potentialEnergy = new PotentialEnergy(molecularAssembly);
            molecularAssembly.setPotential(potentialEnergy);
        }
        xrayEnergy = xraystructure.getXRayEnergy();
        if (xrayEnergy == null) {
            xrayEnergy = new XRayEnergy(xraystructure, nxyz, nb, nocc,
                    refinementMode);
            xraystructure.setXRayEnergy(xrayEnergy);
        }
    }

    /**
     * Implementation of the {@link Optimizable} interface for the RefinementEnergy.
     *
     * @param x
     * @param g
     * @return
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double e = 0.0;
        Arrays.fill(g, 0.0);
        switch (refinementMode) {
            case COORDINATES:
                // Compute the chemical energy and gradient.
                if (gChemical == null || gChemical.length != nxyz) {
                    gChemical = new double[nxyz];
                }
                e = potentialEnergy.energyAndGradient(x, gChemical);

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != nxyz) {
                    gXray = new double[nxyz];
                }
                e += weight * xrayEnergy.energyAndGradient(x, gXray);

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nxyz; i++) {
                    g[i] = gChemical[i] + weight * gXray[i];
                }
                break;
            case BFACTORS:
                // Compute the X-ray target energy and gradient.
                e = xrayEnergy.energyAndGradient(x, g);
                break;
            case COORDINATES_AND_BFACTORS:
                // Compute the chemical energy and gradient.
                if (gChemical == null || gChemical.length != n) {
                    gChemical = new double[n];
                }
                e = potentialEnergy.energyAndGradient(x, gChemical);

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != n) {
                    gXray = new double[n];
                }
                e += weight * xrayEnergy.energyAndGradient(x, gXray);

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nxyz; i++) {
                    g[i] = gChemical[i] + weight * gXray[i];
                }
                // bfactors, occ
                if (n > nxyz) {
                    for (int i = nxyz; i < n; i++) {
                        g[i] = weight * gXray[i];
                    }
                }
                break;
            default:
                String message = "Unknown refinment mode.";
                logger.log(Level.SEVERE, message);
        }
        return e;
    }

    @Override
    public void setOptimizationScaling(double[] scaling) {
        switch (refinementMode) {
            case COORDINATES:
            case BFACTORS:
            case COORDINATES_AND_BFACTORS:
                potentialEnergy.setOptimizationScaling(scaling);
                xrayEnergy.setOptimizationScaling(scaling);
                break;
            default:
                String message = "Unknown refinement mode.";
                logger.log(Level.SEVERE, message);
        }
    }

    @Override
    public double[] getOptimizationScaling() {
        switch (refinementMode) {
            case COORDINATES:
            case BFACTORS:
            case COORDINATES_AND_BFACTORS:
                return potentialEnergy.getOptimizationScaling();
            default:
                String message = "Unknown refinement mode.";
                logger.log(Level.SEVERE, message);
        }
        return null;
    }
}
