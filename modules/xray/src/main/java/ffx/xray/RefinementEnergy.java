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

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.MolecularAssembly;
import ffx.xray.RefinementMinimize.RefinementMode;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class RefinementEnergy implements Potential {

    private static final Logger logger = Logger.getLogger(RefinementEnergy.class.getName());
    private final MolecularAssembly molecularAssembly[];
    private final XRayStructure xraystructure;
    private final List<Integer> xindex[];
    protected XRayEnergy xrayEnergy;
    private RefinementMode refinementMode;
    private final int nxyz;
    private final int nb;
    private final int nocc;
    private final int n;
    private double weight;
    private double xChemical[][];
    private double gChemical[][];
    private double gXray[];
    protected double[] optimizationScaling = null;

    public RefinementEnergy(MolecularAssembly molecularAssembly[],
            XRayStructure xraystructure, int nxyz, int nb, int nocc,
            RefinementMode refinementmode) {
        this.molecularAssembly = molecularAssembly;
        this.xraystructure = xraystructure;
        this.weight = xraystructure.refinementdata.xweight;
        xindex = xraystructure.xindex;
        this.refinementMode = refinementmode;
        this.nxyz = nxyz;
        this.nb = nb;
        this.nocc = nocc;
        this.n = nxyz + nb + nocc;

        for (MolecularAssembly ma : molecularAssembly) {
            ForceFieldEnergy fe = ma.getPotentialEnergy();
            if (fe == null) {
                fe = new ForceFieldEnergy(ma);
                ma.setPotential(fe);
            }
            fe.setScaling(null);
        }
        if (!xraystructure.scaled) {
            xraystructure.scalebulkfit();
            xraystructure.printstats();
        }
        xrayEnergy = xraystructure.getXRayEnergy();
        if (xrayEnergy == null) {
            xrayEnergy = new XRayEnergy(xraystructure, nxyz, nb, nocc,
                    refinementMode);
            xraystructure.setXRayEnergy(xrayEnergy);
        } else {
            xrayEnergy.setNXYZ(nxyz);
            xrayEnergy.setNB(nb);
            xrayEnergy.setNOcc(nocc);
            xrayEnergy.setRefinementMode(refinementmode);
        }
        xrayEnergy.setScaling(null);

        int assemblysize = molecularAssembly.length;
        xChemical = new double[assemblysize][];
        gChemical = new double[assemblysize][];
        for (int i = 0; i < assemblysize; i++) {
            int len = molecularAssembly[i].getAtomArray().length * 3;
            xChemical[i] = new double[len];
            gChemical[i] = new double[len];
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

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        int assemblysize = molecularAssembly.length;
        switch (refinementMode) {
            case COORDINATES:
                // Compute the chemical energy and gradient.
                for (int i = 0; i < assemblysize; i++) {
                    ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
                    getAssemblyi(i, x, xChemical[i]);
                    double curE = fe.energyAndGradient(xChemical[i], gChemical[i]);
                    e += (curE - e) / (i + 1);
                    setAssemblyi(i, g, gChemical[i]);
                }
                // normalize gradients for multiple-counted atoms
                if (assemblysize > 1) {
                    for (int i = 0; i < nxyz; i++) {
                        g[i] /= assemblysize;
                    }
                }

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != nxyz) {
                    gXray = new double[nxyz];
                }
                e += weight * xrayEnergy.energyAndGradient(x, gXray);

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nxyz; i++) {
                    g[i] += weight * gXray[i];
                }
                break;
            case BFACTORS:
                // Compute the X-ray target energy and gradient.
                e = xrayEnergy.energyAndGradient(x, g);
                break;
            case COORDINATES_AND_BFACTORS:
                // Compute the chemical energy and gradient.
                for (int i = 0; i < assemblysize; i++) {
                    ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
                    getAssemblyi(i, x, xChemical[i]);
                    double curE = fe.energyAndGradient(xChemical[i], gChemical[i]);
                    e += (curE - e) / (i + 1);
                    setAssemblyi(i, g, gChemical[i]);
                }
                // normalize gradients for multiple-counted atoms
                if (assemblysize > 1) {
                    for (int i = 0; i < nxyz; i++) {
                        g[i] /= assemblysize;
                    }
                }

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != n) {
                    gXray = new double[n];
                }
                e += weight * xrayEnergy.energyAndGradient(x, gXray);

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nxyz; i++) {
                    g[i] += weight * gXray[i];
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

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }
        return e;
    }

    public void getAssemblyi(int i, double x[], double xchem[]) {
        assert (x != null && xchem != null);
        for (int j = 0; j < xchem.length; j += 3) {
            int index = (j + 1) / 3;
            int aindex = xindex[i].get(index) * 3;
            xchem[j] = x[aindex];
            xchem[j + 1] = x[aindex + 1];
            xchem[j + 2] = x[aindex + 2];
        }
    }

    public void setAssemblyi(int i, double x[], double xchem[]) {
        assert (x != null && xchem != null);
        for (int j = 0; j < xchem.length; j += 3) {
            int index = (j + 1) / 3;
            int aindex = xindex[i].get(index) * 3;
            x[aindex] += xchem[j];
            x[aindex + 1] += xchem[j + 1];
            x[aindex + 2] += xchem[j + 2];
        }
    }

    @Override
    public void setScaling(double[] scaling) {
        optimizationScaling = scaling;
    }

    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    @Override
    public double[] getMass() {
        return xrayEnergy.getMass();
    }

    @Override
    public int getNumberOfVariables() {
        return n;
    }

    @Override
    public double[] getCoordinates(double[] parameters) {
        return xrayEnergy.getCoordinates(parameters);
    }
}
