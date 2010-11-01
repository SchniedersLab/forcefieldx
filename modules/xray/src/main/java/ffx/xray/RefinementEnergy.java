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

import static ffx.numerics.VectorMath.b2u;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.RefinementMinimize.RefinementMode;
import java.util.ArrayList;
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
    private final RefinementData refinementdata;
    private final Atom[] atomarray;
    private final int nAtoms;
    private final List<Integer> xindex[];
    protected XRayEnergy xrayEnergy;
    private RefinementMode refinementMode;
    protected int nxyz;
    protected int nb;
    protected int nocc;
    private int n;
    private double weight;
    private double xChemical[][];
    private double gChemical[][];
    private double gXray[];
    protected double[] optimizationScaling = null;

    public RefinementEnergy(MolecularAssembly molecularAssembly,
            XRayStructure xraystructure, RefinementMode refinementmode) {
        this(new MolecularAssembly[]{molecularAssembly}, xraystructure,
                refinementmode, null);
    }

    public RefinementEnergy(MolecularAssembly molecularAssembly,
            XRayStructure xraystructure, RefinementMode refinementmode,
            double scaling[]) {
        this(new MolecularAssembly[]{molecularAssembly}, xraystructure,
                refinementmode, scaling);
    }

    public RefinementEnergy(MolecularAssembly molecularAssembly[],
            XRayStructure xraystructure, RefinementMode refinementmode,
            double scaling[]) {
        this.molecularAssembly = molecularAssembly;
        this.refinementdata = xraystructure.refinementdata;
        this.atomarray = xraystructure.atomarray;
        this.nAtoms = atomarray.length;
        this.xindex = xraystructure.xindex;
        this.weight = refinementdata.xweight;
        this.refinementMode = refinementmode;
        this.optimizationScaling = scaling;

        // determine size of fit
        n = nxyz = nb = nocc = 0;
        switch (refinementmode) {
            case COORDINATES:
                nxyz = nAtoms * 3;
                break;
            case COORDINATES_AND_BFACTORS:
            case COORDINATES_AND_OCCUPANCIES:
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                // coordinate params
                nxyz = nAtoms * 3;
            case BFACTORS:
            case OCCUPANCIES:
            case BFACTORS_AND_OCCUPANCIES:
                // bfactor params
                if (refinementmode == RefinementMode.BFACTORS
                        || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                        || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                        || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
                    int resnum = -1;
                    int nres = refinementdata.nresiduebfactor + 1;
                    for (Atom a : atomarray) {
                        // ignore hydrogens!!!
                        if (a.getAtomicNumber() == 1) {
                            continue;
                        }
                        if (a.getAnisou() == null) {
                            if (refinementdata.addanisou) {
                                double anisou[] = new double[6];
                                double u = b2u(a.getTempFactor());
                                anisou[0] = anisou[1] = anisou[2] = u;
                                anisou[3] = anisou[4] = anisou[5] = 0.0;
                                a.setAnisou(anisou);
                                nb += 6;
                            } else if (refinementdata.residuebfactor) {
                                if (resnum != a.getResidueNumber()) {
                                    if (nres >= refinementdata.nresiduebfactor) {
                                        nb++;
                                        nres = 1;
                                    } else {
                                        nres++;
                                    }
                                    resnum = a.getResidueNumber();
                                }
                            } else {
                                nb++;
                            }
                        } else {
                            nb += 6;
                        }
                    }
                    if (refinementdata.residuebfactor) {
                        if (nres < refinementdata.nresiduebfactor) {
                            nb--;
                        }
                    }
                }

                // occupancy params
                if (refinementmode == RefinementMode.OCCUPANCIES
                        || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                        || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                        || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
                    for (ArrayList<Residue> list : xraystructure.altresidues) {
                        nocc += list.size();
                    }
                    for (ArrayList<Molecule> list : xraystructure.altmolecules) {
                        nocc += list.size();
                    }
                }
                break;
        }
        n = nxyz + nb + nocc;

        // initialize force field and Xray energies
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
            case OCCUPANCIES:
            case BFACTORS_AND_OCCUPANCIES:
                // Compute the X-ray target energy and gradient.
                e = xrayEnergy.energyAndGradient(x, g);
                break;
            case COORDINATES_AND_BFACTORS:
            case COORDINATES_AND_OCCUPANCIES:
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
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
                String message = "Unknown refinement mode.";
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

    public double getXWeight() {
        return this.weight;
    }

    public void setXWeight(double weight) {
        this.weight = weight;
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
        return xrayEnergy.getNumberOfVariables();
    }

    @Override
    public double[] getCoordinates(double[] parameters) {
        return xrayEnergy.getCoordinates(parameters);
    }
}
