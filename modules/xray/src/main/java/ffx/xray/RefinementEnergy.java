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

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Thermostat;
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
 * Combine the X-ray target and chemical potential energy using the
 * {@link Potential} interface
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class RefinementEnergy implements Potential, AlgorithmListener {

    private static final Logger logger = Logger.getLogger(RefinementEnergy.class.getName());
    private final MolecularAssembly molecularAssembly[];
    private final DataContainer data;
    private final RefinementModel refinementmodel;
    private final Atom[] atomarray;
    private final int nAtoms;
    private final List<Integer> xindex[];
    protected Potential dataEnergy;
    protected Thermostat thermostat;
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

    /**
     * constructor for energy
     *
     * @param data input {@link DiffractionData data} for refinement
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for refinement
     */
    public RefinementEnergy(DataContainer data, RefinementMode refinementmode) {
        this(data, refinementmode, null);
    }

    /**
     * constructor for energy
     *
     * @param data input {@link DiffractionData data} for refinement
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for refinement
     * @param scaling scaling of refinement parameters
     */
    public RefinementEnergy(DataContainer data,
            RefinementMode refinementmode, double scaling[]) {
        this.molecularAssembly = data.getMolecularAssembly();
        this.data = data;
        this.refinementmodel = data.getRefinementModel();
        this.atomarray = data.getAtomArray();
        this.nAtoms = atomarray.length;
        this.xindex = refinementmodel.xindex;
        this.weight = data.getWeight();
        this.refinementMode = refinementmode;
        this.optimizationScaling = scaling;
        this.thermostat = null;

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
                if (data instanceof DiffractionData) {
                    DiffractionData diffractiondata = (DiffractionData) data;
                    // bfactor params
                    if (refinementmode == RefinementMode.BFACTORS
                            || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                            || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS
                            || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
                        int resnum = -1;
                        int nres = diffractiondata.nresiduebfactor + 1;
                        for (Atom a : atomarray) {
                            // ignore hydrogens!!!
                            if (a.getAtomicNumber() == 1) {
                                continue;
                            }
                            if (a.getAnisou() == null) {
                                if (diffractiondata.addanisou) {
                                    double anisou[] = new double[6];
                                    double u = b2u(a.getTempFactor());
                                    anisou[0] = anisou[1] = anisou[2] = u;
                                    anisou[3] = anisou[4] = anisou[5] = 0.0;
                                    a.setAnisou(anisou);
                                    nb += 6;
                                } else if (diffractiondata.residuebfactor) {
                                    if (resnum != a.getResidueNumber()) {
                                        if (nres >= diffractiondata.nresiduebfactor) {
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
                        if (diffractiondata.residuebfactor) {
                            if (nres < diffractiondata.nresiduebfactor) {
                                nb--;
                            }
                        }
                    }

                    // occupancy params
                    if (refinementmode == RefinementMode.OCCUPANCIES
                            || refinementmode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                            || refinementmode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                            || refinementmode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
                        for (ArrayList<Residue> list : refinementmodel.altresidues) {
                            nocc += list.size();
                        }
                        for (ArrayList<Molecule> list : refinementmodel.altmolecules) {
                            nocc += list.size();
                        }
                    }
                } else {
                    logger.severe("Refinement method not supported for this data type!");
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

        if (data instanceof DiffractionData) {
            DiffractionData diffractiondata = (DiffractionData) data;
            if (!diffractiondata.scaled[0]) {
                diffractiondata.printStats();
            }

            dataEnergy = new XRayEnergy(diffractiondata, nxyz, nb, nocc,
                    refinementMode);
            dataEnergy.setScaling(null);
        } else if (data instanceof RealSpaceData) {
            RealSpaceData realspacedata = (RealSpaceData) data;
            dataEnergy = new RealSpaceEnergy(realspacedata, nxyz, 0, 0,
                    refinementMode);
            dataEnergy.setScaling(null);
        }

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
     * Implementation of the {@link Potential} interface for the RefinementEnergy.
     *
     * @param x input params
     * @param g output gradients
     * @return energy
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double e = 0.0;
        Arrays.fill(g, 0.0);

        double ktscale = 1.0;
        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getCurrentTemperture() * Thermostat.kB);
        }

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
                e *= ktscale;
                // normalize gradients for multiple-counted atoms
                if (assemblysize > 1) {
                    for (int i = 0; i < nxyz; i++) {
                        g[i] /= assemblysize;
                    }
                }
                for (int i = 0; i < nxyz; i++) {
                    g[i] *= ktscale;
                }

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != nxyz) {
                    gXray = new double[nxyz];
                }
                e += weight * dataEnergy.energyAndGradient(x, gXray);

                /*
                double normchem = 0.0;
                for (int i = 0; i < nxyz; i++) {
                normchem += g[i] * g[i];
                }
                normchem = Math.sqrt(normchem) / nxyz;
                double normxray = 0.0;
                for (int i = 0; i < nxyz; i++) {
                normxray += gXray[i] * gXray[i];
                }
                normxray = Math.sqrt(normxray) / nxyz;
                System.out.println("chem: " + normchem + " xray: " + normxray + " weight wa: " + normchem / normxray);
                 */

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nxyz; i++) {
                    g[i] += weight * gXray[i];
                }
                break;
            case BFACTORS:
            case OCCUPANCIES:
            case BFACTORS_AND_OCCUPANCIES:
                // Compute the X-ray target energy and gradient.
                e = dataEnergy.energyAndGradient(x, g);
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
                e += weight * dataEnergy.energyAndGradient(x, gXray);

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

    /**
     * get the molecular assembly associated with index i of n, put in xchem
     *
     * @param i the desired molecular assembly index to "set" xchem to
     * @param x all parameters
     * @param xchem the xchem parameters for the particular molecular assembly
     * that will be passed to {@link ForceFieldEnergy}
     */
    public void getAssemblyi(int i, double x[], double xchem[]) {
        assert (x != null && xchem != null);
        for (int j = 0; j < xchem.length; j += 3) {
            int index = j / 3;
            int aindex = xindex[i].get(index) * 3;
            xchem[j] = x[aindex];
            xchem[j + 1] = x[aindex + 1];
            xchem[j + 2] = x[aindex + 2];
        }
    }

    /**
     * get the molecular assembly associated with index i of n, put in x
     *
     * @param i the desired molecular assembly index to "set" x to
     * @param x all parameters
     * @param xchem the xchem parameters for the particular molecular assembly
     * that will be passed to {@link ForceFieldEnergy}
     */
    public void setAssemblyi(int i, double x[], double xchem[]) {
        assert (x != null && xchem != null);
        for (int j = 0; j < xchem.length; j += 3) {
            int index = j / 3;
            int aindex = xindex[i].get(index) * 3;
            x[aindex] += xchem[j];
            x[aindex + 1] += xchem[j + 1];
            x[aindex + 2] += xchem[j + 2];
        }
    }

    /**
     * get the current dataset weight
     *
     * @return weight wA
     */
    public double getXWeight() {
        return this.weight;
    }

    /**
     * set the current dataset weight
     *
     * @param weight requested weight wA
     */
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
        return dataEnergy.getMass();
    }

    @Override
    public int getNumberOfVariables() {
        return dataEnergy.getNumberOfVariables();
    }

    @Override
    public double[] getCoordinates(double[] parameters) {
        return dataEnergy.getCoordinates(parameters);
    }

    @Override
    public boolean algorithmUpdate(MolecularAssembly active) {
        double ktscale = 1.0;
        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getCurrentTemperture() * Thermostat.kB);
        }
        logger.info("kTscale: " + ktscale);
        logger.info(data.printEnergyUpdate());

        return true;
    }

    // this should probably be part of the potential class
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    // this should probably be part of the potential class
    public Thermostat getThermostat() {
        return thermostat;
    }
}
