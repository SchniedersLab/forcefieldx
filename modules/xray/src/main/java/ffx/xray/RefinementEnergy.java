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
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Thermostat;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.numerics.VectorMath.b2u;
import static ffx.xray.RefinementMinimize.RefinementMode.BFACTORS;
import static ffx.xray.RefinementMinimize.RefinementMode.BFACTORS_AND_OCCUPANCIES;
import static ffx.xray.RefinementMinimize.RefinementMode.COORDINATES;
import static ffx.xray.RefinementMinimize.RefinementMode.COORDINATES_AND_BFACTORS;
import static ffx.xray.RefinementMinimize.RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES;
import static ffx.xray.RefinementMinimize.RefinementMode.COORDINATES_AND_OCCUPANCIES;
import static ffx.xray.RefinementMinimize.RefinementMode.OCCUPANCIES;

/**
 * Combine the X-ray target and chemical potential energy using the
 * {@link Potential} interface
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 * @since 1.0
 *
 */
public class RefinementEnergy implements LambdaInterface, Potential, AlgorithmListener {

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
    private double ktscale;
    private double xChemical[][];
    private double gChemical[][];
    private double gXray[];
    private double totalEnergy;
    protected double[] optimizationScaling = null;

    /**
     * constructor for energy
     *
     * @param data input {@link DiffractionData data} for refinement
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for
     * refinement
     */
    public RefinementEnergy(DataContainer data, RefinementMode refinementmode) {
        this(data, refinementmode, null);
    }

    /**
     * constructor for energy
     *
     * @param data input {@link DiffractionData data} for refinement
     * @param refinementmode {@link RefinementMinimize.RefinementMode} for
     * refinement
     * @param scaling scaling of refinement parameters
     */
    public RefinementEnergy(DataContainer data,
            RefinementMode refinementmode, double scaling[]) {
        this.molecularAssembly = data.getMolecularAssembly();
        this.data = data;
        this.refinementmodel = data.getRefinementModel();
        this.atomarray = data.getAtomArray();
        this.nAtoms = atomarray.length;
        this.xindex = refinementmodel.xIndex;
        this.refinementMode = refinementmode;
        this.optimizationScaling = scaling;
        this.thermostat = null;
        this.ktscale = 1.0;

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
                        for (ArrayList<Residue> list : refinementmodel.altResidues) {
                            nocc += list.size();
                        }
                        for (ArrayList<Molecule> list : refinementmodel.altMolecules) {
                            nocc += list.size();
                        }
                    }
                } else {
                    logger.severe(" Refinement method not supported for this data type!");
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

    public Potential getDataEnergy() {
        return dataEnergy;
    }

    @Override
    public double energy(double[] x) {
        double weight = data.getWeight();
        double e = 0.0;

        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
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
                    double curE = fe.energy(xChemical[i]);
                    e += (curE - e) / (i + 1);
                }
                double chemE = e;

                e = chemE * ktscale;

                // Compute the X-ray target energy.
                double xE = dataEnergy.energy(x);
                e += weight * xE;
                break;
            case BFACTORS:
            case OCCUPANCIES:
            case BFACTORS_AND_OCCUPANCIES:
                // Compute the X-ray target energy and gradient.
                e = dataEnergy.energy(x);
                break;
            case COORDINATES_AND_BFACTORS:
            case COORDINATES_AND_OCCUPANCIES:
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                // Compute the chemical energy and gradient.
                for (int i = 0; i < assemblysize; i++) {
                    ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
                    getAssemblyi(i, x, xChemical[i]);
                    double curE = fe.energy(xChemical[i]);
                    e += (curE - e) / (i + 1);
                }
                e += weight * dataEnergy.energy(x);
                break;
            default:
                String message = "Unknown refinement mode.";
                logger.log(Level.SEVERE, message);
        }

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
            }
        }
        totalEnergy = e;
        return e;
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of the {@link Potential} interface for the
     * RefinementEnergy.
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double weight = data.getWeight();
        double e = 0.0;
        fill(g, 0.0);

        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
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
                double chemE = e;
                //System.out.println("chem E: " + e + " scaled chem E: " + ktscale * e);

                e = chemE * ktscale;
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
                double xE = dataEnergy.energyAndGradient(x, gXray);
                //System.out.println("Xray E: " + xE + " scaled Xray E: " + weight * xE);
                e += weight * xE;

                /*
                 * double normchem = 0.0; for (int i = 0; i < nxyz; i++) {
                 * normchem += g[i] * g[i]; } normchem = Math.sqrt(normchem) /
                 * nxyz; double normxray = 0.0; for (int i = 0; i < nxyz; i++) {
                 * normxray += gXray[i] * gXray[i]; } normxray =
                 * Math.sqrt(normxray) / nxyz; System.out.println("chem: " +
                 * normchem + " xray: " + normxray + " weight wa: " + normchem /
                 * normxray);
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
        totalEnergy = e;
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
        return data.getWeight();
    }

    /**
     * get the current kT scaling weight
     *
     * @return kT scale
     */
    public double getKTScale() {
        return this.ktscale;
    }

    /**
     * set the current kT scaling weight
     *
     * @param ktscale requested kT scale
     */
    public void setKTScale(double ktscale) {
        this.ktscale = ktscale;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double[] scaling) {
        optimizationScaling = scaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        return dataEnergy.getMass();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        return dataEnergy.getNumberOfVariables();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        return dataEnergy.getCoordinates(parameters);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean algorithmUpdate(MolecularAssembly active) {
        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }
        logger.info(" kTscale: " + ktscale);
        logger.info(data.printEnergyUpdate());

        return true;
    }

    // this should probably be part of the potential class
    /**
     * <p>
     * Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link ffx.algorithms.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    // this should probably be part of the potential class
    /**
     * <p>
     * Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        for (MolecularAssembly ma : molecularAssembly) {
            ForceFieldEnergy fe = ma.getPotentialEnergy();
            fe.setLambda(lambda);
        }

        if (data instanceof DiffractionData) {
            XRayEnergy xrayenergy = (XRayEnergy) dataEnergy;
            xrayenergy.setLambda(lambda);
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realspaceenergy = (RealSpaceEnergy) dataEnergy;
            realspaceenergy.setLambda(lambda);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getLambda() {
        double l = 1.0;
        if (data instanceof DiffractionData) {
            XRayEnergy xrayenergy = (XRayEnergy) dataEnergy;
            l = xrayenergy.getLambda();
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realspaceenergy = (RealSpaceEnergy) dataEnergy;
            l = realspaceenergy.getLambda();
        }

        return l;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        double e = 0.0;

        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }

        int assemblysize = molecularAssembly.length;
        // Compute the chemical energy and gradient.
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
            double curE = fe.getdEdL();
            e += (curE - e) / (i + 1);
        }
        e *= ktscale;

        if (data instanceof DiffractionData) {
            XRayEnergy xrayenergy = (XRayEnergy) dataEnergy;
            e += xrayenergy.getdEdL();
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realspaceenergy = (RealSpaceEnergy) dataEnergy;
            e += realspaceenergy.getdEdL();
        }

        return e;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        double e = 0.0;

        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }

        int assemblysize = molecularAssembly.length;
        // Compute the chemical energy and gradient.
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
            double curE = fe.getd2EdL2();
            e += (curE - e) / (i + 1);
        }
        e *= ktscale;

        return e;
    }

    /*
     * FIXME: needs to handle multiple conformations
     */
    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        double weight = data.getWeight();

        if (thermostat != null) {
            ktscale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }

        int assemblysize = molecularAssembly.length;
        // Compute the chemical energy and gradient.
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
            fe.getdEdXdL(gChemical[i]);
        }
        for (int i = 0; i < assemblysize; i++) {
            for (int j = 0; j < nxyz; j++) {
                gradient[j] += gChemical[i][j];
            }
        }
        // normalize gradients for multiple-counted atoms
        if (assemblysize > 1) {
            for (int i = 0; i < nxyz; i++) {
                gradient[i] /= assemblysize;
            }
        }
        for (int i = 0; i < nxyz; i++) {
            gradient[i] *= ktscale;
        }

        // clear gradients for X-ray calculation
        for (Atom a : atomarray) {
            a.setXYZGradient(0.0, 0.0, 0.0);
        }

        // Compute the X-ray target energy and gradient.
        if (gXray == null || gXray.length != nxyz) {
            gXray = new double[nxyz];
        }
        if (data instanceof DiffractionData) {
            XRayEnergy xrayenergy = (XRayEnergy) dataEnergy;
            xrayenergy.getdEdXdL(gXray);
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realspaceenergy = (RealSpaceEnergy) dataEnergy;
            realspaceenergy.getdEdXdL(gXray);
        }

        // Add the chemical and X-ray gradients.
        for (int i = 0; i < nxyz; i++) {
            gradient[i] += weight * gXray[i];
        }
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return dataEnergy.getVariableTypes();
    }

    @Override
    public void setEnergyTermState(STATE state) {
        int assemblysize = molecularAssembly.length;
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy fe = molecularAssembly[i].getPotentialEnergy();
            fe.setEnergyTermState(state);
        }
        dataEnergy.setEnergyTermState(state);
    }
}
