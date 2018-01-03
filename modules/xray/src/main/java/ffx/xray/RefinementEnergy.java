/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.xray;

import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Thermostat;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Descriptions;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import ffx.realspace.RealSpaceData;
import ffx.realspace.RealSpaceEnergy;
import ffx.xray.RefinementMinimize.RefinementMode;
import static ffx.numerics.VectorMath.b2u;

/**
 * Combine the X-ray target and chemical potential energy using the
 * {@link Potential} interface
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class RefinementEnergy implements LambdaInterface, Potential, AlgorithmListener {

    private static final Logger logger = Logger.getLogger(RefinementEnergy.class.getName());

    /**
     * MolecularAssembly instances being refined.
     */
    private final MolecularAssembly molecularAssemblies[];
    /**
     * Compute fast varying forces, slowly varying forces, or both.
     */
    private STATE state = STATE.BOTH;
    /**
     * Container to huge experimental data.
     */
    private final DataContainer data;
    /**
     * The refinement model.
     */
    private final RefinementModel refinementModel;
    /**
     * An array of atoms being refined.
     */
    private final Atom[] atomArray;
    /**
     * The number of atoms being refined.
     */
    private final int nAtoms;
    /**
     * An array of active atoms.
     */
    private final Atom[] activeAtomArray;
    /**
     * The number of active atoms.
     */
    private final int nActive;
    /**
     * An array of XYZIndex values.
     */
    private final List<Integer> xIndex[];
    /**
     * The Potential based on experimental data.
     */
    protected Potential dataEnergy;
    /**
     * A thermostat instance.
     */
    protected Thermostat thermostat;
    /**
     * The refinement mode being used.
     */
    private RefinementMode refinementMode;
    /**
     * The number of XYZ coordinates being refined.
     */
    protected int nXYZ;
    /**
     * The number of b-factor parameters being refined.
     */
    protected int nBFactor;
    /**
     * The number of occupancy parameters being refined.
     */
    protected int nOccupancy;
    /**
     * The total number of parameters being refined.
     */
    private int n;
    /**
     * The kT scale factor.
     */
    private double kTScale;
    /**
     * Atomic coordinates for computing the chemical energy.
     */
    private double xChemical[][];
    /**
     * Array for storing chemical gradient.
     */
    private double gChemical[][];
    /**
     * Array for storing the experimental gradient.
     */
    private double gXray[];
    /**
     * Total potential energy.
     */
    private double totalEnergy;
    /**
     * Optimization scale factors.
     */
    protected double[] optimizationScaling = null;
    /**
     * If true, collect lambda derivatives.
     */
    protected boolean lambdaTerm;
    /**
     * Print a file if there is an error in the energy.
     */
    private boolean printOnFailure;

    /**
     * RefinementEnergy Constructor.
     *
     * @param data input {@link DiffractionData data} for refinement
     * @param refinementMode {@link RefinementMinimize.RefinementMode} for
     * refinement
     */
    public RefinementEnergy(DataContainer data, RefinementMode refinementMode) {
        this(data, refinementMode, null);
    }

    /**
     * RefinementEnergy Constructor.
     *
     * @param data input {@link DiffractionData data} for refinement
     * @param refinementMode {@link RefinementMinimize.RefinementMode} for
     * refinement
     * @param optimizationScaling scaling of refinement parameters
     */
    public RefinementEnergy(DataContainer data,
            RefinementMode refinementMode, double optimizationScaling[]) {

        this.data = data;
        this.refinementMode = refinementMode;
        this.optimizationScaling = optimizationScaling;
        molecularAssemblies = data.getMolecularAssemblies();
        refinementModel = data.getRefinementModel();
        atomArray = data.getAtomArray();
        nAtoms = atomArray.length;

        // Determine if lambda derivatives are needed.
        ForceField forceField = molecularAssemblies[0].getForceField();
        lambdaTerm = forceField.getBoolean(ForceFieldBoolean.LAMBDATERM, false);
        printOnFailure = forceField.getBoolean(ForceFieldBoolean.PRINT_ON_FAILURE, true);

        // Fill an active atom array.
        int count = 0;
        int nUse = 0;
        for (Atom a : atomArray) {
            if (a.isActive()) {
                count++;
            }
            if (a.getUse()) {
                nUse++;
            }
        }
        nActive = count;
        activeAtomArray = new Atom[count];
        count = 0;
        for (Atom a : atomArray) {
            if (a.isActive()) {
                activeAtomArray[count++] = a;
            }
        }

        xIndex = refinementModel.getxIndex();
        thermostat = null;
        kTScale = 1.0;

        // determine size of fit
        n = nXYZ = nBFactor = nOccupancy = 0;
        switch (refinementMode) {
            case COORDINATES:
                nXYZ = nActive * 3;
                break;
            case COORDINATES_AND_BFACTORS:
            case COORDINATES_AND_OCCUPANCIES:
            case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES:
                // coordinate params
                nXYZ = nActive * 3;
            case BFACTORS:
            case OCCUPANCIES:
            case BFACTORS_AND_OCCUPANCIES:
                if (data instanceof DiffractionData) {
                    DiffractionData diffractionData = (DiffractionData) data;
                    // bfactor params
                    if (refinementMode == RefinementMode.BFACTORS
                            || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                            || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                            || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
                        int resnum = -1;
                        int nres = diffractionData.getnResidueBFactor() + 1;
                        for (Atom a : atomArray) {
                            // Ignore hydrogens and atoms that are not active.
                            if (a.getAtomicNumber() == 1 || !a.isActive()) {
                                continue;
                            }
                            if (a.getAnisou(null) == null) {
                                if (diffractionData.isAddAnisou()) {
                                    logger.info(format(" Adding ANISOU to %s.", a.describe(Descriptions.Resnum_Name)));
                                    double anisou[] = new double[6];
                                    double u = b2u(a.getTempFactor());
                                    anisou[0] = anisou[1] = anisou[2] = u;
                                    anisou[3] = anisou[4] = anisou[5] = 0.0;
                                    a.setAnisou(anisou);
                                    nBFactor += 6;
                                } else if (diffractionData.isResidueBFactor()) {
                                    if (resnum != a.getResidueNumber()) {
                                        if (nres >= diffractionData.getnResidueBFactor()) {
                                            nBFactor++;
                                            nres = 1;
                                        } else {
                                            nres++;
                                        }
                                        resnum = a.getResidueNumber();
                                    }
                                } else {
                                    nBFactor++;
                                }
                            } else {
                                nBFactor += 6;
                            }
                        }
                        if (diffractionData.isResidueBFactor()) {
                            if (nres < diffractionData.getnResidueBFactor()) {
                                nBFactor--;
                            }
                        }
                    }

                    // occupancy params
                    if (refinementMode == RefinementMode.OCCUPANCIES
                            || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                            || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                            || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
                        for (ArrayList<Residue> list : refinementModel.getAltResidues()) {
                            nOccupancy += list.size();
                        }
                        for (ArrayList<Molecule> list : refinementModel.getAltMolecules()) {
                            nOccupancy += list.size();
                        }
                        if (nActive != nAtoms) {
                            logger.severe(" Occupancy refinement is not supported with inactive atoms.");
                        }
                    }
                } else {
                    logger.severe(" Refinement method not supported for this data type!");
                }
                break;
        }

        logger.info(String.format("\n RefinementEnergy\n  Number of atoms:\t\t%d\n  Atoms being used:  \t\t%d\n  Active atoms: \t\t%d",
                nAtoms, nUse, nActive));

        n = nXYZ + nBFactor + nOccupancy;
        logger.info(String.format("  Number of variables:\t\t%d (nXYZ %d, nB %d, nOcc %d)\n",
                n, nXYZ, nBFactor, nOccupancy));

        // initialize force field and Xray energies
        for (MolecularAssembly molecularAssembly : molecularAssemblies) {
            ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
            if (forceFieldEnergy == null) {
                forceFieldEnergy = ForceFieldEnergy.energyFactory(molecularAssembly);
                molecularAssembly.setPotential(forceFieldEnergy);
            }
            forceFieldEnergy.setScaling(null);
        }

        if (data instanceof DiffractionData) {
            DiffractionData diffractionData = (DiffractionData) data;
            if (!diffractionData.getScaled()[0]) {
                diffractionData.printStats();
            }
            dataEnergy = new XRayEnergy(diffractionData, nXYZ, nBFactor, nOccupancy,
                    refinementMode);
            dataEnergy.setScaling(null);
        } else if (data instanceof RealSpaceData) {
            RealSpaceData realSpaceData = (RealSpaceData) data;
            dataEnergy = new RealSpaceEnergy(realSpaceData, nXYZ, 0, 0,
                    refinementMode);
            dataEnergy.setScaling(null);
        }

        int assemblySize = molecularAssemblies.length;
        xChemical = new double[assemblySize][];
        gChemical = new double[assemblySize][];
        for (int i = 0; i < assemblySize; i++) {
            int len = molecularAssemblies[i].getActiveAtomArray().length * 3;
            xChemical[i] = new double[len];
            gChemical[i] = new double[len];
        }
    }

    public Potential getDataEnergy() {
        return dataEnergy;
    }

    public Atom[] getActiveAtoms() {
        return activeAtomArray;
    }

    /**
     * Sets the printOnFailure flag; if override is true, over-rides any
     * existing property. Essentially sets the default value of printOnFailure
     * for an algorithm. For example, rotamer optimization will generally run
     * into force field issues in the normal course of execution as it tries
     * unphysical self and pair configurations, so the algorithm should not
     * print out a large number of error PDBs.
     *
     * @param onFail To set
     * @param override Override properties
     */
    public void setPrintOnFailure(boolean onFail, boolean override) {
        if (override) {
            // Ignore any pre-existing value
            printOnFailure = onFail;
        } else {
            try {
                molecularAssemblies[0].getForceField().getBoolean(ForceFieldBoolean.PRINT_ON_FAILURE);
                /*
                 * If the call was successful, the property was explicitly set
                 * somewhere and should be kept. If an exception was thrown, the
                 * property was never set explicitly, so over-write.
                 */
            } catch (Exception ex) {
                printOnFailure = onFail;
            }
        }
    }

    public boolean getPrintOnFailure() {
        return printOnFailure;
    }

    @Override
    public double energy(double[] x) {
        double weight = data.getWeight();
        double e = 0.0;

        if (thermostat != null) {
            kTScale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        int assemblysize = molecularAssemblies.length;
        switch (refinementMode) {
            case COORDINATES:
                // Compute the chemical energy.
                for (int i = 0; i < assemblysize; i++) {
                    try {
                        ForceFieldEnergy fe = molecularAssemblies[i].getPotentialEnergy();
                        getAssemblyi(i, x, xChemical[i]);
                        double curE = fe.energy(xChemical[i]);
                        e += (curE - e) / (i + 1);
                    } catch (EnergyException ex) {
                        ex.printStackTrace();
                        if (printOnFailure) {
                            String timeString = LocalDateTime.now().format(DateTimeFormatter.
                                    ofPattern("yyyy_MM_dd-HH_mm_ss"));

                            String filename = String.format("%s-ERROR-%s.pdb",
                                    FilenameUtils.removeExtension(molecularAssemblies[i].getFile().getName()),
                                    timeString);

                            PotentialsFunctions ef = new PotentialsUtils();
                            filename = ef.versionFile(filename);
                            logger.info(String.format(" Writing on-error snapshot to file %s", filename));
                            ef.saveAsPDB(molecularAssemblies[i], new File(filename));
                        }
                        if (ex.doCauseSevere()) {
                            ex.printStackTrace();
                            logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
                        } else {
                            ex.printStackTrace();
                            throw ex; // Rethrow exception
                        }

                        return 0; // Should ordinarily be unreachable.
                    }
                }
                double chemE = e;

                e = chemE * kTScale;

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
                    try {
                        ForceFieldEnergy fe = molecularAssemblies[i].getPotentialEnergy();
                        getAssemblyi(i, x, xChemical[i]);
                        double curE = fe.energy(xChemical[i]);
                        e += (curE - e) / (i + 1);
                    } catch (EnergyException ex) {
                        ex.printStackTrace();
                        if (printOnFailure) {
                            String timeString = LocalDateTime.now().format(DateTimeFormatter.
                                    ofPattern("yyyy_MM_dd-HH_mm_ss"));

                            String filename = String.format("%s-ERROR-%s.pdb",
                                    FilenameUtils.removeExtension(molecularAssemblies[i].getFile().getName()),
                                    timeString);

                            PotentialsFunctions ef = new PotentialsUtils();
                            filename = ef.versionFile(filename);
                            logger.info(String.format(" Writing on-error snapshot to file %s", filename));
                            ef.saveAsPDB(molecularAssemblies[i], new File(filename));
                        }
                        if (ex.doCauseSevere()) {
                            ex.printStackTrace();
                            logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
                        } else {
                            ex.printStackTrace();
                            throw ex; // Rethrow exception
                        }

                        return 0; // Should ordinarily be unreachable.
                    }
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
            kTScale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }

        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        int assemblysize = molecularAssemblies.length;
        switch (refinementMode) {
            case COORDINATES:
                // Compute the chemical energy and gradient.
                for (int i = 0; i < assemblysize; i++) {
                    try {
                        ForceFieldEnergy fe = molecularAssemblies[i].getPotentialEnergy();
                        getAssemblyi(i, x, xChemical[i]);
                        double curE = fe.energyAndGradient(xChemical[i], gChemical[i]);
                        e += (curE - e) / (i + 1);
                        setAssemblyi(i, g, gChemical[i]);
                    } catch (EnergyException ex) {
                        ex.printStackTrace();
                        if (printOnFailure) {
                            String timeString = LocalDateTime.now().format(DateTimeFormatter.
                                    ofPattern("yyyy_MM_dd-HH_mm_ss"));

                            String filename = String.format("%s-ERROR-%s.pdb",
                                    FilenameUtils.removeExtension(molecularAssemblies[i].getFile().getName()),
                                    timeString);

                            PotentialsFunctions ef = new PotentialsUtils();
                            filename = ef.versionFile(filename);
                            logger.info(String.format(" Writing on-error snapshot to file %s", filename));
                            ef.saveAsPDB(molecularAssemblies[i], new File(filename));
                        }
                        if (ex.doCauseSevere()) {
                            ex.printStackTrace();
                            logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
                        } else {
                            ex.printStackTrace();
                            throw ex; // Rethrow exception
                        }

                        return 0; // Should ordinarily be unreachable.
                    }
                }
                double chemE = e;

                e = chemE * kTScale;
                // normalize gradients for multiple-counted atoms
                if (assemblysize > 1) {
                    for (int i = 0; i < nXYZ; i++) {
                        g[i] /= assemblysize;
                    }
                }
                for (int i = 0; i < nXYZ; i++) {
                    g[i] *= kTScale;
                }

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != nXYZ) {
                    gXray = new double[nXYZ];
                }
                double xE = dataEnergy.energyAndGradient(x, gXray);
                //System.out.println("Xray E: " + xE + " scaled Xray E: " + weight * xE);
                e += weight * xE;

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nXYZ; i++) {
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
                    try {
                        ForceFieldEnergy fe = molecularAssemblies[i].getPotentialEnergy();
                        getAssemblyi(i, x, xChemical[i]);
                        double curE = fe.energyAndGradient(xChemical[i], gChemical[i]);
                        e += (curE - e) / (i + 1);
                        setAssemblyi(i, g, gChemical[i]);
                    } catch (EnergyException ex) {
                        ex.printStackTrace();
                        if (printOnFailure) {
                            String timeString = LocalDateTime.now().format(DateTimeFormatter.
                                    ofPattern("yyyy_MM_dd-HH_mm_ss"));

                            String filename = String.format("%s-ERROR-%s.pdb",
                                    FilenameUtils.removeExtension(molecularAssemblies[i].getFile().getName()),
                                    timeString);

                            PotentialsFunctions ef = new PotentialsUtils();
                            filename = ef.versionFile(filename);
                            logger.info(String.format(" Writing on-error snapshot to file %s", filename));
                            ef.saveAsPDB(molecularAssemblies[i], new File(filename));
                        }
                        if (ex.doCauseSevere()) {
                            ex.printStackTrace();
                            logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
                        } else {
                            ex.printStackTrace();
                            throw ex; // Rethrow exception
                        }

                        return 0; // Should ordinarily be unreachable.
                    }
                }
                // normalize gradients for multiple-counted atoms
                if (assemblysize > 1) {
                    for (int i = 0; i < nXYZ; i++) {
                        g[i] /= assemblysize;
                    }
                }

                // Compute the X-ray target energy and gradient.
                if (gXray == null || gXray.length != n) {
                    gXray = new double[n];
                }
                e += weight * dataEnergy.energyAndGradient(x, gXray);

                // Add the chemical and X-ray gradients.
                for (int i = 0; i < nXYZ; i++) {
                    g[i] += weight * gXray[i];
                }

                // bfactors, occ
                if (n > nXYZ) {
                    for (int i = nXYZ; i < n; i++) {
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
     * Get the MolecularAssembly associated with index i of n; put in xChem.
     *
     * @param i The desired MolecularAssembly index for xChem.
     * @param x All parameters.
     * @param xChem The xChem parameters for the particular MolecularAssembly
     * that will be passed to {@link ForceFieldEnergy}.
     */
    public void getAssemblyi(int i, double x[], double xChem[]) {
        assert (x != null && xChem != null);
        for (int j = 0; j < xChem.length; j += 3) {
            int index = j / 3;
            int aindex = xIndex[i].get(index) * 3;
            xChem[j] = x[aindex];
            xChem[j + 1] = x[aindex + 1];
            xChem[j + 2] = x[aindex + 2];
        }
    }

    /**
     * Set the MolecularAssembly associated with index i of n; put in x.
     *
     * @param i the desired MolecularAssembly index for "setting" x.
     * @param x All parameters.
     * @param xChem The xChem parameters for the particular MolecularAssembly
     * that will be passed to {@link ForceFieldEnergy}.
     */
    public void setAssemblyi(int i, double x[], double xChem[]) {
        assert (x != null && xChem != null);
        for (int j = 0; j < xChem.length; j += 3) {
            int index = j / 3;
            int aindex = xIndex[i].get(index) * 3;
            x[aindex] += xChem[j];
            x[aindex + 1] += xChem[j + 1];
            x[aindex + 2] += xChem[j + 2];
        }
    }

    /**
     * Get the current data weight (wA).
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
        return this.kTScale;
    }

    /**
     * set the current kT scaling weight
     *
     * @param ktscale requested kT scale
     */
    public void setKTScale(double ktscale) {
        this.kTScale = ktscale;
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
            kTScale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }
        logger.info(" kTscale: " + kTScale);
        logger.info(data.printEnergyUpdate());
        return true;
    }

    /**
     * <p>
     * Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link ffx.algorithms.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

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
        for (MolecularAssembly molecularAssembly : molecularAssemblies) {
            ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
            forceFieldEnergy.setLambda(lambda);
        }
        if (data instanceof DiffractionData) {
            XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
            xRayEnergy.setLambda(lambda);
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
            realSpaceEnergy.setLambda(lambda);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getLambda() {
        double lambda = 1.0;
        if (data instanceof DiffractionData) {
            XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
            lambda = xRayEnergy.getLambda();
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
            lambda = realSpaceEnergy.getLambda();
        }
        return lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        double dEdL = 0.0;
        if (thermostat != null) {
            kTScale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }
        int assemblysize = molecularAssemblies.length;
        /**
         * Compute the chemical energy and gradient.
         */
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy forceFieldEnergy = molecularAssemblies[i].getPotentialEnergy();
            double curdEdL = forceFieldEnergy.getdEdL();
            dEdL += (curdEdL - dEdL) / (i + 1);
        }
        dEdL *= kTScale;
        double weight = data.getWeight();
        if (data instanceof DiffractionData) {
            XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
            dEdL += weight * xRayEnergy.getdEdL();
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
            dEdL += weight * realSpaceEnergy.getdEdL();
        }
        return dEdL;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        double d2EdL2 = 0.0;
        if (thermostat != null) {
            kTScale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }
        int assemblysize = molecularAssemblies.length;
        /**
         * Compute the chemical energy and gradient.
         */
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy forceFieldEnergy = molecularAssemblies[i].getPotentialEnergy();
            double curE = forceFieldEnergy.getd2EdL2();
            d2EdL2 += (curE - d2EdL2) / (i + 1);
        }
        d2EdL2 *= kTScale;

        /**
         * No 2nd derivative for scattering term.
         */
        return d2EdL2;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        double weight = data.getWeight();
        if (thermostat != null) {
            kTScale = Thermostat.convert / (thermostat.getTargetTemperature() * Thermostat.kB);
        }
        int assemblysize = molecularAssemblies.length;

        /**
         * Compute the chemical energy and gradient.
         */
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy forcefieldEnergy = molecularAssemblies[i].getPotentialEnergy();
            Arrays.fill(gChemical[i], 0.0);
            forcefieldEnergy.getdEdXdL(gChemical[i]);
        }
        for (int i = 0; i < assemblysize; i++) {
            for (int j = 0; j < nXYZ; j++) {
                gradient[j] += gChemical[i][j];
            }
        }

        /**
         * Normalize gradients for multiple-counted atoms.
         */
        if (assemblysize > 1) {
            for (int i = 0; i < nXYZ; i++) {
                gradient[i] /= assemblysize;
            }
        }
        for (int i = 0; i < nXYZ; i++) {
            gradient[i] *= kTScale;
        }

        /**
         * Compute the X-ray target energy and gradient.
         */
        if (gXray == null || gXray.length != nXYZ) {
            gXray = new double[nXYZ];
        } else {
            for (int j = 0; j < nXYZ; j++) {
                gXray[j] = 0.0;
            }
        }
        if (data instanceof DiffractionData) {
            XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
            xRayEnergy.getdEdXdL(gXray);
        } else if (data instanceof RealSpaceData) {
            RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
            realSpaceEnergy.getdEdXdL(gXray);
        }

        // Add the chemical and X-ray gradients.
        for (int i = 0; i < nXYZ; i++) {
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
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        int assemblysize = molecularAssemblies.length;
        for (int i = 0; i < assemblysize; i++) {
            ForceFieldEnergy fe = molecularAssemblies[i].getPotentialEnergy();
            fe.setEnergyTermState(state);
        }
        dataEnergy.setEnergyTermState(state);
    }

    @Override
    public void setVelocity(double[] velocity) {
        if (this.nBFactor > 0 || this.nOccupancy > 0) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        if (velocity == null) {
            return;
        }
        int index = 0;
        double vel[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atomArray[i].isActive()) {
                vel[0] = velocity[index++];
                vel[1] = velocity[index++];
                vel[2] = velocity[index++];
                atomArray[i].setVelocity(vel);
            }
        }
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        if (this.nBFactor > 0 || this.nOccupancy > 0) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        if (acceleration == null) {
            return;
        }
        int index = 0;
        double accel[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atomArray[i].isActive()) {
                accel[0] = acceleration[index++];
                accel[1] = acceleration[index++];
                accel[2] = acceleration[index++];
                atomArray[i].setAcceleration(accel);
            }
        }
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        if (this.nBFactor > 0 || this.nOccupancy > 0) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        if (previousAcceleration == null) {
            return;
        }
        int index = 0;
        double prev[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atomArray[i].isActive()) {
                prev[0] = previousAcceleration[index++];
                prev[1] = previousAcceleration[index++];
                prev[2] = previousAcceleration[index++];
                atomArray[i].setPreviousAcceleration(prev);
            }
        }
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        if (this.nBFactor > 0 || this.nOccupancy > 0) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        int n = getNumberOfVariables();
        if (velocity == null || velocity.length < n) {
            velocity = new double[n];
        }
        int index = 0;
        double v[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atomArray[i];
            if (a.isActive()) {
                a.getVelocity(v);
                velocity[index++] = v[0];
                velocity[index++] = v[1];
                velocity[index++] = v[2];
            }
        }
        return velocity;
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        if (this.nBFactor > 0 || this.nOccupancy > 0) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        int n = getNumberOfVariables();
        if (acceleration == null || acceleration.length < n) {
            acceleration = new double[n];
        }
        int index = 0;
        double a[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atomArray[i].isActive()) {
                atomArray[i].getAcceleration(a);
                acceleration[index++] = a[0];
                acceleration[index++] = a[1];
                acceleration[index++] = a[2];
            }
        }
        return acceleration;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        if (this.nBFactor > 0 || this.nOccupancy > 0) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        int n = getNumberOfVariables();
        if (previousAcceleration == null || previousAcceleration.length < n) {
            previousAcceleration = new double[n];
        }
        int index = 0;
        double a[] = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atomArray[i].isActive()) {
                atomArray[i].getPreviousAcceleration(a);
                previousAcceleration[index++] = a[0];
                previousAcceleration[index++] = a[1];
                previousAcceleration[index++] = a[2];
            }
        }
        return previousAcceleration;
    }
}
