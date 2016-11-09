/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;

import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.algorithms.Thermostat.convert;
import static ffx.algorithms.Thermostat.kB;
import static ffx.numerics.VectorMath.b2u;
import static ffx.numerics.VectorMath.determinant3;
import static ffx.numerics.VectorMath.u2b;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class XRayEnergy implements LambdaInterface, Potential {

    private static final Logger logger = Logger.getLogger(XRayEnergy.class.getName());
    private static final double kBkcal = kB / convert;
    private static final double eightPI2 = 8.0 * PI * PI;
    private static final double eightPI23 = eightPI2 * eightPI2 * eightPI2;

    private final DiffractionData diffractionData;
    private final RefinementModel refinementModel;
    private final Atom atomArray[];
    private final Atom activeAtomArray[];
    private final int nAtoms;
    private final int nActiveAtoms;
    private int nXYZ;
    private int nB;
    private int nOCC;

    private RefinementMode refinementMode;
    private boolean refineXYZ = false;
    private boolean refineOCC = false;
    private boolean refineB = false;
    private boolean xrayTerms = true;
    private boolean restraintTerms = true;
    protected double[] optimizationScaling = null;

    private final double bMass;
    private final double kTbNonzero;
    private final double kTbSimWeight;
    private final double occMass;
    private final double temperature = 50.0;

    protected double lambda = 1.0;
    private boolean lambdaTerm = false;
    private double totalEnergy;
    private double dEdL;
    private double d2UdL2 = 0.0;
    private double g2[];
    private double dUdXdL[];
    private STATE state = STATE.BOTH;

    /**
     * Diffraction data energy target
     *
     * @param diffractionData {@link DiffractionData} object to associate with
     * the target
     * @param nXYZ number of xyz parameters
     * @param nB number of b factor parameters
     * @param nOCC number of occupancy parameters
     * @param refinementMode the {@link RefinementMinimize.RefinementMode} type
     * of refinement requested
     */
    public XRayEnergy(DiffractionData diffractionData, int nXYZ, int nB, int nOCC,
            RefinementMode refinementMode) {
        this.diffractionData = diffractionData;
        this.refinementModel = diffractionData.getRefinementModel();
        this.refinementMode = refinementMode;
        this.atomArray = refinementModel.getTotalAtomArray();
        this.nAtoms = atomArray.length;
        this.nXYZ = nXYZ;
        this.nB = nB;
        this.nOCC = nOCC;

        bMass = diffractionData.getbMass();
        kTbNonzero = 0.5 * kBkcal * temperature * diffractionData.getbNonZeroWeight();
        kTbSimWeight = kBkcal * temperature * diffractionData.getbSimWeight();
        occMass = diffractionData.getOccMass();

        ForceField forceField = diffractionData.getAssembly()[0].getForceField();
        lambdaTerm = forceField.getBoolean(ForceFieldBoolean.LAMBDATERM, false);

        // Fill an active atom array.
        int count = 0;
        for (Atom a : atomArray) {
            if (a.isActive()) {
                count++;
            }
        }
        nActiveAtoms = count;
        activeAtomArray = new Atom[count];
        count = 0;
        for (Atom a : atomArray) {
            if (a.isActive()) {
                activeAtomArray[count++] = a;
            }
        }

        dUdXdL = new double[count * 3];
        g2 = new double[count * 3];

        setRefinementBooleans();

        if (refineB) {
            logger.info(" B-Factor Refinement Parameters");
            logger.info("  Temperature:                 " + temperature);
            logger.info("  Non-zero restraint weight:   " + diffractionData.getbNonZeroWeight());
            logger.info("  Similarity restraint weight: " + diffractionData.getbSimWeight());
        }

        logger.info(String.format(" XRayEnergy variables:  %d (nXYZ %d, nB %d, nOcc %d)\n",
                nXYZ + nB + nOCC, nXYZ, nB, nOCC));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {
        double e = 0.0;
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        if (refineXYZ) {
            // update coordinates
            diffractionData.setFFTCoordinates(x);
        }
        if (refineB) {
            // update B factors
            setBFactors(x);
        }
        if (refineOCC) {
            // update occupancies
            setOccupancies(x);
        }

        if (xrayTerms) {

            if (lambdaTerm) {
                diffractionData.setLambdaTerm(false);
            }

            // Compute new structure factors.
            diffractionData.computeAtomicDensity();
            // Compute crystal likelihood.
            e = diffractionData.computeLikelihood();

            if (lambdaTerm) {

                // Turn off all atoms scaled by lambda.
                diffractionData.setLambdaTerm(true);

                // Compute new structure factors.
                diffractionData.computeAtomicDensity();

                // Compute crystal likelihood.
                double e2 = diffractionData.computeLikelihood();

                dEdL = e - e2;

                e = lambda * e + (1.0 - lambda) * e2;

                diffractionData.setLambdaTerm(false);
            }

        }

        if (restraintTerms) {
            if (refineB) {
                // add B restraints
                e += getBFactorRestraints();
            }
        }

        /**
         * Scale the coordinates and gradients.
         */
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
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double e = 0.0;
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        if (refineXYZ) {
            for (Atom a : activeAtomArray) {
                a.setXYZGradient(0.0, 0.0, 0.0);
                a.setLambdaXYZGradient(0.0, 0.0, 0.0);
            }

            // update coordinates
            diffractionData.setFFTCoordinates(x);
        }

        if (refineB) {
            for (Atom a : activeAtomArray) {
                a.setTempFactorGradient(0.0);
                if (a.getAnisou(null) != null) {
                    if (a.getAnisouGradient(null) == null) {
                        double ganisou[] = new double[6];
                        a.setAnisouGradient(ganisou);
                    } else {
                        double ganisou[] = a.getAnisouGradient(null);
                        ganisou[0] = ganisou[1] = ganisou[2] = 0.0;
                        ganisou[3] = ganisou[4] = ganisou[5] = 0.0;
                        a.setAnisouGradient(ganisou);
                    }
                }
            }

            // update B factors
            setBFactors(x);
        }

        if (refineOCC) {
            for (Atom a : activeAtomArray) {
                a.setOccupancyGradient(0.0);
            }

            // update occupancies
            setOccupancies(x);
        }

        if (xrayTerms) {

            if (lambdaTerm) {
                diffractionData.setLambdaTerm(false);
            }

            // compute new structure factors
            diffractionData.computeAtomicDensity();

            // compute crystal likelihood
            e = diffractionData.computeLikelihood();

            // compute the crystal gradients
            diffractionData.computeAtomicGradients(refinementMode);

            if (refineXYZ) {
                // pack gradients into gradient array
                getXYZGradients(g);
            }

            if (lambdaTerm) {

                int n = dUdXdL.length;
                System.arraycopy(g, 0, dUdXdL, 0, n);

                for (Atom a : activeAtomArray) {
                    a.setXYZGradient(0.0, 0.0, 0.0);
                    a.setLambdaXYZGradient(0.0, 0.0, 0.0);
                }

                // Turn off all atoms scaled by lambda.
                diffractionData.setLambdaTerm(true);

                // Compute new structure factors.
                diffractionData.computeAtomicDensity();

                // Compute crystal likelihood.
                double e2 = diffractionData.computeLikelihood();

                // compute the crystal gradients
                diffractionData.computeAtomicGradients(refinementMode);

                dEdL = e - e2;
                e = lambda * e + (1.0 - lambda) * e2;
                getXYZGradients(g2);
                for (int i = 0; i < g.length; i++) {
                    dUdXdL[i] -= g2[i];
                    g[i] = lambda * g[i] + (1.0 - lambda) * g2[i];
                }

                diffractionData.setLambdaTerm(false);
            }

        }

        if (restraintTerms) {
            if (refineB) {
                // add B restraints
                e += getBFactorRestraints();

                // pack gradients into gradient array
                getBFactorGradients(g);
            }

            if (refineOCC) {
                // pack gradients into gradient array
                getOccupancyGradients(g);
            }
        }

        /**
         * Scale the coordinates and gradients.
         */
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
     * <p>
     * Getter for the field <code>refinementMode</code>.</p>
     *
     * @return a {@link ffx.xray.RefinementMinimize.RefinementMode} object.
     */
    public RefinementMode getRefinementMode() {
        return refinementMode;
    }

    /**
     * <p>
     * Setter for the field <code>refinementMode</code>.</p>
     *
     * @param refinementmode a
     * {@link ffx.xray.RefinementMinimize.RefinementMode} object.
     */
    public void setRefinementMode(RefinementMode refinementmode) {
        this.refinementMode = refinementmode;
        setRefinementBooleans();
    }

    /**
     * if the refinement mode has changed, this should be called to update which
     * parameters are being fit
     */
    private void setRefinementBooleans() {
        // reset, if previously set
        refineXYZ = false;
        refineB = false;
        refineOCC = false;

        if (refinementMode == RefinementMode.COORDINATES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineXYZ = true;
        }

        if (refinementMode == RefinementMode.BFACTORS
                || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineB = true;
        }

        if (refinementMode == RefinementMode.OCCUPANCIES
                || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineOCC = true;
        }
    }

    /**
     * Get the number of xyz parameters being fit.
     *
     * @return the number of xyz parameters
     */
    public int getNXYZ() {
        return nXYZ;
    }

    /**
     * set the number of xyz parameters
     *
     * @param nXYZ requested number of xyz parameters
     */
    public void setNXYZ(int nXYZ) {
        this.nXYZ = nXYZ;
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
     * set the number of B factor parameters
     *
     * @param nB requested number of B factor parameters
     */
    public void setNB(int nB) {
        this.nB = nB;
    }

    /**
     * get the number of occupancy parameters being fit
     *
     * @return the number of occupancy parameters
     */
    public int getNOcc() {
        return nOCC;
    }

    /**
     * set the number of occupancy parameters
     *
     * @param nOCC requested number of occupancy parameters
     */
    public void setNOcc(int nOCC) {
        this.nOCC = nOCC;
    }

    /**
     * fill gradient array with B factor gradients
     *
     * @param g array to add gradients to
     */
    public void getBFactorGradients(double g[]) {
        assert (g != null);
        double grad[] = null;
        int index = nXYZ;
        int resnum = -1;
        int nres = diffractionData.getnResidueBFactor() + 1;
        for (Atom a : activeAtomArray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou(null) != null) {
                grad = a.getAnisouGradient(grad);
                g[index++] = grad[0];
                g[index++] = grad[1];
                g[index++] = grad[2];
                g[index++] = grad[3];
                g[index++] = grad[4];
                g[index++] = grad[5];
            } else if (diffractionData.isResidueBFactor()) {
                if (resnum != a.getResidueNumber()) {
                    if (nres >= diffractionData.getnResidueBFactor()) {
                        if (resnum > -1
                                && index < nXYZ + nB - 1) {
                            index++;
                        }
                        nres = 1;
                    } else {
                        nres++;
                    }
                    g[index] += a.getTempFactorGradient();
                    resnum = a.getResidueNumber();
                } else {
                    g[index] += a.getTempFactorGradient();
                }
            } else {
                g[index++] = a.getTempFactorGradient();
            }
        }
    }

    /**
     * Fill gradient array with occupancy gradients.
     *
     * @param g array to add gradients to
     */
    public void getOccupancyGradients(double g[]) {
        double ave;
        int index = nXYZ + nB;

        // First: Alternate Residues
        for (ArrayList<Residue> list : refinementModel.getAltResidues()) {
            ave = 0.0;
            for (Residue r : list) {
                for (Atom a : r.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        ave += a.getOccupancyGradient();
                    }
                }
            }
            ave /= list.size();
            for (Residue r : list) {
                for (Atom a : r.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        g[index] += a.getOccupancyGradient();
                    }
                }
                if (list.size() > 1) {
                    g[index] -= ave;
                }
                index++;
            }
        }

        // Now the molecules (HETATMs).
        for (ArrayList<Molecule> list : refinementModel.getAltMolecules()) {
            ave = 0.0;
            for (Molecule m : list) {
                for (Atom a : m.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        ave += a.getOccupancyGradient();
                    }
                }
            }
            ave /= list.size();
            for (Molecule m : list) {
                for (Atom a : m.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        g[index] += a.getOccupancyGradient();
                    }
                }
                if (list.size() > 1) {
                    g[index] -= ave;
                }
                index++;
            }
        }
    }

    /**
     * Fill gradient array with atomic coordinate partial derivatives.
     *
     * @param g gradient array
     */
    public void getXYZGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : activeAtomArray) {
            a.getXYZGradient(grad);
            g[index++] = grad[0];
            g[index++] = grad[1];
            g[index++] = grad[2];
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double x[]) {
        assert (x != null);
        double xyz[] = new double[3];
        int index = 0;
        fill(x, 0.0);

        if (refineXYZ) {
            for (Atom a : activeAtomArray) {
                a.getXYZ(xyz);
                x[index++] = xyz[0];
                x[index++] = xyz[1];
                x[index++] = xyz[2];
            }
        }

        if (refineB) {
            double anisou[] = null;
            int resnum = -1;
            int nat = 0;
            int nres = diffractionData.getnResidueBFactor() + 1;
            for (Atom a : activeAtomArray) {
                // ignore hydrogens!!!
                if (a.getAtomicNumber() == 1) {
                    continue;
                }
                if (a.getAnisou(null) != null) {
                    anisou = a.getAnisou(anisou);
                    x[index++] = anisou[0];
                    x[index++] = anisou[1];
                    x[index++] = anisou[2];
                    x[index++] = anisou[3];
                    x[index++] = anisou[4];
                    x[index++] = anisou[5];
                } else if (diffractionData.isResidueBFactor()) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= diffractionData.getnResidueBFactor()) {
                            if (resnum > -1
                                    && index < nXYZ + nB - 1) {
                                x[index] /= nat;
                                index++;
                            }
                            nat = 1;
                            nres = 1;
                        } else {
                            nres++;
                            nat++;
                        }
                        x[index] += a.getTempFactor();
                        resnum = a.getResidueNumber();
                    } else {
                        x[index] += a.getTempFactor();
                        nat++;
                    }
                } else {
                    x[index++] = a.getTempFactor();
                }
            }

            if (diffractionData.isResidueBFactor()) {
                if (nat > 1) {
                    x[index] /= nat;
                }
            }
        }

        if (refineOCC) {
            for (ArrayList<Residue> list : refinementModel.getAltResidues()) {
                for (Residue r : list) {
                    for (Atom a : r.getAtomList()) {
                        if (a.getOccupancy() < 1.0) {
                            x[index++] = a.getOccupancy();
                            break;
                        }
                    }
                }
            }
            for (ArrayList<Molecule> list : refinementModel.getAltMolecules()) {
                for (Molecule m : list) {
                    for (Atom a : m.getAtomList()) {
                        if (a.getOccupancy() < 1.0) {
                            x[index++] = a.getOccupancy();
                            break;
                        }
                    }
                }
            }
        }

        return x;
    }

    /**
     * set atomic B factors based on current position
     *
     * @param x current parameters to set B factors with
     */
    public void setBFactors(double x[]) {
        double tmpanisou[] = new double[6];
        int index = nXYZ;
        int nneg = 0;
        int resnum = -1;
        int nres = diffractionData.getnResidueBFactor() + 1;
        for (Atom a : activeAtomArray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou(null) == null) {
                double biso = x[index];
                if (diffractionData.isResidueBFactor()) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= diffractionData.getnResidueBFactor()) {
                            if (resnum > -1
                                    && index < nXYZ + nB - 1) {
                                index++;
                                biso = x[index];
                            }
                            nres = 1;
                        } else {
                            nres++;
                        }
                        resnum = a.getResidueNumber();
                    }
                } else {
                    index++;
                }

                if (biso > 0.0) {
                    a.setTempFactor(biso);
                } else {
                    nneg++;
                    a.setTempFactor(0.01);
                    if (nneg < 5) {
                        logger.info(" Isotropic atom: " + a.toString() + " negative B factor");
                    }
                }
            } else {
                double anisou[] = a.getAnisou(null);
                tmpanisou[0] = x[index++];
                tmpanisou[1] = x[index++];
                tmpanisou[2] = x[index++];
                tmpanisou[3] = x[index++];
                tmpanisou[4] = x[index++];
                tmpanisou[5] = x[index++];
                double det = determinant3(tmpanisou);
                if (det > 0.0) {
                    System.arraycopy(tmpanisou, 0, anisou, 0, 6);
                    a.setAnisou(anisou);
                    det = Math.pow(det, 0.3333);
                    a.setTempFactor(u2b(det));
                } else {
                    nneg++;
                    a.setTempFactor(0.01);
                    anisou[0] = anisou[1] = anisou[2] = b2u(0.01);
                    anisou[3] = anisou[4] = anisou[5] = 0.0;
                    a.setAnisou(anisou);
                    if (nneg < 5) {
                        logger.info(" Anisotropic atom: " + a.toString() + " negative ANISOU");
                    }
                }
            }
        }

        if (nneg > 0) {
            logger.info(" " + nneg + " of " + nAtoms
                    + " atoms with negative B factors! Attempting to correct.\n  (If this problem persists, increase bsimweight)");
            /*
             * if (nneg > 50){ kTbsim *= 2.0; logger.info("excessive number of
             * negative Bs, increasing similarity restraint: " + kTbsim); }
             */
        }

        // set hydrogen based on bonded atom
        for (Atom a : activeAtomArray) {
            if (a.getAtomicNumber() == 1) {
                Atom b = a.getBonds().get(0).get1_2(a);
                a.setTempFactor(b.getTempFactor());
            }
        }
    }

    /**
     * set atomic xyz coordinates based on current position
     *
     * @param x current parameters to set coordinates with
     */
    public void setCoordinates(double x[]) {
        assert (x != null);
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : activeAtomArray) {
            xyz[0] = x[index++];
            xyz[1] = x[index++];
            xyz[2] = x[index++];
            a.moveTo(xyz);
        }
    }

    /**
     * set atom occupancies based on current position
     *
     * @param x current parameters to set occupancies with
     */
    public void setOccupancies(double x[]) {
        double occ = 0.0;
        int index = nXYZ + nB;
        for (ArrayList<Residue> list : refinementModel.getAltResidues()) {
            for (Residue r : list) {
                occ = x[index++];
                for (Atom a : r.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        a.setOccupancy(occ);
                    }
                }
            }
        }
        for (ArrayList<Molecule> list : refinementModel.getAltMolecules()) {
            for (Molecule m : list) {
                occ = x[index++];
                for (Atom a : m.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        a.setOccupancy(occ);
                    }
                }
            }
        }
    }

    /**
     * determine similarity and non-zero B factor restraints (done independently
     * of getBFactorGradients), affects atomic gradients
     *
     * @return energy of the restraint
     */
    public double getBFactorRestraints() {
        Atom a1, a2;
        double b1, b2, bdiff;
        double anisou1[] = null;
        double anisou2[] = null;
        double gradb;
        double det1, det2;
        double gradu[] = new double[6];
        double e = 0.0;

        for (Atom a : activeAtomArray) {
            double biso = a.getTempFactor();
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }

            if (a.getAnisou(null) == null) {
                // isotropic B restraint

                // non-zero restraint: -kTln[Z], Z is ADP partition function
                e += -3.0 * kTbNonzero * Math.log(biso / (4.0 * Math.PI));
                gradb = -3.0 * kTbNonzero / biso;
                a.addToTempFactorGradient(gradb);

                // similarity harmonic restraint
                ArrayList<Bond> bonds = a.getBonds();
                for (Bond b : bonds) {
                    if (a.compareTo(b.getAtom(0)) == 0) {
                        a1 = b.getAtom(0);
                        a2 = b.getAtom(1);
                    } else {
                        a1 = b.getAtom(1);
                        a2 = b.getAtom(0);
                    }
                    if (a2.getAtomicNumber() == 1) {
                        continue;
                    }
                    if (a2.xyzIndex < a1.xyzIndex) {
                        continue;
                    }
                    if (!a1.getAltLoc().equals(' ')
                            && !a1.getAltLoc().equals('A')
                            && a2.getAltLoc().equals(' ')) {
                        continue;
                    }

                    b1 = a1.getTempFactor();
                    b2 = a2.getTempFactor();
                    // transform B similarity restraints to U scale
                    bdiff = b2u(b1 - b2);
                    e += kTbSimWeight * Math.pow(bdiff, 2.0);
                    gradb = 2.0 * kTbSimWeight * bdiff;
                    a1.addToTempFactorGradient(gradb);
                    a2.addToTempFactorGradient(-gradb);
                }
            } else {
                // anisotropic B restraint
                anisou1 = a.getAnisou(anisou1);
                det1 = determinant3(anisou1);

                // non-zero restraint: -kTln[Z], Z is ADP partition function
                e += u2b(-kTbNonzero * Math.log(det1 * eightPI2 * Math.PI));
                gradu[0] = u2b(-kTbNonzero * ((anisou1[1] * anisou1[2] - anisou1[5] * anisou1[5]) / det1));
                gradu[1] = u2b(-kTbNonzero * ((anisou1[0] * anisou1[2] - anisou1[4] * anisou1[4]) / det1));
                gradu[2] = u2b(-kTbNonzero * ((anisou1[0] * anisou1[1] - anisou1[3] * anisou1[3]) / det1));
                gradu[3] = u2b(-kTbNonzero * ((2.0 * (anisou1[4] * anisou1[5] - anisou1[3] * anisou1[2])) / det1));
                gradu[4] = u2b(-kTbNonzero * ((2.0 * (anisou1[3] * anisou1[5] - anisou1[4] * anisou1[1])) / det1));
                gradu[5] = u2b(-kTbNonzero * ((2.0 * (anisou1[3] * anisou1[4] - anisou1[5] * anisou1[0])) / det1));
                a.addToAnisouGradient(gradu);

                // similarity harmonic restraint based on determinants
                ArrayList<Bond> bonds = a.getBonds();
                for (Bond b : bonds) {
                    if (a.compareTo(b.getAtom(0)) == 0) {
                        a1 = b.getAtom(0);
                        a2 = b.getAtom(1);
                    } else {
                        a1 = b.getAtom(1);
                        a2 = b.getAtom(0);
                    }
                    if (a2.getAtomicNumber() == 1) {
                        continue;
                    }
                    if (a2.xyzIndex < a1.xyzIndex) {
                        continue;
                    }
                    if (a2.getAnisou(null) == null) {
                        continue;
                    }
                    if (!a1.getAltLoc().equals(' ')
                            && !a1.getAltLoc().equals('A')
                            && a2.getAltLoc().equals(' ')) {
                        continue;
                    }

                    anisou2 = a2.getAnisou(anisou2);
                    det2 = determinant3(anisou2);
                    bdiff = det1 - det2;
                    e += eightPI23 * kTbSimWeight * Math.pow(bdiff, 2.0);
                    gradb = eightPI23 * 2.0 * kTbSimWeight * bdiff;

                    // parent atom
                    gradu[0] = gradb * (anisou1[1] * anisou1[2] - anisou1[5] * anisou1[5]);
                    gradu[1] = gradb * (anisou1[0] * anisou1[2] - anisou1[4] * anisou1[4]);
                    gradu[2] = gradb * (anisou1[0] * anisou1[1] - anisou1[3] * anisou1[3]);
                    gradu[3] = gradb * (2.0 * (anisou1[4] * anisou1[5] - anisou1[3] * anisou1[2]));
                    gradu[4] = gradb * (2.0 * (anisou1[3] * anisou1[5] - anisou1[4] * anisou1[1]));
                    gradu[5] = gradb * (2.0 * (anisou1[3] * anisou1[4] - anisou1[5] * anisou1[0]));
                    a1.addToAnisouGradient(gradu);

                    // bonded atom
                    gradu[0] = gradb * (anisou2[5] * anisou2[5] - anisou2[1] * anisou2[2]);
                    gradu[1] = gradb * (anisou2[4] * anisou2[4] - anisou2[0] * anisou2[2]);
                    gradu[2] = gradb * (anisou2[3] * anisou2[3] - anisou2[0] * anisou2[1]);
                    gradu[3] = gradb * (2.0 * (anisou2[3] * anisou2[2] - anisou2[4] * anisou2[5]));
                    gradu[4] = gradb * (2.0 * (anisou2[4] * anisou2[1] - anisou2[3] * anisou2[5]));
                    gradu[5] = gradb * (2.0 * (anisou2[5] * anisou2[0] - anisou2[3] * anisou2[4]));
                    a2.addToAnisouGradient(gradu);
                }
            }
        }
        return e;
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
        double mass[] = new double[nXYZ + nB + nOCC];
        int i = 0;
        if (refineXYZ) {
            for (Atom a : activeAtomArray) {
                double m = a.getMass();
                mass[i++] = m;
                mass[i++] = m;
                mass[i++] = m;
            }
        }

        if (refineB) {
            for (int j = i; j < nXYZ + nB; i++, j++) {
                mass[j] = bMass;
            }
        }

        if (refineOCC) {
            for (int j = i; j < nXYZ + nB + nOCC; i++, j++) {
                mass[j] = occMass;
            }
        }
        return mass;
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
        return nXYZ + nB + nOCC;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
        } else {
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getLambda() {
        return lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        return dEdL;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        int n = dUdXdL.length;
        System.arraycopy(dUdXdL, 0, gradient, 0, n);
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        VARIABLE_TYPE vtypes[] = new VARIABLE_TYPE[nXYZ + nB + nOCC];
        int i = 0;
        if (refineXYZ) {
            for (Atom a : activeAtomArray) {
                vtypes[i++] = VARIABLE_TYPE.X;
                vtypes[i++] = VARIABLE_TYPE.Y;
                vtypes[i++] = VARIABLE_TYPE.Z;
            }
        }

        if (refineB) {
            for (int j = i; j < nXYZ + nB; i++, j++) {
                vtypes[j] = VARIABLE_TYPE.OTHER;
            }
        }

        if (refineOCC) {
            for (int j = i; j < nXYZ + nB + nOCC; i++, j++) {
                vtypes[j] = VARIABLE_TYPE.OTHER;
            }
        }
        return vtypes;
    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    /*
     * RESPA setup
     */
    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        switch (state) {
            case FAST:
                xrayTerms = false;
                restraintTerms = true;
                break;
            case SLOW:
                xrayTerms = true;
                restraintTerms = false;
                break;
            default:
                xrayTerms = true;
                restraintTerms = true;
        }
    }

    @Override
    public void setVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
