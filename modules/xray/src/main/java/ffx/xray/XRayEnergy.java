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
import java.util.Arrays;
import java.util.logging.Logger;

import ffx.numerics.Potential;
import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.algorithms.Thermostat.convert;
import static ffx.algorithms.Thermostat.kB;
import static ffx.numerics.VectorMath.*;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 * @since 1.0
 *
 */
public class XRayEnergy implements LambdaInterface, Potential {

    private static final Logger logger = Logger.getLogger(XRayEnergy.class.getName());
    private static final double kBkcal = kB / convert;
    private static final double eightpi2 = 8.0 * Math.PI * Math.PI;
    private static final double eightpi23 = eightpi2 * eightpi2 * eightpi2;
    private final DiffractionData diffractiondata;
    private final RefinementModel refinementmodel;
    private final Atom atomarray[];
    private final int nAtoms;
    private int nxyz;
    private int nb;
    private int nocc;
    private RefinementMode refinementMode;
    private boolean refinexyz = false;
    private boolean refineocc = false;
    private boolean refineb = false;
    private boolean xrayterms = true;
    private boolean restraintterms = true;
    protected double[] optimizationScaling = null;
    private double bmass;
    private double kTbnonzero;
    private double kTbsim;
    private double occmass;
    private double temp = 50.0;
    protected double lambda = 1.0;
    private double totalEnergy;

    /**
     * Diffraction data energy target
     *
     * @param diffractiondata {@link DiffractionData} object to associate with
     * the target
     * @param nxyz number of xyz parameters
     * @param nb number of b factor parameters
     * @param nocc number of occupancy parameters
     * @param refinementmode the {@link RefinementMinimize.RefinementMode} type
     * of refinement requested
     */
    public XRayEnergy(DiffractionData diffractiondata, int nxyz, int nb, int nocc,
            RefinementMode refinementmode) {
        this.diffractiondata = diffractiondata;
        this.refinementmodel = diffractiondata.getRefinementModel();
        this.refinementMode = refinementmode;
        this.atomarray = refinementmodel.atomarray;
        this.nAtoms = atomarray.length;
        this.nxyz = nxyz;
        this.nb = nb;
        this.nocc = nocc;

        bmass = diffractiondata.bmass;
        kTbnonzero = 0.5 * kBkcal * temp * diffractiondata.bnonzeroweight;
        kTbsim = kBkcal * temp * diffractiondata.bsimweight;
        occmass = diffractiondata.occmass;

        setRefinementBooleans();

        if (refineb) {
            logger.info(" B-Factor Refinement Parameters");
            logger.info(" Temperature:                 " + temp);
            logger.info(" Non-zero restraint weight:   " + diffractiondata.bnonzeroweight);
            logger.info(" Similarity restraint weight: " + diffractiondata.bsimweight);
        }
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

        if (refinexyz) {
            // update coordinates
            diffractiondata.setFFTCoordinates(x);
        }
        if (refineb) {
            // update B factors
            setBFactors(x);
        }
        if (refineocc) {
            // update occupancies
            setOccupancies(x);
        }

        if (xrayterms) {
            // compute new structure factors
            diffractiondata.computeAtomicDensity();

            // compute crystal likelihood
            e = diffractiondata.computeLikelihood();
        }

        if (restraintterms) {
            if (refineb) {
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

        if (refinexyz) {
            for (Atom a : atomarray) {
                a.setXYZGradient(0.0, 0.0, 0.0);
            }

            // update coordinates
            diffractiondata.setFFTCoordinates(x);
        }

        if (refineb) {
            for (Atom a : atomarray) {
                a.setTempFactorGradient(0.0);
                if (a.getAnisou() != null) {
                    if (a.getAnisouGradient() == null) {
                        double ganisou[] = new double[6];
                        a.setAnisouGradient(ganisou);
                    } else {
                        double ganisou[] = a.getAnisouGradient();
                        ganisou[0] = ganisou[1] = ganisou[2] = 0.0;
                        ganisou[3] = ganisou[4] = ganisou[5] = 0.0;
                    }
                }
            }

            // update B factors
            setBFactors(x);
        }

        if (refineocc) {
            for (Atom a : atomarray) {
                a.setOccupancyGradient(0.0);
            }

            // update occupancies
            setOccupancies(x);
        }

        if (xrayterms) {
            // compute new structure factors
            diffractiondata.computeAtomicDensity();

            // compute crystal likelihood
            e = diffractiondata.computeLikelihood();

            // compute the crystal gradients
            diffractiondata.computeAtomicGradients(refinementMode);

            if (refinexyz) {
                // pack gradients into gradient array
                getXYZGradients(g);
            }
        }

        if (restraintterms) {
            if (refineb) {
                // add B restraints
                e += getBFactorRestraints();

                // pack gradients into gradient array
                getBFactorGradients(g);
            }

            if (refineocc) {
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
     * <p>Getter for the field
     * <code>refinementMode</code>.</p>
     *
     * @return a {@link ffx.xray.RefinementMinimize.RefinementMode} object.
     */
    public RefinementMode getRefinementMode() {
        return refinementMode;
    }

    /**
     * <p>Setter for the field
     * <code>refinementMode</code>.</p>
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
        refinexyz = false;
        refineb = false;
        refineocc = false;

        if (refinementMode == RefinementMode.COORDINATES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refinexyz = true;
        }

        if (refinementMode == RefinementMode.BFACTORS
                || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineb = true;
        }

        if (refinementMode == RefinementMode.OCCUPANCIES
                || refinementMode == RefinementMode.BFACTORS_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineocc = true;
        }
    }

    /**
     * get the number of xyz parameters being fit
     *
     * @return the number of xyz parameters
     */
    public int getNXYZ() {
        return nxyz;
    }

    /**
     * set the number of xyz parameters
     *
     * @param nxyz requested number of xyz parameters
     */
    public void setNXYZ(int nxyz) {
        this.nxyz = nxyz;
    }

    /**
     * get the number of B factor parameters being fit
     *
     * @return the number of B factor parameters
     */
    public int getNB() {
        return nb;
    }

    /**
     * set the number of B factor parameters
     *
     * @param nb requested number of B factor parameters
     */
    public void setNB(int nb) {
        this.nb = nb;
    }

    /**
     * get the number of occupancy parameters being fit
     *
     * @return the number of occupancy parameters
     */
    public int getNOcc() {
        return nocc;
    }

    /**
     * set the number of occupancy parameters
     *
     * @param nocc requested number of occupancy parameters
     */
    public void setNOcc(int nocc) {
        this.nocc = nocc;
    }

    /**
     * fill gradient array with B factor gradients
     *
     * @param g array to add gradients to
     */
    public void getBFactorGradients(double g[]) {
        assert (g != null);
        double grad[];
        int index = nxyz;
        int resnum = -1;
        int nres = diffractiondata.nresiduebfactor + 1;
        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() != null) {
                grad = a.getAnisouGradient();
                g[index++] = grad[0];
                g[index++] = grad[1];
                g[index++] = grad[2];
                g[index++] = grad[3];
                g[index++] = grad[4];
                g[index++] = grad[5];
            } else if (diffractiondata.residuebfactor) {
                if (resnum != a.getResidueNumber()) {
                    if (nres >= diffractiondata.nresiduebfactor) {
                        if (resnum > -1
                                && index < nxyz + nb - 1) {
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
     * fill gradient array with occupancy gradients
     *
     * @param g array to add gradients to
     */
    public void getOccupancyGradients(double g[]) {
        double ave;
        int index = nxyz + nb;

        // first: alternate residues
        for (ArrayList<Residue> list : refinementmodel.altresidues) {
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

        // now the molecules (HETATMs)
        for (ArrayList<Molecule> list : refinementmodel.altmolecules) {
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
     * fill gradient array with xyz gradients
     *
     * @param g array to add gradients to
     */
    public void getXYZGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : atomarray) {
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
        Arrays.fill(x, 0.0);

        if (refinexyz) {
            for (Atom a : atomarray) {
                a.getXYZ(xyz);
                x[index++] = xyz[0];
                x[index++] = xyz[1];
                x[index++] = xyz[2];
            }
        }

        if (refineb) {
            double anisou[];
            int resnum = -1;
            int nat = 0;
            int nres = diffractiondata.nresiduebfactor + 1;
            for (Atom a : atomarray) {
                // ignore hydrogens!!!
                if (a.getAtomicNumber() == 1) {
                    continue;
                }
                if (a.getAnisou() != null) {
                    anisou = a.getAnisou();
                    x[index++] = anisou[0];
                    x[index++] = anisou[1];
                    x[index++] = anisou[2];
                    x[index++] = anisou[3];
                    x[index++] = anisou[4];
                    x[index++] = anisou[5];
                } else if (diffractiondata.residuebfactor) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= diffractiondata.nresiduebfactor) {
                            if (resnum > -1
                                    && index < nxyz + nb - 1) {
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

            if (diffractiondata.residuebfactor) {
                if (nat > 1) {
                    x[index] /= nat;
                }
            }
        }

        if (refineocc) {
            for (ArrayList<Residue> list : refinementmodel.altresidues) {
                for (Residue r : list) {
                    for (Atom a : r.getAtomList()) {
                        if (a.getOccupancy() < 1.0) {
                            x[index++] = a.getOccupancy();
                            break;
                        }
                    }
                }
            }
            for (ArrayList<Molecule> list : refinementmodel.altmolecules) {
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
        int index = nxyz;
        int nneg = 0;
        int resnum = -1;
        int nres = diffractiondata.nresiduebfactor + 1;
        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() == null) {
                double biso = x[index];
                if (diffractiondata.residuebfactor) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= diffractiondata.nresiduebfactor) {
                            if (resnum > -1
                                    && index < nxyz + nb - 1) {
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
                        logger.info("isotropic atom: " + a.toString() + " negative B factor");
                    }
                }
            } else {
                double anisou[] = a.getAnisou();
                tmpanisou[0] = x[index++];
                tmpanisou[1] = x[index++];
                tmpanisou[2] = x[index++];
                tmpanisou[3] = x[index++];
                tmpanisou[4] = x[index++];
                tmpanisou[5] = x[index++];
                double det = determinant3(tmpanisou);
                if (det > 0.0) {
                    System.arraycopy(tmpanisou, 0, anisou, 0, 6);
                    det = Math.pow(det, 0.3333);
                    a.setTempFactor(u2b(det));
                } else {
                    nneg++;
                    a.setTempFactor(0.01);
                    anisou[0] = anisou[1] = anisou[2] = b2u(0.01);
                    anisou[3] = anisou[4] = anisou[5] = 0.0;
                    if (nneg < 5) {
                        logger.info("anisotropic atom: " + a.toString() + " negative ANISOU");
                    }
                }
            }
        }

        if (nneg > 0) {
            logger.info(nneg + " of " + nAtoms
                    + " atoms with negative B factors! Attempting to correct.\n  (If this problem persists, increase bsimweight)");
            /*
             * if (nneg > 50){ kTbsim *= 2.0; logger.info("excessive number of
             * negative Bs, increasing similarity restraint: " + kTbsim); }
             */
        }

        // set hydrogen based on bonded atom
        for (Atom a : atomarray) {
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
        int n = getNumberOfVariables();
        assert (x != null);
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : atomarray) {
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
        int index = nxyz + nb;
        for (ArrayList<Residue> list : refinementmodel.altresidues) {
            for (Residue r : list) {
                occ = x[index++];
                for (Atom a : r.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        a.setOccupancy(occ);
                    }
                }
            }
        }
        for (ArrayList<Molecule> list : refinementmodel.altmolecules) {
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
        double anisou1[];
        double anisou2[];
        double gradb;
        double det1, det2;
        double gradu[] = new double[6];
        double e = 0.0;

        for (Atom a : atomarray) {
            double biso = a.getTempFactor();
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }

            if (a.getAnisou() == null) {
                // isotropic B restraint

                // non-zero restraint: -kTln[Z], Z is ADP partition function
                e += -3.0 * kTbnonzero * Math.log(biso / (4.0 * Math.PI));
                gradb = -3.0 * kTbnonzero / biso;
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
                    e += kTbsim * Math.pow(bdiff, 2.0);
                    gradb = 2.0 * kTbsim * bdiff;
                    a1.addToTempFactorGradient(gradb);
                    a2.addToTempFactorGradient(-gradb);
                }
            } else {
                // anisotropic B restraint
                anisou1 = a.getAnisou();
                det1 = determinant3(anisou1);

                // non-zero restraint: -kTln[Z], Z is ADP partition function
                e += u2b(-kTbnonzero * Math.log(det1 * eightpi2 * Math.PI));
                gradu[0] = u2b(-kTbnonzero * ((anisou1[1] * anisou1[2] - anisou1[5] * anisou1[5]) / det1));
                gradu[1] = u2b(-kTbnonzero * ((anisou1[0] * anisou1[2] - anisou1[4] * anisou1[4]) / det1));
                gradu[2] = u2b(-kTbnonzero * ((anisou1[0] * anisou1[1] - anisou1[3] * anisou1[3]) / det1));
                gradu[3] = u2b(-kTbnonzero * ((2.0 * (anisou1[4] * anisou1[5] - anisou1[3] * anisou1[2])) / det1));
                gradu[4] = u2b(-kTbnonzero * ((2.0 * (anisou1[3] * anisou1[5] - anisou1[4] * anisou1[1])) / det1));
                gradu[5] = u2b(-kTbnonzero * ((2.0 * (anisou1[3] * anisou1[4] - anisou1[5] * anisou1[0])) / det1));
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
                    if (a2.getAnisou() == null) {
                        continue;
                    }
                    if (!a1.getAltLoc().equals(' ')
                            && !a1.getAltLoc().equals('A')
                            && a2.getAltLoc().equals(' ')) {
                        continue;
                    }

                    anisou2 = a2.getAnisou();
                    det2 = determinant3(anisou2);
                    bdiff = det1 - det2;
                    e += eightpi23 * kTbsim * Math.pow(bdiff, 2.0);
                    gradb = eightpi23 * 2.0 * kTbsim * bdiff;

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
        double mass[] = new double[nxyz + nb + nocc];
        int i = 0;
        if (refinexyz) {
            for (Atom a : atomarray) {
                double m = a.getMass();
                mass[i++] = m;
                mass[i++] = m;
                mass[i++] = m;
            }
        }

        if (refineb) {
            for (int j = i; j < nxyz + nb; i++, j++) {
                mass[j] = bmass;
            }
        }

        if (refineocc) {
            for (int j = i; j < nxyz + nb + nocc; i++, j++) {
                mass[j] = occmass;
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
        return nxyz + nb + nocc;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            diffractiondata.setLambda(lambda);
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
        diffractiondata.setLambda(1.0);
        // compute new structure factors
        diffractiondata.computeAtomicDensity();

        // compute crystal likelihood
        double e = diffractiondata.computeLikelihood();
        diffractiondata.setLambda(lambda);

        return e;
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
        // compute the crystal gradients
        diffractiondata.computeAtomicGradients(refinementMode);

        // pack gradients into gradient array
        getXYZGradients(gradient);
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        VARIABLE_TYPE vtypes[] = new VARIABLE_TYPE[nxyz + nb + nocc];
        int i = 0;
        if (refinexyz) {
            for (Atom a : atomarray) {
                vtypes[i++] = VARIABLE_TYPE.X;
                vtypes[i++] = VARIABLE_TYPE.Y;
                vtypes[i++] = VARIABLE_TYPE.Z;
            }
        }

        if (refineb) {
            for (int j = i; j < nxyz + nb; i++, j++) {
                vtypes[j] = VARIABLE_TYPE.OTHER;
            }
        }

        if (refineocc) {
            for (int j = i; j < nxyz + nb + nocc; i++, j++) {
                vtypes[j] = VARIABLE_TYPE.OTHER;
            }
        }
        return vtypes;
    }

    /*
     * RESPA setup
     */
    @Override
    public void setEnergyTermState(STATE state) {
        switch (state) {
            case FAST:
                xrayterms = false;
                restraintterms = true;
                break;
            case SLOW:
                xrayterms = true;
                restraintterms = false;
                break;
            default:
                xrayterms = true;
                restraintterms = true;
        }
    }
}
