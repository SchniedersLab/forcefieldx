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

import static ffx.algorithms.Thermostat.kB;
import static ffx.numerics.VectorMath.determinant3;
import static ffx.numerics.VectorMath.b2u;
import static ffx.numerics.VectorMath.u2b;

import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.RefinementMinimize.RefinementMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class XRayEnergy implements Potential {

    private static final Logger logger = Logger.getLogger(XRayEnergy.class.getName());
    private final DiffractionData diffractiondata;
    private final RefinementData refinementdata;
    private final Atom atomarray[];
    private final int nAtoms;
    private int nxyz;
    private int nb;
    private int nocc;
    private RefinementMode refinementMode;
    private boolean refinexyz = false;
    private boolean refineocc = false;
    private boolean refineb = false;
    protected double[] optimizationScaling = null;
    private double bmass;
    private double temp = 0.1;
    private double kTbnonzero;
    private double kTbsim;
    private double occmass;

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
        this.refinementdata = diffractiondata.refinementdata[0];
        this.refinementMode = refinementmode;
        this.atomarray = diffractiondata.atomarray;
        this.nAtoms = atomarray.length;
        this.nxyz = nxyz;
        this.nb = nb;
        this.nocc = nocc;

        bmass = refinementdata.bmass;
        kTbnonzero = 1.5 * kB * temp * refinementdata.bnonzeroweight;
        kTbsim = 0.01 * kB * temp * refinementdata.bsimweight;
        occmass = refinementdata.occmass;

        setRefinementBooleans();

        if (refineb) {
            logger.info("total B non-zero restraint weight: " + temp * refinementdata.bnonzeroweight);
            logger.info("total B similarity restraint weight: " + temp * refinementdata.bsimweight);
        }
    }

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
        return e;
    }

    public RefinementMode getRefinementMode() {
        return refinementMode;
    }

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
     * @return the number of xyz parameters
     */
    public int getNXYZ() {
        return nxyz;
    }

    /**
     * set the number of xyz parameters
     * @param nxyz requested number of xyz parameters
     */
    public void setNXYZ(int nxyz) {
        this.nxyz = nxyz;
    }

    /**
     * get the number of B factor parameters being fit
     * @return the number of B factor parameters
     */
    public int getNB() {
        return nb;
    }

    /**
     * set the number of B factor parameters
     * @param nb requested number of B factor parameters
     */
    public void setNB(int nb) {
        this.nb = nb;
    }

    /**
     * get the number of occupancy parameters being fit
     * @return the number of occupancy parameters
     */
    public int getNOcc() {
        return nocc;
    }

    /**
     * set the number of occupancy parameters
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
        int nres = refinementdata.nresiduebfactor + 1;
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
            } else if (refinementdata.residuebfactor) {
                if (resnum != a.getResidueNumber()) {
                    if (nres >= refinementdata.nresiduebfactor) {
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
        for (ArrayList<Residue> list : diffractiondata.altresidues) {
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
        for (ArrayList<Molecule> list : diffractiondata.altmolecules) {
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
            int nres = refinementdata.nresiduebfactor + 1;
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
                } else if (refinementdata.residuebfactor) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= refinementdata.nresiduebfactor) {
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

            if (refinementdata.residuebfactor) {
                if (nat > 1) {
                    x[index] /= nat;
                }
            }
        }

        if (refineocc) {
            for (ArrayList<Residue> list : diffractiondata.altresidues) {
                for (Residue r : list) {
                    for (Atom a : r.getAtomList()) {
                        if (a.getOccupancy() < 1.0) {
                            x[index++] = a.getOccupancy();
                            break;
                        }
                    }
                }
            }
            for (ArrayList<Molecule> list : diffractiondata.altmolecules) {
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
        int nres = refinementdata.nresiduebfactor + 1;
        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() == null) {
                double biso = x[index];
                if (refinementdata.residuebfactor) {
                    if (resnum != a.getResidueNumber()) {
                        if (nres >= refinementdata.nresiduebfactor) {
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
                    a.setTempFactor(0.1);
                    anisou[0] = anisou[1] = anisou[2] = b2u(0.1);
                    anisou[3] = anisou[4] = anisou[5] = 0.0;
                    if (nneg < 5) {
                        logger.info("anisotropic atom: " + a.toString() + " negative ANISOU");
                    }
                }
            }
        }

        if (nneg > 0) {
            logger.info(nneg + " of " + nAtoms
                    + " atoms with negative B factors! Attempting to correct....");
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
        for (ArrayList<Residue> list : diffractiondata.altresidues) {
            for (Residue r : list) {
                occ = x[index++];
                for (Atom a : r.getAtomList()) {
                    if (a.getOccupancy() < 1.0) {
                        a.setOccupancy(occ);
                    }
                }
            }
        }
        for (ArrayList<Molecule> list : diffractiondata.altmolecules) {
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
        double c = kTbnonzero * Math.log(4.0 * Math.PI);
        double pi6 = 512.0 * Math.pow(Math.PI, 6.0);
        double anisou[];
        double banisou1[] = new double[6];
        double banisou2[] = new double[6];
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

                // non-zero restraint
                e += -kTbnonzero * Math.log(Math.pow(biso, 3.0)) + c;
                gradb = -3.0 * kTbnonzero / biso;
                a.addToTempFactorGradient(gradb);

                // similarity / harmonic restraint
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
                    bdiff = b1 - b2;
                    e += kTbsim * Math.pow(bdiff, 2.0);
                    gradb = 2.0 * kTbsim * bdiff;
                    a1.addToTempFactorGradient(gradb);
                    a2.addToTempFactorGradient(-gradb);
                }
            } else {
                // anisotropic B restraint
                anisou = a.getAnisou();
                for (int i = 0; i < 6; i++) {
                    banisou1[i] = u2b(anisou[i]);
                }
                det1 = determinant3(banisou1);

                // non-zero restraint
                e += -kTbnonzero * Math.log(det1) + c;
                gradu[0] = -kTbnonzero * ((pi6 * (-banisou1[0] * banisou1[0] + banisou1[1] * banisou1[2])) / det1);
                gradu[1] = -kTbnonzero * ((pi6 * (-banisou1[4] * banisou1[4] + banisou1[0] * banisou1[2])) / det1);
                gradu[2] = -kTbnonzero * ((pi6 * (-banisou1[3] * banisou1[3] + banisou1[0] * banisou1[1])) / det1);
                gradu[3] = -kTbnonzero * ((2.0 * pi6 * (banisou1[4] * banisou1[5] - banisou1[3] * banisou1[2])) / det1);
                gradu[4] = -kTbnonzero * ((2.0 * pi6 * (banisou1[3] * banisou1[5] - banisou1[4] * banisou1[1])) / det1);
                gradu[5] = -kTbnonzero * ((2.0 * pi6 * (banisou1[3] * banisou1[4] - banisou1[5] * banisou1[0])) / det1);
                for (int i = 0; i < 6; i++) {
                    gradu[i] = b2u(gradu[i]);
                }
                a.addToAnisouGradient(gradu);

                // similarity / harmonic restraint
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

                    anisou = a2.getAnisou();
                    for (int i = 0; i < 6; i++) {
                        banisou2[i] = u2b(anisou[i]);
                    }
                    det2 = determinant3(banisou2);
                    bdiff = det1 - det2;
                    e += kTbsim * Math.pow(bdiff, 2.0);
                    gradb = 2.0 * kTbsim * bdiff;
                    gradu[0] = gradb * (banisou1[1] * banisou1[2] - banisou1[5] * banisou1[5]);
                    gradu[1] = gradb * (banisou1[0] * banisou1[2] - banisou1[4] * banisou1[4]);
                    gradu[2] = gradb * (banisou1[0] * banisou1[1] - banisou1[3] * banisou1[3]);
                    gradu[3] = gradb * (2.0 * banisou1[4] * banisou1[5] - banisou1[2] * banisou1[3]);
                    gradu[4] = gradb * (2.0 * banisou1[3] * banisou1[5] - banisou1[1] * banisou1[4]);
                    gradu[5] = gradb * (2.0 * banisou1[3] * banisou1[4] - banisou1[0] * banisou1[5]);
                    for (int i = 0; i < 6; i++) {
                        gradu[i] = b2u(gradu[i]);
                    }
                    a1.addToAnisouGradient(gradu);
                    gradu[0] = gradb * (banisou2[5] * banisou2[5] - banisou2[1] * banisou2[2]);
                    gradu[1] = gradb * (banisou2[4] * banisou2[4] - banisou2[0] * banisou2[2]);
                    gradu[2] = gradb * (banisou2[3] * banisou2[3] - banisou2[0] * banisou2[1]);
                    gradu[3] = gradb * (2.0 * banisou2[2] * banisou2[3] - banisou2[4] * banisou2[5]);
                    gradu[4] = gradb * (2.0 * banisou2[1] * banisou2[4] - banisou2[3] * banisou2[5]);
                    gradu[5] = gradb * (2.0 * banisou2[0] * banisou2[5] - banisou2[3] * banisou2[4]);
                    for (int i = 0; i < 6; i++) {
                        gradu[i] = b2u(gradu[i]);
                    }
                    a2.addToAnisouGradient(gradu);
                }
            }
        }
        return e;
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

    @Override
    public int getNumberOfVariables() {
        return nxyz + nb + nocc;
    }
}
