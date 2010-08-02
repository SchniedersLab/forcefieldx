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
import ffx.xray.RefinementMinimize.RefinementMode;
import java.util.ArrayList;
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
public class XRayEnergy implements Potential {

    private static final Logger logger = Logger.getLogger(XRayEnergy.class.getName());
    private final XRayStructure xraystructure;
    private final CrystalReciprocalSpace crs_fc;
    private final CrystalReciprocalSpace crs_fs;
    private final RefinementData refinementdata;
    private final Atom atomarray[];
    private final int nAtoms;
    private int nxyz;
    private int nb;
    private int nocc;
    private RefinementMode refinementMode;
    protected double[] optimizationScaling = null;
    private double bmass;
    private double temp = 0.1;
    private double kTbnonzero;
    private double kTbsim;

    public XRayEnergy(XRayStructure xraystructure, int nxyz, int nb, int nocc,
            RefinementMode refinementmode) {
        this.xraystructure = xraystructure;
        this.refinementdata = xraystructure.refinementdata;
        this.crs_fc = refinementdata.crs_fc;
        this.crs_fs = refinementdata.crs_fs;
        this.refinementMode = refinementmode;
        this.atomarray = xraystructure.atomarray;
        this.nAtoms = atomarray.length;
        this.nxyz = nxyz;
        this.nb = nb;
        this.nocc = nocc;

        bmass = refinementdata.bmass;
        kTbnonzero = 1.5 * kB * temp * refinementdata.bnonzeroweight;
        kTbsim = kB * temp * refinementdata.bsimweight;

        if (refinementMode == RefinementMode.BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS) {
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
        switch (refinementMode) {
            case COORDINATES:
                for (Atom a : atomarray) {
                    a.setXYZGradient(0.0, 0.0, 0.0);
                }
                // update coordinates
                crs_fc.setCoordinates(x);
                crs_fs.setCoordinates(x);

                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);

                // compute crystal likelihood
                e = xraystructure.sigmaaminimize.calculateLikelihood();

                // compute the crystal gradients
                crs_fc.computeAtomicGradients(refinementdata.dfc,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);
                crs_fs.computeAtomicGradients(refinementdata.dfs,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);

                // pack gradients into gradient array
                getXYZGradients(g);

                break;
            case BFACTORS:
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

                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);

                // compute crystal likelihood
                e = xraystructure.sigmaaminimize.calculateLikelihood();

                // compute the crystal gradients
                crs_fc.computeAtomicGradients(refinementdata.dfc,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);
                crs_fs.computeAtomicGradients(refinementdata.dfs,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);

                // add B restraints
                e += getBFactorRestraints();

                // pack gradients into gradient array
                getBFactorGradients(g);

                break;
            case COORDINATES_AND_BFACTORS:
                for (Atom a : atomarray) {
                    a.setXYZGradient(0.0, 0.0, 0.0);
                    a.setTempFactorGradient(0.0);
                    if (a.getAnisouGradient() == null) {
                        double ganisou[] = new double[6];
                        a.setAnisouGradient(ganisou);
                    } else {
                        double ganisou[] = a.getAnisouGradient();
                        ganisou[0] = ganisou[1] = ganisou[2] = 0.0;
                        ganisou[3] = ganisou[4] = ganisou[5] = 0.0;
                    }
                }
                // update coordinates
                crs_fc.setCoordinates(x);
                crs_fs.setCoordinates(x);
                // update B factors
                setBFactors(x);

                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);

                // compute crystal likelihood
                e = xraystructure.sigmaaminimize.calculateLikelihood();

                // compute the crystal gradients
                crs_fc.computeAtomicGradients(refinementdata.dfc,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);
                crs_fs.computeAtomicGradients(refinementdata.dfs,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);

                // add B restraints
                e += getBFactorRestraints();

                // pack gradients into gradient array
                getXYZGradients(g);
                getBFactorGradients(g);

                break;
            default:
                String message = "refinement mode not implemented.";
                logger.log(Level.SEVERE, message);
                break;
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
    }

    public int getNXYZ() {
        return nxyz;
    }

    public void setNXYZ(int nxyz) {
        this.nxyz = nxyz;
    }

    public int getNB() {
        return nb;
    }

    public void setNB(int nb) {
        this.nb = nb;
    }

    public int getNOcc() {
        return nocc;
    }

    public void setNOcc(int nocc) {
        this.nocc = nocc;
    }

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

        if (refinementMode == RefinementMode.COORDINATES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS) {
            for (Atom a : atomarray) {
                a.getXYZ(xyz);
                x[index++] = xyz[0];
                x[index++] = xyz[1];
                x[index++] = xyz[2];
            }
        }

        if (refinementMode == RefinementMode.BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS) {
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

        return x;
    }

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
                e += -kTbnonzero * Math.log(Math.pow(biso, 3.0)) + c;
                gradb = -3.0 * kTbnonzero / biso;
                a.addToTempFactorGradient(gradb);

                if (refinementdata.bsimweight > 0.0) {
                    ArrayList<Bond> bonds = a.getBonds();
                    for (Bond b : bonds) {
                        if (a.compareTo(b.getAtom(0)) == 0) {
                            a1 = b.getAtom(0);
                            a2 = b.getAtom(1);
                        } else {
                            a1 = b.getAtom(1);
                            a2 = b.getAtom(0);
                        }
                        // ignore hydrogens!!!
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
                }
            } else {
                anisou = a.getAnisou();
                for (int i = 0; i < 6; i++) {
                    banisou1[i] = u2b(anisou[i]);
                }
                det1 = determinant3(banisou1);
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

                if (refinementdata.bsimweight > 0.0) {
                    ArrayList<Bond> bonds = a.getBonds();
                    for (Bond b : bonds) {
                        if (a.compareTo(b.getAtom(0)) == 0) {
                            a1 = b.getAtom(0);
                            a2 = b.getAtom(1);
                        } else {
                            a1 = b.getAtom(1);
                            a2 = b.getAtom(0);
                        }
                        // ignore hydrogens!!!
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
        if (refinementMode == RefinementMode.COORDINATES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS) {
            for (Atom a : atomarray) {
                double m = a.getMass();
                mass[i++] = m;
                mass[i++] = m;
                mass[i++] = m;
            }
        }

        if (refinementMode == RefinementMode.BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS) {
            for (int j = i; j < nxyz + nb + nocc; j++) {
                mass[j] = bmass;
            }
        }
        return mass;
    }

    @Override
    public int getNumberOfVariables() {
        return nxyz + nb + nocc;
    }
}
