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
import ffx.xray.RefinementMinimize.RefinementMode;
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
    private double kT32;

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
        kT32 = 1.5 * kB * temp * refinementdata.bresweight;

        logger.info("total B restraint weight: " + temp * refinementdata.bresweight);
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

                // compute the crystal gradients (requires inverse FFT)
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
                setBFactors(x, 0);

                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);

                // compute crystal likelihood
                e = xraystructure.sigmaaminimize.calculateLikelihood();

                // compute the crystal gradients (requires inverse FFT)
                crs_fc.computeAtomicGradients(refinementdata.dfc,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);
                crs_fs.computeAtomicGradients(refinementdata.dfs,
                        refinementdata.freer, refinementdata.rfreeflag,
                        refinementMode);

                // add B restraints
                e += getBFactorRestraints();

                // pack gradients into gradient array
                getBFactorGradients(g, 0);

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
                setBFactors(x, nxyz);

                // compute new structure factors
                crs_fc.computeDensity(refinementdata.fc);
                crs_fs.computeDensity(refinementdata.fs);

                // compute crystal likelihood
                e = xraystructure.sigmaaminimize.calculateLikelihood();

                // compute the crystal gradients (requires inverse FFT)
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
                getBFactorGradients(g, nxyz);

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

    public void getBFactorGradients(double g[], int offset) {
        assert (g != null);
        double grad[];
        int index = offset;
        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() == null) {
                g[index++] = a.getTempFactorGradient();
            } else {
                grad = a.getAnisouGradient();
                g[index++] = grad[0];
                g[index++] = grad[1];
                g[index++] = grad[2];
                g[index++] = grad[3];
                g[index++] = grad[4];
                g[index++] = grad[5];
            }
        }
    }

    public void getBFactors(double x[], int offset) {
        double anisou[];
        int index = offset;
        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() == null) {
                x[index++] = a.getTempFactor();
            } else {
                anisou = a.getAnisou();
                x[index++] = anisou[0];
                x[index++] = anisou[1];
                x[index++] = anisou[2];
                x[index++] = anisou[3];
                x[index++] = anisou[4];
                x[index++] = anisou[5];
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
        int n = getNumberOfVariables();
        assert (x != null);
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : atomarray) {
            a.getXYZ(xyz);
            x[index++] = xyz[0];
            x[index++] = xyz[1];
            x[index++] = xyz[2];
        }

        return x;
    }

    public void setBFactors(double x[], int offset) {
        double tmpanisou[] = new double[6];
        int index = offset;
        for (Atom a : atomarray) {
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() == null) {
                a.setTempFactor(x[index++]);
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
                    for (int i = 0; i < 6; i++) {
                        anisou[i] = tmpanisou[i];
                    }
                    det = Math.pow(det, 0.3333);
                    a.setTempFactor(u2b(det));
                }
            }
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
            xyz[2] = x[index++];
            xyz[2] = x[index++];
            a.moveTo(xyz);
        }
    }

    public double getBFactorRestraints() {
        double c = kT32 * Math.log(4.0 * Math.PI);
        double pi6 = 512.0 * Math.pow(Math.PI, 6.0);
        double anisou[];
        double banisou[] = new double[6];
        double gradb;
        double det;
        double gradu[] = new double[6];
        double e = 0.0;

        for (Atom a : atomarray) {
            double biso = a.getTempFactor();
            // ignore hydrogens!!!
            if (a.getAtomicNumber() == 1) {
                continue;
            }
            if (a.getAnisou() == null) {
                e += -kT32 * Math.log(Math.pow(biso, 3.0)) + c;
                gradb = -3.0 * kT32 / biso;
                a.addToTempFactorGradient(gradb);
            } else {
                anisou = a.getAnisou();
                for (int i = 0; i < 6; i++) {
                    banisou[i] = u2b(anisou[i]);
                }
                det = determinant3(banisou);
                e += -kT32 * Math.log(det) + c;
                gradu[0] = -kT32 * ((pi6 * (-banisou[0] * banisou[0] + banisou[1] * banisou[2])) / det);
                gradu[1] = -kT32 * ((pi6 * (-banisou[4] * banisou[4] + banisou[0] * banisou[2])) / det);
                gradu[2] = -kT32 * ((pi6 * (-banisou[3] * banisou[3] + banisou[0] * banisou[1])) / det);
                gradu[3] = -kT32 * ((2.0 * pi6 * (banisou[4] * banisou[5] - banisou[3] * banisou[2])) / det);
                gradu[4] = -kT32 * ((2.0 * pi6 * (banisou[3] * banisou[5] - banisou[4] * banisou[1])) / det);
                gradu[5] = -kT32 * ((2.0 * pi6 * (banisou[3] * banisou[4] - banisou[5] * banisou[0])) / det);
                for (int i = 0; i < 6; i++) {
                    gradu[i] = b2u(gradu[i]);
                }
                a.addToAnisouGradient(gradu);
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
            for (Atom a : atomarray) {
                // ignore hydrogens!!!
                if (a.getAtomicNumber() == 1) {
                    continue;
                }
                if (a.getAnisou() == null) {
                    mass[i++] = bmass;
                } else {
                    mass[i++] = bmass;
                    mass[i++] = bmass;
                    mass[i++] = bmass;
                    mass[i++] = bmass;
                    mass[i++] = bmass;
                    mass[i++] = bmass;
                }
            }
        }
        return mass;
    }

    @Override
    public int getNumberOfVariables() {
        return nxyz + nb + nocc;
    }
}
