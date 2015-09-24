/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.potential;

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;

import static ffx.crystal.SpaceGroup.CrystalSystem.CUBIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.HEXAGONAL;
import static ffx.crystal.SpaceGroup.CrystalSystem.MONOCLINIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.ORTHORHOMBIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.TETRAGONAL;
import static ffx.crystal.SpaceGroup.CrystalSystem.TRICLINIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.TRIGONAL;

/**
 * This class computes the energy and Cartesian coordinate gradient, plus finite
 * difference derivatives of lattice parameters.
 *
 * @author Jooyeon Park
 */
public class XtalEnergy implements Potential {

    /**
     * The logger.
     */
    private static final Logger logger = Logger.getLogger(XtalEnergy.class.getName());

    private final ForceFieldEnergy forceFieldEnergy;
    private final Atom activeAtoms[];
    private final int nActive;

    private final double xyz[];
    private final double gr[];
    private final int nParams;

    private final Crystal crystal;
    private final VARIABLE_TYPE type[];
    private final double mass[];
    private final Crystal unitCell;
    private double scaling[];
    private double totalEnergy;

    public XtalEnergy(ForceFieldEnergy forceFieldEnergy, MolecularAssembly molecularAssembly) {
        this.forceFieldEnergy = forceFieldEnergy;
        Atom atoms[] = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;

        int n = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                n++;
            }
        }
        nActive = n;

        activeAtoms = new Atom[nActive];
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                activeAtoms[index++] = a;
            }
        }

        nParams = 3 * nActive + 6;
        crystal = forceFieldEnergy.getCrystal();
        unitCell = crystal.getUnitCell();
        xyz = new double[3 * nActive];
        gr = new double[3 * nActive];
        type = new VARIABLE_TYPE[nParams];
        mass = new double[nParams];

        index = 0;
        for (int i = 0; i < nActive; i++) {
            double m = activeAtoms[i].getMass();
            mass[index] = m;
            mass[index + 1] = m;
            mass[index + 2] = m;
            type[index] = VARIABLE_TYPE.X;
            type[index + 1] = VARIABLE_TYPE.Y;
            type[index + 2] = VARIABLE_TYPE.Z;
            index += 3;
        }
        for (int i = nActive * 3; i < nActive * 3 + 6; i++) {
            mass[i] = 1.0;
            type[i] = VARIABLE_TYPE.OTHER;
        }
    }

    @Override
    public double energy(double[] x) {
        /**
         * Un-scale coordinates if applicable.
         */
        if (scaling != null) {
            for (int i = 0; i < nParams; i++) {
                x[i] /= scaling[i];
            }
        }

        /**
         * Set atomic coordinates & lattice parameters.
         */
        setCoordinates(x);

        totalEnergy = forceFieldEnergy.energy(false, false);

        /**
         * Scale coordinates if applicable.
         */
        if (scaling != null) {
            for (int i = 0; i < nParams; i++) {
                x[i] *= scaling[i];
            }
        }

        return totalEnergy;

    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {
        /**
         * Un-scale coordinates if applicable.
         */
        if (scaling != null) {
            for (int i = 0; i < nParams; i++) {
                x[i] /= scaling[i];
            }
        }

        /**
         * Set atomic coordinates & lattice parameters.
         */
        setCoordinates(x);

        /**
         * Calculate system energy and Cartesian coordinate gradient.
         */
        double e = forceFieldEnergy.energyAndGradient(xyz, gr);

        /**
         * Both coordinates and gradient are scaled if applicable.
         */
        packGradient(x, g);

        /**
         * Calculate finite-difference partial derivatives of lattice
         * parameters.
         */
        unitCellParameterDerivatives(x, g);

        totalEnergy = e;

        return totalEnergy;

    }

    /**
     *
     * @param x
     * @param g
     */
    private void unitCellParameterDerivatives(double x[], double g[]) {

        double eps = 1.0e-5;
        double deps = Math.toDegrees(eps);

        double a = unitCell.a;
        double b = unitCell.b;
        double c = unitCell.c;
        double alpha = unitCell.alpha;
        double beta = unitCell.beta;
        double gamma = unitCell.gamma;

        int index = 3 * nActive;
        switch (crystal.spaceGroup.crystalSystem) {
            case TRICLINIC:
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, deps);
                index++;
                g[index] = finiteDifference(x, index, deps);
                index++;
                g[index] = finiteDifference(x, index, deps);
                break;
            case MONOCLINIC:
                // alpha == gamma == 90.0
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = 0.0;
                index++;
                g[index] = finiteDifference(x, index, deps);
                index++;
                g[index] = 0.0;
                break;
            case ORTHORHOMBIC:
                // alpha == beta == gamma == 90.0
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                break;
            case TETRAGONAL:
                // a == b && alpha == beta == gamma == 90.0
                g[index] = finiteDifference2(x, index, index + 1, eps);
                index++;
                g[index] = g[index - 1];
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                break;
            case TRIGONAL:
                if (a == b && b == c && alpha == beta && beta == gamma) {
                    // Rhombohedral axes, primitive cell.
                    g[index] = finiteDifference3(x, index, index + 1, index + 2, eps);
                    index++;
                    g[index] = g[index - 1];
                    index++;
                    g[index] = g[index - 2];
                    index++;
                    g[index] = finiteDifference3(x, index, index + 1, index + 2, deps);
                    index++;
                    g[index] = g[index - 1];
                    index++;
                    g[index] = g[index - 2];

                } else if (a == b && alpha == 90.0 && beta == 90.0 && gamma == 120.0) {
                    // Hexagonal axes, triple obverse cell.
                    g[index] = finiteDifference2(x, index, index + 1, eps);
                    index++;
                    g[index] = g[index - 1];
                    index++;
                    g[index] = finiteDifference(x, index, eps);
                    index++;
                    g[index] = 0.0;
                    index++;
                    g[index] = 0.0;
                    index++;
                    g[index] = 0.0;

                }
                break;
            case HEXAGONAL:
                // a == b && alpha == beta == 90.0 && gamma == 120.0
                g[index] = finiteDifference2(x, index, index + 1, eps);
                index++;
                g[index] = g[index - 1];
                index++;
                g[index] = finiteDifference(x, index, eps);
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                break;
            case CUBIC:
                // a == b == c && alpha == beta == gamma == 90.0
                g[index] = finiteDifference3(x, index, index + 1, index + 2, eps);
                index++;
                g[index] = g[index - 1];
                index++;
                g[index] = g[index - 2];
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                index++;
                g[index] = 0.0;
                break;
        }
        /**
         * Scale finite-difference partial derivatives of lattice parameters.
         */
        if (scaling != null) {
            index = 3 * nActive;
            g[index] /= scaling[index];
            index++;
            g[index] /= scaling[index];
            index++;
            g[index] /= scaling[index];
            index++;
            g[index] /= scaling[index];
            index++;
            g[index] /= scaling[index];
            index++;
            g[index] /= scaling[index];
        }
    }

    /**
     * Calculate finite-difference derivative for any parameter.
     */
    private double finiteDifference(double[] x, int index, double eps) {
        double scale = 1.0;
        if (scaling != null) {
            scale = scaling[index];
        }
        double xoriginal = x[index];
        double param = x[index] / scale;

        x[index] = (param + eps / 2.0) * scale;
        double ePlus = energy(x);
        x[index] = (param - eps / 2.0) * scale;
        double eMinus = energy(x);

        x[index] = xoriginal;

        return (ePlus - eMinus) / eps;
    }

    /**
     * @param x
     * @param index1
     * @param index2
     * @param eps
     *
     * @return finite-difference derivative.
     */
    private double finiteDifference2(double[] x, int index1, int index2, double eps) {
        double scale1 = 1.0;
        double scale2 = 1.0;

        if (scaling != null) {
            scale1 = scaling[index1];
            scale2 = scaling[index2];
        }

        double param1 = x[index1] / scale1;
        double param2 = x[index2] / scale2;

        x[index1] = (param1 + eps / 2.0) * scale1;
        x[index2] = (param2 + eps / 2.0) * scale2;
        double ePlus = energy(x);
        x[index1] = (param1 - eps / 2.0) * scale1;
        x[index2] = (param2 - eps / 2.0) * scale2;
        double eMinus = energy(x);

        x[index1] = param1 * scale1;
        x[index2] = param2 * scale2;

        return (ePlus - eMinus) / eps;
    }

    /**
     * @param x
     * @param index1
     * @param index2
     * @param index3
     * @param eps
     * @return finite-difference derivative.
     */
    private double finiteDifference3(double[] x, int index1, int index2, int index3, double eps) {
        double scale1 = 1.0;
        double scale2 = 1.0;
        double scale3 = 1.0;

        if (scaling != null) {
            scale1 = scaling[index1];
            scale2 = scaling[index2];
            scale3 = scaling[index3];
        }

        double param1 = x[index1] / scale1;
        double param2 = x[index2] / scale2;
        double param3 = x[index3] / scale3;

        x[index1] = (param1 + eps / 2.0) * scale1;
        x[index2] = (param2 + eps / 2.0) * scale2;
        x[index3] = (param3 + eps / 2.0) * scale3;
        double ePlus = energy(x);
        x[index1] = (param1 - eps / 2.0) * scale1;
        x[index2] = (param2 - eps / 2.0) * scale2;
        x[index3] = (param3 - eps / 2.0) * scale3;
        double eMinus = energy(x);

        x[index1] = param1 * scale1;
        x[index2] = param2 * scale2;
        x[index2] = param2 * scale2;

        return (ePlus - eMinus) / eps;
    }

    /**
     * Apply scaling for the optimizer if applicable.
     *
     * @param x
     * @param g
     */
    private void packGradient(double x[], double g[]) {
        // Scale fractional coordinates and gradient.
        if (scaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                g[i] /= scaling[i];
                x[i] *= scaling[i];
            }
        }

    }

    /**
     * Sets atomic coordinates and lattice parameters.
     *
     * @param x First 3*nActive parameters are coordinates, next 6 are x
     * parameters.
     */
    private void setCoordinates(double x[]) {
        assert (x != null);

        int index = nActive * 3;
        double a = x[index];
        double b = x[index + 1];
        double c = x[index + 2];
        double alpha = x[index + 3];
        double beta = x[index + 4];
        double gamma = x[index + 5];

        switch (crystal.spaceGroup.crystalSystem) {
            case TRICLINIC:
                break;
            case MONOCLINIC:
                // alpha == gamma == 90.0
                alpha = 90.0;
                gamma = 90.0;
                break;
            case ORTHORHOMBIC:
                // alpha == beta == gamma == 90.0
                alpha = 90.0;
                beta = 90.0;
                gamma = 90.0;
                break;
            case TETRAGONAL:
                // a == b && alpha == beta == gamma == 90.0
                double temp = (a + b) / 2.0;
                a = temp;
                b = temp;
                alpha = 90.0;
                beta = 90.0;
                gamma = 90.0;
                break;
            case TRIGONAL:
                if (a == b && b == c && alpha == beta && beta == gamma) {
                    temp = (a + b + c) / 3.0;
                    a = temp;
                    b = temp;
                    c = temp;
                    temp = (alpha + beta + gamma) / 3.0;
                    alpha = temp;
                    beta = temp;
                    gamma = temp;
                } else if (a == b && alpha == 90.0 && beta == 90.0 && gamma == 120.0) {
                    // Hexagonal axes, triple obverse cell.
                    temp = (a + b) / 2.0;
                    a = temp;
                    b = temp;
                    alpha = 90.0;
                    beta = 90.0;
                    gamma = 120.0;
                }
                break;
            case HEXAGONAL:
                // a == b && alpha == beta == 90.0 && gamma == 120.0
                temp = (a + b) / 2.0;
                a = temp;
                b = temp;
                alpha = 90.0;
                beta = 90.0;
                gamma = 120.0;
                break;
            case CUBIC:
                // a == b == c && alpha == beta == gamma == 90.0
                temp = (a + b + c) / 3.0;
                a = temp;
                b = temp;
                c = temp;
                alpha = 90.0;
                beta = 90.0;
                gamma = 90.0;
                break;
        }
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
        forceFieldEnergy.setCrystal(crystal);

        index = 0;
        for (int i = 0; i < nActive; i++) {
            Atom atom = activeAtoms[i];
            double xx = x[index];
            double yy = x[index + 1];
            double zz = x[index + 2];
            xyz[index] = xx;
            xyz[index + 1] = yy;
            xyz[index + 2] = zz;
            index += 3;
            atom.moveTo(xx, yy, zz);
        }
    }

    @Override
    public void setScaling(double[] scaling) {
        this.scaling = scaling;
    }

    @Override
    public double[] getScaling() {
        return scaling;
    }

    @Override
    public double[] getCoordinates(double[] x) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        int index = 0;
        for (int i = 0; i < nActive; i++) {
            Atom a = activeAtoms[i];
            x[index] = a.getX();
            index++;
            x[index] = a.getY();
            index++;
            x[index] = a.getZ();
            index++;
        }
        x[index] = unitCell.a;
        index++;
        x[index] = unitCell.b;
        index++;
        x[index] = unitCell.c;
        index++;
        x[index] = unitCell.alpha;
        index++;
        x[index] = unitCell.beta;
        index++;
        x[index] = unitCell.gamma;
        return x;
    }

    @Override
    public double[] getMass() {
        return mass;
    }

    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    @Override
    public int getNumberOfVariables() {
        return nParams;
    }

    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return type;
    }

    @Override
    public void setEnergyTermState(STATE state) {
        forceFieldEnergy.setEnergyTermState(state);
    }
    
    @Override
    public void reInit() {
        throw new UnsupportedOperationException(String.format(" No reInit method defined for %s", XtalEnergy.class.toString()));
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
