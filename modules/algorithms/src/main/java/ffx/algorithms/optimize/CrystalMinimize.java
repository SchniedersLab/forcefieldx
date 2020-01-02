//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.optimize;

import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.numerics.optimization.OptimizationListener;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.MolecularAssembly.FractionalMode;
import ffx.potential.XtalEnergy;
import ffx.utilities.Constants;

/**
 * Minimize the energy of a system to an RMS gradient per atom convergence criteria.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class CrystalMinimize extends Minimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(CrystalMinimize.class.getName());

    private final ForceFieldEnergy forceFieldEnergy;
    private Crystal crystal;
    private Crystal unitCell;

    /**
     * <p>
     * Constructor for Minimize.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param xtalEnergy        a {@link ffx.potential.XtalEnergy} object.
     * @param algorithmListener a {@link ffx.algorithms.AlgorithmListener} object.
     */
    public CrystalMinimize(MolecularAssembly molecularAssembly, XtalEnergy xtalEnergy,
                           AlgorithmListener algorithmListener) {
        super(molecularAssembly, xtalEnergy, algorithmListener);
        crystal = molecularAssembly.getCrystal();
        unitCell = crystal.getUnitCell();
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();

        for (int i = 0; i < n - 6; i += 3) {
            scaling[i] = 12.0 * crystal.a;
            scaling[i + 1] = 12.0 * crystal.b;
            scaling[i + 2] = 12.0 * crystal.c;
        }
        scaling[n - 6] = 4.0 * sqrt(crystal.a);
        scaling[n - 5] = 4.0 * sqrt(crystal.b);
        scaling[n - 4] = 4.0 * sqrt(crystal.c);
        scaling[n - 3] = 0.02 * sqrt(crystal.volume);
        scaling[n - 2] = 0.02 * sqrt(crystal.volume);
        scaling[n - 1] = 0.02 * sqrt(crystal.volume);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param m             The number of previous steps used to estimate the Hessian.
     * @param maxIterations The maximum number of minimization steps.
     * @param eps           The minimization convergence criteria.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(int m, double eps, int maxIterations) {
        super.minimize(m, eps, maxIterations);

        crystal = molecularAssembly.getCrystal();
        logger.info("\n Final lattice parameters" + crystal);

        return potential;
    }

    /**
     * Print out the partial derivatives of the Energy with respect to
     * components of the 3 vectors that define the primitive cell.
     */
    public void printTensor() {
        computeStressTensor(true);
    }

    /**
     * <p>computeStressTensor.</p>
     *
     * @param verbose a boolean.
     */
    public void computeElasticityTensor(boolean verbose) {

        double[] xOrig = new double[forceFieldEnergy.getNumberOfVariables()];
        forceFieldEnergy.getCoordinates(xOrig);

        forceFieldEnergy.energy(xOrig, verbose);

        FractionalMode currentFractionalMode = molecularAssembly.getFractionalMode();

        molecularAssembly.setFractionalMode(FractionalMode.MOLECULE);
        molecularAssembly.computeFractionalCoordinates();

        // Conversion from kcal/mol/Ang^3 to Atm.
        double PRESCON = 6.85684112e4;

        // Conversion from Atm to GPa.
        double GPA = 101325 * 1.0e-9;

        double volume = crystal.getUnitCell().volume / crystal.getUnitCell().spaceGroup.getNumberOfSymOps();

        boolean currentCheckRestrictions = crystal.getCheckRestrictions();

        crystal.setCheckRestrictions(false);

        logger.info(" The Strain and Stress Tensor code is under development.");

        double delta = 5.0e-3;
        double eps = 2.0e-3;

        double c11 = 0.0, c22 = 0.0, c33 = 0.0;
        double c12 = 0.0, c13 = 0.0, c23 = 0.0;
        double c14 = 0.0, c15 = 0.0, c16 = 0.0;
        double c24 = 0.0, c25 = 0.0, c26 = 0.0;
        double c34 = 0.0, c35 = 0.0, c36 = 0.0;
        double c44 = 0.0, c55 = 0.0, c66 = 0.0;
        double c45 = 0.0, c46 = 0.0, c56 = 0.0;

        // 1 -> XX
        // 2 -> YY
        // 3 -> ZZ
        // 4 -> YZ
        // 5 -> XZ
        // 6 -> XY

        double x[] = new double[forceFieldEnergy.getNumberOfVariables()];
        arraycopy(xOrig, 0, x, 0, x.length);

        switch (crystal.spaceGroup.crystalSystem) {
            case CUBIC:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;
                double cavg = (c11 + c22 + c33) / 3.0;
                c11 = c22 = c33 = cavg;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c12 + c13 + c23) / 3.0;
                c12 = c13 = c23 = cavg;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                c66 = dE2dA2(6, 6, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c44 + c55 + c66) / 3.0;
                c44 = c55 = c66 = cavg;
                break;
            case HEXAGONAL:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c11 + c22) / 2.0;
                c11 = c22 = cavg;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c13 + c23) / 2.0;
                c13 = c23 = cavg;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c44 + c55) / 2.0;
                c44 = c55 = cavg;
                c66 = 0.5 * (c11 - c12);
                break;
            case TRIGONAL:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c11 + c22) / 2.0;
                c11 = c22 = cavg;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c13 + c23) / 2.0;
                c13 = c23 = cavg;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c44 + c55) / 2.0;
                c44 = c55 = cavg;
                c66 = 0.5 * (c11 - c12);

                // c14 and c24 should be equal and opposite sign.
                c14 = dE2dA2(1, 4, delta, x, eps) * PRESCON * GPA / volume;
                c24 = dE2dA2(2, 4, delta, x, eps) * PRESCON * GPA / volume;
                c56 = c14;

                String pg = crystal.getUnitCell().spaceGroup.pointGroupName;
                if (pg.equalsIgnoreCase("PG3")
                        || pg.equalsIgnoreCase("PG3bar")) {
                    // Should be equal in magnitude and opposite sign.
                    c15 = dE2dA2(1, 5, delta, x, eps) * PRESCON * GPA / volume;
                    c25 = dE2dA2(2, 5, delta, x, eps) * PRESCON * GPA / volume;
                    c46 = c25;
                }
                break;
            case TETRAGONAL:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c11 + c22) / 2.0;
                c11 = c22 = cavg;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c13 + c23) / 2.0;
                c13 = c23 = cavg;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                cavg = (c44 + c55) / 2.0;
                c44 = c55 = cavg;
                c66 = dE2dA2(6, 6, delta, x, eps) * PRESCON * GPA / volume;
                pg = crystal.getUnitCell().spaceGroup.pointGroupName;
                if (pg.equalsIgnoreCase("PG4")
                        || pg.equalsIgnoreCase("PG4bar")
                        || pg.equalsIgnoreCase("PG4/m")) {
                    // Should be equal in magnitude and opposite sign.
                    c16 = dE2dA2(1, 6, delta, x, eps) * PRESCON * GPA / volume;
                    c26 = dE2dA2(2, 6, delta, x, eps) * PRESCON * GPA / volume;
                }
                break;
            case ORTHORHOMBIC:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                c66 = dE2dA2(6, 6, delta, x, eps) * PRESCON * GPA / volume;
                break;
            case MONOCLINIC:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;

                c15 = dE2dA2(1, 5, delta, x, eps) * PRESCON * GPA / volume;
                c25 = dE2dA2(2, 5, delta, x, eps) * PRESCON * GPA / volume;
                c35 = dE2dA2(3, 5, delta, x, eps) * PRESCON * GPA / volume;
                c46 = dE2dA2(4, 6, delta, x, eps) * PRESCON * GPA / volume;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                c66 = dE2dA2(6, 6, delta, x, eps) * PRESCON * GPA / volume;
                break;
            case TRICLINIC:
            default:
                c11 = dE2dA2(1, 1, delta, x, eps) * PRESCON * GPA / volume;
                c22 = dE2dA2(2, 2, delta, x, eps) * PRESCON * GPA / volume;
                c33 = dE2dA2(3, 3, delta, x, eps) * PRESCON * GPA / volume;

                c12 = dE2dA2(1, 2, delta, x, eps) * PRESCON * GPA / volume;
                c13 = dE2dA2(1, 3, delta, x, eps) * PRESCON * GPA / volume;
                c23 = dE2dA2(2, 3, delta, x, eps) * PRESCON * GPA / volume;

                c14 = dE2dA2(1, 4, delta, x, eps) * PRESCON * GPA / volume;
                c15 = dE2dA2(1, 5, delta, x, eps) * PRESCON * GPA / volume;
                c16 = dE2dA2(1, 6, delta, x, eps) * PRESCON * GPA / volume;

                c24 = dE2dA2(2, 4, delta, x, eps) * PRESCON * GPA / volume;
                c25 = dE2dA2(2, 5, delta, x, eps) * PRESCON * GPA / volume;
                c26 = dE2dA2(2, 6, delta, x, eps) * PRESCON * GPA / volume;

                c34 = dE2dA2(3, 4, delta, x, eps) * PRESCON * GPA / volume;
                c35 = dE2dA2(3, 5, delta, x, eps) * PRESCON * GPA / volume;
                c36 = dE2dA2(3, 6, delta, x, eps) * PRESCON * GPA / volume;

                c44 = dE2dA2(4, 4, delta, x, eps) * PRESCON * GPA / volume;
                c55 = dE2dA2(5, 5, delta, x, eps) * PRESCON * GPA / volume;
                c66 = dE2dA2(6, 6, delta, x, eps) * PRESCON * GPA / volume;

                c45 = dE2dA2(4, 5, delta, x, eps) * PRESCON * GPA / volume;
                c46 = dE2dA2(4, 6, delta, x, eps) * PRESCON * GPA / volume;
                c56 = dE2dA2(5, 6, delta, x, eps) * PRESCON * GPA / volume;
        }


        if (eps > 0) {
            logger.info(format(" Elasticity Tensor using minimization and FD step size of %6.3e (GPa) =", delta));
        } else {
            logger.info(format(" Elasticity Tensor using rigid body fractional coordinates and FD step size of %6.3e (GPa) =", delta));
        }
        logger.info(format(" [ %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f ]", c11, c12, c13, c14, c15, c16));
        logger.info(format(" [ %12s %12.3f %12.3f %12.3f %12.3f %12.3f ]", "", c22, c23, c24, c25, c26));
        logger.info(format(" [ %12s %12s %12.3f %12.3f %12.3f %12.3f ]", "", "", c33, c34, c35, c36));
        logger.info(format(" [ %12s %12s %12s %12.3f %12.3f %12.3f ]", "", "", "", c44, c45, c46));
        logger.info(format(" [ %12s %12s %12s %12s %12.3f %12.3f ]", "", "", "", "", c55, c56));
        logger.info(format(" [ %12s %12s %12s %12s %12s %12.3f ]", "", "", "", "", "", c66));

        // forceFieldEnergy.energy(xOrig, verbose);
        molecularAssembly.setFractionalMode(currentFractionalMode);
        crystal.setCheckRestrictions(currentCheckRestrictions);

    }

    /**
     * Apply a strain delta.
     *
     * @param voight Voight index
     * @param delta  Strain delta
     * @param strain Strain matrix
     */
    private void applyStrain(int voight, double delta, double strain[][]) {
        switch (voight) {
            case 1: // XX
                strain[0][0] += delta;
                //strain[1][1] -= 1.0 * delta / 3.0;
                //strain[2][2] -= 1.0 * delta / 3.0;
                break;
            case 2: // YY
                //strain[0][0] -= 1.0 * delta / 3.0;
                strain[1][1] += delta;
                //strain[2][2] -= 1.0 * delta / 3.0;
                break;
            case 3: // ZZ
                //strain[0][0] -= 1.0 * delta / 3.0;
                //strain[1][1] -= 1.0 * delta / 3.0;
                strain[2][2] += delta;
                break;
            case 4: // YZ
                strain[1][2] += delta / 2.0;
                strain[2][1] += delta / 2.0;
                break;
            case 5: // XZ
                strain[0][2] += delta / 2.0;
                strain[2][0] += delta / 2.0;
                break;
            case 6: // XY
                strain[0][1] += delta / 2.0;
                strain[1][0] += delta / 2.0;
                break;
        }
    }

    private void computeStressTensor(boolean verbose) {

        double x[] = new double[forceFieldEnergy.getNumberOfVariables()];
        forceFieldEnergy.getCoordinates(x);

        forceFieldEnergy.energy(x, verbose);

        FractionalMode currentFractionalMode = molecularAssembly.getFractionalMode();

        molecularAssembly.setFractionalMode(FractionalMode.ATOM);
        molecularAssembly.computeFractionalCoordinates();

        double delta = 1.0e-5;
        double[][] dEdV = new double[3][3];

        dEdV[0][0] = dEdA(0, 0, delta, x);
        dEdV[0][1] = dEdA(0, 1, delta, x);
        dEdV[0][2] = dEdA(0, 2, delta, x);
        dEdV[1][1] = dEdA(1, 1, delta, x);
        dEdV[1][2] = dEdA(1, 2, delta, x);
        dEdV[2][2] = dEdA(2, 2, delta, x);

        logger.info("\n The Stress Tensor code is under development.");

        logger.info("\n dE/dLvec Derivatives: ");
        logger.info(format(" [ %16.8f %16.8f %16.8f ]", dEdV[0][0], dEdV[0][1], dEdV[0][2]));
        logger.info(format(" [ %16.8f %16.8f %16.8f ]", dEdV[1][0], dEdV[1][1], dEdV[1][2]));
        logger.info(format(" [ %16.8f %16.8f %16.8f ]", dEdV[2][0], dEdV[2][1], dEdV[2][2]));

        dEdV[1][0] = dEdV[0][1];
        dEdV[2][0] = dEdV[0][2];
        dEdV[2][1] = dEdV[1][2];

        // Compute a numerical virial.
        double[][] virial = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j <= i; j++) {
                virial[j][i] = 0.0;
                for (int k = 0; k < 3; k++) {
                    virial[j][i] = virial[j][i] + dEdV[k][j] * unitCell.Ai[i][k];
                }
                virial[i][j] = virial[j][i];
            }
        }

        logger.info("\n Numerical Virial: ");
        logger.info(format(" [ %16.8f %16.8f %16.8f ]", virial[0][0], virial[0][1], virial[0][2]));
        logger.info(format(" [ %16.8f %16.8f %16.8f ]", virial[1][0], virial[1][1], virial[1][2]));
        logger.info(format(" [ %16.8f %16.8f %16.8f ]", virial[2][0], virial[2][1], virial[2][2]));


        double dedv = (virial[0][0] + virial[1][1] + virial[2][2]) / (3 * unitCell.volume);

        // Conversion from kcal/mol/Ang^3 to Atm.
        double PRESCON = 6.85684112e4;
        double temperature = 298.15;
        int nAtoms = molecularAssembly.getAtomArray().length;
        double pressure = PRESCON * (nAtoms * Constants.R * temperature / unitCell.volume - dedv);

        logger.info(format("\n Pressure estimate (%6.2f Kelvin): %12.8f (atm)", temperature, pressure));

        molecularAssembly.setFractionalMode(currentFractionalMode);
    }

    private double dEdA(int ii, int jj, double delta, double x[]) {

        // Store current unit cell parameters.
        double a = unitCell.a;
        double b = unitCell.b;
        double c = unitCell.c;
        double alpha = unitCell.alpha;
        double beta = unitCell.beta;
        double gamma = unitCell.gamma;

        double cellVectors[][] = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cellVectors[i][j] = unitCell.Ai[i][j];
            }
        }

        cellVectors[ii][jj] += delta;
        try {
            if (!crystal.setCellVectors(cellVectors)) {
                return 0.0;
            }
        } catch (Exception e) {
            return 0.0;
        }

        // Update the ForceFieldEnergy with the new boundary conditions.
        forceFieldEnergy.setCrystal(crystal);

        // Update and retrieve coordinates.
        molecularAssembly.moveToFractionalCoordinates();
        forceFieldEnergy.getCoordinates(x);

        // Compute the energy.
        double eplus = forceFieldEnergy.energy(x);

        cellVectors[ii][jj] -= 2 * delta;
        try {
            if (!crystal.setCellVectors(cellVectors)) {
                return 0.0;
            }
        } catch (Exception e) {
            return 0.0;
        }

        // Update the ForceFieldEnergy with the new boundary conditions.
        forceFieldEnergy.setCrystal(crystal);

        // Update and retrieve coordinates.
        molecularAssembly.moveToFractionalCoordinates();
        forceFieldEnergy.getCoordinates(x);

        // Compute the energy.
        double eminus = forceFieldEnergy.energy(x);

        // Revert to the original boundary conditions.
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
        forceFieldEnergy.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();

        return (eplus - eminus) / (2 * delta);
    }


    private double dE2dA2(int voight1, int voight2, double delta, double x[], double eps) {

        // Store current unit cell parameters.
        double a = unitCell.a;
        double b = unitCell.b;
        double c = unitCell.c;
        double alpha = unitCell.alpha;
        double beta = unitCell.beta;
        double gamma = unitCell.gamma;

        // F(1,1)
        double[][] dStrain = new double[3][3];
        applyStrain(voight1, delta, dStrain);
        applyStrain(voight2, delta, dStrain);
        try {
            if (!crystal.perturbCellVectors(dStrain)) {
                logger.info(" Crytsal method perturbCellVectors returned false.");
                return 0.0;
            }
        } catch (Exception e) {
            logger.info(" Exception from Crytsal method perturbCellVectors.");
            return 0.0;
        }
        forceFieldEnergy.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();
        if (eps > 0.0) {
            Minimize minimize = new Minimize(molecularAssembly, forceFieldEnergy, null);
            minimize.minimize(eps);
        }
        forceFieldEnergy.getCoordinates(x);
        double e11 = forceFieldEnergy.energy(x);

        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);

        // F(-1,-1)
        dStrain = new double[3][3];
        applyStrain(voight1, -delta, dStrain);
        applyStrain(voight2, -delta, dStrain);
        try {
            if (!crystal.perturbCellVectors(dStrain)) {
                logger.info(" Crytsal method perturbCellVectors returned false.");
                return 0.0;
            }
        } catch (Exception e) {
            logger.info(" Exception from Crytsal method perturbCellVectors.");
            return 0.0;
        }
        forceFieldEnergy.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();
        if (eps > 0.0) {
            Minimize minimize = new Minimize(molecularAssembly, forceFieldEnergy, null);
            minimize.minimize(eps);
        }
        forceFieldEnergy.getCoordinates(x);
        double em1m1 = forceFieldEnergy.energy(x);

        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);

        // F(1,-1)
        dStrain = new double[3][3];
        applyStrain(voight1, delta, dStrain);
        applyStrain(voight2, -delta, dStrain);
        try {
            if (!crystal.perturbCellVectors(dStrain)) {
                logger.info(" Crytsal method perturbCellVectors returned false.");
                return 0.0;
            }
        } catch (Exception e) {
            logger.info(" Exception from Crytsal method perturbCellVectors.");
            return 0.0;
        }
        forceFieldEnergy.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();
        if (eps > 0) {
            Minimize minimize = new Minimize(molecularAssembly, forceFieldEnergy, null);
            minimize.minimize(eps);
        }
        forceFieldEnergy.getCoordinates(x);
        double e1m1 = forceFieldEnergy.energy(x);

        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);

        // F(-1,1)
        dStrain = new double[3][3];
        applyStrain(voight1, -delta, dStrain);
        applyStrain(voight2, delta, dStrain);
        try {
            if (!crystal.perturbCellVectors(dStrain)) {
                logger.info(" Crytsal method perturbCellVectors returned false.");
                return 0.0;
            }
        } catch (Exception e) {
            logger.info(" Exception from Crytsal method perturbCellVectors.");
            return 0.0;
        }
        forceFieldEnergy.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();
        if (eps > 0) {
            Minimize minimize = new Minimize(molecularAssembly, forceFieldEnergy, null);
            minimize.minimize(eps);
        }
        forceFieldEnergy.getCoordinates(x);
        double em11 = forceFieldEnergy.energy(x);

        // Change back to original unit cell parameters.
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
        forceFieldEnergy.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();
        forceFieldEnergy.getCoordinates(x);

        return (e11 - e1m1 - em11 + em1m1) / (4 * delta * delta);
    }

}
