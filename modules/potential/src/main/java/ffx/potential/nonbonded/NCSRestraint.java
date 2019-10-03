//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.nonbonded;

import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.pow;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.ForceField;
import static ffx.numerics.math.VectorMath.rsq;

/**
 * Given unit cell parameters and symmetry operators, NCS copies are restrained
 * to the asymmetric unit atoms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class NCSRestraint implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(NCSRestraint.class.getName());
    private Crystal ncsCrystal;
    private SpaceGroup spaceGroup;
    private final Atom[] atoms;
    private int nAtoms;
    private final double forceConstant = 10.0;
    private final double[][] transOp = new double[3][3];
    private int nSymm;
    private final double[] a1 = new double[3];
    private final double[] a2 = new double[3];
    private final double[] dx = new double[3];
    private double lambda = 1.0;
    private final double lambdaExp = 2.0;
    private double lambdaPow = pow(lambda, lambdaExp);
    private double dLambdaPow = lambdaExp * pow(lambda, lambdaExp - 1.0);
    private double d2LambdaPow = lambdaExp * (lambdaExp - 1.0) * pow(lambda, lambdaExp - 2.0);
    private double dEdL = 0.0;
    private double d2EdL2 = 0.0;
    private final double[] lambdaGradient;
    private boolean lambdaTerm;

    /**
     * This NCSRestraint is based on the unit cell parameters and symmetry
     * operators of the supplied crystal.
     *
     * @param atoms      the Atom array to construct this NCSRestraint from.
     * @param forceField the ForceField parameters
     * @param crystal    the Crystal specifies symmetry and PBCs.
     */
    public NCSRestraint(Atom[] atoms, ForceField forceField, Crystal crystal) {
        this.ncsCrystal = crystal;
        this.atoms = atoms;
        nAtoms = atoms.length;
        spaceGroup = this.ncsCrystal.spaceGroup;
        nSymm = spaceGroup.getNumberOfSymOps();
        assert (nAtoms % nSymm == 0);
        lambdaTerm = forceField.getBoolean("LAMBDATERM", false);
        if (lambdaTerm) {
            lambdaGradient = new double[nAtoms * 3];
        } else {
            lambdaGradient = null;
            this.lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
        }

        logger.info(format("\n NCS Restraint%s", crystal.toString()));

    }

    /**
     * <p>residual.</p>
     *
     * @param gradient a boolean.
     * @param print    a boolean.
     * @return a double.
     */
    public double residual(boolean gradient, boolean print) {

        // Check that the number of atom is multiple of the number of symmetry operators.
        if (nAtoms % nSymm != 0) {
            return 0;
        }

        if (lambdaTerm) {
            dEdL = 0.0;
            d2EdL2 = 0.0;
            fill(lambdaGradient, 0.0);
        }

        double residual = 0.0;
        int nAsymmAtoms = nAtoms / nSymm;
        double fx2 = forceConstant * 2.0;
        for (int i = 1; i < nSymm; i++) {
            SymOp symOp = spaceGroup.getSymOp(i);
            ncsCrystal.getTransformationOperator(symOp, transOp);
            int offset = nAsymmAtoms * i;
            for (int j = 0; j < nAsymmAtoms; j++) {
                // Reference atom of the asymmetric unit
                Atom atom1 = atoms[j];
                atom1.getXYZ(a1);

                if (atom1.isHydrogen()) {
                    continue;
                }

                // Apply the symmetry operator to superpose the reference atom with its symmetry mate.
                ncsCrystal.applySymOp(a1, a1, symOp);

                // Symmetry mate atom.
                Atom atom2 = atoms[offset + j];
                atom2.getXYZ(a2);
                // Compute their separation.
                dx[0] = a1[0] - a2[0];
                dx[1] = a1[1] - a2[1];
                dx[2] = a1[2] - a2[2];
                // Apply the minimum image convention.
                double r2 = rsq(dx);
                // double r2 = ncsCrystal.image(dx);
                //logger.info(String.format(" %d %16.8f", j, Math.sqrt(r2)));
                residual += r2;
                if (gradient || lambdaTerm) {
                    final double dedx = dx[0] * fx2;
                    final double dedy = dx[1] * fx2;
                    final double dedz = dx[2] * fx2;
                    atom2.addToXYZGradient(-lambdaPow * dedx, -lambdaPow * dedy, -lambdaPow * dedz);

                    // Rotate the force on the reference atom back into the asymmetric unit.
                    final double dedx1 = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
                    final double dedy1 = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
                    final double dedz1 = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];
                    atom1.addToXYZGradient(lambdaPow * dedx1, lambdaPow * dedy1, lambdaPow * dedz1);
                    if (lambdaTerm) {
                        int j3 = j * 3;
                        lambdaGradient[j3] += dLambdaPow * dedx1;
                        lambdaGradient[j3 + 1] += dLambdaPow * dedy1;
                        lambdaGradient[j3 + 2] += dLambdaPow * dedz1;
                        int oj3 = (offset + j) * 3;
                        lambdaGradient[oj3] -= dLambdaPow * dedx;
                        lambdaGradient[oj3 + 1] -= dLambdaPow * dedy;
                        lambdaGradient[oj3 + 2] -= dLambdaPow * dedz;
                    }
                }
            }
        }
        if (lambdaTerm) {
            dEdL = dLambdaPow * forceConstant * residual;
            d2EdL2 = d2LambdaPow * forceConstant * residual;
        }
        return forceConstant * residual * lambdaPow;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = lambda;

            double lambdaWindow = 1.0;

            if (this.lambda < lambdaWindow) {
                double dldgl = 1.0 / lambdaWindow;
                double l = dldgl * this.lambda;
                double l2 = l * l;
                double l3 = l2 * l;
                double l4 = l2 * l2;
                double l5 = l4 * l;
                double c3 = 10.0;
                double c4 = -15.0;
                double c5 = 6.0;
                double threeC3 = 3.0 * c3;
                double sixC3 = 6.0 * c3;
                double fourC4 = 4.0 * c4;
                double twelveC4 = 12.0 * c4;
                double fiveC5 = 5.0 * c5;
                double twentyC5 = 20.0 * c5;
                lambdaPow = c3 * l3 + c4 * l4 + c5 * l5;
                dLambdaPow = (threeC3 * l2 + fourC4 * l3 + fiveC5 * l4) * dldgl;
                d2LambdaPow = (sixC3 * l + twelveC4 * l2 + twentyC5 * l3) * dldgl * dldgl;
            } else {
                lambdaPow = 1.0;
                dLambdaPow = 0.0;
                d2LambdaPow = 0.0;
            }
        } else {
            this.lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
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
        if (lambdaTerm) {
            return dEdL;
        } else {
            return 0.0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        if (lambdaTerm) {
            return d2EdL2;
        } else {
            return 0.0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        if (lambdaTerm) {
            int n3 = nAtoms * 3;
            for (int i = 0; i < n3; i++) {
                gradient[i] += lambdaGradient[i];
            }
        }
    }
}
