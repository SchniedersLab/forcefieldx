/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.nonbonded;

import java.util.Arrays;
import java.util.logging.Logger;

import static java.lang.Math.pow;

import ffx.crystal.Crystal;
import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;

/**
 * Restrain atoms to their initial coordinates.
 *
 * @author Michael J. Schnieders
 */
public class CoordRestraint implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(CoordRestraint.class.getName());
    private MolecularAssembly molecularAssembly = null;
    private Crystal crystal = null;
    private final Atom atoms[];
    private final double initialCoordinates[][];
    private int nAtoms = 0;
    /**
     * Force constant in Kcal/mole/Angstrom.
     */
    private final double forceConstant = 10.0;
    private final double a1[] = new double[3];
    private final double dx[] = new double[3];
    private double lambda = 1.0;
    private final double lambdaExp = 2.0;
    private double lambdaPow = pow(lambda, lambdaExp);
    private double dLambdaPow = lambdaExp * pow(lambda, lambdaExp - 1.0);
    private double d2LambdaPow = lambdaExp * (lambdaExp - 1.0) * pow(lambda, lambdaExp - 2.0);
    private double dEdL = 0.0;
    private double d2EdL2 = 0.0;
    private final double lambdaGradient[];
    private boolean lambdaTerm = false;

    /**
     * This NCSRestraint is based on the unit cell parameters and symmetry
     * operators of the supplied crystal.
     *
     * @param molecularAssembly
     * @param crystal
     */
    public CoordRestraint(MolecularAssembly molecularAssembly, Crystal crystal) {
        this.molecularAssembly = molecularAssembly;
        this.crystal = crystal;
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        ForceField forceField = molecularAssembly.getForceField();
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);
        if (lambdaTerm) {
            lambdaGradient = new double[nAtoms * 3];
        } else {
            lambdaGradient = null;
        }

        initialCoordinates = new double[3][nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            initialCoordinates[0][i] = a.getX();
            initialCoordinates[1][i] = a.getY();
            initialCoordinates[2][i] = a.getZ();
        }

        logger.info("\n Coordinate Restraint Activated");
    }

    public double residual(boolean gradient, boolean print) {

        if (lambdaTerm) {
            dEdL = 0.0;
            d2EdL2 = 0.0;
            Arrays.fill(lambdaGradient, 0.0);
        }

        double residual = 0.0;
        double fx2 = forceConstant * 2.0;
        for (int j = 0; j < nAtoms; j++) {
            // Current atomic coordinates.
            Atom atom1 = atoms[j];
            atom1.getXYZ(a1);
            // Compute their separation.
            dx[0] = a1[0] - initialCoordinates[0][j];
            dx[1] = a1[1] - initialCoordinates[1][j];
            dx[2] = a1[2] - initialCoordinates[2][j];
            // Apply the minimum image convention.
            double r2 = crystal.image(dx);
            //logger.info(String.format(" %d %16.8f", j, Math.sqrt(r2)));
            residual += r2;
            if (gradient || lambdaTerm) {
                final double dedx = dx[0] * fx2;
                final double dedy = dx[1] * fx2;
                final double dedz = dx[2] * fx2;
                atom1.addToXYZGradient(lambdaPow * dedx, lambdaPow * dedy, lambdaPow * dedz);
                if (lambdaTerm) {
                    int j3 = j * 3;
                    lambdaGradient[j3] += dLambdaPow * dedx;
                    lambdaGradient[j3 + 1] += dLambdaPow * dedy;
                    lambdaGradient[j3 + 2] += dLambdaPow * dedz;
                }
            }
        }
        if (lambdaTerm) {
            dEdL = dLambdaPow * forceConstant * residual;
            d2EdL2 = d2LambdaPow * forceConstant * residual;
        }
        return forceConstant * residual * lambdaPow;
    }

    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = lambda;
            lambdaPow = pow(lambda, lambdaExp);
            dLambdaPow = lambdaExp * pow(lambda, lambdaExp - 1.0);
            d2LambdaPow = lambdaExp * (lambdaExp - 1.0) * pow(lambda, lambdaExp - 2.0);
        } else {
            this.lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        return dEdL;
    }

    @Override
    public double getd2EdL2() {
        return d2EdL2;
    }

    @Override
    public void getdEdXdL(double[] gradient) {
        int n3 = nAtoms * 3;
        for (int i = 0; i < n3; i++) {
            gradient[i] += lambdaGradient[i];
        }
    }
}
