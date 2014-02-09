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
    private final double lambdaExp = 1.0;
    private double lambdaPow = pow(lambda, lambdaExp);
    private double dLambdaPow = lambdaExp * pow(lambda, lambdaExp - 1.0);
    private double d2LambdaPow = 0;
    private final double lambdaWindow = 1.0;
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
        this.crystal = crystal.getUnitCell();
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        ForceField forceField = molecularAssembly.getForceField();

        lambdaTerm = false;
        //lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);

        if (lambdaTerm) {
            lambdaGradient = new double[nAtoms * 3];
        } else {
            lambdaGradient = null;
            this.lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
        }

        logger.info("\n Coordinate Restraint:");

        initialCoordinates = new double[3][nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            initialCoordinates[0][i] = a.getX();
            initialCoordinates[1][i] = a.getY();
            initialCoordinates[2][i] = a.getZ();
            a.print();
        }

    }

    public double residual(boolean gradient, boolean print) {

        if (lambdaTerm) {
            dEdL = 0.0;
            d2EdL2 = 0.0;
            Arrays.fill(lambdaGradient, 0.0);
        }

        double residual = 0.0;
        double fx2 = forceConstant * 2.0;
        for (int i = 0; i < nAtoms; i++) {
            // Current atomic coordinates.
            Atom atom = atoms[i];
            atom.getXYZ(a1);
            // Compute their separation.
            dx[0] = a1[0] - initialCoordinates[0][i];
            dx[1] = a1[1] - initialCoordinates[1][i];
            dx[2] = a1[2] - initialCoordinates[2][i];
            // Apply the minimum image convention.
            double r2 = crystal.image(dx);
            //logger.info(String.format(" %d %16.8f", j, Math.sqrt(r2)));
            residual += r2;
            if (gradient || lambdaTerm) {
                final double dedx = dx[0] * fx2;
                final double dedy = dx[1] * fx2;
                final double dedz = dx[2] * fx2;
                if (gradient) {
                    atom.addToXYZGradient(lambdaPow * dedx, lambdaPow * dedy, lambdaPow * dedz);
                }
                if (lambdaTerm) {
                    int j3 = i * 3;
                    lambdaGradient[j3] = dLambdaPow * dedx;
                    lambdaGradient[j3 + 1] = dLambdaPow * dedy;
                    lambdaGradient[j3 + 2] = dLambdaPow * dedz;
                }
            }
        }
        if (lambdaTerm) {
            dEdL = dLambdaPow * forceConstant * residual;
            d2EdL2 = d2LambdaPow * forceConstant * residual;
        }

        //logger.info(String.format(" Restraint Energy %16.8f", forceConstant * residual * lambdaPow));
        return forceConstant * residual * lambdaPow;
    }

    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = 1.0 - lambda;
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
                dLambdaPow = -(threeC3 * l2 + fourC4 * l3 + fiveC5 * l4) * dldgl;
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

    public void setLambdaTerm(boolean lambdaTerm) {
        this.lambdaTerm = lambdaTerm;
        setLambda(lambda);
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        if (lambdaTerm) {
            return dEdL;
        } else {
            return 0.0;
        }
    }

    @Override
    public double getd2EdL2() {
        if (lambdaTerm) {
            return d2EdL2;
        } else {
            return 0.0;
        }
    }

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
