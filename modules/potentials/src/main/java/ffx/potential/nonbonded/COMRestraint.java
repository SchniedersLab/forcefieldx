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
package ffx.potential.nonbonded;

import java.util.Arrays;
import java.util.logging.Logger;

import static java.lang.Math.pow;

import ffx.crystal.Crystal;
import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;

import static ffx.numerics.VectorMath.rsq;

/**
 * Restrain molecules to their center of mass.
 *
 * @author Julia Park
 */
public class COMRestraint implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(COMRestraint.class.getName());
    private final Atom atoms[];
    private int nAtoms = 0;
    /**
     * Force constant in Kcal/mole/Angstrom.
     */
    private final double forceConstant;
    private final double initialCOM[] = new double[3];
    private final double currentCOM[] = new double[3];
    private final double dx[] = new double[3];
    private final double dcomdx[];

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
    public COMRestraint(MolecularAssembly molecularAssembly, Crystal crystal) {
        ForceField forceField = molecularAssembly.getForceField();
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);

        if (lambdaTerm) {
            lambdaGradient = new double[nAtoms * 3];
        } else {
            lambdaGradient = null;
            lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
        }
        dcomdx = new double[nAtoms];
        forceConstant = forceField.getDouble(ForceField.ForceFieldDouble.RESTRAINT_K, 10.0);
        boolean gradient = false;
        computeCOM(initialCOM, gradient);
        logger.info("\n COM restraint initialized");
    }

    public double residual(boolean gradient, boolean print) {
        if (lambdaTerm) {
            dEdL = 0.0;
            d2EdL2 = 0.0;
            Arrays.fill(lambdaGradient, 0.0);
        }
        double residual = 0.0;
        double fx2 = forceConstant * 2.0;
        boolean computedcomdx = true;
        Arrays.fill(currentCOM, 0.0);
        computeCOM(currentCOM, computedcomdx);
        dx[0] = currentCOM[0] - initialCOM[0];
        dx[1] = currentCOM[1] - initialCOM[1];
        dx[2] = currentCOM[2] - initialCOM[2];
        double r2 = rsq(dx);
        residual = r2;
        for (int i = 0; i < nAtoms; i++) {
            if (gradient || lambdaTerm) {
                final double dedx = dx[0] * fx2 * dcomdx[i];
                final double dedy = dx[1] * fx2 * dcomdx[i];
                final double dedz = dx[2] * fx2 * dcomdx[i];
                // Current atomic coordinates.
                Atom atom = atoms[i];
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
        return forceConstant * residual * lambdaPow;
    }

    private void computeCOM(double[] com, boolean derivative) {
        double totalMass = 0.0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            double mass = a.getMass();
            com[0] = a.getX() * mass;
            com[1] = a.getY() * mass;
            com[2] = a.getZ() * mass;
            totalMass += mass;
        }
        com[0] /= totalMass;
        com[1] /= totalMass;
        com[2] /= totalMass;
        if (derivative) {
            for (int i = 0; i < nAtoms; i++) {
                Atom a = atoms[i];
                dcomdx[i] = a.getMass() / totalMass;
            }
        }
    }

//    private int countMolecules() {
//        int count = 0;
//        // Move polymers togethers.
//        Polymer polymers[] = molecularAssembly.getChains();
//        if (polymers != null && polymers.length > 0) {
//            count++;
//        }
//        List<Molecule> molecules = molecularAssembly.getMolecules();
//        if (molecules != null) {
//            count += molecules.size();
//        }
//        List<MSNode> waters = molecularAssembly.getWaters();
//        if (waters != null) {
//            count += waters.size();
//        }
//        List<MSNode> ions = molecularAssembly.getIons();
//        if (ions != null) {
//            count += ions.size();
//        }
//        return count;
//    }
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = lambda;
            if (this.lambda <= lambdaWindow) {
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
