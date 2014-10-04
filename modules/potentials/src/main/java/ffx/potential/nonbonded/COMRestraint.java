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

import java.util.List;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.pow;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
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
    private final int nAtoms;
    private MolecularAssembly molecularAssembly = null;
    private Crystal crystal = null;
    private final int nMolecules;
    /**
     * Force constant in Kcal/mole/Angstrom.
     */
    private final double forceConstant;
    private final double initialCOM[][];
    private double currentCOM[][];
    private double dx[] = new double[3];
    private double dcomdx[];

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
        this.molecularAssembly = molecularAssembly;
        this.crystal = crystal.getUnitCell();
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        nMolecules = countMolecules();
        initialCOM = new double[3][nMolecules];
        currentCOM = new double[3][nMolecules];
        ForceField forceField = molecularAssembly.getForceField();

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

        computeCOM(initialCOM, nMolecules);

        logger.info("\n COM restraint initialized");
    }

    public double residual(boolean gradient, boolean print) {
        if (lambdaTerm) {
            dEdL = 0.0;
            d2EdL2 = 0.0;
            fill(lambdaGradient, 0.0);
        }
        double residual = 0.0;
        double fx2 = forceConstant * 2.0;
        //boolean computedcomdx = true;
        //fill(currentCOM, 0.0);
        computeCOM(currentCOM, nMolecules);
        computedcomdx();
        for (int i = 0; i < nMolecules; i++) {
            dx[0] = currentCOM[0][i] - initialCOM[0][i];
            dx[1] = currentCOM[1][i] - initialCOM[1][i];
            dx[2] = currentCOM[2][i] - initialCOM[2][i];

            double r2 = rsq(dx);
            residual += r2;
            for (int j = 0; j < nAtoms; j++) {
                if (gradient || lambdaTerm) {
                    final double dedx = dx[0] * fx2 * dcomdx[j];
                    final double dedy = dx[1] * fx2 * dcomdx[j];
                    final double dedz = dx[2] * fx2 * dcomdx[j];
                    // Current atomic coordinates.
                    Atom atom = atoms[j];
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
        }
        if (lambdaTerm) {
            dEdL = dLambdaPow * forceConstant * residual;
            d2EdL2 = d2LambdaPow * forceConstant * residual;
        }
        return forceConstant * residual * lambdaPow;
    }

    private void computeCOM(double[][] com, int nMolecules) {
        int i = 0;
        while (i < nMolecules) {
            Polymer polymers[] = molecularAssembly.getChains();
            if (polymers != null && polymers.length > 0) {
                // Find the center of mass
                for (Polymer polymer : polymers) {
                    List<Atom> list = polymer.getAtomList();
                    com[0][i] = 0.0;
                    com[1][i] = 0.0;
                    com[2][i] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0][i] += atom.getX() * m;
                        com[1][i] += atom.getY() * m;
                        com[2][i] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0][i] /= totalMass;
                    com[1][i] /= totalMass;
                    com[2][i] /= totalMass;
                    i++;
                }
            }

            // Loop over each molecule
            List<Molecule> molecules = molecularAssembly.getMolecules();
            for (MSNode molecule : molecules) {
                List<Atom> list = molecule.getAtomList();
                // Find the center of mass
                com[0][i] = 0.0;
                com[1][i] = 0.0;
                com[2][i] = 0.0;
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    com[0][i] += atom.getX() * m;
                    com[1][i] += atom.getY() * m;
                    com[2][i] += atom.getZ() * m;
                    totalMass += m;
                }
                com[0][i] /= totalMass;
                com[1][i] /= totalMass;
                com[2][i] /= totalMass;
                i++;
            }

            // Loop over each water
            List<MSNode> waters = molecularAssembly.getWaters();
            for (MSNode water : waters) {
                List<Atom> list = water.getAtomList();
                // Find the center of mass
                com[0][i] = 0.0;
                com[1][i] = 0.0;
                com[2][i] = 0.0;
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    com[0][i] += atom.getX() * m;
                    com[1][i] += atom.getY() * m;
                    com[2][i] += atom.getZ() * m;
                    totalMass += m;
                }
                com[0][i] /= totalMass;
                com[1][i] /= totalMass;
                com[2][i] /= totalMass;
                i++;
            }

            // Loop over each ion
            List<MSNode> ions = molecularAssembly.getIons();
            for (MSNode ion : ions) {
                List<Atom> list = ion.getAtomList();
                // Find the center of mass
                com[0][i] = 0.0;
                com[1][i] = 0.0;
                com[2][i] = 0.0;
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    com[0][i] += atom.getX() * m;
                    com[1][i] += atom.getY() * m;
                    com[2][i] += atom.getZ() * m;
                    totalMass += m;
                }
                com[0][i] /= totalMass;
                com[1][i] /= totalMass;
                com[2][i] /= totalMass;
                i++;
            }
        }
    }

    private void computedcomdx() {
//        double totalMass = 0.0;
        int i = 0;
        while (i < nAtoms) {
            Polymer polymers[] = molecularAssembly.getChains();
            if (polymers != null && polymers.length > 0) {
                for (Polymer polymer : polymers) {
                    List<Atom> list = polymer.getAtomList();
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        totalMass += m;
                    }
                    for (Atom atom : list) {
                        dcomdx[i] = atom.getMass();
                        dcomdx[i] /= totalMass;
                        i++;
                    }
                }
            }

            // Loop over each molecule
            List<Molecule> molecules = molecularAssembly.getMolecules();
            for (MSNode molecule : molecules) {
                List<Atom> list = molecule.getAtomList();
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    totalMass += m;
                }
                for (Atom atom : list) {
                    dcomdx[i] = atom.getMass();
                    dcomdx[i] /= totalMass;
                    i++;
                }
            }

            // Loop over each water
            List<MSNode> waters = molecularAssembly.getWaters();
            for (MSNode water : waters) {
                List<Atom> list = water.getAtomList();
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    totalMass += m;
                }
                for (Atom atom : list) {
                    dcomdx[i] = atom.getMass();
                    dcomdx[i] /= totalMass;
                    i++;
                }
            }

            // Loop over each ion
            List<MSNode> ions = molecularAssembly.getIons();
            for (MSNode ion : ions) {
                List<Atom> list = ion.getAtomList();
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    totalMass += m;
                }
                for (Atom atom : list) {
                    dcomdx[i] = atom.getMass();
                    dcomdx[i] /= totalMass;
                    i++;
                }
            }
        }
//        for (int i = 0; i < nAtoms; i++) {
//            Atom a = atoms[i];
//            dcomdx[j] = a.getMass() / totalMass;
//        }

    }

    private int countMolecules() {
        int count = 0;
        // Move polymers togethers.
        Polymer polymers[] = molecularAssembly.getChains();
        if (polymers != null && polymers.length > 0) {
            count += polymers.length;
        }
        List<Molecule> molecules = molecularAssembly.getMolecules();
        if (molecules != null) {
            count += molecules.size();
        }
        List<MSNode> waters = molecularAssembly.getWaters();
        if (waters != null) {
            count += waters.size();
        }
        List<MSNode> ions = molecularAssembly.getIons();
        if (ions != null) {
            count += ions.size();
        }
        return count;
    }

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
