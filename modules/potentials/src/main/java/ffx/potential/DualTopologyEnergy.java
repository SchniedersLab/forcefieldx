/**
 * Title: Force Field X
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
package ffx.potential;

import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.pow;
import static java.util.Arrays.fill;

import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;

/**
 * Compute the potential energy and derivatives for a dual-topology AMOEBA
 * system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class DualTopologyEnergy implements Potential, LambdaInterface {

    /**
     * Logger for the DualTopologyEnergy class.
     */
    private static final Logger logger = Logger.getLogger(DualTopologyEnergy.class.getName());
    /**
     * Topology 1 number of atoms.
     */
    private final int nAtoms1;
    /**
     * Topology 2 number of atoms.
     */
    private final int nAtoms2;
    /**
     * Topology 1 number of softcore atoms.
     */
    private final int nSoftCore1;
    /**
     * Topology 2 number of softcore atoms.
     */
    private final int nSoftCore2;
    /**
     * Shared atoms between topologies 1 and 2.
     */
    private final int nShared;
    /**
     * Total number of softcore and shared atoms: nTotal = nShared + nSoftcore1
     * + nSoftcore2
     */
    private final int nTotal;
    /**
     * Total number of variables: nVariables = nTotal * 3;
     */
    private final int nVariables;
    /**
     * Current potential energy of topology 1 (kcal/mol).
     */
    private double energy1 = 0;
    /**
     * Current potential energy of topology 2 (kcal/mol).
     */
    private double energy2 = 0;
    /**
     * Include a valence restaint energy for atoms being "disappeared."
     */
    private boolean doValenceRestraint1 = true;
    /**
     * Include a valence restaint energy for atoms being "disappeared."
     */
    private boolean doValenceRestraint2 = true;
    /**
     * Softcore restraint energy of topology 1 (kcal/mol).
     */
    private double restraintEnergy1 = 0;
    /**
     * Softcore restraint energy of topology 2 (kcal/mol).
     */
    private double restraintEnergy2 = 0;
    /**
     * Scale factor applied to the energy of Topology 1. One intended use is
     * relative crystal deposition free energy where topology 1 and topology 2
     * have different numbers of molecules in the asymmetric unit.
     */
    private double scaleEnergy1 = 1.0;
    /**
     * Scale factor applied to the energy of Topology 2.
     */
    private double scaleEnergy2 = 1.0;
    /**
     * Total energy of the dual topology, including lambda scaling.
     */
    private double totalEnergy = 0;
    /**
     * Current lambda value.
     */
    private double lambda = 1.0;
    /**
     * Current lambda value minus one.
     */
    private double oneMinusLambda = 0.0;
    /**
     * Lambda raised to the power of lambdaExponent: lambda^lambdaExponent
     */
    private double lambdaPow = 1.0;
    /**
     * One minus Lambda raised to the power of lambdaExponent:
     * (1-lambda)^lambdaExponent
     */
    private double oneMinusLambdaPow = 0.0;
    /**
     * First derivative with respect to lambda of lambda^lambdaExponent
     * lambdaExponent*lambda^(lambdaExponent-1)
     */
    private double dLambdaPow = 0.0;
    /**
     * First derivative with respect to lambda of (1-lambda)^lambdaExponent
     * -lambdaExponent*(one-lambda)^(lambdaExponent-1)
     */
    private double dOneMinusLambdaPow = 0.0;
    /**
     * Second derivative with respect to lambda of lambda^lambdaExponent
     * lambdaExponent*(lambdaExponent-1)*lambda^(lambdaExponent-2)
     */
    private double d2LambdaPow = 0.0;
    /**
     * Second derivative with respect to lambda of (1-lambda)^lambdaExponent
     * lambdaExponent*(lambdaExponent-1)*(1-lambda)^(lambdaExponent-2)
     */
    private double d2OneMinusLambdaPow = 0.0;
    /**
     * Lambda exponent that controls the thermodynamic path between topologies.
     */
    private final double lambdaExponent = 1.0;
    /**
     * Atom array for topology 1.
     */
    private final Atom[] atoms1;
    /**
     * Atom array for topology 2.
     */
    private final Atom[] atoms2;
    /**
     * Mass array for shared and softcore atoms.
     */
    private final double mass[];
    /**
     * VARIABLE_TYPE array for shared and softcore atoms.
     */
    private final VARIABLE_TYPE variableTypes[];
    /**
     * Scaling array for shared and softcore atoms.
     */
    private double scaling[] = null;
    /**
     * Topology 1 coordinates.
     */
    private final double x1[];
    /**
     * Topology 2 coordinates.
     */
    private final double x2[];
    /**
     * Topology 1 coordinate gradient.
     */
    private final double g1[];
    /**
     * Topology 2 coordinate gradient.
     */
    private final double g2[];
    /**
     * Topology 1 restraint gradient (just bonded terms of softcore atoms).
     */
    private final double rg1[];
    /**
     * Topology 2 restraint gradient (just bonded terms of softcore atoms).
     */
    private final double rg2[];
    /**
     * Topology 1 derivative of the coordinate gradient with respect to lambda.
     */
    private final double gl1[];
    /**
     * Topology 2 derivative of the coordinate gradient with respect to lambda.
     */
    private final double gl2[];
    /**
     * Topology 1 ForceFieldEnergy.
     */
    private final ForceFieldEnergy forceFieldEnergy1;
    /**
     * Topology 2 ForceFieldEnergy.
     */
    private final ForceFieldEnergy forceFieldEnergy2;

    public DualTopologyEnergy(MolecularAssembly topology1, MolecularAssembly topology2) {
        forceFieldEnergy1 = topology1.getPotentialEnergy();
        forceFieldEnergy2 = topology2.getPotentialEnergy();
        atoms1 = topology1.getAtomArray();
        atoms2 = topology2.getAtomArray();
        nAtoms1 = atoms1.length;
        nAtoms2 = atoms2.length;
        x1 = new double[nAtoms1 * 3];
        x2 = new double[nAtoms2 * 3];
        g1 = new double[nAtoms1 * 3];
        g2 = new double[nAtoms2 * 3];
        rg1 = new double[nAtoms1 * 3];
        rg2 = new double[nAtoms2 * 3];
        gl1 = new double[nAtoms1 * 3];
        gl2 = new double[nAtoms2 * 3];

        ForceField forceField1 = topology1.getForceField();
        this.doValenceRestraint1 = forceField1.getBoolean(
                ForceField.ForceFieldBoolean.LAMBDA_VALENCE_RESTRAINTS, true);
        ForceField forceField2 = topology2.getForceField();
        this.doValenceRestraint2 = forceField2.getBoolean(
                ForceField.ForceFieldBoolean.LAMBDA_VALENCE_RESTRAINTS, true);

        /**
         * Check that all atoms that are not undergoing alchemy are common to
         * both topologies.
         */
        int atomCount1 = 0;
        int atomCount2 = 0;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a1 = atoms1[i];
            if (!a1.applyLambda()) {
                atomCount1++;
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a2 = atoms2[i];
            if (!a2.applyLambda()) {
                atomCount2++;
            }
        }

        assert (atomCount1 == atomCount2);
        nShared = atomCount1;
        nSoftCore1 = nAtoms1 - nShared;
        nSoftCore2 = nAtoms2 - nShared;
        nTotal = nShared + nSoftCore1 + nSoftCore2;
        nVariables = 3 * nTotal;

        /**
         * Check that all Dual-Topology atoms start with identical coordinates.
         */
        int i1 = 0;
        int i2 = 0;
        for (int i = 0; i < nShared; i++) {
            Atom a1 = atoms1[i1++];
            while (a1.applyLambda()) {
                a1 = atoms1[i1++];
            }
            Atom a2 = atoms2[i2++];
            while (a2.applyLambda()) {
                a2 = atoms2[i2++];
            }
            assert (a1.getX() == a2.getX());
            assert (a1.getY() == a2.getY());
            assert (a1.getZ() == a2.getZ());
        }

        /**
         * All variables are coordinates.
         */
        int index = 0;
        variableTypes = new VARIABLE_TYPE[nVariables];
        for (int i = 0; i < nTotal; i++) {
            variableTypes[index++] = VARIABLE_TYPE.X;
            variableTypes[index++] = VARIABLE_TYPE.Y;
            variableTypes[index++] = VARIABLE_TYPE.Z;
        }

        /**
         * Fill the mass array.
         */
        int commonIndex = 0;
        int softcoreIndex = 3 * nShared;
        mass = new double[nVariables];
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            double m = a.getMass();
            if (!a.applyLambda()) {
                mass[commonIndex++] = m;
                mass[commonIndex++] = m;
                mass[commonIndex++] = m;
            } else {
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (a.applyLambda()) {
                double m = a.getMass();
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
            }
        }
    }

    @Override
    public double energy(double[] x) {
        /**
         * Update the coordinates of both topologies.
         */
        unpackCoordinates(x);

        /**
         * Compute the energy of topology 1.
         */
        energy1 = forceFieldEnergy1.energy(x1);
        if (doValenceRestraint1) {
            forceFieldEnergy1.setLambdaBondedTerms(true);
            restraintEnergy1 = forceFieldEnergy1.energy(x1);
            forceFieldEnergy1.setLambdaBondedTerms(false);
        } else {
            restraintEnergy1 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n",
                    scaleEnergy1 * lambdaPow * energy1, scaleEnergy1 * oneMinusLambdaPow * restraintEnergy1));
        }

        /**
         * Compute the energy of topology 2.
         */
        energy2 = forceFieldEnergy2.energy(x2);
        if (doValenceRestraint2) {
            forceFieldEnergy2.setLambdaBondedTerms(true);
            restraintEnergy2 = forceFieldEnergy2.energy(x2);
            forceFieldEnergy2.setLambdaBondedTerms(false);
        } else {
            restraintEnergy2 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 2 Energy & Restraints: %15.8f %15.8f\n",
                    scaleEnergy2 * oneMinusLambdaPow * energy2, scaleEnergy2 * lambdaPow * restraintEnergy2));
        }

        /**
         * Apply the dual-topology scaling for the total energy.
         */
        totalEnergy = scaleEnergy1 * (lambdaPow * energy1 + oneMinusLambdaPow * restraintEnergy1)
                + scaleEnergy2 * (oneMinusLambdaPow * energy2 + lambdaPow * restraintEnergy2);

        /**
         * Rescale the coordinates.
         */
        packCoordinates(x);

        return totalEnergy;
    }

    /**
     * The coordinate and gradient arrays are unpacked/packed based on the dual
     * topology.
     *
     * @param x
     * @param g
     * @return The DualTopologyEnergy total energy.
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {

        /**
         * Update the coordinates of both topologies.
         */
        unpackCoordinates(x);

        /**
         * Compute the energy and gradient of topology 1.
         */
        energy1 = forceFieldEnergy1.energyAndGradient(x1, g1);
        if (doValenceRestraint1) {
            forceFieldEnergy1.setLambdaBondedTerms(true);
            restraintEnergy1 = forceFieldEnergy1.energyAndGradient(x1, rg1);
            forceFieldEnergy1.setLambdaBondedTerms(false);
        } else {
            restraintEnergy1 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n",
                    scaleEnergy1 * lambdaPow * energy1, scaleEnergy1 * oneMinusLambdaPow * restraintEnergy1));
        }

        /**
         * Compute the energy and gradient of topology 2.
         */
        energy2 = forceFieldEnergy2.energyAndGradient(x2, g2);
        if (doValenceRestraint2) {
            forceFieldEnergy2.setLambdaBondedTerms(true);
            restraintEnergy2 = forceFieldEnergy2.energyAndGradient(x2, rg2);
            forceFieldEnergy2.setLambdaBondedTerms(false);
        } else {
            restraintEnergy2 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 2 Energy & Restraints: %15.8f %15.8f\n",
                    scaleEnergy2 * oneMinusLambdaPow * energy2, scaleEnergy2 * lambdaPow * restraintEnergy2));
        }

        /**
         * Apply the dual-topology scaling for the total energy.
         */
        totalEnergy = scaleEnergy1 * (lambdaPow * energy1 + oneMinusLambdaPow * restraintEnergy1)
                + scaleEnergy2 * (oneMinusLambdaPow * energy2 + lambdaPow * restraintEnergy2);

        /**
         * Scale and pack the gradient.
         */
        packGradient(x, g);

        return totalEnergy;
    }

    @Override
    public void setScaling(double[] scaling) {
        this.scaling = scaling;
    }

    @Override
    public double[] getScaling() {
        return scaling;
    }

    private void packCoordinates(double x[]) {
        if (scaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= scaling[i];
            }
        }
    }

    private void packGradient(double x[], double g[]) {
        if (g == null) {
            g = new double[nVariables];
        }
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        /**
         * Coordinate Gradient from Topology 1.
         */
        int index = 0;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            if (!a.applyLambda()) {
                g[indexCommon++] = scaleEnergy1 * (lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++]);
                g[indexCommon++] = scaleEnergy1 * (lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++]);
                g[indexCommon++] = scaleEnergy1 * (lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++]);
            } else {
                g[indexUnique++] = scaleEnergy1 * (lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++]);
                g[indexUnique++] = scaleEnergy1 * (lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++]);
                g[indexUnique++] = scaleEnergy1 * (lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++]);
            }
        }
        /**
         * Coordinate Gradient from Topology 2.
         */
        indexCommon = 0;
        index = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (!a.applyLambda()) {
                g[indexCommon++] += scaleEnergy2 * (oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++]);
                g[indexCommon++] += scaleEnergy2 * (oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++]);
                g[indexCommon++] += scaleEnergy2 * (oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++]);
            } else {
                g[indexUnique++] = scaleEnergy2 * (oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++]);
                g[indexUnique++] = scaleEnergy2 * (oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++]);
                g[indexUnique++] = scaleEnergy2 * (oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++]);
            }
        }

        if (scaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= scaling[i];
                g[i] /= scaling[i];
            }
        }
    }

    private void unpackCoordinates(double x[]) {

        /**
         * Unscale the coordinates.
         */
        if (scaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= scaling[i];
            }
        }

        int index = 0;
        int indexCommon = 0;
        int indexUnique = 3 * nShared;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            if (!a.applyLambda()) {
                x1[index++] = x[indexCommon++];
                x1[index++] = x[indexCommon++];
                x1[index++] = x[indexCommon++];
            } else {
                x1[index++] = x[indexUnique++];
                x1[index++] = x[indexUnique++];
                x1[index++] = x[indexUnique++];
            }
        }

        index = 0;
        indexCommon = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (!a.applyLambda()) {
                x2[index++] = x[indexCommon++];
                x2[index++] = x[indexCommon++];
                x2[index++] = x[indexCommon++];
            } else {
                x2[index++] = x[indexUnique++];
                x2[index++] = x[indexUnique++];
                x2[index++] = x[indexUnique++];
            }
        }
    }

    @Override
    public double[] getCoordinates(double[] x) {
        if (x == null) {
            x = new double[nVariables];
        }
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            if (!a.applyLambda()) {
                x[indexCommon++] = a.getX();
                x[indexCommon++] = a.getY();
                x[indexCommon++] = a.getZ();
            } else {
                x[indexUnique++] = a.getX();
                x[indexUnique++] = a.getY();
                x[indexUnique++] = a.getZ();
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (a.applyLambda()) {
                x[indexUnique++] = a.getX();
                x[indexUnique++] = a.getY();
                x[indexUnique++] = a.getZ();
            }
        }
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
        return nVariables;
    }

    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return variableTypes;
    }

    @Override
    public void setEnergyTermState(STATE state) {
        forceFieldEnergy1.setEnergyTermState(state);
        forceFieldEnergy2.setEnergyTermState(state);
    }

    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            oneMinusLambda = 1.0 - lambda;
            forceFieldEnergy1.setLambda(lambda);
            forceFieldEnergy2.setLambda(oneMinusLambda);

            lambdaPow = pow(lambda, lambdaExponent);
            dLambdaPow = lambdaExponent * pow(lambda, lambdaExponent - 1.0);
            if (lambdaExponent >= 2.0) {
                d2LambdaPow = lambdaExponent * (lambdaExponent - 1.0) * pow(lambda, lambdaExponent - 2.0);
            } else {
                d2LambdaPow = 0.0;
            }

            oneMinusLambdaPow = pow(oneMinusLambda, lambdaExponent);
            dOneMinusLambdaPow = -lambdaExponent * pow(oneMinusLambda, lambdaExponent - 1.0);
            if (lambdaExponent >= 2.0) {
                d2OneMinusLambdaPow = lambdaExponent * (lambdaExponent - 1.0) * pow(oneMinusLambda, lambdaExponent - 2.0);
            } else {
                d2OneMinusLambdaPow = 0.0;
            }
        } else {
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.severe(message);
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        double dEdL1 = scaleEnergy1 * (lambdaPow * forceFieldEnergy1.getdEdL()
                + dLambdaPow * energy1 + dOneMinusLambdaPow * restraintEnergy1);
        double dEdL2 = scaleEnergy2 * (-oneMinusLambdaPow * forceFieldEnergy2.getdEdL()
                + dOneMinusLambdaPow * energy2 + dLambdaPow * restraintEnergy2);
        return dEdL1 + dEdL2;
    }

    @Override
    public double getd2EdL2() {
        double d2EdL2_1 = scaleEnergy1 * (lambdaPow * forceFieldEnergy1.getd2EdL2()
                + 2.0 * dLambdaPow * forceFieldEnergy1.getdEdL()
                + d2LambdaPow * energy1 + d2OneMinusLambdaPow * restraintEnergy1);
        double d2EdL2_2 = scaleEnergy2 * (oneMinusLambdaPow * forceFieldEnergy2.getd2EdL2()
                - 2.0 * dOneMinusLambdaPow * forceFieldEnergy2.getdEdL()
                + d2OneMinusLambdaPow * energy2 + d2LambdaPow * restraintEnergy2);
        return d2EdL2_1 + d2EdL2_2;
    }

    @Override
    public void getdEdXdL(double[] g) {
        if (g == null) {
            g = new double[nVariables];
        }

        fill(gl1, 0.0);
        fill(gl2, 0.0);

        forceFieldEnergy1.getdEdXdL(gl1);
        forceFieldEnergy2.getdEdXdL(gl2);

        int index = 0;
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        /**
         * Coordinate Gradient from Topology 1.
         */
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            if (!a.applyLambda()) {
                g[indexCommon++] = scaleEnergy1 * (lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + dOneMinusLambdaPow * rg1[index++]);
                g[indexCommon++] = scaleEnergy1 * (lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + dOneMinusLambdaPow * rg1[index++]);
                g[indexCommon++] = scaleEnergy1 * (lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + dOneMinusLambdaPow * rg1[index++]);
            } else {
                g[indexUnique++] = scaleEnergy1 * (lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + dOneMinusLambdaPow * rg1[index++]);
                g[indexUnique++] = scaleEnergy1 * (lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + dOneMinusLambdaPow * rg1[index++]);
                g[indexUnique++] = scaleEnergy1 * (lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + dOneMinusLambdaPow * rg1[index++]);
            }
        }

        /**
         * Coordinate Gradient from Topology 2.
         */
        index = 0;
        indexCommon = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (!a.applyLambda()) {
                g[indexCommon++] += scaleEnergy2 * (-oneMinusLambdaPow * gl2[index]
                        + dOneMinusLambdaPow * g2[index] + dLambdaPow * rg2[index++]);
                g[indexCommon++] += scaleEnergy2 * (-oneMinusLambdaPow * gl2[index]
                        + dOneMinusLambdaPow * g2[index] + dLambdaPow * rg2[index++]);
                g[indexCommon++] += scaleEnergy2 * (-oneMinusLambdaPow * gl2[index]
                        + dOneMinusLambdaPow * g2[index] + dLambdaPow * rg2[index++]);
            } else {
                g[indexUnique++] = scaleEnergy2 * (-oneMinusLambdaPow * gl2[index]
                        + dOneMinusLambdaPow * g2[index] + dLambdaPow * rg2[index++]);
                g[indexUnique++] = scaleEnergy2 * (-oneMinusLambdaPow * gl2[index]
                        + dOneMinusLambdaPow * g2[index] + dLambdaPow * rg2[index++]);
                g[indexUnique++] = scaleEnergy2 * (-oneMinusLambdaPow * gl2[index]
                        + dOneMinusLambdaPow * g2[index] + dLambdaPow * rg2[index++]);
            }
        }
    }

    /**
     * Scale the energy of Topology 1.
     *
     * @param scaleEnergy1 the scale factor to apply to the energy of topology
     * 1.
     */
    public void setScaleEnergy1(double scaleEnergy1) {
        this.scaleEnergy1 = scaleEnergy1;
    }

    /**
     * Scale the energy of Topology 2.
     *
     * @param scaleEnergy2 the scale factor to apply to the energy of topology
     * 2.
     */
    public void setScaleEnergy2(double scaleEnergy2) {
        this.scaleEnergy2 = scaleEnergy2;
    }
}
