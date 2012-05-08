/**
 * Title: Force Field X Description: Force Field X - Software for Molecular
 * Biophysics. Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import java.util.logging.Logger;

import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;

/**
 * Compute the potential energy and derivatives for a dual-topology AMOEBA
 * system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class DualTopology implements Potential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(DualTopology.class.getName());
    private final ForceFieldEnergy forceFieldEnergy1;
    private final ForceFieldEnergy forceFieldEnergy2;
    private double energy1 = 0;
    private double energy2 = 0;
    private double totalEnergy = 0;
    private final Atom[] atoms1;
    private final Atom[] atoms2;
    private final int nAtoms1, nAtoms2;
    private final int nUnique1, nUnique2, nCommon, nTotal, nVariables;
    private final VARIABLE_TYPE variableTypes[];
    private final double mass[];
    private final double x1[], x2[], g1[], g2[], gl1[], gl2[];
    private double scaling[] = null;
    private double lambda = 1.0;
    private double oneMinusLambda = 0.0;

    public DualTopology(MolecularAssembly topology1, MolecularAssembly topology2) {
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
        gl1 = new double[nAtoms1 * 3];
        gl2 = new double[nAtoms2 * 3];

        /**
         * Check that all atoms that are not undergoing alchemy are common to
         * both topologies.
         */
        int atomCount1 = 0;
        int atomCount2 = 0;
        for (int i = 0; i < nAtoms1; i++) {
            if (!atoms1[i].applyLambda()) {
                atomCount1++;
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            if (!atoms2[i].applyLambda()) {
                atomCount2++;
            }
        }
        assert (atomCount1 == atomCount2);
        nCommon = atomCount1;
        nUnique1 = nAtoms1 - nCommon;
        nUnique2 = nAtoms2 - nCommon;
        nTotal = nCommon + nUnique1 + nUnique2;
        nVariables = 3 * nTotal;

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
        int uniqueIndex = 3 * nCommon;
        mass = new double[nVariables];
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            double m = a.getMass();
            if (!a.applyLambda()) {
                mass[commonIndex++] = m;
                mass[commonIndex++] = m;
                mass[commonIndex++] = m;
            } else {
                mass[uniqueIndex++] = m;
                mass[uniqueIndex++] = m;
                mass[uniqueIndex++] = m;
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (a.applyLambda()) {
                double m = a.getMass();
                mass[uniqueIndex++] = m;
                mass[uniqueIndex++] = m;
                mass[uniqueIndex++] = m;
            }
        }
    }

    /**
     * The coordinates and gradient arrays are unpacked/packed based on the dual
     * topology.
     *
     * @param x
     * @param g
     * @return
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

        /**
         * Compute the energy and gradient of topology 2.
         */
        energy2 = forceFieldEnergy2.energyAndGradient(x2, g2);

        /**
         * Scale and pack the gradient.
         */
        packGradient(g);

        /**
         * Apply the dual-topology scaling for the total energy.
         */
        totalEnergy = lambda * energy1 + oneMinusLambda * energy2;

        logger.fine(String.format(" E1: %15.8f E2: %15.8f Etot: %15.8f", energy1, energy2, totalEnergy));

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

    private void packGradient(double g[]) {
        if (g == null) {
            g = new double[nVariables];
        }
        int indexCommon = 0;
        int indexUnique = nCommon * 3;
        double grad[] = new double[3];
        /**
         * Coordinate Gradient from Topology 1.
         */
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            a.getXYZGradient(grad);
            if (!a.applyLambda()) {
                g[indexCommon++] = lambda * grad[0];
                g[indexCommon++] = lambda * grad[1];
                g[indexCommon++] = lambda * grad[2];
            } else {
                g[indexUnique++] = lambda * grad[0];
                g[indexUnique++] = lambda * grad[1];
                g[indexUnique++] = lambda * grad[2];
            }
        }
        /**
         * Coordinate Gradient from Topology 2.
         */
        indexCommon = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            a.getXYZGradient(grad);
            if (!a.applyLambda()) {
                g[indexCommon++] += oneMinusLambda * grad[0];
                g[indexCommon++] += oneMinusLambda * grad[1];
                g[indexCommon++] += oneMinusLambda * grad[2];
            } else {
                g[indexUnique++] = oneMinusLambda * grad[0];
                g[indexUnique++] = oneMinusLambda * grad[1];
                g[indexUnique++] = oneMinusLambda * grad[2];
            }
        }
    }

    private void unpackCoordinates(double x[]) {
        int index = 0;
        int indexCommon = 0;
        int indexUnique = 3 * nCommon;
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
        int indexUnique = nCommon * 3;
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
                a.print();
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            if (a.applyLambda()) {
                x[indexUnique++] = a.getX();
                x[indexUnique++] = a.getY();
                x[indexUnique++] = a.getZ();
                a.print();
            }
        }
        return x;
    }

    @Override
    public double[] getMass() {
        return mass;
    }

    @Override
    public double getTotal() {
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
        double dEdL1 = lambda * forceFieldEnergy1.getdEdL() + energy1;
        double dEdL2 = oneMinusLambda * forceFieldEnergy2.getdEdL() - energy2;
        return dEdL1 + dEdL2;
    }

    @Override
    public double getd2EdL2() {
        double d2EdL2_1 = lambda * forceFieldEnergy1.getd2EdL2() + forceFieldEnergy1.getdEdL();
        double d2EdL2_2 = oneMinusLambda * forceFieldEnergy2.getd2EdL2() - forceFieldEnergy2.getdEdL();
        return d2EdL2_1 + d2EdL2_2;
    }

    @Override
    public void getdEdXdL(double[] g) {
        if (g == null) {
            g = new double[nVariables];
        }

        forceFieldEnergy1.getdEdXdL(gl1);
        forceFieldEnergy2.getdEdXdL(gl2);

        int index = 0;
        int indexCommon = 0;
        int indexUnique = nCommon * 3;
        double grad[] = new double[3];
        /**
         * Coordinate Gradient from Topology 1.
         */
        for (int i = 0; i < nAtoms1; i++) {
            Atom a = atoms1[i];
            a.getXYZGradient(grad);
            if (!a.applyLambda()) {
                g[indexCommon++] = lambda * gl1[index++] + grad[0];
                g[indexCommon++] = lambda * gl1[index++] + grad[1];
                g[indexCommon++] = lambda * gl1[index++] + grad[2];
            } else {
                g[indexUnique++] = lambda * gl1[index++] + grad[0];
                g[indexUnique++] = lambda * gl1[index++] + grad[1];
                g[indexUnique++] = lambda * gl1[index++] + grad[2];
            }
        }

        /**
         * Coordinate Gradient from Topology 2.
         */
        index = 0;
        indexCommon = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a = atoms2[i];
            a.getXYZGradient(grad);
            if (!a.applyLambda()) {
                g[indexCommon++] += (oneMinusLambda * gl2[index++] - grad[0]);
                g[indexCommon++] += (oneMinusLambda * gl2[index++] - grad[1]);
                g[indexCommon++] += (oneMinusLambda * gl2[index++] - grad[2]);
            } else {
                g[indexUnique++] = oneMinusLambda * gl2[index++] - grad[0];
                g[indexUnique++] = oneMinusLambda * gl2[index++] - grad[1];
                g[indexUnique++] = oneMinusLambda * gl2[index++] - grad[2];
            }
        }
    }
}
