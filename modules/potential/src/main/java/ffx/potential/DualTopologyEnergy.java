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

/*import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;*/
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.pow;

import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.ForceField;

/**
 * Compute the potential energy and derivatives for a dual-topology AMOEBA
 * system.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
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
     * End-state restraint energy of topology 1 (kcal/mol).
     */
    private double restraintEnergy1 = 0;
    /**
     * End-state restraint energy of topology 2 (kcal/mol).
     */
    private double restraintEnergy2 = 0;
    /**
     * dEdL of topology 1 (kcal/mol).
     */
    private double dEdL_1 = 0;
    /**
     * dEdL of topology 2 (kcal/mol).
     */
    private double dEdL_2 = 0;
    /**
     * End-state restraint dEdL of topology 1 (kcal/mol).
     */
    private double restraintdEdL_1 = 0;
    /**
     * End-state restraint dEdL of topology 2 (kcal/mol).
     */
    private double restraintdEdL_2 = 0;
    /**
     * d2EdL2 of topology 1 (kcal/mol).
     */
    private double d2EdL2_1 = 0;
    /**
     * d2EdL2 of topology 2 (kcal/mol).
     */
    private double d2EdL2_2 = 0;
    /**
     * End-state restraint d2EdL2 of topology 1 (kcal/mol).
     */
    private double restraintd2EdL2_1 = 0;
    /**
     * End-state restraint d2EdL2 of topology 2 (kcal/mol).
     */
    private double restraintd2EdL2_2 = 0;
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
     * Topology 1 restraint gradient for end state bonded terms.
     */
    private final double rg1[];
    /**
     * Topology 2 restraint gradient for end state bonded terms.
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
     * Topology 1 derivative of the coordinate gradient with respect to lambda
     * for end state bonded terms
     */
    private final double rgl1[];
    /**
     * Topology 2 derivative of the coordinate gradient with respect to lambda
     * for end state bonded terms
     */
    private final double rgl2[];
    /**
     * Topology 1 Potential.
     */
    private final Potential potential1;
    /**
     * Topology 2 Potential.
     */
    private final Potential potential2;
    /**
     * Topology 1 LambdaInterface.
     */
    private final LambdaInterface lambdaInterface1;
    /**
     * Topology 2 LambdaInterface.
     */
    private final LambdaInterface lambdaInterface2;
    /**
     * Topology 1 ForceFieldEnergy.
     */
    private final ForceFieldEnergy forceFieldEnergy1;
    /**
     * Topology 2 ForceFieldEnergy.
     */
    private final ForceFieldEnergy forceFieldEnergy2;

    /*private final ParallelTeam dualTopologyTeam;
    private final EnergyRegion energyRegion;*/

    private final int nActive1;
    private final int nActive2;
    private final Atom activeAtoms1[];
    private final Atom activeAtoms2[];

    public DualTopologyEnergy(Potential topology1, Atom atoms1[],
            Potential topology2, Atom atoms2[]) {
        potential1 = topology1;
        potential2 = topology2;
        lambdaInterface1 = (LambdaInterface) potential1;
        lambdaInterface2 = (LambdaInterface) potential2;
        this.atoms1 = atoms1;
        this.atoms2 = atoms2;
        nAtoms1 = atoms1.length;
        nAtoms2 = atoms2.length;
        forceFieldEnergy1 = null;
        forceFieldEnergy2 = null;
        doValenceRestraint1 = false;
        doValenceRestraint2 = false;

        /**
         * Check that all atoms that are not undergoing alchemy are common to
         * both topologies.
         */
        int shared1 = 0;
        int shared2 = 0;
        int activeCount1 = 0;
        int activeCount2 = 0;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a1 = atoms1[i];
            if (a1.isActive()) {
                activeCount1++;
                if (!a1.applyLambda()) {
                    shared1++;
                }
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a2 = atoms2[i];
            if (a2.isActive()) {
                activeCount2++;
                if (!a2.applyLambda()) {
                    shared2++;
                }
            }
        }
        nActive1 = activeCount1;
        nActive2 = activeCount2;
        activeAtoms1 = new Atom[activeCount1];
        activeAtoms2 = new Atom[activeCount2];
        int index = 0;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a1 = atoms1[i];
            if (a1.isActive()) {
                activeAtoms1[index++] = a1;
            }
        }
        index = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a2 = atoms2[i];
            if (a2.isActive()) {
                activeAtoms2[index++] = a2;
            }
        }

        assert (shared1 == shared2);
        nShared = shared1;
        nSoftCore1 = nActive1 - nShared;
        nSoftCore2 = nActive2 - nShared;
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
         * Allocate memory for coordinates and derivatives.
         */
        x1 = new double[nActive1 * 3];
        x2 = new double[nActive2 * 3];
        g1 = new double[nActive1 * 3];
        g2 = new double[nActive2 * 3];
        rg1 = new double[nActive1 * 3];
        rg2 = new double[nActive2 * 3];
        gl1 = new double[nActive1 * 3];
        gl2 = new double[nActive2 * 3];
        rgl1 = new double[nActive1 * 3];
        rgl2 = new double[nActive2 * 3];

        /**
         * All variables are coordinates.
         */
        index = 0;
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
        for (int i = 0; i < nActive1; i++) {
            Atom a = activeAtoms1[i];
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
        for (int i = 0; i < nActive2; i++) {
            Atom a = activeAtoms2[i];
            if (a.applyLambda()) {
                double m = a.getMass();
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
            }
        }
        /*dualTopologyTeam = new ParallelTeam(2);
        // Will probably not be an issue, even with only one thread allocated.
        energyRegion = new EnergyRegion();*/
    }

    public DualTopologyEnergy(MolecularAssembly topology1, MolecularAssembly topology2) {
        forceFieldEnergy1 = topology1.getPotentialEnergy();
        forceFieldEnergy2 = topology2.getPotentialEnergy();
        potential1 = forceFieldEnergy1;
        potential2 = forceFieldEnergy2;
        lambdaInterface1 = forceFieldEnergy1;
        lambdaInterface2 = forceFieldEnergy2;
        atoms1 = topology1.getAtomArray();
        atoms2 = topology2.getAtomArray();
        nAtoms1 = atoms1.length;
        nAtoms2 = atoms2.length;

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
        int shared1 = 0;
        int shared2 = 0;
        int activeCount1 = 0;
        int activeCount2 = 0;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a1 = atoms1[i];
            if (a1.isActive()) {
                activeCount1++;
                if (!a1.applyLambda()) {
                    shared1++;
                }
            }
        }
        for (int i = 0; i < nAtoms2; i++) {
            Atom a2 = atoms2[i];
            if (a2.isActive()) {
                activeCount2++;
                if (!a2.applyLambda()) {
                    shared2++;
                }
            }
        }
        nActive1 = activeCount1;
        nActive2 = activeCount2;
        activeAtoms1 = new Atom[activeCount1];
        activeAtoms2 = new Atom[activeCount2];
        int index = 0;
        for (int i = 0; i < nAtoms1; i++) {
            Atom a1 = atoms1[i];
            if (a1.isActive()) {
                activeAtoms1[index++] = a1;
            }
        }
        index = 0;
        for (int i = 0; i < nAtoms2; i++) {
            Atom a2 = atoms2[i];
            if (a2.isActive()) {
                activeAtoms2[index++] = a2;
            }
        }

        assert (shared1 == shared2);
        nShared = shared1;
        nSoftCore1 = nActive1 - nShared;
        nSoftCore2 = nActive2 - nShared;
        nTotal = nShared + nSoftCore1 + nSoftCore2;
        nVariables = 3 * nTotal;

        /**
         * Allocate memory for coordinates and derivatives.
         */
        x1 = new double[nActive1 * 3];
        x2 = new double[nActive2 * 3];
        g1 = new double[nActive1 * 3];
        g2 = new double[nActive2 * 3];
        rg1 = new double[nActive1 * 3];
        rg2 = new double[nActive2 * 3];
        gl1 = new double[nActive1 * 3];
        gl2 = new double[nActive2 * 3];
        rgl1 = new double[nActive1 * 3];
        rgl2 = new double[nActive2 * 3];

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
        index = 0;
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
        for (int i = 0; i < nActive1; i++) {
            Atom a = activeAtoms1[i];
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
        for (int i = 0; i < nActive2; i++) {
            Atom a = activeAtoms2[i];
            if (a.applyLambda()) {
                double m = a.getMass();
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
                mass[softcoreIndex++] = m;
            }
        }
        /*dualTopologyTeam = new ParallelTeam(2);
        // Will probably not be an issue, even with only one thread allocated.
        energyRegion = new EnergyRegion();*/
    }

    @Override
    public double energy(double[] x) {
        /**
         * Update the coordinates of both topologies.
         */
        unpackCoordinates(x);

        /*
         try {
         dualTopologyTeam.execute(energyRegion);
         } catch (Exception e) {
         String message = "Fatal exception computing the Dual-Topology Energy.\n";
         logger.log(Level.SEVERE, message, e);
         } */
        /**
         * Compute the energy of topology 1.
         */
        energy1 = potential1.energy(x1);
        if (doValenceRestraint1 && potential1 instanceof ForceFieldEnergy) {
            ForceFieldEnergy ffE1 = (ForceFieldEnergy) potential1;
            ffE1.setLambdaBondedTerms(true);
            restraintEnergy1 = potential1.energy(x1);
            ffE1.setLambdaBondedTerms(false);
        } else {
            restraintEnergy1 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n",
                    lambdaPow * energy1, oneMinusLambdaPow * restraintEnergy1));
        }
        /**
         * Compute the energy of topology 2.
         */

        energy2 = potential2.energy(x2);
        if (doValenceRestraint2 && potential2 instanceof ForceFieldEnergy) {
            ForceFieldEnergy ffE2 = (ForceFieldEnergy) potential2;
            ffE2.setLambdaBondedTerms(true);
            restraintEnergy2 = potential2.energy(x2);
            ffE2.setLambdaBondedTerms(false);
        } else {
            restraintEnergy2 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 2 Energy & Restraints: %15.8f %15.8f\n",
                    oneMinusLambdaPow * energy2, lambdaPow * restraintEnergy2));
        }
        /**
         * Apply the dual-topology scaling for the total energy.
         */
        totalEnergy = lambdaPow * energy1 + oneMinusLambdaPow * restraintEnergy1
                + oneMinusLambdaPow * energy2 + lambdaPow * restraintEnergy2;

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
     * @param x the coordinate array.
     * @param g the gradient array.
     * @return the DualTopologyEnergy total energy.
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {

        /**
         * Update the coordinates of both topologies.
         */
        unpackCoordinates(x);

        /*
         try {
         energyRegion.setGradient(true);
         dualTopologyTeam.execute(energyRegion);
         } catch (Exception e) {
         String message = "Fatal exception computing the Dual-Topology Energy.\n";
         logger.log(Level.SEVERE, message, e);
         } */
        /**
         * Initialize dUdXdL arrays.
         */
        fill(gl1, 0.0);
        fill(gl2, 0.0);
        fill(rgl1, 0.0);
        fill(rgl2, 0.0);

        /**
         * Compute the energy and gradient of topology 1.
         */
        energy1 = potential1.energyAndGradient(x1, g1);
        dEdL_1 = lambdaInterface1.getdEdL();
        d2EdL2_1 = lambdaInterface1.getd2EdL2();
        lambdaInterface1.getdEdXdL(gl1);
        if (doValenceRestraint1) {
            forceFieldEnergy1.setLambdaBondedTerms(true);
            restraintEnergy1 = forceFieldEnergy1.energyAndGradient(x1, rg1);
            restraintdEdL_1 = forceFieldEnergy1.getdEdL();
            restraintd2EdL2_1 = forceFieldEnergy1.getd2EdL2();
            forceFieldEnergy1.getdEdXdL(rgl1);
            forceFieldEnergy1.setLambdaBondedTerms(false);
        } else {
            restraintEnergy1 = 0.0;
            restraintdEdL_1 = 0.0;
            restraintd2EdL2_1 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n",
                    lambdaPow * energy1, oneMinusLambdaPow * restraintEnergy1));
        }

        /**
         * Compute the energy and gradient of topology 2.
         */
        energy2 = potential2.energyAndGradient(x2, g2);
        dEdL_2 = -lambdaInterface2.getdEdL();
        d2EdL2_2 = lambdaInterface2.getd2EdL2();
        lambdaInterface2.getdEdXdL(gl2);
        if (doValenceRestraint2) {
            forceFieldEnergy2.setLambdaBondedTerms(true);
            restraintEnergy2 = forceFieldEnergy2.energyAndGradient(x2, rg2);
            restraintdEdL_2 = -forceFieldEnergy2.getdEdL();
            restraintd2EdL2_2 = forceFieldEnergy2.getd2EdL2();
            forceFieldEnergy2.getdEdXdL(rgl2);
            forceFieldEnergy2.setLambdaBondedTerms(false);
        } else {
            restraintEnergy2 = 0.0;
            restraintdEdL_2 = 0.0;
            restraintd2EdL2_2 = 0.0;
        }
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Topology 2 Energy & Restraints: %15.8f %15.8f\n",
                    oneMinusLambdaPow * energy2, lambdaPow * restraintEnergy2));
        }

        /**
         * Apply the dual-topology scaling for the total energy.
         */
        totalEnergy = lambdaPow * energy1 + oneMinusLambdaPow * restraintEnergy1
                + oneMinusLambdaPow * energy2 + lambdaPow * restraintEnergy2;

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
        for (int i = 0; i < nActive1; i++) {
            Atom a = activeAtoms1[i];
            if (!a.applyLambda()) {
                g[indexCommon++] = lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++];
                g[indexCommon++] = lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++];
                g[indexCommon++] = lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++];
            } else {
                g[indexUnique++] = lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++];
                g[indexUnique++] = lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++];
                g[indexUnique++] = lambdaPow * g1[index] + oneMinusLambdaPow * rg1[index++];
            }
        }
        /**
         * Coordinate Gradient from Topology 2.
         */
        indexCommon = 0;
        index = 0;
        for (int i = 0; i < nActive2; i++) {
            Atom a = activeAtoms2[i];
            if (!a.applyLambda()) {
                g[indexCommon++] += oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++];
                g[indexCommon++] += oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++];
                g[indexCommon++] += oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++];
            } else {
                g[indexUnique++] = oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++];
                g[indexUnique++] = oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++];
                g[indexUnique++] = oneMinusLambdaPow * g2[index] + lambdaPow * rg2[index++];
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
    
    @Override
    public void reInit() {
        throw new UnsupportedOperationException(String.format(" No reInit method defined for %s", DualTopologyEnergy.class.toString()));
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
        for (int i = 0; i < nActive1; i++) {
            Atom a = activeAtoms1[i];
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
        for (int i = 0; i < nActive2; i++) {
            Atom a = activeAtoms2[i];
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
        for (int i = 0; i < nActive1; i++) {
            Atom a = activeAtoms1[i];

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
        for (int i = 0; i < nActive2; i++) {
            Atom a = activeAtoms2[i];
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
        potential1.setEnergyTermState(state);
        potential2.setEnergyTermState(state);
    }

    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            oneMinusLambda = 1.0 - lambda;
            lambdaInterface1.setLambda(lambda);
            lambdaInterface2.setLambda(oneMinusLambda);

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
        double e1 = lambdaPow * dEdL_1 + dLambdaPow * energy1
                + oneMinusLambdaPow * restraintdEdL_1 + dOneMinusLambdaPow * restraintEnergy1;
        double e2 = oneMinusLambdaPow * dEdL_2 + dOneMinusLambdaPow * energy2
                + lambdaPow * restraintdEdL_2 + dLambdaPow * restraintEnergy2;
        return e1 + e2;
    }

    @Override
    public double getd2EdL2() {
        double e1 = lambdaPow * d2EdL2_1 + 2.0 * dLambdaPow * dEdL_1 + d2LambdaPow * energy1
                + oneMinusLambdaPow * restraintd2EdL2_1 + 2.0 * dOneMinusLambdaPow * restraintdEdL_1
                + d2OneMinusLambdaPow * restraintEnergy1;
        double e2 = oneMinusLambdaPow * d2EdL2_2 + 2.0 * dOneMinusLambdaPow * dEdL_2
                + d2OneMinusLambdaPow * energy2
                + lambdaPow * restraintd2EdL2_2 + 2.0 * dLambdaPow * restraintdEdL_2 + d2LambdaPow * restraintEnergy2;
        return e1 + e2;
    }

    @Override
    public void getdEdXdL(double[] g) {
        if (g == null) {
            g = new double[nVariables];
        }

        int index = 0;
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        /**
         * Coordinate Gradient from Topology 1.
         */
        for (int i = 0; i < nActive1; i++) {
            Atom a = activeAtoms1[i];
            if (!a.applyLambda()) {
                g[indexCommon++] = lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + oneMinusLambdaPow * rgl1[index] + dOneMinusLambdaPow * rg1[index++];
                g[indexCommon++] = lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + oneMinusLambdaPow * rgl1[index] + dOneMinusLambdaPow * rg1[index++];
                g[indexCommon++] = lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + oneMinusLambdaPow * rgl1[index] + dOneMinusLambdaPow * rg1[index++];
            } else {
                g[indexUnique++] = lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + oneMinusLambdaPow * rgl1[index] + dOneMinusLambdaPow * rg1[index++];
                g[indexUnique++] = lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + oneMinusLambdaPow * rgl1[index] + dOneMinusLambdaPow * rg1[index++];
                g[indexUnique++] = lambdaPow * gl1[index] + dLambdaPow * g1[index]
                        + oneMinusLambdaPow * rgl1[index] + dOneMinusLambdaPow * rg1[index++];
            }
        }

        /**
         * Coordinate Gradient from Topology 2.
         */
        index = 0;
        indexCommon = 0;
        for (int i = 0; i < nActive2; i++) {
            Atom a = activeAtoms2[i];
            if (a.isActive()) {
                if (!a.applyLambda()) {
                    g[indexCommon++] += (-oneMinusLambdaPow * gl2[index] + dOneMinusLambdaPow * g2[index]
                            - lambdaPow * rgl2[index] + dLambdaPow * rg2[index++]);
                    g[indexCommon++] += (-oneMinusLambdaPow * gl2[index] + dOneMinusLambdaPow * g2[index]
                            - lambdaPow * rgl2[index] + dLambdaPow * rg2[index++]);
                    g[indexCommon++] += (-oneMinusLambdaPow * gl2[index] + dOneMinusLambdaPow * g2[index]
                            - lambdaPow * rgl2[index] + dLambdaPow * rg2[index++]);
                } else {
                    g[indexUnique++] = (-oneMinusLambdaPow * gl2[index] + dOneMinusLambdaPow * g2[index]
                            - lambdaPow * rgl2[index] + dLambdaPow * rg2[index++]);
                    g[indexUnique++] = (-oneMinusLambdaPow * gl2[index] + dOneMinusLambdaPow * g2[index]
                            - lambdaPow * rgl2[index] + dLambdaPow * rg2[index++]);
                    g[indexUnique++] = (-oneMinusLambdaPow * gl2[index] + dOneMinusLambdaPow * g2[index]
                            - lambdaPow * rgl2[index] + dLambdaPow * rg2[index++]);
                }
            }
        }
    }

    @Override
    public void setVelocity(double[] velocity) {
        double vel[] = new double[3];
        int indexCommon = 0;
        int indexUnique = 3 * nShared;
        for (int i = 0; i < nActive1; i++) {
            Atom atom = activeAtoms1[i];
            if (!atom.applyLambda()) {
                vel[0] = velocity[indexCommon++];
                vel[1] = velocity[indexCommon++];
                vel[2] = velocity[indexCommon++];
            } else {
                vel[0] = velocity[indexUnique++];
                vel[1] = velocity[indexUnique++];
                vel[2] = velocity[indexUnique++];
            }
            atom.setVelocity(vel);
        }

        indexCommon = 0;
        for (int i = 0; i < nActive2; i++) {
            Atom atom = activeAtoms2[i];
            if (!atom.applyLambda()) {
                vel[0] = velocity[indexCommon++];
                vel[1] = velocity[indexCommon++];
                vel[2] = velocity[indexCommon++];
            } else {
                vel[0] = velocity[indexUnique++];
                vel[1] = velocity[indexUnique++];
                vel[2] = velocity[indexUnique++];
            }
            atom.setVelocity(vel);
        }
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        double accel[] = new double[3];
        int indexCommon = 0;
        int indexUnique = 3 * nShared;
        for (int i = 0; i < nActive1; i++) {
            Atom atom = activeAtoms1[i];
            if (!atom.applyLambda()) {
                accel[0] = acceleration[indexCommon++];
                accel[1] = acceleration[indexCommon++];
                accel[2] = acceleration[indexCommon++];
            } else {
                accel[0] = acceleration[indexUnique++];
                accel[1] = acceleration[indexUnique++];
                accel[2] = acceleration[indexUnique++];
            }
            atom.setAcceleration(accel);
        }
        indexCommon = 0;
        for (int i = 0; i < nActive2; i++) {
            Atom atom = activeAtoms2[i];
            if (!atom.applyLambda()) {
                accel[0] = acceleration[indexCommon++];
                accel[1] = acceleration[indexCommon++];
                accel[2] = acceleration[indexCommon++];
            } else {
                accel[0] = acceleration[indexUnique++];
                accel[1] = acceleration[indexUnique++];
                accel[2] = acceleration[indexUnique++];
            }
            atom.setAcceleration(accel);
        }
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        double prev[] = new double[3];
        int indexCommon = 0;
        int indexUnique = 3 * nShared;
        for (int i = 0; i < nActive1; i++) {
            Atom atom = activeAtoms1[i];
            if (!atom.applyLambda()) {
                prev[0] = previousAcceleration[indexCommon++];
                prev[1] = previousAcceleration[indexCommon++];
                prev[2] = previousAcceleration[indexCommon++];
            } else {
                prev[0] = previousAcceleration[indexUnique++];
                prev[1] = previousAcceleration[indexUnique++];
                prev[2] = previousAcceleration[indexUnique++];
            }
            atom.setPreviousAcceleration(prev);
        }
        indexCommon = 0;
        for (int i = 0; i < nActive2; i++) {
            Atom atom = activeAtoms2[i];
            if (!atom.applyLambda()) {
                prev[0] = previousAcceleration[indexCommon++];
                prev[1] = previousAcceleration[indexCommon++];
                prev[2] = previousAcceleration[indexCommon++];
            } else {
                prev[0] = previousAcceleration[indexUnique++];
                prev[1] = previousAcceleration[indexUnique++];
                prev[2] = previousAcceleration[indexUnique++];
            }
            atom.setPreviousAcceleration(prev);
        }
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        if (velocity == null || velocity.length < nVariables) {
            velocity = new double[nVariables];
        }
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        double vel[] = new double[3];
        for (int i = 0; i < nActive1; i++) {
            Atom atom = activeAtoms1[i];
            atom.getVelocity(vel);
            if (!atom.applyLambda()) {
                velocity[indexCommon++] = vel[0];
                velocity[indexCommon++] = vel[1];
                velocity[indexCommon++] = vel[2];
            } else {
                velocity[indexUnique++] = vel[0];
                velocity[indexUnique++] = vel[1];
                velocity[indexUnique++] = vel[2];
            }

        }
        for (int i = 0; i < nActive2; i++) {
            Atom atom = activeAtoms2[i];
            if (atom.applyLambda()) {
                atom.getVelocity(vel);
                velocity[indexUnique++] = vel[0];
                velocity[indexUnique++] = vel[1];
                velocity[indexUnique++] = vel[2];
            }
        }

        return velocity;
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        if (acceleration == null || acceleration.length < nVariables) {
            acceleration = new double[nVariables];
        }
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        double accel[] = new double[3];
        for (int i = 0; i < nActive1; i++) {
            Atom atom = activeAtoms1[i];
            atom.getAcceleration(accel);
            if (!atom.applyLambda()) {
                acceleration[indexCommon++] = accel[0];
                acceleration[indexCommon++] = accel[1];
                acceleration[indexCommon++] = accel[2];
            } else {
                acceleration[indexUnique++] = accel[0];
                acceleration[indexUnique++] = accel[1];
                acceleration[indexUnique++] = accel[2];
            }

        }
        for (int i = 0; i < nActive2; i++) {
            Atom atom = activeAtoms2[i];
            if (atom.applyLambda()) {
                atom.getAcceleration(accel);
                acceleration[indexUnique++] = accel[0];
                acceleration[indexUnique++] = accel[1];
                acceleration[indexUnique++] = accel[2];
            }
        }

        return acceleration;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        if (previousAcceleration == null || previousAcceleration.length < nVariables) {
            previousAcceleration = new double[nVariables];
        }
        int indexCommon = 0;
        int indexUnique = nShared * 3;
        double prev[] = new double[3];
        for (int i = 0; i < nActive1; i++) {
            Atom atom = activeAtoms1[i];
            atom.getPreviousAcceleration(prev);
            if (!atom.applyLambda()) {
                previousAcceleration[indexCommon++] = prev[0];
                previousAcceleration[indexCommon++] = prev[1];
                previousAcceleration[indexCommon++] = prev[2];
            } else {
                previousAcceleration[indexUnique++] = prev[0];
                previousAcceleration[indexUnique++] = prev[1];
                previousAcceleration[indexUnique++] = prev[2];
            }

        }
        for (int i = 0; i < nActive2; i++) {
            Atom atom = activeAtoms2[i];
            if (atom.applyLambda()) {
                atom.getPreviousAcceleration(prev);
                previousAcceleration[indexUnique++] = prev[0];
                previousAcceleration[indexUnique++] = prev[1];
                previousAcceleration[indexUnique++] = prev[2];
            }
        }

        return previousAcceleration;
    }

    /*private class EnergyRegion extends ParallelRegion {

        private final EnergySection1 e1 = new EnergySection1();
        private final EnergySection2 e2 = new EnergySection2();
        boolean gradient = false;

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        @Override
        public void run() throws Exception {
            try {
                execute(e1, e2);
                gradient = false;
            } catch (Exception e) {
                String message = "Fatal exception computing the Dual-Topology Energy.\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class EnergySection1 extends ParallelSection {

            @Override
            public void run() throws Exception {
                if (!gradient) {
                    energy1 = potential1.energy(x1);
                    if (doValenceRestraint1 && potential1 instanceof ForceFieldEnergy) {
                        ForceFieldEnergy ffE1 = (ForceFieldEnergy) potential1;
                        ffE1.setLambdaBondedTerms(true);
                        restraintEnergy1 = potential1.energy(x1);
                        ffE1.setLambdaBondedTerms(false);
                    } else {
                        restraintEnergy1 = 0.0;
                    }
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info(String.format(" Topology 1 Energy & Restraints: %15.8f %15.8f\n",
                                lambdaPow * energy1, oneMinusLambdaPow * restraintEnergy1));
                    }
                } else {
                    fill(gl1, 0.0);
                    fill(rgl1, 0.0);
                    energy1 = potential1.energyAndGradient(x1, g1);
                    dEdL_1 = lambdaInterface1.getdEdL();
                    d2EdL2_1 = lambdaInterface1.getd2EdL2();
                    lambdaInterface1.getdEdXdL(gl1);
                    if (doValenceRestraint1) {
                        forceFieldEnergy1.setLambdaBondedTerms(true);
                        restraintEnergy1 = forceFieldEnergy1.energyAndGradient(x1, rg1);
                        restraintdEdL_1 = forceFieldEnergy1.getdEdL();
                        restraintd2EdL2_1 = forceFieldEnergy1.getd2EdL2();
                        forceFieldEnergy1.getdEdXdL(rgl1);
                        forceFieldEnergy1.setLambdaBondedTerms(false);
                    } else {
                        restraintEnergy1 = 0.0;
                        restraintdEdL_1 = 0.0;
                        restraintd2EdL2_1 = 0.0;
                    }
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info(String.format(" Topology 1 Gradient & Restraints: %15.8f %15.8f\n",
                                lambdaPow * energy1, oneMinusLambdaPow * restraintEnergy1));
                    }
                }
            }
        }

        private class EnergySection2 extends ParallelSection {

            @Override
            public void run() throws Exception {
                if (!gradient) {
                    energy2 = potential2.energy(x2);
                    if (doValenceRestraint2 && potential2 instanceof ForceFieldEnergy) {
                        ForceFieldEnergy ffE2 = (ForceFieldEnergy) potential2;
                        ffE2.setLambdaBondedTerms(true);
                        restraintEnergy2 = potential2.energy(x2);
                        ffE2.setLambdaBondedTerms(false);
                    } else {
                        restraintEnergy2 = 0.0;
                    }
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info(String.format(" Topology 2 Energy & Restraints: %15.8f %15.8f\n",
                                oneMinusLambdaPow * energy2, lambdaPow * restraintEnergy2));
                    }
                } else {
                    fill(gl2, 0.0);
                    fill(rgl2, 0.0);
                    energy2 = potential2.energyAndGradient(x2, g2);
                    dEdL_2 = -lambdaInterface2.getdEdL();
                    d2EdL2_2 = lambdaInterface2.getd2EdL2();
                    lambdaInterface2.getdEdXdL(gl2);

                    if (doValenceRestraint2) {
                        forceFieldEnergy2.setLambdaBondedTerms(true);
                        restraintEnergy2 = forceFieldEnergy2.energyAndGradient(x2, rg2);
                        restraintdEdL_2 = -forceFieldEnergy2.getdEdL();
                        restraintd2EdL2_2 = forceFieldEnergy2.getd2EdL2();
                        forceFieldEnergy2.getdEdXdL(rgl2);
                        forceFieldEnergy2.setLambdaBondedTerms(false);
                    } else {
                        restraintEnergy2 = 0.0;
                        restraintdEdL_2 = 0.0;
                        restraintd2EdL2_2 = 0.0;
                    }
                    if (logger.isLoggable(Level.INFO)) {
                        logger.info(String.format(" Topology 2 Gradient & Restraints: %15.8f %15.8f\n",
                                oneMinusLambdaPow * energy2, lambdaPow * restraintEnergy2));
                    }
                }
            }
        }
    }*/
}
