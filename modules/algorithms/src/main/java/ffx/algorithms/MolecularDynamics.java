/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.algorithms;

import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Random;

import ffx.potential.PotentialEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;

/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class MolecularDynamics implements Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    /**
     * Boltzmann constant in units of g*Ang**2/ps**2/mole/K
     */
    private static final double kB = 0.83144725;
    /**
     * Conversion from kcal/mole to g*Ang**2/ps**2.
     */
    private static final double convert = 4.1840e2;
    private final int n;
    private final int dof;
    private final double[] x;
    private final double[] xPrevious;
    private final double[] v;
    private final double[] a;
    private final double[] aPrevious;
    private final double[] grad;
    private final double[] mass;
    private double targetTemp;
    private double currentTemp;
    private double kinetic;
    private double potential;
    private double total;
    private double dt;
    private final MolecularAssembly molecularAssembly;
    private final PotentialEnergy potentialEnergy;
    private final Atom atoms[];
    private AlgorithmListener algorithmListener;
    private boolean done;
    private boolean terminate;

    public MolecularDynamics(MolecularAssembly assembly, AlgorithmListener listener) {
        this.molecularAssembly = assembly;
        this.algorithmListener = listener;
        if (molecularAssembly.getPotentialEnergy() == null) {
            potentialEnergy = new PotentialEnergy(molecularAssembly);
            molecularAssembly.setPotential(potentialEnergy);
        } else {
            potentialEnergy = molecularAssembly.getPotentialEnergy();
        }
        ArrayList<Atom> atomList = molecularAssembly.getAtomList();
        n = atomList.size();
        atoms = atomList.toArray(new Atom[n]);
        Arrays.sort(atoms);

        mass = new double[n];
        dof = n * 3;
        x = new double[dof];
        xPrevious = new double[dof];
        v = new double[dof];
        a = new double[dof];
        aPrevious = new double[dof];
        grad = new double[dof];
    }

    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double temperature, final boolean initVelocities) {
        terminate = false;
        done = false;
        /**
         * Convert the time step from femtoseconds to picoseconds.
         */
        dt = timeStep * 1.0e-3;

        int printFrequency = 1;
        if (printInterval > dt) {
            printFrequency = (int) (printInterval / dt);
        }

        /**
         * Set the target temperature.
         */
        if (temperature > 0.0) {
            targetTemp = temperature;
        } else {
            targetTemp = 0.0;
        }

        /**
         * Initialize atomic positions and masses.
         */
        int j = 0;
        for (int i = 0; i < n; i++) {
            Atom atom = atoms[i];
            mass[i] = atom.getMass();
            double xyz[] = atom.getXYZ();
            x[j++] = xyz[0];
            x[j++] = xyz[1];
            x[j++] = xyz[2];
        }

        /**
         * Initialize velocities.
         */
        if (targetTemp != 0.0 && initVelocities) {
            maxwell();
        } else {
            for (int i = 0; i < dof; i++) {
                v[i] = 0.0;
            }
        }

        if (targetTemp > 0.0) {
            kinetic = kineticEnergy();
        } else {
            kinetic = 0.0;
        }

        potentialEnergy.setOptimizationScaling(null);
        potential = potentialEnergy.energyAndGradient(x, grad);
        total = kinetic + potential;

        /**
         * Initialize current and previous accelerations.
         */
        int index = 0;
        for (int i = 0; i < n; i++) {
            double m = mass[i];
            for (j = 0; j < 3; j++, index++) {
                a[index] = -convert * grad[index] / m;
                aPrevious[index] = a[index];
            }
        }

        final double kT = kB * targetTemp;
        logger.info(String.format("\n Initial kT per degree of freedom:  %7.3f\n", convert * kinetic / (dof * kT)));
        logger.info(String.format("\n   Step      Kinetic    Potential        Total     Temp     Time"));
        logger.info(String.format("       %13.4f%13.4f%13.4f %8.2f ", kinetic, potential, total, currentTemp));

        /**
         * Integrate Newton's equations of motion for the requested number of
         * steps, unless early termination is requested.
         */
        long time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {
            beeman(dt);
            kinetic = kineticEnergy();
            total = kinetic + potential;
            if (step % printFrequency == 0) {
                time = System.nanoTime() - time;
                logger.info(String.format(" %6d%13.4f%13.4f%13.4f%9.2f%9.3f", step, kinetic, potential, total, currentTemp, time * 1.0e-9));
                time = System.nanoTime();
            }
            if (algorithmListener != null && step % printFrequency == 0) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            if (terminate) {
                logger.info(String.format("\n Terminating after %6d time steps.\n", step));
                done = true;
                break;
            }
        }

        if (!terminate) {
            logger.info(String.format("\n Final kT per degree of freedom:    %7.3f\n", convert * kinetic / (dof * kT)));
            logger.info(String.format(" Completed %6d time steps.\n", nSteps));
        }
    }

    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Exception terminating minimization.\n", e);
                }
            }
        }
    }

    /**
     * Integrate Newton's equations of motion using a Beeman multistep recursion
     * formula; the actual coefficients are Brooks' "Better Beeman" values.
     *
     * @since 1.0
     */
    public void beeman(final double dt) {
        final double dt_8 = 0.125 * dt;
        final double dt2_8 = dt * dt_8;

        /**
         * Store the current atom positions, then find new atom positions
         * and half-step velocities via Beeman recusion.
         */
        for (int i = 0; i < dof; i++) {
            xPrevious[i] = x[i];
            double temp = 5.0 * a[i] - aPrevious[i];
            x[i] += v[i] * dt + temp * dt2_8;
            v[i] += temp * dt_8;
        }

        /**
         * Compute the potential energy and gradients.
         */
        potential = potentialEnergy.energyAndGradient(x, grad);

        /**
         * Use Newton's second law to get the next acceleration and find
         * the full-step velocities using the Beeman recusion.
         */
        for (int index = 0, i = 0; i < n; i++) {
            double m = mass[i];
            for (int j = 0; j < 3; j++, index++) {
                aPrevious[index] = a[index];
                a[index] = -convert * grad[index] / m;
                v[index] += (3.0 * a[index] + aPrevious[index]) * dt_8;
            }
        }
    }

    /**
     * Initialize velocities from a Maxwell-Boltzmann distribution of momenta.
     * The varience of each independent momentum component is kT * mass.
     */
    private void maxwell() {
        Random random = new Random();
        final double kT = kB * targetTemp;
        int index = 0;
        for (int i = 0; i < n; i++) {
            double m = mass[i];
            double variance = m * kT;
            double sd = sqrt(variance);
            v[index++] = random.nextGaussian() * sd / m;
            v[index++] = random.nextGaussian() * sd / m;
            v[index++] = random.nextGaussian() * sd / m;
        }
    }

    /**
     * Compute the kinetic energy of the system.
     */
    private double kineticEnergy() {
        double e = 0.0;
        int index = 0;
        for (int i = 0; i < n; i++) {
            double velocity = v[index++];
            double v2 = velocity * velocity;
            velocity = v[index++];
            v2 += velocity * velocity;
            velocity = v[index++];
            v2 += velocity * velocity;
            e += mass[i] * v2;
        }
        e *= 0.5 / convert;
        currentTemp = convert * e / (0.5 * kB * dof);
        return e;
    }
}
