/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

/**
 * The abstract Thermostat class implements methods common to all thermostats
 * for initalizing velocities from a Maxwell-Boltzmann distribution and
 * computing the instantaneous temperature. Abstract methods are declared for
 * half-step and full-step modification of velocities thermostat implementions.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public abstract class Thermostat {

    private static final Logger logger = Logger.getLogger(Thermostat.class.getName());
    /**
     * Boltzmann constant in units of g*Ang**2/ps**2/mole/K
     */
    public static final double kB = 0.83144725;
    /**
     * Conversion from kcal/mole to g*Ang**2/ps**2.
     */
    public static final double convert = 4.1840e2;
    
    /**
     * Gas constant (in Kcal/mole/Kelvin).
     */
    public static final double R = 1.9872066e-3;
    
    public enum Thermostats {

        ADIABATIC, BERENDSEN, BUSSI;
    };
    protected double targetTemperature;
    protected double currentTemperature;
    protected double kT;
    protected double kineticEnergy;
    protected int dof;
    protected double x[];
    protected double v[];
    protected double mass[];
    protected double totalMass;
    protected final double centerOfMass[] = new double[3];
    protected final double linearMomentum[] = new double[3];
    protected final double angularMomentum[] = new double[3];
    protected Random random;
    protected Thermostats name;

    public Thermostat(int dof, double x[], double v[], double mass[], double t) {
        this.dof = dof;
        this.x = x;
        this.v = v;
        this.mass = mass;
        assert (x.length == dof);
        assert (v.length == dof);
        assert (mass.length == dof);
        random = new Random(0);
        setTargetTemperature(t);
    }

    public void setRandomSeed(long seed) {
        random.setSeed(seed);
    }
    
    protected void log(Level level) {
        if (logger.isLoggable(level)) {
            logger.log(level, String.format("\nThermostat target temperature %6.2.", targetTemperature));
            logger.log(level, String.format("\nkT per degree of freedom:  %7.3f\n", convert * kineticEnergy / (dof * kT)));
        }
    }

    public double getCurrentTemperture() {
        return currentTemperature;
    }

    public double getKineticEnergy() {
        return kineticEnergy;
    }

    public double getTargetTemperature() {
        return targetTemperature;
    }

    /**
     * Set the target temperature.
     *
     * @param t Target temperature must be greater than absolute zero.
     *
     * @since 1.0
     */
    public final void setTargetTemperature(double t) {
        /**
         * Obey the Third Law of Thermodynamics.
         */
        assert (t > 0.0);
        targetTemperature = t;
        kT = t * kB;
    }

    /**
     * Reset velocities from a Maxwell-Boltzmann distribution of momenta.
     * The varience of each independent momentum component is kT * mass.
     */
    public void maxwell() {

        for (int i = 0; i < dof; i++) {
            double m = mass[i];
            v[i++] = random.nextGaussian() * Math.sqrt(m * kT) / m;
        }
        centerOfMassMotion(true, true);
        kineticEnergy();
    }

    protected void centerOfMassMotion(boolean remove, boolean print) {
        totalMass = 0.0;
        for (int i = 0; i < 3; i++) {
            centerOfMass[i] = 0.0;
            linearMomentum[i] = 0.0;
            angularMomentum[i] = 0.0;
        }

        for (int index = 0, i = 0; i < dof/3; i++) {
            double m = mass[index];
            double xx = x[index];
            double vx = v[index++];
            double yy = x[index];
            double vy = v[index++];
            double zz = x[index];
            double vz = v[index++];
            totalMass += m;
            centerOfMass[0] += xx * m;
            centerOfMass[1] += yy * m;
            centerOfMass[2] += zz * m;
            linearMomentum[0] += vx * m;
            linearMomentum[1] += vy * m;
            linearMomentum[2] += vz * m;
            angularMomentum[0] += (yy * vz - zz * vy) * m;
            angularMomentum[1] += (zz * vx - xx * vz) * m;
            angularMomentum[2] += (xx * vy - yy * vx) * m;
        }

        angularMomentum[0] -= (centerOfMass[1] * linearMomentum[2] - centerOfMass[2] * linearMomentum[1]) / totalMass;
        angularMomentum[1] -= (centerOfMass[2] * linearMomentum[0] - centerOfMass[0] * linearMomentum[2]) / totalMass;
        angularMomentum[2] -= (centerOfMass[0] * linearMomentum[1] - centerOfMass[1] * linearMomentum[0]) / totalMass;
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;
        linearMomentum[0] /= totalMass;
        linearMomentum[1] /= totalMass;
        linearMomentum[2] /= totalMass;

        if (print) {
            logger.info(String.format(" Center of Mass   (%12.3f,%12.3f,%12.3f)", centerOfMass[0], centerOfMass[1], centerOfMass[2]));
            logger.info(String.format(" Linear Momemtum  (%12.3f,%12.3f,%12.3f)", linearMomentum[0], linearMomentum[1], linearMomentum[2]));
            logger.info(String.format(" Angular Momemtum (%12.3f,%12.3f,%12.3f)", angularMomentum[0], angularMomentum[1], angularMomentum[2]));
        }

        if (remove) {
            removeCenterOfMassMotion(print);
            centerOfMassMotion(false, print);
        }
    }

    /**
     * Remove center of mass translational and rotational velocity by
     * inverting the moment of inertia tensor.
     */
    private void removeCenterOfMassMotion(boolean print) {
        double xx = 0.0;
        double yy = 0.0;
        double zz = 0.0;
        double xy = 0.0;
        double xz = 0.0;
        double yz = 0.0;
        for (int index = 0, i = 0; i < dof/3; i++) {
            double m = mass[index];
            double xi = x[index++] - centerOfMass[0];
            double yi = x[index++] - centerOfMass[1];
            double zi = x[index++] - centerOfMass[2];
            xx += xi * xi * m;
            yy += yi * yi * m;
            zz += zi * zi * m;
            xy += xi * yi * m;
            xz += xi * zi * m;
            yz += yi * zi * m;
        }

        RealMatrix inertia = new Array2DRowRealMatrix(3, 3);
        inertia.setEntry(0, 0, yy + zz);
        inertia.setEntry(1, 0, -xy);
        inertia.setEntry(2, 0, -xz);
        inertia.setEntry(0, 1, -xy);
        inertia.setEntry(1, 1, xx + zz);
        inertia.setEntry(2, 1, -yz);
        inertia.setEntry(0, 2, -xz);
        inertia.setEntry(1, 2, -yz);
        inertia.setEntry(2, 2, xx + yy);
        inertia = new LUDecompositionImpl(inertia).getSolver().getInverse();
        xx = inertia.getEntry(0, 0);
        yy = inertia.getEntry(1, 1);
        zz = inertia.getEntry(2, 2);
        xy = inertia.getEntry(0, 1);
        xz = inertia.getEntry(0, 2);
        yz = inertia.getEntry(1, 2);
        double ox = angularMomentum[0] * xx + angularMomentum[1] * xy + angularMomentum[2] * xz;
        double oy = angularMomentum[0] * xy + angularMomentum[1] * yy + angularMomentum[2] * yz;
        double oz = angularMomentum[0] * xz + angularMomentum[1] * yz + angularMomentum[2] * zz;
        /**
         * Remove center of mass translational momentum.
         */
        for (int index = 0, i = 0; i < dof/3; i++) {
            v[index++] -= linearMomentum[0];
            v[index++] -= linearMomentum[1];
            v[index++] -= linearMomentum[2];
        }

        /**
         * Only remove center of mass rotational momentum for 
         * non-periodic systems.
         */
        if (false) {
            for (int index = 0, i = 0; i < dof/3; i++) {
                double xi = x[index++] - centerOfMass[0];
                double yi = x[index++] - centerOfMass[1];
                double zi = x[index] - centerOfMass[2];
                index -= 2;
                v[index++] += (-oy * zi + oz * yi);
                v[index++] += (-oz * xi + ox * zi);
                v[index++] += (-ox * yi + oy * xi);
            }
        }
        /**
         * Update the degrees of freedom.
         */
        if (print) {
            logger.info(String.format(" Center of mass motion removed.\n"
                                      + " Total degress of freedom %d - 3 = %d", dof, dof-3));
        }
    }

    /**
     * Compute the current temperature and kinetic energy of the system.
     */
    protected void kineticEnergy() {
        double e = 0.0;
        for (int i = 0; i < dof; i++) {
            double velocity = v[i];
            double v2 = velocity * velocity;
            e += mass[i] * v2;
        }
        currentTemperature = e / (kB * (dof - 3));
        e *= 0.5 / convert;
        kineticEnergy = e;
    }

    public abstract void halfStep(double dt);

    public abstract void fullStep(double dt);
}
