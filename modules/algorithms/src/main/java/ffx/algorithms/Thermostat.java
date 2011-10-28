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
import static java.lang.Math.sqrt;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import ffx.numerics.Potential.VARIABLE_TYPE;

/**
 * The abstract Thermostat class implements methods common to all thermostats
 * for initalizing velocities from a Maxwell-Boltzmann distribution and
 * computing the instantaneous temperature. Abstract methods are declared for
 * half-step and full-step modification of velocities thermostat implementions.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
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
    protected double currentKineticEnergy;
    protected int nVariables;
    protected int dof;
    protected boolean removingCenterOfMassMotion;
    protected double x[];
    protected double v[];
    protected double mass[];
    protected VARIABLE_TYPE type[];
    protected double totalMass;
    protected final double centerOfMass[] = new double[3];
    protected final double linearMomentum[] = new double[3];
    protected final double angularMomentum[] = new double[3];
    protected Random random;
    protected Thermostats name;

    /**
     * <p>Constructor for Thermostat.</p>
     *
     * @param nVariables a int.
     * @param x an array of double.
     * @param v an array of double.
     * @param mass an array of double.
     * @param t a double.
     */
    public Thermostat(int nVariables, double x[], double v[], double mass[],
                      VARIABLE_TYPE type[], double t) {
        assert (nVariables > 3);

        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.mass = mass;
        this.type = type;
        assert (x.length == nVariables);
        assert (v.length == nVariables);
        assert (mass.length == nVariables);
        assert (type.length == nVariables);
        random = new Random();
        setTargetTemperature(t);
        /**
         * Set the degrees of freedom to nVariables - 3 because we will
         * remove center of mass motion.
         */
        removingCenterOfMassMotion = true;
        dof = nVariables - 3;
    }

    /**
     * If center of mass motion is being removed, then the mean
     * kinetic energy of the system will be 3 * kT/2 less than if
     * center of mass motion is allowed.
     *
     * @param remove <code>true</code> if center of mass motion is being removed.
     */
    public void removingCenterOfMassMotion(boolean remove) {
        removingCenterOfMassMotion = remove;
        if (removingCenterOfMassMotion) {
            dof = nVariables - 3;
        } else {
            dof = nVariables;
        }
    }

    public boolean removingCOM() {
        return removingCenterOfMassMotion;
    }

    /**
     * <p>The setRandomSeed method is used to initialize the Random number
     * generator to the same starting state, such that separate runs
     * produce the same Maxwell-Boltzmann initial velocities.
     * same </p>
     *
     * @param seed The seed.
     */
    public void setRandomSeed(long seed) {
        random.setSeed(seed);
    }

    /**
     * <p>Log the target temperature and current number of kT per
     * degree of freedom (should be 0.5 kT at equilibrium).</p>
     *
     * @param level a {@link java.util.logging.Level} object.
     */
    protected void log(Level level) {
        if (logger.isLoggable(level)) {
            logger.log(level, "\n" + toString());
            logger.log(level, String.format(" Target temperature:           %7.2f Kelvin", targetTemperature));
            logger.log(level, String.format(" Current temperature:          %7.2f Kelvin", currentTemperature));
            logger.log(level, String.format(" Number of variables:          %7d", nVariables));
            logger.log(level, String.format(" Number of degrees of freedom: %7d", dof));
            logger.log(level, String.format(" kT per degree of freedom:     %7.2f", convert * currentKineticEnergy / (dof * kT)));
        }
    }

    /**
     * <p>getCurrentTemperture</p>
     *
     * @return a double.
     */
    public double getCurrentTemperture() {
        return currentTemperature;
    }

    /**
     * <p>Getter for the field <code>kineticEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getKineticEnergy() {
        return currentKineticEnergy;
    }

    /**
     * <p>Getter for the field <code>targetTemperature</code>.</p>
     *
     * @return a double.
     */
    public double getTargetTemperature() {
        return targetTemperature;
    }

    /**
     * Set the target temperature.
     *
     * @param t Target temperature must be greater than absolute zero.
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
     * The variance of each independent momentum component is kT * mass.
     */
    public void maxwell(double targetTemperature) {
        for (int i = 0; i < nVariables; i++) {
            double m = mass[i];
            v[i] = random.nextGaussian() * sqrt(kB * targetTemperature / m);
        }

        /**
         * Remove the center of mass motion.
         */
        if (removingCenterOfMassMotion) {
            centerOfMassMotion(true, true);
        }

        /**
         * Find the current kinetic energy and temperature.
         */
        kineticEnergy();

        /**
         * The current temperature will deviate slightly from the target
         * temperature if the center of mass motion was removed and/or
         * due to finite system size.
         *
         * Scale the velocities to reach the target temperature.
         */
        double scale = Math.sqrt(targetTemperature / currentTemperature);
        for (int i = 0; i < nVariables; i++) {
            v[i] *= scale;
        }

        /**
         * Update the kinetic energy and current temperature.
         */
        kineticEnergy();

        log(Level.INFO);
    }

    public void maxwell() {
        maxwell(targetTemperature);
    }

    /**
     * <p>centerOfMassMotion</p>
     *
     * @param remove a boolean.
     * @param print a boolean.
     */
    protected void centerOfMassMotion(boolean remove, boolean print) {
        totalMass = 0.0;
        for (int i = 0; i < 3; i++) {
            centerOfMass[i] = 0.0;
            linearMomentum[i] = 0.0;
            angularMomentum[i] = 0.0;
        }

        int index = 0;
        while (index < nVariables) {
            if (type[index] == VARIABLE_TYPE.OTHER) {
                index++;
                continue;
            }
            assert (type[index] == VARIABLE_TYPE.X);
            double m = mass[index];
            double xx = x[index];
            double vx = v[index++];
            assert (type[index] == VARIABLE_TYPE.Y);
            double yy = x[index];
            double vy = v[index++];
            assert (type[index] == VARIABLE_TYPE.Z);
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
            logger.info(String.format("\n Center of Mass   (%12.3f,%12.3f,%12.3f)", centerOfMass[0], centerOfMass[1], centerOfMass[2]));
            logger.info(String.format(" Linear Momentum  (%12.3f,%12.3f,%12.3f)", linearMomentum[0], linearMomentum[1], linearMomentum[2]));
            logger.info(String.format(" Angular Momentum (%12.3f,%12.3f,%12.3f)", angularMomentum[0], angularMomentum[1], angularMomentum[2]));
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
        int index = 0;
        while (index < nVariables) {
            if (type[index] == VARIABLE_TYPE.OTHER) {
                index++;
                continue;
            }
            double m = mass[index];
            assert (type[index] == VARIABLE_TYPE.X);
            double xi = x[index++] - centerOfMass[0];
            assert (type[index] == VARIABLE_TYPE.Y);
            double yi = x[index++] - centerOfMass[1];
            assert (type[index] == VARIABLE_TYPE.Z);
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
        index = 0;
        while (index < nVariables) {
            if (type[index] == VARIABLE_TYPE.OTHER) {
                index++;
                continue;
            }
            v[index++] -= linearMomentum[0];
            v[index++] -= linearMomentum[1];
            v[index++] -= linearMomentum[2];
        }

        /**
         * Only remove center of mass rotational momentum for
         * non-periodic systems.
         */
        if (false) {
            index = 0;
            while (index < nVariables) {
                if (type[index] == VARIABLE_TYPE.OTHER) {
                    index++;
                    continue;
                }
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
            logger.info(String.format(" Center of mass motion removed."));
        }
    }

    /**
     * Compute the current temperature and kinetic energy of the system.
     */
    protected void kineticEnergy() {
        double e = 0.0;
        for (int i = 0; i < nVariables; i++) {
            double velocity = v[i];
            double v2 = velocity * velocity;
            e += mass[i] * v2;
        }
        currentTemperature = e / (kB * dof);
        e *= 0.5 / convert;
        currentKineticEnergy = e;
    }

    /**
     * <p>halfStep</p>
     *
     * @param dt a double.
     */
    public abstract void halfStep(double dt);

    /**
     * <p>fullStep</p>
     *
     * @param dt a double.
     */
    public abstract void fullStep(double dt);
}
