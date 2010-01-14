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

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

import no.uib.cipr.matrix.Matrix;

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

    private static final Logger logger = Logger.getLogger(ffx.algorithms.Thermostat.class.getName());
    /**
     * Boltzmann constant in units of g*Ang**2/ps**2/mole/K
     */
    public static final double kB = 0.83144725;
    /**
     * Conversion from kcal/mole to g*Ang**2/ps**2.
     */
    public static final double convert = 4.1840e2;

    public enum Thermostats {

        BERENDSEN, BUSSI, ISOTHERMAL
    };
    protected double targetTemperature;
    protected double currentTemperature;
    protected double kT;
    protected double kineticEnergy;
    protected int n;
    protected int dof;
    protected double x[];
    protected double v[];
    protected double mass[];
    protected Thermostats name;

    public Thermostat(int n, double x[], double v[], double mass[], double t) {
        this.n = n;
        this.dof = n * 3;
        this.x = x;
        this.v = v;
        this.mass = mass;
        assert (v.length == dof);
        assert (mass.length == n);
        setTargetTemperature(t);
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
    public void setTargetTemperature(double t) {
        /**
         * Obey the Third Law of Thermodynamics.
         */
        assert (t > 0.0);
        targetTemperature = t;
        kT = t * kB;
    }

    /**
     * Rest velocities from a Maxwell-Boltzmann distribution of momenta.
     * The varience of each independent momentum component is kT * mass.
     */
    public void maxwell() {
        Random random = new Random();
        for (int index = 0, i = 0; i < n; i++) {
            double m = mass[i];
            double variance = m * kT;
            double sd = Math.sqrt(variance);
            v[index++] = random.nextGaussian() * sd / m;
            v[index++] = random.nextGaussian() * sd / m;
            v[index++] = random.nextGaussian() * sd / m;
        }
        centerOfMassMotion(true);
        kineticEnergy();
    }
    protected double mTot;
    protected final double com[] = new double[3];
    protected final double lm[] = new double[3];
    protected final double am[] = new double[3];

    protected void centerOfMassMotion(boolean remove) {
        mTot = 0.0;
        for (int i = 0; i < 3; i++) {
            com[i] = 0.0;
            lm[i] = 0.0;
            am[i] = 0.0;
        }

        for (int index = 0, i = 0; i < n; i++) {
            double m = mass[i];
            double xx = x[index];
            double vx = v[index++];
            double yy = x[index];
            double vy = v[index++];
            double zz = x[index];
            double vz = v[index++];
            mTot += m;
            com[0] += xx * m;
            com[1] += yy * m;
            com[2] += zz * m;
            lm[0] += vx * m;
            lm[1] += vy * m;
            lm[2] += vz * m;
            am[0] += (yy * vz - zz * vy) * m;
            am[1] += (zz * vx - xx * vz) * m;
            am[2] += (xx * vy - yy * vx) * m;
        }

        am[0] -= (com[1] * lm[2] - com[2] * lm[1]) / mTot;
        am[1] -= (com[2] * lm[0] - com[0] * lm[2]) / mTot;
        am[2] -= (com[0] * lm[1] - com[1] * lm[0]) / mTot;
        com[0] /= mTot;
        com[1] /= mTot;
        com[2] /= mTot;
        lm[0] /= mTot;
        lm[1] /= mTot;
        lm[2] /= mTot;

        logger.info(String.format(" Center of Mass   (%12.3f,%12.3f,%12.3f)", com[0], com[1], com[2]));
        logger.info(String.format(" Linear Momemtum  (%12.3f,%12.3f,%12.3f)", lm[0], lm[1], lm[2]));
        logger.info(String.format(" Angular Momemtum (%12.3f,%12.3f,%12.3f)", am[0], am[1], am[2]));

        if (remove) {
            removeCenterOfMassMotion();
            centerOfMassMotion(false);
        }
    }

    /**
     * Remove center of mass translational and rotational velocity by
     * inverting the moment of intertia tensor.
     */
    private void removeCenterOfMassMotion() {
        double xx = 0.0;
        double yy = 0.0;
        double zz = 0.0;
        double xy = 0.0;
        double xz = 0.0;
        double yz = 0.0;
        for (int index = 0, i = 0; i < n; i++) {
            double m = mass[i];
            double xi = x[index++] - com[0];
            double yi = x[index++] - com[1];
            double zi = x[index++] - com[2];
            xx += xi * xi * m;
            yy += yi * yi * m;
            zz += zi * zi * m;
            xy += xi * yi * m;
            xz += xi * zi * m;
            yz += yi * zi * m;
        }

        DenseMatrix inertia = new DenseMatrix(3, 3);
        inertia.set(0, 0, yy + zz);
        inertia.set(1, 0, -xy);
        inertia.set(2, 0, -xz);
        inertia.set(0, 1, -xy);
        inertia.set(1, 1, xx + zz);
        inertia.set(2, 1, -yz);
        inertia.set(0, 2, -xz);
        inertia.set(1, 2, -yz);
        inertia.set(2, 2, xx + yy);

        DenseMatrix I = Matrices.identity(3);
        DenseMatrix AI = I.copy();
        inertia.solve(I, AI);

        xx = AI.get(0, 0);
        yy = AI.get(1, 1);
        zz = AI.get(2, 2);
        xy = AI.get(0, 1);
        xz = AI.get(0, 2);
        yz = AI.get(1, 2);
        double ox = am[0]*xx + am[1]*xy + am[2]*xz;
        double oy = am[0]*xy + am[1]*yy + am[2]*yz;
        double oz = am[0]*xz + am[1]*yz + am[2]*zz;
        /**
         * Remove center of mass translational and rotational velocity.
         */
        for (int index = 0, i = 0; i < n; i++) {
            double xi = x[index++] - com[0];
            double yi = x[index++] - com[1];
            double zi = x[index] - com[2];
            index -= 2;
            v[index++] += (-lm[0] - oy*zi + oz*yi);
            v[index++] += (-lm[1] - oz*xi + ox*zi);
            v[index++] += (-lm[2] - ox*yi + oy*xi);
        }
        /**
         * Update the degrees of freedom.
         */
        dof = 3 * n - 6;
        logger.info(String.format(" Center of mass motion removed.\n" +
                " Total degress of freedom %d - 6 = %d", 3*n, dof));
    }

    /**
     * Compute the current temperature and kinetic energy of the system.
     */
    protected void kineticEnergy() {
        double e = 0.0;
        for (int index = 0, i = 0; i < n; i++) {
            double velocity = v[index++];
            double v2 = velocity * velocity;
            velocity = v[index++];
            v2 += velocity * velocity;
            velocity = v[index++];
            v2 += velocity * velocity;
            e += mass[i] * v2;
        }
        e *= 0.5 / convert;
        currentTemperature = convert * e / (0.5 * kB * dof);
        kineticEnergy = e;
    }

    public abstract void halfStep(double dt);

    public abstract void fullStep(double dt);
}
