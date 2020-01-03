//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.dynamics.thermostats;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.Constraint;
import ffx.numerics.Potential.VARIABLE_TYPE;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;

/**
 * The abstract Thermostat class implements methods common to all thermostats
 * for initializing velocities from a Maxwell-Boltzmann distribution and
 * computing the instantaneous temperature. Abstract methods are declared for
 * half-step and full-step modification of velocities for thermostat
 * implementations.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class Thermostat {

    private static final Logger logger = Logger.getLogger(Thermostat.class.getName());

    /**
     * Parse a string into a Thermostat enumeration.
     *
     * @param str Thermostat String.
     * @return An instance of the ThermostatEnum.
     */
    public static ThermostatEnum parseThermostat(String str) {
        try {
            return ThermostatEnum.valueOf(str.toUpperCase());
        } catch (Exception e) {
            logger.info(format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
            return ThermostatEnum.BERENDSEN;
        }
    }

    /**
     * The identity of this Thermostat.
     */
    protected ThermostatEnum name;
    /**
     * The target temperature that this thermostat should maintain.
     */
    double targetTemperature;
    /**
     * The current temperature of the degrees of freedom.
     */
    double currentTemperature;
    /**
     * The value of kT in kcal/mol at the target temperature.
     */
    protected double kT;
    /**
     * The current kinetic energy of the system.
     */
    private double currentKineticEnergy;
    /**
     * Number of variables.
     */
    protected int nVariables;
    /**
     * Number of degrees of freedom, which can be less than the number of
     * variables. For example, removing translational motion removes 3 degrees
     * of freedom.
     */
    int dof;
    /**
     * Current values of variables.
     */
    protected double[] x;
    /**
     * Current velocity of the variables.
     */
    protected double[] v;
    /**
     * Mass for each variable.
     */
    protected double[] mass;
    /**
     * The type of each variable.
     */
    protected VARIABLE_TYPE[] type;
    /**
     * The center of mass coordinates.
     */
    private final double[] centerOfMass = new double[3];
    /**
     * The linear momentum.
     */
    private final double[] linearMomentum = new double[3];
    /**
     * The angular momentum.
     */
    private final double[] angularMomentum = new double[3];
    /**
     * Flag to indicate that center of mass motion should be removed.
     */
    private boolean removeCenterOfMassMotion;
    /**
     * The random number generator that the Thermostat will use to initialize velocities.
     */
    protected Random random;
    /**
     * Reduce logging.
     */
    private boolean quiet = false;
    /**
     * Any geometric constraints to apply during integration.
     */
    protected List<Constraint> constraints;
    /**
     * Number of degrees of freedom removed by constraints.
     */
    private int constrainedDoF;

    /**
     * <p>
     * Constructor for Thermostat.</p>
     *
     * @param n                 Number of degrees of freedom.
     * @param x                 Atomic coordinates.
     * @param v                 Velocities.
     * @param mass              Mass of each degrees of freedom.
     * @param type              the VARIABLE_TYPE of each variable.
     * @param targetTemperature a double.
     */
    public Thermostat(int n, double[] x, double[] v, double[] mass,
                      VARIABLE_TYPE[] type, double targetTemperature) {
        this(n, x, v, mass, type, targetTemperature, new ArrayList<>());
    }

    public Thermostat(int n, double[] x, double[] v, double[] mass,
                      VARIABLE_TYPE[] type, double targetTemperature, List<Constraint> constraints) {
        assert (n > 3);

        this.nVariables = n;
        this.x = x;
        this.v = v;
        this.mass = mass;
        this.type = type;
        assert (x.length == nVariables);
        assert (v.length == nVariables);
        assert (mass.length == nVariables);
        assert (type.length == nVariables);
        random = new Random();
        setTargetTemperature(targetTemperature);

        this.constraints = new ArrayList<>(constraints);
        // Not every type of constraint constrains just one DoF.
        // SETTLE constraints, for example, constrain three.
        constrainedDoF = constraints.stream().
                mapToInt(Constraint::getNumDegreesFrozen).
                sum();

        // Set the degrees of freedom to nVariables - 3 because we will remove center of mass motion.
        removeCenterOfMassMotion = true;
        dof = nVariables - 3 - constrainedDoF;

        // Update the kinetic energy.
        kineticEnergy();
    }

    /**
     * To allow chemical perturbations during MD.
     *
     * @param n                        Number of degrees of freedom.
     * @param x                        Atomic coordinates.
     * @param v                        Velocities.
     * @param mass                     Mass of each degrees of freedom.
     * @param type                     the VARIABLE_TYPE of each variable.
     * @param removeCenterOfMassMotion a boolean.
     */
    public void setNumberOfVariables(int n, double[] x, double[] v, double[] mass, VARIABLE_TYPE[] type,
                                     boolean removeCenterOfMassMotion) {
        this.nVariables = n;
        this.x = x;
        this.v = v;
        this.mass = mass;
        this.type = type;
        assert (x.length == nVariables);
        assert (v.length == nVariables);
        assert (mass.length == nVariables);
        assert (type.length == nVariables);
        setRemoveCenterOfMassMotion(removeCenterOfMassMotion);
    }

    /**
     * If center of mass motion is being removed, then the mean kinetic energy
     * of the system will be 3 * kT/2 less than if center of mass motion is
     * allowed.
     *
     * @param remove <code>true</code> if center of mass motion is being
     *               removed.
     */
    public void setRemoveCenterOfMassMotion(boolean remove) {
        removeCenterOfMassMotion = remove;
        if (removeCenterOfMassMotion) {
            dof = nVariables - 3 - constrainedDoF;
        } else {
            dof = nVariables - constrainedDoF;
        }
    }

    /**
     * <p>Getter for the field <code>removeCenterOfMassMotion</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getRemoveCenterOfMassMotion() {
        return removeCenterOfMassMotion;
    }

    /**
     * <p>Setter for the field <code>quiet</code>.</p>
     *
     * @param quiet a boolean.
     */
    public void setQuiet(boolean quiet) {
        this.quiet = quiet;
    }

    /**
     * <p>
     * The setRandomSeed method is used to initialize the Random number
     * generator to the same starting state, such that separate runs produce the
     * same Maxwell-Boltzmann initial velocities. same </p>
     *
     * @param seed The seed.
     */
    public void setRandomSeed(long seed) {
        random.setSeed(seed);
    }

    /**
     * <p>
     * Log the target temperature and current number of kT per degree of freedom
     * (should be 0.5 kT at equilibrium).</p>
     *
     * @param level a {@link java.util.logging.Level} object.
     */
    protected void log(Level level) {
        if (logger.isLoggable(level) && !quiet) {
            StringBuilder sb = new StringBuilder("\n");
            sb.append(toString()).append("\n");
            sb.append(format(" Target temperature:           %7.2f Kelvin\n", targetTemperature));
            sb.append(format(" Current temperature:          %7.2f Kelvin\n", currentTemperature));
            sb.append(format(" Number of variables:          %7d\n", nVariables));
            sb.append(format(" Number of degrees of freedom: %7d\n", dof));
            sb.append(format(" Kinetic Energy:               %7.2f\n", currentKineticEnergy));
            sb.append(format(" kT per degree of freedom:     %7.2f\n", KCAL_TO_GRAM_ANG2_PER_PS2 * currentKineticEnergy / (dof * kT)));
            logger.log(level, sb.toString());
        }
    }

    /**
     * <p>
     * getCurrentTemperature</p>
     *
     * @return a double.
     */
    public double getCurrentTemperature() {
        return currentTemperature;
    }

    /**
     * <p>
     * Getter for the field <code>kineticEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getKineticEnergy() {
        return currentKineticEnergy;
    }

    /**
     * <p>
     * Getter for the field <code>targetTemperature</code>.</p>
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
    public void setTargetTemperature(double t) {
        // Obey the Third Law of Thermodynamics.
        assert (t > 0.0);
        targetTemperature = t;
        kT = t * kB;
    }

    /**
     * Return 3 velocities from a Maxwell-Boltzmann distribution of momenta. The
     * variance of each independent momentum component is kT * mass.
     *
     * @param mass The mass for the degrees of freedom.
     * @return three velocity components.
     */
    public double[] maxwellIndividual(double mass) {
        double[] vv = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = random.nextGaussian() * sqrt(kB * targetTemperature / mass);
        }
        return vv;
    }

    /**
     * Reset velocities from a Maxwell-Boltzmann distribution of momenta based
     * on the supplied target temperature. The variance of each independent
     * momentum component is kT * mass.
     *
     * @param targetTemperature the target Temperature for the Maxwell
     *                          distribution.
     */
    public void maxwell(double targetTemperature) {

        setTargetTemperature(targetTemperature);

        for (int i = 0; i < nVariables; i++) {
            double m = mass[i];
            v[i] = random.nextGaussian() * sqrt(kB * targetTemperature / m);
        }

        // Remove the center of mass motion.
        if (removeCenterOfMassMotion) {
            centerOfMassMotion(true, !quiet);
        }

        // Find the current kinetic energy and temperature.
        kineticEnergy();

        /*
          The current temperature will deviate slightly from the target
          temperature if the center of mass motion was removed and/or due to
          finite system size.

          Scale the velocities to enforce the target temperature.
         */
        double scale = sqrt(targetTemperature / currentTemperature);
        for (int i = 0; i < nVariables; i++) {
            v[i] *= scale;
        }

        // Update the kinetic energy and current temperature.
        kineticEnergy();

        log(Level.INFO);
    }

    /**
     * Reset velocities from a Maxwell-Boltzmann distribution based on the
     * current target temperature of thermostat.
     */
    public void maxwell() {
        maxwell(targetTemperature);
    }

    /**
     * Compute the center of mass, linear momentum and angular momentum.
     *
     * @param remove If true, the center of mass motion will be removed.
     * @param print  If true, the center of mass and momenta will be printed.
     */
    public void centerOfMassMotion(boolean remove, boolean print) {
        double totalMass = 0.0;
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
            StringBuilder sb = new StringBuilder(format(" Center of Mass   (%12.3f,%12.3f,%12.3f)\n",
                    centerOfMass[0], centerOfMass[1], centerOfMass[2]));
            sb.append(format(" Linear Momentum  (%12.3f,%12.3f,%12.3f)\n",
                    linearMomentum[0], linearMomentum[1], linearMomentum[2]));
            sb.append(format(" Angular Momentum (%12.3f,%12.3f,%12.3f)",
                    angularMomentum[0], angularMomentum[1], angularMomentum[2]));
            logger.info(sb.toString());
        }

        if (remove) {
            removeCenterOfMassMotion(print);
            centerOfMassMotion(false, print);
        }
    }

    /**
     * Remove center of mass translational and rotational velocity by inverting
     * the moment of inertia tensor.
     *
     * @param print If true, log removal of center of mass motion.
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

        /*
         * RealMatrix inertia = new Array2DRowRealMatrix(3, 3);
         * inertia.setEntry(0, 0, yy + zz); inertia.setEntry(1, 0, -xy);
         * inertia.setEntry(2, 0, -xz); inertia.setEntry(0, 1, -xy);
         * inertia.setEntry(1, 1, xx + zz); inertia.setEntry(2, 1, -yz);
         * inertia.setEntry(0, 2, -xz); inertia.setEntry(1, 2, -yz);
         * inertia.setEntry(2, 2, xx + yy); inertia = new
         * LUDecomposition(inertia).getSolver().getInverse(); xx =
         * inertia.getEntry(0, 0); yy = inertia.getEntry(1, 1); zz =
         * inertia.getEntry(2, 2); xy = inertia.getEntry(0, 1); xz =
         * inertia.getEntry(0, 2); yz = inertia.getEntry(1, 2); double ox =
         * angularMomentum[0] * xx + angularMomentum[1] * xy +
         * angularMomentum[2] * xz; double oy = angularMomentum[0] * xy +
         * angularMomentum[1] * yy + angularMomentum[2] * yz; double oz =
         * angularMomentum[0] * xz + angularMomentum[1] * yz +
         * angularMomentum[2] * zz;
         */

        // Remove center of mass translational momentum.
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

        /*
         * Only remove center of mass rotational momentum for non-periodic
         * systems.
         */
        /*
         * if (false) { index = 0; while (index < nVariables) { if (type[index]
         * == VARIABLE_TYPE.OTHER) { index++; continue; } double xi = x[index++]
         * - centerOfMass[0]; double yi = x[index++] - centerOfMass[1]; double
         * zi = x[index] - centerOfMass[2]; index -= 2; v[index++] += (-oy * zi
         * + oz * yi); v[index++] += (-oz * xi + ox * zi); v[index++] += (-ox *
         * yi + oy * xi); } }
         */
        if (print) {
            logger.info(" Center of mass motion removed.");
        }
    }

    /**
     * Compute the current temperature and kinetic energy of the system.
     */
    public final void kineticEnergy() {
        double e = 0.0;
        for (int i = 0; i < nVariables; i++) {
            double velocity = v[i];
            double v2 = velocity * velocity;
            e += mass[i] * v2;
        }
        currentTemperature = e / (kB * dof);
        e *= 0.5 / KCAL_TO_GRAM_ANG2_PER_PS2;
        currentKineticEnergy = e;
    }

    /**
     * The half-step temperature correction.
     *
     * @param dt a double.
     */
    public abstract void halfStep(double dt);

    /**
     * The full-step temperature correction.
     *
     * @param dt a double.
     */
    public abstract void fullStep(double dt);
}
