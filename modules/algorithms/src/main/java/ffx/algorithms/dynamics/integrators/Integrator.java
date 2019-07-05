//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.dynamics.integrators;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import static java.lang.System.arraycopy;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.numerics.Constraint;

/**
 * The Integrator class is responsible for propagation of degrees of freedom
 * through time. Implementations must define their behavior at pre-force and
 * post-force evaluation time points.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class Integrator {

    private static final Logger logger = Logger.getLogger(Integrator.class.getName());

    /**
     * Parse an integrator String into an instance of the IntegratorEnum enum.
     *
     * @param str Integrator string.
     * @return Integrator enum.
     */
    public static IntegratorEnum parseIntegrator(String str) {
        try {
            return IntegratorEnum.valueOf(str.toUpperCase().replaceAll("\\s+", ""));
        } catch (Exception e) {
            logger.info(String.format(" Could not parse %s as an integrator; defaulting to Verlet.", str));
            return IntegratorEnum.VERLET;
        }
    }

    /**
     * Number of variables.
     */
    protected int nVariables;
    /**
     * Mass of each degree of freedom.
     */
    protected double[] mass;
    /**
     * Coordinates for each degree of freedom in Angstroms.
     */
    protected double[] x;
    /**
     * Velocity of each degree of freedom in Angstroms per picosecond.
     */
    protected double[] v;
    /**
     * Acceleration of each degree of freedom in Angstroms per square picosecond.
     */
    protected double[] a;
    /**
     * Acceleration of each degree of freedom for the previous step.
     */
    double[] aPrevious;

    /**
     * Time step (psec).
     */
    protected double dt;
    /**
     * Half the time step (psec).
     */
    double dt_2;
    /**
     * Any geometric constraints to apply during integration.
     */
    protected List<Constraint> constraints = new ArrayList<>();
    /**
     * If there are constraints present.
     */
    protected boolean useConstraints = false;
    /**
     * Numerical tolerance (as a fraction of bond length) permitted for numerical solutions to constraints.
     */
    protected final double constraintTolerance = ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;

    /**
     * Constructor for Integrator.
     *
     * @param nVariables number of Variables.
     * @param x          Cartesian coordinates (Angstroms).
     * @param v          Velocities.
     * @param a          Accelerations.
     * @param aPrevious  Previous Accelerations.
     * @param mass       Mass.
     */
    public Integrator(int nVariables, double[] x, double[] v, double[] a,
                      double[] aPrevious, double[] mass) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.aPrevious = aPrevious;
        this.mass = mass;
        dt = 1.0e-3;
        dt_2 = dt / 2.0;
    }

    /**
     * Constructor for Integrator that do not use previous accelerations.
     *
     * @param nVariables number of Variables.
     * @param x          Cartesian coordinates (Angstroms).
     * @param v          Velocities.
     * @param a          Accelerations.
     * @param mass       Mass.
     */
    public Integrator(int nVariables, double[] x, double[] v, double[] a, double[] mass) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.mass = mass;
        this.aPrevious = new double[nVariables];
        dt = 1.0e-3;
    }

    /**
     * Update the integrator to be consistent with chemical perturbations.
     *
     * @param nVariables the number of variables being integrated.
     * @param x          the current value of each variable.
     * @param v          the current velocity of each variable.
     * @param a          the current acceleration of each variable.
     * @param aPrevious  the previous acceleration of each variable.
     * @param mass       the mass for each variable.
     */
    public void setNumberOfVariables(int nVariables, double[] x, double[] v,
                                     double[] a, double[] aPrevious, double[] mass) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.aPrevious = aPrevious;
        this.mass = mass;
    }

    /**
     * Update the integrator to be consistent with chemical perturbations.
     *
     * @param nVariables the number of variables being integrated.
     * @param x          the current value of each variable.
     * @param v          the current velocity of each variable.
     * @param a          the current acceleration of each variable.
     * @param mass       the mass for each variable.
     */
    public void setNumberOfVariables(int nVariables, double[] x, double[] v,
                                     double[] a, double[] mass) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.mass = mass;
    }

    /**
     * Copy acceleration to previous acceleration.
     */
    public void copyAccelerationToPrevious() {
        if (aPrevious == null || aPrevious.length < a.length) {
            aPrevious = new double[a.length];
        }
        arraycopy(a, 0, aPrevious, 0, nVariables);
    }

    /**
     * Get the time step.
     *
     * @return the time step (psec).
     */
    public double getTimeStep() {
        return dt;
    }

    /**
     * Adds a set of Constraints that this Integrator must respect.
     *
     * @param addedConstraints Constraints to add.
     */
    public void addConstraints(List<Constraint> addedConstraints) {
        constraints.addAll(addedConstraints);
        useConstraints = true;
    }

    /**
     * Returns a copy of the list of Constraints.
     *
     * @return All Constraints this Integrator respects.
     */
    public List<Constraint> getConstraints() {
        return new ArrayList<>(constraints);
    }

    public void removeConstraint(Constraint constraint) {
        constraints.remove(constraint);
        useConstraints = !constraints.isEmpty();
    }

    public void removeConstraints(Collection<Constraint> toRemove) {
        constraints = constraints.stream().
                filter((Constraint c) -> !toRemove.contains(c)).
                collect(Collectors.toList());
        useConstraints = !constraints.isEmpty();
    }

    /**
     * Set the time step.
     *
     * @param dt the time step (psec).
     */
    abstract public void setTimeStep(double dt);

    /**
     * Integrator pre-force evaluation operation.
     *
     * @param potential the Potential this integrator operates on.
     */
    abstract public void preForce(Potential potential);

    /**
     * Integrator post-force evaluation operation.
     *
     * @param gradient the gradient for the post-force operation.
     */
    abstract public void postForce(double[] gradient);

}
