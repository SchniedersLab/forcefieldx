// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
// ******************************************************************************
package ffx.algorithms.dynamics.integrators;

import ffx.potential.SystemState;
import ffx.numerics.Constraint;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * The Integrator class is responsible for propagation of degrees of freedom through time.
 * Implementations must define their behavior at pre-force and post-force evaluation time points.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class Integrator {

  private static final Logger logger = Logger.getLogger(Integrator.class.getName());
  /**
   * Numerical tolerance (as a fraction of bond length) permitted for numerical solutions to
   * constraints.
   */
  protected final double constraintTolerance = ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;

  /**
   * The MDState class contains the current state of the Molecular Dynamics simulation.
   */
  protected final SystemState state;
  /** Time step (psec). */
  protected double dt;
  /** Any geometric constraints to apply during integration. */
  protected List<Constraint> constraints = new ArrayList<>();
  /** If there are constraints present. */
  protected boolean useConstraints = false;
  /** Half the time step (psec). */
  double dt_2;

  /**
   * Constructor for Integrator.
   *
   * @param state The MD state to operate on.
   */
  public Integrator(SystemState state) {
    this.state = state;
    dt = 1.0e-3;
    dt_2 = dt / 2.0;
  }


  /**
   * Parse an integrator String into an instance of the IntegratorEnum enum.
   *
   * @param str Integrator string.
   * @return Integrator enum.
   */
  public static IntegratorEnum parseIntegrator(String str) {
    try {
      String integrator = str.toUpperCase().replaceAll("\\s+", "");
      integrator = integrator.replaceAll("-", "_");
      return IntegratorEnum.valueOf(integrator);
    } catch (Exception e) {
      logger.info(String.format(" Could not parse %s as an integrator; defaulting to Verlet.", str));
      return IntegratorEnum.VERLET;
    }
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

  /** Copy acceleration to previous acceleration. */
  public void copyAccelerationToPrevious() {
    state.copyAccelerationsToPrevious();
  }

  /**
   * Returns a copy of the list of Constraints.
   *
   * @return All Constraints this Integrator respects.
   */
  public List<Constraint> getConstraints() {
    return new ArrayList<>(constraints);
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
   * Set the time step.
   *
   * @param dt the time step (psec).
   */
  public abstract void setTimeStep(double dt);

  /**
   * Integrator post-force evaluation operation.
   *
   * @param gradient the gradient for the post-force operation.
   */
  public abstract void postForce(double[] gradient);

  /**
   * Integrator pre-force evaluation operation.
   *
   * @param potential the Potential this integrator operates on.
   */
  public abstract void preForce(Potential potential);

  public void removeConstraint(Constraint constraint) {
    constraints.remove(constraint);
    useConstraints = !constraints.isEmpty();
  }

  public void removeConstraints(Collection<Constraint> toRemove) {
    constraints = constraints.stream().filter((Constraint c) -> !toRemove.contains(c))
        .collect(Collectors.toList());
    useConstraints = !constraints.isEmpty();
  }

}
