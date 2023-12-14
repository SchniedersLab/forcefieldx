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
package ffx.potential.openmm;

import com.sun.jna.ptr.PointerByReference;
import ffx.potential.parameters.ForceField;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addComputePerDof;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addGlobalVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addPerDofVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addUpdateContextState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VerletIntegrator_create;
import static ffx.utilities.Constants.KCAL_TO_KJ;
import static ffx.utilities.Constants.R;
import static java.lang.Math.exp;
import static java.lang.String.format;
import static java.util.Arrays.copyOfRange;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Create and manage an OpenMM Integrator.
 *
 * <p>An Integrator defines a method for simulating a System by integrating the equations of
 * motion.
 *
 * <p>Each Integrator object is bound to a particular Context which it integrates. This connection
 * is specified by passing the Integrator as an argument to the constructor of the Context.
 */
class OpenMMIntegrator {

  private static final Logger logger = Logger.getLogger(OpenMMIntegrator.class.getName());

  /**
   * Constraint tolerance as a fraction of the constrained bond length.
   */
  private final double constraintTolerance;
  /**
   * Langevin friction coefficient.
   */
  private final double frictionCoefficient;
  /**
   * OpenMM Integrator pointer.
   */
  private PointerByReference integratorPointer = null;

  /**
   * Create an Integrator instance.
   *
   * @param forceField          the ForceField instance containing integrator parameters.
   * @param constraintTolerance The integrator constraint tolerance.
   */
  public OpenMMIntegrator(ForceField forceField, double constraintTolerance) {
    this.constraintTolerance = constraintTolerance;
    frictionCoefficient = forceField.getDouble("FRICTION_COEFF", 91.0);
  }

  /**
   * Use the integrator to step forward.
   *
   * @param steps The number of steps to take.
   */
  public void step(int steps) {
    OpenMM_Integrator_step(integratorPointer, steps);
  }

  /**
   * Return a reference to the integrator.
   *
   * @return Integrator reference.
   */
  public PointerByReference getIntegratorPointer() {
    return integratorPointer;
  }

  /**
   * Create a integrator.
   *
   * @param integratorString Name of the integrator to use.
   * @param timeStep         Time step (psec).
   * @param temperature      Target temperature (kelvin).
   * @param openMMSystem     OpenMM System.
   * @return Integrator reference.
   */
  public PointerByReference createIntegrator(String integratorString, double timeStep, double temperature, OpenMMSystem openMMSystem) {
    switch (integratorString) {
      default -> createVerletIntegrator(timeStep);
      case "LANGEVIN" -> {
        CompositeConfiguration properties = openMMSystem.getForceField().getProperties();
        int seed = 0;
        if (properties.containsKey("integrator-seed")) {
          seed = properties.getInt("integrator-seed", 0);
        }
        createLangevinIntegrator(timeStep, temperature, frictionCoefficient, seed);
      }
      case "MTS" -> createCustomMTSIntegrator(timeStep, openMMSystem);
      case "LANGEVIN-MTS" ->
          createCustomMTSLangevinIntegrator(timeStep, temperature, frictionCoefficient, openMMSystem);
    }

    return integratorPointer;
  }

  /**
   * Create a Langevin integrator.
   *
   * @param dt            Time step (psec).
   * @param temperature   Temperature (K).
   * @param frictionCoeff Frictional coefficient.
   * @param seed          Random number seed.
   */
  public void createLangevinIntegrator(double dt, double temperature, double frictionCoeff, int seed) {
    free();
    integratorPointer = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, dt);
    OpenMM_LangevinIntegrator_setRandomNumberSeed(integratorPointer, seed);
    OpenMM_Integrator_setConstraintTolerance(integratorPointer, constraintTolerance);
    logger.info("  Langevin Integrator");
    logger.info(format("  Target Temperature:   %6.2f (K)", temperature));
    logger.info(format("  Friction Coefficient: %6.2f (1/psec)", frictionCoeff));
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
  }

  /**
   * Create a Custom MTS Integrator.
   *
   * @param dt The outer time step (psec).
   */
  public void createCustomMTSIntegrator(double dt, OpenMMSystem openMMSystem) {
    createCustomIntegrator(dt);

    int n = 4;
    // Force group 1 contains slowly varying forces.
    // Force group 0 contains the fast varying forces.
    int[] forceGroups = {1, 0};
    // There will be 1 force evaluation per outer step, and 4 per inner step.
    int[] subSteps = {1, 4};
    if (openMMSystem.hasAmoebaCavitationForce()) {
      n = 8;
      // Force group 2 contains the cavitation force.
      // Force group 1 contains slowly varying forces.
      // Force group 0 contains the fast varying forces.
      forceGroups = new int[]{2, 1, 0};
      // There will be 1 force evaluation per outer step.
      // There will be 2 force evaluations per middle step.
      // There will be 8 force evaluations per inner step.
      subSteps = new int[]{1, 2, 8};
    }

    OpenMM_CustomIntegrator_addPerDofVariable(integratorPointer, "x1", 0.0);
    OpenMM_CustomIntegrator_addUpdateContextState(integratorPointer);
    createMTSSubStep(1, forceGroups, subSteps);
    OpenMM_CustomIntegrator_addConstrainVelocities(integratorPointer);
    logger.info("  Custom MTS Integrator");
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
    logger.info(format("  Inner Time step:      %6.2f (fsec)", dt / n * 1000));
    logger.info(format("  Friction Coefficient: %6.2f", frictionCoefficient));
  }

  /**
   * Create substeps for the MTS CustomIntegrator.
   *
   * @param parentSubsteps The number of substeps for the previous force group.
   * @param forceGroups    The force groups to be evaluated.
   * @param subSteps       The number of substeps for each force group.
   */
  public void createMTSSubStep(int parentSubsteps, int[] forceGroups, int[] subSteps) {
    int forceGroup = forceGroups[0];
    int steps = subSteps[0];
    int stepsPerParentStep = steps / parentSubsteps;
    if (stepsPerParentStep < 1 || steps % parentSubsteps != 0) {
      throw new IllegalArgumentException(
          "The number for substeps for each group must be a multiple of the number for the previous group");
    }
    if (forceGroup < 0 || forceGroup > 31) {
      throw new IllegalArgumentException("Force group must be between 0 and 31");
    }
    for (int i = 0; i < stepsPerParentStep; i++) {
      OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v",
          "v+0.5*(dt/" + steps + ")*f" + forceGroup + "/m");
      if (forceGroups.length == 1) {
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x", "x+(dt/" + steps + ")*v");
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x1", "x");
        OpenMM_CustomIntegrator_addConstrainPositions(integratorPointer);
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v",
            "v+(x-x1)/(dt/" + steps + ")");
        OpenMM_CustomIntegrator_addConstrainVelocities(integratorPointer);
      } else {
        createMTSSubStep(steps, copyOfRange(forceGroups, 1, forceGroups.length),
            copyOfRange(subSteps, 1, subSteps.length));
      }
      OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v",
          "v+0.5*(dt/" + steps + ")*f" + forceGroup + "/m");
    }
  }

  /**
   * Create a Custom MTS Langevin integrator.
   *
   * @param dt            The outer time step (psec).
   * @param temperature   The target temperature (K).
   * @param frictionCoeff The friction coefficient (1/psec).
   */
  public void createCustomMTSLangevinIntegrator(double dt, double temperature,
                                                double frictionCoeff, OpenMMSystem openMMSystem) {
    createCustomIntegrator(dt);

    int n = 4;
    // Force group 1 contains slowly varying forces.
    // Force group 0 contains the fast varying forces.
    int[] forceGroups = {1, 0};
    // There will be 1 force evaluation per outer step, and 4 per inner step.
    int[] subSteps = {1, 4};
    if (openMMSystem.hasAmoebaCavitationForce()) {
      n = 8;
      // Force group 2 contains the cavitation force.
      // Force group 1 contains slowly varying forces.
      // Force group 0 contains the fast varying forces.
      forceGroups = new int[]{2, 1, 0};
      // There will be 1 force evaluation per outer step.
      // There will be 2 force evaluations per middle step.
      // There will be 8 force evaluations per inner step.
      subSteps = new int[]{1, 2, 8};
    }

    OpenMM_CustomIntegrator_addGlobalVariable(integratorPointer, "a",
        exp(-frictionCoeff * dt / n));
    OpenMM_CustomIntegrator_addGlobalVariable(integratorPointer, "b",
        sqrt(1.0 - exp(-2.0 * frictionCoeff * dt / n)));
    OpenMM_CustomIntegrator_addGlobalVariable(integratorPointer, "kT",
        R * temperature * KCAL_TO_KJ);
    OpenMM_CustomIntegrator_addPerDofVariable(integratorPointer, "x1", 0.0);
    StringBuilder sb = new StringBuilder(" Update Context State\n");
    OpenMM_CustomIntegrator_addUpdateContextState(integratorPointer);

    createMTSLangevinSubStep(1, forceGroups, subSteps, sb);
    // Log the substeps.
    logger.finest(" Langevin-MTS steps:" + sb);
    OpenMM_CustomIntegrator_addConstrainVelocities(integratorPointer);
    logger.info("  Custom MTS Langevin Integrator");
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
    logger.info(format("  Inner Time step:      %6.2f (fsec)", dt / n * 1000));
    logger.info(format("  Friction Coefficient: %6.2f (1/psec)", frictionCoeff));
  }

  /**
   * Create substeps for the MTS Langevin CustomIntegrator.
   *
   * @param parentSubsteps The number of substeps for the previous force group.
   * @param forceGroups    The force groups to be evaluated.
   * @param subSteps       The number of substeps for each force group.
   */
  public void createMTSLangevinSubStep(int parentSubsteps, int[] forceGroups, int[] subSteps,
                                       StringBuilder sb) {
    int forceGroup = forceGroups[0];
    int steps = subSteps[0];
    int stepsPerParentStep = steps / parentSubsteps;
    if (stepsPerParentStep < 1 || steps % parentSubsteps != 0) {
      throw new IllegalArgumentException(
          "The number for substeps for each group must be a multiple of the number for the previous group");
    }
    if (forceGroup < 0 || forceGroup > 31) {
      throw new IllegalArgumentException("Force group must be between 0 and 31");
    }

    sb.append(" Force Group: ").append(forceGroup).append(" ForceGroup length: ")
        .append(forceGroups.length).append(" Steps: ").append(steps)
        .append(" Step Per Parent Step: ").append(stepsPerParentStep).append(" Parent Sub Steps: ")
        .append(parentSubsteps).append("\n");

    for (int i = 0; i < stepsPerParentStep; i++) {
      String step = "v+0.5*(dt/" + steps + ")*f" + forceGroup + "/m";
      sb.append(" v = ").append(step).append("\n");
      OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", step);
      // String step;
      if (forceGroups.length == 1) {
        step = "x+(dt/" + 2 * steps + ")*v";
        sb.append(" x = ").append(step).append("\n");
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x", step);
        step = "a*v + b*sqrt(kT/m)*gaussian";
        sb.append(" v = ").append(step).append("\n");
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", step);
        step = "x+(dt/" + 2 * steps + ")*v";
        sb.append(" x = ").append(step).append("\n");
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x", step);
        step = "x";
        sb.append(" x1 = ").append(step).append("\n");
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x1", step);
        sb.append(" Constrain Positions\n");
        OpenMM_CustomIntegrator_addConstrainPositions(integratorPointer);
        step = "v+(x-x1)/(dt/" + steps + ")";
        sb.append(" v = ").append(step).append("\n");
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", step);
        sb.append(" Constrain Velocities\n");
        OpenMM_CustomIntegrator_addConstrainVelocities(integratorPointer);
      } else {
        createMTSLangevinSubStep(steps, copyOfRange(forceGroups, 1, forceGroups.length),
            copyOfRange(subSteps, 1, subSteps.length), sb);
      }
      step = "v+0.5*(dt/" + steps + ")*f" + forceGroup + "/m";
      sb.append(" v = ").append(step).append("\n");
      OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", step);
    }
  }

  /**
   * Create a Verlet integrator.
   *
   * @param dt Time step (psec).
   */
  public void createVerletIntegrator(double dt) {
    free();
    integratorPointer = OpenMM_VerletIntegrator_create(dt);
    OpenMM_Integrator_setConstraintTolerance(integratorPointer, constraintTolerance);
    logger.info("  Verlet Integrator");
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
  }

  /**
   * Create a Custom integrator.
   *
   * @param dt Time step (psec).
   */
  public void createCustomIntegrator(double dt) {
    free();
    integratorPointer = OpenMM_CustomIntegrator_create(dt);
    OpenMM_Integrator_setConstraintTolerance(integratorPointer, constraintTolerance);
  }

  /**
   * Destroy the integrator instance.
   */
  public void free() {
    if (integratorPointer != null) {
      logger.fine(" Free OpenMM Integrator.");
      OpenMM_Integrator_destroy(integratorPointer);
      logger.fine(" Free OpenMM Integrator completed.");
      integratorPointer = null;
    }
  }
}
