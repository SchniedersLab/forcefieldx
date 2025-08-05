// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import ffx.openmm.Integrator;
import ffx.openmm.LangevinIntegrator;
import ffx.openmm.VerletIntegrator;
import ffx.potential.parameters.ForceField;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.logging.Logger;

import static ffx.potential.ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;
import static java.lang.String.format;

/**
 * Create and manage an OpenMM Integrator.
 *
 * <p>An Integrator defines a method for simulating a System by integrating the equations of
 * motion.
 *
 * <p>Each Integrator object is bound to a particular Context which it integrates. This connection
 * is specified by passing the Integrator as an argument to the constructor of the Context.
 */
public class OpenMMIntegrator {

  private static final Logger logger = Logger.getLogger(OpenMMIntegrator.class.getName());

  /**
   * Constraint tolerance as a fraction of the constrained bond length.
   */
  private final static double constraintTolerance = DEFAULT_CONSTRAINT_TOLERANCE;

  /**
   * Prevent instantiation.
   */
  private OpenMMIntegrator() {

  }

  /**
   * Create a integrator.
   *
   * @param name Name of the integrator to use.
   * @param timeStep         Time step (psec).
   * @param temperature      Target temperature (kelvin).
   * @param openMMSystem     OpenMM System.
   * @return Integrator reference.
   */
  public static Integrator createIntegrator(String name, double timeStep, double temperature, OpenMMSystem openMMSystem) {
    switch (name) {
      default -> {
        return createVerletIntegrator(timeStep);
      }
      case "LANGEVIN" -> {
        return createLangevinIntegrator(timeStep, temperature, openMMSystem.getForceField());
      }
      case "MTS" -> {
        return createCustomMTSIntegrator(timeStep, openMMSystem);
      }
      case "LANGEVIN-MTS" -> {
        return createCustomMTSLangevinIntegrator(timeStep, temperature, openMMSystem);
      }
    }
  }

  /**
   * Create a Langevin integrator.
   *
   * @param dt          Time step (psec).
   * @param temperature Temperature (K).
   * @param forceField  Force field.
   */
  public static LangevinIntegrator createLangevinIntegrator(double dt, double temperature, ForceField forceField) {
    CompositeConfiguration properties = forceField.getProperties();
    int seed = 0;
    if (properties.containsKey("integrator-seed")) {
      seed = properties.getInt("integrator-seed", 0);
    }
    double frictionCoeff = forceField.getDouble("FRICTION_COEFF", 91.0);
    LangevinIntegrator langevinIntegrator = new LangevinIntegrator(dt, temperature, frictionCoeff);
    langevinIntegrator.setRandomNumberSeed(seed);
    langevinIntegrator.setConstraintTolerance(constraintTolerance);
    logger.info("  Langevin Integrator");
    logger.info(format("  Target Temperature:   %6.2f (K)", temperature));
    logger.info(format("  Friction Coefficient: %6.2f (1/psec)", frictionCoeff));
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
    return langevinIntegrator;
  }

  /**
   * Create a Custom MTS Integrator.
   *
   * @param dt The outer time step (psec).
   */
  public static CustomMTSIntegrator createCustomMTSIntegrator(double dt, OpenMMSystem openMMSystem) {
    return new CustomMTSIntegrator(dt, constraintTolerance, openMMSystem.hasAmoebaCavitationForce());
  }


  /**
   * Create a Custom MTS Langevin integrator.
   *
   * @param dt           The outer time step (psec).
   * @param temperature  The target temperature (K).
   * @param openMMSystem OpenMM System.
   */
  public static CustomMTSLangevinIntegrator createCustomMTSLangevinIntegrator(double dt, double temperature, OpenMMSystem openMMSystem) {

    ForceField forceField = openMMSystem.getForceField();
    double frictionCoeff = forceField.getDouble("FRICTION_COEFF", 91.0);
    CustomMTSLangevinIntegrator customMTSLangevinIntegrator =
        new CustomMTSLangevinIntegrator(dt, temperature, frictionCoeff, openMMSystem.hasAmoebaCavitationForce());
    return customMTSLangevinIntegrator;
  }

  /**
   * Create a Verlet integrator.
   *
   * @param dt Time step (psec).
   */
  public static VerletIntegrator createVerletIntegrator(double dt) {
    VerletIntegrator verletIntegrator = new VerletIntegrator(dt);
    verletIntegrator.setConstraintTolerance(constraintTolerance);
    logger.info("\n  Verlet Integrator");
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
    return verletIntegrator;
  }

}
