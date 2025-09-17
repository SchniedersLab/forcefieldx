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
package ffx.openmm;

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_addException;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getNumExceptions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_getUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_setExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_setUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_GayBerneForce_usesPeriodicBoundaryConditions;

/**
 * This class implements the Gay-Berne anisotropic potential.  This is similar to a Lennard-Jones potential,
 * but it represents the particles as ellipsoids rather than point particles.  In addition to the standard
 * sigma and epsilon parameters, each particle has three widths sx, sy, and sz that give the diameter of the
 * ellipsoid along each axis.  It also has three scale factors ex, ey, and ez that scale the strength
 * of the interaction along each axis.  You can think of this force as a Lennard-Jones interaction computed
 * based on the distance between the nearest points on two ellipsoids.  The scale factors act as multipliers
 * for epsilon along each axis, so the strength of the interaction along the ellipsoid's x axis is multiplied by
 * ex, and likewise for the other axes.  If two particles each have all their widths set to sigma and all their
 * scale factors set to 1, the interaction simplifies to a standard Lennard-Jones force between point particles.
 * <p>
 * The orientation of a particle's ellipsoid is determined based on the positions of two other particles.
 * The vector to the first particle sets the direction of the x axis.  The vector to the second particle
 * (after subtracting out any x component) sets the direction of the y axis.  If the ellipsoid is axially
 * symmetric (sy=sz and ey=ez), you can omit the second particle and define only an x axis direction.
 * If the ellipsoid is a sphere (all three widths and all three scale factors are equal), both particles
 * can be omitted.
 * <p>
 * To determine the values of sigma and epsilon for an interaction, this class uses Lorentz-Berthelot
 * combining rules: it takes the arithmetic mean of the sigmas and the geometric mean of the epsilons for
 * the two interacting particles.  You also can specify "exceptions", particular pairs of particles for
 * which different values should be used.
 * <p>
 * To use this class, create a GayBerneForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define parameters must be exactly
 * equal to the number of particles in the System, or else an exception will be thrown when you try to
 * create a Context.  After a particle has been added, you can modify its force field parameters by calling
 * setParticleParameters().  This will have no effect on Contexts that already exist unless you call
 * updateParametersInContext().
 * <p>
 * When using a cutoff, by default interactions are sharply truncated at the cutoff distance.  Optionally
 * you can instead use a switching function to make the interaction smoothly go to zero over a finite
 * distance range.  To enable this, call setUseSwitchingFunction().  You must also call setSwitchingDistance()
 * to specify the distance at which the interaction should begin to decrease.  The switching distance must be
 * less than the cutoff distance.
 */
public class GayBerneForce extends Force {

  /**
   * Create a new GayBerneForce.
   */
  public GayBerneForce() {
    super(OpenMM_GayBerneForce_create());
  }

  /**
   * Add an exception to the force.
   *
   * @param particle1 The index of the first particle.
   * @param particle2 The index of the second particle.
   * @param sigma     The sigma parameter for the exception.
   * @param epsilon   The epsilon parameter for the exception.
   * @param replace   Whether to replace an existing exception.
   * @return The index of the exception that was added.
   */
  public int addException(int particle1, int particle2, double sigma, double epsilon, int replace) {
    return OpenMM_GayBerneForce_addException(pointer, particle1, particle2, sigma, epsilon, replace);
  }

  /**
   * Add a particle to the force.
   *
   * @param sigma     The sigma parameter.
   * @param epsilon   The epsilon parameter.
   * @param xparticle The x-axis particle type.
   * @param yparticle The y-axis particle type.
   * @param ex        The x-axis shape parameter.
   * @param ey        The y-axis shape parameter.
   * @param ez        The z-axis shape parameter.
   * @param sx        The x-axis strength parameter.
   * @param sy        The y-axis strength parameter.
   * @param sz        The z-axis strength parameter.
   * @return The index of the particle that was added.
   */
  public int addParticle(double sigma, double epsilon, int xparticle, int yparticle,
                         double ex, double ey, double ez, double sx, double sy, double sz) {
    return OpenMM_GayBerneForce_addParticle(pointer, sigma, epsilon, xparticle, yparticle,
        ex, ey, ez, sx, sy, sz);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_GayBerneForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance, measured in nm.
   */
  public double getCutoffDistance() {
    return OpenMM_GayBerneForce_getCutoffDistance(pointer);
  }

  /**
   * Get the parameters for an exception.
   *
   * @param index     The index of the exception.
   * @param particle1 The index of the first particle (output).
   * @param particle2 The index of the second particle (output).
   * @param sigma     The sigma parameter for the exception (output).
   * @param epsilon   The epsilon parameter for the exception (output).
   */
  public void getExceptionParameters(int index, IntByReference particle1, IntByReference particle2,
                                     DoubleByReference sigma, DoubleByReference epsilon) {
    OpenMM_GayBerneForce_getExceptionParameters(pointer, index, particle1, particle2, sigma, epsilon);
  }

  /**
   * Get the parameters for an exception.
   *
   * @param index     The index of the exception.
   * @param particle1 The index of the first particle (output).
   * @param particle2 The index of the second particle (output).
   * @param sigma     The sigma parameter for the exception (output).
   * @param epsilon   The epsilon parameter for the exception (output).
   */
  public void getExceptionParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                     DoubleBuffer sigma, DoubleBuffer epsilon) {
    OpenMM_GayBerneForce_getExceptionParameters(pointer, index, particle1, particle2, sigma, epsilon);
  }

  /**
   * Get the nonbonded method.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_GayBerneForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of exceptions.
   *
   * @return The number of exceptions.
   */
  public int getNumExceptions() {
    return OpenMM_GayBerneForce_getNumExceptions(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_GayBerneForce_getNumParticles(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index     The index of the particle.
   * @param sigma     The sigma parameter (output).
   * @param epsilon   The epsilon parameter (output).
   * @param xparticle The x-axis particle type (output).
   * @param yparticle The y-axis particle type (output).
   * @param ex        The x-axis shape parameter (output).
   * @param ey        The y-axis shape parameter (output).
   * @param ez        The z-axis shape parameter (output).
   * @param sx        The x-axis strength parameter (output).
   * @param sy        The y-axis strength parameter (output).
   * @param sz        The z-axis strength parameter (output).
   */
  public void getParticleParameters(int index, DoubleByReference sigma, DoubleByReference epsilon,
                                    IntByReference xparticle, IntByReference yparticle,
                                    DoubleByReference ex, DoubleByReference ey, DoubleByReference ez,
                                    DoubleByReference sx, DoubleByReference sy, DoubleByReference sz) {
    OpenMM_GayBerneForce_getParticleParameters(pointer, index, sigma, epsilon, xparticle, yparticle,
        ex, ey, ez, sx, sy, sz);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index     The index of the particle.
   * @param sigma     The sigma parameter (output).
   * @param epsilon   The epsilon parameter (output).
   * @param xparticle The x-axis particle type (output).
   * @param yparticle The y-axis particle type (output).
   * @param ex        The x-axis shape parameter (output).
   * @param ey        The y-axis shape parameter (output).
   * @param ez        The z-axis shape parameter (output).
   * @param sx        The x-axis strength parameter (output).
   * @param sy        The y-axis strength parameter (output).
   * @param sz        The z-axis strength parameter (output).
   */
  public void getParticleParameters(int index, DoubleBuffer sigma, DoubleBuffer epsilon,
                                    IntBuffer xparticle, IntBuffer yparticle,
                                    DoubleBuffer ex, DoubleBuffer ey, DoubleBuffer ez,
                                    DoubleBuffer sx, DoubleBuffer sy, DoubleBuffer sz) {
    OpenMM_GayBerneForce_getParticleParameters(pointer, index, sigma, epsilon, xparticle, yparticle,
        ex, ey, ez, sx, sy, sz);
  }

  /**
   * Get the switching distance.
   *
   * @return The switching distance, measured in nm.
   */
  public double getSwitchingDistance() {
    return OpenMM_GayBerneForce_getSwitchingDistance(pointer);
  }

  /**
   * Get whether a switching function is used.
   *
   * @return 1 if a switching function is used, 0 otherwise.
   */
  public int getUseSwitchingFunction() {
    return OpenMM_GayBerneForce_getUseSwitchingFunction(pointer);
  }

  /**
   * Set the cutoff distance.
   *
   * @param distance The cutoff distance, measured in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_GayBerneForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the parameters for an exception.
   *
   * @param index     The index of the exception.
   * @param particle1 The index of the first particle.
   * @param particle2 The index of the second particle.
   * @param sigma     The sigma parameter for the exception.
   * @param epsilon   The epsilon parameter for the exception.
   */
  public void setExceptionParameters(int index, int particle1, int particle2, double sigma, double epsilon) {
    OpenMM_GayBerneForce_setExceptionParameters(pointer, index, particle1, particle2, sigma, epsilon);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_GayBerneForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index     The index of the particle.
   * @param sigma     The sigma parameter.
   * @param epsilon   The epsilon parameter.
   * @param xparticle The x-axis particle type.
   * @param yparticle The y-axis particle type.
   * @param ex        The x-axis shape parameter.
   * @param ey        The y-axis shape parameter.
   * @param ez        The z-axis shape parameter.
   * @param sx        The x-axis strength parameter.
   * @param sy        The y-axis strength parameter.
   * @param sz        The z-axis strength parameter.
   */
  public void setParticleParameters(int index, double sigma, double epsilon, int xparticle, int yparticle,
                                    double ex, double ey, double ez, double sx, double sy, double sz) {
    OpenMM_GayBerneForce_setParticleParameters(pointer, index, sigma, epsilon, xparticle, yparticle,
        ex, ey, ez, sx, sy, sz);
  }

  /**
   * Set the switching distance.
   *
   * @param distance The switching distance, measured in nm.
   */
  public void setSwitchingDistance(double distance) {
    OpenMM_GayBerneForce_setSwitchingDistance(pointer, distance);
  }

  /**
   * Set whether to use a switching function.
   *
   * @param use 1 to use a switching function, 0 otherwise.
   */
  public void setUseSwitchingFunction(int use) {
    OpenMM_GayBerneForce_setUseSwitchingFunction(pointer, use);
  }

  /**
   * Update the parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_GayBerneForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_GayBerneForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}