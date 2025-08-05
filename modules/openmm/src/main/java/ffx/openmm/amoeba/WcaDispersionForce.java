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
package ffx.openmm.amoeba;

import com.sun.jna.ptr.DoubleByReference;
import ffx.openmm.Context;
import ffx.openmm.Force;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getAwater;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getDispoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getEpsh;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getEpso;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getRminh;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getRmino;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getShctd;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_getSlevy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setAwater;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setDispoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setEpsh;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setEpso;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setRminh;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setRmino;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setShctd;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setSlevy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * Weeks-Chandler-Andersen Dispersion Force.
 */
public class WcaDispersionForce extends Force {

  /**
   * Create a new Amoeba WCA dispersion force.
   */
  public WcaDispersionForce() {
    super(OpenMM_AmoebaWcaDispersionForce_create());
  }

  /**
   * Add a particle to the force field term.
   *
   * @param radius  The radius of the particle.
   * @param epsilon The well depth of the particle.
   */
  public void addParticle(double radius, double epsilon) {
    OpenMM_AmoebaWcaDispersionForce_addParticle(pointer, radius, epsilon);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaWcaDispersionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the water density parameter.
   *
   * @return The water density parameter.
   */
  public double getAwater() {
    return OpenMM_AmoebaWcaDispersionForce_getAwater(pointer);
  }

  /**
   * Get the dispersion offset.
   *
   * @return The dispersion offset.
   */
  public double getDispoff() {
    return OpenMM_AmoebaWcaDispersionForce_getDispoff(pointer);
  }

  /**
   * Get the water hydrogen epsilon parameter.
   *
   * @return The water hydrogen epsilon parameter.
   */
  public double getEpsh() {
    return OpenMM_AmoebaWcaDispersionForce_getEpsh(pointer);
  }

  /**
   * Get the water oxygen epsilon parameter.
   *
   * @return The water oxygen epsilon parameter.
   */
  public double getEpso() {
    return OpenMM_AmoebaWcaDispersionForce_getEpso(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_AmoebaWcaDispersionForce_getNumParticles(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index   The index of the particle.
   * @param radius  The radius of the particle (output).
   * @param epsilon The well depth of the particle (output).
   */
  public void getParticleParameters(int index, DoubleByReference radius, DoubleByReference epsilon) {
    OpenMM_AmoebaWcaDispersionForce_getParticleParameters(pointer, index, radius, epsilon);
  }

  /**
   * Get the water hydrogen radius parameter.
   *
   * @return The water hydrogen radius parameter.
   */
  public double getRminh() {
    return OpenMM_AmoebaWcaDispersionForce_getRminh(pointer);
  }

  /**
   * Get the water oxygen radius parameter.
   *
   * @return The water oxygen radius parameter.
   */
  public double getRmino() {
    return OpenMM_AmoebaWcaDispersionForce_getRmino(pointer);
  }

  /**
   * Get the overlap factor.
   *
   * @return The overlap factor.
   */
  public double getShctd() {
    return OpenMM_AmoebaWcaDispersionForce_getShctd(pointer);
  }

  /**
   * Get the Levy parameter.
   *
   * @return The Levy parameter.
   */
  public double getSlevy() {
    return OpenMM_AmoebaWcaDispersionForce_getSlevy(pointer);
  }

  /**
   * Set the water density parameter.
   *
   * @param awater The water density parameter.
   */
  public void setAwater(double awater) {
    OpenMM_AmoebaWcaDispersionForce_setAwater(pointer, awater);
  }

  /**
   * Set the dispersion offset.
   *
   * @param dispoff The dispersion offset.
   */
  public void setDispoff(double dispoff) {
    OpenMM_AmoebaWcaDispersionForce_setDispoff(pointer, dispoff);
  }

  /**
   * Set the water hydrogen epsilon parameter.
   *
   * @param epsh The water hydrogen epsilon parameter.
   */
  public void setEpsh(double epsh) {
    OpenMM_AmoebaWcaDispersionForce_setEpsh(pointer, epsh);
  }

  /**
   * Set the water oxygen epsilon parameter.
   *
   * @param epso The water oxygen epsilon parameter.
   */
  public void setEpso(double epso) {
    OpenMM_AmoebaWcaDispersionForce_setEpso(pointer, epso);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index   The index of the particle to set.
   * @param radius  The radius of the particle.
   * @param epsilon The well depth of the particle.
   */
  public void setParticleParameters(int index, double radius, double epsilon) {
    OpenMM_AmoebaWcaDispersionForce_setParticleParameters(pointer, index, radius, epsilon);
  }

  /**
   * Set the water hydrogen radius parameter.
   *
   * @param rminh The water hydrogen radius parameter.
   */
  public void setRminh(double rminh) {
    OpenMM_AmoebaWcaDispersionForce_setRminh(pointer, rminh);
  }

  /**
   * Set the water oxygen radius parameter.
   *
   * @param rmino The water oxygen radius parameter.
   */
  public void setRmino(double rmino) {
    OpenMM_AmoebaWcaDispersionForce_setRmino(pointer, rmino);
  }

  /**
   * Set the overlap factor.
   *
   * @param shctd The overlap factor.
   */
  public void setShctd(double shctd) {
    OpenMM_AmoebaWcaDispersionForce_setShctd(pointer, shctd);
  }

  /**
   * Set the Levy parameter.
   *
   * @param slevy The Levy parameter.
   */
  public void setSlevy(double slevy) {
    OpenMM_AmoebaWcaDispersionForce_setSlevy(pointer, slevy);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The OpenMM context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}