// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import ffx.openmm.Context;
import ffx.openmm.Force;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_destroy;
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
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;

/**
 * Weeks-Chandler-Andersen Dispersion Force.
 */
public class WcaDispersionForce extends Force {

  /**
   * Create a new Amoeba WCA dispersion force.
   */
  public WcaDispersionForce() {
    pointer = OpenMM_AmoebaWcaDispersionForce_create();
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
   * Set the water oxygen epsilon parameter.
   *
   * @param epso The water oxygen epsilon parameter.
   */
  public void setEpso(double epso) {
    OpenMM_AmoebaWcaDispersionForce_setEpso(pointer, epso * OpenMM_KJPerKcal);
  }

  /**
   * Set the water hydrogen epsilon parameter.
   *
   * @param epsh The water hydrogen epsilon parameter.
   */
  public void setEpsh(double epsh) {
    OpenMM_AmoebaWcaDispersionForce_setEpsh(pointer, epsh * OpenMM_KJPerKcal);
  }

  /**
   * Set the water oxygen radius parameter.
   *
   * @param rmino The water oxygen radius parameter.
   */
  public void setRmino(double rmino) {
    OpenMM_AmoebaWcaDispersionForce_setRmino(pointer, rmino * OpenMM_NmPerAngstrom);
  }

  /**
   * Set the water hydrogen radius parameter.
   *
   * @param rminh The water hydrogen radius parameter.
   */
  public void setRminh(double rminh) {
    OpenMM_AmoebaWcaDispersionForce_setRminh(pointer, rminh * OpenMM_NmPerAngstrom);
  }

  /**
   * Set the dispersion offset.
   *
   * @param dispoff The dispersion offset.
   */
  public void setDispoff(double dispoff) {
    OpenMM_AmoebaWcaDispersionForce_setDispoff(pointer, dispoff * OpenMM_NmPerAngstrom);
  }

  /**
   * Set the water density parameter.
   *
   * @param awater The water density parameter.
   */
  public void setAwater(double awater) {
    OpenMM_AmoebaWcaDispersionForce_setAwater(pointer, awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
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
   * Set the overlap factor.
   *
   * @param shctd The overlap factor.
   */
  public void setShctd(double shctd) {
    OpenMM_AmoebaWcaDispersionForce_setShctd(pointer, shctd);
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
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaWcaDispersionForce_destroy(pointer);
      pointer = null;
    }
  }
}
