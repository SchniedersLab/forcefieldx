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

import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.parameters.VDWType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_create;
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
import static java.lang.String.format;

public class AmoebaWcaDispersionForce extends OpenMMForce {

  private static final Logger logger = Logger.getLogger(AmoebaGeneralizedKirkwoodForce.class.getName());

  /**
   * Create a new Amoeba WCA dispersion force.
   *
   * @param openMMEnergy The OpenMM energy term.
   */
  public AmoebaWcaDispersionForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      return;
    }
    DispersionRegion dispersionRegion = gk.getDispersionRegion();
    if (dispersionRegion == null) {
      return;
    }

    double epso = 0.1100;
    double epsh = 0.0135;
    double rmino = 1.7025;
    double rminh = 1.3275;
    double awater = 0.033428;
    double slevy = 1.0;
    double dispoff = dispersionRegion.getDispersionOffset();
    double shctd = dispersionRegion.getDispersionOverlapFactor();

    VanDerWaals vdW = openMMEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    double radScale = 1.0;
    if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    forcePointer = OpenMM_AmoebaWcaDispersionForce_create();

    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    for (Atom atom : atoms) {
      VDWType vdwType = atom.getVDWType();
      double radius = vdwType.radius;
      double eps = vdwType.wellDepth;
      addParticle(OpenMM_NmPerAngstrom * radius * radScale, OpenMM_KJPerKcal * eps);
    }

    setEpso(epso * OpenMM_KJPerKcal);
    setEpsh(epsh * OpenMM_KJPerKcal);
    setRmino(rmino * OpenMM_NmPerAngstrom);
    setRminh(rminh * OpenMM_NmPerAngstrom);
    setDispoff(dispoff * OpenMM_NmPerAngstrom);
    setAwater(awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
    setSlevy(slevy);
    setShctd(shctd);

    int forceGroup = openMMEnergy.getSystem().getForceField().getInteger("GK_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  WCA dispersion force \t\t\t%d", forceGroup));
  }

  /**
   * Convenience method to construct an AMOEBA WCA Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the WCA information.
   * @return An AMOEBA WCA Force, or null if there are no WCA interactions.
   */
  public static OpenMMForce constructForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      return null;
    }
    DispersionRegion dispersionRegion = gk.getDispersionRegion();
    if (dispersionRegion == null) {
      return null;
    }
    return new AmoebaWcaDispersionForce(openMMEnergy);
  }

  /**
   * Update the WCA force.
   *
   * @param atoms        The atoms to update.
   * @param openMMEnergy The OpenMM energy term.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    double radScale = 1.0;
    if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    double lambdaElec = openMMEnergy.getSystem().getLambdaElec();

    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      double useFactor = 1.0;
      if (!atom.getUse()) {
        useFactor = 0.0;
      }

      // Scale all implicit solvent terms with the square of electrostatics lambda
      // (so dUdisp / dL is 0 at lambdaElec = 0).
      double lambdaScale = lambdaElec * lambdaElec;
      if (!atom.applyLambda()) {
        lambdaScale = 1.0;
      }
      useFactor *= lambdaScale;

      VDWType vdwType = atom.getVDWType();
      double radius = vdwType.radius;
      double eps = vdwType.wellDepth;
      setParticleParameters(index, OpenMM_NmPerAngstrom * radius * radScale, OpenMM_KJPerKcal * eps * useFactor);
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Add a particle to the force field term.
   *
   * @param radius  The radius of the particle.
   * @param epsilon The well depth of the particle.
   */
  public void addParticle(double radius, double epsilon) {
    OpenMM_AmoebaWcaDispersionForce_addParticle(forcePointer, radius, epsilon);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index   The index of the particle to set.
   * @param radius  The radius of the particle.
   * @param epsilon The well depth of the particle.
   */
  public void setParticleParameters(int index, double radius, double epsilon) {
    OpenMM_AmoebaWcaDispersionForce_setParticleParameters(forcePointer, index, radius, epsilon);
  }


  /**
   * Set the water oxygen epsilon parameter.
   *
   * @param epso The water oxygen epsilon parameter.
   */
  public void setEpso(double epso) {
    OpenMM_AmoebaWcaDispersionForce_setEpso(forcePointer, epso * OpenMM_KJPerKcal);
  }

  /**
   * Set the water hydrogen epsilon parameter.
   *
   * @param epsh The water hydrogen epsilon parameter.
   */
  public void setEpsh(double epsh) {
    OpenMM_AmoebaWcaDispersionForce_setEpsh(forcePointer, epsh * OpenMM_KJPerKcal);
  }

  /**
   * Set the water oxygen radius parameter.
   *
   * @param rmino The water oxygen radius parameter.
   */
  public void setRmino(double rmino) {
    OpenMM_AmoebaWcaDispersionForce_setRmino(forcePointer, rmino * OpenMM_NmPerAngstrom);
  }

  /**
   * Set the water hydrogen radius parameter.
   *
   * @param rminh The water hydrogen radius parameter.
   */
  public void setRminh(double rminh) {
    OpenMM_AmoebaWcaDispersionForce_setRminh(forcePointer, rminh * OpenMM_NmPerAngstrom);
  }

  /**
   * Set the dispersion offset.
   *
   * @param dispoff The dispersion offset.
   */
  public void setDispoff(double dispoff) {
    OpenMM_AmoebaWcaDispersionForce_setDispoff(forcePointer, dispoff * OpenMM_NmPerAngstrom);
  }

  /**
   * Set the water density parameter.
   *
   * @param awater The water density parameter.
   */
  public void setAwater(double awater) {
    OpenMM_AmoebaWcaDispersionForce_setAwater(forcePointer, awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
  }

  /**
   * Set the Levy parameter.
   *
   * @param slevy The Levy parameter.
   */
  public void setSlevy(double slevy) {
    OpenMM_AmoebaWcaDispersionForce_setSlevy(forcePointer, slevy);
  }

  /**
   * Set the overlap factor.
   *
   * @param shctd The overlap factor.
   */
  public void setShctd(double shctd) {
    OpenMM_AmoebaWcaDispersionForce_setShctd(forcePointer, shctd);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The OpenMM context.
   */
  public void updateParametersInContext(OpenMMContext context) {
    if (context.hasContextPointer()) {
      OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(forcePointer, context.getContextPointer());
    }
  }
}
