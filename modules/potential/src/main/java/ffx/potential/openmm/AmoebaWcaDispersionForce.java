// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import ffx.openmm.Force;
import ffx.openmm.amoeba.WcaDispersionForce;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.parameters.VDWType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;

public class AmoebaWcaDispersionForce extends WcaDispersionForce {

  private static final Logger logger = Logger.getLogger(AmoebaGeneralizedKirkwoodForce.class.getName());

  /**
   * Create a new Amoeba WCA dispersion force.
   *
   * @param openMMEnergy The OpenMM energy term.
   */
  public AmoebaWcaDispersionForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      destroy();
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
    if (vdwForm.radiusSize == VDWType.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

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

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("GK_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  WCA dispersion force \t\t\t%d", forceGroup));
  }

  /**
   * Convenience method to construct an AMOEBA WCA Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the WCA information.
   * @return An AMOEBA WCA Force, or null if there are no WCA interactions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
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
    if (vdwForm.radiusSize == VDWType.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    double lambdaElec = openMMEnergy.getPmeNode().getAlchemicalParameters().permLambda;

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
}
