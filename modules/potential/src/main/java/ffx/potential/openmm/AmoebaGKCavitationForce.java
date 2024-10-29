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
package ffx.potential.openmm;

import ffx.openmm.Force;
import ffx.openmm.amoeba.GKCavitationForce;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.implicit.ChandlerCavitation;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.nonbonded.implicit.GaussVol;

import java.util.logging.Level;
import java.util.logging.Logger;

// import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_NonbondedMethod.OpenMM_AmoebaGKCavitationForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static java.lang.String.format;

/**
 * AmoebaCavitationForce.
 */
public class AmoebaGKCavitationForce extends GKCavitationForce {

  private static final Logger logger = Logger.getLogger(AmoebaGKCavitationForce.class.getName());

  /**
   * Constructor.
   *
   * @param openMMEnergy OpenMM energy.
   */
  public AmoebaGKCavitationForce(OpenMMEnergy openMMEnergy) {
    logger.severe(" The AmoebaGKCavitationForce is not currently supported.");
    // TODO: Implement the AmoebaGKCavitationForce as a plugin.

    GeneralizedKirkwood generalizedKirkwood = openMMEnergy.getGK();
    if (generalizedKirkwood == null) {
      destroy();
      return;
    }
    ChandlerCavitation chandlerCavitation = generalizedKirkwood.getChandlerCavitation();
    if (chandlerCavitation == null) {
      destroy();
      return;
    }
    GaussVol gaussVol = chandlerCavitation.getGaussVol();
    if (gaussVol == null) {
      destroy();
      return;
    }

    double surfaceTension = chandlerCavitation.getSurfaceTension()
        * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom / OpenMM_NmPerAngstrom;
    double[] rad = gaussVol.getRadii();

    int index = 0;
    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    for (Atom atom : atoms) {
      int isHydrogen = OpenMM_False;
      double radius = rad[index++];
      if (atom.isHydrogen()) {
        isHydrogen = OpenMM_True;
        radius = 0.0;
      }
      addParticle(radius * OpenMM_NmPerAngstrom, surfaceTension, isHydrogen);
    }

    // TODO: Uncomment this when the AmoebaGKCavitationForce plugin is ready.
    // setNonbondedMethod(OpenMM_AmoebaGKCavitationForce_NoCutoff);

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("GK_FORCE_GROUP", 2);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  GaussVol cavitation force \t\t%d", forceGroup));
  }

  /**
   * Convenience method to construct an AMOEBA Cavitation Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the cavitation information.
   * @return An AMOEBA Cavitation Force, or null if there are no cavitation interactions.
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
    return new AmoebaGKCavitationForce(openMMEnergy);
  }

  /**
   * Update the Cavitation force.
   *
   * @param atoms        The atoms to update.
   * @param openMMEnergy The OpenMM energy term.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood generalizedKirkwood = openMMEnergy.getGK();
    if (generalizedKirkwood == null) {
      return;
    }
    ChandlerCavitation chandlerCavitation = generalizedKirkwood.getChandlerCavitation();
    if (chandlerCavitation == null) {
      return;
    }
    GaussVol gaussVol = chandlerCavitation.getGaussVol();
    if (gaussVol == null) {
      return;
    }

    double surfaceTension = chandlerCavitation.getSurfaceTension()
        * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom / OpenMM_NmPerAngstrom;

    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    double lambdaElec = pme.getAlchemicalParameters().permLambda;

    // Changing cavitation radii is not supported.
    // for (int i=0; i<nAtoms; i++) {
    //  gaussVol.updateAtom(i);
    // }
    double[] rad = gaussVol.getRadii();

    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      double useFactor = 1.0;
      if (!atom.getUse()) {
        useFactor = 0.0;
      }
      // Scale all implicit solvent terms with the square of electrostatics lambda
      // (so dUcav / dL is 0 at lambdaElec = 0).
      double lambdaScale = lambdaElec * lambdaElec;
      if (!atom.applyLambda()) {
        lambdaScale = 1.0;
      }
      useFactor *= lambdaScale;

      double radius = rad[index];
      int isHydrogen = OpenMM_False;
      if (atom.isHydrogen()) {
        isHydrogen = OpenMM_True;
        radius = 0.0;
      }

      setParticleParameters(index, radius * OpenMM_NmPerAngstrom, surfaceTension * useFactor, isHydrogen);
    }
    updateParametersInContext(openMMEnergy.getContext());
  }
}
