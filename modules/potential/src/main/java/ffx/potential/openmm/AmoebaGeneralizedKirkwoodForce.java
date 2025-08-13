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

import ffx.openmm.Force;
import ffx.openmm.amoeba.GeneralizedKirkwoodForce;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static java.lang.String.format;

public class AmoebaGeneralizedKirkwoodForce extends GeneralizedKirkwoodForce {

  private static final Logger logger = Logger.getLogger(AmoebaGeneralizedKirkwoodForce.class.getName());

  public AmoebaGeneralizedKirkwoodForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      destroy();
      return;
    }

    setSolventDielectric(gk.getSolventPermittivity());
    double soluteDielectric = gk.getSolutePermittivity();
    if (soluteDielectric != 1.0) {
      logger.severe(" Solute dielectric is not 1.0, which is not supported by OpenMM.");
    }
    setSoluteDielectric(1.0);
    setDielectricOffset(gk.getDescreenOffset() * OpenMM_NmPerAngstrom);


    boolean usePerfectRadii = gk.getUsePerfectRadii();
    double perfectRadiiScale = 1.0;
    if (usePerfectRadii) {
      // No descreening when using perfect radii (OpenMM will just load the base radii).
      perfectRadiiScale = 0.0;
    }

    // Turn on tanh rescaling only when not using perfect radii.
    int tanhRescale = 0;
    if (gk.getTanhCorrection() && !usePerfectRadii) {
      tanhRescale = 1;
    }
    double[] betas = gk.getTanhBetas();
    setTanhRescaling(tanhRescale);
    setTanhParameters(betas[0], betas[1], betas[2]);

    double[] baseRadius = gk.getBaseRadii();
    if (usePerfectRadii) {
      baseRadius = gk.getPerfectRadii();
    }

    double[] overlapScale = gk.getOverlapScale();
    double[] descreenRadius = gk.getDescreenRadii();
    double[] neckFactor = gk.getNeckScale();

    if (!usePerfectRadii && logger.isLoggable(Level.FINE)) {
      logger.fine("   GK Base Radii  Descreen Radius  Overlap Scale  Overlap");
    }

    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      MultipoleType multipoleType = atoms[i].getMultipoleType();
      double base = baseRadius[i] * OpenMM_NmPerAngstrom;
      double descreen = descreenRadius[i] * OpenMM_NmPerAngstrom * perfectRadiiScale;
      double overlap = overlapScale[i] * perfectRadiiScale;
      double neck = neckFactor[i] * perfectRadiiScale;
      addParticle(multipoleType.charge, base, overlap, descreen, neck);
      if (!usePerfectRadii && logger.isLoggable(Level.FINE)) {
        logger.fine(format("   %s %8.6f %8.6f %5.3f", atoms[i].toString(), baseRadius[i], descreenRadius[i], overlapScale[i]));
      }
    }

    setProbeRadius(gk.getProbeRadius() * OpenMM_NmPerAngstrom);

    GeneralizedKirkwood.NonPolarModel nonpolar = gk.getNonPolarModel();
    switch (nonpolar) {
      case BORN_CAV_DISP, BORN_SOLV -> {
        // Configure a Born Radii based surface area term.
        double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm;
        setIncludeCavityTerm(OpenMM_True);
        setSurfaceAreaFactor(-surfaceTension);
      }
      // Other models do not use a Born Radii based surface area term.
      default -> setIncludeCavityTerm(OpenMM_False);
    }

    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Generalized Kirkwood force \t\t%d", forceGroup));
  }

  /**
   * Convenience method to construct an AMOEBA Generalized Kirkwood Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the GK information.
   * @return An AMOEBA GK Force, or null if there are no GK interactions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      return null;
    }
    return new AmoebaGeneralizedKirkwoodForce(openMMEnergy);
  }


  /**
   * Update the force.
   *
   * @param atoms        The atoms to update.
   * @param openMMEnergy The OpenMM energy.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null || pointer == null) {
      return;
    }

    // Update the GK solute parameters.
    int nAtoms = openMMEnergy.getMolecularAssembly().getAtomArray().length;
    for (int i = 0; i < nAtoms; i++) {
      gk.udpateSoluteParameters(i);
    }

    boolean usePerfectRadii = gk.getUsePerfectRadii();
    double perfectRadiiScale = 1.0;
    if (usePerfectRadii) {
      // No descreening when using perfect radii (OpenMM will just load the base radii).
      perfectRadiiScale = 0.0;
    }

    double[] baseRadii = gk.getBaseRadii();
    if (usePerfectRadii) {
      baseRadii = gk.getPerfectRadii();
    }
    double[] overlapScale = gk.getOverlapScale();
    double[] descreenRadius = gk.getDescreenRadii();
    double[] neckFactors = gk.getNeckScale();

    boolean nea = gk.getNativeEnvironmentApproximation();
    double lambdaElec = openMMEnergy.getPmeNode().getAlchemicalParameters().permLambda;

    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      double chargeUseFactor = 1.0;
      if (!atom.getUse() || !atom.getElectrostatics()) {
        chargeUseFactor = 0.0;
      }

      double lambdaScale = lambdaElec;
      if (!atom.applyLambda()) {
        lambdaScale = 1.0;
      }

      double baseSize = baseRadii[index] * OpenMM_NmPerAngstrom;
      double descreenSize = descreenRadius[index] * OpenMM_NmPerAngstrom * perfectRadiiScale;

      chargeUseFactor *= lambdaScale;
      double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;
      overlapScaleUseFactor = overlapScaleUseFactor * perfectRadiiScale;
      double overlap = overlapScale[index] * overlapScaleUseFactor;
      double neckFactor = neckFactors[index] * overlapScaleUseFactor;

      MultipoleType multipoleType = atom.getMultipoleType();
      setParticleParameters(index, multipoleType.charge * chargeUseFactor, baseSize, overlap, descreenSize, neckFactor);
    }

    // OpenMM Bug: Surface Area is not Updated by "updateParametersInContext"
    GeneralizedKirkwood.NonPolarModel nonpolar = gk.getNonPolarModel();
    switch (nonpolar) {
      case BORN_CAV_DISP, BORN_SOLV -> {
        // Configure a Born Radii based surface area term.
        double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm;
        setSurfaceAreaFactor(-surfaceTension);
      }
    }

    updateParametersInContext(openMMEnergy.getContext());
  }
}
