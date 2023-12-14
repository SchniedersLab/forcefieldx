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

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.parameters.MultipoleType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePair;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePairNoExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_SingleParticle;
import static java.lang.String.format;

/**
 * FixedChargeGBForce.
 */
public class FixedChargeGBForce extends OpenMMCustomGBForce {

  private static final Logger logger = Logger.getLogger(FixedChargeGBForce.class.getName());

  /**
   * FixedChargeGBForce constructor.
   *
   * @param openMMEnergy OpenMM Energy that contains the GK instance.
   */
  public FixedChargeGBForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      destroy();
      return;
    }

    double sTens = 0.0;
    if (gk.getNonPolarModel() == GeneralizedKirkwood.NonPolarModel.BORN_SOLV
        || gk.getNonPolarModel() == GeneralizedKirkwood.NonPolarModel.BORN_CAV_DISP) {
      sTens = gk.getSurfaceTension();
      sTens *= OpenMM_KJPerKcal;
      sTens *= 100.0; // 100 square Angstroms per square nanometer.
      // logger.info(String.format(" FFX surface tension: %9.5g kcal/mol/Ang^2", sTens));
      // logger.info(String.format(" OpenMM surface tension: %9.5g kJ/mol/nm^2", sTens));
    }

    addPerParticleParameter("q");
    addPerParticleParameter("radius");
    addPerParticleParameter("scale");
    addPerParticleParameter("surfaceTension");
    addGlobalParameter("solventDielectric", gk.getSolventPermittivity());
    addGlobalParameter("soluteDielectric", 1.0);
    addGlobalParameter("dOffset", gk.getDielecOffset() * OpenMM_NmPerAngstrom); // Factor of 0.1 for Ang to nm.
    addGlobalParameter("probeRadius", gk.getProbeRadius() * OpenMM_NmPerAngstrom);

    addComputedValue("I", """
        0.5*((1/L^3-1/U^3)/3.0+(1/U^4-1/L^4)/8.0*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);
        U=r+sr2;
        C=2/3*(1/or1^3-1/L^3)*step(sr2-r-or1);
        L = step(sr2 - r1r)*sr2mr + (1 - step(sr2 - r1r))*L;
        sr2mr = sr2 - r;
        r1r = radius1 + r;
        L = step(r1sr2 - r)*radius1 + (1 - step(r1sr2 - r))*L;
        r1sr2 = radius1 + sr2;
        L = r - sr2;
        sr2 = scale2 * radius2;
        or1 = radius1;
        or2 = radius2;
        """, OpenMM_CustomGBForce_ParticlePairNoExclusions);

    addComputedValue("B", """
        step(BB-radius)*BB + (1 - step(BB-radius))*radius;
        BB = 1 / ( (3.0*III)^(1.0/3.0) );
        III = step(II)*II + (1 - step(II))*1.0e-9/3.0;
        II = maxI - I;
        maxI = 1/(3.0*radius^3);
        """, OpenMM_CustomGBForce_SingleParticle);

    addEnergyTerm(
        "surfaceTension*(radius+probeRadius+dOffset)^2*((radius+dOffset)/B)^6/6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B",
        OpenMM_CustomGBForce_SingleParticle);

    // Particle pair term is the generalized Born cross term.
    addEnergyTerm("""
        -138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;
        f=sqrt(r^2+B1*B2*exp(-r^2/(2.455*B1*B2)));
        """, OpenMM_CustomGBForce_ParticlePair);

    double[] baseRadii = gk.getBaseRadii();
    double[] overlapScale = gk.getOverlapScale();
    OpenMMDoubleArray doubleArray = new OpenMMDoubleArray(0);
    MolecularAssembly molecularAssembly = openMMEnergy.getMolecularAssembly();
    Atom[] atoms = molecularAssembly.getAtomArray();
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      MultipoleType multipoleType = atoms[i].getMultipoleType();
      doubleArray.append(multipoleType.charge);
      doubleArray.append(OpenMM_NmPerAngstrom * baseRadii[i]);
      doubleArray.append(overlapScale[i]);
      doubleArray.append(sTens);
      addParticle(doubleArray);
      doubleArray.resize(0);
    }
    doubleArray.destroy();

    double cut = gk.getCutoff();
    setCutoffDistance(OpenMM_NmPerAngstrom * cut);

    int forceGroup = molecularAssembly.getForceField().getInteger("GK_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Custom generalized Born force \t%d", forceGroup));
  }

  /**
   * Construct a GB force.
   *
   * @param openMMEnergy OpenMM Energy that contains the GK instance.
   * @return The GB force.
   */
  public static OpenMMForce constructForce(OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
      return null;
    }
    return new FixedChargeGBForce(openMMEnergy);
  }

  /**
   * Update the GB force.
   *
   * @param atoms        The atoms to update.
   * @param openMMEnergy OpenMM Energy that contains the GK instance.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    double[] baseRadii = gk.getBaseRadii();
    double[] overlapScale = gk.getOverlapScale();
    boolean nea = gk.getNativeEnvironmentApproximation();

    double sTens = 0.0;
    if (gk.getNonPolarModel() == GeneralizedKirkwood.NonPolarModel.BORN_SOLV
        || gk.getNonPolarModel() == GeneralizedKirkwood.NonPolarModel.BORN_CAV_DISP) {
      sTens = gk.getSurfaceTension();
      sTens *= OpenMM_KJPerKcal;
      sTens *= 100.0; // 100 square Angstroms per square nanometer.
    }

    OpenMMDoubleArray parameters = new OpenMMDoubleArray(0);
    double lambdaElec = openMMEnergy.getSystem().getLambdaElec();
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

      chargeUseFactor *= lambdaScale;
      MultipoleType multipoleType = atom.getMultipoleType();
      double charge = multipoleType.charge * chargeUseFactor;
      double surfaceTension = sTens * chargeUseFactor;

      double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;
      double oScale = overlapScale[index] * overlapScaleUseFactor;
      double baseRadius = baseRadii[index];

      parameters.append(charge);
      parameters.append(OpenMM_NmPerAngstrom * baseRadius);
      parameters.append(oScale);
      parameters.append(surfaceTension);
      setParticleParameters(index, parameters);
      parameters.resize(0);
    }
    parameters.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }

}
