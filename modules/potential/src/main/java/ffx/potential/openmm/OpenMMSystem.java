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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMMLibrary;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.nonbonded.implicit.ChandlerCavitation;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setDefaultCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePair;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePairNoExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_SingleParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addComputedValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyTerm;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addInteractionGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_resize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setForceGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntSet_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntSet_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntSet_insert;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addConstraint;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getNumConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setDefaultPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setParticleMass;
import static ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.VDW_TYPE.LENNARD_JONES;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * Create and manage an OpenMM System.
 *
 * <p>The definition of a System involves four elements:
 *
 * <p>The particles and constraints are defined directly by the System object, while forces are
 * defined by objects that extend the Force class. After creating a System, call addParticle() once
 * for each particle, addConstraint() for each constraint, and addForce() for each Force.
 *
 * <p>In addition, particles may be designated as "virtual sites". These are particles whose
 * positions are computed automatically based on the positions of other particles. To define a
 * virtual site, call setVirtualSite(), passing in a VirtualSite object that defines the rules for
 * computing its position.
 */
public class OpenMMSystem {

  private static final Logger logger = Logger.getLogger(OpenMMSystem.class.getName());

  private static final double DEFAULT_MELD_SCALE_FACTOR = -1.0;
  /**
   * The ForceFieldEnergyOpenMM instance.
   */
  private final OpenMMEnergy openMMEnergy;
  /**
   * The OpenMMContext instance.
   */
  private final OpenMMContext openMMContext;
  private final double meldScaleFactor;
  /**
   * Andersen thermostat collision frequency.
   */
  private final double collisionFreq;
  /**
   * When using MELD, our goal will be to scale down the potential by this factor. A negative value
   * indicates we're not using MELD.
   */
  private final boolean useMeld;
  /**
   * The Force Field in use.
   */
  private final ForceField forceField;
  /**
   * Array of atoms in the system.
   */
  private final Atom[] atoms;
  /**
   * Number of atoms in the system.
   */
  private final int nAtoms;
  /**
   * OpenMM System.
   */
  private PointerByReference system;
  /**
   * Barostat to be added if NPT (isothermal-isobaric) dynamics is requested.
   */
  private PointerByReference ommBarostat = null;
  /**
   * OpenMM center-of-mass motion remover.
   */
  private PointerByReference commRemover = null;
  /**
   * This flag indicates bonded force constants and equilibria are updated (e.g. during ManyBody
   * titration).
   */
  private boolean updateBondedTerms = false;
  /**
   * OpenMM Custom Bond Force
   */
  private BondForce bondForce = null;
  /**
   * OpenMM Custom Angle Force
   */
  private AngleForce angleForce = null;
  /**
   * OpenMM Custom In-Plane Angle Force
   */
  private InPlaneAngleForce inPlaneAngleForce = null;
  /**
   * OpenMM Custom Stretch-Bend Force
   */
  private StretchBendForce stretchBendForce = null;
  /**
   * OpenMM Custom Urey-Bradley Force
   */
  private UreyBradleyForce ureyBradleyForce = null;
  /**
   * OpenMM Custom Out-of-Plane Bend Force
   */
  private OutOfPlaneBendForce outOfPlaneBendForce = null;
  /**
   * OpenMM Custom Pi-Torsion Force
   */
  private PiOrbitalTorsionForce piOrbitalTorsionForce = null;
  /**
   * OpenMM AMOEBA Torsion Force.
   */
  private TorsionForce torsionForce = null;
  /**
   * OpenMM Improper Torsion Force.
   */
  private ImproperTorsionForce improperTorsionForce = null;
  /**
   * OpenMM Torsion-Torsion Force.
   */
  private AmoebaTorsionTorsionForce amoebaTorsionTorsionForce = null;
  /**
   * OpenMM Stretch-Torsion Force.
   */
  private StretchTorsionForce stretchTorsionForce = null;
  /**
   * OpenMM Angle-Torsion Force.
   */
  private AngleTorsionForce angleTorsionForce = null;
  /**
   * OpenMM Restraint-Torsion Force.
   */
  private RestrainTorsionsForce restrainTorsionsForce = null;
  /**
   * OpenMM Restrain-Position Force.
   */
  private RestrainPositionsForce restrainPositionsForce = null;
  /**
   * OpenMM Restrain-Groups Force.
   */
  private RestrainGroupsForce restrainGroupsForce = null;

  /**
   * OpenMM AMOEBA van der Waals Force.
   */
  private AmoebaVDWForce amoebaVDWForce = null;
  /**
   * OpenMM AMOEBA Multipole Force.
   */
  private AmoebaMultipoleForce amoebaMultipoleForce = null;
  /**
   * OpenMM Generalized Kirkwood Force.
   */
  private AmoebaGeneralizedKirkwoodForce amoebaGeneralizedKirkwoodForce = null;
  /**
   * OpenMM AMOEBA WCA Dispersion Force.
   */
  private AmoebaWcaDispersionForce amoebaWcaDispersionForce = null;
  /**
   * OpenMM AMOEBA WCA Cavitation Force.
   */
  private AmoebaCavitationForce amoebaCavitationForce = null;
  /**
   * OpenMM Custom GB Force.
   */
  private PointerByReference customGBForce = null;
  /**
   * OpenMM Fixed Charge Non-Bonded Force.
   */
  private FixedChargeNonbondedForce fixedChargeNonBondedForce = null;
  /**
   * Fixed charge softcore vdW force boolean.
   */
  private boolean softcoreCreated = false;
  /**
   * Lambda flag to indicate control of electrostatic scaling. If both elec and vdW are being
   * scaled, then vdW is scaled first, followed by elec.
   */
  private boolean elecLambdaTerm;
  /**
   * Lambda flag to indicate control of vdW scaling. If both elec and vdW are being scaled, then
   * vdW is scaled first, followed by elec.
   */
  private boolean vdwLambdaTerm;
  /**
   * Lambda flag to indicate control of torsional force constants (L=0 corresponds to torsions
   * being off, and L=1 to torsions at full strength).
   */
  private boolean torsionLambdaTerm;
  /**
   * Value of the van der Waals lambda state variable.
   */
  private double lambdaVDW = 1.0;
  /**
   * Value of the electrostatics lambda state variable.
   */
  private double lambdaElec = 1.0;
  /**
   * Value of the electrostatics lambda state variable.
   */
  private double lambdaTorsion = 1.0;
  /**
   * The lambda value that defines when the electrostatics will start to turn on for full path
   * non-bonded term scaling.
   *
   * <p>A value of 0.6 works well for Chloride ion solvation, which is a difficult case due to the
   * ion having a formal negative charge and a large polarizability.
   */
  private double electrostaticStart = 0.6;
  /**
   * Electrostatics lambda is raised to this power.
   */
  private double electrostaticLambdaPower;
  /**
   * van der Waals softcore alpha.
   */
  private double vdWSoftcoreAlpha = 0.25;
  /**
   * OpenMM thermostat. Currently, an Andersen thermostat is supported.
   */
  private PointerByReference ommThermostat = null;
  /**
   * van der Waals softcore beta.
   */
  private double vdwSoftcorePower = 3.0;
  /**
   * Torsional lambda power.
   */
  private double torsionalLambdaPower = 2.0;

  /**
   * OpenMMSystem constructor.
   *
   * @param openMMEnergy ForceFieldEnergyOpenMM instance.
   */
  public OpenMMSystem(OpenMMEnergy openMMEnergy) {
    this.openMMEnergy = openMMEnergy;
    this.openMMContext = openMMEnergy.getContext();

    // Create the OpenMM System
    system = OpenMM_System_create();
    logger.info("\n System created");

    MolecularAssembly molecularAssembly = openMMEnergy.getMolecularAssembly();
    forceField = molecularAssembly.getForceField();
    atoms = molecularAssembly.getAtomArray();
    nAtoms = atoms.length;

    // Load atoms.
    try {
      addAtoms();
    } catch (Exception e) {
      logger.severe(" Atom without mass encountered.");
    }

    // Check for MELD use. If we're using MELD, set all lambda terms to true.
    meldScaleFactor = forceField.getDouble("MELD_SCALE_FACTOR", DEFAULT_MELD_SCALE_FACTOR);
    if (meldScaleFactor <= 1.0 && meldScaleFactor > 0.0) {
      useMeld = true;
      elecLambdaTerm = true;
      vdwLambdaTerm = true;
      torsionLambdaTerm = true;
    } else {
      useMeld = false;
      elecLambdaTerm = false;
      vdwLambdaTerm = false;
      torsionLambdaTerm = false;
    }

    // Read alchemical information -- this needs to be done before creating forces.
    elecLambdaTerm = forceField.getBoolean("ELEC_LAMBDATERM", elecLambdaTerm);
    vdwLambdaTerm = forceField.getBoolean("VDW_LAMBDATERM", vdwLambdaTerm);
    torsionLambdaTerm = forceField.getBoolean("TORSION_LAMBDATERM", torsionLambdaTerm);

    if (!forceField.getBoolean("LAMBDATERM", false)) {
      openMMEnergy.setLambdaTerm(elecLambdaTerm || vdwLambdaTerm || torsionLambdaTerm);
    } else {
      openMMEnergy.setLambdaTerm(true);
    }

    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW != null) {
      vdWSoftcoreAlpha = vdW.getAlpha();
      vdwSoftcorePower = (int) vdW.getBeta();
    }

    electrostaticStart = forceField.getDouble("PERMANENT_LAMBDA_START", electrostaticStart);
    if (electrostaticStart > 1.0) {
      electrostaticStart = 1.0;
    } else if (electrostaticStart < 0.0) {
      electrostaticStart = 0.0;
    }
    electrostaticLambdaPower = forceField.getDouble("PERMANENT_LAMBDA_EXPONENT", 2.0);

    if (useMeld) {
      // lambda path starts at 0.0
      openMMEnergy.setLambdaStart(0.0);
      // electrostaticStart is ignored for MELD.
      electrostaticStart = 0.0;
      // electrostaticLambdaPower is ignored for MELD.
      electrostaticLambdaPower = 1.0;
      // vdW is linearly scaled for MELD.
      vdwSoftcorePower = 1;
      // No softcore offset for MELD.
      vdWSoftcoreAlpha = 0.0;
      // Torsions are linearly scaled for MELD.
      torsionalLambdaPower = 1.0;
      // Only need single-sided dU/dL
      openMMEnergy.setTwoSidedFiniteDifference(false);
    }

    collisionFreq = forceField.getDouble("COLLISION_FREQ", 0.1);
  }

  /**
   * Add forces to the system.
   */
  public void addForces() {

    // Set up rigid constraints.
    // These flags need to be set before bonds and angles are created below.
    boolean rigidHydrogen = forceField.getBoolean("RIGID_HYDROGEN", false);
    boolean rigidBonds = forceField.getBoolean("RIGID_BONDS", false);
    boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    if (rigidHydrogen) {
      addHydrogenConstraints();
    }
    if (rigidBonds) {
      addUpBondConstraints();
    }
    if (rigidHydrogenAngles) {
      setUpHydrogenAngleConstraints();
    }

    logger.info("\n Bonded Terms\n");

    // Add Bond Force.
    if (rigidBonds) {
      logger.info(" Not creating AmoebaBondForce because bonds are constrained.");
    } else {
      bondForce = (BondForce) BondForce.constructForce(openMMEnergy);
      addForce(bondForce);
    }

    // Add Angle Force.
    angleForce = (AngleForce) AngleForce.constructForce(openMMEnergy);
    addForce(angleForce);

    // Add In-Plane Angle Force.
    inPlaneAngleForce = (InPlaneAngleForce) InPlaneAngleForce.constructForce(openMMEnergy);
    addForce(inPlaneAngleForce);

    // Add Stretch-Bend Force.
    stretchBendForce = (StretchBendForce) StretchBendForce.constructForce(openMMEnergy);
    addForce(stretchBendForce);

    // Add Urey-Bradley Force.
    ureyBradleyForce = (UreyBradleyForce) UreyBradleyForce.constructForce(openMMEnergy);
    addForce(ureyBradleyForce);

    // Out-of Plane Bend Force.
    outOfPlaneBendForce = (OutOfPlaneBendForce) OutOfPlaneBendForce.constructForce(openMMEnergy);
    addForce(outOfPlaneBendForce);

    // Add Pi-Torsion Force.
    piOrbitalTorsionForce = (PiOrbitalTorsionForce) PiOrbitalTorsionForce.constructForce(openMMEnergy);
    addForce(piOrbitalTorsionForce);

    // Add Torsion Force.
    torsionForce = (TorsionForce) TorsionForce.constructForce(openMMEnergy);
    addForce(torsionForce);

    // Add Improper Torsion Force.
    improperTorsionForce = (ImproperTorsionForce) ImproperTorsionForce.constructForce(openMMEnergy);
    addForce(improperTorsionForce);

    // Add Torsion-Torsion Force.
    amoebaTorsionTorsionForce = (AmoebaTorsionTorsionForce) AmoebaTorsionTorsionForce.constructForce(openMMEnergy);
    addForce(amoebaTorsionTorsionForce);

    // Add stretch-torsion coupling terms.
    stretchTorsionForce = (StretchTorsionForce) StretchTorsionForce.constructForce(openMMEnergy);
    addForce(stretchTorsionForce);

    // Add angle-torsion coupling terms.
    angleTorsionForce = (AngleTorsionForce) AngleTorsionForce.constructForce(openMMEnergy);
    addForce(angleTorsionForce);

    // Add Restraint-Torsions
    restrainTorsionsForce = (RestrainTorsionsForce) RestrainTorsionsForce.constructForce(openMMEnergy);
    addForce(restrainTorsionsForce);

    // Add Restrain-Position force.
    restrainPositionsForce = (RestrainPositionsForce) RestrainPositionsForce.constructForce(openMMEnergy);
    addForce(restrainPositionsForce);

    // Add a Restrain-Bond force for each functional form.
    for (BondType.BondFunction function : BondType.BondFunction.values()) {
      RestrainBondsForce restrainBondsForce = (RestrainBondsForce) RestrainBondsForce.constructForce(function, openMMEnergy);
      addForce(restrainBondsForce);
    }

    // Add Restrain-Groups force.
    restrainGroupsForce = (RestrainGroupsForce) RestrainGroupsForce.constructForce(openMMEnergy);
    addForce(restrainGroupsForce);

    setDefaultPeriodicBoxVectors();

    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW != null) {
      logger.info("\n Non-Bonded Terms\n");
      VanDerWaalsForm vdwForm = vdW.getVDWForm();
      if (vdwForm.vdwType == LENNARD_JONES) {
        fixedChargeNonBondedForce = (FixedChargeNonbondedForce) FixedChargeNonbondedForce.constructForce(openMMEnergy);
        addForce(fixedChargeNonBondedForce);
        if (fixedChargeNonBondedForce != null) {
          GeneralizedKirkwood gk = openMMEnergy.getGK();
          if (gk != null) {
            addCustomGBForce();
          }
        }
      } else {
        // Add vdW Force.
        amoebaVDWForce = (AmoebaVDWForce) AmoebaVDWForce.constructForce(openMMEnergy);
        addForce(amoebaVDWForce);

        // Add Multipole Force.
        amoebaMultipoleForce = (AmoebaMultipoleForce) AmoebaMultipoleForce.constructForce(openMMEnergy);
        addForce(amoebaMultipoleForce);

        // Add Generalized Kirkwood Force.
        GeneralizedKirkwood gk = openMMEnergy.getGK();
        if (gk != null) {
          amoebaGeneralizedKirkwoodForce = (AmoebaGeneralizedKirkwoodForce) AmoebaGeneralizedKirkwoodForce.constructForce(openMMEnergy);
          addForce(amoebaGeneralizedKirkwoodForce);

          // Add WCA Dispersion Force.
          DispersionRegion dispersionRegion = gk.getDispersionRegion();
          if (dispersionRegion != null) {
            amoebaWcaDispersionForce = (AmoebaWcaDispersionForce) AmoebaWcaDispersionForce.constructForce(openMMEnergy);
            addForce(amoebaWcaDispersionForce);
          }

          // Add a GaussVol Cavitation Force.
          ChandlerCavitation chandlerCavitation = gk.getChandlerCavitation();
          if (chandlerCavitation != null && chandlerCavitation.getGaussVol() != null) {
            amoebaCavitationForce = (AmoebaCavitationForce) AmoebaCavitationForce.constructForce(openMMEnergy);
            addForce(amoebaCavitationForce);
          }
        }
      }
    }

    if (openMMEnergy.getLambdaTerm()) {
      logger.info(format("\n Lambda path start:              %6.3f", openMMEnergy.getLambdaStart()));
      logger.info(format(" Lambda scales torsions:          %s", torsionLambdaTerm));
      if (torsionLambdaTerm) {
        logger.info(format(" torsion lambda power:           %6.3f", torsionalLambdaPower));
      }
      logger.info(format(" Lambda scales vdW interactions:  %s", vdwLambdaTerm));
      if (vdwLambdaTerm) {
        logger.info(format(" van Der Waals alpha:            %6.3f", vdWSoftcoreAlpha));
        logger.info(format(" van Der Waals lambda power:     %6.3f", vdwSoftcorePower));
      }
      logger.info(format(" Lambda scales electrostatics:    %s", elecLambdaTerm));

      if (elecLambdaTerm) {
        logger.info(format(" Electrostatics start:           %6.3f", electrostaticStart));
        logger.info(format(" Electrostatics lambda power:    %6.3f", electrostaticLambdaPower));
      }
      logger.info(format(" Using Meld:                      %s", useMeld));
      if (useMeld) {
        logger.info(format(" Meld scale factor:              %6.3f", meldScaleFactor));
      }
    }
  }

  public void addForce(OpenMMForce force) {
    if (force != null) {
      OpenMM_System_addForce(system, force.forcePointer);
    }
  }

  /**
   * Add an Andersen thermostat to the system.
   *
   * @param targetTemp Target temperature in Kelvins.
   */
  public void addAndersenThermostatForce(double targetTemp) {
    addAndersenThermostatForce(targetTemp, collisionFreq);
  }

  /**
   * Add an Andersen thermostat to the system.
   *
   * @param targetTemp    Target temperature in Kelvins.
   * @param collisionFreq Collision frequency in 1/psec.
   */
  public void addAndersenThermostatForce(double targetTemp, double collisionFreq) {
    if (ommThermostat == null) {
      ommThermostat = OpenMM_AndersenThermostat_create(targetTemp, collisionFreq);
      OpenMM_System_addForce(system, ommThermostat);
      logger.info("\n Adding an Andersen thermostat");
      logger.info(format("  Target Temperature:   %6.2f (K)", targetTemp));
      logger.info(format("  Collision Frequency:  %6.2f (1/psec)", collisionFreq));
    } else {
      OpenMM_AndersenThermostat_setDefaultTemperature(ommThermostat, targetTemp);
      OpenMM_AndersenThermostat_setDefaultCollisionFrequency(ommThermostat, collisionFreq);
      logger.fine(" Updated the Andersen thermostat");
      logger.fine(format("  Target Temperature:   %6.2f (K)", targetTemp));
      logger.fine(format("  Collision Frequency:  %6.2f (1/psec)", collisionFreq));
    }
  }

  /**
   * Adds a force that removes center-of-mass motion.
   */
  public void addCOMMRemoverForce() {
    int frequency = 100;
    if (commRemover == null) {
      commRemover = OpenMM_CMMotionRemover_create(frequency);
      OpenMM_System_addForce(system, commRemover);
      logger.info("\n Adding a center of mass motion remover");
      logger.info(format("  Frequency:            %6d", frequency));
    }
  }

  /**
   * Add a Monte Carlo Barostat to the system.
   *
   * @param targetPressure The target pressure (in atm).
   * @param targetTemp     The target temperature.
   * @param frequency      The frequency to apply the barostat.
   */
  public void addMonteCarloBarostatForce(double targetPressure, double targetTemp, int frequency) {
    if (ommBarostat == null) {
      double pressureInBar = targetPressure * Constants.ATM_TO_BAR;
      ommBarostat = OpenMM_MonteCarloBarostat_create(pressureInBar, targetTemp, frequency);
      CompositeConfiguration properties = openMMEnergy.getMolecularAssembly().getProperties();
      if (properties.containsKey("barostat-seed")) {
        int randomSeed = properties.getInt("barostat-seed", 0);
        logger.info(format(" Setting random seed %d for Monte Carlo Barostat", randomSeed));
        OpenMM_MonteCarloBarostat_setRandomNumberSeed(ommBarostat, randomSeed);
      }
      OpenMM_System_addForce(system, ommBarostat);
      logger.info("\n Adding a Monte Carlo barostat");
      logger.info(format("  Target Pressure:      %6.2f (atm)", targetPressure));
      logger.info(format("  Target Temperature:   %6.2f (K)", targetTemp));
      logger.info(format("  MC Move Frequency:    %6d", frequency));
    } else {
      logger.fine("\n Updating the Monte Carlo barostat");
      logger.fine(format("  Target Pressure:      %6.2f (atm)", targetPressure));
      logger.fine(format("  Target Temperature:   %6.2f (K)", targetTemp));
      logger.fine(format("  MC Move Frequency:    %6d", frequency));
    }
  }

  /**
   * Calculate the number of degrees of freedom.
   *
   * @return Number of degrees of freedom.
   */
  public int calculateDegreesOfFreedom() {
    // Begin from the 3 times the number of active atoms.
    int dof = openMMEnergy.getNumberOfVariables();
    // Remove OpenMM constraints.
    dof = dof - OpenMM_System_getNumConstraints(system);
    // Remove center of mass motion.
    if (commRemover != null) {
      dof -= 3;
    }
    return dof;
  }

  public double getTemperature(double kineticEnergy) {
    double dof = calculateDegreesOfFreedom();
    return 2.0 * kineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);
  }

  /**
   * Get the ForceField in use.
   *
   * @return ForceField instance.
   */
  public ForceField getForceField() {
    return forceField;
  }

  /**
   * Destroy the system.
   */
  public void free() {
    if (system != null) {
      logger.fine(" Free OpenMM system.");
      OpenMM_System_destroy(system);
      logger.fine(" Free OpenMM system completed.");
      system = null;
    }
  }

  /**
   * Print current lambda values.
   */
  public void printLambdaValues() {
    logger.info(format("\n Lambda Values\n Torsion: %6.3f vdW: %6.3f Elec: %6.3f ", lambdaTorsion,
        lambdaVDW, lambdaElec));
  }

  /**
   * Set the overall lambda value for the system.
   *
   * @param lambda Current lambda value.
   */
  public void setLambda(double lambda) {

    // Initially set all lambda values to 1.0.
    lambdaTorsion = 1.0;

    // Applied to softcore vdW forces.
    lambdaVDW = 1.0;

    // Applied to normal electrostatic parameters for alchemical atoms.
    lambdaElec = 1.0;

    if (torsionLambdaTerm) {
      // Multiply torsional potentials by L^2 (dU/dL = 0 at L=0).
      lambdaTorsion = pow(lambda, torsionalLambdaPower);
      if (useMeld) {
        lambdaTorsion = meldScaleFactor + lambda * (1.0 - meldScaleFactor);
      }
    }

    if (elecLambdaTerm && vdwLambdaTerm) {
      // Lambda effects both vdW and electrostatics.
      if (lambda < electrostaticStart) {
        // Begin turning vdW on with electrostatics off.
        lambdaElec = 0.0;
      } else {
        // Turn electrostatics on during the latter part of the path.
        double elecWindow = 1.0 - electrostaticStart;
        lambdaElec = (lambda - electrostaticStart) / elecWindow;
        lambdaElec = pow(lambdaElec, electrostaticLambdaPower);
      }
      lambdaVDW = lambda;
      if (useMeld) {
        lambdaElec = sqrt(meldScaleFactor + lambda * (1.0 - meldScaleFactor));
        lambdaVDW = meldScaleFactor + lambda * (1.0 - meldScaleFactor);
      }
    } else if (vdwLambdaTerm) {
      // Lambda effects vdW, with electrostatics turned off.
      lambdaElec = 0.0;
      lambdaVDW = lambda;
      if (useMeld) {
        lambdaVDW = meldScaleFactor + lambda * (1.0 - meldScaleFactor);
      }

    } else if (elecLambdaTerm) {
      // Lambda effects electrostatics, but not vdW.
      lambdaElec = lambda;
      if (useMeld) {
        lambdaElec = sqrt(meldScaleFactor + lambda * (1.0 - meldScaleFactor));
      }
    }
  }

  public boolean getVdwLambdaTerm() {
    return vdwLambdaTerm;
  }

  public void setUpdateBondedTerms(boolean updateBondedTerms) {
    this.updateBondedTerms = updateBondedTerms;
  }

  /**
   * Return a reference to the System.
   *
   * @return System referenece.
   */
  PointerByReference getSystem() {
    return system;
  }

  /**
   * Set the default values of the vectors defining the axes of the periodic box (measured in nm).
   *
   * <p>Any newly created Context will have its box vectors set to these. They will affect any
   * Force added to the System that uses periodic boundary conditions.
   *
   * <p>Triclinic boxes are supported, but the vectors must satisfy certain requirements. In
   * particular, a must point in the x direction, b must point "mostly" in the y direction, and c
   * must point "mostly" in the z direction. See the documentation for details.
   */
  private void setDefaultPeriodicBoxVectors() {
    Crystal crystal = openMMEnergy.getCrystal();
    if (!crystal.aperiodic()) {
      OpenMM_Vec3 a = new OpenMM_Vec3();
      OpenMM_Vec3 b = new OpenMM_Vec3();
      OpenMM_Vec3 c = new OpenMM_Vec3();
      double[][] Ai = crystal.Ai;
      a.x = Ai[0][0] * OpenMM_NmPerAngstrom;
      a.y = Ai[0][1] * OpenMM_NmPerAngstrom;
      a.z = Ai[0][2] * OpenMM_NmPerAngstrom;
      b.x = Ai[1][0] * OpenMM_NmPerAngstrom;
      b.y = Ai[1][1] * OpenMM_NmPerAngstrom;
      b.z = Ai[1][2] * OpenMM_NmPerAngstrom;
      c.x = Ai[2][0] * OpenMM_NmPerAngstrom;
      c.y = Ai[2][1] * OpenMM_NmPerAngstrom;
      c.z = Ai[2][2] * OpenMM_NmPerAngstrom;
      OpenMM_System_setDefaultPeriodicBoxVectors(system, a, b, c);
    }
  }

  public double getLambdaElec() {
    return lambdaElec;
  }

  /**
   * Update parameters if the Use flags and/or Lambda value has changed.
   *
   * @param atoms Atoms in this list are considered.
   */
  public void updateParameters(@Nullable Atom[] atoms) {
    if (vdwLambdaTerm) {
      if (fixedChargeNonBondedForce != null) {
        if (!softcoreCreated) {
          addCustomNonbondedSoftcoreForce();
          // Re-initialize the context.
          openMMContext.reinitialize();
          softcoreCreated = true;
        }
        openMMContext.setParameter("vdw_lambda", lambdaVDW);
      } else if (amoebaVDWForce != null) {
        openMMContext.setParameter("AmoebaVdwLambda", lambdaVDW);
        if (softcoreCreated) {
          // Avoid any updateParametersInContext calls if vdwLambdaTerm is true, but not other
          // alchemical terms.
          if (!torsionLambdaTerm && !elecLambdaTerm) {
            return;
          }
        } else {
          softcoreCreated = true;
        }
      }
    }

    // Note Stretch-Torsion and Angle-Torsion terms (for nucleic acids)
    // and Torsion-Torsion terms (for protein backbones) are not udpated yet.
    if (updateBondedTerms) {
      if (bondForce != null) {
        bondForce.updateForce(openMMEnergy);
      }
      if (angleForce != null) {
        angleForce.updateForce(openMMEnergy);
      }
      if (inPlaneAngleForce != null) {
        inPlaneAngleForce.updateForce(openMMEnergy);
      }
      if (stretchBendForce != null) {
        stretchBendForce.updateForce(openMMEnergy);
      }
      if (ureyBradleyForce != null) {
        ureyBradleyForce.updateForce(openMMEnergy);
      }
      if (outOfPlaneBendForce != null) {
        outOfPlaneBendForce.updateForce(openMMEnergy);
      }
      if (piOrbitalTorsionForce != null) {
        piOrbitalTorsionForce.updateForce(openMMEnergy);
      }
    }

    if (torsionLambdaTerm || updateBondedTerms) {
      if (torsionForce != null) {
        torsionForce.setLambdaTorsion(lambdaTorsion);
        torsionForce.updateForce(openMMEnergy);
      }
      if (improperTorsionForce != null) {
        improperTorsionForce.setLambdaTorsion(lambdaTorsion);
        improperTorsionForce.updateForce(openMMEnergy);
      }
    }

    if (restrainTorsionsForce != null) {
      restrainTorsionsForce.updateForce(openMMEnergy);
    }

    if (atoms == null || atoms.length == 0) {
      return;
    }

    // Update fixed charge non-bonded parameters.
    if (fixedChargeNonBondedForce != null) {
      fixedChargeNonBondedForce.updateForce(atoms, openMMEnergy);
    }

    // Update fixed charge GB parameters.
    if (customGBForce != null) {
      updateCustomGBForce(atoms);
    }

    // Update AMOEBA vdW parameters.
    if (amoebaVDWForce != null) {
      amoebaVDWForce.updateForce(atoms, openMMEnergy);
    }

    // Update AMOEBA polarizable multipole parameters.
    if (amoebaMultipoleForce != null) {
      amoebaMultipoleForce.updateForce(atoms, openMMEnergy);
    }

    // Update GK force.
    if (amoebaGeneralizedKirkwoodForce != null) {
      amoebaGeneralizedKirkwoodForce.updateForce(atoms, openMMEnergy);
    }

    // Update WCA Force.
    if (amoebaWcaDispersionForce != null) {
      amoebaWcaDispersionForce.updateForce(atoms, openMMEnergy);
    }

    // Update WCA Force.
    if (amoebaCavitationForce != null) {
      amoebaCavitationForce.updateForce(atoms, openMMEnergy);
    }
  }

  /**
   * Adds atoms from the molecular assembly to the OpenMM System and reports to the user the number
   * of particles added.
   */
  private void addAtoms() throws Exception {
    double totalMass = 0.0;
    for (Atom atom : atoms) {
      double mass = atom.getMass();
      totalMass += mass;
      if (mass < 0.0) {
        throw new Exception(" Atom with mass less than 0.");
      }
      if (mass == 0.0) {
        logger.info(format(" Atom %s has zero mass.", atom));
      }
      OpenMM_System_addParticle(system, mass);
    }
    logger.log(Level.INFO, format("  Atoms \t\t%6d", atoms.length));
    logger.log(Level.INFO, format("  Mass  \t\t%12.3f", totalMass));
  }

  /**
   * This method sets the mass of inactive atoms to zero.
   */
  public void updateAtomMass() {
    int index = 0;
    for (Atom atom : atoms) {
      double mass = 0.0;
      if (atom.isActive()) {
        mass = atom.getMass();
      }
      OpenMM_System_setParticleMass(system, index++, mass);
    }
  }

  /**
   * 1. Handle interactions between non-alchemical atoms with our default OpenMM NonBondedForce.
   * Note that alchemical atoms must have eps=0 to turn them off in this force.
   * <p>
   * 2. Handle interactions between alchemical atoms and mixed non-alchemical <-> alchemical
   * interactions with an OpenMM CustomNonBondedForce.
   */
  private void addCustomNonbondedSoftcoreForce() {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      return;
    }

    /*
     Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
     for epsilon is supported.
    */
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    if (vdwForm.vdwType != LENNARD_JONES || vdwForm.radiusRule != ARITHMETIC
        || vdwForm.epsilonRule != GEOMETRIC) {
      logger.info(format(" VDW Type:         %s", vdwForm.vdwType));
      logger.info(format(" VDW Radius Rule:  %s", vdwForm.radiusRule));
      logger.info(format(" VDW Epsilon Rule: %s", vdwForm.epsilonRule));
      logger.log(Level.SEVERE, " Unsupported van der Waals functional form.");
      return;
    }

    // Sterics mixing rules.
    String stericsMixingRules = " epsilon = sqrt(epsilon1*epsilon2);";
    stericsMixingRules += " rmin = 0.5 * (sigma1 + sigma2) * 1.122462048309372981;";

    // Softcore Lennard-Jones, with a form equivalent to that used in FFX VanDerWaals class.
    String stericsEnergyExpression = "(vdw_lambda^beta)*epsilon*x*(x-2.0);";
    // Effective softcore distance for sterics.
    stericsEnergyExpression += " x = 1.0 / (alpha*(1.0-vdw_lambda)^2.0 + (r/rmin)^6.0);";
    // Define energy expression for sterics.
    String energyExpression = stericsEnergyExpression + stericsMixingRules;

    PointerByReference fixedChargeSoftcore = OpenMM_CustomNonbondedForce_create(energyExpression);

    // Get the Alpha and Beta constants from the VanDerWaals instance.
    double alpha = vdWSoftcoreAlpha;
    double beta = vdwSoftcorePower;

    OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "vdw_lambda", 1.0);
    OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "alpha", alpha);
    OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "beta", beta);
    OpenMM_CustomNonbondedForce_addPerParticleParameter(fixedChargeSoftcore, "sigma");
    OpenMM_CustomNonbondedForce_addPerParticleParameter(fixedChargeSoftcore, "epsilon");

    // Add particles.
    PointerByReference alchemicalGroup = OpenMM_IntSet_create();
    PointerByReference nonAlchemicalGroup = OpenMM_IntSet_create();
    DoubleByReference charge = new DoubleByReference();
    DoubleByReference sigma = new DoubleByReference();
    DoubleByReference eps = new DoubleByReference();

    int index = 0;
    for (Atom atom : atoms) {
      if (atom.applyLambda()) {
        OpenMM_IntSet_insert(alchemicalGroup, index);
      } else {
        OpenMM_IntSet_insert(nonAlchemicalGroup, index);
      }

      fixedChargeNonBondedForce.getParticleParameters(index, charge, sigma, eps);
      double sigmaValue = sigma.getValue();
      double epsValue = eps.getValue();

      // Handle cases where sigma is 0.0; for example Amber99 tyrosine hydrogen atoms.
      if (sigmaValue == 0.0) {
        sigmaValue = 1.0;
        epsValue = 0.0;
      }

      PointerByReference particleParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(particleParameters, sigmaValue);
      OpenMM_DoubleArray_append(particleParameters, epsValue);
      OpenMM_CustomNonbondedForce_addParticle(fixedChargeSoftcore, particleParameters);
      OpenMM_DoubleArray_destroy(particleParameters);

      index++;
    }

    OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup, alchemicalGroup);
    OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup, nonAlchemicalGroup);
    OpenMM_IntSet_destroy(alchemicalGroup);
    OpenMM_IntSet_destroy(nonAlchemicalGroup);

    Crystal crystal = openMMEnergy.getCrystal();
    if (crystal.aperiodic()) {
      OpenMM_CustomNonbondedForce_setNonbondedMethod(fixedChargeSoftcore,
          OpenMMLibrary.OpenMM_CustomNonbondedForce_NonbondedMethod.OpenMM_CustomNonbondedForce_NoCutoff);
    } else {
      OpenMM_CustomNonbondedForce_setNonbondedMethod(fixedChargeSoftcore,
          OpenMMLibrary.OpenMM_CustomNonbondedForce_NonbondedMethod.OpenMM_CustomNonbondedForce_CutoffPeriodic);
    }

    NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
    double off = nonbondedCutoff.off;
    double cut = nonbondedCutoff.cut;
    if (cut == off) {
      logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
      if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
        logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
        cut = 1E40;
      } else {
        logger.info(format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut, off));
        cut *= 0.99;
      }
    }

    OpenMM_CustomNonbondedForce_setCutoffDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * off);
    OpenMM_CustomNonbondedForce_setUseSwitchingFunction(fixedChargeSoftcore, OpenMM_True);
    OpenMM_CustomNonbondedForce_setSwitchingDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * cut);

    // Add energy parameter derivative
    // OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(fixedChargeSoftcore,
    // "vdw_lambda");

    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);

    OpenMM_Force_setForceGroup(fixedChargeSoftcore, forceGroup);
    OpenMM_System_addForce(system, fixedChargeSoftcore);

    // Alchemical with Alchemical could be either softcore or normal interactions (softcore here).
    PointerByReference alchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(stericsEnergyExpression);

    // Non-Alchemical with Alchemical is essentially always softcore.
    PointerByReference nonAlchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(stericsEnergyExpression);

    // Currently both are treated the same (so we could condense the code below).
    OpenMM_CustomBondForce_addPerBondParameter(alchemicalAlchemicalStericsForce, "rmin");
    OpenMM_CustomBondForce_addPerBondParameter(alchemicalAlchemicalStericsForce, "epsilon");
    OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "vdw_lambda", 1.0);
    OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "alpha", alpha);
    OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "beta", beta);

    OpenMM_CustomBondForce_addPerBondParameter(nonAlchemicalAlchemicalStericsForce, "rmin");
    OpenMM_CustomBondForce_addPerBondParameter(nonAlchemicalAlchemicalStericsForce, "epsilon");
    OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "vdw_lambda",
        1.0);
    OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "alpha", alpha);
    OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "beta", beta);

    int range = fixedChargeNonBondedForce.getNumExceptions();

    IntByReference atomi = new IntByReference();
    IntByReference atomj = new IntByReference();
    int[][] torsionMask = vdW.getMask14();

    for (int i = 0; i < range; i++) {
      fixedChargeNonBondedForce.getExceptionParameters(i, atomi, atomj, charge, sigma, eps);

      // Omit both Exclusions (1-2, 1-3) and Exceptions (scaled 1-4) from the
      // CustomNonbondedForce.
      OpenMM_CustomNonbondedForce_addExclusion(fixedChargeSoftcore, atomi.getValue(),
          atomj.getValue());

      // Deal with scaled 1-4 torsions using the CustomBondForce
      int[] maskI = torsionMask[atomi.getValue()];
      int jID = atomj.getValue();
      boolean epsException = false;

      for (int mask : maskI) {
        if (mask == jID) {
          epsException = true;
          break;
        }
      }

      if (epsException) {
        Atom atom1 = atoms[atomi.getValue()];
        Atom atom2 = atoms[atomj.getValue()];

        boolean bothAlchemical = false;
        boolean oneAlchemical = false;

        if (atom1.applyLambda() && atom2.applyLambda()) {
          bothAlchemical = true;
        } else if ((atom1.applyLambda() && !atom2.applyLambda()) || (!atom1.applyLambda()
            && atom2.applyLambda())) {
          oneAlchemical = true;
        }

        if (bothAlchemical) {
          PointerByReference bondParameters = OpenMM_DoubleArray_create(0);
          OpenMM_DoubleArray_append(bondParameters, sigma.getValue() * 1.122462048309372981);
          OpenMM_DoubleArray_append(bondParameters, eps.getValue());
          OpenMM_CustomBondForce_addBond(alchemicalAlchemicalStericsForce, atomi.getValue(),
              atomj.getValue(), bondParameters);
          OpenMM_DoubleArray_destroy(bondParameters);
        } else if (oneAlchemical) {
          PointerByReference bondParameters = OpenMM_DoubleArray_create(0);
          OpenMM_DoubleArray_append(bondParameters, sigma.getValue() * 1.122462048309372981);
          OpenMM_DoubleArray_append(bondParameters, eps.getValue());
          OpenMM_CustomBondForce_addBond(nonAlchemicalAlchemicalStericsForce, atomi.getValue(),
              atomj.getValue(), bondParameters);
          OpenMM_DoubleArray_destroy(bondParameters);
        }
      }
    }

    OpenMM_Force_setForceGroup(alchemicalAlchemicalStericsForce, forceGroup);
    OpenMM_Force_setForceGroup(nonAlchemicalAlchemicalStericsForce, forceGroup);

    // OpenMM_CustomBondForce_addEnergyParameterDerivative(alchemicalAlchemicalStericsForce,
    // "vdw_lambda");
    // OpenMM_CustomBondForce_addEnergyParameterDerivative(nonAlchemicalAlchemicalStericsForce,
    // "vdw_lambda");

    OpenMM_System_addForce(system, alchemicalAlchemicalStericsForce);
    OpenMM_System_addForce(system, nonAlchemicalAlchemicalStericsForce);

    logger.log(Level.INFO, format("  Added fixed charge softcore force \t%d", forceGroup));
    logger.log(Level.INFO, format("   Alpha = %8.6f and beta = %8.6f", alpha, beta));
  }

  /**
   * Add a custom GB force to the OpenMM System.
   */
  private void addCustomGBForce() {
    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk == null) {
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

    customGBForce = OpenMM_CustomGBForce_create();
    OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "q");
    OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "radius");
    OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "scale");
    OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "surfaceTension");

    OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "solventDielectric",
        gk.getSolventPermittivity());
    OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "soluteDielectric", 1.0);
    OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "dOffset",
        gk.getDielecOffset() * OpenMM_NmPerAngstrom); // Factor of 0.1 for Ang to nm.
    OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "probeRadius",
        gk.getProbeRadius() * OpenMM_NmPerAngstrom);

    OpenMM_CustomGBForce_addComputedValue(customGBForce, "I",
        // "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
        // "step(r+sr2-or1)*0.5*((1/L^3-1/U^3)/3+(1/U^4-1/L^4)/8*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);"
        "0.5*((1/L^3-1/U^3)/3.0+(1/U^4-1/L^4)/8.0*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);"
            + "U=r+sr2;"
            // + "C=2*(1/or1-1/L)*step(sr2-r-or1);"
            + "C=2/3*(1/or1^3-1/L^3)*step(sr2-r-or1);"
            // + "L=step(or1-D)*or1 + (1-step(or1-D))*D;"
            // + "D=step(r-sr2)*(r-sr2) + (1-step(r-sr2))*(sr2-r);"
            + "L = step(sr2 - r1r)*sr2mr + (1 - step(sr2 - r1r))*L;" + "sr2mr = sr2 - r;"
            + "r1r = radius1 + r;" + "L = step(r1sr2 - r)*radius1 + (1 - step(r1sr2 - r))*L;"
            + "r1sr2 = radius1 + sr2;" + "L = r - sr2;" + "sr2 = scale2 * radius2;"
            + "or1 = radius1; or2 = radius2", OpenMM_CustomGBForce_ParticlePairNoExclusions);

    OpenMM_CustomGBForce_addComputedValue(customGBForce, "B",
        // "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
        // "psi=I*or; or=radius-0.009"
        "step(BB-radius)*BB + (1 - step(BB-radius))*radius;" + "BB = 1 / ( (3.0*III)^(1.0/3.0) );"
            + "III = step(II)*II + (1 - step(II))*1.0e-9/3.0;" + "II = maxI - I;"
            + "maxI = 1/(3.0*radius^3)", OpenMM_CustomGBForce_SingleParticle);

    OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
        "surfaceTension*(radius+probeRadius+dOffset)^2*((radius+dOffset)/B)^6/6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B",
        OpenMM_CustomGBForce_SingleParticle);

    // Particle pair term is the generalized Born cross term.
    OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
        "-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
            + "f=sqrt(r^2+B1*B2*exp(-r^2/(2.455*B1*B2)))", OpenMM_CustomGBForce_ParticlePair);

    double[] baseRadii = gk.getBaseRadii();
    double[] overlapScale = gk.getOverlapScale();
    PointerByReference doubleArray = OpenMM_DoubleArray_create(0);
    for (int i = 0; i < nAtoms; i++) {
      MultipoleType multipoleType = atoms[i].getMultipoleType();
      OpenMM_DoubleArray_append(doubleArray, multipoleType.charge);
      OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadii[i]);
      OpenMM_DoubleArray_append(doubleArray, overlapScale[i]);
      OpenMM_DoubleArray_append(doubleArray, sTens);

      OpenMM_CustomGBForce_addParticle(customGBForce, doubleArray);
      OpenMM_DoubleArray_resize(doubleArray, 0);
    }
    OpenMM_DoubleArray_destroy(doubleArray);

    double cut = gk.getCutoff();
    OpenMM_CustomGBForce_setCutoffDistance(customGBForce, cut);

    int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
    OpenMM_Force_setForceGroup(customGBForce, forceGroup);
    OpenMM_System_addForce(system, customGBForce);

    logger.log(Level.INFO, format("  Custom generalized Born force \t%d", forceGroup));
  }

  /**
   * Updates the custom GB force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateCustomGBForce(Atom[] atoms) {
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

    PointerByReference doubleArray = OpenMM_DoubleArray_create(0);
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

      OpenMM_DoubleArray_append(doubleArray, charge);
      OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadius);
      OpenMM_DoubleArray_append(doubleArray, oScale);
      OpenMM_DoubleArray_append(doubleArray, surfaceTension);
      OpenMM_CustomGBForce_setParticleParameters(customGBForce, index, doubleArray);

      // Reset the double array for the next atom.
      OpenMM_DoubleArray_resize(doubleArray, 0);
    }
    OpenMM_DoubleArray_destroy(doubleArray);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomGBForce_updateParametersInContext(customGBForce, openMMContext.getContextPointer());
    }
  }

  public boolean hasAmoebaCavitationForce() {
    return amoebaCavitationForce != null;
  }

  /**
   * Add a constraint to every bond.
   */
  private void addUpBondConstraints() {
    Bond[] bonds = openMMEnergy.getBonds();
    if (bonds == null || bonds.length < 1) {
      return;
    }

    logger.info(" Adding constraints for all bonds.");
    for (Bond bond : bonds) {
      Atom atom1 = bond.getAtom(0);
      Atom atom2 = bond.getAtom(1);
      int iAtom1 = atom1.getXyzIndex() - 1;
      int iAtom2 = atom2.getXyzIndex() - 1;
      OpenMM_System_addConstraint(system, iAtom1, iAtom2,
          bond.bondType.distance * OpenMM_NmPerAngstrom);
    }
  }

  /**
   * Add a constraint to every bond that includes a hydrogen atom.
   */
  private void addHydrogenConstraints() {
    Bond[] bonds = openMMEnergy.getBonds();
    if (bonds == null || bonds.length < 1) {
      return;
    }

    logger.info(" Adding constraints for hydrogen bonds.");
    for (Bond bond : bonds) {
      Atom atom1 = bond.getAtom(0);
      Atom atom2 = bond.getAtom(1);
      if (atom1.isHydrogen() || atom2.isHydrogen()) {
        BondType bondType = bond.bondType;
        int iAtom1 = atom1.getXyzIndex() - 1;
        int iAtom2 = atom2.getXyzIndex() - 1;
        OpenMM_System_addConstraint(system, iAtom1, iAtom2,
            bondType.distance * OpenMM_NmPerAngstrom);
      }
    }
  }

  /**
   * Add a constraint to every angle that includes two hydrogen atoms.
   */
  private void setUpHydrogenAngleConstraints() {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return;
    }

    logger.info(" Adding hydrogen angle constraints.");

    for (Angle angle : angles) {
      if (isHydrogenAngle(angle)) {
        Atom atom1 = angle.getAtom(0);
        Atom atom3 = angle.getAtom(2);

        // Calculate a "false bond" length between atoms 1 and 3 to constrain the angle using the
        // law of cosines.
        Bond bond1 = angle.getBond(0);
        double distance1 = bond1.bondType.distance;

        Bond bond2 = angle.getBond(1);
        double distance2 = bond2.bondType.distance;

        // Equilibrium angle value in degrees.
        double angleVal = angle.angleType.angle[angle.nh];

        // Law of cosines.
        double falseBondLength = sqrt(
            distance1 * distance1 + distance2 * distance2 - 2.0 * distance1 * distance2 * cos(
                toRadians(angleVal)));

        int iAtom1 = atom1.getXyzIndex() - 1;
        int iAtom3 = atom3.getXyzIndex() - 1;
        OpenMM_System_addConstraint(system, iAtom1, iAtom3,
            falseBondLength * OpenMM_NmPerAngstrom);
      }
    }
  }

  /**
   * Check to see if an angle is a hydrogen angle. This method only returns true for hydrogen
   * angles that are less than 160 degrees.
   *
   * @param angle Angle to check.
   * @return boolean indicating whether an angle is a hydrogen angle that is less than 160 degrees.
   */
  private boolean isHydrogenAngle(Angle angle) {
    if (angle.containsHydrogen()) {
      // Equilibrium angle value in degrees.
      double angleVal = angle.angleType.angle[angle.nh];
      // Make sure angle is less than 160 degrees because greater than 160 degrees will not be
      // constrained
      // well using the law of cosines.
      if (angleVal < 160.0) {
        Atom atom1 = angle.getAtom(0);
        Atom atom2 = angle.getAtom(1);
        Atom atom3 = angle.getAtom(2);
        // Setting constraints only on angles where atom1 or atom3 is a hydrogen while atom2 is
        // not a hydrogen.
        return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
      }
    }
    return false;
  }

}
