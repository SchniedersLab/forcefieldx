package ffx.potential.openmm;

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMMLibrary;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.RestraintTorsion;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.RestrainGroups;
import ffx.potential.nonbonded.RestrainPosition;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.nonbonded.implicit.ChandlerCavitation;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.nonbonded.implicit.GaussVol;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.nonbonded.pme.SCFAlgorithm;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ImproperTorsionType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiOrbitalTorsionType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWPairType;
import ffx.potential.parameters.VDWType;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_set;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_NonbondedMethod.OpenMM_AmoebaGKCavitationForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent12;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent13;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent14;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent15;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent11;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_Bisector;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_NoAxisType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ThreeFold;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZBisect;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZOnly;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZThenX;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_PME;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Direct;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Extrapolated;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Mutual;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_addMultipole;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setAEwald;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCovalentMap;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMultipoleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPmeGridDimensions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPolarizationType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_AlchemicalMethod.OpenMM_AmoebaVdwForce_Decouple;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_CutoffPeriodic;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticleType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addTypePair;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setAlchemicalMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleExclusions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcoreAlpha;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcorePower;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_updateParametersInContext;
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
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static edu.uiowa.jopenmm.OpenMMLibrary.*;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePair;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePairNoExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_SingleParticle;
import static ffx.potential.nonbonded.GeneralizedKirkwood.NonPolarModel.GAUSS_DISP;
import static ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_SIZE.RADIUS;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_TYPE.R_MIN;
import static ffx.potential.nonbonded.VanDerWaalsForm.VDW_TYPE.LENNARD_JONES;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
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
  ForceField forceField;
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
   * If true, all torsions are treated as 6-fold, and all angles are treated as possibly changing
   * between normal and in-plane types.
   */
  private final boolean manyBodyTitration;
  /**
   * OpenMM Custom Bond Force
   */
  private PointerByReference bondForce = null;
  /**
   * OpenMM Custom Angle Force
   */
  private PointerByReference angleForce = null;
  /**
   * OpenMM Custom Stretch-Bend Force
   */
  private PointerByReference stretchBendForce = null;
  /**
   * OpenMM Custom In-Plane Angle Force
   */
  private PointerByReference inPlaneAngleForce = null;
  /**
   * OpenMM Custom Urey-Bradley Force
   */
  private PointerByReference ureyBradleyForce = null;
  /**
   * OpenMM Custom Out-of-Plane Bend Force
   */
  private PointerByReference outOfPlaneBendForce = null;
  /**
   * OpenMM Custom Pi-Torsion Force
   */
  private PointerByReference piTorsionForce = null;
  /**
   * OpenMM AMOEBA Torsion Force.
   */
  private PointerByReference torsionForce = null;
  private PointerByReference[] restraintTorsions = null;
  /**
   * OpenMM Improper Torsion Force.
   */
  private PointerByReference improperTorsionForce = null;
  /**
   * OpenMM AMOEBA van der Waals Force.
   */
  private PointerByReference amoebaVDWForce = null;
  /**
   * OpenMM AMOEBA Multipole Force.
   */
  private PointerByReference amoebaMultipoleForce = null;
  /**
   * OpenMM Generalized Kirkwood Force.
   */
  private PointerByReference amoebaGeneralizedKirkwoodForce = null;
  /**
   * OpenMM AMOEBA WCA Dispersion Force.
   */
  private PointerByReference amoebaWcaDispersionForce = null;
  /**
   * OpenMM AMOEBA WCA Cavitation Force.
   */
  private PointerByReference amoebaCavitationForce = null;
  /**
   * OpenMM Custom GB Force.
   */
  private PointerByReference customGBForce = null;
  /**
   * OpenMM Fixed Charge Non-Bonded Force.
   */
  private PointerByReference fixedChargeNonBondedForce = null;
  /**
   * Fixed charge softcore vdW force boolean.
   */
  private boolean softcoreCreated = false;
  /**
   * Boolean array, holds charge exclusion list.
   */
  private boolean[] chargeExclusion;
  /**
   * Boolean array, holds van Der Waals exclusion list.
   */
  private boolean[] vdWExclusion;
  /**
   * Double array, holds charge quantity value for exceptions.
   */
  private double[] exceptionChargeProd;
  /**
   * Double array, holds epsilon quantity value for exceptions.
   */
  private double[] exceptionEps;
  /**
   * A map from vdW class values to OpenMM vdW types.
   */
  private Map<Integer, Integer> vdwClassToOpenMMType;
  /**
   * A class for a special vdW type that specifies zero energy (eps = 0.0; sigma = 1.0) for use
   * with the FFX "use" flag (e.g. use = false should give zero vdW energy for a many-body
   * expansion).
   */
  private int vdWClassForNoInteraction;
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

    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);

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

    // Set up rigid constraints. These flags need to be set before bonds and angles are created
    // below.
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
      addBondForce();
    }

    // Add Angle Force.
    addAngleForce();
    addInPlaneAngleForce();

    // Add Stretch-Bend Force.
    addStretchBendForce();

    // Add Urey-Bradley Force.
    addUreyBradleyForce();

    // Out-of Plane Bend Force.
    addOutOfPlaneBendForce();

    // Add Torsion Force.
    addTorsionForce();

    // Add Improper Torsion Force.
    addImproperTorsionForce();

    // Add Pi-Torsion Force.
    addPiTorsionForce();

    // Add Torsion-Torsion Force.
    addTorsionTorsionForce();

    // Add coordinate restraints.
    addRestrainPositionForce();

    // Add bond restraints.
    addRestraintBondForce();

    // Add Restrain Groups.
    addRestrainGroupsForce();

    // Add stretch-torsion coupling terms.
    addStretchTorsionForce();

    // Add angle-torsion coupling terms.
    addAngleTorsionForce();

    setDefaultPeriodicBoxVectors();

    addRestraintTorsions();

    if (vdW != null) {
      logger.info("\n Non-Bonded Terms\n");
      VanDerWaalsForm vdwForm = vdW.getVDWForm();
      if (vdwForm.vdwType == LENNARD_JONES) {
        addFixedChargeNonBondedForce();
      } else {
        // Add vdW Force.
        addAmoebaVDWForce();

        // Add Multipole Force.
        addAmoebaMultipoleForce();
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
          openMMContext.reinitContext();
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
        updateBondForce();
      }
      if (angleForce != null) {
        updateAngleForce();
      }
      if (stretchBendForce != null) {
        updateStretchBendForce();
      }
      if (inPlaneAngleForce != null) {
        updateInPlaneAngleForce();
      }
      if (ureyBradleyForce != null) {
        updateUreyBradleyForce();
      }
      if (outOfPlaneBendForce != null) {
        updateOutOfPlaneBendForce();
      }
      if (piTorsionForce != null) {
        updatePiTorsionForce();
      }
    }

    if (torsionLambdaTerm || updateBondedTerms) {
      if (torsionForce != null) {
        updateTorsionForce();
      }
      if (improperTorsionForce != null) {
        updateImproperTorsionForce();
      }
    }

    if (restraintTorsions != null && restraintTorsions.length > 0) {
      updateRestraintTorsions();
    }

    if (atoms == null || atoms.length == 0) {
      return;
    }

    // Update fixed charge non-bonded parameters.
    if (fixedChargeNonBondedForce != null) {
      updateFixedChargeNonBondedForce(atoms);
    }

    // Update fixed charge GB parameters.
    if (customGBForce != null) {
      updateCustomGBForce(atoms);
    }

    // Update AMOEBA vdW parameters.
    if (amoebaVDWForce != null) {
      updateAmoebaVDWForce(atoms);
    }

    // Update AMOEBA polarizable multipole parameters.
    if (amoebaMultipoleForce != null) {
      updateAmoebaMultipoleForce(atoms);
    }

    // Update GK force.
    if (amoebaGeneralizedKirkwoodForce != null) {
      updateGeneralizedKirkwoodForce(atoms);
    }

    // Update WCA Force.
    if (amoebaWcaDispersionForce != null) {
      updateWCAForce(atoms);
    }

    // Update WCA Force.
    if (amoebaCavitationForce != null) {
      updateCavitationForce(atoms);
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
   * Add a bond force to the OpenMM System.
   */
  private void addBondForce() {
    Bond[] bonds = openMMEnergy.getBonds();
    if (bonds == null || bonds.length < 1) {
      return;
    }

    String energy;
    BondType bondType = bonds[0].getBondType();
    if (bondType.bondFunction == BondType.BondFunction.QUARTIC) {
      energy = format("k*(d^2 + %.15g*d^3 + %.15g*d^4); d=r-r0",
          bondType.cubic / OpenMM_NmPerAngstrom,
          bondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
    } else {
      energy = "k*(d^2); d=r-r0";
    }
    bondForce = OpenMM_CustomBondForce_create(energy);
    OpenMM_CustomBondForce_addPerBondParameter(bondForce, "r0");
    OpenMM_CustomBondForce_addPerBondParameter(bondForce, "k");
    OpenMM_Force_setName(bondForce, "AmoebaBond");

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    for (Bond bond : bonds) {
      int i1 = bond.getAtom(0).getXyzIndex() - 1;
      int i2 = bond.getAtom(1).getXyzIndex() - 1;
      bondType = bond.bondType;
      double r0 = bondType.distance * OpenMM_NmPerAngstrom;
      double k = kParameterConversion * bondType.forceConstant * bond.bondType.bondUnit;
      OpenMM_DoubleArray_append(parameters, r0);
      OpenMM_DoubleArray_append(parameters, k);
      OpenMM_CustomBondForce_addBond(bondForce, i1, i2, parameters);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_DoubleArray_destroy(parameters);

    int forceGroup = forceField.getInteger("BOND_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(bondForce, forceGroup);
    OpenMM_System_addForce(system, bondForce);
    logger.log(Level.INFO, format("  Bonds \t\t%6d\t\t%1d", bonds.length, forceGroup));
  }

  /**
   * Update an existing bond force for the OpenMM System.
   */
  private void updateBondForce() {
    Bond[] bonds = openMMEnergy.getBonds();
    if (bonds == null || bonds.length < 1) {
      return;
    }

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int index = 0;
    for (Bond bond : bonds) {
      int i1 = bond.getAtom(0).getXyzIndex() - 1;
      int i2 = bond.getAtom(1).getXyzIndex() - 1;
      BondType bondType = bond.bondType;
      double r0 = bondType.distance * OpenMM_NmPerAngstrom;
      double k = kParameterConversion * bondType.forceConstant * bondType.bondUnit;
      OpenMM_DoubleArray_append(parameters, r0);
      OpenMM_DoubleArray_append(parameters, k);
      OpenMM_CustomBondForce_setBondParameters(bondForce, index++, i1, i2, parameters);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_DoubleArray_destroy(parameters);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomBondForce_updateParametersInContext(bondForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add an angle force to the OpenMM System.
   */
  private void addAngleForce() {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return;
    }
    boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    String energy;
    AngleType angleType = angles[0].angleType;
    if (angleType.angleFunction == AngleType.AngleFunction.SEXTIC) {
      energy = format(
          "k*(d^2 + %.15g*d^3 + %.15g*d^4 + %.15g*d^5 + %.15g*d^6); d=%.15g*theta-theta0",
          angleType.cubic, angleType.quartic, angleType.pentic, angleType.sextic, 180.0 / PI);
    } else {
      energy = format("k*(d^2); d=%.15g*theta-theta0", 180.0 / PI);
    }
    angleForce = OpenMM_CustomAngleForce_create(energy);
    OpenMM_CustomAngleForce_addPerAngleParameter(angleForce, "theta0");
    OpenMM_CustomAngleForce_addPerAngleParameter(angleForce, "k");
    OpenMM_Force_setName(angleForce, "Angle");

    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int angleCount = 0;
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;

      if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
        // Skip In-Plane angles unless this is ManyBody Titration.
      } else if (isHydrogenAngle(angle) && rigidHydrogenAngles) {
        logger.log(Level.INFO, " Constrained angle %s was not added the AngleForce.", angle);
      } else {
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        if (angleMode == AngleType.AngleMode.IN_PLANE) {
          // This is a place-holder Angle, in case the In-Plane Angle is swtiched to a
          // Normal Angle during in the udpateAngleForce.
          k = 0.0;
        }
        OpenMM_DoubleArray_append(parameters, theta0);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomAngleForce_addAngle(angleForce, i1, i2, i3, parameters);
        angleCount++;
        OpenMM_DoubleArray_resize(parameters, 0);
      }
    }
    OpenMM_DoubleArray_destroy(parameters);

    int forceGroup = forceField.getInteger("ANGLE_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(angleForce, forceGroup);
    OpenMM_System_addForce(system, angleForce);
    logger.log(Level.INFO, format("  Angles \t\t%6d\t\t%1d", angleCount, forceGroup));
  }

  /**
   * Update the angle force.
   */
  private void updateAngleForce() {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return;
    }
    boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int index = 0;
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;
      if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
        // Skip In-Plane angles unless this is ManyBody Titration.
      }
      // Update angles that do not involve rigid hydrogen atoms.
      else if (!rigidHydrogenAngles || !isHydrogenAngle(angle)) {
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        if (angleMode == AngleType.AngleMode.IN_PLANE) {
          // Zero the force constant for In-Plane Angles.
          k = 0.0;
        }
        OpenMM_DoubleArray_append(parameters, theta0);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomAngleForce_setAngleParameters(angleForce, index++, i1, i2, i3, parameters);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
    }
    OpenMM_DoubleArray_destroy(parameters);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomAngleForce_updateParametersInContext(angleForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add an in-plane angle force to the OpenMM System.
   */
  private void addInPlaneAngleForce() {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return;
    }

    AngleType angleType = angles[0].angleType;
    String energy = format("""
        k*(d^2 + %.15g*d^3 + %.15g*d^4 + %.15g*d^5 + %.15g*d^6);
        d=theta-theta0;
        theta = %.15g*pointangle(x1, y1, z1, projx, projy, projz, x3, y3, z3);
        projx = x2-nx*dot;
        projy = y2-ny*dot;
        projz = z2-nz*dot;
        dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3);
        nx = px/norm;
        ny = py/norm;
        nz = pz/norm;
        norm = sqrt(px*px + py*py + pz*pz);
        px = (d1y*d2z-d1z*d2y);
        py = (d1z*d2x-d1x*d2z);
        pz = (d1x*d2y-d1y*d2x);
        d1x = x1-x4;
        d1y = y1-y4;
        d1z = z1-z4;
        d2x = x3-x4;
        d2y = y3-y4;
        d2z = z3-z4;
        """, angleType.cubic, angleType.quartic, angleType.pentic, angleType.sextic, 180.0 / PI);
    inPlaneAngleForce = OpenMM_CustomCompoundBondForce_create(4, energy);
    OpenMM_CustomCompoundBondForce_addPerBondParameter(inPlaneAngleForce, "theta0");
    OpenMM_CustomCompoundBondForce_addPerBondParameter(inPlaneAngleForce, "k");
    OpenMM_Force_setName(inPlaneAngleForce, "InPlaneAngle");

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;

      if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
        // Skip Normal angles unless this is ManyBody Titration.
      } else {
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        int i4 = 0;
        if (angleMode == AngleType.AngleMode.NORMAL) {
          // This is a place-holder Angle, in case the Normal Angle is switched to an
          // In-Plane Angle during in the udpateInPlaneAngleForce.
          k = 0.0;
          Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
          if (fourthAtom != null) {
            i4 = fourthAtom.getXyzIndex() - 1;
          } else {
            while (i1 == i4 || i2 == i4 || i3 == i4) {
              i4++;
            }
          }
        } else {
          i4 = angle.getAtom4().getXyzIndex() - 1;
        }
        OpenMM_IntArray_append(particles, i1);
        OpenMM_IntArray_append(particles, i2);
        OpenMM_IntArray_append(particles, i3);
        OpenMM_IntArray_append(particles, i4);
        OpenMM_DoubleArray_append(parameters, theta0);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomCompoundBondForce_addBond(inPlaneAngleForce, particles, parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    int forceGroup = forceField.getInteger("IN_PLANE_ANGLE_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(inPlaneAngleForce, forceGroup);
    OpenMM_System_addForce(system, inPlaneAngleForce);
    logger.log(Level.INFO, format("  In-Plane Angles \t%6d\t\t%1d", angles.length, forceGroup));
  }

  /**
   * Update the in-plane angle force.
   */
  private void updateInPlaneAngleForce() {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return;
    }
    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int index = 0;
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;
      if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
        // Skip Normal angles unless this is ManyBody Titration.
      } else {
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        // There is no 4th atom for normal angles, so set the index to first atom.
        int i4 = 0;
        if (angleMode == AngleType.AngleMode.NORMAL) {
          // Zero the force constant for Normal Angles.
          k = 0.0;
          Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
          if (fourthAtom != null) {
            i4 = fourthAtom.getXyzIndex() - 1;
          } else {
            while (i1 == i4 || i2 == i4 || i3 == i4) {
              i4++;
            }
          }
        } else {
          i4 = angle.getAtom4().getXyzIndex() - 1;
        }
        OpenMM_IntArray_append(particles, i1);
        OpenMM_IntArray_append(particles, i2);
        OpenMM_IntArray_append(particles, i3);
        OpenMM_IntArray_append(particles, i4);
        OpenMM_DoubleArray_append(parameters, theta0);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomCompoundBondForce_setBondParameters(inPlaneAngleForce, index++, particles,
            parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomCompoundBondForce_updateParametersInContext(inPlaneAngleForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add the Urey-Bradley force to the OpenMM System.
   */
  private void addUreyBradleyForce() {
    UreyBradley[] ureyBradleys = openMMEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return;
    }

    ureyBradleyForce = OpenMM_HarmonicBondForce_create();
    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    for (UreyBradley ureyBradley : ureyBradleys) {
      int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
      int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
      UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
      double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
      // The implementation of UreyBradley in FFX & Tinker: k x^2
      // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
      double k =
          2.0 * ureyBradleyType.forceConstant * ureyBradleyType.ureyUnit * kParameterConversion;
      OpenMM_HarmonicBondForce_addBond(ureyBradleyForce, i1, i2, length, k);
    }

    int forceGroup = forceField.getInteger("UREY_BRADLEY_FORCE", 0);
    OpenMM_Force_setForceGroup(ureyBradleyForce, forceGroup);
    OpenMM_System_addForce(system, ureyBradleyForce);
    logger.log(Level.INFO,
        format("  Urey-Bradleys \t%6d\t\t%1d", ureyBradleys.length, forceGroup));
  }

  /**
   * Update the Urey-Bradley force.
   */
  private void updateUreyBradleyForce() {
    UreyBradley[] ureyBradleys = openMMEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return;
    }

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    int index = 0;
    for (UreyBradley ureyBradley : ureyBradleys) {
      int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
      int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
      UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
      double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
      // The implementation of UreyBradley in FFX & Tinker: k x^2
      // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
      double k =
          2.0 * ureyBradleyType.forceConstant * ureyBradleyType.ureyUnit * kParameterConversion;
      OpenMM_HarmonicBondForce_setBondParameters(ureyBradleyForce, index++, i1, i2, length, k);
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_HarmonicBondForce_updateParametersInContext(ureyBradleyForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add an out-of-plane bend force to the OpenMM System.
   */
  private void addOutOfPlaneBendForce() {
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return;
    }

    OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBends[0].outOfPlaneBendType;
    String energy = format(
        "k*(theta^2 + %.15g*theta^3 + %.15g*theta^4 + %.15g*theta^5 + %.15g*theta^6); "
            + "theta = %.15g*pointangle(x2, y2, z2, x4, y4, z4, projx, projy, projz); "
            + "projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot; "
            + "dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3); "
            + "nx = px/norm; ny = py/norm; nz = pz/norm; " + "norm = sqrt(px*px + py*py + pz*pz); "
            + "px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x); "
            + "d1x = x1-x4; d1y = y1-y4; d1z = z1-z4; " + "d2x = x3-x4; d2y = y3-y4; d2z = z3-z4",
        outOfPlaneBendType.cubic, outOfPlaneBendType.quartic, outOfPlaneBendType.pentic,
        outOfPlaneBendType.sextic, 180.0 / PI);
    outOfPlaneBendForce = OpenMM_CustomCompoundBondForce_create(4, energy);
    OpenMM_CustomCompoundBondForce_addPerBondParameter(outOfPlaneBendForce, "k");
    OpenMM_Force_setName(outOfPlaneBendForce, "OutOfPlaneBend");

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
      outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
      int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
      int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
      int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
      int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
      double k =
          OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * outOfPlaneBendType.opBendUnit;
      OpenMM_IntArray_append(particles, i1);
      OpenMM_IntArray_append(particles, i2);
      OpenMM_IntArray_append(particles, i3);
      OpenMM_IntArray_append(particles, i4);
      OpenMM_DoubleArray_append(parameters, k);
      OpenMM_CustomCompoundBondForce_addBond(outOfPlaneBendForce, particles, parameters);
      OpenMM_IntArray_resize(particles, 0);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);
    int forceGroup = forceField.getInteger("OUT_OF_PLANE_BEND_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(outOfPlaneBendForce, forceGroup);
    OpenMM_System_addForce(system, outOfPlaneBendForce);
    logger.log(Level.INFO,
        format("  Out-of-Plane Bends \t%6d\t\t%1d", outOfPlaneBends.length, forceGroup));
  }

  /**
   * Update the Out-of-Plane bend force.
   */
  private void updateOutOfPlaneBendForce() {
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return;
    }

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int index = 0;
    for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
      OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
      int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
      int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
      int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
      int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
      double k =
          OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * outOfPlaneBendType.opBendUnit;
      OpenMM_IntArray_append(particles, i1);
      OpenMM_IntArray_append(particles, i2);
      OpenMM_IntArray_append(particles, i3);
      OpenMM_IntArray_append(particles, i4);
      OpenMM_DoubleArray_append(parameters, k);
      OpenMM_CustomCompoundBondForce_setBondParameters(outOfPlaneBendForce, index++, particles,
          parameters);
      OpenMM_IntArray_resize(particles, 0);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomCompoundBondForce_updateParametersInContext(outOfPlaneBendForce, openMMContext.getContextPointer());
    }

  }

  /**
   * Add a stretch-bend force to the OpenMM System.
   */
  private void addStretchBendForce() {
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return;
    }

    String energy = format(
        "(k1*(distance(p1,p2)-r12) + k2*(distance(p2,p3)-r23))*(%.15g*(angle(p1,p2,p3)-theta0))",
        180.0 / PI);
    stretchBendForce = OpenMM_CustomCompoundBondForce_create(3, energy);
    OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "r12");
    OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "r23");
    OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "theta0");
    OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "k1");
    OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "k2");
    OpenMM_Force_setName(stretchBendForce, "AmoebaStretchBend");

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    for (StretchBend stretchBend : stretchBends) {
      int i1 = stretchBend.getAtom(0).getXyzIndex() - 1;
      int i2 = stretchBend.getAtom(1).getXyzIndex() - 1;
      int i3 = stretchBend.getAtom(2).getXyzIndex() - 1;
      double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
      double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
      double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
      double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      OpenMM_IntArray_append(particles, i1);
      OpenMM_IntArray_append(particles, i2);
      OpenMM_IntArray_append(particles, i3);
      OpenMM_DoubleArray_append(parameters, r12);
      OpenMM_DoubleArray_append(parameters, r23);
      OpenMM_DoubleArray_append(parameters, theta0);
      OpenMM_DoubleArray_append(parameters, k1);
      OpenMM_DoubleArray_append(parameters, k2);
      OpenMM_CustomCompoundBondForce_addBond(stretchBendForce, particles, parameters);
      OpenMM_IntArray_resize(particles, 0);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    int forceGroup = forceField.getInteger("STRETCH_BEND_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(stretchBendForce, forceGroup);
    OpenMM_System_addForce(system, stretchBendForce);
    logger.log(Level.INFO,
        format("  Stretch-Bends \t%6d\t\t%1d", stretchBends.length, forceGroup));
  }

  /**
   * Update the Stretch-Bend force.
   */
  private void updateStretchBendForce() {
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return;
    }

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int index = 0;
    for (StretchBend stretchBend : stretchBends) {
      int i1 = stretchBend.getAtom(0).getXyzIndex() - 1;
      int i2 = stretchBend.getAtom(1).getXyzIndex() - 1;
      int i3 = stretchBend.getAtom(2).getXyzIndex() - 1;
      double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
      double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
      double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
      double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      OpenMM_IntArray_append(particles, i1);
      OpenMM_IntArray_append(particles, i2);
      OpenMM_IntArray_append(particles, i3);
      OpenMM_DoubleArray_append(parameters, r12);
      OpenMM_DoubleArray_append(parameters, r23);
      OpenMM_DoubleArray_append(parameters, theta0);
      OpenMM_DoubleArray_append(parameters, k1);
      OpenMM_DoubleArray_append(parameters, k2);
      OpenMM_CustomCompoundBondForce_setBondParameters(stretchBendForce, index++, particles,
          parameters);
      OpenMM_IntArray_resize(particles, 0);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomCompoundBondForce_updateParametersInContext(stretchBendForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add a torsion force to the OpenMM System.
   */
  private void addTorsionForce() {
    Torsion[] torsions = openMMEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return;
    }

    torsionForce = OpenMM_PeriodicTorsionForce_create();
    for (Torsion torsion : torsions) {
      int a1 = torsion.getAtom(0).getXyzIndex() - 1;
      int a2 = torsion.getAtom(1).getXyzIndex() - 1;
      int a3 = torsion.getAtom(2).getXyzIndex() - 1;
      int a4 = torsion.getAtom(3).getXyzIndex() - 1;
      TorsionType torsionType = torsion.torsionType;
      int nTerms = torsionType.phase.length;
      for (int j = 0; j < nTerms; j++) {
        OpenMM_PeriodicTorsionForce_addTorsion(torsionForce, a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree,
            OpenMM_KJPerKcal * torsionType.torsionUnit * torsionType.amplitude[j]);
      }
      // Enforce 6-fold torsions since TorsionType instances can have different lengths
      // when side-chain protonation changes.
      if (manyBodyTitration) {
        for (int j = nTerms; j < 6; j++) {
          OpenMM_PeriodicTorsionForce_addTorsion(torsionForce, a1, a2, a3, a4, j + 1, 0.0, 0.0);
        }
      }
    }
    int fGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);

    OpenMM_Force_setForceGroup(torsionForce, fGroup);
    OpenMM_System_addForce(system, torsionForce);

    logger.log(Level.INFO, format("  Torsions \t\t%6d\t\t%1d", torsions.length, fGroup));
  }

  /**
   * Update the Torsion force.
   */
  private void updateTorsionForce() {
    // Check if this system has torsions.
    Torsion[] torsions = openMMEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return;
    }

    int index = 0;
    for (Torsion torsion : torsions) {
      TorsionType torsionType = torsion.torsionType;
      int nTerms = torsionType.phase.length;
      int a1 = torsion.getAtom(0).getXyzIndex() - 1;
      int a2 = torsion.getAtom(1).getXyzIndex() - 1;
      int a3 = torsion.getAtom(2).getXyzIndex() - 1;
      int a4 = torsion.getAtom(3).getXyzIndex() - 1;
      for (int j = 0; j < nTerms; j++) {
        double forceConstant =
            OpenMM_KJPerKcal * torsionType.torsionUnit * torsionType.amplitude[j] * lambdaTorsion;
        OpenMM_PeriodicTorsionForce_setTorsionParameters(torsionForce, index++, a1, a2, a3, a4,
            j + 1, torsionType.phase[j] * OpenMM_RadiansPerDegree, forceConstant);
      }
      // Enforce 6-fold torsions since TorsionType instances can have different lengths
      // when side-chain protonation changes.
      if (manyBodyTitration) {
        for (int j = nTerms; j < 6; j++) {
          OpenMM_PeriodicTorsionForce_setTorsionParameters(torsionForce, index++, a1, a2, a3, a4,
              j + 1, 0.0, 0.0);
        }
      }
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_PeriodicTorsionForce_updateParametersInContext(torsionForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add an improper-torsion force to the OpenMM System.
   */
  private void addImproperTorsionForce() {
    ImproperTorsion[] improperTorsions = openMMEnergy.getImproperTorsions();
    if (improperTorsions == null || improperTorsions.length < 1) {
      return;
    }

    improperTorsionForce = OpenMM_PeriodicTorsionForce_create();
    for (ImproperTorsion improperTorsion : improperTorsions) {
      int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
      ImproperTorsionType improperTorsionType = improperTorsion.improperType;
      OpenMM_PeriodicTorsionForce_addTorsion(improperTorsionForce, a1, a2, a3, a4,
          improperTorsionType.periodicity, improperTorsionType.phase * OpenMM_RadiansPerDegree,
          OpenMM_KJPerKcal * improperTorsion.improperType.impTorUnit * improperTorsion.scaleFactor
              * improperTorsionType.k);
    }

    int forceGroup = forceField.getInteger("IMPROPER_TORSION_FORCE_GROUP", 0);

    OpenMM_Force_setForceGroup(improperTorsionForce, forceGroup);
    OpenMM_System_addForce(system, improperTorsionForce);

    logger.log(Level.INFO,
        format("  Improper Torsions \t%6d\t\t%1d", improperTorsions.length, forceGroup));
  }

  /**
   * Update the Improper Torsion force.
   */
  private void updateImproperTorsionForce() {
    ImproperTorsion[] improperTorsions = openMMEnergy.getImproperTorsions();
    if (improperTorsions == null || improperTorsions.length < 1) {
      return;
    }

    int nImproperTorsions = improperTorsions.length;
    for (int i = 0; i < nImproperTorsions; i++) {
      ImproperTorsion improperTorsion = improperTorsions[i];
      int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
      ImproperTorsionType improperTorsionType = improperTorsion.improperType;
      double forceConstant =
          OpenMM_KJPerKcal * improperTorsion.improperType.impTorUnit * improperTorsion.scaleFactor
              * improperTorsionType.k * lambdaTorsion;
      OpenMM_PeriodicTorsionForce_setTorsionParameters(improperTorsionForce, i, a1, a2, a3, a4,
          improperTorsionType.periodicity, improperTorsionType.phase * OpenMM_RadiansPerDegree,
          forceConstant);
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_PeriodicTorsionForce_updateParametersInContext(improperTorsionForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add a Pi-Torsion force to the OpenMM System.
   */
  private void addPiTorsionForce() {
    PiOrbitalTorsion[] piOrbitalTorsions = openMMEnergy.getPiOrbitalTorsions();
    if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
      return;
    }

    String energy = "2*k*sin(phi)^2;"
        + "phi = pointdihedral(x3+c1x, y3+c1y, z3+c1z, x3, y3, z3, x4, y4, z4, x4+c2x, y4+c2y, z4+c2z); "
        + "c1x = (d14y*d24z-d14z*d24y); c1y = (d14z*d24x-d14x*d24z); c1z = (d14x*d24y-d14y*d24x); "
        + "c2x = (d53y*d63z-d53z*d63y); c2y = (d53z*d63x-d53x*d63z); c2z = (d53x*d63y-d53y*d63x); "
        + "d14x = x1-x4; d14y = y1-y4; d14z = z1-z4; "
        + "d24x = x2-x4; d24y = y2-y4; d24z = z2-z4; "
        + "d53x = x5-x3; d53y = y5-y3; d53z = z5-z3; "
        + "d63x = x6-x3; d63y = y6-y3; d63z = z6-z3";
    piTorsionForce = OpenMM_CustomCompoundBondForce_create(6, energy);
    OpenMM_CustomCompoundBondForce_addPerBondParameter(piTorsionForce, "k");
    OpenMM_Force_setName(piTorsionForce, "PiTorsion");

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
      int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
      int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
      int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
      PiOrbitalTorsionType type = piOrbitalTorsion.piOrbitalTorsionType;
      double k =
          OpenMM_KJPerKcal * type.forceConstant * piOrbitalTorsion.piOrbitalTorsionType.piTorsUnit;
      OpenMM_IntArray_append(particles, a1);
      OpenMM_IntArray_append(particles, a2);
      OpenMM_IntArray_append(particles, a3);
      OpenMM_IntArray_append(particles, a4);
      OpenMM_IntArray_append(particles, a5);
      OpenMM_IntArray_append(particles, a6);
      OpenMM_DoubleArray_append(parameters, k);
      OpenMM_CustomCompoundBondForce_addBond(piTorsionForce, particles, parameters);
      OpenMM_IntArray_resize(particles, 0);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    int forceGroup = forceField.getInteger("PI_ORBITAL_TORSION_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(piTorsionForce, forceGroup);
    OpenMM_System_addForce(system, piTorsionForce);
    logger.log(Level.INFO,
        format("  Pi-Orbital Torsions  \t%6d\t\t%1d", piOrbitalTorsions.length, forceGroup));
  }

  /**
   * Update the Pi-Torsion force.
   */
  private void updatePiTorsionForce() {
    PiOrbitalTorsion[] piOrbitalTorsions = openMMEnergy.getPiOrbitalTorsions();
    if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
      return;
    }

    PointerByReference particles = OpenMM_IntArray_create(0);
    PointerByReference parameters = OpenMM_DoubleArray_create(0);
    int index = 0;
    for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
      int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
      int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
      int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
      PiOrbitalTorsionType type = piOrbitalTorsion.piOrbitalTorsionType;
      double k =
          OpenMM_KJPerKcal * type.forceConstant * piOrbitalTorsion.piOrbitalTorsionType.piTorsUnit;
      OpenMM_IntArray_append(particles, a1);
      OpenMM_IntArray_append(particles, a2);
      OpenMM_IntArray_append(particles, a3);
      OpenMM_IntArray_append(particles, a4);
      OpenMM_IntArray_append(particles, a5);
      OpenMM_IntArray_append(particles, a6);
      OpenMM_DoubleArray_append(parameters, k);
      OpenMM_CustomCompoundBondForce_setBondParameters(piTorsionForce, index++, particles,
          parameters);
      OpenMM_IntArray_resize(particles, 0);
      OpenMM_DoubleArray_resize(parameters, 0);
    }
    OpenMM_IntArray_destroy(particles);
    OpenMM_DoubleArray_destroy(parameters);

    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomCompoundBondForce_updateParametersInContext(piTorsionForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add a Torsion-Torsion force to the OpenMM System.
   */
  private void addTorsionTorsionForce() {
    TorsionTorsion[] torsionTorsions = openMMEnergy.getTorsionTorsions();
    if (torsionTorsions == null || torsionTorsions.length < 1) {
      return;
    }

    // Load the torsion-torsions.
    int nTypes = 0;
    LinkedHashMap<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();
    PointerByReference amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create();
    for (TorsionTorsion torsionTorsion : torsionTorsions) {
      int ia = torsionTorsion.getAtom(0).getXyzIndex() - 1;
      int ib = torsionTorsion.getAtom(1).getXyzIndex() - 1;
      int ic = torsionTorsion.getAtom(2).getXyzIndex() - 1;
      int id = torsionTorsion.getAtom(3).getXyzIndex() - 1;
      int ie = torsionTorsion.getAtom(4).getXyzIndex() - 1;

      TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
      String key = torsionTorsionType.getKey();

      // Check if the TorTor parameters have already been added to the Hash.
      int gridIndex = 0;
      if (torTorTypes.containsKey(key)) {

        // If the TorTor has been added, get its (ordered) index in the Hash.
        int index = 0;
        for (String entry : torTorTypes.keySet()) {
          if (entry.equalsIgnoreCase(key)) {
            gridIndex = index;
            break;
          } else {
            index++;
          }
        }
      } else {
        // Add the new TorTor.
        torTorTypes.put(key, torsionTorsionType);
        gridIndex = nTypes;
        nTypes++;
      }

      Atom atom = torsionTorsion.getChiralAtom();
      int iChiral = -1;
      if (atom != null) {
        iChiral = atom.getXyzIndex() - 1;
      }
      OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(amoebaTorsionTorsionForce, ia, ib, ic, id,
          ie, iChiral, gridIndex);
    }

    // Load the Torsion-Torsion parameters.
    PointerByReference values = OpenMM_DoubleArray_create(6);
    int gridIndex = 0;
    for (String key : torTorTypes.keySet()) {
      TorsionTorsionType torTorType = torTorTypes.get(key);
      int nx = torTorType.nx;
      int ny = torTorType.ny;
      double[] tx = torTorType.tx;
      double[] ty = torTorType.ty;
      double[] f = torTorType.energy;
      double[] dx = torTorType.dx;
      double[] dy = torTorType.dy;
      double[] dxy = torTorType.dxy;

      // Create the 3D grid.
      PointerByReference grid3D = OpenMM_3D_DoubleArray_create(nx, ny, 6);
      int xIndex = 0;
      int yIndex = 0;
      for (int j = 0; j < nx * ny; j++) {
        int addIndex = 0;
        OpenMM_DoubleArray_set(values, addIndex++, tx[xIndex]);
        OpenMM_DoubleArray_set(values, addIndex++, ty[yIndex]);
        OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * f[j]);
        OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dx[j]);
        OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dy[j]);
        OpenMM_DoubleArray_set(values, addIndex, OpenMM_KJPerKcal * dxy[j]);
        OpenMM_3D_DoubleArray_set(grid3D, yIndex, xIndex, values);
        xIndex++;
        if (xIndex == nx) {
          xIndex = 0;
          yIndex++;
        }
      }
      OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(amoebaTorsionTorsionForce,
          gridIndex++, grid3D);
      OpenMM_3D_DoubleArray_destroy(grid3D);
    }
    OpenMM_DoubleArray_destroy(values);

    int forceGroup = forceField.getInteger("TORSION_TORSION_FORCE_GROUP", 0);

    OpenMM_Force_setForceGroup(amoebaTorsionTorsionForce, forceGroup);
    OpenMM_System_addForce(system, amoebaTorsionTorsionForce);
    logger.log(Level.INFO,
        format("  Torsion-Torsions  \t%6d\t\t%1d", torsionTorsions.length, forceGroup));
  }

  /**
   * Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic acid model).
   */
  private void addStretchTorsionForce() {
    StretchTorsion[] stretchTorsions = openMMEnergy.getStretchTorsions();
    if (stretchTorsions == null || stretchTorsions.length < 1) {
      return;
    }

    PointerByReference stretchTorsionForce = OpenMM_CustomCompoundBondForce_create(4,
        StretchTorsion.stretchTorsionForm());
    OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi1", 0);
    OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi2", Math.PI);
    OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi3", 0);

    for (int m = 1; m < 4; m++) {
      for (int n = 1; n < 4; n++) {
        OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchTorsionForce,
            String.format("k%d%d", m, n));
      }
    }

    for (int m = 1; m < 4; m++) {
      OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchTorsionForce,
          String.format("b%d", m));
    }

    final double unitConv = OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;

    for (StretchTorsion strTors : stretchTorsions) {
      double[] constants = strTors.getConstants();
      PointerByReference strTorsParams = OpenMM_DoubleArray_create(0);
      for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
          int index = (3 * m) + n;
          double kmn = constants[index] * unitConv;
          OpenMM_DoubleArray_append(strTorsParams, kmn);
        }
      }

      OpenMM_DoubleArray_append(strTorsParams, strTors.bondType1.distance * OpenMM_NmPerAngstrom);
      OpenMM_DoubleArray_append(strTorsParams, strTors.bondType2.distance * OpenMM_NmPerAngstrom);
      OpenMM_DoubleArray_append(strTorsParams, strTors.bondType3.distance * OpenMM_NmPerAngstrom);

      PointerByReference strTorsParticles = OpenMM_IntArray_create(0);
      Atom[] atoms = strTors.getAtomArray(true);
      for (int i = 0; i < 4; i++) {
        OpenMM_IntArray_append(strTorsParticles, atoms[i].getXyzIndex() - 1);
      }

      OpenMM_CustomCompoundBondForce_addBond(stretchTorsionForce, strTorsParticles, strTorsParams);
      OpenMM_DoubleArray_destroy(strTorsParams);
      OpenMM_IntArray_destroy(strTorsParticles);
    }

    int forceGroup = forceField.getInteger("STRETCH_TORSION_FORCE_GROUP", 0);

    OpenMM_Force_setForceGroup(stretchTorsionForce, forceGroup);
    OpenMM_System_addForce(system, stretchTorsionForce);

    logger.log(Level.INFO,
        format("  Stretch-Torsions  \t%6d\t\t%1d", stretchTorsions.length, forceGroup));
  }

  /**
   * Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic acid model).
   */
  private void addAngleTorsionForce() {
    AngleTorsion[] angleTorsions = openMMEnergy.getAngleTorsions();
    if (angleTorsions == null || angleTorsions.length < 1) {
      return;
    }

    PointerByReference angleTorsionForce = OpenMM_CustomCompoundBondForce_create(4,
        AngleTorsion.angleTorsionForm());
    OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi1", 0);
    OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi2", Math.PI);
    OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi3", 0);

    for (int m = 1; m < 3; m++) {
      for (int n = 1; n < 4; n++) {
        OpenMM_CustomCompoundBondForce_addPerBondParameter(angleTorsionForce,
            format("k%d%d", m, n));
      }
    }

    for (int m = 1; m < 3; m++) {
      OpenMM_CustomCompoundBondForce_addPerBondParameter(angleTorsionForce, format("a%d", m));
    }

    for (AngleTorsion angleTorsion : angleTorsions) {
      double[] constants = angleTorsion.getConstants();
      PointerByReference atorsParams = OpenMM_DoubleArray_create(0);
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 3; n++) {
          int index = (3 * m) + n;
          double kmn = constants[index] * OpenMM_KJPerKcal;
          OpenMM_DoubleArray_append(atorsParams, kmn);
        }
      }

      Atom[] atoms = angleTorsion.getAtomArray(true);

      // One thing that concerns me is whether it's correct to get angle[0] instead of angle[num
      // hydrogens].
      // This is the way it is in FFX, but that may be a bug.

      OpenMM_DoubleArray_append(atorsParams,
          angleTorsion.angleType1.angle[0] * OpenMM_RadiansPerDegree);
      OpenMM_DoubleArray_append(atorsParams,
          angleTorsion.angleType2.angle[0] * OpenMM_RadiansPerDegree);

      PointerByReference atorsParticles = OpenMM_IntArray_create(0);
      for (int i = 0; i < 4; i++) {
        OpenMM_IntArray_append(atorsParticles, atoms[i].getXyzIndex() - 1);
      }

      OpenMM_CustomCompoundBondForce_addBond(angleTorsionForce, atorsParticles, atorsParams);
      OpenMM_DoubleArray_destroy(atorsParams);
      OpenMM_IntArray_destroy(atorsParticles);
    }

    int forceGroup = forceField.getInteger("ANGLE_TORSION_FORCE_GROUP", 0);

    OpenMM_Force_setForceGroup(angleTorsionForce, forceGroup);
    OpenMM_System_addForce(system, angleTorsionForce);

    logger.log(Level.INFO,
        format("  Angle-Torsions  \t%6d\t\t%1d", angleTorsions.length, forceGroup));
  }

  private void addRestraintTorsions() {
    List<RestraintTorsion> rTors = openMMEnergy.getRestraintTorsions();
    if (rTors != null && !rTors.isEmpty()) {
      int nRT = rTors.size();
      restraintTorsions = new PointerByReference[nRT];
      for (int i = 0; i < nRT; i++) {
        PointerByReference rtOMM = OpenMM_PeriodicTorsionForce_create();
        RestraintTorsion rt = rTors.get(i);
        int a1 = rt.getAtom(0).getXyzIndex() - 1;
        int a2 = rt.getAtom(1).getXyzIndex() - 1;
        int a3 = rt.getAtom(2).getXyzIndex() - 1;
        int a4 = rt.getAtom(3).getXyzIndex() - 1;
        int nTerms = rt.torsionType.terms;
        for (int j = 0; j < nTerms; j++) {
          OpenMM_PeriodicTorsionForce_addTorsion(rtOMM, a1, a2, a3, a4, j + 1,
              rt.torsionType.phase[j] * OpenMM_RadiansPerDegree,
              OpenMM_KJPerKcal * rt.units * rt.torsionType.amplitude[j]);
        }
        int fGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);

        OpenMM_Force_setForceGroup(rtOMM, fGroup);
        OpenMM_System_addForce(system, rtOMM);
        restraintTorsions[i] = rtOMM;
      }
      logger.info(format(" Added %d restraint torsions to OpenMM.", nRT));
    }
  }

  private void updateRestraintTorsions() {
    int nRT = restraintTorsions.length;
    if (nRT == 0) {
      return;
    }
    List<RestraintTorsion> rTors = openMMEnergy.getRestraintTorsions();
    double lambda = openMMEnergy.getLambda();
    for (int i = 0; i < nRT; i++) {
      RestraintTorsion rt = rTors.get(i);
      PointerByReference rtOMM = restraintTorsions[i];
      // Only update parameters if torsions are being scaled by lambda.
      if (rt.applyLambda()) {
        int index = 0;
        TorsionType torsionType = rt.torsionType;
        int nTerms = torsionType.phase.length;
        int a1 = rt.getAtom(0).getXyzIndex() - 1;
        int a2 = rt.getAtom(1).getXyzIndex() - 1;
        int a3 = rt.getAtom(2).getXyzIndex() - 1;
        int a4 = rt.getAtom(3).getXyzIndex() - 1;
        for (int j = 0; j < nTerms; j++) {
          double forceConstant =
              OpenMM_KJPerKcal * rt.units * torsionType.amplitude[j] * rt.mapLambda(lambda);
          OpenMM_PeriodicTorsionForce_setTorsionParameters(rtOMM, index++, a1, a2, a3, a4, j + 1,
              torsionType.phase[j] * OpenMM_RadiansPerDegree, forceConstant);
        }
      }

      if (openMMContext.hasContextPointer()) {
        OpenMM_PeriodicTorsionForce_updateParametersInContext(rtOMM, openMMContext.getContextPointer());
      }
    }
  }

  /**
   * Uses arithmetic mean to define sigma and geometric mean for epsilon.
   */
  private void addFixedChargeNonBondedForce() {
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

    fixedChargeNonBondedForce = OpenMM_NonbondedForce_create();

    // OpenMM vdW force requires a diameter (i.e. not radius).
    double radScale = 1.0;
    if (vdwForm.radiusSize == RADIUS) {
      radScale = 2.0;
    }

    // OpenMM vdw force requires atomic sigma values (i.e. not r-min).
    if (vdwForm.radiusType == R_MIN) {
      radScale /= 1.122462048309372981;
    }

    // Add particles.
    for (Atom atom : atoms) {
      VDWType vdwType = atom.getVDWType();
      double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
      double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
      double charge = 0.0;
      MultipoleType multipoleType = atom.getMultipoleType();
      if (multipoleType != null && atom.getElectrostatics()) {
        charge = multipoleType.charge;
      }
      OpenMM_NonbondedForce_addParticle(fixedChargeNonBondedForce, charge, sigma, eps);
    }

    // Define 1-4 scale factors.
    double lj14Scale = vdwForm.getScale14();
    double coulomb14Scale = 1.0 / 1.2;
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme != null) {
      coulomb14Scale = pme.getScale14();
    }
    Bond[] bonds = openMMEnergy.getBonds();
    if (bonds != null && bonds.length > 0) {
      PointerByReference bondArray = OpenMM_BondArray_create(0);
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getXyzIndex() - 1;
        int i2 = bond.getAtom(1).getXyzIndex() - 1;
        OpenMM_BondArray_append(bondArray, i1, i2);
      }
      OpenMM_NonbondedForce_createExceptionsFromBonds(fixedChargeNonBondedForce, bondArray,
          coulomb14Scale, lj14Scale);
      OpenMM_BondArray_destroy(bondArray);

      int num = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);
      chargeExclusion = new boolean[num];
      vdWExclusion = new boolean[num];
      exceptionChargeProd = new double[num];
      exceptionEps = new double[num];

      IntByReference particle1 = new IntByReference();
      IntByReference particle2 = new IntByReference();
      DoubleByReference chargeProd = new DoubleByReference();
      DoubleByReference sigma = new DoubleByReference();
      DoubleByReference eps = new DoubleByReference();

      for (int i = 0; i < num; i++) {
        OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i, particle1,
            particle2, chargeProd, sigma, eps);
        if (abs(chargeProd.getValue()) > 0.0) {
          chargeExclusion[i] = false;
          exceptionChargeProd[i] = chargeProd.getValue();
        } else {
          exceptionChargeProd[i] = 0.0;
          chargeExclusion[i] = true;
        }
        if (abs(eps.getValue()) > 0.0) {
          vdWExclusion[i] = false;
          exceptionEps[i] = eps.getValue();
        } else {
          vdWExclusion[i] = true;
          exceptionEps[i] = 0.0;
        }
      }
    }

    Crystal crystal = openMMEnergy.getCrystal();
    if (crystal.aperiodic()) {
      OpenMM_NonbondedForce_setNonbondedMethod(fixedChargeNonBondedForce,
          OpenMMLibrary.OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_NoCutoff);
    } else {
      OpenMM_NonbondedForce_setNonbondedMethod(fixedChargeNonBondedForce,
          OpenMMLibrary.OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_PME);

      if (pme != null) {
        // Units of the Ewald coefficient are A^-1; Multiply by AngstromsPerNM to convert to
        // (Nm^-1).
        double aEwald = OpenMM_AngstromsPerNm * pme.getEwaldCoefficient();
        int nx = pme.getReciprocalSpace().getXDim();
        int ny = pme.getReciprocalSpace().getYDim();
        int nz = pme.getReciprocalSpace().getZDim();
        OpenMM_NonbondedForce_setPMEParameters(fixedChargeNonBondedForce, aEwald, nx, ny, nz);
      }

      NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
      double off = nonbondedCutoff.off;
      double cut = nonbondedCutoff.cut;
      OpenMM_NonbondedForce_setCutoffDistance(fixedChargeNonBondedForce,
          OpenMM_NmPerAngstrom * off);
      OpenMM_NonbondedForce_setUseSwitchingFunction(fixedChargeNonBondedForce, OpenMM_True);
      if (cut == off) {
        logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
        if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
          logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
          cut = 1E40;
        } else {
          logger.info(String.format(
              " Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut,
              off));
          cut *= 0.99;
        }
      }
      OpenMM_NonbondedForce_setSwitchingDistance(fixedChargeNonBondedForce,
          OpenMM_NmPerAngstrom * cut);
    }

    OpenMM_NonbondedForce_setUseDispersionCorrection(fixedChargeNonBondedForce, OpenMM_False);

    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
    int pmeGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
    if (forceGroup != pmeGroup) {
      logger.severe(String.format(" ERROR: VDW-FORCE-GROUP is %d while PME-FORCE-GROUP is %d. "
              + "This is invalid for fixed-charge force fields with combined nonbonded forces.",
          forceGroup, pmeGroup));
    }

    OpenMM_Force_setForceGroup(fixedChargeNonBondedForce, forceGroup);
    OpenMM_System_addForce(system, fixedChargeNonBondedForce);

    logger.log(Level.INFO, format("  Fixed charge non-bonded force \t%1d", forceGroup));

    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk != null) {
      addCustomGBForce();
    }
  }

  /**
   * Updates the fixed-charge non-bonded force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateFixedChargeNonBondedForce(Atom[] atoms) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    // Only 6-12 LJ with arithmetic mean to define sigma and geometric mean for epsilon is
    // supported.
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    if (vdwForm.vdwType != LENNARD_JONES || vdwForm.radiusRule != ARITHMETIC
        || vdwForm.epsilonRule != GEOMETRIC) {
      logger.log(Level.SEVERE, " Unsupported van der Waals functional form.");
      return;
    }

    // OpenMM vdW force requires a diameter (i.e. not radius).
    double radScale = 1.0;
    if (vdwForm.radiusSize == RADIUS) {
      radScale = 2.0;
    }

    // OpenMM vdw force requires atomic sigma values (i.e. not r-min).
    if (vdwForm.radiusType == R_MIN) {
      radScale /= 1.122462048309372981;
    }

    // Update parameters.
    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      boolean applyLambda = atom.applyLambda();

      double charge = Double.MIN_VALUE;
      MultipoleType multipoleType = atom.getMultipoleType();
      if (multipoleType != null && atom.getElectrostatics()) {
        charge = multipoleType.charge;
      }

      VDWType vdwType = atom.getVDWType();
      double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
      double eps = OpenMM_KJPerKcal * vdwType.wellDepth;

      if (applyLambda) {
        // If we're using vdwLambdaTerm, this atom's vdW interactions are handled by the Custom
        // Non-Bonded force.
        if (vdwLambdaTerm) {
          eps = 0.0;
        }
        // Always scale the charge by lambdaElec
        charge *= lambdaElec;
      }

      if (!atom.getUse()) {
        eps = 0.0;
        charge = 0.0;
      }

      OpenMM_NonbondedForce_setParticleParameters(fixedChargeNonBondedForce, index, charge, sigma,
          eps);
    }

    // Update Exceptions.
    IntByReference particle1 = new IntByReference();
    IntByReference particle2 = new IntByReference();
    DoubleByReference chargeProd = new DoubleByReference();
    DoubleByReference sigma = new DoubleByReference();
    DoubleByReference eps = new DoubleByReference();

    int numExceptions = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);
    for (int i = 0; i < numExceptions; i++) {

      // Only update exceptions.
      if (chargeExclusion[i] && vdWExclusion[i]) {
        continue;
      }

      OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i, particle1,
          particle2, chargeProd, sigma, eps);

      int i1 = particle1.getValue();
      int i2 = particle2.getValue();

      double qq = exceptionChargeProd[i];
      double epsilon = exceptionEps[i];

      Atom atom1 = atoms[i1];
      Atom atom2 = atoms[i2];

      /*
      Note that the minimum epsilon value cannot be zero, or OpenMM may
      report an error that the number of Exceptions has changed.
      */
      double minEpsilon = 1.0e-12;
      double lambdaValue = lambdaElec;

      if (lambdaValue < minEpsilon) {
        lambdaValue = minEpsilon;
      }

      if (atom1.applyLambda()) {
        qq *= lambdaValue;
        if (vdwLambdaTerm) {
          epsilon = minEpsilon;
        }
      }
      if (atom2.applyLambda()) {
        qq *= lambdaValue;
        if (vdwLambdaTerm) {
          epsilon = minEpsilon;
        }
      }
      if (!atom1.getUse() || !atom2.getUse()) {
        qq = minEpsilon;
        epsilon = minEpsilon;
      }
      OpenMM_NonbondedForce_setExceptionParameters(fixedChargeNonBondedForce, i, i1, i2, qq,
          sigma.getValue(), epsilon);
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_NonbondedForce_updateParametersInContext(fixedChargeNonBondedForce, openMMContext.getContextPointer());
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

      OpenMM_NonbondedForce_getParticleParameters(fixedChargeNonBondedForce, index, charge, sigma,
          eps);
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

    OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup,
        alchemicalGroup);
    OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup,
        nonAlchemicalGroup);
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
        logger.info(
            format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut,
                off));
        cut *= 0.99;
      }
    }

    OpenMM_CustomNonbondedForce_setCutoffDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * off);
    OpenMM_CustomNonbondedForce_setUseSwitchingFunction(fixedChargeSoftcore, OpenMM_True);
    OpenMM_CustomNonbondedForce_setSwitchingDistance(fixedChargeSoftcore,
        OpenMM_NmPerAngstrom * cut);

    // Add energy parameter derivative
    // OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(fixedChargeSoftcore,
    // "vdw_lambda");

    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);

    OpenMM_Force_setForceGroup(fixedChargeSoftcore, forceGroup);
    OpenMM_System_addForce(system, fixedChargeSoftcore);

    // Alchemical with Alchemical could be either softcore or normal interactions (softcore here).
    PointerByReference alchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(
        stericsEnergyExpression);

    // Non-Alchemical with Alchemical is essentially always softcore.
    PointerByReference nonAlchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(
        stericsEnergyExpression);

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

    int range = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);

    IntByReference atomi = new IntByReference();
    IntByReference atomj = new IntByReference();
    int[][] torsionMask = vdW.getMask14();

    for (int i = 0; i < range; i++) {
      OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i, atomi, atomj,
          charge, sigma, eps);

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

  /**
   * Add an AMOEBA van der Waals force to the OpenMM System.
   */
  private void addAmoebaVDWForce() {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      return;
    }

    amoebaVDWForce = OpenMM_AmoebaVdwForce_create();

    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
    Crystal crystal = openMMEnergy.getCrystal();

    double radScale = 1.0;
    if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    /*
      The vdW class used to specify no vdW interactions for an atom will be Zero
      if all atom classes are greater than zero. Otherwise:
      vdWClassForNoInteraction = min(atomClass) - 1
     */
    vdWClassForNoInteraction = 0;
    // Add vdW parameters to the force and record their type.
    vdwClassToOpenMMType = new HashMap<>();

    Map<String, VDWType> vdwTypes = forceField.getVDWTypes();

    for (VDWType vdwType : vdwTypes.values()) {
      int atomClass = vdwType.atomClass;
      if (!vdwClassToOpenMMType.containsKey(atomClass)) {
        double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
        double rad = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
        int type = OpenMM_AmoebaVdwForce_addParticleType(amoebaVDWForce, rad, eps);
        vdwClassToOpenMMType.put(atomClass, type);
        if (atomClass <= vdWClassForNoInteraction) {
          vdWClassForNoInteraction = atomClass - 1;
        }
      }
    }

    // Add a special vdW type for zero vdW energy and forces (e.g. to support the FFX "use" flag).
    int type = OpenMM_AmoebaVdwForce_addParticleType(amoebaVDWForce, OpenMM_NmPerAngstrom, 0.0);
    vdwClassToOpenMMType.put(vdWClassForNoInteraction, type);

    Map<String, VDWPairType> vdwPairTypeMap = forceField.getVDWPairTypes();
    for (VDWPairType vdwPairType : vdwPairTypeMap.values()) {
      int c1 = vdwPairType.atomClasses[0];
      int c2 = vdwPairType.atomClasses[1];
      int type1 = vdwClassToOpenMMType.get(c1);
      int type2 = vdwClassToOpenMMType.get(c2);
      double rMin = vdwPairType.radius * OpenMM_NmPerAngstrom;
      double eps = vdwPairType.wellDepth * OpenMM_KJPerKcal;
      OpenMM_AmoebaVdwForce_addTypePair(amoebaVDWForce, type1, type2, rMin, eps);
      OpenMM_AmoebaVdwForce_addTypePair(amoebaVDWForce, type2, type1, rMin, eps);
    }

    ExtendedSystem extendedSystem = vdW.getExtendedSystem();
    double[] vdwPrefactorAndDerivs = new double[3];

    int[] ired = vdW.getReductionIndex();
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      VDWType vdwType = atom.getVDWType();
      int atomClass = vdwType.atomClass;
      type = vdwClassToOpenMMType.get(atomClass);
      int isAlchemical = atom.applyLambda() ? 1 : 0;
      double scaleFactor = 1.0;
      if (extendedSystem != null) {
        extendedSystem.getVdwPrefactor(i, vdwPrefactorAndDerivs);
        scaleFactor = vdwPrefactorAndDerivs[0];
      }

      OpenMM_AmoebaVdwForce_addParticle_1(amoebaVDWForce, ired[i], type, vdwType.reductionFactor,
          isAlchemical, scaleFactor);
    }

    double cutoff = nonbondedCutoff.off * OpenMM_NmPerAngstrom;
    OpenMM_AmoebaVdwForce_setCutoffDistance(amoebaVDWForce, cutoff);

    if (vdW.getDoLongRangeCorrection()) {
      OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_True);
    } else {
      OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_False);
    }

    if (crystal.aperiodic()) {
      OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce, OpenMM_AmoebaVdwForce_NoCutoff);
    } else {
      OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
          OpenMM_AmoebaVdwForce_CutoffPeriodic);
    }

    if (vdwLambdaTerm) {
      OpenMM_AmoebaVdwForce_setAlchemicalMethod(amoebaVDWForce, OpenMM_AmoebaVdwForce_Decouple);
      OpenMM_AmoebaVdwForce_setSoftcoreAlpha(amoebaVDWForce, vdWSoftcoreAlpha);
      OpenMM_AmoebaVdwForce_setSoftcorePower(amoebaVDWForce, (int) vdwSoftcorePower);
    }

    int[][] bondMask = vdW.getMask12();
    int[][] angleMask = vdW.getMask13();

    // Create exclusion lists.
    PointerByReference exclusions = OpenMM_IntArray_create(0);
    for (int i = 0; i < nAtoms; i++) {
      OpenMM_IntArray_append(exclusions, i);
      final int[] bondMaski = bondMask[i];
      for (int value : bondMaski) {
        OpenMM_IntArray_append(exclusions, value);
      }
      final int[] angleMaski = angleMask[i];
      for (int value : angleMaski) {
        OpenMM_IntArray_append(exclusions, value);
      }
      OpenMM_AmoebaVdwForce_setParticleExclusions(amoebaVDWForce, i, exclusions);
      OpenMM_IntArray_resize(exclusions, 0);
    }
    OpenMM_IntArray_destroy(exclusions);

    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
    OpenMM_Force_setForceGroup(amoebaVDWForce, forceGroup);
    OpenMM_System_addForce(system, amoebaVDWForce);

    logger.log(Level.INFO, format("  AMOEBA van der Waals force \t\t%d", forceGroup));
  }

  /**
   * Updates the AMOEBA van der Waals force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateAmoebaVDWForce(Atom[] atoms) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    double radScale = 1.0;
    if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    ExtendedSystem extendedSystem = vdW.getExtendedSystem();
    double[] vdwPrefactorAndDerivs = new double[3];

    int[] ired = vdW.getReductionIndex();
    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      VDWType vdwType = atom.getVDWType();

      // Get the OpenMM index for this vdW type.
      int type = vdwClassToOpenMMType.get(vdwType.atomClass);
      if (!atom.getUse()) {
        // Get the OpenMM index for a special vdW type that has no interactions.
        type = vdwClassToOpenMMType.get(vdWClassForNoInteraction);
      }
      int isAlchemical = atom.applyLambda() ? 1 : 0;
      double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
      double rad = OpenMM_NmPerAngstrom * vdwType.radius * radScale;

      double scaleFactor = 1.0;
      if (extendedSystem != null) {
        extendedSystem.getVdwPrefactor(index, vdwPrefactorAndDerivs);
        scaleFactor = vdwPrefactorAndDerivs[0];
      }

      OpenMM_AmoebaVdwForce_setParticleParameters(amoebaVDWForce, index, ired[index], rad, eps,
          vdwType.reductionFactor, isAlchemical, type, scaleFactor);
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_AmoebaVdwForce_updateParametersInContext(amoebaVDWForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add an AMOEBA polarizable multipole force to the OpenMM System.
   */
  private void addAmoebaMultipoleForce() {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme == null) {
      return;
    }

    int[][] axisAtom = pme.getAxisAtoms();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();

    double polarScale = 1.0;
    SCFAlgorithm scfAlgorithm = null;

    if (pme.getPolarizationType() != Polarization.MUTUAL) {
      OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce,
          OpenMM_AmoebaMultipoleForce_Direct);
      if (pme.getPolarizationType() == Polarization.NONE) {
        polarScale = 0.0;
      }
    } else {
      String algorithm = forceField.getString("SCF_ALGORITHM", "CG");
      try {
        algorithm = algorithm.replaceAll("-", "_").toUpperCase();
        scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
      } catch (Exception e) {
        scfAlgorithm = SCFAlgorithm.CG;
      }

      if (scfAlgorithm == SCFAlgorithm.EPT) {
        /*
         * Citation:
         * Simmonett, A. C.;  Pickard, F. C. t.;  Shao, Y.;  Cheatham, T. E., 3rd; Brooks, B. R., Efficient treatment of induced dipoles. The Journal of chemical physics 2015, 143 (7), 074115-074115.
         */
        OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce,
            OpenMM_AmoebaMultipoleForce_Extrapolated);
        PointerByReference exptCoefficients = OpenMM_DoubleArray_create(4);
        OpenMM_DoubleArray_set(exptCoefficients, 0, -0.154);
        OpenMM_DoubleArray_set(exptCoefficients, 1, 0.017);
        OpenMM_DoubleArray_set(exptCoefficients, 2, 0.657);
        OpenMM_DoubleArray_set(exptCoefficients, 3, 0.475);
        // PointerByReference exptCoefficients = OpenMM_DoubleArray_create(1);
        // OpenMM_DoubleArray_set(exptCoefficients, 0, 1.044);
        OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(amoebaMultipoleForce,
            exptCoefficients);
        OpenMM_DoubleArray_destroy(exptCoefficients);
      } else {
        OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce,
            OpenMM_AmoebaMultipoleForce_Mutual);
      }
    }

    PointerByReference dipoles = OpenMM_DoubleArray_create(3);
    PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      MultipoleType multipoleType = pme.getMultipoleType(i);
      PolarizeType polarType = pme.getPolarizeType(i);

      // Define the frame definition.
      int axisType = switch (multipoleType.frameDefinition) {
        case NONE -> OpenMM_AmoebaMultipoleForce_NoAxisType;
        case ZONLY -> OpenMM_AmoebaMultipoleForce_ZOnly;
        case ZTHENX -> OpenMM_AmoebaMultipoleForce_ZThenX;
        case BISECTOR -> OpenMM_AmoebaMultipoleForce_Bisector;
        case ZTHENBISECTOR -> OpenMM_AmoebaMultipoleForce_ZBisect;
        case THREEFOLD -> OpenMM_AmoebaMultipoleForce_ThreeFold;
      };

      double useFactor = 1.0;
      if (!atoms[i].getUse() || !atoms[i].getElectrostatics()) {
        useFactor = 0.0;
      }

      double lambdaScale = lambdaElec; // Should be 1.0 at this point.
      if (!atom.applyLambda()) {
        lambdaScale = 1.0;
      }

      useFactor *= lambdaScale;

      // Load local multipole coefficients.
      for (int j = 0; j < 3; j++) {
        OpenMM_DoubleArray_set(dipoles, j,
            multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          OpenMM_DoubleArray_set(quadrupoles, l++,
              multipoleType.quadrupole[j][k] * quadrupoleConversion * useFactor / 3.0);
        }
      }

      int zaxis = -1;
      int xaxis = -1;
      int yaxis = -1;
      int[] refAtoms = axisAtom[i];
      if (refAtoms != null) {
        zaxis = refAtoms[0];
        if (refAtoms.length > 1) {
          xaxis = refAtoms[1];
          if (refAtoms.length > 2) {
            yaxis = refAtoms[2];
          }
        }
      } else {
        axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
      }

      double charge = multipoleType.charge * useFactor;

      // Add the multipole.
      OpenMM_AmoebaMultipoleForce_addMultipole(amoebaMultipoleForce, charge, dipoles, quadrupoles,
          axisType, zaxis, xaxis, yaxis, polarType.thole,
          polarType.pdamp * dampingFactorConversion,
          polarType.polarizability * polarityConversion * polarScale);
    }
    OpenMM_DoubleArray_destroy(dipoles);
    OpenMM_DoubleArray_destroy(quadrupoles);

    Crystal crystal = openMMEnergy.getCrystal();
    if (!crystal.aperiodic()) {
      OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce,
          OpenMM_AmoebaMultipoleForce_PME);
      OpenMM_AmoebaMultipoleForce_setCutoffDistance(amoebaMultipoleForce,
          pme.getEwaldCutoff() * OpenMM_NmPerAngstrom);
      OpenMM_AmoebaMultipoleForce_setAEwald(amoebaMultipoleForce,
          pme.getEwaldCoefficient() / OpenMM_NmPerAngstrom);

      double ewaldTolerance = 1.0e-04;
      OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(amoebaMultipoleForce, ewaldTolerance);

      PointerByReference gridDimensions = OpenMM_IntArray_create(3);
      ReciprocalSpace recip = pme.getReciprocalSpace();
      OpenMM_IntArray_set(gridDimensions, 0, recip.getXDim());
      OpenMM_IntArray_set(gridDimensions, 1, recip.getYDim());
      OpenMM_IntArray_set(gridDimensions, 2, recip.getZDim());
      OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(amoebaMultipoleForce, gridDimensions);
      OpenMM_IntArray_destroy(gridDimensions);
    } else {
      OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce,
          OpenMM_AmoebaMultipoleForce_NoCutoff);
    }

    OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(amoebaMultipoleForce, 500);
    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(amoebaMultipoleForce,
        pme.getPolarEps());

    int[][] ip11 = pme.getPolarization11();

    PointerByReference covalentMap = OpenMM_IntArray_create(0);
    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];

      // 1-2 Mask
      OpenMM_IntArray_resize(covalentMap, 0);
      for (Atom ak : ai.get12List()) {
        OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
      }
      OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
          OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);

      // 1-3 Mask
      OpenMM_IntArray_resize(covalentMap, 0);
      for (Atom ak : ai.get13List()) {
        OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
      }
      OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
          OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);

      // 1-4 Mask
      OpenMM_IntArray_resize(covalentMap, 0);
      for (Atom ak : ai.get14List()) {
        OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
      }
      OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
          OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);

      // 1-5 Mask
      OpenMM_IntArray_resize(covalentMap, 0);
      for (Atom ak : ai.get15List()) {
        OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
      }
      OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
          OpenMM_AmoebaMultipoleForce_Covalent15, covalentMap);

      // 1-1 Polarization Groups.
      OpenMM_IntArray_resize(covalentMap, 0);
      for (int j = 0; j < ip11[i].length; j++) {
        OpenMM_IntArray_append(covalentMap, ip11[i][j]);
      }
      OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
          OpenMM_AmoebaMultipoleForce_PolarizationCovalent11, covalentMap);

      // AMOEBA does not scale between 1-2, 1-3, etc. polarization groups.
    }

    OpenMM_IntArray_destroy(covalentMap);

    int forceGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
    OpenMM_Force_setForceGroup(amoebaMultipoleForce, forceGroup);
    OpenMM_System_addForce(system, amoebaMultipoleForce);

    logger.log(Level.INFO, format("  AMOEBA polarizable multipole force \t%d", forceGroup));

    GeneralizedKirkwood gk = openMMEnergy.getGK();
    if (gk != null) {
      addGeneralizedKirkwoodForce();
    }

    if (scfAlgorithm == SCFAlgorithm.EPT) {
      logger.info("   Using extrapolated perturbation theory for polarization energy.");
    }
  }

  /**
   * Updates the Amoeba electrostatic multipolar force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateAmoebaMultipoleForce(Atom[] atoms) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    double polarScale = 1.0;
    if (pme.getPolarizationType() == Polarization.NONE) {
      polarScale = 0.0;
    }

    PointerByReference dipoles = OpenMM_DoubleArray_create(3);
    PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      MultipoleType multipoleType = pme.getMultipoleType(index);
      PolarizeType polarizeType = pme.getPolarizeType(index);
      int[] axisAtoms = atom.getAxisAtomIndices();

      double useFactor = 1.0;
      if (!atom.getUse() || !atom.getElectrostatics()) {
        useFactor = 0.0;
      }

      double lambdaScale = lambdaElec;
      if (!atom.applyLambda()) {
        lambdaScale = 1.0;
      }
      useFactor *= lambdaScale;

      // Define the frame definition.
      int axisType = switch (multipoleType.frameDefinition) {
        case NONE -> OpenMM_AmoebaMultipoleForce_NoAxisType;
        case ZONLY -> OpenMM_AmoebaMultipoleForce_ZOnly;
        case ZTHENX -> OpenMM_AmoebaMultipoleForce_ZThenX;
        case BISECTOR -> OpenMM_AmoebaMultipoleForce_Bisector;
        case ZTHENBISECTOR -> OpenMM_AmoebaMultipoleForce_ZBisect;
        case THREEFOLD -> OpenMM_AmoebaMultipoleForce_ThreeFold;
      };

      // Load local multipole coefficients.
      for (int j = 0; j < 3; j++) {
        OpenMM_DoubleArray_set(dipoles, j,
            multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          OpenMM_DoubleArray_set(quadrupoles, l++,
              multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * useFactor);
        }
      }

      int zaxis = 1;
      int xaxis = 1;
      int yaxis = 1;

      if (axisAtoms != null) {
        zaxis = axisAtoms[0];
        if (axisAtoms.length > 1) {
          xaxis = axisAtoms[1];
          if (axisAtoms.length > 2) {
            yaxis = axisAtoms[2];
          }
        }
      } else {
        axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
      }

      // Set the multipole parameters.
      OpenMM_AmoebaMultipoleForce_setMultipoleParameters(amoebaMultipoleForce, index,
          multipoleType.charge * useFactor, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
          polarizeType.thole, polarizeType.pdamp * dampingFactorConversion,
          polarizeType.polarizability * polarityConversion * polarScale * useFactor);
    }

    OpenMM_DoubleArray_destroy(dipoles);
    OpenMM_DoubleArray_destroy(quadrupoles);

    if (openMMContext.hasContextPointer()) {
      OpenMM_AmoebaMultipoleForce_updateParametersInContext(amoebaMultipoleForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add a Generalized Kirkwood force to the OpenMM System.
   */
  private void addGeneralizedKirkwoodForce() {
    GeneralizedKirkwood gk = openMMEnergy.getGK();

    amoebaGeneralizedKirkwoodForce = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(amoebaGeneralizedKirkwoodForce,
        gk.getSolventPermittivity());
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(amoebaGeneralizedKirkwoodForce, 1.0);
    OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(amoebaGeneralizedKirkwoodForce,
        gk.getDescreenOffset() * OpenMM_NmPerAngstrom);

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
    OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(amoebaGeneralizedKirkwoodForce,
        tanhRescale);
    OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(amoebaGeneralizedKirkwoodForce,
        betas[0], betas[1], betas[2]);

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

    for (int i = 0; i < nAtoms; i++) {
      MultipoleType multipoleType = atoms[i].getMultipoleType();
      double base = baseRadius[i] * OpenMM_NmPerAngstrom;
      double descreen = descreenRadius[i] * OpenMM_NmPerAngstrom * perfectRadiiScale;
      double overlap = overlapScale[i] * perfectRadiiScale;
      double neck = neckFactor[i] * perfectRadiiScale;
      OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(amoebaGeneralizedKirkwoodForce,
          multipoleType.charge, base, overlap, descreen, neck);

      if (!usePerfectRadii && logger.isLoggable(Level.FINE)) {
        logger.fine(format("   %s %8.6f %8.6f %5.3f", atoms[i].toString(), baseRadius[i],
            descreenRadius[i], overlapScale[i]));
      }
    }

    OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(amoebaGeneralizedKirkwoodForce,
        gk.getProbeRadius() * OpenMM_NmPerAngstrom);

    GeneralizedKirkwood.NonPolarModel nonpolar = gk.getNonPolarModel();
    switch (nonpolar) {
      default -> {
        // Configure a Born Radii based surface area term.
        double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal * OpenMM_AngstromsPerNm
            * OpenMM_AngstromsPerNm;
        OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(amoebaGeneralizedKirkwoodForce,
            OpenMM_True);
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(amoebaGeneralizedKirkwoodForce,
            -surfaceTension);
      }
      case CAV, CAV_DISP, GAUSS_DISP, SEV_DISP, HYDROPHOBIC_PMF, NONE ->
        // This NonPolar model does not use a Born Radii based surface area term.
          OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(
              amoebaGeneralizedKirkwoodForce, OpenMM_False);
    }

    int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
    OpenMM_Force_setForceGroup(amoebaGeneralizedKirkwoodForce, forceGroup);
    OpenMM_System_addForce(system, amoebaGeneralizedKirkwoodForce);

    logger.log(Level.INFO, format("  Generalized Kirkwood force \t\t%d", forceGroup));

    // Add dispersion
    switch (nonpolar) {
      case CAV_DISP, GAUSS_DISP, SEV_DISP, BORN_CAV_DISP -> addWCAForce();
      default -> {
        // WCA force is not being used.
      }
    }

    // Add cavitation
    if (nonpolar == GAUSS_DISP) {
      addCavitationForce();
    }
  }

  /**
   * Updates the AMOEBA Generalized Kirkwood force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateGeneralizedKirkwoodForce(Atom[] atoms) {
    GeneralizedKirkwood gk = openMMEnergy.getGK();

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
      OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters_1(amoebaGeneralizedKirkwoodForce,
          index, multipoleType.charge * chargeUseFactor, baseSize, overlap, descreenSize,
          neckFactor);
    }

    //        OpenMM Bug: Surface Area is not Updated by "updateParametersInContext"
    //
    //        NonPolar nonpolar = gk.getNonPolarModel();
    //        switch (nonpolar) {
    //            case BORN_SOLV:
    //            case BORN_CAV_DISP:
    //            default:
    //                // Configure a Born Radii based surface area term.
    //                double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal
    //                        * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm * lambdaElec;
    //
    // OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(amoebaGeneralizedKirkwoodForce,
    // -surfaceTension);
    //                break;
    //            case CAV:
    //            case CAV_DISP:
    //            case HYDROPHOBIC_PMF:
    //            case NONE:
    //                // This NonPolar model does not use a Born Radii based surface area term.
    //                break;
    //        }

    if (openMMContext.hasContextPointer()) {
      OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(amoebaGeneralizedKirkwoodForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add a nonpolar Weeks-Chandler-Andersen dispersion force to the OpenMM System.
   */
  private void addWCAForce() {

    GeneralizedKirkwood gk = openMMEnergy.getGK();
    DispersionRegion dispersionRegion = gk.getDispersionRegion();

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

    amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create();

    for (Atom atom : atoms) {
      VDWType vdwType = atom.getVDWType();
      double radius = vdwType.radius;
      double eps = vdwType.wellDepth;
      OpenMM_AmoebaWcaDispersionForce_addParticle(amoebaWcaDispersionForce,
          OpenMM_NmPerAngstrom * radius * radScale, OpenMM_KJPerKcal * eps);
    }

    OpenMM_AmoebaWcaDispersionForce_setEpso(amoebaWcaDispersionForce, epso * OpenMM_KJPerKcal);
    OpenMM_AmoebaWcaDispersionForce_setEpsh(amoebaWcaDispersionForce, epsh * OpenMM_KJPerKcal);
    OpenMM_AmoebaWcaDispersionForce_setRmino(amoebaWcaDispersionForce,
        rmino * OpenMM_NmPerAngstrom);
    OpenMM_AmoebaWcaDispersionForce_setRminh(amoebaWcaDispersionForce,
        rminh * OpenMM_NmPerAngstrom);
    OpenMM_AmoebaWcaDispersionForce_setDispoff(amoebaWcaDispersionForce,
        dispoff * OpenMM_NmPerAngstrom);
    OpenMM_AmoebaWcaDispersionForce_setAwater(amoebaWcaDispersionForce,
        awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
    OpenMM_AmoebaWcaDispersionForce_setSlevy(amoebaWcaDispersionForce, slevy);
    OpenMM_AmoebaWcaDispersionForce_setShctd(amoebaWcaDispersionForce, shctd);

    int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);

    OpenMM_Force_setForceGroup(amoebaWcaDispersionForce, forceGroup);
    OpenMM_System_addForce(system, amoebaWcaDispersionForce);

    logger.log(Level.INFO, format("  WCA dispersion force \t\t\t%d", forceGroup));
  }

  /**
   * Updates the WCA force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateWCAForce(Atom[] atoms) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    double radScale = 1.0;
    if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

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
      OpenMM_AmoebaWcaDispersionForce_setParticleParameters(amoebaWcaDispersionForce, index,
          OpenMM_NmPerAngstrom * radius * radScale, OpenMM_KJPerKcal * eps * useFactor);
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(amoebaWcaDispersionForce, openMMContext.getContextPointer());
    }
  }

  /**
   * Add a GaussVol cavitation force to the OpenMM System.
   */
  private void addCavitationForce() {

    GeneralizedKirkwood generalizedKirkwood = openMMEnergy.getGK();
    ChandlerCavitation chandlerCavitation = generalizedKirkwood.getChandlerCavitation();
    GaussVol gaussVol = chandlerCavitation.getGaussVol();
    if (gaussVol == null) {
      return;
    }

    amoebaCavitationForce = OpenMM_AmoebaGKCavitationForce_create();
    double surfaceTension =
        chandlerCavitation.getSurfaceTension() * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom
            / OpenMM_NmPerAngstrom;
    double[] rad = gaussVol.getRadii();

    int index = 0;
    for (Atom atom : atoms) {
      int isHydrogen = OpenMM_False;
      double radius = rad[index++];
      if (atom.isHydrogen()) {
        isHydrogen = OpenMM_True;
        radius = 0.0;
      }
      OpenMM_AmoebaGKCavitationForce_addParticle(amoebaCavitationForce,
          radius * OpenMM_NmPerAngstrom, surfaceTension, isHydrogen);
    }

    OpenMM_AmoebaGKCavitationForce_setNonbondedMethod(amoebaCavitationForce,
        OpenMM_AmoebaGKCavitationForce_NoCutoff);

    int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 2);
    OpenMM_Force_setForceGroup(amoebaCavitationForce, forceGroup);
    OpenMM_System_addForce(system, amoebaCavitationForce);

    logger.log(Level.INFO, format("  GaussVol cavitation force \t\t%d", forceGroup));
  }

  /**
   * Updates the Cavitation force for changes in Use flags or Lambda.
   *
   * @param atoms Array of atoms to update.
   */
  private void updateCavitationForce(Atom[] atoms) {
    GeneralizedKirkwood generalizedKirkwood = openMMEnergy.getGK();
    ChandlerCavitation chandlerCavitation = generalizedKirkwood.getChandlerCavitation();
    GaussVol gaussVol = chandlerCavitation.getGaussVol();
    if (gaussVol == null) {
      return;
    }

    double surfaceTension =
        chandlerCavitation.getSurfaceTension() * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom
            / OpenMM_NmPerAngstrom;

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

      OpenMM_AmoebaGKCavitationForce_setParticleParameters(amoebaCavitationForce, index,
          radius * OpenMM_NmPerAngstrom, surfaceTension * useFactor, isHydrogen);
    }

    if (openMMContext.hasContextPointer()) {
      OpenMM_AmoebaGKCavitationForce_updateParametersInContext(amoebaCavitationForce, openMMContext.getContextPointer());
    }
  }

  public boolean hasAmoebaCavitationForce() {
    return amoebaCavitationForce != null;
  }

  /**
   * Adds harmonic restraints (CoordRestraint objects) to OpenMM as a custom external force.
   */
  private void addRestrainPositionForce() {
    List<RestrainPosition> restrainPositionList = openMMEnergy.getRestrainPositions();
    if (restrainPositionList == null || restrainPositionList.isEmpty()) {
      return;
    }

    int nRestraints = restrainPositionList.size();
    int forceGroup = forceField.getInteger("RESTRAIN_POSITION_FORCE_GROUP", 0);

    for (RestrainPosition restrainPosition : openMMEnergy.getRestrainPositions()) {
      double forceConstant = restrainPosition.getForceConstant();
      forceConstant *= OpenMM_KJPerKcal;
      forceConstant *= (OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm);
      Atom[] restAtoms = restrainPosition.getAtoms();
      int nRestAts = restrainPosition.getNumAtoms();
      double[][] oCoords = restrainPosition.getEquilibriumCoordinates();
      for (int i = 0; i < nRestAts; i++) {
        oCoords[i][0] *= OpenMM_NmPerAngstrom;
        oCoords[i][1] *= OpenMM_NmPerAngstrom;
        oCoords[i][2] *= OpenMM_NmPerAngstrom;
      }

      PointerByReference theRestraint = OpenMM_CustomExternalForce_create("k*periodicdistance(x,y,z,x0,y0,z0)^2");
      OpenMM_CustomExternalForce_addGlobalParameter(theRestraint, "k", forceConstant);
      OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "x0");
      OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "y0");
      OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "z0");

      PointerByReference xyzOrigArray = OpenMM_DoubleArray_create(3);
      for (int i = 0; i < nRestAts; i++) {
        int ommIndex = restAtoms[i].getXyzIndex() - 1;
        for (int j = 0; j < 3; j++) {
          OpenMM_DoubleArray_set(xyzOrigArray, j, oCoords[i][j]);
        }
        OpenMM_CustomExternalForce_addParticle(theRestraint, ommIndex, xyzOrigArray);
      }
      OpenMM_DoubleArray_destroy(xyzOrigArray);
      OpenMM_Force_setForceGroup(theRestraint, forceGroup);
      OpenMM_System_addForce(system, theRestraint);
    }

    logger.log(Level.INFO, format("  Restrain Positions\t%6d\t%d", nRestraints, forceGroup));
  }

  /**
   * Adds restraint bonds, if any.
   */
  private void addRestraintBondForce() {
    List<RestraintBond> restraintBonds = openMMEnergy.getRestraintBonds();
    if (restraintBonds == null || restraintBonds.isEmpty()) {
      return;
    }

    int forceGroup = forceField.getInteger("BOND_RESTRAINT_FORCE_GROUP", 0);

    // OpenMM's HarmonicBondForce class uses k, not 1/2*k as does FFX.
    double kParameterConversion =
        2.0 * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    // Map from bond functional forms to the restraint-bonds using that functional form.
    Map<BondType.BondFunction, PointerByReference> restraintForces = new HashMap<>();

    for (RestraintBond restraintBond : openMMEnergy.getRestraintBonds()) {
      PointerByReference theForce;
      BondType bondType = restraintBond.bondType;
      BondType.BondFunction bondFunction = bondType.bondFunction;
      if (restraintForces.containsKey(bondFunction)) {
        theForce = restraintForces.get(bondFunction);
      } else {
        theForce = OpenMM_CustomBondForce_create(bondFunction.toMathematicalForm());
        OpenMM_CustomBondForce_addPerBondParameter(theForce, "k");
        OpenMM_CustomBondForce_addPerBondParameter(theForce, "r0");
        if (bondFunction.hasFlatBottom()) {
          OpenMM_CustomBondForce_addPerBondParameter(theForce, "fb");
        }

        switch (bondFunction) {
          case QUARTIC, FLAT_BOTTOM_QUARTIC -> {
            OpenMM_CustomBondForce_addGlobalParameter(theForce, "cubic",
                bondType.cubic / OpenMM_NmPerAngstrom);
            OpenMM_CustomBondForce_addGlobalParameter(theForce, "quartic",
                bondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
          }
          default -> {
            // Do nothing.
          }
        }

        OpenMM_Force_setForceGroup(theForce, forceGroup);
        OpenMM_System_addForce(system, theForce);
        restraintForces.put(bondFunction, theForce);
      }

      double forceConst = bondType.forceConstant * bondType.bondUnit * kParameterConversion;
      double equilDist = bondType.distance * OpenMM_NmPerAngstrom;
      Atom[] ats = restraintBond.getAtomArray();
      int at1 = ats[0].getXyzIndex() - 1;
      int at2 = ats[1].getXyzIndex() - 1;

      PointerByReference bondParams = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(bondParams, forceConst);
      OpenMM_DoubleArray_append(bondParams, equilDist);
      if (bondFunction.hasFlatBottom()) {
        OpenMM_DoubleArray_append(bondParams, bondType.flatBottomRadius * OpenMM_NmPerAngstrom);
      }
      OpenMM_CustomBondForce_addBond(theForce, at1, at2, bondParams);
      OpenMM_DoubleArray_destroy(bondParams);
    }

    logger.log(Level.INFO,
        format("  Restraint bonds force \t%6d\t%d", restraintBonds.size(), forceGroup));
  }

  private void addRestrainGroupsForce() {
    RestrainGroups restrainGroups = openMMEnergy.getRestrainGroups();
    if (restrainGroups == null) {
      return;
    }

    // In the expression below, u and l are the upper and lower threshold
    PointerByReference force = OpenMM_CustomCentroidBondForce_create(2,
        "step(distance(g1,g2)-u)*k*(distance(g1,g2)-u)^2+step(l-distance(g1,g2))*k*(distance(g1,g2)-l)^2");
    OpenMM_CustomCentroidBondForce_addPerBondParameter(force, "k");
    OpenMM_CustomCentroidBondForce_addPerBondParameter(force, "l");
    OpenMM_CustomCentroidBondForce_addPerBondParameter(force, "u");

    // Create the Restrain Groups.
    int nGroups = restrainGroups.getNumberOfGroups();
    for (int j = 0; j < nGroups; j++) {
      PointerByReference group = OpenMM_IntArray_create(0);
      PointerByReference weight = OpenMM_DoubleArray_create(0);
      int[] groupMembers = restrainGroups.getGroupMembers(j);
      for (int i : groupMembers) {
        OpenMM_IntArray_append(group, i);
        OpenMM_DoubleArray_append(weight, atoms[i].getMass());
      }
      OpenMM_CustomCentroidBondForce_addGroup(force, group, weight);
      OpenMM_IntArray_destroy(group);
      OpenMM_DoubleArray_destroy(weight);
    }

    // Add the restraints between groups.
    double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    int nRestraints = restrainGroups.getNumberOfRestraints();
    int[] group1 = restrainGroups.getGroup1();
    int[] group2 = restrainGroups.getGroup2();
    double[] forceConstants = restrainGroups.getForceConstants();
    double[] smallerDistance = restrainGroups.getSmallerDistance();
    double[] largerDistance = restrainGroups.getLargerDistance();
    for (int i = 0; i < nRestraints; i++) {
      PointerByReference group = OpenMM_IntArray_create(0);
      OpenMM_IntArray_append(group, group1[i]);
      OpenMM_IntArray_append(group, group2[i]);
      PointerByReference params = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(params, forceConstants[i] * convert);
      OpenMM_DoubleArray_append(params, smallerDistance[i] * OpenMM_NmPerAngstrom);
      OpenMM_DoubleArray_append(params, largerDistance[i] * OpenMM_NmPerAngstrom);
      OpenMM_CustomCentroidBondForce_addBond(force, group, params);
      OpenMM_IntArray_destroy(group);
      OpenMM_DoubleArray_destroy(params);
    }

    // Add the constraint force.
    int forceGroup = forceField.getInteger("RESTRAIN_GROUPS_FORCE_GROUP", 0);
    OpenMM_Force_setForceGroup(force, forceGroup);

    if (openMMEnergy.getCrystal().aperiodic()) {
      OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(force, OpenMM_False);
    } else {
      OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(force, OpenMM_True);
    }
    OpenMM_System_addForce(system, force);
    logger.log(Level.INFO, format("  Restrain Groups \t%6d\t\t%1d", nRestraints, forceGroup));
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
