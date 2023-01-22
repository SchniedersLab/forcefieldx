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
package ffx.potential.nonbonded;

import static ffx.numerics.atomic.AtomicDoubleArray.atomicDoubleArrayFactory;
import static ffx.potential.nonbonded.implicit.DispersionRegion.DEFAULT_DISPERSION_OFFSET;
import static ffx.potential.parameters.ForceField.toEnumForm;
import static ffx.potential.parameters.SoluteType.setSoluteRadii;
import static ffx.utilities.Constants.dWater;
import static ffx.utilities.KeywordGroup.ImplicitSolvent;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.max;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.nonbonded.implicit.BornGradRegion;
import ffx.potential.nonbonded.implicit.BornRadiiRegion;
import ffx.potential.nonbonded.implicit.BornTanhRescaling;
import ffx.potential.nonbonded.implicit.ChandlerCavitation;
import ffx.potential.nonbonded.implicit.ConnollyRegion;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.nonbonded.implicit.GKEnergyRegion;
import ffx.potential.nonbonded.implicit.GaussVol;
import ffx.potential.nonbonded.implicit.InducedGKFieldRegion;
import ffx.potential.nonbonded.implicit.InitializationRegion;
import ffx.potential.nonbonded.implicit.PermanentGKFieldRegion;
import ffx.potential.nonbonded.implicit.SurfaceAreaRegion;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.SoluteType;
import ffx.potential.parameters.SoluteType.SOLUTE_RADII_TYPE;
import ffx.utilities.FFXKeyword;

import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This Generalized Kirkwood class implements GK for the AMOEBA polarizable atomic multipole force
 * field in parallel using a {@link ffx.potential.nonbonded.NeighborList}.
 *
 * @author Michael J. Schnieders<br> derived from:<br> TINKER code by Michael J. Schnieders and Jay
 *     W. Ponder<br>
 * @see <a href="http://dx.doi.org/10.1021/ct7001336" target="_blank">M. J. Schnieders and J. W.
 *     Ponder, Polarizable atomic multipole solutes in a generalized Kirkwood continuum, Journal of
 *     Chemical Theory and Computation 2007, 3, (6), 2083-2097.</a><br>
 * @since 1.0
 */
public class GeneralizedKirkwood implements LambdaInterface {

  /**
   * Default dielectric offset
   */
  public static final double DEFAULT_DIELECTRIC_OFFSET = 0.09;

  private static final Logger logger = Logger.getLogger(GeneralizedKirkwood.class.getName());
  /**
   * Default Bondi scale factor.
   */
  private static final double DEFAULT_SOLUTE_SCALE = 1.0;

  /**
   * Dielectric offset from:
   *
   * <p>W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson, "A Semianalytical Treatment of
   * Solvation for Molecular Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129 (1990)
   */
  private final double dOffset = DEFAULT_DIELECTRIC_OFFSET;
  /**
   * Force field in use.
   */
  private final ForceField forceField;
  /**
   * Treatment of polarization.
   */
  private final Polarization polarization;
  /**
   * Particle mesh Ewald instance, which contains variables such as expanded coordinates and
   * multipoles in the global frame that GK uses.
   */
  private final ParticleMeshEwald particleMeshEwald;
  /**
   * Parallel team object for shared memory parallelization.
   */
  private final ParallelTeam parallelTeam;
  /**
   * Initialize GK output variables.
   */
  private final InitializationRegion initializationRegion;
  /**
   * Parallel computation of Born Radii.
   */
  private final BornRadiiRegion bornRadiiRegion;
  /**
   * Parallel computation of the Permanent GK Field.
   */
  private final PermanentGKFieldRegion permanentGKFieldRegion;
  /**
   * Parallel computation of the Induced GK Field.
   */
  private final InducedGKFieldRegion inducedGKFieldRegion;
  /**
   * Parallel computation of the GK continuum electrostatics energy.
   */
  private final GKEnergyRegion gkEnergyRegion;
  /**
   * Parallel computation of Born radii chain rule term.
   */
  private final BornGradRegion bornGradRegion;
  /**
   * Parallel computation of Dispersion energy.
   */
  private final DispersionRegion dispersionRegion;
  /**
   * Parallel computation of Cavitation.
   */
  private final SurfaceAreaRegion surfaceAreaRegion;
  /**
   * Volume to Surface Area Cavitation Dependence.
   */
  private final ChandlerCavitation chandlerCavitation;
  /**
   * Maps radii overrides (by AtomType) specified from the command line. e.g.
   * -DradiiOverride=134r1.20,135r1.20 sets atom types 134,135 to Bondi=1.20
   */
  private final HashMap<Integer, Double> radiiOverride = new HashMap<>();
  /**
   * Conversion from electron**2/Ang to kcal/mole.
   */
  public double electric;
  /**
   * Empirical scaling of the Bondi radii.
   */
  private final double bondiScale;

  /**
   * Array of Atoms being considered.
   */
  private Atom[] atoms;
  /**
   * Number of atoms.
   */
  private int nAtoms;
  /**
   * Cartesian coordinates of each atom.
   */
  private double[] x, y, z;
  /**
   * Base radius of each atom.
   */
  private double[] baseRadius;
  /**
   * Descreen radius of each atom.
   */
  private double[] descreenRadius;
  /**
   * Overlap scale factor for each atom, when using the Hawkins, Cramer & Truhlar pairwise
   * descreening algorithm.
   *
   * <p>G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized Models of Aqueous Free Energies
   * of Solvation Based on Pairwise Descreening of Solute Atomic Charges from a Dielectric Medium",
   * J. Phys. Chem., 100, 19824-19839 (1996).
   */
  private double[] overlapScale;
  /**
   * Sneck scaling parameter for each atom. Set based on maximum Sneck scaling parameter and number
   * of bound non-hydrogen atoms
   */
  private double[] neckScale;
  /**
   * If true, the descreening integral includes overlaps with the volume of the descreened atom
   */
  private final boolean perfectHCTScale;

  /**
   * The requested permittivity for the solvent.
   */
  @FFXKeyword(name = "solvent-dielectric", keywordGroup = ImplicitSolvent, defaultValue = "78.3",
      description = "The dielectric constant used for the solvent in generalized Kirkwood calculations."
          + "The default of 78.3 corresponds to water.")
  private final double solventDielectric;
  /**
   * The requested permittivity for the solute.
   */
  @FFXKeyword(name = "solute-dielectric", keywordGroup = ImplicitSolvent, defaultValue = "1.0",
      description = "The dielectric constant used for the solute(s) in generalized Kirkwood calculations."
          + "The default of 1.0 is consistent with all solute dielectric response arising from either"
          + "polarization via induced dipoles and/or permanent dipole realignment.")
  private final double soluteDielectric;

  /**
   * Default overlap scale factor for the Hawkins, Cramer & Truhlar pairwise descreening algorithm:
   * 0.69 New default overlap scale factor set during implicit solvent model optimization: 0.72
   */
  private static final double DEFAULT_HCT_SCALE = 0.72;

  /**
   * Base overlap HCT overlap scale factor.
   */
  @FFXKeyword(name = "hct-scale", keywordGroup = ImplicitSolvent, defaultValue = "0.72",
      description = "The default overlap scale factor for Hawkins-Cramer-Truhlar pairwise descreening.")
  private final double hctScale;

  /**
   * If true, HCT overlap scale factors are element-specific
   */
  @FFXKeyword(name = "element-hct-scale", clazz = Boolean.class, keywordGroup = ImplicitSolvent, defaultValue = "false",
      description = "Flag to turn on element specific overlap scale factors for Hawkins-Cramer-Truhlar pairwise descreening.")
  private final boolean elementHCTScale;

  /**
   * Element-specific HCT overlap scale factors
   */
  private final HashMap<Integer, Double> elementHCTScaleFactors;

  /**
   * If true, the descreening size of atoms is based on their force field vdW radius
   */
  @FFXKeyword(name = "descreen-vdw", keywordGroup = ImplicitSolvent, defaultValue = "true",
      description = "If true, the descreening size of each atom is based on its force field van der Waals radius.")
  private final boolean descreenVDW;

  /**
   * If true, hydrogen atoms displace solvent during the pairwise descreening integral.
   */
  @FFXKeyword(name = "descreen-hydrogen", keywordGroup = ImplicitSolvent, defaultValue = "false",
      description = "If true, hydrogen atoms are contribute to the pairwise descreening integrals.")
  private final boolean descreenHydrogen;

  private static final double DEFAULT_DESCREEN_OFFSET = 0.0;
  /**
   * Offset applied to the pairwise descreening integral to improve stability at small separation.
   */
  @FFXKeyword(name = "descreen-offset", keywordGroup = ImplicitSolvent, defaultValue = "0.0",
      description = "Offset applied to the pairwise descreening integral to improve stability at small separation.")
  private final double descreenOffset;

  /**
   * Apply a neck correction during descreening.
   */
  @FFXKeyword(name = "neck-correction", keywordGroup = ImplicitSolvent, defaultValue = "false",
      description = "Apply a neck correction during descreening.")
  private final boolean neckCorrection;

  /**
   * Default Sneck scaling factor from Aguilar/Onufriev 2010
   */
  private static final double DEFAULT_NECK_SCALE = 0.6784;
  /**
   * Maximum Sneck scaling parameter value
   */
  @FFXKeyword(name = "neck-scale", keywordGroup = ImplicitSolvent, defaultValue = "",
      description = "The overlap scale factor to use during the descreening neck correction.")
  private double sneck;

  /**
   * Use the Corrigan et al. chemically aware neck correction; atoms with more heavy atom bonds are
   * less capable of forming interstitial necks.
   */
  @FFXKeyword(name = "chemically-aware-neck-scale", keywordGroup = ImplicitSolvent, defaultValue = "true",
      description = "If the neck descreening correction is being used, apply a smaller overlap scale"
          + "factors as the number of bonded heavy atoms increases.")
  private final boolean chemicallyAwareSneck;

  /**
   * If true, the descreening integral includes the tanh correction to better approximate molecular
   * surface
   */
  @FFXKeyword(name = "tanh-correction", keywordGroup = ImplicitSolvent, defaultValue = "false",
      description = "If the neck descreening correction is being used, apply a smaller overlap scale"
          + "factors as the number of bonded heavy atoms increases.")
  private final boolean tanhCorrection;

  /**
   * Default value of beta0 for tanh scaling
   */
  public static final double DEFAULT_TANH_BETA0 = 0.9563;
  /**
   * Default value of beta1 for tanh scaling
   */
  public static final double DEFAULT_TANH_BETA1 = 0.2578;
  /**
   * Default value of beta2 for tanh scaling
   */
  public static final double DEFAULT_TANH_BETA2 = 0.0810;

  /**
   * The coefficient beta0 for tanh rescaling of descreening integrals.
   */
  @FFXKeyword(name = "tanh-beta0", keywordGroup = ImplicitSolvent, defaultValue = "0.9563",
      description = "The coefficient beta0 for tanh rescaling of descreening integrals.")
  private double beta0;

  /**
   * The coefficient beta1 for tanh rescaling of descreening integrals.
   */
  @FFXKeyword(name = "tanh-beta1", keywordGroup = ImplicitSolvent, defaultValue = "0.2578",
      description = "The coefficient beta1 for tanh rescaling of descreening integrals.")
  private double beta1;

  /**
   * The coefficient beta2 for tanh rescaling of descreening integrals.
   */
  @FFXKeyword(name = "tanh-beta2", keywordGroup = ImplicitSolvent, defaultValue = "0.0810",
      description = "The coefficient beta2 for tanh rescaling of descreening integrals.")
  private double beta2;

  /**
   * Default constant for the Generalized Kirkwood cross-term.
   */
  public static final double DEFAULT_GKC = 2.455;
  /**
   * The Generalized Kirkwood cross-term parameter.
   */
  @FFXKeyword(name = "gkc", keywordGroup = ImplicitSolvent, defaultValue = "2.455",
      description = "The Generalized Kirkwood cross-term parameter.")
  public final double gkc;

  private static final NonPolarModel DEFAULT_NONPOLAR_MODEL = NonPolarModel.CAV_DISP;

  /**
   * Treatment of non-polar interactions.
   */
  @FFXKeyword(name = "nonpolar-model", clazz = String.class,
      keywordGroup = ImplicitSolvent, defaultValue = "cav-disp",
      description = "[CAV / CAV-DISP / GAUSS-DISP / SEV-DISP / NONE ] "
          + "The non-polar contribution to the implicit solvent.")
  private final NonPolarModel nonPolarModel;

  /**
   * Default surface tension for non-polar models without an explicit dispersion term. This is lower
   * than CAVDISP, since the favorable dispersion term is implicitly included.
   */
  private static final double DEFAULT_CAVONLY_SURFACE_TENSION = 0.0049;

  /**
   * Water probe radius.
   */
  public final double probe;

  /**
   * Default solvent pressure for apolar models with an explicit volume term.
   *
   * <p>From work by Chandler et al., the following relationship for cavitation free energy as a
   * function of spherical particle size was found: Cross-Over = 3 * S.T. / S.P.
   *
   * <p>A S.P. of 0.0334 kcal/mol/A^3 was obtained using explicit AMOEBA water simulations and
   * solvent excluded volumes.
   *
   * <p>A S.P. of 0.0343 kcal/mol/A^3 is obtained assuming a macroscopic surface tension of 0.103
   * kcal/mol/A^3 and a cross-over of 9.0 (i.e. S.P. = 3 * S.T. / Cross-Over)
   *
   * <p>Both values are in reasonably good agreement, and 0.0334 is chosen as our default.
   */
  public static final double DEFAULT_SOLVENT_PRESSURE = 0.0334;

  /**
   * Default surface tension for apolar models with an explicit dispersion term.
   *
   * <p>Experimental value: 0.103 kcal/mol/Ang^2
   */
  public static final double DEFAULT_CAVDISP_SURFACE_TENSION = 0.103;

  /**
   * Using a S.P. of 0.0334 kcal/mol/A^3, and a limiting surface tension of 0.103 kcal/mol/A^2, the
   * cross-over point is 9.2515 A.
   *
   * <p>Using a S.P. of 0.0334 kcal/mol/A^3, and a limiting surface tension of 0.08 (i.e. 80% of the
   * experimentally observed surface tension of 0.103 kcal/mol/A^2), we derive a cross-over of:
   *
   * <p>9.251 A = 3 * 0.103 kcal/mol/A^2 / 0.0334 kcal/mol/A^3.
   */
  public static final double DEFAULT_CROSSOVER =
      3.0 * DEFAULT_CAVDISP_SURFACE_TENSION / DEFAULT_SOLVENT_PRESSURE;

  /**
   * Cavitation surface tension coefficient (kcal/mol/A^2).
   */
  @FFXKeyword(name = "surface-tension", keywordGroup = ImplicitSolvent, defaultValue = "0.103",
      description = "The cavitation surface tension coefficient (kcal/mol/A^2).")
  private final double surfaceTension;

  /**
   * Cavitation solvent pressure coefficient (kcal/mol/A^3).
   */
  @FFXKeyword(name = "solvent-pressure", keywordGroup = ImplicitSolvent, defaultValue = "0.0334",
      description = "The solvent pressure for nonpolar models with an explicit volume term (kcal/mol/A^3).")
  private final double solventPressue;

  /**
   * The base radii to use for GK.
   */
  @FFXKeyword(name = "gk-radius", clazz = String.class, keywordGroup = ImplicitSolvent, defaultValue = "solute",
      description = "[SOLUTE / VDW / CONSENSUS] "
          + "The base atomic radii to use for generalized Kirkwood calculations. The default is to use solute radii, "
          + "which were fit to experimental solvation free energy differences. Alternatively, force field"
          + "specific van der Waals radii (vdw) or consensus Bondi radii (consensus) can be chosen.")
  private SOLUTE_RADII_TYPE soluteRadiiType;

  /**
   * Born radius of each atom.
   */
  private double[] born;
  /**
   * Flag to indicate if an atom should be included.
   */
  private boolean[] use = null;
  /**
   * Periodic boundary conditions and symmetry.
   */
  private Crystal crystal;
  /**
   * Atomic coordinates for each symmetry operator.
   */
  private double[][][] sXYZ;
  /**
   * Multipole moments for each symmetry operator.
   */
  private double[][][] globalMultipole;
  /**
   * Induced dipoles for each symmetry operator.
   */
  private double[][][] inducedDipole;
  /**
   * Induced dipole chain rule terms for each symmetry operator.
   */
  private double[][][] inducedDipoleCR;
  /**
   * AtomicDoubleArray implementation to use.
   */
  private AtomicDoubleArrayImpl atomicDoubleArrayImpl;
  /**
   * Atomic Gradient array.
   */
  private AtomicDoubleArray3D grad;
  /**
   * Atomic Torque array.
   */
  private AtomicDoubleArray3D torque;
  /**
   * Atomic Born radii chain-rule array.
   */
  private AtomicDoubleArray bornRadiiChainRule;
  /**
   * Atomic GK field array.
   */
  private final AtomicDoubleArray3D fieldGK;
  /**
   * Atomic GK field chain-rule array.
   */
  private final AtomicDoubleArray3D fieldGKCR;
  /**
   * Neighbor lists for each atom and symmetry operator.
   */
  private int[][][] neighborLists;
  /**
   * This field is because re-initializing the force field resizes some arrays but not others; that
   * second category must, when called on, be resized not to the current number of atoms but to the
   * maximum number of atoms (and thus to the size of the first category of arrays).
   */
  private int maxNumAtoms;
  /**
   * GK cut-off distance.
   */
  private double cutoff;
  /**
   * GK cut-off distance squared.
   */
  private double cut2;
  /**
   * Boolean flag to indicate GK will be scaled by the lambda state variable.
   */
  private boolean lambdaTerm;
  /**
   * The current value of the lambda state variable.
   */
  private double lambda = 1.0;
  /**
   * lPow equals lambda^polarizationLambdaExponent, where polarizationLambdaExponent is also used by
   * PME.
   */
  private double lPow = 1.0;
  /**
   * First derivative of lPow with respect to l.
   */
  private double dlPow = 0.0;
  /**
   * Second derivative of lPow with respect to l.
   */
  private double dl2Pow = 0.0;
  /**
   * Total Solvation Energy.
   */
  private double solvationEnergy = 0.0;
  /**
   * Electrostatic Solvation Energy.
   */
  private double gkEnergy = 0.0;
  /**
   * Electrostatic Solvation Energy from Permanent Multipoles.
   */
  private double gkPermanentEnergy = 0.0;
  /**
   * Electrostatic Solvation Energy from Induced Dipoles.
   */
  private double gkPolarizationEnergy = 0.0;
  /**
   * Dispersion Solvation Energy.
   */
  private double dispersionEnergy = 0.0;
  /**
   * Cavitation Solvation Energy.
   */
  private double cavitationEnergy = 0.0;
  /**
   * Time to compute GK electrostatics.
   */
  private long gkTime = 0;
  /**
   * Time to compute Dispersion energy.
   */
  private long dispersionTime = 0;
  /**
   * Time to compute Cavitation energy.
   */
  private long cavitationTime = 0;
  /**
   * Forces all atoms to be considered during Born radius updates.
   */
  private final boolean nativeEnvironmentApproximation;
  /**
   * Flag to turn on use of perfect Born radii.
   */
  private final boolean perfectRadii;

  /**
   * Constructor for GeneralizedKirkwood.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param particleMeshEwald a {@link ParticleMeshEwald} object.
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
   */
  public GeneralizedKirkwood(ForceField forceField, Atom[] atoms,
      ParticleMeshEwald particleMeshEwald, Crystal crystal, ParallelTeam parallelTeam,
      double electric, double gkCutoff) {
    this.forceField = forceField;
    this.atoms = atoms;
    this.particleMeshEwald = particleMeshEwald;
    this.crystal = crystal;
    this.parallelTeam = parallelTeam;
    this.electric = electric;
    this.cutoff = gkCutoff;
    nAtoms = atoms.length;
    maxNumAtoms = nAtoms;
    polarization = particleMeshEwald.polarization;

    // Set the Kirkwood multipolar reaction field constants for solvent.
    solventDielectric = forceField.getDouble("SOLVENT_DIELECTRIC", dWater);
    // Set the Kirkwood multipolar reaction field constants for solute.
    soluteDielectric = forceField.getDouble("SOLUTE_DIELECTRIC", 1.0);

    // Define how force arrays will be accumulated.
    String value = forceField.getString("ARRAY_REDUCTION", "MULTI");
    try {
      atomicDoubleArrayImpl = AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
    } catch (Exception e) {
      atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
      logger.info(
          format(" Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value,
              atomicDoubleArrayImpl));
    }

    String gkRadius = forceField.getString("GK_RADIUS", "SOLUTE");
    try {
      soluteRadiiType = SOLUTE_RADII_TYPE.valueOf(gkRadius.trim().toUpperCase());
    } catch (Exception e) {
      soluteRadiiType = SOLUTE_RADII_TYPE.SOLUTE;
    }

    // Define default Bondi scale factor, and HCT overlap scale factors.
    if (soluteRadiiType != SOLUTE_RADII_TYPE.SOLUTE) {
      bondiScale = forceField.getDouble("SOLUTE_SCALE", DEFAULT_SOLUTE_SCALE);
    } else {
      // Default Bondi scale factor for Solute Radii is 1.0.
      bondiScale = forceField.getDouble("SOLUTE_SCALE", 1.0);
    }

    hctScale = forceField.getDouble("HCT_SCALE", DEFAULT_HCT_SCALE);
    elementHCTScale = forceField.getBoolean("ELEMENT_HCT_SCALE", false);
    descreenVDW = forceField.getBoolean("DESCREEN_VDW", true);
    descreenHydrogen = forceField.getBoolean("DESCREEN_HYDROGEN", false);
    if (descreenVDW && !descreenHydrogen) {
      perfectHCTScale = forceField.getBoolean("PERFECT_HCT_SCALE", false);
    } else {
      perfectHCTScale = false;
    }
    descreenOffset = forceField.getDouble("DESCREEN_OFFSET", DEFAULT_DESCREEN_OFFSET);
    // If true, the descreening integral includes the neck correction to better approximate molecular surface.
    neckCorrection = forceField.getBoolean("NECK_CORRECTION", false);
    double sn = forceField.getDouble("SNECK", DEFAULT_NECK_SCALE);
    sneck = forceField.getDouble("NECK_SCALE", sn);
    chemicallyAwareSneck = forceField.getBoolean("CHEMICALLY_AWARE_SNECK", true);
    tanhCorrection = forceField.getBoolean("TANH_CORRECTION", false);
    double b0 = forceField.getDouble("BETA0", DEFAULT_TANH_BETA0);
    double b1 = forceField.getDouble("BETA1", DEFAULT_TANH_BETA1);
    double b2 = forceField.getDouble("BETA2", DEFAULT_TANH_BETA2);
    beta0 = forceField.getDouble("TANH_BETA0", b0);
    beta1 = forceField.getDouble("TANH_BETA1", b1);
    beta2 = forceField.getDouble("TANH_BETA2", b2);
    
    BornTanhRescaling.setBeta0(beta0);
    BornTanhRescaling.setBeta1(beta1);
    BornTanhRescaling.setBeta2(beta2);

    // Default overlap element specific scale factors for the Hawkins, Cramer & Truhlar pairwise descreening algorithm.
    HashMap<Integer, Double> DEFAULT_HCT_ELEMENTS = new HashMap<>();
    // Fit default values from Corrigan et. al. interstitial spaces work
    DEFAULT_HCT_ELEMENTS.put(1,0.7200);
    DEFAULT_HCT_ELEMENTS.put(6,0.6950);
    DEFAULT_HCT_ELEMENTS.put(7,0.7673);
    DEFAULT_HCT_ELEMENTS.put(8,0.7965);
    DEFAULT_HCT_ELEMENTS.put(15,0.6117);
    DEFAULT_HCT_ELEMENTS.put(16,0.7204);

    // Add default values for all elements
    elementHCTScaleFactors = new HashMap<>();
    elementHCTScaleFactors.put(1, DEFAULT_HCT_ELEMENTS.get(1));
    elementHCTScaleFactors.put(6, DEFAULT_HCT_ELEMENTS.get(6));
    elementHCTScaleFactors.put(7, DEFAULT_HCT_ELEMENTS.get(7));
    elementHCTScaleFactors.put(8, DEFAULT_HCT_ELEMENTS.get(8));
    elementHCTScaleFactors.put(15, DEFAULT_HCT_ELEMENTS.get(15));
    elementHCTScaleFactors.put(16, DEFAULT_HCT_ELEMENTS.get(16));

    // Fill elementHCTScaleFactors array based on input hct-element property flags
    // Overwrite defaults with input values, add new values if input elements aren't represented
    // Input lines read in as "hct-element atomicNumber hctValue"
    // So "hct-element 1 0.7200" would set the HCT scaling factor value for hydrogen (atomic number 1) to 0.7200
    String[] tmpElemScaleFactors = forceField.getProperties().getStringArray("hct-element");
    for(String elemScale : tmpElemScaleFactors){
      String[] singleScale = elemScale.trim().split(" +");
      if(elementHCTScaleFactors.containsKey(Integer.parseInt(singleScale[0]))) {
        elementHCTScaleFactors.replace(Integer.parseInt(singleScale[0]), Double.parseDouble(singleScale[1]));
      } else {
        elementHCTScaleFactors.put(Integer.parseInt(singleScale[0]), Double.parseDouble(singleScale[1]));
      }
    }

    perfectRadii = forceField.getBoolean("PERFECT_RADII", false);

    // Process any radii override values.
    String radiiProp = forceField.getString("GK_RADIIOVERRIDE", null);
    if (radiiProp != null) {
      String[] tokens = radiiProp.split("A");
      for (String token : tokens) {
        if (!token.contains("r")) {
          logger.severe("Invalid radius override.");
        }
        int separator = token.indexOf("r");
        int type = Integer.parseInt(token.substring(0, separator));
        double factor = Double.parseDouble(token.substring(separator + 1));
        logger.info(format(" (GK) Scaling AtomType %d with bondi factor %.2f", type, factor));
        radiiOverride.put(type, factor);
      }
    }

    NonPolarModel npModel;
    try {
      String cavModel = forceField.getString("CAVMODEL", DEFAULT_NONPOLAR_MODEL.name());
      cavModel = forceField.getString("NONPOLAR_MODEL", cavModel);
      npModel = getNonPolarModel(cavModel.toUpperCase());
    } catch (Exception ex) {
      npModel = NonPolarModel.NONE;
      logger.warning(format(" Error parsing non-polar model (set to NONE) %s", ex));
    }
    nonPolarModel = npModel;

    int threadCount = parallelTeam.getThreadCount();
    fieldGK = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);
    fieldGKCR = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);

    nativeEnvironmentApproximation =
        forceField.getBoolean("NATIVE_ENVIRONMENT_APPROXIMATION", false);
    probe = forceField.getDouble("PROBE_RADIUS", 1.4);
    cut2 = cutoff * cutoff;
    lambdaTerm = forceField.getBoolean("GK_LAMBDATERM", forceField.getBoolean("LAMBDATERM", false));

    // If polarization lambda exponent is set to 0.0, then we're running
    // Dual-Topology and the GK energy will be scaled with the overall system lambda value.
    double polLambdaExp = forceField.getDouble("POLARIZATION_LAMBDA_EXPONENT", 3.0);
    if (polLambdaExp == 0.0) {
      lambdaTerm = false;
      logger.info(" GK lambda term set to false.");
    }

    // If PME includes polarization and is a function of lambda, GK must also.
    if (!lambdaTerm
        && particleMeshEwald.getPolarizationType() != Polarization.NONE) {
      if (forceField.getBoolean("ELEC_LAMBDATERM",
          forceField.getBoolean("LAMBDATERM", false))) {
        logger.info(" If PME includes polarization and is a function of lambda, GK must also.");
        lambdaTerm = true;
      }
    }

    initAtomArrays();

    double tensionDefault;
    switch (nonPolarModel) {
      case CAV:
        tensionDefault = DEFAULT_CAVONLY_SURFACE_TENSION;
        surfaceAreaRegion = new SurfaceAreaRegion(atoms, x, y, z, use,
            neighborLists, grad, threadCount, probe, tensionDefault);
        dispersionRegion = null;
        chandlerCavitation = null;
        break;
      case CAV_DISP:
        tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
        surfaceAreaRegion = new SurfaceAreaRegion(atoms, x, y, z, use,
            neighborLists, grad, threadCount, probe, tensionDefault);
        dispersionRegion = new DispersionRegion(threadCount, atoms, forceField);
        chandlerCavitation = null;
        break;
      case SEV_DISP:
        tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
        double[] radii = new double[nAtoms];
        int index = 0;
        for (Atom atom : atoms) {
          radii[index] = atom.getVDWType().radius / 2.0;
          index++;
        }
        ConnollyRegion connollyRegion = new ConnollyRegion(atoms, radii, threadCount);
        double wiggle = forceField.getDouble("WIGGLE", ConnollyRegion.DEFAULT_WIGGLE);
        connollyRegion.setWiggle(wiggle);
        chandlerCavitation = new ChandlerCavitation(atoms, connollyRegion, forceField);
        dispersionRegion = new DispersionRegion(threadCount, atoms, forceField);
        surfaceAreaRegion = null;
        break;
      case GAUSS_DISP:
        dispersionRegion = new DispersionRegion(threadCount, atoms, forceField);
        surfaceAreaRegion = null;
        tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
        GaussVol gaussVol = new GaussVol(atoms, forceField, parallelTeam);
        chandlerCavitation = new ChandlerCavitation(atoms, gaussVol, forceField);
        break;
      case BORN_CAV_DISP:
        tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
        surfaceAreaRegion = null;
        chandlerCavitation = null;
        dispersionRegion = new DispersionRegion(threadCount, atoms, forceField);
        break;
      case HYDROPHOBIC_PMF:
      case BORN_SOLV:
      case NONE:
      default:
        tensionDefault = DEFAULT_CAVONLY_SURFACE_TENSION;
        surfaceAreaRegion = null;
        dispersionRegion = null;
        chandlerCavitation = null;
        break;
    }

    gkc = forceField.getDouble("GKC", DEFAULT_GKC);

    surfaceTension = forceField.getDouble("SURFACE_TENSION", tensionDefault);
    solventPressue = forceField.getDouble("SOLVENT_PRESSURE", DEFAULT_SOLVENT_PRESSURE);
    // Volume to surface area cross-over point (A).
    double crossOver = forceField.getDouble("CROSS_OVER", DEFAULT_CROSSOVER);
    if (chandlerCavitation != null) {
      // Set the cross-over first.
      chandlerCavitation.setCrossOver(crossOver);
      // Surface tension and solvent pressure will over-write cross-over if it's not appropriate.
      chandlerCavitation.setSurfaceTension(surfaceTension);
      chandlerCavitation.setSolventPressure(solventPressue);
      crossOver = chandlerCavitation.getCrossOver();
    }
    if (surfaceAreaRegion != null) {
      surfaceAreaRegion.setSurfaceTension(surfaceTension);
    }
    double dispersionOffset = forceField.getDouble("DISPERSION_OFFSET", DEFAULT_DISPERSION_OFFSET);
    if (dispersionRegion != null) {
      dispersionRegion.setDispersionOffest(dispersionOffset);
    }

    initializationRegion = new InitializationRegion(threadCount);
    bornRadiiRegion = new BornRadiiRegion(threadCount, nAtoms, forceField, neckCorrection,
        tanhCorrection, perfectHCTScale);
    permanentGKFieldRegion = new PermanentGKFieldRegion(threadCount, soluteDielectric, solventDielectric, gkc);
    inducedGKFieldRegion = new InducedGKFieldRegion(threadCount, soluteDielectric, solventDielectric, gkc);
    if (!perfectRadii) {
      bornGradRegion = new BornGradRegion(threadCount, neckCorrection, tanhCorrection, perfectHCTScale);
    } else {
      // No Born chain-rule terms when using Perfect Born Radii.
      bornGradRegion = null;
    }
    boolean gkQI = forceField.getBoolean("GK_QI", false);
    gkEnergyRegion = new GKEnergyRegion(threadCount, polarization, nonPolarModel, surfaceTension,
        probe, electric, soluteDielectric, solventDielectric, gkc, gkQI);

    logger.info("  Continuum Solvation ");
    logger.info(format("   Radii:                              %8s", soluteRadiiType));
    if (cutoff != Double.POSITIVE_INFINITY) {
      logger.info(format("   Generalized Kirkwood Cut-Off:       %8.4f (A)", cutoff));
    } else {
      logger.info("   Generalized Kirkwood Cut-Off:           NONE");
    }
    logger.info(format("   GKC:                                %8.4f", gkc));
    logger.info(format("   Solvent Dielectric:                 %8.4f", solventDielectric));
    logger.info(format("   Solute Dielectric:                  %8.4f", soluteDielectric));

    if (perfectRadii) {
      logger.info(format("   Use Perfect Born Radii:             %8B", perfectRadii));
    } else {
      logger.info(format("   Descreen with vdW Radii:            %8B", descreenVDW));
      logger.info(format("   Descreen with Hydrogen Atoms:       %8B", descreenHydrogen));
      logger.info(format("   Descreen Offset:                    %8.4f", descreenOffset));
      if (neckCorrection) {
        logger.info(format("   Use Neck Correction:                %8B", neckCorrection));
        logger.info(format("   Sneck Scale Factor:                 %8.4f", sneck));
        logger.info(format("   Chemically Aware Sneck:             %8B", chemicallyAwareSneck));
      }
      if (tanhCorrection) {
        logger.info(format("   Use Tanh Correction                 %8B", tanhCorrection));
        logger.info(format("    Beta0:                             %8.4f", beta0));
        logger.info(format("    Beta1:                             %8.4f", beta1));
        logger.info(format("    Beta2:                             %8.4f", beta2));
      }
      if (perfectHCTScale) {
        logger.info(format("   GaussVol HCT Scale Factors:         %8B", perfectHCTScale));
      }
      logger.info(format("   General HCT Scale Factor:           %8.4f",
          forceField.getDouble("HCT-SCALE", DEFAULT_HCT_SCALE)));
      if (elementHCTScale) {
        logger.info(format("   Element-Specific HCT Scale Factors: %8B", elementHCTScale));
        Integer[] elementHCTkeyset = elementHCTScaleFactors.keySet().toArray(new Integer[0]);
        Arrays.sort(elementHCTkeyset);
        for(Integer key : elementHCTkeyset){
          logger.info(format("    HCT-Element # %d:                   %8.4f",key,elementHCTScaleFactors.get(key)));
        }
      }
    }

    logger.info(
        format("   Non-Polar Model:                  %10s",
            nonPolarModel.toString().replace('_', '-')));

    if (nonPolarModel.equals(NonPolarModel.GAUSS_DISP)) {
      logger.info(
          format("    GaussVol Radii Offset:               %2.4f",
              forceField.getDouble("GAUSSVOL_RADII_OFFSET", 0.0)));
      logger.info(
          format("    GaussVol Radii Scale:                %2.4f",
              forceField.getDouble("GAUSSVOL_RADII_SCALE", 1.0)));
    }

    if (dispersionRegion != null) {
      logger.info(
          format(
              "   Dispersion Integral Offset:         %8.4f (A)",
              dispersionRegion.getDispersionOffset()));
    }

    if (surfaceAreaRegion != null) {
      logger.info(format("   Cavitation Probe Radius:            %8.4f (A)", probe));
      logger.info(
          format("   Cavitation Surface Tension:         %8.4f (Kcal/mol/A^2)", surfaceTension));
    } else if (chandlerCavitation != null) {
      logger.info(
          format("   Cavitation Solvent Pressure:        %8.4f (Kcal/mol/A^3)", solventPressue));
      logger.info(
          format("   Cavitation Surface Tension:         %8.4f (Kcal/mol/A^2)", surfaceTension));
      logger.info(format("   Cavitation Cross-Over Radius:       %8.4f (A)", crossOver));
    }

    // Print out all Base Radii
    if (logger.isLoggable(Level.FINE)) {
      logger.fine("      Base Radii  Descreen Radius  Overlap Scale  Neck Scale");
      for (int i = 0; i < nAtoms; i++) {
        logger.info(
            format("   %s %8.6f %8.6f %5.3f %5.3f",
                atoms[i].toString(), baseRadius[i], descreenRadius[i], overlapScale[i],
                neckScale[i]));
      }
    }
  }

  /**
   * GK is using perfect radii where available.
   *
   * @return True if using perfect radii.
   */
  public boolean getUsePerfectRadii() {
    return perfectRadii;
  }

  /**
   * Return perfect Born radii read in as keywords, or base radii if perfect radii are not
   * available.
   *
   * @return Array of perfect Born radii.
   */
  public double[] getPerfectRadii() {
    bornRadiiRegion.init(atoms, crystal, sXYZ, neighborLists, baseRadius, descreenRadius,
        overlapScale, neckScale, descreenOffset, use, cut2, nativeEnvironmentApproximation, born);
    return bornRadiiRegion.getPerfectRadii();
  }

  /**
   * getNonPolarModel.
   *
   * @param nonpolarModel a {@link java.lang.String} object.
   * @return a {@link NonPolarModel} object.
   */
  public static NonPolarModel getNonPolarModel(String nonpolarModel) {
    try {
      return NonPolarModel.valueOf(toEnumForm(nonpolarModel));
    } catch (IllegalArgumentException ex) {
      logger.warning(" Unrecognized nonpolar model requested; defaulting to NONE.");
      return NonPolarModel.NONE;
    }
  }

  /**
   * computeBornRadii
   */
  public void computeBornRadii() {
    // The solute radii can change during titration based rotamer optimization.
    applySoluteRadii();

    // Update tanh correction parameters.
    if (tanhCorrection) {
      BornTanhRescaling.setBeta0(beta0);
      BornTanhRescaling.setBeta1(beta1);
      BornTanhRescaling.setBeta2(beta2);
    }

    try {
      bornRadiiRegion.init(
          atoms,
          crystal,
          sXYZ,
          neighborLists,
          baseRadius,
          descreenRadius,
          overlapScale,
          neckScale,
          descreenOffset,
          use,
          cut2,
          nativeEnvironmentApproximation,
          born);
      parallelTeam.execute(bornRadiiRegion);
    } catch (Exception e) {
      String message = "Fatal exception computing Born radii.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * computeInducedGKField
   */
  public void computeInducedGKField() {
    try {
      fieldGK.reset(parallelTeam);
      fieldGKCR.reset(parallelTeam);
      inducedGKFieldRegion.init(atoms, inducedDipole, inducedDipoleCR, crystal, sXYZ,
          neighborLists, use, cut2, born, fieldGK, fieldGKCR);
      parallelTeam.execute(inducedGKFieldRegion);
      fieldGK.reduce(parallelTeam);
      fieldGKCR.reduce(parallelTeam);
    } catch (Exception e) {
      String message = "Fatal exception computing induced GK field.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * computePermanentGKField
   */
  public void computePermanentGKField() {
    try {
      fieldGK.reset(parallelTeam);
      permanentGKFieldRegion.init(
          atoms, globalMultipole, crystal, sXYZ, neighborLists, use, cut2, born, fieldGK);
      parallelTeam.execute(permanentGKFieldRegion);
      fieldGK.reduce(parallelTeam);
    } catch (Exception e) {
      String message = "Fatal exception computing permanent GK field.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * getBaseRadii.
   *
   * @return an array of {@link double} objects.
   */
  public double[] getBaseRadii() {
    return baseRadius;
  }

  /**
   * getDescreenRadii.
   *
   * @return an array of {@link double} objects.
   */
  public double[] getDescreenRadii() {
    return descreenRadius;
  }

  /**
   * Returns the cavitation component of the solvation energy.
   *
   * @return Cavitation energy
   */
  public double getCavitationEnergy() {
    return cavitationEnergy;
  }

  /**
   * Return the GaussVol instance.
   *
   * @return GaussVol.
   */
  public ChandlerCavitation getChandlerCavitation() {
    return chandlerCavitation;
  }

  /**
   * Getter for the field <code>cutoff</code>.
   *
   * @return a double.
   */
  public double getCutoff() {
    return cutoff;
  }

  /**
   * Setter for the field <code>cutoff</code>.
   *
   * @param cutoff a double.
   */
  public void setCutoff(double cutoff) {
    this.cutoff = cutoff;
    this.cut2 = cutoff * cutoff;
  }

  /**
   * Returns the dielectric offset (in Angstroms).
   *
   * @return Currently: 0.09 Angstroms.
   */
  public double getDielecOffset() {
    return dOffset;
  }

  /**
   * Return the descreening dielectric offset.
   *
   * @return The offset (A).
   */
  public double getDescreenOffset() {
    return descreenOffset;
  }

  /**
   * Returns the dispersion component of the solvation energy.
   *
   * @return Dispersion energy
   */
  public double getDispersionEnergy() {
    return dispersionEnergy;
  }

  public DispersionRegion getDispersionRegion() {
    return dispersionRegion;
  }

  public AtomicDoubleArray3D getFieldGK() {
    return fieldGK;
  }

  public AtomicDoubleArray3D getFieldGKCR() {
    return fieldGKCR;
  }

  public SurfaceAreaRegion getSurfaceAreaRegion() {
    return surfaceAreaRegion;
  }

  /**
   * Returns the GK component of the solvation energy.
   *
   * @return GK electrostatic energy
   */
  public double getGeneralizedKirkwoordEnergy() {
    return gkEnergy;
  }

  /**
   * Returns the GK component of the solvation energy.
   *
   * @return GK electrostatic energy
   */
  public double getGeneralizedKirkwoordPermanentEnergy() {
    return gkPermanentEnergy;
  }

  /**
   * Returns the GK component of the solvation energy.
   *
   * @return GK electrostatic energy
   */
  public double getGeneralizedKirkwoordPolariztionEnergy() {
    return gkPolarizationEnergy;
  }

  public AtomicDoubleArray3D getGrad() {
    return grad;
  }

  /**
   * getInteractions
   *
   * @return a int.
   */
  public int getInteractions() {
    return gkEnergyRegion.getInteractions();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getLambda() {
    return lambda;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Updates the value of lPow.
   */
  @Override
  public void setLambda(double lambda) {
    if (lambdaTerm) {
      this.lambda = lambda;
    } else {
      // If the lambdaTerm flag is false, lambda is set to one.
      this.lambda = 1.0;
      lPow = 1.0;
      dlPow = 0.0;
      dl2Pow = 0.0;
    }
  }

  /**
   * Checks whether GK uses the Native Environment Approximation.
   *
   * <p>This (previously known as born-use-all) is useful for rotamer optimization under continuum
   * solvent. If a large number of sidechain atoms are completely removed from the GK/GB calculation,
   * the remaining sidechains are overly solvated. The NEA says "we will keep all sidechains not
   * under optimization in some default state and let them contribute to Born radii calculations, but
   * still exclude their solvation energy components."
   *
   * @return Whether the NEA is in use.
   */
  public boolean getNativeEnvironmentApproximation() {
    return nativeEnvironmentApproximation;
  }

  /**
   * getNonPolarModel.
   *
   * @return a {@link NonPolarModel} object.
   */
  public NonPolarModel getNonPolarModel() {
    return nonPolarModel;
  }

  /**
   * Getter for the field <code>overlapScale</code>.
   *
   * @return an array of {@link double} objects.
   */
  public double[] getOverlapScale() {
    return overlapScale;
  }

  /**
   * Getter for the field <code>neckScale</code>.
   *
   * @return an array of {@link double} objects.
   */
  public double[] getNeckScale() {
    return neckScale;
  }

  public boolean getTanhCorrection() {
    return tanhCorrection;
  }

  /**
   * Returns the probe radius (typically 1.4 Angstroms).
   *
   * @return Radius of the solvent probe.
   */
  public double getProbeRadius() {
    return probe;
  }

  /**
   * Returns the solvent relative permittivity (typically 78.3).
   *
   * @return Relative permittivity of the solvent.
   */
  public double getSolventPermittivity() {
    return solventDielectric;
  }

  /**
   * Returns the solvent relative permittivity (typically 1.0).
   *
   * @return Relative permittivity of the solute.
   */
  public double getSolutePermittivity() {
    return soluteDielectric;
  }

  /**
   * Getter for the field <code>surfaceTension</code>.
   *
   * @return a double.
   */
  public double getSurfaceTension() {
    return surfaceTension;
  }

  public AtomicDoubleArray3D getTorque() {
    return torque;
  }

  /**
   * {@inheritDoc}
   *
   * <p>The 2nd derivative is 0.0. (U=Lambda*Egk, dU/dL=Egk, d2U/dL2=0.0)
   */
  @Override
  public double getd2EdL2() {
    if (lambdaTerm) {
      return dl2Pow * solvationEnergy;
    } else {
      return 0.0;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    if (lambdaTerm) {
      return dlPow * solvationEnergy;
    }
    return 0.0;
  }

  /**
   * {@inheritDoc}
   *
   * <p>These contributions are already aggregated into the arrays used by PME.
   */
  @Override
  public void getdEdXdL(double[] gradient) {
    // This contribution is collected by GeneralizedKirkwood.reduce
  }

  public void init() {
    initializationRegion.init(this, atoms, lambdaTerm, grad, torque, bornRadiiChainRule);
    initializationRegion.executeWith(parallelTeam);
  }

  public void reduce(
      AtomicDoubleArray3D g,
      AtomicDoubleArray3D t,
      AtomicDoubleArray3D lg,
      AtomicDoubleArray3D lt) {
    grad.reduce(parallelTeam);
    torque.reduce(parallelTeam);
    for (int i = 0; i < nAtoms; i++) {
      g.add(0, i, lPow * grad.getX(i), lPow * grad.getY(i), lPow * grad.getZ(i));
      t.add(0, i, lPow * torque.getX(i), lPow * torque.getY(i), lPow * torque.getZ(i));
      if (lambdaTerm) {
        lg.add(0, i, dlPow * grad.getX(i), dlPow * grad.getY(i), dlPow * grad.getZ(i));
        lt.add(0, i, dlPow * torque.getX(i), dlPow * torque.getY(i), dlPow * torque.getZ(i));
      }
    }
  }

  public AtomicDoubleArray getSelfEnergy() {
    // Init and execute gkEnergyRegion (?)
    return gkEnergyRegion.getSelfEnergy();
  }

  public double[] getBorn() {
    return bornRadiiRegion.getBorn();
  }

  /**
   * Setter for element-specific HCT overlap scale factors
   *
   * @param elementHCT HashMap containing element name keys and scale factor values
   */
  public void setElementHCTScaleFactors(HashMap<Integer, Double> elementHCT) {
    for (Integer element : elementHCT.keySet()) {
      if (elementHCTScaleFactors.containsKey(element)) {
        elementHCTScaleFactors.replace(element, elementHCT.get(element));
      } else {
        logger.info("No HCT scale factor set for element: " + element);
      }
    }
    initAtomArrays();
  }

  /**
   * Setter for the field <code>atoms</code>.
   *
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   */
  public void setAtoms(Atom[] atoms) {
    this.atoms = atoms;
    nAtoms = atoms.length;
    maxNumAtoms = max(nAtoms, maxNumAtoms);
    initAtomArrays();
  }

  /**
   * Setter for the field <code>crystal</code>.
   *
   * @param crystal a {@link ffx.crystal.Crystal} object.
   */
  public void setCrystal(Crystal crystal) {
    this.crystal = crystal;
  }

  /**
   * setNeighborList.
   *
   * @param neighbors an array of {@link int} objects.
   */
  public void setNeighborList(int[][][] neighbors) {
    this.neighborLists = neighbors;
  }

  /**
   * Setter for the field <code>use</code>.
   *
   * @param use an array of {@link boolean} objects.
   */
  public void setUse(boolean[] use) {
    this.use = use;
  }

  /**
   * solvationEnergy
   *
   * @param gradient a boolean.
   * @param print a boolean.
   * @return a double.
   */
  public double solvationEnergy(boolean gradient, boolean print) {
    return solvationEnergy(0.0, gradient, print);
  }

  /**
   * solvationEnergy
   *
   * @param gkInducedCorrectionEnergy GK vacuum to SCRF polarization energy cost.
   * @param gradient a boolean.
   * @param print a boolean.
   * @return a double.
   */
  public double solvationEnergy(double gkInducedCorrectionEnergy, boolean gradient, boolean print) {
    cavitationEnergy = 0.0;
    dispersionEnergy = 0.0;
    gkEnergy = 0.0;

    try {
      // Find the GK energy.
      gkTime = -System.nanoTime();
      gkEnergyRegion.init(
          atoms,
          globalMultipole,
          inducedDipole,
          inducedDipoleCR,
          crystal,
          sXYZ,
          neighborLists,
          use,
          cut2,
          baseRadius,
          born,
          gradient,
          parallelTeam,
          grad,
          torque,
          bornRadiiChainRule);
      parallelTeam.execute(gkEnergyRegion);
      gkTime += System.nanoTime();

      // Find the nonpolar energy.
      switch (nonPolarModel) {
        case CAV:
          cavitationTime = -System.nanoTime();
          parallelTeam.execute(surfaceAreaRegion);
          cavitationEnergy = surfaceAreaRegion.getEnergy();
          cavitationTime += System.nanoTime();
          break;
        case CAV_DISP:
          dispersionTime = -System.nanoTime();
          dispersionRegion.init(atoms, crystal, use, neighborLists, x, y, z, cut2, gradient, grad);
          parallelTeam.execute(dispersionRegion);
          dispersionEnergy = dispersionRegion.getEnergy();
          dispersionTime += System.nanoTime();
          cavitationTime = -System.nanoTime();
          parallelTeam.execute(surfaceAreaRegion);
          cavitationEnergy = surfaceAreaRegion.getEnergy();
          cavitationTime += System.nanoTime();
          break;
        case SEV_DISP:
          dispersionTime = -System.nanoTime();
          dispersionRegion.init(atoms, crystal, use, neighborLists, x, y, z, cut2, gradient, grad);
          parallelTeam.execute(dispersionRegion);
          dispersionEnergy = dispersionRegion.getEnergy();
          dispersionTime += System.nanoTime();
          cavitationTime = -System.nanoTime();
          cavitationTime = -System.nanoTime();
          double[][] positions = new double[nAtoms][3];
          for (int i = 0; i < nAtoms; i++) {
            positions[i][0] = atoms[i].getX();
            positions[i][1] = atoms[i].getY();
            positions[i][2] = atoms[i].getZ();
          }
          chandlerCavitation.energyAndGradient(positions, grad);
          cavitationEnergy = chandlerCavitation.getEnergy();
          cavitationTime += System.nanoTime();
          break;
        case GAUSS_DISP:
          dispersionTime = -System.nanoTime();
          dispersionRegion.init(atoms, crystal, use, neighborLists, x, y, z, cut2, gradient, grad);
          parallelTeam.execute(dispersionRegion);
          dispersionEnergy = dispersionRegion.getEnergy();
          dispersionTime += System.nanoTime();
          cavitationTime = -System.nanoTime();
          positions = new double[nAtoms][3];
          for (int i = 0; i < nAtoms; i++) {
            positions[i][0] = atoms[i].getX();
            positions[i][1] = atoms[i].getY();
            positions[i][2] = atoms[i].getZ();
          }
          chandlerCavitation.energyAndGradient(positions, grad);
          cavitationEnergy = chandlerCavitation.getEnergy();
          cavitationTime += System.nanoTime();
          break;
        case BORN_CAV_DISP:
          dispersionTime = -System.nanoTime();
          dispersionRegion.init(atoms, crystal, use, neighborLists, x, y, z, cut2, gradient, grad);
          parallelTeam.execute(dispersionRegion);
          dispersionEnergy = dispersionRegion.getEnergy();
          dispersionTime += System.nanoTime();
          break;
        case HYDROPHOBIC_PMF:
        case BORN_SOLV:
        case NONE:
        default:
          break;
      }
    } catch (Exception e) {
      String message = "Fatal exception computing the continuum solvation energy.";
      logger.log(Level.SEVERE, message, e);
    }

    // Compute the Born radii chain rule term.
    if (gradient && !perfectRadii) {
      try {
        gkTime -= System.nanoTime();
        bornGradRegion.init(
            atoms, crystal, sXYZ, neighborLists,
            baseRadius, descreenRadius, overlapScale, neckScale, descreenOffset,
            bornRadiiRegion.getUnscaledBornIntegral(), use, cut2,
            nativeEnvironmentApproximation, born, grad, bornRadiiChainRule);
        bornGradRegion.executeWith(parallelTeam);
        gkTime += System.nanoTime();
      } catch (Exception e) {
        String message = "Fatal exception computing Born radii chain rule term.";
        logger.log(Level.SEVERE, message, e);
      }
    }

    gkEnergy = gkEnergyRegion.getEnergy() + gkInducedCorrectionEnergy;
    gkPermanentEnergy = gkEnergyRegion.getPermanentEnergy();
    gkPolarizationEnergy = gkEnergy - gkPermanentEnergy;

    // The following expression is equivalent to the former.
    // gkPolarizationEnergy = gkEnergyRegion.getPolarizationEnergy() + gkInducedCorrectionEnergy;

    // Solvation energy is the sum of cavitation, dispersion and GK
    solvationEnergy = cavitationEnergy + dispersionEnergy + gkEnergy;

    if (print) {
      logger.info(format(" Generalized Kirkwood%16.8f %10.3f", gkEnergy, gkTime * 1e-9));
      switch (nonPolarModel) {
        case CAV:
          logger.info(
              format(
                  " Cavitation          %16.8f %10.3f", cavitationEnergy, cavitationTime * 1e-9));
          break;
        case CAV_DISP:
        case SEV_DISP:
        case GAUSS_DISP:
          logger.info(
              format(
                  " Cavitation          %16.8f %10.3f", cavitationEnergy, cavitationTime * 1e-9));
          logger.info(
              format(
                  " Dispersion          %16.8f %10.3f", dispersionEnergy, dispersionTime * 1e-9));
          break;
        case BORN_CAV_DISP:
          logger.info(
              format(
                  " Dispersion          %16.8f %10.3f", dispersionEnergy, dispersionTime * 1e-9));
          break;
        case HYDROPHOBIC_PMF:
        case BORN_SOLV:
        case NONE:
        default:
          break;
      }
    }

    if (lambdaTerm) {
      return lPow * solvationEnergy;
    } else {
      return solvationEnergy;
    }
  }

  private void initAtomArrays() {
    sXYZ = particleMeshEwald.coordinates;

    x = sXYZ[0][0];
    y = sXYZ[0][1];
    z = sXYZ[0][2];

    globalMultipole = particleMeshEwald.globalMultipole;
    inducedDipole = particleMeshEwald.inducedDipole;
    inducedDipoleCR = particleMeshEwald.inducedDipoleCR;
    neighborLists = particleMeshEwald.neighborLists;

    if (grad == null) {
      int threadCount = parallelTeam.getThreadCount();
      grad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);
      torque = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);
      bornRadiiChainRule = atomicDoubleArrayFactory(atomicDoubleArrayImpl, threadCount, nAtoms);
    } else {
      grad.alloc(nAtoms);
      torque.alloc(nAtoms);
      bornRadiiChainRule.alloc(nAtoms);
    }

    fieldGK.alloc(nAtoms);
    fieldGKCR.alloc(nAtoms);

    if (baseRadius == null || baseRadius.length < nAtoms) {
      baseRadius = new double[nAtoms];
      descreenRadius = new double[nAtoms];
      overlapScale = new double[nAtoms];
      neckScale = new double[nAtoms];
      born = new double[nAtoms];
      use = new boolean[nAtoms];
    }

    fill(use, true);

    setSoluteRadii(forceField, atoms, soluteRadiiType);
    applySoluteRadii();

    if (dispersionRegion != null) {
      dispersionRegion.allocate(atoms);
    }

    if (surfaceAreaRegion != null) {
      surfaceAreaRegion.init();
    }

    if (chandlerCavitation != null) {
      GaussVol gaussVol = chandlerCavitation.getGaussVol();
      try {
        gaussVol.allocate(atoms);
      } catch (Exception e) {
        logger.severe(" " + e);
      }
    }
  }

  /**
   * Update GK solute parameters for a given atom. This should only be called after each atom is
   * assigned a "SoluteType".
   *
   * @param i The atom to update.
   */
  public void udpateSoluteParameters(int i) {
    Atom atom = atoms[i];
    SoluteType soluteType = atom.getSoluteType();
    if (soluteType == null) {
      logger.severe(format(" No SoluteType for atom %s.", atom));
      return;
    }

    // Assign the base radius.
    baseRadius[i] = soluteType.gkDiameter * 0.5 * bondiScale;
    // Assign a default overlap scale factor.
    overlapScale[i] = hctScale;
    // Use element specific HCT scaling factors.
    if (elementHCTScale) {
      int atomicNumber = atom.getAtomicNumber();
      if (elementHCTScaleFactors.containsKey(atomicNumber)) {
        overlapScale[i] = elementHCTScaleFactors.get(atomicNumber);
      } else {
        logger.fine(format(" No element-specific HCT scale factor for %s", atom));
        overlapScale[i] = hctScale;
      }
    }
    // Assign the default descreen radius to equal the base radius.
    descreenRadius[i] = baseRadius[i];

    // Handle radii override values.
    AtomType atomType = atom.getAtomType();
    if (radiiOverride.containsKey(atomType.type)) {
      double override = radiiOverride.get(atomType.type);
      // Remove default bondiFactor, and apply override.
      baseRadius[i] = baseRadius[i] * override / bondiScale;
      logger.fine(format(
          " Scaling %s (atom type %d) to %7.4f (Bondi factor %7.4f)",
          atom, atomType.type, baseRadius[i], override));
      descreenRadius[i] = baseRadius[i];
    }

    // Apply the descreenWithVDW flag.
    if (descreenVDW) {
      descreenRadius[i] = atom.getVDWType().radius / 2.0;
    }

    // Apply the descreenWithHydrogen flag.
    if (!descreenHydrogen && atom.getAtomicNumber() == 1) {
      overlapScale[i] = 0.0;
    }

    // Set Sneck scaling parameters based on atom type and number of bound non-hydrogen atoms.
    // The Sneck values for hydrogen atoms are controlled by their heavy atom.
    if (!atom.isHydrogen()) {
      // If the overlap scale factor is zero, then so is the neck overlap.
      if (!neckCorrection || overlapScale[i] == 0.0) {
        neckScale[i] = 0.0;
      } else {
        if (chemicallyAwareSneck) {
          // Determine number of bound heavy atoms for each non-hydrogen atom if chemically aware Sneck is being used
          int numBoundHeavyAtoms = 0;
          for (Atom boundAtom : atom.get12List()) {
            if (!boundAtom.isHydrogen()) {
              numBoundHeavyAtoms++;
            }
          }
          // Use this number to determine Sneck scaling parameter
          if (numBoundHeavyAtoms == 0) {
            // Sneck for lone ions or molecules like methane, which are not descreened by any other atoms
            neckScale[i] = 1.0;
          } else {
            neckScale[i] = atom.getSoluteType().sneck * (5.0 - numBoundHeavyAtoms) / 4.0;
          }
        } else {
          // Non-chemically aware Sneck - set neckScale to the max (input) Sneck value for all non-hydrogen atoms
          neckScale[i] = atom.getSoluteType().sneck;
        }
      }

      // Set hydrogen atom neck scaling factors to match those of their heavy atom binding partner.
      // By default, hydrogen atoms don't descreen other atoms. However, when they're descreened,
      // their contribution to the mixed neck value should matches their heavy atom.
      for (Atom boundAtom : atom.get12List()) {
        if (boundAtom.isHydrogen()) {
          int hydrogenIndex = boundAtom.getIndex();
          neckScale[hydrogenIndex - 1] = neckScale[i];
        }
      }
    }
  }

  /**
   * Apply solute radii definitions used to calculate Born radii.
   */
  public void applySoluteRadii() {
    // Set base radii, descreen radii, HCT overlap scale factor and neck scale factor.
    for (int i = 0; i < nAtoms; i++) {
      udpateSoluteParameters(i);
    }

    // Compute "perfect" HCT scale factors.
    if (perfectHCTScale) {
      double[][] coords = new double[nAtoms][3];
      int index = 0;
      for (Atom atom : atoms) {
        coords[index][0] = atom.getX();
        coords[index][1] = atom.getY();
        coords[index][2] = atom.getZ();
        index++;
      }
      GaussVol gaussVol = new GaussVol(atoms, forceField, parallelTeam);
      gaussVol.computeVolumeAndSA(coords);
      double[] selfVolumesFractions = gaussVol.getSelfVolumeFractions();
      for (int i = 0; i < nAtoms; i++) {
        // Use the self volume fractions, plus add the GK overlap scale.
        overlapScale[i] = selfVolumesFractions[i] * hctScale;

        // Apply the descreenWithHydrogen flag.
        if (!descreenHydrogen && atoms[i].getAtomicNumber() == 1) {
          overlapScale[i] = 0.0;
        }
      }
    }
  }

  public void setSneck(double sneck_input) {
    this.sneck = sneck_input;
    initAtomArrays();
  }

  public double[] getTanhBetas() {
    double[] betas = {beta0, beta1, beta2};
    return betas;
  }

  public void setTanhBetas(double[] betas) {
    this.beta0 = betas[0];
    this.beta1 = betas[1];
    this.beta2 = betas[2];
  }

  void setLambdaFunction(double lPow, double dlPow, double dl2Pow) {
    if (lambdaTerm) {
      this.lPow = lPow;
      this.dlPow = dlPow;
      this.dl2Pow = dl2Pow;
    } else {
      // If the lambdaTerm flag is false, lambda must be set to one.
      this.lambda = 1.0;
      this.lPow = 1.0;
      this.dlPow = 0.0;
      this.dl2Pow = 0.0;
    }
  }

  public enum NonPolarModel {
    CAV,
    CAV_DISP,
    SEV_DISP,
    GAUSS_DISP,
    HYDROPHOBIC_PMF,
    BORN_CAV_DISP,
    BORN_SOLV,
    NONE
  }
}
