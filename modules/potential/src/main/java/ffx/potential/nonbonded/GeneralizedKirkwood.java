// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import static ffx.potential.parameters.SoluteRadii.applyGKRadii;
import static ffx.utilities.Constants.DEFAULT_ELECTRIC;
import static ffx.utilities.Constants.dWater;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.pow;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.nonbonded.implicit.BornGradRegion;
import ffx.potential.nonbonded.implicit.BornRadiiRegion;
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
import ffx.potential.parameters.SoluteRadii;
import uk.ac.manchester.tornado.api.type.annotations.Atomic;

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
   * Default offset applied to radii for use with Gaussian Volumes to correct for not including
   * hydrogen atoms.
   */
  public static final double DEFAULT_GAUSSVOL_RADII_OFFSET = 0.0;
  /** Default dielectric offset */
  public static final double DEFAULT_DIELECTRIC_OFFSET = 0.09;
  /** Default constant for the Generalized Kirkwood cross-term. */
  public static final double DEFAULT_GKC = 2.455;

  private static final Logger logger = Logger.getLogger(GeneralizedKirkwood.class.getName());
  /** Default Bondi scale factor. */
  private static final double DEFAULT_SOLUTE_SCALE = 1.0;
  /**
   * Default overlap scale factor for the Hawkins, Cramer & Truhlar pairwise descreening algorithm: 0.69
   * New default overlap scale factor set during implicit solvent model optimization: 0.72
   */
  private static final double DEFAULT_HCT_SCALE = 0.72;
  /**
   * Default Sneck scaling factor from Aguilar/Onufriev 2010
   */
  private static final double DEFAULT_SNECK = 0.4058;
  /**
   * Default surface tension for apolar models without an explicit dispersion term. This is lower
   * than CAVDISP, since the favorable dispersion term is implicitly included.
   */
  private static final double DEFAULT_CAV_SURFACE_TENSION = 0.0049;
  /** Water probe radius. */
  public final double probe;
  /** Cavitation surface tension coefficient (kcal/mol/A^2). */
  private final double surfaceTension;
  /** Cavitation solvent pressure coefficient (kcal/mol/A^3). */
  private final double solventPressue;
  /**
   * Dielectric offset from:
   *
   * <p>W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson, "A Semianalytical Treatment of
   * Solvation for Molecular Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129 (1990)
   */
  private final double dOffset = DEFAULT_DIELECTRIC_OFFSET;
  /** Force field in use. */
  private final ForceField forceField;
  /** Treatment of polarization. */
  private final Polarization polarization;
  /** Treatment of non-polar interactions. */
  private final NonPolar nonPolar;
  /**
   * Particle mesh Ewald instance, which contains variables such as expanded coordinates and
   * multipoles in the global frame that GK uses.
   */
  private final ParticleMeshEwald particleMeshEwald;
  /** Parallel team object for shared memory parallelization. */
  private final ParallelTeam parallelTeam;
  /** Initialize GK output variables. */
  private final InitializationRegion initializationRegion;
  /** Parallel computation of Born Radii. */
  private final BornRadiiRegion bornRadiiRegion;
  /** Parallel computation of the Permanent GK Field. */
  private final PermanentGKFieldRegion permanentGKFieldRegion;
  /** Parallel computation of the Induced GK Field. */
  private final InducedGKFieldRegion inducedGKFieldRegion;
  /** Parallel computation of the GK continuum electrostatics energy. */
  private final GKEnergyRegion gkEnergyRegion;
  /** Parallel computation of Born radii chain rule term. */
  private final BornGradRegion bornGradRegion;
  /** Parallel computation of Dispersion energy. */
  private final DispersionRegion dispersionRegion;
  /** Parallel computation of Cavitation. */
  private final SurfaceAreaRegion surfaceAreaRegion;
  /** Volume to Surface Area Cavitation Dependence. */
  private final ChandlerCavitation chandlerCavitation;
  /**
   * Maps radii overrides (by AtomType) specified from the command line. e.g.
   * -DradiiOverride=134r1.20,135r1.20 sets atom types 134,135 to Bondi=1.20
   */
  private final HashMap<Integer, Double> radiiOverride = new HashMap<>();
  /** Conversion from electron**2/Ang to kcal/mole. */
  public double electric;
  /** Empirical scaling of the Bondi radii. */
  private double bondiScale;
  /** Volume to surface area cross-over point (A). */
  private double crossOver;
  /** GaussVol radii offset. */
  private double offset;
  /** The requested permittivity. */
  private double epsilon;
  /** Array of Atoms being considered. */
  private Atom[] atoms;
  /** Number of atoms. */
  private int nAtoms;
  /** Cartesian coordinates of each atom. */
  private double[] x, y, z;
  /** Base radius of each atom. */
  private double[] baseRadius;
  /** Descreen radius of each atom. */
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
   * Sneck scaling parameter for each atom. Set based on maximum Sneck scaling parameter and
   * number of bound non-hydrogen atoms
   */
  private double[] neckScale;
  /** If true, the descreening size of atoms is based on their force field vdW radius */
  private final boolean descreenWithVDW;
  /** If true, hydrogen atoms displace solvent */
  private final boolean descreenWithHydrogen;
  /** If true, the descreening integral includes overlaps with the volume of the descreened atom */
  private boolean perfectHCTScale;
  /** If true, the descreening integral includes the neck correction to better approximate molecular surface */
  private boolean neckCorrection;
  /** Maximum Sneck scaling parameter value */
  private double sneck;
  /** Base overlap scale factor. */
  private double gkOverlapScale;
  /** If true, HCT overlap scale factors are element-specific */
  private boolean elementHCTScale;
  /** Element-specific HCT overlap scale factors */
  private HashMap<String,Double> elementHCTScaleFactors;
  private double hct_n;
  private double hct_c;
  private double hct_h;
  private double hct_o;
  private double hct_p;
  private double hct_s;
  /** Born radius of each atom. */
  private double[] born;
  private double[] bornAfterTanh;
  /** Flag to indicate if an atom should be included. */
  private boolean[] use = null;
  /** Periodic boundary conditions and symmetry. */
  private Crystal crystal;
  /** Atomic coordinates for each symmetry operator. */
  private double[][][] sXYZ;
  /** Multipole moments for each symmetry operator. */
  private double[][][] globalMultipole;
  /** Induced dipoles for each symmetry operator. */
  private double[][][] inducedDipole;
  /** Induced dipole chain rule terms for each symmetry operator. */
  private double[][][] inducedDipoleCR;
  /** AtomicDoubleArray implementation to use. */
  private AtomicDoubleArrayImpl atomicDoubleArrayImpl;
  /** Atomic Gradient array. */
  private AtomicDoubleArray3D grad;
  /** Atomic Torque array. */
  private AtomicDoubleArray3D torque;
  /** Atomic Born radii chain-rule array. */
  private AtomicDoubleArray bornRadiiChainRule;
  /** Atomic GK field array. */
  private AtomicDoubleArray3D fieldGK;
  /** Atomic GK field chain-rule array. */
  private AtomicDoubleArray3D fieldGKCR;
  /** Neighbor lists for each atom and symmetry operator. */
  private int[][][] neighborLists;
  /**
   * This field is because re-initializing the force field resizes some arrays but not others; that
   * second category must, when called on, be resized not to the current number of atoms but to the
   * maximum number of atoms (and thus to the size of the first category of arrays).
   */
  private int maxNumAtoms;
  /** GK cut-off distance. */
  private double cutoff;
  /** GK cut-off distance squared. */
  private double cut2;
  /** Boolean flag to indicate GK will be scaled by the lambda state variable. */
  private boolean lambdaTerm;
  /** The current value of the lambda state variable. */
  private double lambda = 1.0;
  /**
   * lPow equals lambda^polarizationLambdaExponent, where polarizationLambdaExponent is also used by
   * PME.
   */
  private double lPow = 1.0;
  /** First derivative of lPow with respect to l. */
  private double dlPow = 0.0;
  /** Second derivative of lPow with respect to l. */
  private double dl2Pow = 0.0;
  /** Total Solvation Energy. */
  private double solvationEnergy = 0.0;
  /** Electrostatic Solvation Energy. */
  private double gkEnergy = 0.0;
  /** Dispersion Solvation Energy. */
  private double dispersionEnergy = 0.0;
  /** Cavitation Solvation Energy. */
  private double cavitationEnergy = 0.0;
  /** Time to compute GK electrostatics. */
  private long gkTime = 0;
  /** Time to compute Dispersion energy. */
  private long dispersionTime = 0;
  /** Time to compute Cavitation energy. */
  private long cavitationTime = 0;
  /** If true, prevents Born radii from updating. */
  private boolean fixedRadii = false;
  /** Forces all atoms to be considered during Born radius updates. */
  private boolean nativeEnvironmentApproximation;

  /**
   * Constructor for GeneralizedKirkwood.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param particleMeshEwald a {@link ffx.potential.nonbonded.ParticleMeshEwald} object.
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
   */
  public GeneralizedKirkwood(
      ForceField forceField,
      Atom[] atoms,
      ParticleMeshEwald particleMeshEwald,
      Crystal crystal,
      ParallelTeam parallelTeam) {

    this.forceField = forceField;
    this.atoms = atoms;
    this.particleMeshEwald = particleMeshEwald;
    this.crystal = crystal;
    this.parallelTeam = parallelTeam;
    nAtoms = atoms.length;
    maxNumAtoms = nAtoms;
    polarization = particleMeshEwald.polarization;

    // Set the conversion from electron**2/Ang to kcal/mole
    electric = forceField.getDouble("ELECTRIC", DEFAULT_ELECTRIC);

    // Set the Kirkwood multipolar reaction field constants.
    epsilon = forceField.getDouble("GK_EPSILON", dWater);

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

    // Define default Bondi scale factor, and HCT overlap scale factors.
    bondiScale = forceField.getDouble("SOLUTE_SCALE", DEFAULT_SOLUTE_SCALE);
    gkOverlapScale = forceField.getDouble("HCT_SCALE", DEFAULT_HCT_SCALE);
    elementHCTScale = forceField.getBoolean("ELEMENT_HCT_SCALE", false);
    descreenWithVDW = forceField.getBoolean("DESCREEN_VDW", true);
    descreenWithHydrogen = forceField.getBoolean("DESCREEN_HYDROGEN", false);
    if (descreenWithVDW && !descreenWithHydrogen) {
      perfectHCTScale = forceField.getBoolean("PERFECT_HCT_SCALE", false);
    } else {
      perfectHCTScale = false;
    }
    neckCorrection = forceField.getBoolean("NECK_CORRECTION",false);
    sneck = forceField.getDouble("SNECK",DEFAULT_SNECK);
    elementHCTScaleFactors = new HashMap<>();
    hct_n = forceField.getDouble("HCT_N",DEFAULT_HCT_SCALE);
    hct_c = forceField.getDouble("HCT_C",DEFAULT_HCT_SCALE);
    hct_h = forceField.getDouble("HCT_H",DEFAULT_HCT_SCALE);
    hct_o = forceField.getDouble("HCT_O",DEFAULT_HCT_SCALE);
    hct_p = forceField.getDouble("HCT_P",DEFAULT_HCT_SCALE);
    hct_s = forceField.getDouble("HCT_S",DEFAULT_HCT_SCALE);
    // Add default values for all elements
    elementHCTScaleFactors.put("N",hct_n);
    elementHCTScaleFactors.put("C",hct_c);
    elementHCTScaleFactors.put("H",hct_h);
    elementHCTScaleFactors.put("O",hct_o);
    elementHCTScaleFactors.put("P",hct_p);
    elementHCTScaleFactors.put("S",hct_s);

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

    NonPolar nonpolarModel;
    try {
      String cavModel = forceField.getString("CAVMODEL", "CAV_DISP").toUpperCase();
      nonpolarModel = getNonPolarModel(cavModel);
    } catch (Exception ex) {
      nonpolarModel = NonPolar.NONE;
      logger.warning(format(" Error parsing non-polar model (set to NONE) %s", ex.toString()));
    }
    nonPolar = nonpolarModel;

    int threadCount = parallelTeam.getThreadCount();
    fieldGK = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);
    fieldGKCR = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, threadCount);

    nativeEnvironmentApproximation =
        forceField.getBoolean("NATIVE_ENVIRONMENT_APPROXIMATION", false);
    probe = forceField.getDouble("PROBE_RADIUS", 1.4);
    cutoff = forceField.getDouble("GK_CUTOFF", particleMeshEwald.getEwaldCutoff());
    cut2 = cutoff * cutoff;
    lambdaTerm = forceField.getBoolean("GK_LAMBDATERM", forceField.getBoolean("LAMBDATERM", false));

    /*
     If polarization lambda exponent is set to 0.0, then we're running
     Dual-Topology and the GK energy will be scaled with the overall system lambda value.
    */
    double polLambdaExp = forceField.getDouble("POLARIZATION_LAMBDA_EXPONENT", 3.0);
    if (polLambdaExp == 0.0) {
      lambdaTerm = false;
      logger.info(" GK lambda term set to false.");
    }

    // If PME includes polarization and is a function of lambda, GK must also.
    if (!lambdaTerm && particleMeshEwald.getPolarizationType() != Polarization.NONE) {
      if (forceField.getBoolean("ELEC_LAMBDATERM", forceField.getBoolean("LAMBDATERM", false))) {
        logger.info(" If PME includes polarization and is a function of lambda, GK must also.");
        lambdaTerm = true;
      }
    }

    initAtomArrays();

    double tensionDefault;
    switch (nonPolar) {
      case CAV:
        tensionDefault = DEFAULT_CAV_SURFACE_TENSION;
        surfaceAreaRegion =
            new SurfaceAreaRegion(
                atoms, x, y, z, use, neighborLists, grad, threadCount, probe, tensionDefault);
        dispersionRegion = null;
        chandlerCavitation = null;
        break;
      case CAV_DISP:
        tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
        surfaceAreaRegion =
            new SurfaceAreaRegion(
                atoms, x, y, z, use, neighborLists, grad, threadCount, probe, tensionDefault);
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
        boolean[] isHydrogen = new boolean[nAtoms];
        radii = new double[nAtoms];
        double[] volume = new double[nAtoms];
        double[] gamma = new double[nAtoms];
        double fourThirdsPI = 4.0 / 3.0 * PI;
        offset = forceField.getDouble("GAUSSVOL_RADII_OFFSET", DEFAULT_GAUSSVOL_RADII_OFFSET);
        index = 0;
        for (Atom atom : atoms) {
          isHydrogen[index] = atom.isHydrogen();
          radii[index] = atom.getVDWType().radius / 2.0;
          radii[index] += offset;
          volume[index] = fourThirdsPI * pow(radii[index], 3);
          gamma[index] = 1.0;
          index++;
        }
        GaussVol gaussVol = new GaussVol(nAtoms, radii, volume, gamma, isHydrogen, parallelTeam);
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
        tensionDefault = DEFAULT_CAV_SURFACE_TENSION;
        surfaceAreaRegion = null;
        dispersionRegion = null;
        chandlerCavitation = null;
        break;
    }

    surfaceTension = forceField.getDouble("SURFACE_TENSION", tensionDefault);
    solventPressue = forceField.getDouble("SOLVENT_PRESSURE", DEFAULT_SOLVENT_PRESSURE);
    crossOver = forceField.getDouble("CROSS_OVER", DEFAULT_CROSSOVER);
    if (chandlerCavitation != null) {
      // Set the cross-over first.
      chandlerCavitation.setCrossOver(crossOver);
      // Surface tension and solvent pressure will over-write cross-over if its not appropriate.
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
    bornRadiiRegion = new BornRadiiRegion(threadCount, forceField, perfectHCTScale);
    permanentGKFieldRegion = new PermanentGKFieldRegion(threadCount, forceField);
    inducedGKFieldRegion = new InducedGKFieldRegion(threadCount, forceField);
    bornGradRegion = new BornGradRegion(threadCount, perfectHCTScale, neckCorrection, bornRadiiRegion.getTanhCorrectionBoolean());
    gkEnergyRegion =
        new GKEnergyRegion(threadCount, forceField, polarization, nonPolar, surfaceTension, probe);

    logger.info("  Continuum Solvation ");
    logger.info(format("   Generalized Kirkwood Cut-Off:       %8.3f (A)", cutoff));
    logger.info(format("   Solvent Dielectric:                 %8.3f", epsilon));
    logger.info(format("   Descreen with vdW Radii:            %8B", descreenWithVDW));
    logger.info(format("   Descreen with Hydrogen Atoms:       %8B", descreenWithHydrogen));
    logger.info(format("   Use Neck Correction:                %8B", neckCorrection));
    logger.info(format("   Sneck:                              %8.4f",sneck));
    logger.info(format("   Use Tanh Correction                 %8B",bornRadiiRegion.getTanhCorrectionBoolean()));
    if(bornRadiiRegion.getTanhCorrectionBoolean()){
      logger.info(format("    Beta0:                             %8.4f",bornRadiiRegion.getBeta0()));
      logger.info(format("    Beta1:                             %8.4f",bornRadiiRegion.getBeta1()));
      logger.info(format("    Beta2:                             %8.4f",bornRadiiRegion.getBeta2()));
    }
    logger.info(format("   GaussVol HCT Scale Factors:         %8B", perfectHCTScale));
    logger.info(format("   Element-Specific HCT Scale Factors: %8B", elementHCTScale));
    if(elementHCTScale){
      logger.info(format("    HCT-H:                             %8.4f",hct_h));
      logger.info(format("    HCT-C:                             %8.4f",hct_c));
      logger.info(format("    HCT-N:                             %8.4f",hct_n));
      logger.info(format("    HCT-O:                             %8.4f",hct_o));
      logger.info(format("    HCT-P:                             %8.4f",hct_p));
      logger.info(format("    HCT-S:                             %8.4f",hct_s));
    }
    logger.info(format("   HCT Scale Factor:                   %8.4f",forceField.getDouble("HCT-SCALE",DEFAULT_HCT_SCALE)));
    logger.info(format("   GKC:                                %8.3f",forceField.getDouble("GKC",DEFAULT_GKC)));
    SoluteRadii.logRadiiSource(forceField);
    logger.info(
        format("   Non-Polar Model:                  %10s", nonPolar.toString().replace('_', '-')));

    if (dispersionRegion != null) {
      logger.info(
          format(
              "   Dispersion Integral Offset:         %8.3f (A)",
              dispersionRegion.getDispersionOffset()));
    }

    if (surfaceAreaRegion != null) {
      logger.info(format("   Cavitation Probe Radius:            %8.3f (A)", probe));
      logger.info(
          format("   Cavitation Surface Tension:         %8.3f (Kcal/mol/A^2)", surfaceTension));
    } else if (chandlerCavitation != null) {
      logger.info(
          format("   Cavitation Solvent Pressure:        %8.3f (Kcal/mol/A^3)", solventPressue));
      logger.info(
          format("   Cavitation Surface Tension:         %8.3f (Kcal/mol/A^2)", surfaceTension));
      logger.info(format("   Cavitation Cross-Over Radius:       %8.3f (A)", crossOver));
    }

    // Print out all Base Radii
    if (logger.isLoggable(Level.FINE)) {
      logger.fine("   GK Base Radii  Descreen Radius  Overlap Scale");
      for (int i = 0; i < nAtoms; i++) {
        logger.info(
            format("   %s %8.6f %8.6f %5.3f",
                atoms[i].toString(), baseRadius[i], descreenRadius[i], overlapScale[i]));
      }
    }
  }

  /**
   * getNonPolarModel.
   *
   * @param nonpolarModel a {@link java.lang.String} object.
   * @return a {@link ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar} object.
   */
  public static NonPolar getNonPolarModel(String nonpolarModel) {
    try {
      return NonPolar.valueOf(toEnumForm(nonpolarModel));
    } catch (IllegalArgumentException ex) {
      logger.warning(" Unrecognized nonpolar model requested; defaulting to NONE.");
      return NonPolar.NONE;
    }
  }

  /** computeBornRadii */
  public void computeBornRadii() {
    // Born radii are fixed.
    if (fixedRadii) {
      return;
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

  /** computeInducedGKField */
  public void computeInducedGKField() {
    try {
      fieldGK.reset(parallelTeam, 0, nAtoms - 1);
      fieldGKCR.reset(parallelTeam, 0, nAtoms - 1);
      inducedGKFieldRegion.init(
          atoms,
          inducedDipole,
          inducedDipoleCR,
          crystal,
          sXYZ,
          neighborLists,
          use,
          cut2,
          born,
          fieldGK,
          fieldGKCR);
      parallelTeam.execute(inducedGKFieldRegion);
      fieldGK.reduce(parallelTeam, 0, nAtoms - 1);
      fieldGKCR.reduce(parallelTeam, 0, nAtoms - 1);
    } catch (Exception e) {
      String message = "Fatal exception computing induced GK field.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /** computePermanentGKField */
  public void computePermanentGKField() {
    try {
      fieldGK.reset(parallelTeam, 0, nAtoms - 1);
      permanentGKFieldRegion.init(
          atoms, globalMultipole, crystal, sXYZ, neighborLists, use, cut2, born, fieldGK);
      parallelTeam.execute(permanentGKFieldRegion);
      fieldGK.reduce(parallelTeam, 0, nAtoms - 1);
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

  /**
   * Returns the GK component of the solvation energy.
   *
   * @return GK electrostatic energy
   */
  public double getGeneralizedKirkwoordEnergy() {
    return gkEnergy;
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

  /** {@inheritDoc} */
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
   * @return a {@link ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar} object.
   */
  public NonPolar getNonPolarModel() {
    return nonPolar;
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
   * @return Relative permittivity of solvent.
   */
  public double getSolventPermittivity() {
    return epsilon;
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

  /** {@inheritDoc} */
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
    initializationRegion.init(atoms, lambdaTerm, grad, torque, bornRadiiChainRule);
    initializationRegion.executeWith(parallelTeam);
  }

  public void reduce(
      AtomicDoubleArray3D g,
      AtomicDoubleArray3D t,
      AtomicDoubleArray3D lg,
      AtomicDoubleArray3D lt) {
    grad.reduce(parallelTeam, 0, nAtoms - 1);
    torque.reduce(parallelTeam, 0, nAtoms - 1);
    for (int i = 0; i < nAtoms; i++) {
      g.add(0, i, lPow * grad.getX(i), lPow * grad.getY(i), lPow * grad.getZ(i));
      t.add(0, i, lPow * torque.getX(i), lPow * torque.getY(i), lPow * torque.getZ(i));
      if (lambdaTerm) {
        lg.add(0, i, dlPow * grad.getX(i), dlPow * grad.getY(i), dlPow * grad.getZ(i));
        lt.add(0, i, dlPow * torque.getX(i), dlPow * torque.getY(i), dlPow * torque.getZ(i));
      }
    }
  }

  public AtomicDoubleArray getSelfEnergy(){
    // Init and execute gkEnergyRegion (?)
    return gkEnergyRegion.getSelfEnergy();
  }

  public double[] getBorn(){
    return bornRadiiRegion.getBorn();
  }

  /**
   * Setter for element-specific HCT overlap scale factors
   * @param elementHCT HashMap containing element name keys and scale factor values
   */
  public void setElementHCTScaleFactors(HashMap<String,Double> elementHCT){
    for(String element : elementHCT.keySet()){
      if(elementHCTScaleFactors.containsKey(element.toUpperCase())){
        elementHCTScaleFactors.replace(element.toUpperCase(), elementHCT.get(element));
      } else {
        logger.info("No HCT scale factor set for element: "+element);
      }
    }
    initAtomArrays();
    //logger.info("Set HCT Element Scale Factors to: "+elementHCT.get("N")+" "+elementHCT.get("C")+" "+
      //      elementHCT.get("H")+" "+elementHCT.get("O")+" "+elementHCT.get("P")+" "+elementHCT.get("S"));
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
   * Update the GaussVol radii offset value.
   *
   * @param offset Offset to apply to all radii (A).
   */
  public void setGaussVolRadiiOffset(double offset) {
    this.offset = offset;
    GaussVol gaussVol = chandlerCavitation.getGaussVol();
    if (gaussVol != null) {
      double[] radii = new double[nAtoms];
      double[] volume = new double[nAtoms];
      double fourThirdsPI = 4.0 / 3.0 * PI;
      int index = 0;
      for (Atom atom : atoms) {
        radii[index] = atom.getVDWType().radius / 2.0 + offset;
        volume[index] = fourThirdsPI * pow(radii[index], 3);
        index++;
      }
      try {
        gaussVol.setRadiiAndVolumes(radii, volume);
      } catch (Exception e) {
        logger.severe(" Exception creating GaussVol: " + e.toString());
      }
    }
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
   * @param gkPolarizationEnergy GK vacuum to SCRF polarization energy cost.
   * @param gradient a boolean.
   * @param print a boolean.
   * @return a double.
   */
  public double solvationEnergy(double gkPolarizationEnergy, boolean gradient, boolean print) {

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
      switch (nonPolar) {
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
    if (gradient) {
      try {
        gkTime -= System.nanoTime();
        bornGradRegion.init(
            atoms,
            crystal,
            sXYZ,
            neighborLists,
            baseRadius,
            descreenRadius,
            overlapScale,
            neckScale,
            use,
            cut2,
            nativeEnvironmentApproximation,
            born,
            grad,
            bornRadiiChainRule);
        bornGradRegion.executeWith(parallelTeam);
        gkTime += System.nanoTime();
      } catch (Exception e) {
        String message = "Fatal exception computing Born radii chain rule term.";
        logger.log(Level.SEVERE, message, e);
      }
    }

    gkEnergy = gkEnergyRegion.getEnergy() + gkPolarizationEnergy;

    // Solvation energy is the sum of cavitation, dispersion and GK
    solvationEnergy = cavitationEnergy + dispersionEnergy + gkEnergy;

    if (print) {
      logger.info(format(" Generalized Kirkwood%16.8f %10.3f", gkEnergy, gkTime * 1e-9));
      switch (nonPolar) {
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
    if (fixedRadii) {
      fixedRadii = false;
    }

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

    applyGKRadii(forceField, bondiScale, atoms, baseRadius);

    // Set default HCT overlap scale factor.
    for (int i = 0; i < nAtoms; i++) {
      overlapScale[i] = gkOverlapScale;
      // Use element specific HCT scaling factors
      if(elementHCTScale){
        Atom atom = atoms[i];
        String atomName = atom.getName();
        char elementChar = (atomName.charAt(0));
        String elementName = Character.toString(elementChar);
        boolean notIon = true;
        if(atomName.contains("-")||atomName.contains("+")){
          logger.info("Ion found: "+atomName);
          notIon = false;
        }
        if(elementHCTScaleFactors.get(elementName) != null && notIon) {
          overlapScale[i] = elementHCTScaleFactors.get(elementName);
        } else{
          logger.info("Element-specific HCT scale factor not found for atom "+i+" of element "+atomName);
          overlapScale[i] = gkOverlapScale;
        }
      }
      descreenRadius[i] = baseRadius[i];
    }

    // Compute "perfect" HCT scale factors.
    if (perfectHCTScale) {
      boolean[] isHydrogen = new boolean[nAtoms];
      double[] radii = new double[nAtoms];
      double[] volume = new double[nAtoms];
      double[] gamma = new double[nAtoms];
      double[][] coords = new double[nAtoms][3];
      double fourThirdsPI = 4.0 / 3.0 * PI;
      int index = 0;
      for (Atom atom : atoms) {
        isHydrogen[index] = atom.isHydrogen();
        radii[index] = atom.getVDWType().radius / 2.0;
        volume[index] = fourThirdsPI * pow(radii[index], 3);
        gamma[index] = 1.0;
        coords[index][0] = atom.getX();
        coords[index][1] = atom.getY();
        coords[index][2] = atom.getZ();
        index++;
      }

      GaussVol gaussVol = new GaussVol(nAtoms, radii, volume, gamma, isHydrogen, parallelTeam);
      gaussVol.computeVolumeAndSA(coords);
      double[] selfVolumesFractions = gaussVol.getSelfVolumeFractions();
      for (int i=0; i<nAtoms; i++) {
        // Use the self volume fractions, plus add the GK overlap scale.
        overlapScale[i] = selfVolumesFractions[i] * gkOverlapScale;
      }
    }

    // Set up HCT overlap scale factors and any requested radii overrides.
    for (int i = 0; i < nAtoms; i++) {
      AtomType atomType = atoms[i].getAtomType();
      if (radiiOverride.containsKey(atomType.type)) {
        double override = radiiOverride.get(atomType.type);
        // Remove default bondiFactor, and apply override.
        baseRadius[i] = baseRadius[i] * override / bondiScale;
        logger.info(
            format(
                " Scaling %s (atom type %d) to %7.4f (Bondi factor %7.4f)",
                atoms[i], atomType.type, baseRadius[i], override));
        descreenRadius[i] = baseRadius[i];
      }

      // Apply the descreenWithVDW flag.
      if (descreenWithVDW) {
        descreenRadius[i] = atoms[i].getVDWType().radius / 2.0;

        // The descreening will be with a scaled size = vdW radius * overlapScale.
        // double vdwSize = atoms[i].getVDWType().radius / 2.0 * overlapScale[i];
        // The code will descreen with a fit baseRadius whose size is different from the vdW radius.
        // Define a per atom "overlapScale" to compensate.
        // perAtomScaleFactor * baseRadius = vdwSize
        // perAtomScaleFactor = vdwSize / baseRadius
        // overlapScale[i] = vdwSize / baseRadius[i];
      }

      // Apply the descreenWithHydrogen flag.
      if (!descreenWithHydrogen && atoms[i].getAtomicNumber() == 1) {
        overlapScale[i] = 0.0;
      }
    }

    // Set Sneck scaling parameters based on atom type and number of bound non-hydrogen atoms.
    for(int i = 0; i < nAtoms; i++){
      int numBoundHeavyAtoms = 0;
      // Hydrogen atoms aren't used in descreening.
      if(overlapScale[i] == 0.0){
        neckScale[i] = 0.0;
      } else {
        // Determine number of bound heavy atoms for each non-hydrogen atom.
        for(Atom boundAtom : atoms[i].get12List()){
          if(!boundAtom.isHydrogen()){
            numBoundHeavyAtoms++;
          }
        }
        // Use this number to determine Sneck scaling parameter
        if(numBoundHeavyAtoms == 0){
          // Sneck for lone ions or molecules like methane, which are not descreened by any other atoms
          neckScale[i] = 1.00;
        } else {
          neckScale[i] = sneck*(5.0-numBoundHeavyAtoms)/4.0;
        }
      }
    }

    if (dispersionRegion != null) {
      dispersionRegion.allocate(atoms);
    }

    if (surfaceAreaRegion != null) {
      surfaceAreaRegion.init();
    }

    if (chandlerCavitation != null) {
      GaussVol gaussVol = chandlerCavitation.getGaussVol();
      if (gaussVol != null) {
        boolean[] isHydrogen = new boolean[nAtoms];
        double[] radii = new double[nAtoms];
        double[] volume = new double[nAtoms];
        double[] gamma = new double[nAtoms];
        double fourThirdsPI = 4.0 / 3.0 * PI;
        int index = 0;
        for (Atom atom : atoms) {
          isHydrogen[index] = atom.isHydrogen();
          radii[index] = atom.getVDWType().radius / 2.0;
          radii[index] += offset;
          volume[index] = fourThirdsPI * pow(radii[index], 3);
          gamma[index] = 1.0;
          index++;
        }
        try {
          gaussVol.setGammas(gamma);
          gaussVol.setRadiiAndVolumes(radii, volume);
          gaussVol.setIsHydrogen(isHydrogen);
        } catch (Exception e) {
          logger.severe(" Exception creating GaussVol: " + e.toString());
        }
      }
    }
  }

  public void setSneck(double sneck_input){
    this.sneck = sneck_input;
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

  public enum NonPolar {
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
