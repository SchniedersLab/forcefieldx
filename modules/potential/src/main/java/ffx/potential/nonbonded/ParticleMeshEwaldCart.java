// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import static ffx.numerics.special.Erf.erfc;
import static ffx.potential.parameters.ForceField.ELEC_FORM.PAM;
import static ffx.potential.parameters.ForceField.toEnumForm;
import static ffx.potential.parameters.MultipoleType.assignMultipole;
import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t200;
import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.util.Range;
import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.math.ScalarMath;
import ffx.numerics.multipole.MultipoleTensor;
import ffx.potential.ForceFieldEnergy.Platform;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.ReciprocalSpace.FFTMethod;
import ffx.potential.nonbonded.pme.DirectRegion;
import ffx.potential.nonbonded.pme.ExpandInducedDipolesRegion;
import ffx.potential.nonbonded.pme.InducedDipoleFieldReduceRegion;
import ffx.potential.nonbonded.pme.InducedDipoleFieldRegion;
import ffx.potential.nonbonded.pme.InitializationRegion;
import ffx.potential.nonbonded.pme.OPTRegion;
import ffx.potential.nonbonded.pme.PCGSolver;
import ffx.potential.nonbonded.pme.PermanentFieldRegion;
import ffx.potential.nonbonded.pme.PolarizationEnergyRegion;
import ffx.potential.nonbonded.pme.RealSpaceEnergyRegion;
import ffx.potential.nonbonded.pme.ReciprocalEnergyRegion;
import ffx.potential.nonbonded.pme.ReduceRegion;
import ffx.potential.nonbonded.pme.SORRegion;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ELEC_FORM;
import ffx.potential.parameters.ForceField.ForceFieldType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.utils.EnergyException;
import ffx.utilities.Constants;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;

/**
 * This Particle Mesh Ewald class implements PME for the AMOEBA polarizable mutlipole force field in
 * parallel using a {@link NeighborList} for any {@link Crystal} space group. The real space
 * contribution is contained within this class and the reciprocal space contribution is delegated to
 * the {@link ReciprocalSpace} class.
 *
 * @author Michael J. Schnieders<br> derived from:<br> TINKER code by Jay Ponder, Pengyu Ren and Tom
 *     Darden.<br>
 * @see <br>
 *     <a href="http://dx.doi.org/10.1021/ct300035u" target="_blank"> M. J. Schnieders, J.
 *     Baltrusaitis, Y. Shi, G. Chattree, L. Zheng, W. Yang and P. Ren, The Structure,
 *     Thermodynamics, and Solubility of Organic Crystals from Simulation with a Polarizable Force
 *     Field, Journal of Chemical Theory and Computation 8 (5), 1721-36 (2012)</a>
 * @see <br>
 *     <a href="http://dx.doi.org/10.1021/ct100506d" target="_blank"> M. J. Schnieders, T. D. Fenn
 *     and V. S. Pande, Polarizable atomic multipole X-ray refinement: Particle-mesh Ewald
 *     electrostatics for macromolecular crystals. Journal of Chemical Theory and Computation 7 (4),
 *     1141-56 (2011)</a>
 * @see <br>
 *     <a href="http://dx.doi.org/10.1063/1.1630791" target="_blank"> C. Sagui, L. G. Pedersen, and
 *     T. A. Darden, Journal of Chemical Physics 120 (1), 73 (2004)</a>
 * @see <br>
 *     <a href="http://link.aip.org/link/?JCPSA6/98/10089/1" target="_©lank"> T. Darden, D. York,
 *     and L. Pedersen, Journal of Chemical Physics 98 (12), 10089 (1993)</a>
 * @see <br>
 *     <a href="http://www.ccp5.org" target="_blank"> W. Smith, "Point Multipoles in the Ewald
 *     Summation (Revisited)", CCP5 Newsletter, 46, 18-30 (1998)</a>
 */
@SuppressWarnings("deprecation")
public class ParticleMeshEwaldCart extends ParticleMeshEwald implements LambdaInterface {

  private static final Logger logger = Logger.getLogger(ParticleMeshEwald.class.getName());
  /** Number of unique tensors for given order. */
  private static final int tensorCount = MultipoleTensor.tensorCount(3);
  /** The sqrt of PI. */
  private static final double SQRT_PI = sqrt(Math.PI);
  /**
   * If lambdaTerm is true, some ligand atom interactions with the environment are being turned
   * on/off.
   */
  private final boolean lambdaTerm;

  private final boolean reciprocalSpaceTerm;
  /** Reference to the force field being used. */
  private final ForceField forceField;

  private final double poleps;
  /** Specify an SCF predictor algorithm. */
  private final SCFPredictor scfPredictor;
  private final SCFPredictorParameters scfPredictorParameters;
  private final EwaldParameters ewaldParameters;
  private final ScaleParameters scaleParameters;
  private final AlchemicalParameters alchemicalParameters;
  private final RealSpaceNeighborParameters realSpaceNeighborParameters;
  private final PMETimings pmeTimings;
  /** By default, maxThreads is set to the number of available SMP cores. */
  private final int maxThreads;
  /**
   * The default ParallelTeam encapsulates the maximum number of threads used to parallelize the
   * electrostatics calculation.
   */
  private final ParallelTeam parallelTeam;
  /**
   * The sectionTeam encapsulates 1 or 2 threads.
   *
   * <p>If it contains 1 thread, the real and reciprocal space calculations are done sequentially.
   *
   * <p>If it contains 2 threads, the real and reciprocal space calculations will be done
   * concurrently.
   */
  private final ParallelTeam sectionTeam;
  /**
   * If the real and reciprocal space parts of PME are done sequentially, then the realSpaceTeam is
   * equal parallalTeam.
   *
   * <p>If the real and reciprocal space parts of PME are done concurrently, then the realSpaceTeam
   * will have fewer threads than the default parallelTeam.
   */
  private final ParallelTeam realSpaceTeam;
  /**
   * If the real and reciprocal space parts of PME are done sequentially, then the
   * reciprocalSpaceTeam is equal parallalTeam.
   *
   * <p>If the real and reciprocal space parts of PME are done concurrently, then the
   * reciprocalSpaceTeam will have fewer threads than the default parallelTeam.
   */
  private final ParallelTeam fftTeam;

  private final NeighborList neighborList;
  private final InitializationRegion initializationRegion;
  private final PermanentFieldRegion permanentFieldRegion;
  private final InducedDipoleFieldRegion inducedDipoleFieldRegion;
  private final InducedDipoleFieldReduceRegion inducedDipoleFieldReduceRegion;
  private final ExpandInducedDipolesRegion expandInducedDipolesRegion;
  private final DirectRegion directRegion;
  private final SORRegion sorRegion;
  private final OPTRegion optRegion;
  private final PCGSolver pcgSolver;
  private final PolarizationEnergyRegion polarizationEnergyRegion;
  private final ReciprocalSpace reciprocalSpace;
  private final ReciprocalEnergyRegion reciprocalEnergyRegion;
  private final RealSpaceEnergyRegion realSpaceEnergyRegion;
  private final ReduceRegion reduceRegion;
  private final GeneralizedKirkwood generalizedKirkwood;
  /** Partial derivative with respect to Lambda. */
  private final SharedDouble shareddEdLambda;
  /** Second partial derivative with respect to Lambda. */
  private final SharedDouble sharedd2EdLambda2;
  /** The electrostatics functional form in use. */
  private final ELEC_FORM elecForm;
  /** Unit cell and spacegroup information. */
  private Crystal crystal;
  /** Number of symmetry operators. */
  private int nSymm;
  /** Flag to indicate use of generalized Kirkwood. */
  private boolean generalizedKirkwoodTerm;
  /** If true, compute coordinate gradient. */
  private boolean gradient = false;
  /**
   * If less than 1.0, all atoms with their lambda flags set will have their multipoles and
   * polarizabilities scaled.
   */
  private double lambdaScaleMultipoles = 1.0;
  /** Number of PME multipole interactions. */
  private int interactions;
  /** Number of generalized Kirkwood interactions. */
  private int gkInteractions;
  /** Generalized Kirkwood energy. */
  private double solvationEnergy;
  /** The current LambdaMode of this PME instance (or OFF for no lambda dependence). */
  private LambdaMode lambdaMode = LambdaMode.OFF;
  /** Current state. */
  private double lambda = 1.0;
  /** Permanent multipoles in their local frame. */
  private double[][] localMultipole;

  private MultipoleFrameDefinition[] frame;
  private int[][] axisAtom;
  private double[] ipdamp;
  private double[] thole;
  private double[] polarizability;
  /** 1-2, 1-3, 1-4 and 1-5 connectivity lists. */
  private int[][] mask12;
  private int[][] mask13;
  private int[][] mask14;
  private int[][] mask15;
  /** Flag for ligand atoms. */
  private boolean[] isSoft;
  /** Molecule number for each atom. */
  private int[] molecule;
  /**
   * When computing the polarization energy at Lambda there are 3 pieces.
   *
   * <p>1.) Upol(1) = The polarization energy computed normally (ie. system with ligand).
   *
   * <p>2.) Uenv = The polarization energy of the system without the ligand.
   *
   * <p>3.) Uligand = The polarization energy of the ligand by itself.
   *
   * <p>Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
   *
   * <p>Set the "use" array to true for all atoms for part 1. Set the "use" array to true for all
   * atoms except the ligand for part 2. Set the "use" array to true only for the ligand atoms for
   * part 3.
   *
   * <p>The "use" array can also be employed to turn off atoms for computing the electrostatic
   * energy of sub-structures.
   */
  private boolean[] use;

  private IntegerSchedule permanentSchedule;
  private double[][] cartesianMultipolePhi;
  private double[][] cartesianDipolePhi;
  private double[][] cartesianDipolePhiCR;
  private double[][] vacuumDipolePhi;
  private double[][] vacuumDipolePhiCR;
  /** AtomicDoubleArray implementation to use. */
  private AtomicDoubleArrayImpl atomicDoubleArrayImpl;
  /** Field array for each thread. [threadID][X/Y/Z][atomID] */
  private AtomicDoubleArray3D field;
  /** Chain rule field array for each thread. [threadID][X/Y/Z][atomID] */
  private AtomicDoubleArray3D fieldCR;
  /** Gradient array for each thread. [threadID][X/Y/Z][atomID] */
  private AtomicDoubleArray3D grad;
  /** Torque array for each thread. [threadID][X/Y/Z][atomID] */
  private AtomicDoubleArray3D torque;
  /** Partial derivative of the gradient with respect to Lambda. [threadID][X/Y/Z][atomID] */
  private AtomicDoubleArray3D lambdaGrad;
  /** Partial derivative of the torque with respect to Lambda. [threadID][X/Y/Z][atomID] */
  private AtomicDoubleArray3D lambdaTorque;

  /**
   * ParticleMeshEwald constructor.
   *
   * @param atoms the Atom array to do electrostatics on.
   * @param molecule the molecule number for each atom.
   * @param forceField the ForceField the defines the electrostatics parameters.
   * @param crystal The boundary conditions.
   * @param elecForm The electrostatics functional form.
   * @param neighborList The NeighborList for both van der Waals and PME.
   * @param parallelTeam A ParallelTeam that delegates parallelization.
   */
  public ParticleMeshEwaldCart(
      Atom[] atoms,
      int[] molecule,
      ForceField forceField,
      Crystal crystal,
      NeighborList neighborList,
      ELEC_FORM elecForm,
      ParallelTeam parallelTeam) {
    this.atoms = atoms;
    this.molecule = molecule;
    this.forceField = forceField;
    this.crystal = crystal;
    this.parallelTeam = parallelTeam;
    this.neighborList = neighborList;
    this.elecForm = elecForm;
    neighborLists = neighborList.getNeighborList();
    permanentSchedule = neighborList.getPairwiseSchedule();
    nAtoms = atoms.length;
    nSymm = crystal.spaceGroup.getNumberOfSymOps();
    maxThreads = parallelTeam.getThreadCount();

    electric = forceField.getDouble("ELECTRIC", Constants.DEFAULT_ELECTRIC);
    poleps = forceField.getDouble("POLAR_EPS", 1e-5);

    // If PME-specific lambda term not set, default to force field-wide lambda term.
    lambdaTerm =
        forceField.getBoolean("ELEC_LAMBDATERM", forceField.getBoolean("LAMBDATERM", false));

    CompositeConfiguration properties = forceField.getProperties();
    printInducedDipoles = properties.getBoolean("pme.printInducedDipoles", false);
    noWindowing = properties.getBoolean("pme.noWindowing", false);

    ewaldParameters = new EwaldParameters();
    scaleParameters = new ScaleParameters(forceField);
    reciprocalSpaceTerm = forceField.getBoolean("RECIPTERM", true);

    SCFPredictor scfPredictor;
    try {
      String predictor = forceField.getString("SCF_PREDICTOR", "NONE");
      predictor = predictor.replaceAll("-", "_").toUpperCase();
      scfPredictor = SCFPredictor.valueOf(predictor);
    } catch (Exception e) {
      scfPredictor = SCFPredictor.NONE;
    }
    this.scfPredictor = scfPredictor;
    scfPredictorParameters = new SCFPredictorParameters();
    if (scfPredictor != SCFPredictor.NONE) {
      scfPredictorParameters.init();
    }

    String algorithm = forceField.getString("SCF_ALGORITHM", "CG");
    try {
      algorithm = algorithm.replaceAll("-", "_").toUpperCase();
      scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
    } catch (Exception e) {
      scfAlgorithm = SCFAlgorithm.CG;
    }

    // Define how force arrays will be accumulated.
    String value = forceField.getString("ARRAY_REDUCTION", "MULTI");
    try {
      atomicDoubleArrayImpl = AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
    } catch (Exception e) {
      atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
      logger.info(
          format(
              " Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value, atomicDoubleArrayImpl));
    }
    logger.fine(format("  PME using %s arrays.", atomicDoubleArrayImpl.toString()));

    if (!scfAlgorithm.isSupported(Platform.FFX)) {
      // Can't know a-priori whether this is being constructed under an FFX or OpenMM
      // ForceFieldEnergy, so fine logging.
      logger.fine(
          format(
              " SCF algorithm %s is not supported by FFX reference implementation; falling back to CG!",
              scfAlgorithm));
      scfAlgorithm = SCFAlgorithm.CG;
    }

    pcgSolver = new PCGSolver(maxThreads, poleps, forceField, nAtoms);

    alchemicalParameters = new AlchemicalParameters(forceField, lambdaTerm);

    String polar = forceField.getString("POLARIZATION", "MUTUAL");
    if (elecForm == ELEC_FORM.FIXED_CHARGE) {
      polar = "NONE";
    }
    boolean polarizationTerm = forceField.getBoolean("POLARIZETERM", true);
    if (!polarizationTerm || polar.equalsIgnoreCase("NONE")) {
      polarization = Polarization.NONE;
    } else if (polar.equalsIgnoreCase("DIRECT")) {
      polarization = Polarization.DIRECT;
    } else {
      polarization = Polarization.MUTUAL;
    }

    String temp = forceField.getString("FFT_METHOD", "PJ");
    FFTMethod method;
    try {
      method = ReciprocalSpace.FFTMethod.valueOf(temp.toUpperCase().trim());
    } catch (Exception e) {
      method = ReciprocalSpace.FFTMethod.PJ;
    }
    boolean gpuFFT = method != FFTMethod.PJ;

    if (lambdaTerm) {
      shareddEdLambda = new SharedDouble();
      sharedd2EdLambda2 = new SharedDouble();
    } else {
      shareddEdLambda = null;
      sharedd2EdLambda2 = null;
      lambdaGrad = null;
      lambdaTorque = null;
    }

    directRegion = new DirectRegion(maxThreads);
    sorRegion = new SORRegion(maxThreads, forceField);
    optRegion = new OPTRegion(maxThreads);

    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder();
      sb.append(format("\n Electrostatics       %25s\n", getClass().getSimpleName()));
      sb.append(format("   Polarization:                       %8s\n", polarization.toString()));
      if (polarization == Polarization.MUTUAL) {
        sb.append(format("    SCF Convergence Criteria:         %8.3e\n", poleps));
        sb.append(format("    SCF Predictor:                     %8s\n", this.scfPredictor));
        sb.append(format("    SCF Algorithm:                     %8s\n", scfAlgorithm));
        if (scfAlgorithm == SCFAlgorithm.SOR) {
          sb.append(format("    SOR Parameter:                     %8.3f\n", sorRegion.getSOR()));
        } else {
          sb.append(
              format(
                  "    CG Preconditioner Cut-Off:         %8.3f\n",
                  pcgSolver.preconditionerCutoff));
          sb.append(
              format(
                  "    CG Preconditioner Ewald Coeff.:    %8.3f\n", pcgSolver.preconditionerEwald));
        }
      }
      if (ewaldParameters.aewald > 0.0) {
        sb.append("   Particle-mesh Ewald\n");
        sb.append(format("    Ewald Coefficient:                 %8.3f\n", ewaldParameters.aewald));
        sb.append(format("    Particle Cut-Off:                  %8.3f (A)", ewaldParameters.off));
      } else {
        sb.append(
            format("    Electrostatics Cut-Off:            %8.3f (A)\n", ewaldParameters.off));
      }
      logger.info(sb.toString());
    }

    // Either 1 or 2; see description below.
    int sectionThreads;
    /*
     If real and reciprocal space are done sequentially or OpenCL is used,
     then realSpaceThreads == maxThreads. Otherwise the number of
     realSpaceThreads is set to ffx.realSpaceThreads.
    */
    int realSpaceThreads;
    /*
     If real and reciprocal space are done sequentially then reciprocalThreads == maxThreads.
     If CUDA is used, reciprocalThreads == 1 Otherwise,
     reciprocalThreads = maxThreads - realSpaceThreads
    */
    int reciprocalThreads;

    if (gpuFFT) {
      sectionThreads = 2;
      reciprocalThreads = 1;
      sectionTeam = new ParallelTeam(sectionThreads);
      realSpaceTeam = parallelTeam;
      fftTeam = new ParallelTeam(reciprocalThreads);
    } else {
      boolean concurrent;
      int realThreads = 1;
      try {
        realThreads = forceField.getInteger("PME_REAL_THREADS");
        if (realThreads >= maxThreads || realThreads < 1) {
          throw new Exception("pme-real-threads must be < ffx.nt and greater than 0");
        }
        concurrent = true;
      } catch (Exception e) {
        concurrent = false;
      }
      if (concurrent) {
        sectionThreads = 2;
        realSpaceThreads = realThreads;
        reciprocalThreads = maxThreads - realThreads;
        sectionTeam = new ParallelTeam(sectionThreads);
        realSpaceTeam = new ParallelTeam(realSpaceThreads);
        fftTeam = new ParallelTeam(reciprocalThreads);
      } else {
        // If pme-real-threads is not defined, then do real and reciprocal space parts sequentially.
        sectionThreads = 1;
        sectionTeam = new ParallelTeam(sectionThreads);
        realSpaceTeam = parallelTeam;
        fftTeam = parallelTeam;
      }
    }

    realSpaceNeighborParameters = new RealSpaceNeighborParameters(maxThreads);
    initializationRegion = new InitializationRegion(maxThreads, forceField);
    expandInducedDipolesRegion = new ExpandInducedDipolesRegion(maxThreads);
    initAtomArrays();

    /*
     Note that we always pass on the unit cell crystal to ReciprocalSpace
     instance even if the real space calculations require a
     ReplicatesCrystal.
    */
    if (ewaldParameters.aewald > 0.0 && reciprocalSpaceTerm) {
      reciprocalSpace =
          new ReciprocalSpace(
              this,
              crystal.getUnitCell(),
              forceField,
              atoms,
              ewaldParameters.aewald,
              fftTeam,
              parallelTeam);
      reciprocalEnergyRegion =
          new ReciprocalEnergyRegion(maxThreads, ewaldParameters.aewald, forceField);
    } else {
      reciprocalSpace = null;
      reciprocalEnergyRegion = null;
    }
    permanentFieldRegion = new PermanentFieldRegion(realSpaceTeam, forceField, lambdaTerm);
    inducedDipoleFieldRegion = new InducedDipoleFieldRegion(realSpaceTeam, forceField, lambdaTerm);
    inducedDipoleFieldReduceRegion = new InducedDipoleFieldReduceRegion(maxThreads);
    polarizationEnergyRegion = new PolarizationEnergyRegion(maxThreads, forceField);
    realSpaceEnergyRegion = new RealSpaceEnergyRegion(maxThreads, forceField, elecForm, lambdaTerm);
    reduceRegion = new ReduceRegion(maxThreads, forceField);

    pmeTimings = new PMETimings(maxThreads);

    if (lambdaTerm) {
      logger.info(alchemicalParameters.toString());
    }

    // The GK reaction field is added to the intra-molecular field to give the self-consistent
    // reaction field.
    generalizedKirkwoodTerm = forceField.getBoolean("GKTERM", false);
    if (generalizedKirkwoodTerm || alchemicalParameters.doLigandGKElec) {
      generalizedKirkwood = new GeneralizedKirkwood(forceField, atoms, this, crystal, parallelTeam);
    } else {
      generalizedKirkwood = null;
    }
  }

  /**
   * ewaldCutoff
   *
   * @param coeff a double.
   * @param maxCutoff a double.
   * @param eps a double.
   * @return a double.
   */
  public static double ewaldCutoff(double coeff, double maxCutoff, double eps) {
    // Set the tolerance value; use of 1.0d-8 requires strict convergence of the real Space sum.
    double ratio = erfc(coeff * maxCutoff) / maxCutoff;

    if (ratio > eps) {
      return maxCutoff;
    }

    // Use a binary search to refine the coefficient.
    double xlo = 0.0;
    double xhi = maxCutoff;
    double cutoff = 0.0;
    for (int j = 0; j < 100; j++) {
      cutoff = (xlo + xhi) / 2.0;
      ratio = erfc(coeff * cutoff) / cutoff;
      if (ratio >= eps) {
        xlo = cutoff;
      } else {
        xhi = cutoff;
      }
    }
    return cutoff;
  }

  public void computeInduceDipoleField() {
    expandInducedDipoles();

    if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
      reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
    }
    field.reset(parallelTeam, 0, nAtoms - 1);
    fieldCR.reset(parallelTeam, 0, nAtoms - 1);
    inducedDipoleFieldRegion.init(
        atoms,
        crystal,
        use,
        molecule,
        ipdamp,
        thole,
        coordinates,
        realSpaceNeighborParameters,
        inducedDipole,
        inducedDipoleCR,
        reciprocalSpaceTerm,
        reciprocalSpace,
        lambdaMode,
        ewaldParameters,
        field,
        fieldCR,
        pmeTimings);
    inducedDipoleFieldRegion.executeWith(sectionTeam);
    pmeTimings.realSpaceSCFTotal = inducedDipoleFieldRegion.getRealSpaceSCFTotal();

    if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
      reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
    }

    if (generalizedKirkwoodTerm) {
      pmeTimings.gkEnergyTotal = -System.nanoTime();
      generalizedKirkwood.computeInducedGKField();
      pmeTimings.gkEnergyTotal += System.nanoTime();
      logger.fine(
          format(" Computed GK induced field %8.3f (sec)", pmeTimings.gkEnergyTotal * 1.0e-9));
    }

    inducedDipoleFieldReduceRegion.init(
        atoms,
        inducedDipole,
        inducedDipoleCR,
        generalizedKirkwoodTerm,
        generalizedKirkwood,
        ewaldParameters,
        cartesianDipolePhi,
        cartesianDipolePhiCR,
        field,
        fieldCR);
    inducedDipoleFieldReduceRegion.executeWith(parallelTeam);
  }

  @Override
  public void destroy() {
    if (fftTeam != null) {
      try {
        fftTeam.shutdown();
      } catch (Exception ex) {
        logger.warning(" Exception in shutting down fftTeam");
      }
    }
    if (sectionTeam != null) {
      try {
        sectionTeam.shutdown();
      } catch (Exception ex) {
        logger.warning(" Exception in shutting down sectionTeam");
      }
    }
    if (realSpaceTeam != null) {
      try {
        realSpaceTeam.shutdown();
      } catch (Exception ex) {
        logger.warning(" Exception in shutting down realSpaceTeam");
      }
    }
  }

  /**
   * Calculate the PME electrostatic energy.
   *
   * @param gradient If <code>true</code>, the gradient will be calculated.
   * @param print If <code>true</code>, extra logging is enabled.
   * @return return the total electrostatic energy (permanent + polarization).
   */
  @Override
  public double energy(boolean gradient, boolean print) {

    this.gradient = gradient;

    // Initialize energy variables.
    totalMultipoleEnergy = 0.0;
    permanentMultipoleEnergy = 0.0;
    permanentRealSpaceEnergy = 0.0;
    permanentSelfEnergy = 0.0;
    permanentReciprocalEnergy = 0.0;
    polarizationEnergy = 0.0;
    inducedRealSpaceEnergy = 0.0;
    inducedSelfEnergy = 0.0;
    inducedReciprocalEnergy = 0.0;
    solvationEnergy = 0.0;

    // Initialize number of interactions.
    interactions = 0;
    gkInteractions = 0;

    // Initialize timing variables.
    pmeTimings.init();
    if (reciprocalSpace != null) {
      reciprocalSpace.initTimings();
    }

    // Initialize Lambda variables.
    if (lambdaTerm) {
      shareddEdLambda.set(0.0);
      sharedd2EdLambda2.set(0.0);
    }

    alchemicalParameters.doPermanentRealSpace = true;
    alchemicalParameters.permanentScale = 1.0;
    alchemicalParameters.doPolarization = true;
    alchemicalParameters.polarizationScale = 1.0;

    // Expand coordinates and rotate multipoles into the global frame.
    initializationRegion.init(
        lambdaTerm,
        esvTerm,
        isAtomTitrating,
        lambdaScaleMultipoles,
        atoms,
        coordinates,
        crystal,
        frame,
        axisAtom,
        globalMultipole,
        dMultipoledTirationESV,
        polarizability,
        thole,
        ipdamp,
        use,
        neighborLists,
        realSpaceNeighborParameters.realSpaceLists,
        alchemicalParameters.vaporLists,
        grad,
        torque,
        lambdaGrad,
        lambdaTorque);
    initializationRegion.executeWith(parallelTeam);

    // Initialize GeneralizedKirkwood.
    if (generalizedKirkwoodTerm || alchemicalParameters.doLigandGKElec) {
      generalizedKirkwood.init();
    }

    // Total permanent + polarization energy.
    double energy;
    if (!lambdaTerm) {
      lambdaMode = LambdaMode.OFF;
      energy = computeEnergy(print);
    } else {
      // Condensed phase with all atoms.
      lambdaMode = LambdaMode.CONDENSED;
      energy = condensedEnergy();
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Solvated energy: %20.8f", energy));
      }

      // Condensed phase SCF without ligand atoms.
      lambdaMode = LambdaMode.CONDENSED_NO_LIGAND;
      double temp = energy;
      energy = condensedNoLigandSCF();
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Step 2 energy:   %20.8f", energy - temp));
      }

      // Vapor ligand electrostatics.
      if (alchemicalParameters.doLigandVaporElec) {
        lambdaMode = LambdaMode.VAPOR;
        temp = energy;
        energy = ligandElec();
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Vacuum energy:   %20.8f", energy - temp));
        }
      }
    }

    /*
     Convert torques to gradients on multipole frame defining atoms. Add
     to electrostatic gradient to the total XYZ gradient.
    */
    if (gradient || lambdaTerm) {
      reduceRegion.init(
          lambdaTerm,
          gradient,
          atoms,
          coordinates,
          frame,
          axisAtom,
          grad,
          torque,
          lambdaGrad,
          lambdaTorque);
      reduceRegion.excuteWith(parallelTeam);
    }

    // Log some timings.
    if (logger.isLoggable(Level.FINE)) {
      pmeTimings.printRealSpaceTimings();
      if (ewaldParameters.aewald > 0.0 && reciprocalSpaceTerm) {
        reciprocalSpace.printTimings();
      }
    }

    return permanentMultipoleEnergy + polarizationEnergy;
  }

  public void expandInducedDipoles() {
    if (nSymm > 1) {
      expandInducedDipolesRegion.init(atoms, crystal, inducedDipole, inducedDipoleCR);
      expandInducedDipolesRegion.executeWith(parallelTeam);
    }
  }

  @Override
  public int[][] getAxisAtoms() {
    return axisAtom;
  }

  @Override
  public double getCavitationEnergy() {
    return generalizedKirkwood.getCavitationEnergy();
  }

  @Override
  public double[][][] getCoordinates() {
    return coordinates;
  }

  @Override
  public double getDispersionEnergy() {
    return generalizedKirkwood.getDispersionEnergy();
  }

  @Override
  public ELEC_FORM getElecForm() {
    return elecForm;
  }

  @Override
  public double getEwaldCoefficient() {
    return ewaldParameters.aewald;
  }

  @Override
  public double getEwaldCutoff() {
    return ewaldParameters.off;
  }

  @Override
  public GeneralizedKirkwood getGK() {
    return generalizedKirkwood;
  }

  /**
   * getGeneralizedKirkwoodEnergy.
   *
   * @return a double.
   */
  @Override
  public double getGKEnergy() {
    return generalizedKirkwood.getGeneralizedKirkwoordEnergy();
  }

  /**
   * getGKInteractions
   *
   * @return a int.
   */
  @Override
  public int getGKInteractions() {
    return gkInteractions;
  }

  /**
   * getGradient
   *
   * @param grad an array of double.
   */
  public void getGradients(double[][] grad) {
    double[] x = grad[0];
    double[] y = grad[1];
    double[] z = grad[2];
    for (int i = 0; i < nAtoms; i++) {
      x[i] = this.grad.getX(i);
      y[i] = this.grad.getY(i);
      z[i] = this.grad.getZ(i);
    }
  }

  /**
   * Getter for the field <code>interactions</code>.
   *
   * @return a int.
   */
  @Override
  public int getInteractions() {
    return interactions;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Get the current lambda scale value.
   */
  @Override
  public double getLambda() {
    return lambda;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Set the electrostatic lambda scaling factor.
   */
  @Override
  public void setLambda(double lambda) {
    assert (lambda >= 0.0 && lambda <= 1.0);
    if (!lambdaTerm) {
      return;
    }
    this.lambda = lambda;

    initSoftCoreInit();
    alchemicalParameters.update(lambda);
    if (generalizedKirkwoodTerm) {
      generalizedKirkwood.setLambda(alchemicalParameters.polLambda);
      generalizedKirkwood.setLambdaFunction(
          alchemicalParameters.lPowPol,
          alchemicalParameters.dlPowPol,
          alchemicalParameters.d2lPowPol);
    }
  }

  @Override
  public String getName() {
    return "Cartesian";
  }

  public double getPermanentEnergy() {
    return permanentMultipoleEnergy;
  }

  public double getPermanentRealSpaceEnergy() {
    return permanentRealSpaceEnergy;
  }

  public double getPermanentReciprocalEnergy() {
    return permanentReciprocalEnergy;
  }

  @Override
  public double getPolarEps() {
    return poleps;
  }

  @Override
  public int[][] getPolarization11() {
    return ip11;
  }

  @Override
  public int[][] getPolarization12() {
    return ip12;
  }

  @Override
  public int[][] getPolarization13() {
    return ip13;
  }

  /**
   * Getter for the field <code>polarizationEnergy</code>.
   *
   * @return a double.
   */
  @Override
  public double getPolarizationEnergy() {
    return polarizationEnergy;
  }

  @Override
  public Polarization getPolarizationType() {
    return polarization;
  }

  @Override
  public ReciprocalSpace getReciprocalSpace() {
    return reciprocalSpace;
  }

  @Override
  public double getScale14() {
    return scaleParameters.m14scale;
  }

  /**
   * getGKEnergy
   *
   * @return a double.
   */
  @Override
  public double getSolvationEnergy() {
    return solvationEnergy;
  }

  /** {@inheritDoc} */
  @Override
  public double getd2EdL2() {
    if (sharedd2EdLambda2 == null || !lambdaTerm) {
      return 0.0;
    }
    double d2EdL2 = sharedd2EdLambda2.get();
    if (generalizedKirkwoodTerm || alchemicalParameters.doLigandGKElec) {
      d2EdL2 += generalizedKirkwood.getd2EdL2();
    }
    return d2EdL2;
  }

  /** {@inheritDoc} */
  @Override
  public double getdEdL() {
    if (shareddEdLambda == null || !lambdaTerm) {
      return 0.0;
    }
    double dEdL = shareddEdLambda.get();
    if (generalizedKirkwoodTerm || alchemicalParameters.doLigandGKElec) {
      dEdL += generalizedKirkwood.getdEdL();
    }
    return dEdL;
  }

  /** {@inheritDoc} */
  @Override
  public void getdEdXdL(double[] gradient) {
    if (lambdaGrad == null || !lambdaTerm) {
      return;
    }
    // Note that the Generalized Kirkwood contributions are already in the lambdaGrad array.
    int index = 0;
    for (int i = 0; i < nAtoms; i++) {
      if (atoms[i].isActive()) {
        gradient[index++] += lambdaGrad.getX(i);
        gradient[index++] += lambdaGrad.getY(i);
        gradient[index++] += lambdaGrad.getZ(i);
      }
    }
  }

  @Override
  public void setAtoms(Atom[] atoms, int[] molecule) {
    if (lambdaTerm && atoms.length != nAtoms) {
      logger.severe(" Changing the number of atoms is not compatible with use of Lambda.");
    }
    this.atoms = atoms;
    this.molecule = molecule;
    nAtoms = atoms.length;
    initAtomArrays();

    if (reciprocalSpace != null) {
      reciprocalSpace.setAtoms(atoms);
    }

    if (generalizedKirkwood != null) {
      generalizedKirkwood.setAtoms(atoms);
    }
  }

  @Override
  public void setCrystal(Crystal crystal) {
    // Check if memory allocation is required.
    int nSymmNew = crystal.spaceGroup.getNumberOfSymOps();
    if (nSymm < nSymmNew) {
      coordinates = new double[nSymmNew][3][nAtoms];
      globalMultipole = new double[nSymmNew][nAtoms][10];
      inducedDipole = new double[nSymmNew][nAtoms][3];
      inducedDipoleCR = new double[nSymmNew][nAtoms][3];
      if (generalizedKirkwood != null) {
        generalizedKirkwood.setAtoms(atoms);
      }
      realSpaceNeighborParameters.allocate(nAtoms, nSymmNew);
      pcgSolver.allocateLists(nSymmNew, nAtoms);
    }
    nSymm = nSymmNew;
    neighborLists = neighborList.getNeighborList();
    this.crystal = crystal;
    if (reciprocalSpace != null) {
      reciprocalSpace.setCrystal(crystal.getUnitCell());
    }
  }

  /**
   * Pass in atoms that have been assigned electrostatics from a fixed charge force field.
   *
   * @param atoms An array of atoms.
   */
  @Override
  public void setFixedCharges(Atom[] atoms) {
    for (Atom ai : atoms) {
      if (ai.getResolution() == Resolution.FIXEDCHARGE) {
        int index = ai.getIndex() - 1;
        polarizability[index] = 0.0;
        localMultipole[index][t000] = ai.getMultipoleType().getCharge();
        localMultipole[index][t100] = 0.0;
        localMultipole[index][t010] = 0.0;
        localMultipole[index][t001] = 0.0;
        localMultipole[index][t200] = 0.0;
        localMultipole[index][t020] = 0.0;
        localMultipole[index][t002] = 0.0;
        localMultipole[index][t110] = 0.0;
        localMultipole[index][t011] = 0.0;
        localMultipole[index][t101] = 0.0;
      }
    }
  }

  @Override
  public void setLambdaMultipoleScale(double multipoleScale) {
    lambdaScaleMultipoles = multipoleScale;
  }

  private void initAtomArrays() {
    if (localMultipole == null || localMultipole.length < nAtoms) {
      localMultipole = new double[nAtoms][10];
      frame = new MultipoleFrameDefinition[nAtoms];
      axisAtom = new int[nAtoms][];
      cartesianMultipolePhi = new double[nAtoms][tensorCount];
      directDipole = new double[nAtoms][3];
      directDipoleCR = new double[nAtoms][3];
      vacuumDirectDipole = new double[nAtoms][3];
      vacuumDirectDipoleCR = new double[nAtoms][3];
      if (optRegion != null) {
        int optOrder = optRegion.optOrder;
        optRegion.optDipole = new double[optOrder + 1][nAtoms][3];
        optRegion.optDipoleCR = new double[optOrder + 1][nAtoms][3];
      }
      cartesianDipolePhi = new double[nAtoms][tensorCount];
      cartesianDipolePhiCR = new double[nAtoms][tensorCount];
      vacuumDipolePhi = new double[nAtoms][tensorCount];
      vacuumDipolePhiCR = new double[nAtoms][tensorCount];
      mask12 = new int[nAtoms][];
      mask13 = new int[nAtoms][];
      mask14 = new int[nAtoms][];
      mask15 = new int[nAtoms][];
      ip11 = new int[nAtoms][];
      ip12 = new int[nAtoms][];
      ip13 = new int[nAtoms][];
      thole = new double[nAtoms];
      ipdamp = new double[nAtoms];
      polarizability = new double[nAtoms];

      if (scfAlgorithm == SCFAlgorithm.CG) {
        pcgSolver.allocateVectors(nAtoms);
      }
      pcgSolver.allocateLists(nSymm, nAtoms);

      if (scfPredictor != SCFPredictor.NONE) {
        int predictorOrder = scfPredictorParameters.predictorOrder;
        if (lambdaTerm) {
          scfPredictorParameters.predictorInducedDipole = new double[3][predictorOrder][nAtoms][3];
          scfPredictorParameters.predictorInducedDipoleCR =
              new double[3][predictorOrder][nAtoms][3];
        } else {
          scfPredictorParameters.predictorInducedDipole = new double[1][predictorOrder][nAtoms][3];
          scfPredictorParameters.predictorInducedDipoleCR =
              new double[1][predictorOrder][nAtoms][3];
        }
      }

      // Initialize per-thread memory for collecting the gradient, torque, field and chain-rule
      // field.
      grad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, maxThreads);
      torque = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, maxThreads);
      field = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, maxThreads);
      fieldCR = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, maxThreads);
      lambdaGrad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, maxThreads);
      lambdaTorque = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, maxThreads);
      isSoft = new boolean[nAtoms];
      use = new boolean[nAtoms];

      coordinates = new double[nSymm][3][nAtoms];
      globalMultipole = new double[nSymm][nAtoms][10];
      inducedDipole = new double[nSymm][nAtoms][3];
      inducedDipoleCR = new double[nSymm][nAtoms][3];
      vacuumInducedDipole = new double[nSymm][nAtoms][3];
      vacuumInducedDipoleCR = new double[nSymm][nAtoms][3];

      // The size of reduced neighbor list depends on the size of the real space cutoff.
      realSpaceNeighborParameters.allocate(nAtoms, nSymm);

      // Lambda factors are different for OST and ESV interactions.
      lambdaFactors = new LambdaFactors[maxThreads];
      for (int i = 0; i < maxThreads; i++) {
        if (esvTerm) {
          // Invoked every time through inner loops.
          lambdaFactors[i] = new LambdaFactorsESV();
        } else if (lambdaTerm) {
          // Invoked on calls to setLambda().
          lambdaFactors[i] = new LambdaFactorsOST();
        } else {
          // Invoked never; inoperative defaults.
          lambdaFactors[i] = LambdaDefaults;
        }
      }
    }

    // Initialize the soft core lambda mask to false for all atoms.
    fill(isSoft, false);

    // Initialize the use mask to true for all atoms.
    fill(use, true);

    // Assign multipole parameters.
    assignMultipoles();

    // Assign polarization groups.
    if (elecForm == PAM) {
      assignPolarizationGroups();
    }

    // Fill the thole, inverse polarization damping and polarizability arrays.
    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];

      if (elecForm == PAM) {
        PolarizeType polarizeType = ai.getPolarizeType();
        int index = ai.getIndex() - 1;
        thole[index] = polarizeType.thole;
        ipdamp[index] = polarizeType.pdamp;
        if (!(ipdamp[index] > 0.0)) {
          ipdamp[index] = Double.POSITIVE_INFINITY;
        } else {
          ipdamp[index] = 1.0 / ipdamp[index];
        }
        polarizability[index] = polarizeType.polarizability;
      }

      // Collect 1-2 interactions.
      List<Atom> n12 = ai.get12List();
      mask12[i] = new int[n12.size()];
      int j = 0;
      for (Atom a12 : n12) {
        mask12[i][j++] = a12.getIndex() - 1;
      }

      // Collect 1-3 interactions.
      List<Atom> n13 = ai.get13List();
      mask13[i] = new int[n13.size()];
      j = 0;
      for (Atom a13 : n13) {
        mask13[i][j++] = a13.getIndex() - 1;
      }

      // Collect 1-4 interactions.
      List<Atom> n14 = ai.get14List();
      mask14[i] = new int[n14.size()];
      j = 0;
      for (Atom a14 : n14) {
        mask14[i][j++] = a14.getIndex() - 1;
      }

      // Collect 1-5 interactions.
      List<Atom> n15 = ai.get15List();
      mask15[i] = new int[n15.size()];
      j = 0;
      for (Atom a15 : n15) {
        mask15[i][j++] = a15.getIndex() - 1;
      }
    }
  }

  /** Initialize a boolean array of soft atoms and, if requested, ligand vapor electrostatics. */
  private void initSoftCoreInit() {

    boolean rebuild = false;

    // Initialize a boolean array that marks soft atoms.
    StringBuilder sb = new StringBuilder("\n Softcore Atoms:\n");
    int count = 0;
    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];
      boolean soft = ai.applyLambda();
      if (soft != isSoft[i]) {
        rebuild = true;
      }
      isSoft[i] = soft;
      if (soft) {
        count++;
        sb.append(ai).append("\n");
      }
    }

    // Force rebuild ligand vapor electrostatics are being computed and vaporCrystal is null.
    if (alchemicalParameters.doLigandVaporElec && alchemicalParameters.vaporCrystal == null) {
      rebuild = true;
    }

    if (!rebuild) {
      return;
    }

    if (count > 0 && logger.isLoggable(Level.FINE)) {
      logger.fine(format(" Softcore atom count: %d", count));
      logger.fine(sb.toString());
    }

    // Initialize boundary conditions,
    // an n^2 neighbor list and parallel scheduling for ligand vapor electrostatics.
    if (alchemicalParameters.doLigandVaporElec) {
      double maxr = 10.0;
      for (int i = 0; i < nAtoms; i++) {
        Atom ai = atoms[i];
        if (ai.applyLambda()) {

          // Determine ligand size.
          for (int j = i + 1; j < nAtoms; j++) {
            Atom aj = atoms[j];
            if (aj.applyLambda()) {
              double dx = ai.getX() - aj.getX();
              double dy = ai.getY() - aj.getY();
              double dz = ai.getZ() - aj.getZ();
              double r = sqrt(dx * dx + dy * dy + dz * dz);
              maxr = max(r, maxr);
            }
          }
        }
      }

      double vacuumOff = 2 * maxr;
      alchemicalParameters.vaporCrystal =
          new Crystal(3 * vacuumOff, 3 * vacuumOff, 3 * vacuumOff, 90.0, 90.0, 90.0, "P1");
      alchemicalParameters.vaporCrystal.setAperiodic(true);
      NeighborList vacuumNeighborList =
          new NeighborList(
              null, alchemicalParameters.vaporCrystal, atoms, vacuumOff, 2.0, parallelTeam);
      vacuumNeighborList.setIntermolecular(false, molecule);

      alchemicalParameters.vaporLists = new int[1][nAtoms][];
      double[][] coords = new double[1][nAtoms * 3];
      for (int i = 0; i < nAtoms; i++) {
        coords[0][i * 3] = atoms[i].getX();
        coords[0][i * 3 + 1] = atoms[i].getY();
        coords[0][i * 3 + 2] = atoms[i].getZ();
      }
      boolean print = logger.isLoggable(Level.FINE);
      vacuumNeighborList.buildList(coords, alchemicalParameters.vaporLists, isSoft, true, print);
      alchemicalParameters.vaporPermanentSchedule = vacuumNeighborList.getPairwiseSchedule();
      alchemicalParameters.vaporEwaldSchedule = alchemicalParameters.vaporPermanentSchedule;
      alchemicalParameters.vacuumRanges = new Range[maxThreads];
      vacuumNeighborList.setDisableUpdates(
          forceField.getBoolean("DISABLE_NEIGHBOR_UPDATES", false));
    } else {
      alchemicalParameters.vaporCrystal = null;
      alchemicalParameters.vaporLists = null;
      alchemicalParameters.vaporPermanentSchedule = null;
      alchemicalParameters.vaporEwaldSchedule = null;
      alchemicalParameters.vacuumRanges = null;
    }
  }

  /**
   * 1.) Total system under PBC. A.) Softcore real space for Ligand-Protein and Ligand-Ligand. B.)
   * Reciprocal space scaled by lambda. C.) Polarization scaled by lambda.
   */
  private double condensedEnergy() {
    if (lambda < alchemicalParameters.polLambdaStart) {
      /*
       If the polarization has been completely decoupled, the
       contribution of the complete system is zero.
       We can skip the SCF for part 1 for efficiency.
      */
      alchemicalParameters.polarizationScale = 0.0;
      alchemicalParameters.doPolarization = false;
    } else if (lambda <= alchemicalParameters.polLambdaEnd) {
      alchemicalParameters.polarizationScale = alchemicalParameters.lPowPol;
      alchemicalParameters.doPolarization = true;
    } else {
      alchemicalParameters.polarizationScale = 1.0;
      alchemicalParameters.doPolarization = true;
    }
    alchemicalParameters.doPermanentRealSpace = true;
    alchemicalParameters.permanentScale = alchemicalParameters.lPowPerm;
    alchemicalParameters.dEdLSign = 1.0;

    return computeEnergy(false);
  }

  /**
   * 2.) Condensed phase system without the ligand. A.) No permanent real space electrostatics needs
   * to be calculated because this was handled analytically in step 1.
   *
   * <p>B.) Permanent reciprocal space scaled by (1 - lambda).
   *
   * <p>C.) Polarization scaled by (1 - lambda).
   */
  private double condensedNoLigandSCF() {
    // Turn off the ligand.
    boolean skip = true;
    for (int i = 0; i < nAtoms; i++) {
      if (atoms[i].applyLambda()) {
        use[i] = false;
      } else {
        use[i] = true;
        skip = false;
      }
    }

    // Permanent real space is done for the condensed phase. Scale the reciprocal space part.
    alchemicalParameters.doPermanentRealSpace = false;
    alchemicalParameters.permanentScale = 1.0 - alchemicalParameters.lPowPerm;
    alchemicalParameters.dEdLSign = -1.0;

    // If we are past the end of the polarization lambda window, then only the condensed phase is
    // necessary.
    if (lambda <= alchemicalParameters.polLambdaEnd
        && alchemicalParameters.doNoLigandCondensedSCF) {
      alchemicalParameters.doPolarization = true;
      alchemicalParameters.polarizationScale = 1.0 - alchemicalParameters.lPowPol;
    } else {
      alchemicalParameters.doPolarization = false;
      alchemicalParameters.polarizationScale = 0.0;
    }

    // Turn off GK.
    boolean gkBack = generalizedKirkwoodTerm;
    generalizedKirkwoodTerm = false;

    // If we are disappearing the entire system (ie. a small crystal) then
    // the energy of this step is 0 and we can skip it.
    double energy;
    if (skip) {
      energy = permanentMultipoleEnergy + polarizationEnergy + solvationEnergy;
    } else {
      energy = computeEnergy(false);
      for (int i = 0; i < nAtoms; i++) {
        use[i] = true;
      }
    }

    generalizedKirkwoodTerm = gkBack;

    return energy;
  }

  /**
   * 3.) Aperiodic ligand electrostatics.
   *
   * <p>A.) Real space with an Ewald coefficient of 0.0 (no reciprocal space).
   *
   * <p>B.) Polarization scaled as in Step 2 by (1 - lambda).
   */
  private double ligandElec() {
    for (int i = 0; i < nAtoms; i++) {
      use[i] = atoms[i].applyLambda();
    }

    // Scale the permanent vacuum electrostatics. The softcore alpha is not
    // necessary (nothing in vacuum to collide with).
    alchemicalParameters.doPermanentRealSpace = true;
    alchemicalParameters.permanentScale = 1.0 - alchemicalParameters.lPowPerm;
    alchemicalParameters.dEdLSign = -1.0;
    double lAlphaBack = alchemicalParameters.lAlpha;
    double dlAlphaBack = alchemicalParameters.dlAlpha;
    double d2lAlphaBack = alchemicalParameters.d2lAlpha;
    alchemicalParameters.lAlpha = 0.0;
    alchemicalParameters.dlAlpha = 0.0;
    alchemicalParameters.d2lAlpha = 0.0;

    // If we are past the end of the polarization lambda window, then only
    // the condensed phase is necessary.
    if (lambda <= alchemicalParameters.polLambdaEnd) {
      alchemicalParameters.doPolarization = true;
      alchemicalParameters.polarizationScale = 1.0 - alchemicalParameters.lPowPol;
    } else {
      alchemicalParameters.doPolarization = false;
      alchemicalParameters.polarizationScale = 0.0;
    }

    // Save the current real space PME parameters.
    double offBack = ewaldParameters.off;
    double aewaldBack = ewaldParameters.aewald;
    ewaldParameters.setEwaldParameters(Double.MAX_VALUE, 0.0);
    // Save the current parallelization schedule.
    IntegerSchedule permanentScheduleBack = permanentSchedule;
    IntegerSchedule ewaldScheduleBack = realSpaceNeighborParameters.realSpaceSchedule;
    Range[] rangesBack = realSpaceNeighborParameters.realSpaceRanges;
    permanentSchedule = alchemicalParameters.vaporPermanentSchedule;
    realSpaceNeighborParameters.realSpaceSchedule = alchemicalParameters.vaporEwaldSchedule;
    realSpaceNeighborParameters.realSpaceRanges = alchemicalParameters.vacuumRanges;

    // Use vacuum crystal / vacuum neighborLists.
    Crystal crystalBack = crystal;
    int nSymmBack = nSymm;
    int[][][] listsBack = neighborLists;
    neighborLists = alchemicalParameters.vaporLists;
    crystal = alchemicalParameters.vaporCrystal;
    nSymm = 1;

    // Turn off GK if in use.
    boolean gkBack = generalizedKirkwoodTerm;

    // Turn off Pre-conditioned conjugate gradient SCF solver.
    SCFAlgorithm scfBack = scfAlgorithm;
    scfAlgorithm = SCFAlgorithm.SOR;

    if (alchemicalParameters.doLigandGKElec) {
      generalizedKirkwoodTerm = true;
      generalizedKirkwood.setNeighborList(alchemicalParameters.vaporLists);
      generalizedKirkwood.setLambda(lambda);
      generalizedKirkwood.setCutoff(ewaldParameters.off);
      generalizedKirkwood.setCrystal(alchemicalParameters.vaporCrystal);
      generalizedKirkwood.setLambdaFunction(
          alchemicalParameters.polarizationScale,
          alchemicalParameters.dEdLSign * alchemicalParameters.dlPowPol,
          alchemicalParameters.dEdLSign * alchemicalParameters.d2lPowPol);
    } else {
      generalizedKirkwoodTerm = false;
    }

    double energy = computeEnergy(false);

    // Revert to the saved parameters.
    ewaldParameters.setEwaldParameters(offBack, aewaldBack);
    neighborLists = listsBack;
    crystal = crystalBack;
    nSymm = nSymmBack;
    permanentSchedule = permanentScheduleBack;
    realSpaceNeighborParameters.realSpaceSchedule = ewaldScheduleBack;
    realSpaceNeighborParameters.realSpaceRanges = rangesBack;
    alchemicalParameters.lAlpha = lAlphaBack;
    alchemicalParameters.dlAlpha = dlAlphaBack;
    alchemicalParameters.d2lAlpha = d2lAlphaBack;
    generalizedKirkwoodTerm = gkBack;
    scfAlgorithm = scfBack;

    fill(use, true);

    return energy;
  }

  /**
   * Calculate the PME electrostatic energy for a Lambda state.
   *
   * @param print If <code>true</code>, extra logging is enabled.
   * @return return the total electrostatic energy (permanent + polarization).
   */
  private double computeEnergy(boolean print) {
    // Find the permanent multipole potential, field, etc.
    permanentMultipoleField();

    // Compute Born radii if necessary.
    if (generalizedKirkwoodTerm) {
      pmeTimings.bornRadiiTotal -= System.nanoTime();
      generalizedKirkwood.setUse(use);
      generalizedKirkwood.computeBornRadii();
      pmeTimings.bornRadiiTotal += System.nanoTime();
    }

    // Do the self-consistent field calculation.
    if (polarization != Polarization.NONE && alchemicalParameters.doPolarization) {
      // Compute vacuum dipole moments.
      if (generalizedKirkwoodTerm) {
        generalizedKirkwoodTerm = false;
        // Run the vacuum SCF.
        selfConsistentField(logger.isLoggable(Level.FINE));
        // Store vacuum dipole moments
        saveInducedDipolesToVacuumDipoles();

        generalizedKirkwoodTerm = true;
        // Load the permanent multipole field.
        permanentMultipoleField();
      }

      // Compute induced dipoles.
      selfConsistentField(logger.isLoggable(Level.FINE));

      if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
        if (gradient && polarization == Polarization.DIRECT) {
          reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
          field.reset(parallelTeam, 0, nAtoms - 1);
          fieldCR.reset(parallelTeam, 0, nAtoms - 1);
          inducedDipoleFieldRegion.init(
              atoms,
              crystal,
              use,
              molecule,
              ipdamp,
              thole,
              coordinates,
              realSpaceNeighborParameters,
              inducedDipole,
              inducedDipoleCR,
              reciprocalSpaceTerm,
              reciprocalSpace,
              lambdaMode,
              ewaldParameters,
              field,
              fieldCR,
              pmeTimings);
          inducedDipoleFieldRegion.executeWith(sectionTeam);
          reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
        } else {
          reciprocalSpace.cartToFracInducedDipoles(inducedDipole, inducedDipoleCR);
        }
      }

      if (scfPredictor != SCFPredictor.NONE) {
        scfPredictorParameters.saveMutualInducedDipoles();
      }

      if (printInducedDipoles) {
        StringBuilder sb = new StringBuilder();
        sb.append("     Atom                                         Induced Dipole \n");
        sb.append("    ======                                       ================\n");
        for (int i = 0; i < nAtoms; i++) {
          sb.append(
              format(
                  "%-47s: (%+8.6f %+8.6f %+8.6f)\n",
                  atoms[i],
                  inducedDipole[0][i][0],
                  inducedDipole[0][i][1],
                  inducedDipole[0][i][2]));
        }
        logger.info(sb.toString());
      }
    }

    /*
     Find the total real space energy. This includes:
      1) the permanent multipoles in their own real space potential
      2) the permenant multipoles in their own reciprocal space potential
      3) the permanent multipoles interacting with the induced dipole real space potential
      4) the permanent multipoles interacting with the induced dipole reciprocal space potential
    */

    double[][][] inputDipole = inducedDipole;
    double[][][] inputDipoleCR = inducedDipoleCR;
    double[][] inputPhi = cartesianDipolePhi;
    double[][] inputPhiCR = cartesianDipolePhiCR;

    if (generalizedKirkwoodTerm) {
      // With GK, we need to compute the polarization energy with vacuum induced dipoles.
      inputDipole = vacuumInducedDipole;
      inputDipoleCR = vacuumInducedDipoleCR;
      inputPhi = vacuumDipolePhi;
      inputPhiCR = vacuumDipolePhiCR;
    }

    double eself = 0.0;
    double erecip = 0.0;
    double eselfi = 0.0;
    double erecipi = 0.0;
    polarizationEnergyRegion.setPolarizationEnergy(0.0);
    if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
      reciprocalEnergyRegion.init(
          atoms,
          crystal,
          use,
          globalMultipole,
          cartesianMultipolePhi,
          inputDipole,
          inputDipoleCR,
          inputPhi,
          inputPhiCR,
          reciprocalSpace,
          polarization,
          grad,
          torque,
          lambdaGrad,
          lambdaTorque,
          gradient,
          lambdaTerm,
          shareddEdLambda,
          sharedd2EdLambda2,
          alchemicalParameters);
      reciprocalEnergyRegion.executeWith(parallelTeam);
      eself = reciprocalEnergyRegion.getPermanentSelfEnergy();
      erecip = reciprocalEnergyRegion.getPermanentReciprocalEnergy();
      eselfi = reciprocalEnergyRegion.getInducedDipoleSelfEnergy();
      erecipi = reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
      interactions += nAtoms;
    }

    pmeTimings.realSpaceEnergyTotal -= System.nanoTime();
    realSpaceEnergyRegion.init(
        atoms,
        crystal,
        coordinates,
        frame,
        axisAtom,
        globalMultipole,
        inputDipole,
        inputDipoleCR,
        use,
        molecule,
        ip11,
        mask12,
        mask13,
        mask14,
        mask15,
        isSoft,
        ipdamp,
        thole,
        realSpaceNeighborParameters,
        gradient,
        lambdaTerm,
        lambdaMode,
        polarization,
        ewaldParameters,
        scaleParameters,
        alchemicalParameters,
        pmeTimings.realSpaceEnergyTime,
        // Output
        grad,
        torque,
        lambdaGrad,
        lambdaTorque,
        shareddEdLambda,
        sharedd2EdLambda2);
    realSpaceEnergyRegion.executeWith(parallelTeam);
    double ereal = realSpaceEnergyRegion.getPermanentEnergy();
    double ereali = realSpaceEnergyRegion.getPolarizationEnergy();
    interactions += realSpaceEnergyRegion.getInteractions();
    pmeTimings.realSpaceEnergyTotal += System.nanoTime();

    if (generalizedKirkwoodTerm) {

      // Compute the polarization energy cost to polarize the induced dipoles
      // from vacuum to the SCRF values.
      double eGK = 0.0;

      // Create default alchemical parameters.
      AlchemicalParameters alchemicalParametersGK = new AlchemicalParameters(forceField, false);
      // No permanent multipole contribution.
      alchemicalParametersGK.permanentScale = 0.0;
      alchemicalParametersGK.doPermanentRealSpace = false;
      // Flip the sign on the vacuum polarization energy and derivatives.
      alchemicalParametersGK.polarizationScale = -1.0;

      // Store the derivative and torque contributions in the GK arrays.
      AtomicDoubleArray3D gradGK = generalizedKirkwood.getGrad();
      AtomicDoubleArray3D torqueGK = generalizedKirkwood.getTorque();

      if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
        logger.severe(" Compute reciprocal space contribution to GK.");
        reciprocalEnergyRegion.init(
            atoms,
            crystal,
            use,
            globalMultipole,
            cartesianMultipolePhi,
            vacuumInducedDipole,
            vacuumInducedDipoleCR,
            vacuumDipolePhi,
            vacuumDipolePhiCR,
            reciprocalSpace,
            polarization,
            gradGK,
            torqueGK,
            null,
            null,
            gradient,
            false,
            shareddEdLambda,
            sharedd2EdLambda2,
            alchemicalParametersGK);
        reciprocalEnergyRegion.executeWith(parallelTeam);
        eGK =
            reciprocalEnergyRegion.getInducedDipoleSelfEnergy()
                + reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
      }

      pmeTimings.gkEnergyTotal -= System.nanoTime();
      realSpaceEnergyRegion.init(
          atoms,
          crystal,
          coordinates,
          frame,
          axisAtom,
          globalMultipole,
          vacuumInducedDipole,
          vacuumInducedDipoleCR,
          use,
          molecule,
          ip11,
          mask12,
          mask13,
          mask14,
          mask15,
          isSoft,
          ipdamp,
          thole,
          realSpaceNeighborParameters,
          gradient,
          false,
          lambdaMode,
          polarization,
          ewaldParameters,
          scaleParameters,
          alchemicalParametersGK,
          pmeTimings.realSpaceEnergyTime,
          // Output
          gradGK,
          torqueGK,
          null,
          null,
          shareddEdLambda,
          sharedd2EdLambda2);
      realSpaceEnergyRegion.executeWith(parallelTeam);
      eGK += realSpaceEnergyRegion.getPolarizationEnergy();

      // Normal sign of +1 for the SCRF polarization energy and derivatives.
      alchemicalParametersGK.polarizationScale = 1.0;

      if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
        logger.severe(" Compute reciprocal space contribution to GK.");
        reciprocalEnergyRegion.init(
            atoms,
            crystal,
            use,
            globalMultipole,
            cartesianMultipolePhi,
            inducedDipole,
            inducedDipoleCR,
            cartesianDipolePhi,
            cartesianDipolePhiCR,
            reciprocalSpace,
            polarization,
            gradGK,
            torqueGK,
            null,
            null,
            gradient,
            false,
            shareddEdLambda,
            sharedd2EdLambda2,
            alchemicalParametersGK);
        reciprocalEnergyRegion.executeWith(parallelTeam);
        eGK +=
            reciprocalEnergyRegion.getInducedDipoleSelfEnergy()
                + reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
      }

      pmeTimings.gkEnergyTotal -= System.nanoTime();
      realSpaceEnergyRegion.init(
          atoms,
          crystal,
          coordinates,
          frame,
          axisAtom,
          globalMultipole,
          inducedDipole,
          inducedDipoleCR,
          use,
          molecule,
          ip11,
          mask12,
          mask13,
          mask14,
          mask15,
          isSoft,
          ipdamp,
          thole,
          realSpaceNeighborParameters,
          gradient,
          false,
          lambdaMode,
          polarization,
          ewaldParameters,
          scaleParameters,
          alchemicalParametersGK,
          pmeTimings.realSpaceEnergyTime,
          // Output
          gradGK,
          torqueGK,
          null,
          null,
          shareddEdLambda,
          sharedd2EdLambda2);
      realSpaceEnergyRegion.executeWith(parallelTeam);
      eGK += realSpaceEnergyRegion.getPolarizationEnergy();

      // Compute the solvation free energy.
      solvationEnergy += generalizedKirkwood.solvationEnergy(eGK, gradient, print);
      if (gradient) {
        // Add the GK derivative contributions into the overall derivatives.
        generalizedKirkwood.reduce(grad, torque, lambdaGrad, lambdaTorque);
      }
      gkInteractions += generalizedKirkwood.getInteractions();
      pmeTimings.gkEnergyTotal += System.nanoTime();
    }

    // Collect energy terms.
    permanentRealSpaceEnergy += ereal;
    permanentSelfEnergy += eself;
    permanentReciprocalEnergy += erecip;
    inducedRealSpaceEnergy += ereali;
    inducedSelfEnergy += eselfi;
    inducedReciprocalEnergy += erecipi;
    permanentMultipoleEnergy += eself + erecip + ereal;
    polarizationEnergy += eselfi + erecipi + ereali;
    totalMultipoleEnergy += ereal + eself + erecip + ereali + eselfi + erecipi;

    // Total from single loop polarization energy.
    // polarizationEnergy += polarizationEnergyRegion.getPolarizationEnergy();
    // totalMultipoleEnergy += ereal + eself + erecip +
    // polarizationEnergyRegion.getPolarizationEnergy();

    // Log some info.
    if (logger.isLoggable(Level.FINE) || printDecomposition) {
      StringBuilder sb = new StringBuilder();
      sb.append(format("\n Global Cartesian PME, lambdaMode=%s\n", lambdaMode.toString()));
      sb.append(format(" Multipole Self-Energy:   %16.8f\n", eself));
      sb.append(format(" Multipole Reciprocal:    %16.8f\n", erecip));
      sb.append(format(" Multipole Real Space:    %16.8f\n", ereal));
      sb.append(format(" Polarization Self-Energy:%16.8f\n", eselfi));
      sb.append(format(" Polarization Reciprocal: %16.8f\n", erecipi));
      sb.append(format(" Polarization Real Space: %16.8f\n", ereali));
      if (generalizedKirkwoodTerm) {
        sb.append(format(" Generalized Kirkwood:    %16.8f\n", solvationEnergy));
      }
      logger.info(sb.toString());
    }

    return permanentMultipoleEnergy + polarizationEnergy + solvationEnergy;
  }

  /** Find the permanent multipole potential, field, etc. */
  private void permanentMultipoleField() {
    try {
      // Compute b-Splines and permanent density.
      if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
        reciprocalSpace.computeBSplines();
        reciprocalSpace.splinePermanentMultipoles(globalMultipole, 0, use);
      }

      field.reset(parallelTeam, 0, nAtoms - 1);
      fieldCR.reset(parallelTeam, 0, nAtoms - 1);
      permanentFieldRegion.init(
          atoms,
          crystal,
          coordinates,
          globalMultipole,
          inducedDipole,
          inducedDipoleCR,
          neighborLists,
          scaleParameters,
          use,
          molecule,
          ipdamp,
          thole,
          ip11,
          mask12,
          mask13,
          mask14,
          lambdaMode,
          reciprocalSpaceTerm,
          reciprocalSpace,
          ewaldParameters,
          pcgSolver,
          permanentSchedule,
          realSpaceNeighborParameters,
          field,
          fieldCR,
          pmeTimings);
      // The real space contribution can be calculated at the same time
      // the reciprocal space convolution is being done.
      sectionTeam.execute(permanentFieldRegion);

      pmeTimings.realSpacePermTotal = permanentFieldRegion.getRealSpacePermTotal();

      // Collect the reciprocal space field.
      if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
        reciprocalSpace.computePermanentPhi(cartesianMultipolePhi);
      }
    } catch (RuntimeException e) {
      String message = "Fatal exception computing the permanent multipole field.\n";
      logger.log(Level.WARNING, message, e);
      throw e;
    } catch (Exception e) {
      String message = "Fatal exception computing the permanent multipole field.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  private void saveInducedDipolesToVacuumDipoles() {
    for (int i = 0; i < nAtoms; i++) {
      System.arraycopy(directDipole[i], 0, vacuumDirectDipole[i], 0, 3);
      System.arraycopy(directDipoleCR[i], 0, vacuumDirectDipoleCR[i], 0, 3);
      System.arraycopy(inducedDipole[0][i], 0, vacuumInducedDipole[0][i], 0, 3);
      System.arraycopy(inducedDipoleCR[0][i], 0, vacuumInducedDipoleCR[0][i], 0, 3);
      System.arraycopy(cartesianDipolePhi[i], 0, vacuumDipolePhi[i], 0, tensorCount);
      System.arraycopy(cartesianDipolePhiCR[i], 0, vacuumDipolePhiCR[i], 0, tensorCount);
      if (nSymm > 1) {
        for (int s = 1; s < nSymm; s++) {
          System.arraycopy(inducedDipole[s][i], 0, vacuumInducedDipole[s][i], 0, 3);
          System.arraycopy(inducedDipoleCR[s][i], 0, vacuumInducedDipoleCR[s][i], 0, 3);
        }
      }
    }
  }

  /** Apply the selected polarization model (NONE, Direct or Mutual). */
  private int selfConsistentField(boolean print) {
    if (polarization == Polarization.NONE) {
      return -1;
    }
    long startTime = System.nanoTime();

    // Compute the direct induced dipoles.
    if (generalizedKirkwoodTerm) {
      pmeTimings.gkEnergyTotal = -System.nanoTime();
      generalizedKirkwood.computePermanentGKField();
      pmeTimings.gkEnergyTotal += System.nanoTime();
      logger.fine(
          format(" Computed GK permanent field %8.3f (sec)", pmeTimings.gkEnergyTotal * 1.0e-9));
    }
    directRegion.init(
        atoms,
        polarizability,
        globalMultipole,
        cartesianMultipolePhi,
        field,
        fieldCR,
        generalizedKirkwoodTerm,
        generalizedKirkwood,
        ewaldParameters,
        inducedDipole,
        inducedDipoleCR,
        directDipole,
        directDipoleCR);
    directRegion.executeWith(parallelTeam);

    // Return unless mutual polarization is selected.
    if (polarization != Polarization.MUTUAL) {
      expandInducedDipoles();
      return 0;
    }

    // Predict the current self-consistent induced dipoles using information from previous steps.
    if (scfPredictor != SCFPredictor.NONE) {
      switch (scfPredictor) {
        case ASPC:
          scfPredictorParameters.aspcPredictor();
          break;
        case LS:
          scfPredictorParameters.leastSquaresPredictor();
          break;
        case POLY:
          scfPredictorParameters.polynomialPredictor();
          break;
        default:
          break;
      }
    }

    // Expand the initial induced dipoles to P1 symmetry, if necessary.
    expandInducedDipoles();

    // Converge the self-consistent field.
    int iterations;
    switch (scfAlgorithm) {
      case SOR:
        iterations = scfBySOR(print, startTime);
        break;
      case EPT:
        iterations = scfByEPT(print, startTime);
        break;
      case CG:
      default:
        pcgSolver.init(
            atoms,
            coordinates,
            polarizability,
            ipdamp,
            thole,
            use,
            crystal,
            inducedDipole,
            inducedDipoleCR,
            directDipole,
            directDipoleCR,
            field,
            fieldCR,
            ewaldParameters,
            parallelTeam,
            realSpaceNeighborParameters.realSpaceSchedule,
            pmeTimings.realSpaceSCFTime);
        iterations = pcgSolver.scfByPCG(print, startTime, this);
        break;
    }

    return iterations;
  }

  /** Converge the SCF using Successive Over-Relaxation (SOR). */
  private int scfBySOR(boolean print, long startTime) {
    long directTime = System.nanoTime() - startTime;

    // A request of 0 SCF cycles simplifies mutual polarization to direct polarization.
    StringBuilder sb = null;
    if (print) {
      sb = new StringBuilder("\n Self-Consistent Field\n" + " Iter  RMS Change (Debye)  Time\n");
    }
    int completedSCFCycles = 0;
    int maxSCFCycles = 1000;
    double eps = 100.0;
    double previousEps;
    boolean done = false;
    while (!done) {
      long cycleTime = -System.nanoTime();
      try {
        if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
          reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
        }
        field.reset(parallelTeam, 0, nAtoms - 1);
        fieldCR.reset(parallelTeam, 0, nAtoms - 1);
        inducedDipoleFieldRegion.init(
            atoms,
            crystal,
            use,
            molecule,
            ipdamp,
            thole,
            coordinates,
            realSpaceNeighborParameters,
            inducedDipole,
            inducedDipoleCR,
            reciprocalSpaceTerm,
            reciprocalSpace,
            lambdaMode,
            ewaldParameters,
            field,
            fieldCR,
            pmeTimings);
        inducedDipoleFieldRegion.executeWith(sectionTeam);
        pmeTimings.realSpaceSCFTotal = inducedDipoleFieldRegion.getRealSpaceSCFTotal();
        if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
          reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
        }

        if (generalizedKirkwoodTerm) {
          // GK field.
          pmeTimings.gkEnergyTotal = -System.nanoTime();
          generalizedKirkwood.computeInducedGKField();
          pmeTimings.gkEnergyTotal += System.nanoTime();
          logger.fine(
              format(" Computed GK induced field %8.3f (sec)", pmeTimings.gkEnergyTotal * 1.0e-9));
        }

        sorRegion.init(
            atoms,
            polarizability,
            inducedDipole,
            inducedDipoleCR,
            directDipole,
            directDipoleCR,
            cartesianDipolePhi,
            cartesianDipolePhiCR,
            field,
            fieldCR,
            generalizedKirkwoodTerm,
            generalizedKirkwood,
            ewaldParameters);
        parallelTeam.execute(sorRegion);

        expandInducedDipoles();
      } catch (Exception e) {
        String message = "Exception computing mutual induced dipoles.";
        logger.log(Level.SEVERE, message, e);
      }
      completedSCFCycles++;
      previousEps = eps;
      eps = sorRegion.getEps();
      eps = Constants.ELEC_ANG_TO_DEBYE * sqrt(eps / (double) nAtoms);
      cycleTime += System.nanoTime();
      if (print) {
        sb.append(format(" %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * NS2SEC));
      }

      // If the RMS Debye change increases, fail the SCF process.
      if (eps > previousEps) {
        if (sb != null) {
          logger.warning(sb.toString());
        }
        String message = format(" SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
        throw new EnergyException(message, false);
      }

      // The SCF should converge well before the max iteration check. Otherwise, fail the SCF
      // process.
      if (completedSCFCycles >= maxSCFCycles) {
        if (sb != null) {
          logger.warning(sb.toString());
        }
        String message = format(" Maximum SCF iterations reached: (%d)\n", completedSCFCycles);
        throw new EnergyException(message, false);
      }

      // Check if the convergence criteria has been achieved.
      if (eps < poleps) {
        done = true;
      }
    }
    if (print) {
      sb.append(format(" Direct:                  %7.4f\n", NS2SEC * directTime));
      startTime = System.nanoTime() - startTime;
      sb.append(format(" Total:                   %7.4f", startTime * NS2SEC));
      logger.info(sb.toString());
    }
    return completedSCFCycles;
  }

  /** Set the induced dipoles using extrapolated perturbation theory. */
  private int scfByEPT(boolean print, long startTime) {

    // Zeroth order OPT dipoles are the direct dipoles.
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        optRegion.optDipole[0][i][j] = directDipole[i][j];
        optRegion.optDipoleCR[0][i][j] = directDipoleCR[i][j];
      }
    }

    int optOrder = optRegion.optOrder;
    // Collect OPT dipole contributions from 1st to Nth.
    for (int currentOptOrder = 1; currentOptOrder <= optOrder; currentOptOrder++) {
      try {
        if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
          reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
        }
        field.reset(parallelTeam, 0, nAtoms - 1);
        fieldCR.reset(parallelTeam, 0, nAtoms - 1);
        inducedDipoleFieldRegion.init(
            atoms,
            crystal,
            use,
            molecule,
            ipdamp,
            thole,
            coordinates,
            realSpaceNeighborParameters,
            inducedDipole,
            inducedDipoleCR,
            reciprocalSpaceTerm,
            reciprocalSpace,
            lambdaMode,
            ewaldParameters,
            field,
            fieldCR,
            pmeTimings);
        inducedDipoleFieldRegion.executeWith(sectionTeam);
        pmeTimings.realSpaceSCFTotal = inducedDipoleFieldRegion.getRealSpaceSCFTotal();

        if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
          reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
        }

        if (generalizedKirkwoodTerm) {
          // GK field.
          pmeTimings.gkEnergyTotal = -System.nanoTime();
          generalizedKirkwood.computeInducedGKField();
          pmeTimings.gkEnergyTotal += System.nanoTime();
          logger.fine(
              format(" Computed GK induced field %8.3f (sec)", pmeTimings.gkEnergyTotal * 1.0e-9));
        }

        optRegion.init(
            currentOptOrder,
            atoms,
            polarizability,
            inducedDipole,
            inducedDipoleCR,
            cartesianDipolePhi,
            cartesianDipolePhiCR,
            field,
            fieldCR,
            generalizedKirkwoodTerm,
            generalizedKirkwood,
            ewaldParameters);
        parallelTeam.execute(optRegion);

        expandInducedDipoles();

      } catch (Exception e) {
        String message = "Exception computing opt induced dipoles.";
        logger.log(Level.SEVERE, message, e);
      }
    }

    // Use OPT dipole components to construct Induced Dipoles for the asymmetric unit.
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        inducedDipole[0][i][j] = 0.0;
        inducedDipoleCR[0][i][j] = 0.0;
        double sum = 0.0;
        double sump = 0.0;
        for (int k = 0; k <= optOrder; k++) {
          sum += optRegion.optDipole[k][i][j];
          sump += optRegion.optDipoleCR[k][i][j];
          inducedDipole[0][i][j] += optRegion.optCoefficients[k] * sum;
          inducedDipoleCR[0][i][j] += optRegion.optCoefficients[k] * sump;
        }
      }
    }

    // Expand asymmetric OPT dipoles for crystal symmetry.
    expandInducedDipoles();

    return optOrder;
  }

  /**
   * A precision of 1.0e-8 results in an Ewald coefficient that ensures continuity in the real space
   * gradient, but at the cost of increased amplitudes for high frequency reciprocal space structure
   * factors.
   */
  private double ewaldCoefficient(double cutoff, double precision) {

    double eps = 1.0e-8;
    if (precision < 1.0e-1) {
      eps = precision;
    }

    /*
     * Get an approximate value from cutoff and tolerance.
     */
    double ratio = eps + 1.0;
    double x = 0.5;
    int i = 0;
    // Larger values lead to a more "delta-function-like" Gaussian
    while (ratio >= eps) {
      i++;
      x *= 2.0;
      ratio = erfc(x * cutoff) / cutoff;
    }
    /*
     * Use a binary search to refine the coefficient.
     */
    int k = i + 60;
    double xlo = 0.0;
    double xhi = x;
    for (int j = 0; j < k; j++) {
      x = (xlo + xhi) / 2.0;
      ratio = erfc(x * cutoff) / cutoff;
      if (ratio >= eps) {
        xlo = x;
      } else {
        xhi = x;
      }
    }

    return x;
  }

  /** Given an array of atoms (with atom types), assign multipole types and reference sites. */
  private void assignMultipoles() {
    if (forceField == null) {
      String message = "No force field is defined.\n";
      logger.log(Level.SEVERE, message);
      return;
    }

    if (forceField.getForceFieldTypeCount(ForceFieldType.MULTIPOLE) < 1
        && forceField.getForceFieldTypeCount(ForceFieldType.CHARGE) < 1) {
      String message = "Force field has no permanent electrostatic types.\n";
      logger.log(Level.SEVERE, message);
      return;
    }
    if (nAtoms < 1) {
      String message = "No atoms are defined.\n";
      logger.log(Level.SEVERE, message);
      return;
    }
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];

      if (!assignMultipole(elecForm, atom, forceField, localMultipole[i], i, axisAtom, frame)) {
        logger.info(
            format("No MultipoleType could be assigned:\n %s --> %s", atom, atom.getAtomType()));
        StringBuilder sb = new StringBuilder();
        List<Bond> bonds = atom.getBonds();
        for (Bond b : bonds) {
          Atom a2 = b.get1_2(atom);
          AtomType aType2 = a2.getAtomType();
          sb.append(format("\n  %s --> %s", a2, aType2));
        }
        if (bonds.size() == 1) {
          Atom atom2 = bonds.get(0).get1_2(atom);
          bonds = atom2.getBonds();
          for (Bond b : bonds) {
            Atom a2 = b.get1_2(atom2);
            AtomType aType2 = a2.getAtomType();
            sb.append(format("\n  1-3 %s --> %s", a2, aType2));
          }
        }

        List<MultipoleType> multipoleTypes = forceField.getMultipoleTypes("" + atom.getAtomType().getKey());
        if (multipoleTypes != null || !multipoleTypes.isEmpty()) {
          sb.append("\n Similar Multipole types:");
          for (MultipoleType multipoleType : multipoleTypes) {
            sb.append(format("\n %s", multipoleType));
          }
        }

        logger.log(Level.SEVERE, sb.toString());
      }
    }

    // Check for multipoles that were not assigned correctly.
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < nAtoms; i++) {
      boolean flag = false;
      for (int j = 0; j < 10; j++) {
        if (Double.isNaN(localMultipole[i][j])) {
          flag = true;
          break;
        }
      }
      if (flag) {
        sb.append("\n").append(atoms[i].toString()).append("\n");
        sb.append(format("%d", i + 1));
        for (int j = 0; j < 10; j++) {
          sb.append(format(" %8.3f", localMultipole[i][j]));
        }
        sb.append("\n");
      }
    }
    if (sb.length() > 0) {
      String message = "Fatal exception: Error assigning multipoles. " + sb;
      logger.log(Level.SEVERE, message);
    }
  }

  /**
   * Log the real space electrostatics interaction.
   *
   * @param i Atom i.
   * @param k Atom j.
   * @param r The distance rij.
   * @param eij The interaction energy.
   * @since 1.0
   */
  private void log(int i, int k, double r, double eij) {
    logger.info(
        format(
            "%s %6d-%s %6d-%s %10.4f  %10.4f",
            "ELEC",
            atoms[i].getIndex(),
            atoms[i].getAtomType().name,
            atoms[k].getIndex(),
            atoms[k].getAtomType().name,
            r,
            eij));
  }

  public static class RealSpaceNeighborParameters {

    final int numThreads;
    /**
     * Neighbor lists, without atoms beyond the real space cutoff. [nSymm][nAtoms][nIncludedNeighbors]
     */
    public int[][][] realSpaceLists;
    /** Number of neighboring atoms within the real space cutoff. [nSymm][nAtoms] */
    public int[][] realSpaceCounts;
    /** Optimal pairwise ranges. */
    public Range[] realSpaceRanges;
    /** Pairwise schedule for load balancing. */
    public IntegerSchedule realSpaceSchedule;

    public RealSpaceNeighborParameters(int maxThreads) {
      numThreads = maxThreads;
      realSpaceRanges = new Range[maxThreads];
    }

    public void allocate(int nAtoms, int nSymm) {
      realSpaceSchedule = new PairwiseSchedule(numThreads, nAtoms, realSpaceRanges);
      realSpaceLists = new int[nSymm][nAtoms][];
      realSpaceCounts = new int[nSymm][nAtoms];
    }
  }

  /** Mutable Particle Mesh Ewald constants. */
  public class EwaldParameters {

    public double aewald;
    public double aewald3;
    public double an0;
    public double an1;
    public double an2;
    public double an3;
    public double an4;
    public double an5;
    public double off;
    public double off2;

    public EwaldParameters() {
      double off;
      if (!crystal.aperiodic()) {
        off = forceField.getDouble("EWALD_CUTOFF", PERIODIC_DEFAULT_EWALD_CUTOFF);
      } else {
        off = forceField.getDouble("EWALD_CUTOFF", APERIODIC_DEFAULT_EWALD_CUTOFF);
      }
      double aewald = forceField.getDouble("EWALD_ALPHA", 0.545);
      setEwaldParameters(off, aewald);
    }

    /**
     * Determine the real space Ewald parameters and permanent multipole self energy.
     *
     * @param off Real space cutoff.
     * @param aewald Ewald convergence parameter (0.0 turns off reciprocal space).
     */
    public void setEwaldParameters(double off, double aewald) {
      this.off = off;
      this.aewald = aewald;
      off2 = off * off;
      double alsq2 = 2.0 * aewald * aewald;
      double piEwald = Double.POSITIVE_INFINITY;
      if (aewald > 0.0) {
        piEwald = 1.0 / (SQRT_PI * aewald);
      }
      aewald3 = 4.0 / 3.0 * pow(aewald, 3.0) / SQRT_PI;
      if (aewald > 0.0) {
        an0 = alsq2 * piEwald;
        an1 = alsq2 * an0;
        an2 = alsq2 * an1;
        an3 = alsq2 * an2;
        an4 = alsq2 * an3;
        an5 = alsq2 * an4;
      } else {
        an0 = 0.0;
        an1 = 0.0;
        an2 = 0.0;
        an3 = 0.0;
        an4 = 0.0;
        an5 = 0.0;
      }
    }
  }

  /** Scale factors for group based polarization rules and energy masking rules. */
  public class ScaleParameters {

    /** The interaction energy between 1-2 multipoles is scaled by m12scale. */
    public final double m12scale;
    /** The interaction energy between 1-3 multipoles is scaled by m13scale. */
    public final double m13scale;
    /** The interaction energy between 1-4 multipoles is scaled by m14scale. */
    public final double m14scale;
    /** The interaction energy between 1-5 multipoles is scaled by m15scale. */
    public final double m15scale;

    /**
     * Direct polarization field due to permanent multipoles at polarizable sites within their group
     * are scaled. The scaling is 0.0 in AMOEBA.
     */
    public final double d11scale;
    /**
     * The interaction energy between a permanent multipole and polarizable site that are 1-2 is
     * scaled by p12scale.
     */
    public final double p12scale;
    /**
     * The interaction energy between a permanent multipole and polarizable site that are 1-3 is
     * scaled by p13scale.
     */
    public final double p13scale;

    public final double p14scale;
    public final double p15scale;
    public final double intra14Scale;

    public ScaleParameters(ForceField forceField) {
      if (elecForm == PAM) {
        m12scale = forceField.getDouble("MPOLE_12_SCALE", 0.0);
        m13scale = forceField.getDouble("MPOLE_13_SCALE", 0.0);
        m14scale = forceField.getDouble("MPOLE_14_SCALE", 0.4);
        m15scale = forceField.getDouble("MPOLE_15_SCALE", 0.8);
      } else {
        double mpole14 = forceField.getDouble("CHG_14_SCALE", 2.0);
        mpole14 = 1.0 / mpole14;
        m12scale = forceField.getDouble("MPOLE_12_SCALE", 0.0);
        m13scale = forceField.getDouble("MPOLE_13_SCALE", 0.0);
        m14scale = forceField.getDouble("MPOLE_14_SCALE", mpole14);
        m15scale = forceField.getDouble("MPOLE_15_SCALE", 1.0);
      }
      intra14Scale = forceField.getDouble("POLAR_14_INTRA", 0.5);
      d11scale = forceField.getDouble("DIRECT_11_SCALE", 0.0);
      p12scale = forceField.getDouble("POLAR_12_SCALE", 0.0);
      p13scale = forceField.getDouble("POLAR_13_SCALE", 0.0);
      p14scale = forceField.getDouble("POLAR_14_SCALE", 1.0);
      p15scale = forceField.getDouble("POLAR_15_SCALE", 1.0);
    }
  }

  public class AlchemicalParameters {

    /** Constant α in: r' = sqrt(r^2 + α*(1 - L)^2) */
    public double permLambdaAlpha = 1.0;
    /** Power on L in front of the pairwise multipole potential. */
    public double permLambdaExponent = 3.0;
    /** Begin turning on permanent multipoles at Lambda = 0.4; */
    public double permLambdaStart = 0.4;
    /** Finish turning on permanent multipoles at Lambda = 1.0; */
    public double permLambdaEnd = 1.0;
    /**
     * Start turning on polarization later in the Lambda path to prevent SCF convergence problems
     * when atoms nearly overlap.
     */
    public double polLambdaStart = 0.75;

    public double polLambdaEnd = 1.0;
    /** Power on L in front of the polarization energy. */
    public double polLambdaExponent = 3.0;
    /** Intramolecular electrostatics for the ligand in vapor is included by default. */
    public boolean doLigandVaporElec = true;
    /** Intramolecular electrostatics for the ligand in done in GK implicit solvent. */
    public boolean doLigandGKElec = false;
    /**
     * Condensed phase SCF without the ligand present is included by default. For DualTopologyEnergy
     * calculations it can be turned off.
     */
    public boolean doNoLigandCondensedSCF = true;
    /** lAlpha = α*(1 - L)^2 */
    public double lAlpha = 0.0;

    public double dlAlpha = 0.0;
    public double d2lAlpha = 0.0;
    public double dEdLSign = 1.0;
    /** lPowPerm = L^permanentLambdaExponent */
    public double lPowPerm = 1.0;

    public double dlPowPerm = 0.0;
    public double d2lPowPerm = 0.0;
    public boolean doPermanentRealSpace = true;
    public double permanentScale = 1.0;
    /** lPowPol = L^polarizationLambdaExponent */
    public double lPowPol = 1.0;

    public double dlPowPol = 0.0;
    public double d2lPowPol = 0.0;
    public boolean doPolarization = true;
    /**
     * When computing the polarization energy at L there are 3 pieces.
     *
     * <p>1.) Upol(1) = The polarization energy computed normally (ie. system with ligand).
     *
     * <p>2.) Uenv = The polarization energy of the system without the ligand.
     *
     * <p>3.) Uligand = The polarization energy of the ligand by itself.
     *
     * <p>Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     *
     * <p>Set polarizationScale to L for part 1. Set polarizationScale to (1-L) for parts 2 and 3.
     */
    public double polarizationScale = 1.0;
    /**
     * The polarization Lambda value goes from 0.0 .. 1.0 as the global lambda value varies between
     * polLambdaStart .. polLambadEnd.
     */
    private double polLambda = 1.0;
    /**
     * The permanent Lambda value goes from 0.0 .. 1.0 as the global lambda value varies between
     * permLambdaStart .. permLambdaEnd.
     */
    private double permLambda = 1.0;
    /** Boundary conditions for the vapor end of the alchemical path. */
    private Crystal vaporCrystal = null;

    private int[][][] vaporLists = null;
    private Range[] vacuumRanges = null;
    private IntegerSchedule vaporPermanentSchedule = null;
    private IntegerSchedule vaporEwaldSchedule = null;

    public AlchemicalParameters(ForceField forceField, boolean lambdaTerm) {
      if (lambdaTerm) {
        // Values of PERMANENT_LAMBDA_ALPHA below 2 can lead to unstable  trajectories.
        permLambdaAlpha = forceField.getDouble("PERMANENT_LAMBDA_ALPHA", 2.0);
        if (permLambdaAlpha < 0.0 || permLambdaAlpha > 3.0) {
          logger.warning(
              "Invalid value for permanent-lambda-alpha (<0.0 || >3.0); reverting to 2.0");
          permLambdaAlpha = 2.0;
        }

        /*
         A PERMANENT_LAMBDA_EXPONENT of 2 gives a non-zero d2U/dL2 at the
         beginning of the permanent schedule. Choosing a power of 3 or
         greater ensures a smooth dU/dL and d2U/dL2 over the schedule.

         A value of 0.0 is also admissible for when ExtendedSystem is
         scaling multipoles rather than softcoring them.
        */
        permLambdaExponent = forceField.getDouble("PERMANENT_LAMBDA_EXPONENT", 3.0);
        if (permLambdaExponent < 0.0) {
          logger.warning("Invalid value for permanent-lambda-exponent (<0.0); reverting to 3.0");
          permLambdaExponent = 3.0;
        }

        /*
         A POLARIZATION_LAMBDA_EXPONENT of 2 gives a non-zero d2U/dL2 at
         the beginning of the polarization schedule. Choosing a power of 3
         or greater ensures a smooth dU/dL and d2U/dL2 over the schedule.

         A value of 0.0 is also admissible: when polarization is not being
         softcored but instead scaled, as by ExtendedSystem.
        */
        polLambdaExponent = forceField.getDouble("POLARIZATION_LAMBDA_EXPONENT", 3.0);
        if (polLambdaExponent < 0.0) {
          logger.warning("Invalid value for polarization-lambda-exponent (<0.0); reverting to 3.0");
          polLambdaExponent = 3.0;
        }

        if (noWindowing) {
          permLambdaStart = 0.0;
          polLambdaStart = 0.0;
          permLambdaEnd = 1.0;
          polLambdaEnd = 1.0;
          logger.info(
              "PME-Cart lambda windowing disabled. Permanent and polarization lambda affect entire [0,1].");
        } else {
          // Values of PERMANENT_LAMBDA_START below 0.5 can lead to unstable trajectories.
          permLambdaStart = forceField.getDouble("PERMANENT_LAMBDA_START", 0.4);
          if (permLambdaStart < 0.0 || permLambdaStart > 1.0) {
            logger.warning("Invalid value for perm-lambda-start (<0.0 || >1.0); reverting to 0.4");
            permLambdaStart = 0.4;
          }

          // Values of PERMANENT_LAMBDA_END must be greater than permLambdaStart and <= 1.0.
          permLambdaEnd = forceField.getDouble("PERMANENT_LAMBDA_END", 1.0);
          if (permLambdaEnd < permLambdaStart || permLambdaEnd > 1.0) {
            logger.warning("Invalid value for perm-lambda-end (<start || >1.0); reverting to 1.0");
            permLambdaEnd = 1.0;
          }

          /*
           The POLARIZATION_LAMBDA_START defines the point in the lambda
           schedule when the condensed phase polarization of the ligand
           begins to be turned on. If the condensed phase polarization
           is considered near lambda=0, then SCF convergence is slow,
           even with Thole damping. In addition, 2 (instead of 1)
           condensed phase SCF calculations are necessary from the
           beginning of the window to lambda=1.
          */
          polLambdaStart = forceField.getDouble("POLARIZATION_LAMBDA_START", 0.75);
          if (polLambdaStart < 0.0 || polLambdaStart > 1.0) {
            logger.warning("Invalid value for polarization-lambda-start; reverting to 0.75");
            polLambdaStart = 0.75;
          }

          /*
           The POLARIZATION_LAMBDA_END defines the point in the lambda
           schedule when the condensed phase polarization of ligand has
           been completely turned on. Values other than 1.0 have not been tested.
          */
          polLambdaEnd = forceField.getDouble("POLARIZATION_LAMBDA_END", 1.0);
          if (polLambdaEnd < polLambdaStart || polLambdaEnd > 1.0) {
            logger.warning(
                "Invalid value for polarization-lambda-end (<start || >1.0); reverting to 1.0");
            polLambdaEnd = 1.0;
          }
        }

        // The LAMBDA_VAPOR_ELEC defines if intramolecular electrostatics of the ligand in vapor
        // will be considered.
        doLigandVaporElec = forceField.getBoolean("LIGAND_VAPOR_ELEC", true);
        doLigandGKElec = forceField.getBoolean("LIGAND_GK_ELEC", false);
        doNoLigandCondensedSCF = forceField.getBoolean("NO_LIGAND_CONDENSED_SCF", true);
      }
    }

    public void printLambdaFactors() {
      StringBuilder sb = new StringBuilder();
      sb.append(
          format(
              "  (%4s)  mode:%-20s lambda:%.2f  permScale:%.2f  polScale:%.2f  dEdLSign:%s  doPol:%-5b  doPermRS:%-5b",
              "CART",
              lambdaMode.toString(),
              lambda,
              permanentScale,
              polarizationScale,
              format("%+f", dEdLSign).charAt(0),
              doPolarization,
              doPermanentRealSpace));
      sb.append(
          format(
              "\n    lAlpha:%.2f,%.2f,%.2f  lPowPerm:%.2f,%.2f,%.2f  lPowPol:%.2f,%.2f,%.2f",
              lAlpha,
              dlAlpha,
              d2lAlpha,
              lPowPerm,
              dlPowPerm,
              d2lPowPerm,
              lPowPol,
              dlPowPol,
              d2lPowPol));
      sb.append(
          format(
              "\n    permExp:%.2f  permAlpha:%.2f  permWindow:%.2f,%.2f  polExp:%.2f  polWindow:%.2f,%.2f",
              permLambdaExponent,
              permLambdaAlpha,
              permLambdaStart,
              permLambdaEnd,
              polLambdaExponent,
              polLambdaStart,
              polLambdaEnd));
      logger.info(sb.toString());
    }

    public String toString() {
      StringBuilder sb = new StringBuilder("   Alchemical Parameters\n");
      sb.append(
          format(
              "    Permanent Multipole Range:      %5.3f-%5.3f\n", permLambdaStart, permLambdaEnd));
      sb.append(format("    Permanent Multipole Softcore Alpha:   %5.3f\n", permLambdaAlpha));
      sb.append(format("    Permanent Multipole Lambda Exponent:  %5.3f\n", permLambdaExponent));
      if (polarization != Polarization.NONE) {
        sb.append(format("    Polarization Lambda Exponent:         %5.3f\n", polLambdaExponent));
        sb.append(
            format(
                "    Polarization Range:             %5.3f-%5.3f\n", polLambdaStart, polLambdaEnd));
        sb.append(format("    Condensed SCF Without Ligand:         %B\n", doNoLigandCondensedSCF));
      }
      if (!doLigandGKElec) {
        sb.append(format("    Vapor Electrostatics:                 %B\n", doLigandVaporElec));
      } else {
        sb.append(format("    GK Electrostatics at L=0:             %B\n", doLigandGKElec));
      }
      return sb.toString();
    }

    /*
     * f = sqrt(r^2 + lAlpha)
     *
     * df/dL = -alpha * (1.0 - lambda) / f
     *
     * g = 1 / sqrt(r^2 + lAlpha)
     *
     * dg/dL = alpha * (1.0 - lambda) / (r^2 + lAlpha)^(3/2)
     *
     * define dlAlpha = alpha * 1.0 - lambda)
     *
     * then df/dL = -dlAlpha / f and dg/dL = dlAlpha * g^3
     *
     * Multipoles are turned on from permLambdaStart .. permLambdaEnd.
     *
     * @param lambda
     */
    public void update(double lambda) {

      lPowPerm = 1.0;
      dlPowPerm = 0.0;
      d2lPowPerm = 0.0;
      lAlpha = 0.0;
      dlAlpha = 0.0;
      d2lAlpha = 0.0;
      if (lambda < permLambdaStart) {
        lPowPerm = 0.0;
      } else if (lambda <= permLambdaEnd) {
        double permWindow = permLambdaEnd - permLambdaStart;
        double permLambdaScale = 1.0 / permWindow;
        permLambda = permLambdaScale * (lambda - permLambdaStart);

        lAlpha = permLambdaAlpha * (1.0 - permLambda) * (1.0 - permLambda);
        dlAlpha = permLambdaAlpha * (1.0 - permLambda);
        d2lAlpha = -permLambdaAlpha;

        lPowPerm = pow(permLambda, permLambdaExponent);
        dlPowPerm = permLambdaExponent * pow(permLambda, permLambdaExponent - 1.0);
        d2lPowPerm = 0.0;
        if (permLambdaExponent >= 2.0) {
          d2lPowPerm =
              permLambdaExponent
                  * (permLambdaExponent - 1.0)
                  * pow(permLambda, permLambdaExponent - 2.0);
        }

        dlAlpha *= permLambdaScale;
        d2lAlpha *= (permLambdaScale * permLambdaScale);
        dlPowPerm *= permLambdaScale;
        d2lPowPerm *= (permLambdaScale * permLambdaScale);
      }

      // Polarization is turned on from polarizationLambdaStart .. polarizationLambdaEnd.
      lPowPol = 1.0;
      dlPowPol = 0.0;
      d2lPowPol = 0.0;
      if (lambda < polLambdaStart) {
        lPowPol = 0.0;
      } else if (lambda <= polLambdaEnd) {
        double polWindow = polLambdaEnd - polLambdaStart;
        double polLambdaScale = 1.0 / polWindow;
        polLambda = polLambdaScale * (lambda - polLambdaStart);
        if (polLambdaExponent > 0.0) {
          lPowPol = pow(polLambda, polLambdaExponent);
          if (polLambdaExponent >= 1.0) {
            dlPowPol = polLambdaExponent * pow(polLambda, polLambdaExponent - 1.0);
            if (polLambdaExponent >= 2.0) {
              d2lPowPol =
                  polLambdaExponent
                      * (polLambdaExponent - 1.0)
                      * pow(polLambda, polLambdaExponent - 2.0);
            }
          }
        }
        // Add the chain rule term due to shrinking the lambda range for the polarization energy.
        dlPowPol *= polLambdaScale;
        d2lPowPol *= (polLambdaScale * polLambdaScale);
      }
    }
  }

  public class SCFPredictorParameters {

    /** Induced dipole predictor order. */
    public int predictorOrder;
    /** Induced dipole predictor index. */
    public int predictorStartIndex;
    /** Induced dipole predictor count. */
    public int predictorCount;
    /** Dimensions of [mode][predictorOrder][nAtoms][3] */
    public double[][][][] predictorInducedDipole;
    /** Dimensions of [mode][predictorOrder][nAtoms][3] */
    public double[][][][] predictorInducedDipoleCR;

    public LeastSquaresPredictor leastSquaresPredictor;
    public LevenbergMarquardtOptimizer leastSquaresOptimizer;

    /** Always-stable predictor-corrector for the mutual induced dipoles. */
    public void aspcPredictor() {

      if (predictorCount < 6) {
        return;
      }

      int mode;
      switch (lambdaMode) {
        case OFF:
        case CONDENSED:
          mode = 0;
          break;
        case CONDENSED_NO_LIGAND:
          mode = 1;
          break;
        case VAPOR:
          mode = 2;
          break;
        default:
          mode = 0;
      }

      final double[] aspc = {
          22.0 / 7.0, -55.0 / 14.0, 55.0 / 21.0, -22.0 / 21.0, 5.0 / 21.0, -1.0 / 42.0
      };

      // Initialize a pointer into predictor induced dipole array.
      int index = predictorStartIndex;

      // Expansion loop.
      for (int k = 0; k < 6; k++) {

        // Set the current predictor coefficient.
        double c = aspc[k];
        for (int i = 0; i < nAtoms; i++) {
          for (int j = 0; j < 3; j++) {
            inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
            inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
          }
        }
        index++;
        if (index >= predictorOrder) {
          index = 0;
        }
      }
    }

    public void init() {
      predictorCount = 0;
      int defaultOrder = 6;
      predictorOrder = forceField.getInteger("SCF_PREDICTOR_ORDER", defaultOrder);
      if (scfPredictor == SCFPredictor.LS) {
        leastSquaresPredictor = new LeastSquaresPredictor();
        double eps = 1.0e-4;
        leastSquaresOptimizer =
            new org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer(
                new org.apache.commons.math3.optimization.SimpleVectorValueChecker(eps, eps));
      } else if (scfPredictor == SCFPredictor.ASPC) {
        predictorOrder = 6;
      }
      predictorStartIndex = 0;
    }

    /**
     * The least-squares predictor with induced dipole information from 8-10 previous steps reduces
     * the number SCF iterations by ~50%.
     */
    public void leastSquaresPredictor() {
      if (predictorCount < 2) {
        return;
      }
      try {
        /*
         The Jacobian and target do not change during the LS optimization,
         so it's most efficient to update them once before the
         Least-Squares optimizer starts.
        */
        leastSquaresPredictor.updateJacobianAndTarget();
        int maxEvals = 100;
        fill(leastSquaresPredictor.initialSolution, 0.0);
        leastSquaresPredictor.initialSolution[0] = 1.0;
        org.apache.commons.math3.optimization.PointVectorValuePair optimum =
            leastSquaresOptimizer.optimize(
                maxEvals,
                leastSquaresPredictor,
                leastSquaresPredictor.calculateTarget(),
                leastSquaresPredictor.weights,
                leastSquaresPredictor.initialSolution);
        double[] optimalValues = optimum.getPoint();
        if (logger.isLoggable(Level.FINEST)) {
          logger.finest(format("\n LS RMS:            %10.6f", leastSquaresOptimizer.getRMS()));
          logger.finest(format(" LS Iterations:     %10d", leastSquaresOptimizer.getEvaluations()));
          logger.finest(
              format(" Jacobian Evals:    %10d", leastSquaresOptimizer.getJacobianEvaluations()));
          logger.finest(format(" Chi Square:        %10.6f", leastSquaresOptimizer.getChiSquare()));
          logger.finest(" LS Coefficients");
          for (int i = 0; i < predictorOrder - 1; i++) {
            logger.finest(format(" %2d  %10.6f", i + 1, optimalValues[i]));
          }
        }

        int mode;
        switch (lambdaMode) {
          case OFF:
          case CONDENSED:
            mode = 0;
            break;
          case CONDENSED_NO_LIGAND:
            mode = 1;
            break;
          case VAPOR:
            mode = 2;
            break;
          default:
            mode = 0;
        }

        // Initialize a pointer into predictor induced dipole array.
        int index = predictorStartIndex;

        // Apply the LS coefficients in order to provide an initial guess at the converged induced
        // dipoles.
        for (int k = 0; k < predictorOrder - 1; k++) {

          // Set the current coefficient.
          double c = optimalValues[k];
          for (int i = 0; i < nAtoms; i++) {
            for (int j = 0; j < 3; j++) {
              inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
              inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
            }
          }
          index++;
          if (index >= predictorOrder) {
            index = 0;
          }
        }
      } catch (Exception e) {
        logger.log(Level.WARNING, " Exception computing predictor coefficients", e);
      }
    }

    /** Polynomial predictor for the mutual induced dipoles. */
    public void polynomialPredictor() {

      if (predictorCount == 0) {
        return;
      }

      int mode;
      switch (lambdaMode) {
        case OFF:
        case CONDENSED:
          mode = 0;
          break;
        case CONDENSED_NO_LIGAND:
          mode = 1;
          break;
        case VAPOR:
          mode = 2;
          break;
        default:
          mode = 0;
      }

      // Check the number of previous induced dipole vectors available.
      int n = predictorOrder;
      if (predictorCount < predictorOrder) {
        n = predictorCount;
      }

      // Initialize a pointer into predictor induced dipole array.
      int index = predictorStartIndex;

      // Initialize the sign of the polynomial expansion.
      double sign = -1.0;

      // Expansion loop.
      for (int k = 0; k < n; k++) {

        // Set the current predictor sign and coefficient.
        sign *= -1.0;
        double c = sign * ScalarMath.binomial(n, k);
        for (int i = 0; i < nAtoms; i++) {
          for (int j = 0; j < 3; j++) {
            inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
            inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
          }
        }
        index++;
        if (index >= predictorOrder) {
          index = 0;
        }
      }
    }

    /** Save the current converged mutual induced dipoles. */
    public void saveMutualInducedDipoles() {

      int mode;
      switch (lambdaMode) {
        case OFF:
        case CONDENSED:
          mode = 0;
          break;
        case CONDENSED_NO_LIGAND:
          mode = 1;
          break;
        case VAPOR:
          mode = 2;
          break;
        default:
          mode = 0;
      }

      // Current induced dipoles are saved before those from the previous step.
      predictorStartIndex--;
      if (predictorStartIndex < 0) {
        predictorStartIndex = predictorOrder - 1;
      }

      if (predictorCount < predictorOrder) {
        predictorCount++;
      }

      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          predictorInducedDipole[mode][predictorStartIndex][i][j] =
              inducedDipole[0][i][j] - directDipole[i][j];
          predictorInducedDipoleCR[mode][predictorStartIndex][i][j] =
              inducedDipoleCR[0][i][j] - directDipoleCR[i][j];
        }
      }
    }

    private class LeastSquaresPredictor implements DifferentiableMultivariateVectorFunction {

      double[] weights;
      double[] target;
      double[] values;
      double[][] jacobian;
      double[] initialSolution;
      private final MultivariateMatrixFunction multivariateMatrixFunction = this::jacobian;

      LeastSquaresPredictor() {
        weights = new double[2 * nAtoms * 3];
        target = new double[2 * nAtoms * 3];
        values = new double[2 * nAtoms * 3];
        jacobian = new double[2 * nAtoms * 3][predictorOrder - 1];
        initialSolution = new double[predictorOrder - 1];
        fill(weights, 1.0);
        initialSolution[0] = 1.0;
      }

      @Override
      public MultivariateMatrixFunction jacobian() {
        return multivariateMatrixFunction;
      }

      @Override
      public double[] value(double[] variables) {
        int mode;
        switch (lambdaMode) {
          case OFF:
          case CONDENSED:
            mode = 0;
            break;
          case CONDENSED_NO_LIGAND:
            mode = 1;
            break;
          case VAPOR:
            mode = 2;
            break;
          default:
            mode = 0;
        }

        for (int i = 0; i < nAtoms; i++) {
          int index = 6 * i;
          values[index] = 0;
          values[index + 1] = 0;
          values[index + 2] = 0;
          values[index + 3] = 0;
          values[index + 4] = 0;
          values[index + 5] = 0;
          int pi = predictorStartIndex + 1;
          if (pi >= predictorOrder) {
            pi = 0;
          }
          for (int j = 0; j < predictorOrder - 1; j++) {
            values[index] += variables[j] * predictorInducedDipole[mode][pi][i][0];
            values[index + 1] += variables[j] * predictorInducedDipole[mode][pi][i][1];
            values[index + 2] += variables[j] * predictorInducedDipole[mode][pi][i][2];
            values[index + 3] += variables[j] * predictorInducedDipoleCR[mode][pi][i][0];
            values[index + 4] += variables[j] * predictorInducedDipoleCR[mode][pi][i][1];
            values[index + 5] += variables[j] * predictorInducedDipoleCR[mode][pi][i][2];
            pi++;
            if (pi >= predictorOrder) {
              pi = 0;
            }
          }
        }
        return values;
      }

      double[] calculateTarget() {
        return target;
      }

      void updateJacobianAndTarget() {
        int mode;
        switch (lambdaMode) {
          case OFF:
          case CONDENSED:
            mode = 0;
            break;
          case CONDENSED_NO_LIGAND:
            mode = 1;
            break;
          case VAPOR:
            mode = 2;
            break;
          default:
            mode = 0;
        }

        // Update the target.
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
          target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][0];
          target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][1];
          target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][2];
          target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][0];
          target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][1];
          target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][2];
        }

        // Update the Jacobian.
        index = predictorStartIndex + 1;
        if (index >= predictorOrder) {
          index = 0;
        }
        for (int j = 0; j < predictorOrder - 1; j++) {
          int ji = 0;
          for (int i = 0; i < nAtoms; i++) {
            jacobian[ji++][j] = predictorInducedDipole[mode][index][i][0];
            jacobian[ji++][j] = predictorInducedDipole[mode][index][i][1];
            jacobian[ji++][j] = predictorInducedDipole[mode][index][i][2];
            jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][0];
            jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][1];
            jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][2];
          }
          index++;
          if (index >= predictorOrder) {
            index = 0;
          }
        }
      }

      private double[][] jacobian(double[] variables) {
        return jacobian;
      }
    }
  }

  public class PMETimings {

    public final long[] realSpacePermTime;
    public final long[] realSpaceEnergyTime;
    public final long[] realSpaceSCFTime;
    /** Timing variables. */
    private final int numThreads;

    public long realSpacePermTotal, realSpaceEnergyTotal, realSpaceSCFTotal;
    public long bornRadiiTotal, gkEnergyTotal;

    public PMETimings(int numThreads) {
      this.numThreads = numThreads;
      realSpacePermTime = new long[numThreads];
      realSpaceEnergyTime = new long[numThreads];
      realSpaceSCFTime = new long[numThreads];
    }

    public void init() {
      for (int i = 0; i < numThreads; i++) {
        realSpacePermTime[i] = 0;
        realSpaceEnergyTime[i] = 0;
        realSpaceSCFTime[i] = 0;
      }
      realSpacePermTotal = 0;
      realSpaceEnergyTotal = 0;
      realSpaceSCFTotal = 0;
      bornRadiiTotal = 0;
      gkEnergyTotal = 0;
    }

    public void printRealSpaceTimings() {
      double total = (realSpacePermTotal + realSpaceSCFTotal + realSpaceEnergyTotal) * NS2SEC;
      logger.info(format("\n Real Space: %7.4f (sec)", total));
      logger.info("           Electric Field");
      logger.info(" Thread    Direct  SCF     Energy     Counts");
      long minPerm = Long.MAX_VALUE;
      long maxPerm = 0;
      long minSCF = Long.MAX_VALUE;
      long maxSCF = 0;
      long minEnergy = Long.MAX_VALUE;
      long maxEnergy = 0;
      int minCount = Integer.MAX_VALUE;
      int maxCount = Integer.MIN_VALUE;

      for (int i = 0; i < maxThreads; i++) {
        int count = realSpaceEnergyRegion.getCount(i);
        logger.info(
            format(
                "    %3d   %7.4f %7.4f %7.4f %10d",
                i,
                realSpacePermTime[i] * NS2SEC,
                realSpaceSCFTime[i] * NS2SEC,
                realSpaceEnergyTime[i] * NS2SEC,
                count));
        minPerm = min(realSpacePermTime[i], minPerm);
        maxPerm = max(realSpacePermTime[i], maxPerm);
        minSCF = min(realSpaceSCFTime[i], minSCF);
        maxSCF = max(realSpaceSCFTime[i], maxSCF);
        minEnergy = min(realSpaceEnergyTime[i], minEnergy);
        maxEnergy = max(realSpaceEnergyTime[i], maxEnergy);
        minCount = min(count, minCount);
        maxCount = max(count, maxCount);
      }
      logger.info(
          format(
              " Min      %7.4f %7.4f %7.4f %10d",
              minPerm * NS2SEC, minSCF * NS2SEC, minEnergy * NS2SEC, minCount));
      logger.info(
          format(
              " Max      %7.4f %7.4f %7.4f %10d",
              maxPerm * NS2SEC, maxSCF * NS2SEC, maxEnergy * NS2SEC, maxCount));
      logger.info(
          format(
              " Delta    %7.4f %7.4f %7.4f %10d",
              (maxPerm - minPerm) * NS2SEC,
              (maxSCF - minSCF) * NS2SEC,
              (maxEnergy - minEnergy) * NS2SEC,
              (maxCount - minCount)));
      logger.info(
          format(
              " Actual   %7.4f %7.4f %7.4f %10d",
              realSpacePermTotal * NS2SEC,
              realSpaceSCFTotal * NS2SEC,
              realSpaceEnergyTotal * NS2SEC,
              realSpaceEnergyRegion.getInteractions()));
    }
  }

  /**
   * Flag to indicate use of extended variables that interpolate permanent multipoles and/or atomic
   * polarizability between states (e.g., for protonation and/or tautamers).
   */
  private boolean esvTerm = false;

  /**
   * The ExtendedSystem instance controls use of extended variables.
   */
  private ExtendedSystem extendedSystem = null;

  /**
   * Number of extended system variables.
   */
  private int numESVs = 0;

  /**
   * The partial derivative of the permanent multipole real space energy with respect to each ESV.
   */
  SharedDouble[] dUPermRealdEsv = null;
  /**
   * The partial derivative of the permanent multipole reciprocal space energy with respect to each
   * ESV.
   */
  SharedDouble[] dUPermRecipdEsv = null;
  /**
   * The partial derivative of the permanent multipole self energy with respect to each ESV.
   */
  SharedDouble[] dUPermSelfdEsv = null;

  /**
   * Attach system with extended variable such as titrations.
   *
   * @param system a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public void attachExtendedSystem(ExtendedSystem system) {
    // Set object handles.
    esvTerm = true;
    extendedSystem = system;
    numESVs = extendedSystem.size();

    // Update atoms and reinitialize arrays for consistency with the ExtendedSystem.
    setAtoms(extendedSystem.getExtendedAtoms(), extendedSystem.getExtendedMolecule());

    // Allocate shared derivative storage.
    dUPermRealdEsv = new SharedDouble[numESVs];
    dUPermRecipdEsv = new SharedDouble[numESVs];
    dUPermSelfdEsv = new SharedDouble[numESVs];
    for (int i = 0; i < numESVs; i++) {
      dUPermRealdEsv[i] = new SharedDouble(0.0);
      dUPermRecipdEsv[i] = new SharedDouble(0.0);
      dUPermSelfdEsv[i] = new SharedDouble(0.0);
    }

    updateEsvLambda();
    logger.info(format(" Attached extended system (%d variables) to PME.\n", numESVs));
  }

  /**
   * Flag to indicate if each atom is titrating.
   */
  boolean[] isAtomTitrating = null;
  /**
   * The scalar ESV that is operating on each atom (or 1.0 if the atom is not under ESV control).
   */
  double[] perAtomTitrationESV = null;
  /**
   * The index of the ESV that is operating on each atom (or -1 if the atom is not under ESV
   * control).
   */
  Integer[] perAtomESVIndex = null;

  /**
   * The partial derivative of each multipole with respect to its titration ESV (or 0.0 if the atom
   * is not under titration ESV control).
   */
  double[][][] dMultipoledTirationESV = null;

  /**
   * OST and ESV specific factors that effect real space interactions.
   */
  LambdaFactors[] lambdaFactors = null;

  /** Precalculate ESV factors subsequent to lambda propagation. */
  public void updateEsvLambda() {
    if (!esvTerm) {
      return;
    }

    // Query ExtendedSystem to create local preloads of all lambda quantities.
    numESVs = extendedSystem.size();
    if (perAtomTitrationESV == null || perAtomTitrationESV.length < nAtoms) {
      isAtomTitrating = new boolean[nAtoms];
      perAtomTitrationESV = new double[nAtoms];
      perAtomESVIndex = new Integer[nAtoms];
      fill(isAtomTitrating, false);
      fill(perAtomTitrationESV, 1.0);
      fill(perAtomESVIndex, null);
    }

    // Allocate space for dM/dTitratonESV
    if (dMultipoledTirationESV == null || dMultipoledTirationESV.length != nSymm
        || dMultipoledTirationESV[0].length != nAtoms) {
      dMultipoledTirationESV = new double[nSymm][nAtoms][10];
    }

    for (int i = 0; i < nAtoms; i++) {
      isAtomTitrating[i] = extendedSystem.isExtended(i);
      perAtomTitrationESV[i] = extendedSystem.getLambda(i);
      perAtomESVIndex[i] = extendedSystem.getEsvIndex(i);
      Atom ai = atoms[i];
      if (ai.getPolarizeType() == null) {
        logger.warning("Null polarize type during ESV init.");
        continue;
      }
      polarizability[i] = ai.getScaledPolarizability();
    }
  }

  /**
   * The setFactors(i,k,lambdaMode) method is called every time through the inner PME loops, avoiding
   * an "if (esv)" branch statement.
   * <p>
   * A plain OST run will have an object of type LambdaFactorsOST instead, which contains an empty
   * version of setFactors(i,k,lambdaMode).
   * <p>
   * The OST version instead sets new factors only on lambda updates, in setLambda(i,k).
   * <p>
   * Interactions involving neither lambda receive the (inoperative) defaults below.
   */
  public class LambdaFactors {

    /** lambda * esvLambda[i] * esvLambda[k] */
    protected double lambdaProd = 1.0;
    /** Interatomic buffer distance: alpha*(1-lambda)*(1-lambda). */
    protected double lfAlpha = 0.0;
    /** First lambda derivative of buffer distance. */
    protected double dlfAlpha = 0.0;
    /** Second lambda derivative of buffer distance. */
    protected double d2lfAlpha = 0.0;
    /** Lambda to its permanent exponent. */
    protected double lfPowPerm = 1.0;
    /** First lambda derivative of lPowPerm. */
    protected double dlfPowPerm = 0.0;
    /** Second lambda derivative of lPowPerm. */
    protected double d2lfPowPerm = 0.0;
    /** Lambda to its polarization exponent. */
    protected double lfPowPol = 1.0;
    /** First lambda derivative of lPowPol. */
    protected double dlfPowPol = 0.0;
    /** Second lambda derivative of lPowPol. */
    protected double d2lfPowPol = 0.0;
    /** Derivative of lambdaProduct w.r.t. lambda. */
    protected double dLpdL = 1.0;
    /** Derivative of lambdaProduct w.r.t. esvLambda[i]. */
    protected double dLpdLi = 1.0;
    /** Derivative of lambdaProduct w.r.t. esvLambda[k]. */
    protected double dLpdLk = 1.0;
    /** Atom indices used only by LambdaFactorsESV::print. */
    protected int[] ik = new int[2];

    public void print() {
      StringBuilder sb = new StringBuilder();
      if (this instanceof LambdaFactorsESV) {
        sb.append(format("  (QI-ESV)  i,k:%d,%d", ik[0], ik[1]));
      } else {
        sb.append("  (QI-OST)");
      }
      sb.append(format(
          "  lambda:%.2f  lAlpha:%.2f,%.2f,%.2f  lPowPerm:%.2f,%.2f,%.2f  lPowPol:%.2f,%.2f,%.2f",
          lambdaProd,
          lfAlpha,
          dlfAlpha,
          d2lfAlpha,
          lfPowPerm,
          dlfPowPerm,
          d2lfPowPerm,
          lfPowPol,
          dlfPowPol,
          d2lfPowPol));

      sb.append(format(
          "\n    permExp:%.2f  permAlpha:%.2f  permWindow:%.2f,%.2f  polExp:%.2f  polWindow:%.2f,%.2f",
          alchemicalParameters.permLambdaExponent,
          alchemicalParameters.permLambdaAlpha,
          alchemicalParameters.permLambdaStart,
          alchemicalParameters.permLambdaEnd,
          alchemicalParameters.polLambdaExponent,
          alchemicalParameters.polLambdaStart,
          alchemicalParameters.polLambdaEnd));
      logger.info(sb.toString());
    }

    /**
     * Overriden by the ESV version which updates with every softcore interaction.
     *
     * @param i Atom i index.
     * @param k Atom k index.
     * @param mode the LambdaMode
     */
    public void setFactors(int i, int k, LambdaMode mode) {
      /* no-op */
    }

    /** Overriden by the OST version which updates only during setLambda(). */
    public void setFactors() {
      /* no-op */
    }
  }

  public class LambdaFactorsOST extends LambdaFactors {

    @Override
    public void setFactors() {
      lambdaProd = lambda;
      lfAlpha = alchemicalParameters.lAlpha;
      dlfAlpha = alchemicalParameters.dlAlpha;
      d2lfAlpha = alchemicalParameters.d2lAlpha;
      lfPowPerm = alchemicalParameters.permanentScale;
      dlfPowPerm = alchemicalParameters.dlPowPerm * alchemicalParameters.dEdLSign;
      d2lfPowPerm = alchemicalParameters.d2lPowPerm * alchemicalParameters.dEdLSign;
      lfPowPol = alchemicalParameters.polarizationScale;
      dlfPowPol = alchemicalParameters.dlPowPol * alchemicalParameters.dEdLSign;
      d2lfPowPol = alchemicalParameters.d2lPowPol * alchemicalParameters.dEdLSign;
    }
  }

  public class LambdaFactorsESV extends LambdaFactors {

    @Override
    public void setFactors(int i, int k, LambdaMode mode) {
      logger.info(format("Invoked Qi setFactors() method with i,k=%d,%d", i, k));
      ik[0] = i;
      ik[1] = k;
      final double L = lambda;
      lambdaProd = L * perAtomTitrationESV[i] * perAtomTitrationESV[k];
      final double esvli = perAtomTitrationESV[i];
      final double esvlk = perAtomTitrationESV[k];
      dLpdL = esvli * esvlk;
      dLpdLi = L * esvlk;
      dLpdLk = L * esvli;

      // Permanent Lambda Window (e.g., 0 .. 1).
      double permLambdaExponent = alchemicalParameters.permLambdaExponent;
      double permLambdaStart = alchemicalParameters.permLambdaStart;
      double permLambdaEnd = alchemicalParameters.permLambdaEnd;

      double permWindow = 1.0 / (permLambdaEnd - permLambdaStart);
      double permLambda = (lambdaProd - permLambdaStart) * permWindow;
      lfPowPerm = pow(permLambda, alchemicalParameters.permLambdaExponent);
      dlfPowPerm = (permLambdaExponent < 1)
          ? 0.0 : permLambdaExponent * pow(permLambda, permLambdaExponent - 1) * permWindow;
      d2lfPowPerm = (permLambdaExponent < 2)
          ? 0.0 : permLambdaExponent
          * (permLambdaExponent - 1)
          * pow(permLambda, permLambdaExponent - 2)
          * permWindow
          * permWindow;

      // Polarization Lambda Window (e.g., 0 .. 1).
      double polLambdaExponent = alchemicalParameters.polLambdaExponent;
      double polLambdaStart = alchemicalParameters.polLambdaStart;
      double polLambdaEnd = alchemicalParameters.polLambdaEnd;

      double polWindow = 1.0 / (polLambdaEnd - polLambdaStart);
      double polLambda = (lambdaProd - polLambdaStart) * polWindow;
      lfPowPol = pow(polLambda, polLambdaExponent);
      dlfPowPol = (polLambdaExponent < 1)
          ? 0.0 : polLambdaExponent * pow(polLambda, polLambdaExponent - 1) * polWindow;
      d2lfPowPol = (polLambdaExponent < 2)
          ? 0.0 : polLambdaExponent
          * (polLambdaExponent - 1)
          * pow(polLambda, polLambdaExponent - 2)
          * polWindow
          * polWindow;

      // Permanent Lambda Softcore Alpha.
      double permLambdaAlpha = alchemicalParameters.permLambdaAlpha;
      lfAlpha = permLambdaAlpha * (1.0 - permLambda) * (1.0 - permLambda);
      dlfAlpha = permLambdaAlpha * (1.0 - permLambda) * permWindow;
      d2lfAlpha = -permLambdaAlpha * permWindow * permWindow;

      /*
        Follow the logic of
        1) condensedEnergy
        2) noLigand
        and
        3) ligandElec

        to set permanentScale (lfPowPerm) and polarizationScale (lfPowPol),
        substituting lambdaProduct for lambda.
       */
      switch (mode) {
        case CONDENSED:
          lfPowPerm = 1.0 - lfPowPerm;
          dlfPowPerm = -dlfPowPerm; // handles dEdLSign
          d2lfPowPerm = -d2lfPowPerm; // handles dEdLSign
          if (polarization == Polarization.NONE || lambda < polLambdaStart) {
            lfPowPol = 0.0;
            dlfPowPol = 0.0;
            d2lfPowPol = 0.0;
          } else if (lambda <= polLambdaEnd) {
            // No change necessary.
          } else {
            lfPowPol = 1.0;
            dlfPowPol = 0.0;
            d2lfPowPol = 0.0;
          }
          break;
        case CONDENSED_NO_LIGAND:
          // There's no treatment of polLambdaStart in ::condensedNoLigandSCF?
          lfPowPerm = 1.0 - lfPowPerm;
          dlfPowPerm = -dlfPowPerm;
          d2lfPowPerm = -d2lfPowPerm;
          if (polarization == Polarization.NONE) {
            lfPowPol = 0.0;
            dlfPowPol = 0.0;
            d2lfPowPol = 0.0;
          } else if (lambda <= polLambdaEnd) {
            lfPowPol = 1.0 - lfPowPol;
            dlfPowPol = -dlfPowPol;
            d2lfPowPol = -d2lfPowPol;
          } else {
            lfPowPol = 0.0;
            dlfPowPol = 0.0;
            d2lfPowPol = 0.0;
          }
          break;
        case VAPOR:
          // There's no treatment of polLambdaStart in ::ligandElec?
          lfPowPerm = 1.0 - lfPowPerm;
          dlfPowPerm = -dlfPowPerm;
          d2lfPowPerm = -d2lfPowPerm;
          lfAlpha = 0.0;
          dlfAlpha = 0.0;
          d2lfAlpha = 0.0;
          if (polarization == Polarization.NONE || lambdaProd > polLambdaEnd) {
            lfPowPol = 0.0;
            dlfPowPol = 0.0;
            d2lfPowPol = 0.0;
          } else if (lambdaProd <= polLambdaEnd) {
            lfPowPol = 1.0 - lfPowPol;
            dlfPowPol = -dlfPowPol;
            d2lfPowPol = -d2lfPowPol;
          }
          break;
        case OFF:
        default:
      }
    }
  }

  /**
   * The defaults are effectively final, as the implementation of setFactors in the base class is
   * always a no-op.
   */
  public final LambdaFactors LambdaDefaults = new LambdaFactors();
}
