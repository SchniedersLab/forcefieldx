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

import static ffx.potential.MolecularAssembly.atomIndexing;
import static ffx.potential.nonbonded.pme.EwaldParameters.DEFAULT_EWALD_COEFFICIENT;
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
import static ffx.utilities.Constants.ELEC_ANG_TO_DEBYE;
import static ffx.utilities.Constants.NS2SEC;
import static ffx.utilities.KeywordGroup.ElectrostaticsFunctionalForm;
import static ffx.utilities.KeywordGroup.LocalGeometryFunctionalForm;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static java.util.Collections.sort;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.util.Range;
import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.multipole.MultipoleTensor;
import ffx.potential.ForceFieldEnergy.Platform;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.pme.AlchemicalParameters;
import ffx.potential.nonbonded.pme.DirectRegion;
import ffx.potential.nonbonded.pme.LambdaMode;
import ffx.potential.nonbonded.pme.PMETimings;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.nonbonded.pme.RealSpaceNeighborParameters;
import ffx.potential.nonbonded.pme.SCFAlgorithm;
import ffx.potential.nonbonded.pme.SCFPredictor;
import ffx.potential.nonbonded.pme.SCFPredictorParameters;
import ffx.potential.nonbonded.pme.ScaleParameters;
import ffx.potential.nonbonded.pme.EwaldParameters;
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
import ffx.utilities.FFXKeyword;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

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
public class ParticleMeshEwald implements LambdaInterface {

  private static final Logger logger = Logger.getLogger(ParticleMeshEwald.class.getName());
  /** Number of unique tensors for given order. */
  private static final int tensorCount = MultipoleTensor.tensorCount(3);
  /**
   * If lambdaTerm is true, some ligand atom interactions with the environment are being turned
   * on/off.
   */
  private final boolean lambdaTerm;
  private final boolean reciprocalSpaceTerm;
  /** Reference to the force field being used. */
  private final ForceField forceField;
  private static final double DEFAULT_POLAR_EPS = 1.0e-6;
  @FFXKeyword(name = "polar-eps", keywordGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0e-6",
      description =
          "Sets the convergence criterion applied during computation of self-consistent induced dipoles. "
              + "The calculation is deemed to have converged when the rms change in Debyes in "
              + "the induced dipole at all polarizable sites is less than the value specified with this property. "
              + "The default value in the absence of the keyword is 0.000001 Debyes")
  private final double poleps;
  private final boolean directFallback;
  /** Specify an SCF predictor algorithm. */
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
  /**
   * Polarization modes include "direct", in which induced dipoles do not interact, and "mutual" that
   * converges the self-consistent field to a tolerance specified by the "polar-eps" keyword.
   */
  public Polarization polarization;
  /** Dimensions of [nsymm][xyz][nAtoms]. */
  public double[][][] coordinates;
  /**
   * Neighbor lists, including atoms beyond the real space cutoff. [nsymm][nAtoms][nAllNeighbors]
   */
  public int[][][] neighborLists;
  /** Cartesian multipoles in the global frame with dimensions of [nsymm][nAtoms][10] */
  public double[][][] globalMultipole;
  /** Fractional multipoles in the global frame with dimensions of [nsymm][nAtoms][10] */
  public double[][][] fractionalMultipole;
  /** Dimensions of [nsymm][nAtoms][3] */
  public double[][][] inducedDipole;
  public double[][][] inducedDipoleCR;
  /** Direct induced dipoles. */
  public double[][] directDipole;
  public double[][] directDipoleCR;
  public double[][] directField;
  public double[][] directFieldCR;
  /** Vacuum induced dipoles */
  public double[][][] vacuumInducedDipole;
  public double[][][] vacuumInducedDipoleCR;
  /** Vacuum induced dipoles */
  public double[][] vacuumDirectDipole;
  public double[][] vacuumDirectDipoleCR;
  /**
   * Log the seven components of total electrostatic energy at each evaluation: (Permanent)
   * PermanentRealSpace, PermanentSelf, PermanentRecip (Induced) InducedRealSpace, InducedSelf,
   * InducedRecip, and GeneralizedKirkwood. Self, Recip terms apply only to periodic systems; GK
   * applies only when requested and aperiodic.
   */
  public boolean printDecomposition = false;
  /**
   * Disables windowed lambda ranges by setting permLambdaStart = polLambdaStart = 0.0 and
   * permLambdaEnd = polLambdaEnd = 1.0.
   */
  public boolean noWindowing;

  /** Coulomb constant in units of kcal*Ang/(mol*electron^2) */
  @FFXKeyword(name = "electric", keywordGroup = LocalGeometryFunctionalForm, defaultValue = "332.063713",
      description =
          "Specifies a value for the so-called \"electric constant\" allowing conversion unit of electrostatic potential energy values from electrons^2/Angstrom to kcal/mol. "
              + "Internally, FFX stores a default value for this constant as 332.063713 based on CODATA reference values. "
              + "Since different force fields are intended for use with slightly different values, this keyword allows overriding the default value.")
  public double electric;

  /** An ordered array of atoms in the system. */
  protected Atom[] atoms;
  /** The number of atoms in the system. */
  protected int nAtoms;
  /** Polarization groups. */
  protected int[][] ip11;
  protected int[][] ip12;
  protected int[][] ip13;
  /**
   * Total multipole energy = permanentMultipoleEnergy + polarizationEnergy. <br> This does not
   * include GK.
   */
  protected double totalMultipoleEnergy;
  /**
   * Permanent multipole energy = permanentRealSpaceEnergy + permanentSelfEnergy +
   * permanentReciprocalEnergy.
   */
  protected double permanentMultipoleEnergy;
  protected double permanentRealSpaceEnergy;
  protected double permanentSelfEnergy;
  protected double permanentReciprocalEnergy;
  /** Polarization energy = inducedRealSpaceEnergy + inducedSelfEnergy + inducedReciprocalEnergy. */
  protected double polarizationEnergy;
  protected double inducedRealSpaceEnergy;
  protected double inducedSelfEnergy;
  protected double inducedReciprocalEnergy;
  protected SCFAlgorithm scfAlgorithm;
  /**
   * Log the induced dipole magnitudes and directions. Use the cgo_arrow.py script (available from
   * the wiki) to draw these easily in PyMol.
   */
  protected boolean printInducedDipoles;

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
  private double[][] fracMultipolePhi;
  private double[][] cartesianInducedDipolePhi;
  private double[][] cartesianInducedDipolePhiCR;
  private double[][] fractionalInducedDipolePhi;
  private double[][] fractionalInducedDipolePhiCR;
  private double[][] cartesianVacuumDipolePhi;
  private double[][] cartesianVacuumDipolePhiCR;
  private double[][] fractionalVacuumDipolePhi;
  private double[][] fractionalVacuumDipolePhiCR;
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
   * @param ewaldCutoff The Ewald real space cutoff.
   * @param gkCutoff The generalized Kirkwood cutoff.
   * @param parallelTeam A ParallelTeam that delegates parallelization.
   */
  public ParticleMeshEwald(
      Atom[] atoms,
      int[] molecule,
      ForceField forceField,
      Crystal crystal,
      NeighborList neighborList,
      ELEC_FORM elecForm,
      double ewaldCutoff,
      double gkCutoff,
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
    poleps = forceField.getDouble("POLAR_EPS", DEFAULT_POLAR_EPS);

    // If PME-specific lambda term not set, default to force field-wide lambda term.
    lambdaTerm =
        forceField.getBoolean("ELEC_LAMBDATERM", forceField.getBoolean("LAMBDATERM", false));

    CompositeConfiguration properties = forceField.getProperties();
    printInducedDipoles = properties.getBoolean("pme.printInducedDipoles", false);
    noWindowing = properties.getBoolean("pme.noWindowing", false);

    double aewald = forceField.getDouble("EWALD_ALPHA", DEFAULT_EWALD_COEFFICIENT);
    ewaldParameters = new EwaldParameters(ewaldCutoff, aewald);
    scaleParameters = new ScaleParameters(elecForm, forceField);
    reciprocalSpaceTerm = forceField.getBoolean("RECIPTERM", true);

    SCFPredictor scfPredictor;
    try {
      String predictor = forceField.getString("SCF_PREDICTOR", "NONE");
      predictor = predictor.replaceAll("-", "_").toUpperCase();
      scfPredictor = SCFPredictor.valueOf(predictor);
    } catch (Exception e) {
      scfPredictor = SCFPredictor.NONE;
    }

    scfPredictorParameters = new SCFPredictorParameters(scfPredictor, nAtoms);
    if (scfPredictor != SCFPredictor.NONE) {
      scfPredictorParameters.init(forceField);
    }

    String algorithm = forceField.getString("SCF_ALGORITHM", "CG");
    try {
      algorithm = algorithm.replaceAll("-", "_").toUpperCase();
      scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
    } catch (Exception e) {
      scfAlgorithm = SCFAlgorithm.CG;
    }

    // Fall back to the direct
    directFallback = forceField.getBoolean("DIRECT_SCF_FALLBACK", true);

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

    alchemicalParameters = new AlchemicalParameters(forceField, lambdaTerm, noWindowing,
        polarization);

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
        sb.append(format("    SCF Predictor:                     %8s\n",
            scfPredictorParameters.scfPredictor));
        sb.append(format("    SCF Algorithm:                     %8s\n", scfAlgorithm));
        if (scfAlgorithm == SCFAlgorithm.SOR) {
          sb.append(format("    SOR Parameter:                     %8.3f\n", sorRegion.getSOR()));
        } else {
          sb.append(format("    CG Preconditioner Cut-Off:         %8.3f\n",
              pcgSolver.getPreconditionerCutoff()));
          sb.append(format("    CG Preconditioner Ewald Coeff.:    %8.3f\n",
              pcgSolver.getPreconditionerEwald()));
          sb.append(format("    CG Preconditioner Scale:           %8.3f\n",
              pcgSolver.getPreconditionerScale()));
          sb.append(
              format("    CG Preconditioner Mode:     %15s\n", pcgSolver.getPreconditionerMode()));
        }
      }
      if (ewaldParameters.aewald > 0.0) {
        sb.append("   Particle-mesh Ewald\n");
        sb.append(format("    Ewald Coefficient:                 %8.3f\n", ewaldParameters.aewald));
        sb.append(format("    Particle Cut-Off:                  %8.3f (A)", ewaldParameters.off));
      } else if (ewaldParameters.off != Double.POSITIVE_INFINITY) {
        sb.append(format("    Electrostatics Cut-Off:            %8.3f (A)\n",
            ewaldParameters.off));
      } else {
        sb.append("    Electrostatics Cut-Off:                NONE\n");
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

    realSpaceNeighborParameters = new RealSpaceNeighborParameters(maxThreads);
    initializationRegion = new InitializationRegion(this, maxThreads, forceField);
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
      reciprocalEnergyRegion = new ReciprocalEnergyRegion(maxThreads, ewaldParameters.aewald,
          electric);
    } else {
      reciprocalSpace = null;
      reciprocalEnergyRegion = null;
    }
    permanentFieldRegion = new PermanentFieldRegion(realSpaceTeam, forceField, lambdaTerm);
    inducedDipoleFieldRegion = new InducedDipoleFieldRegion(realSpaceTeam, forceField, lambdaTerm);
    inducedDipoleFieldReduceRegion = new InducedDipoleFieldReduceRegion(maxThreads);
    polarizationEnergyRegion = new PolarizationEnergyRegion(maxThreads, electric);
    realSpaceEnergyRegion = new RealSpaceEnergyRegion(maxThreads, forceField, lambdaTerm, electric);
    reduceRegion = new ReduceRegion(maxThreads, forceField);

    pmeTimings = new PMETimings(maxThreads);

    if (lambdaTerm) {
      logger.info(alchemicalParameters.toString());
    }

    // The GK reaction field is added to the intra-molecular field to give the self-consistent
    // reaction field.
    generalizedKirkwoodTerm = forceField.getBoolean("GKTERM", false);
    if (generalizedKirkwoodTerm || alchemicalParameters.doLigandGKElec) {
      generalizedKirkwood = new GeneralizedKirkwood(forceField, atoms,
          this, crystal, parallelTeam, electric, gkCutoff);
    } else {
      generalizedKirkwood = null;
    }
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
      reciprocalSpace.computeInducedPhi(
          cartesianInducedDipolePhi, cartesianInducedDipolePhiCR,
          fractionalInducedDipolePhi, fractionalInducedDipolePhiCR);
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
        cartesianInducedDipolePhi,
        cartesianInducedDipolePhiCR,
        field,
        fieldCR);
    inducedDipoleFieldReduceRegion.executeWith(parallelTeam);
  }

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

    if (esvTerm) {
      extendedSystem.initEsvPermElec();
      extendedSystem.initEsvIndElec();
    }

    alchemicalParameters.doPermanentRealSpace = true;
    alchemicalParameters.permanentScale = 1.0;
    alchemicalParameters.doPolarization = true;
    alchemicalParameters.polarizationScale = 1.0;

    // Expand coordinates and rotate multipoles into the global frame.
    initializationRegion.init(
        lambdaTerm,
        extendedSystem,
        atoms,
        coordinates,
        crystal,
        frame,
        axisAtom,
        globalMultipole,
        dMultipoledTirationESV,
        dMultipoledTautomerESV,
        polarizability,
        dPolardTitrationESV,
        dPolardTautomerESV,
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
      pmeTimings.printRealSpaceTimings(maxThreads, realSpaceEnergyRegion);
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

  public int[][] getAxisAtoms() {
    return axisAtom;
  }

  public double getCavitationEnergy() {
    return generalizedKirkwood.getCavitationEnergy();
  }

  public double[][][] getCoordinates() {
    return coordinates;
  }

  public double getDispersionEnergy() {
    return generalizedKirkwood.getDispersionEnergy();
  }

  public double getEwaldCoefficient() {
    return ewaldParameters.aewald;
  }

  public double getEwaldCutoff() {
    return ewaldParameters.off;
  }

  public GeneralizedKirkwood getGK() {
    return generalizedKirkwood;
  }

  /**
   * getGeneralizedKirkwoodEnergy.
   *
   * @return a double.
   */
  public double getGKEnergy() {
    return generalizedKirkwood.getGeneralizedKirkwoordEnergy();
  }

  /**
   * getGKInteractions
   *
   * @return a int.
   */
  public int getGKInteractions() {
    return gkInteractions;
  }

  /**
   * Getter for the field <code>interactions</code>.
   *
   * @return a int.
   */
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

  public double getPermanentEnergy() {
    return permanentMultipoleEnergy;
  }

  /**
   * getIndRealEnergy.
   *
   * @return a double.
   */
  public double getIndRealEnergy() {
    return inducedRealSpaceEnergy;
  }

  /**
   * getIndRecipEnergy.
   *
   * @return a double.
   */
  public double getIndRecipEnergy() {
    return inducedReciprocalEnergy;
  }

  /**
   * getIndSelfEnergy.
   *
   * @return a double.
   */
  public double getIndSelfEnergy() {
    return inducedSelfEnergy;
  }

  /**
   * getPermSelfEnergy.
   *
   * @return a double.
   */
  public double getPermSelfEnergy() {
    return permanentSelfEnergy;
  }

  public double getPermRealEnergy() {
    return permanentRealSpaceEnergy;
  }

  public double getPermRecipEnergy() {
    return permanentReciprocalEnergy;
  }

  public double getPolarEps() {
    return poleps;
  }

  public int[][] getPolarization11() {
    return ip11;
  }

  public int[][] getPolarization12() {
    return ip12;
  }

  public int[][] getPolarization13() {
    return ip13;
  }

  /**
   * Getter for the field <code>polarizationEnergy</code>.
   *
   * @return a double.
   */
  public double getPolarizationEnergy() {
    return polarizationEnergy;
  }

  public Polarization getPolarizationType() {
    return polarization;
  }

  public ReciprocalSpace getReciprocalSpace() {
    return reciprocalSpace;
  }

  public double getScale14() {
    return scaleParameters.m14scale;
  }

  /**
   * getGKEnergy
   *
   * @return a double.
   */
  public double getSolvationEnergy() {
    return solvationEnergy;
  }

  /**
   * Get the MultipoleType for Atom i.
   *
   * @param i The atom index.
   * @return The MultipoleType.
   */
  public MultipoleType getMultipoleType(int i) {
    Atom atom = atoms[i];
    MultipoleType multipoleType = atom.getMultipoleType();
    double[] multipole = multipoleType.getMultipole();
    if (esvTerm && extendedSystem.isTitrating(i)) {
      double[] esvMultipole = new double[10];
      System.arraycopy(multipole, 0, esvMultipole, 0, multipole.length);
      double titrationLambda = extendedSystem.getTitrationLambda(i);
      double tautomerLambda = extendedSystem.getTautomerLambda(i);
      esvMultipole = extendedSystem.getTitrationUtils()
          .getMultipole(atom, titrationLambda, tautomerLambda, esvMultipole);
      // Create a new MultipoleType for the tritrating site.
      multipoleType = new MultipoleType(multipoleType, esvMultipole);
    }
    return multipoleType;
  }

  /**
   * Get the PolarizeType for Atom i.
   *
   * @param i The atom index.
   * @return The PolarizeType.
   */
  public PolarizeType getPolarizeType(int i) {
    Atom atom = atoms[i];
    PolarizeType polarizeType = atom.getPolarizeType();
    if (polarizeType != null) {
      if (esvTerm && extendedSystem.isTitrating(i) && (extendedSystem.isTitratingHydrogen(i)
          || extendedSystem.isTitratingSulfur(i))) {
        double titrationLambda = extendedSystem.getTitrationLambda(i);
        double tautomerLambda = extendedSystem.getTautomerLambda(i);
        double esvPolarizability = extendedSystem.getTitrationUtils()
            .getPolarizability(atom, titrationLambda, tautomerLambda, polarizeType.polarizability);
        polarizeType = new PolarizeType(polarizeType, esvPolarizability);
      }
    }
    return polarizeType;
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

  public void setCrystal(Crystal crystal) {
    // Check if memory allocation is required.
    int nSymmNew = crystal.spaceGroup.getNumberOfSymOps();
    if (nSymm < nSymmNew) {
      coordinates = new double[nSymmNew][3][nAtoms];
      globalMultipole = new double[nSymmNew][nAtoms][10];
      fractionalMultipole = new double[nSymmNew][nAtoms][10];
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

  private void initAtomArrays() {
    if (localMultipole == null || localMultipole.length < nAtoms) {
      localMultipole = new double[nAtoms][10];
      frame = new MultipoleFrameDefinition[nAtoms];
      axisAtom = new int[nAtoms][];
      cartesianMultipolePhi = new double[nAtoms][tensorCount];
      fracMultipolePhi = new double[nAtoms][tensorCount];
      directDipole = new double[nAtoms][3];
      directDipoleCR = new double[nAtoms][3];
      directField = new double[nAtoms][3];
      directFieldCR = new double[nAtoms][3];
      vacuumDirectDipole = new double[nAtoms][3];
      vacuumDirectDipoleCR = new double[nAtoms][3];
      if (optRegion != null) {
        int optOrder = optRegion.optOrder;
        optRegion.optDipole = new double[optOrder + 1][nAtoms][3];
        optRegion.optDipoleCR = new double[optOrder + 1][nAtoms][3];
      }
      // Total Induced Dipole Phi (i.e. with GK)
      cartesianInducedDipolePhi = new double[nAtoms][tensorCount];
      cartesianInducedDipolePhiCR = new double[nAtoms][tensorCount];
      fractionalInducedDipolePhi = new double[nAtoms][tensorCount];
      fractionalInducedDipolePhiCR = new double[nAtoms][tensorCount];
      // Vacuum Induced Dipole Phi (i.e. no GK)
      cartesianVacuumDipolePhi = new double[nAtoms][tensorCount];
      cartesianVacuumDipolePhiCR = new double[nAtoms][tensorCount];
      fractionalVacuumDipolePhi = new double[nAtoms][tensorCount];
      fractionalVacuumDipolePhiCR = new double[nAtoms][tensorCount];

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

      if (scfPredictorParameters.scfPredictor != SCFPredictor.NONE) {
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
      fractionalMultipole = new double[nSymm][nAtoms][10];
      inducedDipole = new double[nSymm][nAtoms][3];
      inducedDipoleCR = new double[nSymm][nAtoms][3];
      vacuumInducedDipole = new double[nSymm][nAtoms][3];
      vacuumInducedDipoleCR = new double[nSymm][nAtoms][3];

      // The size of reduced neighbor list depends on the size of the real space cutoff.
      realSpaceNeighborParameters.allocate(nAtoms, nSymm);

      // Lambda factors are different for OST and ESV interactions.
      lambdaFactors = new LambdaFactors[maxThreads];
      for (int i = 0; i < maxThreads; i++) {
        /*if (esvTerm) {
          // Invoked every time through inner loops.
          lambdaFactors[i] = new LambdaFactorsESV();
        } else */
        if (lambdaTerm) {
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
    if (alchemicalParameters.doLigandVaporElec
        && alchemicalParameters.vaporCrystal == null) {
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
              null, alchemicalParameters.vaporCrystal, atoms, vacuumOff, 2.0,
              parallelTeam);
      vacuumNeighborList.setIntermolecular(false, molecule);

      alchemicalParameters.vaporLists = new int[1][nAtoms][];
      double[][] coords = new double[1][nAtoms * 3];
      for (int i = 0; i < nAtoms; i++) {
        coords[0][i * 3] = atoms[i].getX();
        coords[0][i * 3 + 1] = atoms[i].getY();
        coords[0][i * 3 + 2] = atoms[i].getZ();
      }
      boolean print = logger.isLoggable(Level.FINE);
      vacuumNeighborList.buildList(coords, alchemicalParameters.vaporLists, isSoft,
          true, print);
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
    alchemicalParameters.permanentScale =
        1.0 - alchemicalParameters.lPowPerm;
    alchemicalParameters.dEdLSign = -1.0;

    // If we are past the end of the polarization lambda window, then only the condensed phase is
    // necessary.
    if (lambda <= alchemicalParameters.polLambdaEnd
        && alchemicalParameters.doNoLigandCondensedSCF) {
      alchemicalParameters.doPolarization = true;
      alchemicalParameters.polarizationScale =
          1.0 - alchemicalParameters.lPowPol;
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
      Arrays.fill(use, true);
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
    alchemicalParameters.permanentScale =
        1.0 - alchemicalParameters.lPowPerm;
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
      alchemicalParameters.polarizationScale =
          1.0 - alchemicalParameters.lPowPol;
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
          alchemicalParameters.dEdLSign
              * alchemicalParameters.d2lPowPol);
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

      for (int i = 0; i < nAtoms; i++) {
        if (polarization != Polarization.NONE && esvTerm && extendedSystem.isTitrating(i)
            && (extendedSystem.isTitratingHydrogen(i) || extendedSystem.isTitratingSulfur(i))) {
          double dx = field.getX(i);
          double dy = field.getY(i);
          double dz = field.getZ(i);
          double dxCR = fieldCR.getX(i);
          double dyCR = fieldCR.getY(i);
          double dzCR = fieldCR.getZ(i);
          //Add back permanent multipole field to total field for extended system derivatives if mutual polarization is used
          if (polarization == Polarization.MUTUAL) {
            dx += directField[i][0];
            dy += directField[i][1];
            dz += directField[i][2];
            dxCR += directFieldCR[i][0];
            dyCR += directFieldCR[i][1];
            dzCR += directFieldCR[i][2];
          }
          double fix = dx * dPolardTitrationESV[i] * dxCR;
          double fiy = dy * dPolardTitrationESV[i] * dyCR;
          double fiz = dz * dPolardTitrationESV[i] * dzCR;
          double titrdUdL = fix + fiy + fiz;
          double tautdUdL = 0.0;
          if (extendedSystem.isTautomerizing(i)) {
            fix = dx * dPolardTautomerESV[i] * dxCR;
            fiy = dy * dPolardTautomerESV[i] * dyCR;
            fiz = dz * dPolardTautomerESV[i] * dzCR;
            tautdUdL = fix + fiy + fiz;
          }
          //logger.info(format("Index i: %d Polarizability: %6.8f TitrDeriv: %6.8f TautDeriv: %6.8f",
          // i, polarizability[i], dPolardTitrationESV[i], dPolardTautomerESV[i]));
          extendedSystem.addIndElecDeriv(i, titrdUdL * electric * -0.5, tautdUdL * electric * -0.5);
        }
      }
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
          reciprocalSpace.computeInducedPhi(
              cartesianInducedDipolePhi, cartesianInducedDipolePhiCR,
              fractionalInducedDipolePhi, fractionalInducedDipolePhiCR);
        }
      }

      if (scfPredictorParameters.scfPredictor != SCFPredictor.NONE) {
        scfPredictorParameters.saveMutualInducedDipoles(lambdaMode,
            inducedDipole, inducedDipoleCR, directDipole, directDipoleCR);
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

    // With GK, we need to compute the polarization energy with vacuum induced dipoles.
    double[][][] inputDipole = inducedDipole;
    double[][][] inputDipoleCR = inducedDipoleCR;
    double[][] inputPhi = cartesianInducedDipolePhi;
    double[][] inputPhiCR = cartesianInducedDipolePhiCR;
    double[][] fracInputPhi = fractionalInducedDipolePhi;
    double[][] fracInputPhiCR = fractionalInducedDipolePhiCR;
    if (generalizedKirkwoodTerm) {
      inputDipole = vacuumInducedDipole;
      inputDipoleCR = vacuumInducedDipoleCR;
      inputPhi = cartesianVacuumDipolePhi;
      inputPhiCR = cartesianVacuumDipolePhiCR;
      fracInputPhi = fractionalVacuumDipolePhi;
      fracInputPhiCR = fractionalVacuumDipolePhiCR;
    }

    double eself = 0.0;
    double erecip = 0.0;
    double eselfi = 0.0;
    double erecipi = 0.0;
    polarizationEnergyRegion.setPolarizationEnergy(0.0);
    if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
      reciprocalEnergyRegion.init(atoms, crystal, gradient, lambdaTerm, esvTerm, use,
          globalMultipole, fractionalMultipole, dMultipoledTirationESV, dMultipoledTautomerESV,
          cartesianMultipolePhi, fracMultipolePhi,
          polarization, inputDipole, inputDipoleCR, inputPhi, inputPhiCR, fracInputPhi,
          fracInputPhiCR,
          reciprocalSpace, alchemicalParameters, extendedSystem,
          // Output
          grad, torque, lambdaGrad, lambdaTorque, shareddEdLambda, sharedd2EdLambda2);
      reciprocalEnergyRegion.executeWith(parallelTeam);
      eself = reciprocalEnergyRegion.getPermanentSelfEnergy();
      erecip = reciprocalEnergyRegion.getPermanentReciprocalEnergy();
      eselfi = reciprocalEnergyRegion.getInducedDipoleSelfEnergy();
      erecipi = reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
      interactions += nAtoms;
    }

    pmeTimings.realSpaceEnergyTotal -= System.nanoTime();
    realSpaceEnergyRegion.init(atoms, crystal, extendedSystem, esvTerm, coordinates, frame, axisAtom,
        globalMultipole, dMultipoledTirationESV, dMultipoledTautomerESV,
        inputDipole, inputDipoleCR, use, molecule,
        ip11, mask12, mask13, mask14, mask15, isSoft, ipdamp, thole, realSpaceNeighborParameters,
        gradient, lambdaTerm, lambdaMode, polarization,
        ewaldParameters, scaleParameters, alchemicalParameters,
        pmeTimings.realSpaceEnergyTime,
        // Output
        grad, torque, lambdaGrad, lambdaTorque, shareddEdLambda, sharedd2EdLambda2);
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
      AlchemicalParameters alchemicalParametersGK = new AlchemicalParameters(
          forceField, false, noWindowing, polarization);
      // No permanent multipole contribution.
      alchemicalParametersGK.permanentScale = 0.0;
      alchemicalParametersGK.doPermanentRealSpace = false;
      // Flip the sign on the vacuum polarization energy and derivatives.
      alchemicalParametersGK.polarizationScale = -1.0;

      // Store the derivative and torque contributions in the GK arrays.
      AtomicDoubleArray3D gradGK = generalizedKirkwood.getGrad();
      AtomicDoubleArray3D torqueGK = generalizedKirkwood.getTorque();

      if (reciprocalSpaceTerm && ewaldParameters.aewald > 0.0) {
        reciprocalEnergyRegion.init(
            atoms,
            crystal,
            gradient,
            false,
            esvTerm,
            use,
            globalMultipole,
            fractionalMultipole,
            dMultipoledTirationESV,
            dMultipoledTautomerESV,
            cartesianMultipolePhi,
            fracMultipolePhi,
            polarization,
            vacuumInducedDipole,
            vacuumInducedDipoleCR,
            cartesianVacuumDipolePhi,
            cartesianVacuumDipolePhiCR,
            fractionalVacuumDipolePhi,
            fractionalVacuumDipolePhiCR,
            reciprocalSpace,
            alchemicalParametersGK,
            extendedSystem,
            gradGK,
            torqueGK,
            null,
            null,
            shareddEdLambda,
            sharedd2EdLambda2);
        reciprocalEnergyRegion.executeWith(parallelTeam);
        eGK = reciprocalEnergyRegion.getInducedDipoleSelfEnergy()
            + reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
      }

      pmeTimings.gkEnergyTotal -= System.nanoTime();
      realSpaceEnergyRegion.init(
          atoms,
          crystal,
          extendedSystem,
          esvTerm,
          coordinates,
          frame,
          axisAtom,
          globalMultipole,
          dMultipoledTirationESV,
          dMultipoledTautomerESV,
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
        reciprocalEnergyRegion.init(
            atoms,
            crystal,
            gradient,
            false,
            esvTerm,
            use,
            globalMultipole,
            fractionalMultipole,
            dMultipoledTirationESV,
            dMultipoledTautomerESV,
            cartesianMultipolePhi,
            fracMultipolePhi,
            polarization,
            inducedDipole,
            inducedDipoleCR,
            cartesianInducedDipolePhi,
            cartesianInducedDipolePhiCR,
            fractionalInducedDipolePhi,
            fractionalInducedDipolePhiCR,
            reciprocalSpace,
            alchemicalParametersGK,
            extendedSystem,
            gradGK,
            torqueGK,
            null,
            null,
            shareddEdLambda,
            sharedd2EdLambda2);
        reciprocalEnergyRegion.executeWith(parallelTeam);
        eGK += reciprocalEnergyRegion.getInducedDipoleSelfEnergy()
            + reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
      }

      pmeTimings.gkEnergyTotal -= System.nanoTime();
      realSpaceEnergyRegion.init(
          atoms,
          crystal,
          extendedSystem,
          esvTerm,
          coordinates,
          frame,
          axisAtom,
          globalMultipole,
          dMultipoledTirationESV,
          dMultipoledTautomerESV,
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
        reciprocalSpace.splinePermanentMultipoles(globalMultipole, fractionalMultipole, use);
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
        reciprocalSpace.computePermanentPhi(cartesianMultipolePhi, fracMultipolePhi);
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
      System.arraycopy(cartesianInducedDipolePhi[i], 0, cartesianVacuumDipolePhi[i], 0, tensorCount);
      System.arraycopy(cartesianInducedDipolePhiCR[i], 0, cartesianVacuumDipolePhiCR[i], 0,
          tensorCount);
      System.arraycopy(fractionalInducedDipolePhi[i], 0, fractionalVacuumDipolePhi[i], 0,
          tensorCount);
      System.arraycopy(fractionalInducedDipolePhiCR[i], 0, fractionalVacuumDipolePhiCR[i], 0,
          tensorCount);
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
        directDipoleCR,
        directField,
        directFieldCR
    );
    directRegion.executeWith(parallelTeam);

    // Return unless mutual polarization is selected.
    if (polarization != Polarization.MUTUAL) {
      expandInducedDipoles();
      return 0;
    }

    // Predict the current self-consistent induced dipoles using information from previous steps.
    if (scfPredictorParameters.scfPredictor != SCFPredictor.NONE) {
      switch (scfPredictorParameters.scfPredictor) {
        case ASPC:
          scfPredictorParameters.aspcPredictor(lambdaMode, inducedDipole, inducedDipoleCR);
          break;
        case LS:
          scfPredictorParameters.leastSquaresPredictor(lambdaMode, inducedDipole, inducedDipoleCR);
          break;
        case POLY:
          scfPredictorParameters.polynomialPredictor(lambdaMode, inducedDipole, inducedDipoleCR);
          break;
        default:
          break;
      }
    }

    // Expand the initial induced dipoles to P1 symmetry, if necessary.
    expandInducedDipoles();

    // Converge the self-consistent field.
    try {
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
    } catch (EnergyException ex) {
      if (directFallback) {
        // SCF Failure: warn and revert to direct polarization.
        logger.warning(ex.toString());
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
            directDipoleCR,
            directField,
            directFieldCR
        );
        directRegion.executeWith(parallelTeam);
        expandInducedDipoles();
        logger.info(" Direct induced dipoles computed due to SCF failure.");
        return 0;
      } else {
        throw ex;
      }
    }
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
          reciprocalSpace.computeInducedPhi(
              cartesianInducedDipolePhi, cartesianInducedDipolePhiCR,
              fractionalInducedDipolePhi, fractionalInducedDipolePhiCR);
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
            cartesianInducedDipolePhi,
            cartesianInducedDipolePhiCR,
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
        throw new EnergyException(message);
      }

      // The SCF should converge well before the max iteration check. Otherwise, fail the SCF
      // process.
      if (completedSCFCycles >= maxSCFCycles) {
        if (sb != null) {
          logger.warning(sb.toString());
        }
        String message = format(" Maximum SCF iterations reached: (%d)\n", completedSCFCycles);
        throw new EnergyException(message);
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
          reciprocalSpace.computeInducedPhi(
              cartesianInducedDipolePhi, cartesianInducedDipolePhiCR,
              fractionalInducedDipolePhi, fractionalInducedDipolePhiCR);
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
            cartesianInducedDipolePhi,
            cartesianInducedDipolePhiCR,
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
   * Getter for the field <code>totalMultipoleEnergy</code>.
   *
   * @return a double.
   */
  public double getTotalMultipoleEnergy() {
    return permanentMultipoleEnergy + polarizationEnergy;
  }

  /**
   * Setter for the field <code>polarization</code>.
   *
   * @param set a {@link Polarization} object.
   */
  public void setPolarization(Polarization set) {
    this.polarization = set;
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
        for (Bond bond12 : bonds) {
          Atom a2 = bond12.get1_2(atom);
          AtomType aType2 = a2.getAtomType();
          sb.append(format("\n  1-2 %s --> %s", a2, aType2));
        }
        for (Bond bond12 : bonds) {
          Atom atom2 = bond12.get1_2(atom);
          bonds = atom2.getBonds();
          for (Bond bond23 : bonds) {
            Atom a2 = bond23.get1_2(atom2);
            AtomType aType2 = a2.getAtomType();
            sb.append(format("\n  1-3 %s --> %s", a2, aType2));
          }
        }

        List<MultipoleType> multipoleTypes = forceField.getMultipoleTypes(
            "" + atom.getAtomType().getKey());
        if (multipoleTypes != null && !multipoleTypes.isEmpty()) {
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
   * AssignPolarizationGroups.
   */
  protected void assignPolarizationGroups() {
    // Find directly connected group members for each atom.
    List<Integer> group = new ArrayList<>();
    for (int i = 0; i < nAtoms; i++) {
      Atom a = atoms[i];
      if (a.getIndex() - 1 != i) {
        logger.info(format(" PME Index i: %d, %s Index: %d\n Atom: %s",
            i, atomIndexing, a.getIndex(), a));
        logger.severe(" Atom indexing is not consistent in PME.");
      }
    }
    for (Atom ai : atoms) {
      group.clear();
      int index = ai.getIndex() - 1;
      group.add(index);
      PolarizeType polarizeType = ai.getPolarizeType();
      if (polarizeType != null) {
        if (polarizeType.polarizationGroup != null) {
          PolarizeType.growGroup(group, ai);
          sort(group);
          ip11[index] = new int[group.size()];
          int j = 0;
          for (int k : group) {
            ip11[index][j++] = k;
          }
        } else {
          ip11[index] = new int[group.size()];
          int j = 0;
          for (int k : group) {
            ip11[index][j++] = k;
          }
        }
      } else {
        String message = "The polarize keyword was not found for atom "
            + (index + 1) + " with type " + ai.getType();
        logger.severe(message);
      }
    }
    // Find 1-2 group relationships.
    int[] mask = new int[nAtoms];
    List<Integer> list = new ArrayList<>();
    List<Integer> keep = new ArrayList<>();
    for (int i = 0; i < nAtoms; i++) {
      mask[i] = -1;
    }
    for (int i = 0; i < nAtoms; i++) {
      list.clear();
      for (int j : ip11[i]) {
        list.add(j);
        mask[j] = i;
      }
      keep.clear();
      for (int j : list) {
        Atom aj = atoms[j];
        List<Bond> bonds = aj.getBonds();
        for (Bond b : bonds) {
          Atom ak = b.get1_2(aj);
          int k = ak.getIndex() - 1;
          if (mask[k] != i) {
            keep.add(k);
          }
        }
      }
      list.clear();
      for (int j : keep) {
        for (int k : ip11[j]) {
          list.add(k);
        }
      }
      sort(list);
      ip12[i] = new int[list.size()];
      int j = 0;
      for (int k : list) {
        ip12[i][j++] = k;
      }
    }

    // Find 1-3 group relationships.
    for (int i = 0; i < nAtoms; i++) {
      mask[i] = -1;
    }
    for (int i = 0; i < nAtoms; i++) {
      for (int j : ip11[i]) {
        mask[j] = i;
      }
      for (int j : ip12[i]) {
        mask[j] = i;
      }
      list.clear();
      for (int j : ip12[i]) {
        for (int k : ip12[j]) {
          if (mask[k] != i) {
            if (!list.contains(k)) {
              list.add(k);
            }
          }
        }
      }
      ip13[i] = new int[list.size()];
      sort(list);
      int j = 0;
      for (int k : list) {
        ip13[i][j++] = k;
      }
    }
  }

  /**
   * Compute multipole moments for an array of atoms.
   *
   * @param activeAtoms Atom array to consider.
   * @param forceEnergy Force calculation of the electrostatic energy (rotate multipoles, perform
   *     SCF).
   */
  public void computeMoments(Atom[] activeAtoms, boolean forceEnergy) {
    // Zero out total charge, dipole and quadrupole components.
    var netchg = 0.0;
    var netdpl = 0.0;
    var xdpl = 0.0;
    var ydpl = 0.0;
    var zdpl = 0.0;
    var xxqdp = 0.0;
    var xyqdp = 0.0;
    var xzqdp = 0.0;
    var yxqdp = 0.0;
    var yyqdp = 0.0;
    var yzqdp = 0.0;
    var zxqdp = 0.0;
    var zyqdp = 0.0;
    var zzqdp = 0.0;

    // Find the center of mass of the set of active atoms.
    double xmid = 0.0;
    double ymid = 0.0;
    double zmid = 0.0;
    double totalMass = 0;
    for (Atom atom : activeAtoms) {
      var m = atom.getMass();
      totalMass += m;
      xmid = xmid + atom.getX() * m;
      ymid = ymid + atom.getY() * m;
      zmid = zmid + atom.getZ() * m;
    }
    if (totalMass > 0) {
      xmid /= totalMass;
      ymid /= totalMass;
      zmid /= totalMass;
    }
    int n = activeAtoms.length;
    double[] xcm = new double[n];
    double[] ycm = new double[n];
    double[] zcm = new double[n];
    int k = 0;
    for (Atom atom : activeAtoms) {
      xcm[k] = atom.getX() - xmid;
      ycm[k] = atom.getY() - ymid;
      zcm[k] = atom.getZ() - zmid;
      k++;
    }

    if (forceEnergy) {
      energy(false, false);
    }

    // Account for charge, dipoles and induced dipoles.
    k = 0;
    for (Atom atom : activeAtoms) {
      int i = atom.getIndex() - 1;
      double[] globalMultipolei = globalMultipole[0][i];
      double[] inducedDipolei = inducedDipole[0][i];

      var ci = globalMultipolei[t000];
      var dix = globalMultipolei[t100];
      var diy = globalMultipolei[t010];
      var diz = globalMultipolei[t001];
      var uix = inducedDipolei[0];
      var uiy = inducedDipolei[1];
      var uiz = inducedDipolei[2];

      netchg += ci;
      xdpl += xcm[k] * ci + dix + uix;
      ydpl += ycm[k] * ci + diy + uiy;
      zdpl += zcm[k] * ci + diz + uiz;
      xxqdp += xcm[k] * xcm[k] * ci + 2.0 * xcm[k] * (dix + uix);
      xyqdp += xcm[k] * ycm[k] * ci + xcm[k] * (diy + uiy) + ycm[k] * (dix + uix);
      xzqdp += xcm[k] * zcm[k] * ci + xcm[k] * (diz + uiz) + zcm[k] * (dix + uix);
      yxqdp += ycm[k] * xcm[k] * ci + ycm[k] * (dix + uix) + xcm[k] * (diy + uiy);
      yyqdp += ycm[k] * ycm[k] * ci + 2.0 * ycm[k] * (diy + uiy);
      yzqdp += ycm[k] * zcm[k] * ci + ycm[k] * (diz + uiz) + zcm[k] * (diy + uiy);
      zxqdp += zcm[k] * xcm[k] * ci + zcm[k] * (dix + uix) + xcm[k] * (diz + uiz);
      zyqdp += zcm[k] * ycm[k] * ci + zcm[k] * (diy + uiy) + ycm[k] * (diz + uiz);
      zzqdp += zcm[k] * zcm[k] * ci + 2.0 * zcm[k] * (diz + uiz);
      k++;
    }

    // Convert the quadrupole from traced to traceless form.
    var qave = (xxqdp + yyqdp + zzqdp) / 3.0;
    xxqdp = 1.5 * (xxqdp - qave);
    xyqdp = 1.5 * xyqdp;
    xzqdp = 1.5 * xzqdp;
    yxqdp = 1.5 * yxqdp;
    yyqdp = 1.5 * (yyqdp - qave);
    yzqdp = 1.5 * yzqdp;
    zxqdp = 1.5 * zxqdp;
    zyqdp = 1.5 * zyqdp;
    zzqdp = 1.5 * (zzqdp - qave);

    // Add the traceless atomic quadrupoles to total quadrupole.
    for (Atom atom : activeAtoms) {
      int i = atom.getIndex() - 1;
      double[] globalMultipolei = globalMultipole[0][i];
      var qixx = globalMultipolei[t200];
      var qiyy = globalMultipolei[t020];
      var qizz = globalMultipolei[t002];
      var qixy = globalMultipolei[t110];
      var qixz = globalMultipolei[t101];
      var qiyz = globalMultipolei[t011];
      xxqdp += qixx;
      xyqdp += qixy;
      xzqdp += qixz;
      yxqdp += qixy;
      yyqdp += qiyy;
      yzqdp += qiyz;
      zxqdp += qixz;
      zyqdp += qiyz;
      zzqdp += qizz;
    }

    // Convert dipole to Debye and quadrupole to Buckingham.
    xdpl = xdpl * ELEC_ANG_TO_DEBYE;
    ydpl = ydpl * ELEC_ANG_TO_DEBYE;
    zdpl = zdpl * ELEC_ANG_TO_DEBYE;
    xxqdp = xxqdp * ELEC_ANG_TO_DEBYE;
    xyqdp = xyqdp * ELEC_ANG_TO_DEBYE;
    xzqdp = xzqdp * ELEC_ANG_TO_DEBYE;
    yxqdp = yxqdp * ELEC_ANG_TO_DEBYE;
    yyqdp = yyqdp * ELEC_ANG_TO_DEBYE;
    yzqdp = yzqdp * ELEC_ANG_TO_DEBYE;
    zxqdp = zxqdp * ELEC_ANG_TO_DEBYE;
    zyqdp = zyqdp * ELEC_ANG_TO_DEBYE;
    zzqdp = zzqdp * ELEC_ANG_TO_DEBYE;

    // Get dipole magnitude and diagonalize quadrupole tensor.
    netdpl = sqrt(xdpl * xdpl + ydpl * ydpl + zdpl * zdpl);
    double[][] a = new double[3][3];
    a[0][0] = xxqdp;
    a[0][1] = xyqdp;
    a[0][2] = xzqdp;
    a[1][0] = yxqdp;
    a[1][1] = yyqdp;
    a[1][2] = yzqdp;
    a[2][0] = zxqdp;
    a[2][1] = zyqdp;
    a[2][2] = zzqdp;
    EigenDecomposition e = new EigenDecomposition(new Array2DRowRealMatrix(a));
    // Eigenvalues are returned in descending order, but logged below in ascending order.
    var netqdp = e.getRealEigenvalues();

    logger.info("\n Electric Moments\n");
    logger.info(format("  Total Electric Charge:    %13.5f Electrons\n", netchg));
    logger.info(format("  Dipole Moment Magnitude:  %13.5f Debye\n", netdpl));
    logger.info(format("  Dipole X,Y,Z-Components:  %13.5f %13.5f %13.5f\n", xdpl, ydpl, zdpl));
    logger.info(format("  Quadrupole Moment Tensor: %13.5f %13.5f %13.5f", xxqdp, xyqdp, xzqdp));
    logger.info(format("       (Buckinghams)        %13.5f %13.5f %13.5f", yxqdp, yyqdp, yzqdp));
    logger.info(format("                            %13.5f %13.5f %13.5f\n", zxqdp, zyqdp, zzqdp));
    logger.info(
        format(
            "  Principal Axes Quadrupole %13.5f %13.5f %13.5f\n", netqdp[2], netqdp[1], netqdp[0]));
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
   * Attach system with extended variable such as titrations.
   *
   * @param system a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public void attachExtendedSystem(ExtendedSystem system) {
    // Set object handles.
    esvTerm = true;
    extendedSystem = system;
    numESVs = extendedSystem.getNumberOfVariables();

    // Update atoms and reinitialize arrays for consistency with the ExtendedSystem.
    setAtoms(extendedSystem.getExtendedAtoms(), extendedSystem.getExtendedMolecule());
    // Allocate space for dM/dTitratonESV
    if (dMultipoledTirationESV == null || dMultipoledTirationESV.length != nSymm
        || dMultipoledTirationESV[0].length != nAtoms) {
      dMultipoledTirationESV = new double[nSymm][nAtoms][10];
      dMultipoledTautomerESV = new double[nSymm][nAtoms][10];
    }

    if (dPolardTitrationESV == null || dPolardTitrationESV.length != nAtoms) {
      dPolardTitrationESV = new double[nAtoms];
      dPolardTautomerESV = new double[nAtoms];
    }

    //updateEsvLambda();
    logger.info(format(" Attached extended system (%d variables) to PME.\n", numESVs));
  }

  /**
   * The scalar ESV that is operating on each atom (or 1.0 if the atom is not under ESV control).
   */
  double[] perAtomTitrationESV = null;

  /**
   * The partial derivative of each multipole with respect to its titration/tautomer ESV (or 0.0 if
   * the atom is not under titration ESV control).
   */
  double[][][] dMultipoledTirationESV = null;
  double[][][] dMultipoledTautomerESV = null;

  /**
   * The partial derivative of each polarizability with respect to its titration/tautomer ESV (or 0.0
   * if the atom is not under titration ESV control).
   */
  double[] dPolardTitrationESV = null;
  double[] dPolardTautomerESV = null;
  /**
   * OST and ESV specific factors that effect real space interactions.
   */
  LambdaFactors[] lambdaFactors = null;

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
      dlfPowPerm =
          alchemicalParameters.dlPowPerm * alchemicalParameters.dEdLSign;
      d2lfPowPerm = alchemicalParameters.d2lPowPerm
          * alchemicalParameters.dEdLSign;
      lfPowPol = alchemicalParameters.polarizationScale;
      dlfPowPol =
          alchemicalParameters.dlPowPol * alchemicalParameters.dEdLSign;
      d2lfPowPol =
          alchemicalParameters.d2lPowPol * alchemicalParameters.dEdLSign;
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
