/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.optimization.PointVectorValuePair;
import org.apache.commons.math3.optimization.SimpleVectorValueChecker;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.util.Range;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.TensorRecursion;
import ffx.numerics.VectorMath;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Torsion;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.ForceField.ForceFieldType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;

import static ffx.numerics.Erf.erfc;
import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.norm;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.parameters.MultipoleType.ELECTRIC;
import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t003;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t012;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t021;
import static ffx.potential.parameters.MultipoleType.t030;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t102;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t111;
import static ffx.potential.parameters.MultipoleType.t120;
import static ffx.potential.parameters.MultipoleType.t200;
import static ffx.potential.parameters.MultipoleType.t201;
import static ffx.potential.parameters.MultipoleType.t210;
import static ffx.potential.parameters.MultipoleType.t300;

/**
 * This Particle Mesh Ewald class implements PME for the AMOEBA polarizable
 * mutlipole force field in parallel using a {@link NeighborList} for any
 * {@link Crystal} space group. The real space contribution is contained within
 * this class and the reciprocal space contribution is delegated to the
 * {@link ReciprocalSpace} class.
 *
 * @author Michael J. Schnieders<br> derived from:<br> TINKER code by Jay
 * Ponder, Pengyu Ren and Tom Darden.<br>
 *
 * @see <a href="http://dx.doi.org/10.1021/ct300035u" target="_blank"> M. J.
 * Schnieders, J. Baltrusaitis, Y. Shi, G. Chattree, L. Zheng, W. Yang and P.
 * Ren, The Structure, Thermodynamics, and Solubility of Organic Crystals from
 * Simulation with a Polarizable Force Field, Journal of Chemical Theory and
 * Computation 8 (5), 1721-36 (2012)</a>
 *
 * @see <br><a href="http://dx.doi.org/10.1021/ct100506d" target="_blank"> M. J.
 * Schnieders, T. D. Fenn and V. S. Pande, Polarizable atomic multipole X-ray
 * refinement: Particle-mesh Ewald electrostatics for macromolecular crystals.
 * Journal of Chemical Theory and Computation 7 (4), 1141-56 (2011)</a>
 *
 * @see <br><a href="http://dx.doi.org/10.1063/1.1630791" target="_blank"> C.
 * Sagui, L. G. Pedersen, and T. A. Darden, Journal of Chemical Physics 120 (1),
 * 73 (2004)</a>
 *
 * @see <br><a href="http://link.aip.org/link/?JCPSA6/98/10089/1"
 * target="_blank"> T. Darden, D. York, and L. Pedersen, Journal of Chemical
 * Physics 98 (12), 10089 (1993)</a>
 *
 * @see <br><a href="http://www.ccp5.org" target="_blank"> W. Smith, "Point
 * Multipoles in the Ewald Summation (Revisited)", CCP5 Newsletter, 46, 18-30
 * (1998)</a>
 */
public class ParticleMeshEwald implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(ParticleMeshEwald.class.getName());

    public enum ELEC_FORM {

        PAM, FIXED_CHARGE
    }

    /**
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public enum Polarization {

        MUTUAL, DIRECT, NONE
    }

    private Resolution resolution = Resolution.AMOEBA;

    /**
     * Polarization mode.
     */
    protected final Polarization polarization;
    /**
     * Flag to indicate use of generalized Kirkwood.
     */
    private boolean generalizedKirkwoodTerm;
    /**
     * If lambdaTerm is true, some ligand atom interactions with the environment
     * are being turned on/off.
     */
    private final boolean lambdaTerm;
    /**
     * If true, compute coordinate gradient.
     */
    private boolean gradient = false;
    /**
     * If set to false, multipole charges are set to zero (default is true).
     */
    private final boolean useCharges;
    /**
     * If set to false, multipole dipoles are set to zero (default is true).
     */
    private final boolean useDipoles;
    /**
     * If set to false, multipole quadrupoles are set to zero (default is true).
     */
    private final boolean useQuadrupoles;
    /**
     * If set to false, multipoles are fixed in their local frame and torques
     * are zero, which is useful for narrowing down discrepancies between
     * analytic and finite-difference derivatives(default is true).
     */
    private final boolean rotateMultipoles;
    /**
     * Number of PME multipole interactions.
     */
    private int interactions;
    /**
     * Number of generalized Kirkwood interactions.
     */
    private int gkInteractions;
    /**
     * Permanent multipole energy (kcal/mol).
     */
    private double permanentMultipoleEnergy;
    /**
     * Polarization energy (kcal/mol).
     */
    private double polarizationEnergy;
    /**
     * Generalized Kirkwood energy.
     */
    private double generalizedKirkwoodEnergy;
    /**
     * Reference to the force field being used.
     */
    private final ForceField forceField;
    /**
     * Unit cell and spacegroup information.
     */
    private Crystal crystal;
    /**
     * Number of symmetry operators.
     */
    private int nSymm;
    /**
     * An ordered array of atoms in the system.
     */
    private Atom atoms[];
    /**
     * The number of atoms in the system.
     */
    private int nAtoms;
    /**
     * Dimensions of [nsymm][3][nAtoms].
     */
    protected double coordinates[][][];
    /**
     * Neighbor lists, including atoms beyond the real space cutoff.
     * [nsymm][nAtoms][nAllNeighbors]
     */
    protected int neighborLists[][][];
    /**
     * Neighbor lists, without atoms beyond the real space cutoff.
     * [nSymm][nAtoms][nIncludedNeighbors]
     */
    private int[][][] realSpaceLists;
    /**
     * Number of neighboring atoms within the real space cutoff. [nSymm][nAtoms]
     */
    private int[][] realSpaceCounts;
    /**
     * Optimal pairwise ranges.
     */
    private Range realSpaceRanges[];
    /**
     * Pairwise schedule for load balancing.
     */
    private IntegerSchedule realSpaceSchedule;
    /**
     * Neighbor lists, without atoms beyond the preconditioner cutoff.
     * [nSymm][nAtoms][nIncludedNeighbors]
     */
    private int[][][] preconditionerLists;
    /**
     * Number of neighboring atoms within the preconditioner cutoff.
     * [nSymm][nAtoms]
     */
    private int[][] preconditionerCounts;
    private double preconditionerCutoff = 4.5;
    private double preconditionerEwald = 0.0;
    private final int preconditionerListSize = 50;

    /**
     * *************************************************************************
     * Lambda state variables.
     */
    private enum LambdaMode {

        OFF, CONDENSED, CONDENSED_NO_LIGAND, VAPOR
    };
    private LambdaMode lambdaMode = LambdaMode.OFF;
    /**
     * Current state.
     */
    private double lambda = 1.0;
    /**
     * No electrostatics on softcore atoms.
     */
    private boolean noSoftcoreElectrostatics = false;
    /**
     * The polarization Lambda value goes from 0.0 .. 1.0 as the global lambda
     * value varies between polarizationLambdaStart .. 1.0.
     */
    private double polLambda = 1.0;
    /**
     * Constant α in: r' = sqrt(r^2 + α*(1 - L)^2)
     */
    private double permLambdaAlpha = 1.0;
    /**
     * Power on L in front of the pairwise multipole potential.
     */
    private double permLambdaExponent = 1.0;
    /**
     * Start turning on polarization later in the Lambda path to prevent SCF
     * convergence problems when atoms nearly overlap.
     */
    private double polLambdaStart = 0.75;
    private double polLambdaEnd = 1.0;
    /**
     * Power on L in front of the polarization energy.
     */
    private double polLambdaExponent = 3.0;
    /**
     * Intramolecular electrostatics for the ligand in vapor is included by
     * default.
     */
    private boolean doLigandVaporElec = true;
    /**
     * Condensed phase SCF without the ligand present is included by default.
     * For DualTopologyEnergy calculations it can be turned off.
     */
    private boolean doNoLigandCondensedSCF = true;
    /**
     * lAlpha = α*(1 - L)^2
     */
    private double lAlpha = 0.0;
    private double dlAlpha = 0.0;
    private double d2lAlpha = 0.0;
    private double dEdLSign = 1.0;
    /**
     * lPowPerm = L^permanentLambdaExponent
     */
    private double lPowPerm = 1.0;
    private double dlPowPerm = 0.0;
    private double d2lPowPerm = 0.0;
    private boolean doPermanentRealSpace;
    private double permanentScale = 1.0;
    /**
     * lPowPol = L^polarizationLambdaExponent
     */
    private double lPowPol = 1.0;
    private double dlPowPol = 0.0;
    private double d2lPowPol = 0.0;
    private boolean doPolarization;
    /**
     * Specify intermolecularSoftcore.
     */
    private boolean intermolecularSoftcore = false;
    /**
     * Molecule number for each atom.
     */
    private int molecule[] = null;

    /**
     * When computing the polarization energy at L there are 3 pieces.
     *
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand).
     *
     * 2.) Uenv = The polarization energy of the system without the ligand.
     *
     * 3.) Uligand = The polarization energy of the ligand by itself.
     *
     * Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     *
     * Set polarizationScale to L for part 1. Set polarizationScale to (1-L) for
     * parts 2 & 3.
     */
    private double polarizationScale = 1.0;
    /**
     * Flag for ligand atoms.
     */
    private boolean isSoft[];
    /**
     * Flag indicating if softcore variables have been initalized.
     */
    private boolean initSoftCore = false;
    /**
     * When computing the polarization energy at Lambda there are 3 pieces.
     *
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand).
     *
     * 2.) Uenv = The polarization energy of the system without the ligand.
     *
     * 3.) Uligand = The polarization energy of the ligand by itself.
     *
     * Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     *
     * Set the "use" array to true for all atoms for part 1. Set the "use" array
     * to true for all atoms except the ligand for part 2. Set the "use" array
     * to true only for the ligand atoms for part 3.
     *
     * The "use" array can also be employed to turn off atoms for computing the
     * electrostatic energy of sub-structures.
     */
    private boolean use[];
    private Crystal vaporCrystal;
    private int vaporLists[][][];
    private IntegerSchedule vaporPermanentSchedule;
    private IntegerSchedule vaporEwaldSchedule;
    private Range vacuumRanges[];
    /**
     * *************************************************************************
     * Permanent multipole variables.
     */
    /**
     * Permanent multipoles in their local frame.
     */
    private double localMultipole[][];
    private MultipoleType.MultipoleFrameDefinition frame[];
    private int axisAtom[][];
    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    protected double globalMultipole[][][];
    private double cartMultipolePhi[][];
    /**
     * The interaction energy between 1-2 multipoles is scaled by m12scale.
     */
    private final double m12scale;
    /**
     * The interaction energy between 1-3 multipoles is scaled by m13scale.
     */
    private final double m13scale;
    /**
     * The interaction energy between 1-4 multipoles is scaled by m14scale.
     */
    private final double m14scale;
    /**
     * The interaction energy between 1-5 multipoles is scaled by m15scale.
     */
    private final double m15scale;
    /**
     * *************************************************************************
     * Induced dipole variables.
     */
    private final double polsor;
    private final double poleps;
    /**
     * Direct polarization field due to permanent multipoles at polarizable
     * sites within their group are scaled. The scaling is 0.0 in AMOEBA.
     */
    private final double d11scale;
    /**
     * The interaction energy between a permanent multipole and polarizable site
     * that are 1-2 is scaled by p12scale.
     */
    private final double p12scale;
    /**
     * The interaction energy between a permanent multipole and polarizable site
     * that are 1-3 is scaled by p13scale.
     */
    private final double p13scale;
    private double ipdamp[];
    private double thole[];
    private double polarizability[];
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    protected double inducedDipole[][][];
    protected double inducedDipoleCR[][][];

    public enum SCFPredictor {

        NONE, LS, POLY, ASPC
    }
    /**
     * Specify an SCF predictor algorithm.
     */
    private SCFPredictor scfPredictor = SCFPredictor.LS;
    /**
     * Induced dipole predictor order.
     */
    private int predictorOrder;
    /**
     * Induced dipole predictor index.
     */
    private int predictorStartIndex;
    /**
     * Induced dipole predictor count.
     */
    private int predictorCount;
    /**
     * Dimensions of [mode][predictorOrder][nAtoms][3]
     */
    private double predictorInducedDipole[][][][];
    /**
     * Dimensions of [mode][predictorOrder][nAtoms][3]
     */
    private double predictorInducedDipoleCR[][][][];
    private LeastSquaresPredictor leastSquaresPredictor;
    private LevenbergMarquardtOptimizer leastSquaresOptimizer;

    public enum SCFAlgorithm {

        SOR, CG
    }
    private SCFAlgorithm scfAlgorithm = SCFAlgorithm.CG;
    /**
     * Direct induced dipoles.
     */
    private double directDipole[][];
    private double directDipoleCR[][];
    private double cartesianDipolePhi[][];
    private double cartesianDipolePhiCR[][];
    private int ip11[][];
    private int ip12[][];
    private int ip13[][];
    /**
     * *************************************************************************
     * Mutable Particle Mesh Ewald constants.
     */
    private double aewald;
    private double alsq2;
    private double an0;
    private double an1;
    private double an2;
    private double an3;
    private double an4;
    private double an5;
    private double piEwald;
    private double aewald3;
    private double off;
    private double off2;

    /**
     * PCG Variables.
     */
    private double rsd[][];
    private double rsdCR[][];
    private double rsdPre[][];
    private double rsdPreCR[][];
    private double conj[][];
    private double conjCR[][];
    private double vec[][];
    private double vecCR[][];

    /**
     * *************************************************************************
     * Parallel variables.
     */
    /**
     * By default, maxThreads is set to the number of available SMP cores.
     */
    private final int maxThreads;
    /**
     * Either 1 or 2; see description below.
     */
    private final int sectionThreads;
    /**
     * If real and reciprocal space are done sequentially or OpenCL is used,
     * then realSpaceThreads == maxThreads. Otherwise the number of
     * realSpaceThreads is set to ffx.realSpaceThreads.
     */
    private final int realSpaceThreads;
    /**
     * If real and reciprocal space are done sequentially then reciprocalThreads
     * == maxThreads If CUDA is used, reciprocalThreads == 1 Otherwise,
     * reciprocalThreads = maxThreads - realSpaceThreads
     */
    private final int reciprocalThreads;
    /**
     * Gradient array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double grad[][][];
    /**
     * Torque array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double torque[][][];
    /**
     * Field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double field[][][];
    /**
     * Chain rule field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double fieldCR[][][];
    /**
     * Partial derivative of the gradient with respect to Lambda.
     * [threadID][X/Y/Z][atomID]
     */
    private double lambdaGrad[][][];
    /**
     * Partial derivative of the torque with respect to Lambda.
     * [threadID][X/Y/Z][atomID]
     */
    private double lambdaTorque[][][];
    /**
     * Partial derivative with respect to Lambda.
     */
    private final SharedDouble shareddEdLambda;
    /**
     * Second partial derivative with respect to Lambda.
     */
    private final SharedDouble sharedd2EdLambda2;
    /**
     * The default ParallelTeam encapsulates the maximum number of threads used
     * to parallelize the electrostatics calculation.
     */
    private final ParallelTeam parallelTeam;
    /**
     * The sectionTeam encapsulates 1 or 2 threads.
     *
     * If it contains 1 thread, the real and reciprocal space calculations are
     * done sequentially.
     *
     * If it contains 2 threads, the real and reciprocal space calculations will
     * be done concurrently.
     */
    private final ParallelTeam sectionTeam;
    /**
     * If the real and reciprocal space parts of PME are done sequentially, then
     * the realSpaceTeam is equal parallalTeam.
     *
     * If the real and reciprocal space parts of PME are done concurrently, then
     * the realSpaceTeam will have fewer threads than the default parallelTeam.
     */
    private final ParallelTeam realSpaceTeam;
    /**
     * If the real and reciprocal space parts of PME are done sequentially, then
     * the reciprocalSpaceTeam is equal parallalTeam.
     *
     * If the real and reciprocal space parts of PME are done concurrently, then
     * the reciprocalSpaceTeam will have fewer threads than the default
     * parallelTeam.
     */
    private final ParallelTeam fftTeam;
    private final boolean cudaFFT;
    private IntegerSchedule permanentSchedule;
    private NeighborList neighborList;
    private final InitializationRegion initializationRegion;
    private final PermanentFieldRegion permanentFieldRegion;
    private final InducedDipoleFieldRegion inducedDipoleFieldRegion;
    private final ExpandInducedDipolesRegion expandInducedDipolesRegion;
    private final DirectRegion directRegion;
    private final SORRegion sorRegion;
    private final InducedDipolePreconditionerRegion inducedDipolePreconditionerRegion;
    private final PCGRegion pcgRegion;
    private final PCGInitRegion1 pcgInitRegion1;
    private final PCGInitRegion2 pcgInitRegion2;
    private final PCGIterRegion1 pcgIterRegion1;
    private final PCGIterRegion2 pcgIterRegion2;

    private final RealSpaceEnergyRegion realSpaceEnergyRegion;
    private final ReciprocalSpace reciprocalSpace;
    private final ReciprocalEnergyRegion reciprocalEnergyRegion;
    private final ReduceRegion reduceRegion;
    private final GeneralizedKirkwood generalizedKirkwood;
    /**
     * Timing variables.
     */
    private final long realSpacePermTime[];
    private final long realSpaceEnergyTime[];
    private final long realSpaceSCFTime[];
    private long realSpacePermTotal, realSpaceEnergyTotal, realSpaceSCFTotal;
    private long bornRadiiTotal, gkEnergyTotal;
    private ELEC_FORM elecForm = ELEC_FORM.PAM;
    private static final double toSeconds = 1.0e-9;
    /**
     * The sqrt of PI.
     */
    private static final double sqrtPi = sqrt(Math.PI);

    /**
     * ParticleMeshEwald constructor.
     *
     * @param atoms the Atom array to do electrostatics on.
     * @param molecule the molecule number for each atom.
     * @param forceField the ForceField the defines the electrostatics
     * parameters.
     * @param crystal The boundary conditions.
     * @param elecForm The electrostatics functional form.
     * @param neighborList The NeighborList for both van der Waals and PME.
     * @param parallelTeam A ParallelTeam that delegates parallelization.
     */
    public ParticleMeshEwald(Atom atoms[], int molecule[], ForceField forceField,
            Crystal crystal, NeighborList neighborList, ELEC_FORM elecForm, ParallelTeam parallelTeam) {
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

        polsor = forceField.getDouble(ForceFieldDouble.POLAR_SOR, 0.70);
        poleps = forceField.getDouble(ForceFieldDouble.POLAR_EPS, 1e-5);
        if (elecForm == ELEC_FORM.PAM) {
            m12scale = forceField.getDouble(ForceFieldDouble.MPOLE_12_SCALE, 0.0);
            m13scale = forceField.getDouble(ForceFieldDouble.MPOLE_13_SCALE, 0.0);
            m14scale = forceField.getDouble(ForceFieldDouble.MPOLE_14_SCALE, 0.4);
            m15scale = forceField.getDouble(ForceFieldDouble.MPOLE_15_SCALE, 0.8);
        } else {
            m12scale = forceField.getDouble(ForceFieldDouble.MPOLE_12_SCALE, 0.0);
            m13scale = forceField.getDouble(ForceFieldDouble.MPOLE_13_SCALE, 0.0);
            m14scale = forceField.getDouble(ForceFieldDouble.MPOLE_14_SCALE, 0.5);
            m15scale = forceField.getDouble(ForceFieldDouble.MPOLE_15_SCALE, 1.0);
        }
        d11scale = forceField.getDouble(ForceFieldDouble.DIRECT_11_SCALE, 0.0);
        p12scale = forceField.getDouble(ForceFieldDouble.POLAR_12_SCALE, 0.0);
        p13scale = forceField.getDouble(ForceFieldDouble.POLAR_13_SCALE, 0.0);
        useCharges = forceField.getBoolean(ForceFieldBoolean.USE_CHARGES, true);
        useDipoles = forceField.getBoolean(ForceFieldBoolean.USE_DIPOLES, true);
        useQuadrupoles = forceField.getBoolean(ForceFieldBoolean.USE_QUADRUPOLES, true);
        rotateMultipoles = forceField.getBoolean(ForceFieldBoolean.ROTATE_MULTIPOLES, true);
        lambdaTerm = forceField.getBoolean(ForceFieldBoolean.LAMBDATERM, false);

        if (!crystal.aperiodic()) {
            off = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        } else {
            off = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 100.0);
        }
        double ewald_precision = forceField.getDouble(ForceFieldDouble.EWALD_PRECISION, 1.0e-8);
        aewald = forceField.getDouble(ForceFieldDouble.EWALD_ALPHA, ewaldCoefficient(off, ewald_precision));
        setEwaldParameters(off, aewald);

        String predictor = forceField.getString(ForceFieldString.SCF_PREDICTOR, "NONE");
        try {
            predictor = predictor.replaceAll("-", "_").toUpperCase();
            scfPredictor = SCFPredictor.valueOf(predictor);
        } catch (Exception e) {
            scfPredictor = SCFPredictor.NONE;
        }

        if (scfPredictor != SCFPredictor.NONE) {
            predictorCount = 0;
            int defaultOrder = 6;
            predictorOrder = forceField.getInteger(ForceFieldInteger.SCF_PREDICTOR_ORDER, defaultOrder);
            if (scfPredictor == SCFPredictor.LS) {
                leastSquaresPredictor = new LeastSquaresPredictor();
                double eps = 1.0e-4;
                leastSquaresOptimizer = new LevenbergMarquardtOptimizer(new SimpleVectorValueChecker(eps, eps));
            } else if (scfPredictor == SCFPredictor.ASPC) {
                predictorOrder = 6;
            }
            predictorStartIndex = 0;
        }

        String algorithm = forceField.getString(ForceFieldString.SCF_ALGORITHM, "CG");
        try {
            algorithm = algorithm.replaceAll("-", "_").toUpperCase();
            scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
        } catch (Exception e) {
            scfAlgorithm = SCFAlgorithm.CG;
        }

        /**
         * The size of the preconditioner neighbor list depends on the size of
         * the preconditioner cutoff.
         */
        if (scfAlgorithm == SCFAlgorithm.CG) {
            inducedDipolePreconditionerRegion = new InducedDipolePreconditionerRegion(maxThreads);
            pcgRegion = new PCGRegion(maxThreads);
            pcgInitRegion1 = new PCGInitRegion1(maxThreads);
            pcgInitRegion2 = new PCGInitRegion2(maxThreads);
            pcgIterRegion1 = new PCGIterRegion1(maxThreads);
            pcgIterRegion2 = new PCGIterRegion2(maxThreads);
            boolean preconditioner = forceField.getBoolean(ForceFieldBoolean.USE_SCF_PRECONDITIONER, true);
            if (preconditioner) {
                preconditionerCutoff = forceField.getDouble(
                        ForceFieldDouble.CG_PRECONDITIONER_CUTOFF, 4.5);
                preconditionerEwald = forceField.getDouble(
                        ForceFieldDouble.CG_PRECONDITIONER_EWALD, 0.0);
            } else {
                preconditionerCutoff = 0.0;
            }
        } else {
            preconditionerCutoff = 0.0;
            inducedDipolePreconditionerRegion = null;
            pcgRegion = null;
            pcgInitRegion1 = null;
            pcgInitRegion2 = null;
            pcgIterRegion1 = null;
            pcgIterRegion2 = null;
        }

        if (lambdaTerm) {
            /**
             * Values of PERMANENT_LAMBDA_ALPHA below 2 can lead to unstable
             * trajectories.
             */
            permLambdaAlpha = forceField.getDouble(ForceFieldDouble.PERMANENT_LAMBDA_ALPHA, 2.0);
            if (permLambdaAlpha < 0.0 || permLambdaAlpha > 3.0) {
                permLambdaAlpha = 2.0;
            }
            /**
             * A PERMANENT_LAMBDA_EXPONENT of 1 gives linear charging of the
             * permanent electrostatics, which is most efficient. A quadratic
             * schedule (PERMANENT_LAMBDA_EXPONENT) also works, but the dU/dL
             * forces near lambda=1 are may be larger by a factor of 2.
             */
            permLambdaExponent = forceField.getDouble(ForceFieldDouble.PERMANENT_LAMBDA_EXPONENT, 1.0);
            if (permLambdaExponent < 1.0) {
                permLambdaExponent = 1.0;
            }
            /**
             * A POLARIZATION_LAMBDA_EXPONENT of 2 gives a non-zero d2U/dL2 at
             * the beginning of the polarization schedule. Choosing a power of 3
             * or greater ensures a smooth dU/dL and d2U/dL2 over the schedule.
             */
            polLambdaExponent = forceField.getDouble(ForceFieldDouble.POLARIZATION_LAMBDA_EXPONENT, 3.0);
            if (polLambdaExponent < 0.0) {
                polLambdaExponent = 0.0;
            }
            /**
             * The POLARIZATION_LAMBDA_START defines the point in the lambda
             * schedule when the condensed phase polarization of the ligand
             * begins to be turned on. If the condensed phase polarization is
             * considered near lambda=0, then SCF convergence is slow, even with
             * Thole damping. In addition, 2 (instead of 1) condensed phase SCF
             * calculations are necessary from the beginning of the window to
             * lambda=1.
             */
            polLambdaStart = forceField.getDouble(ForceFieldDouble.POLARIZATION_LAMBDA_START, 0.75);
            if (polLambdaStart < 0.0 || polLambdaStart > 0.9) {
                polLambdaStart = 0.75;
            }
            /**
             * The POLARIZATION_LAMBDA_END defines the point in the lambda
             * schedule when the condensed phase polarization of ligand has been
             * completely turned on. Values other than 1.0 have not been tested.
             */
            polLambdaEnd = forceField.getDouble(ForceFieldDouble.POLARIZATION_LAMBDA_END, 1.0);
            if (polLambdaEnd < polLambdaStart
                    || polLambdaEnd > 1.0
                    || polLambdaEnd - polLambdaStart < 0.3) {
                polLambdaEnd = 1.0;
            }

            /**
             * The LAMBDA_VAPOR_ELEC defines if intramolecular electrostatics of
             * the ligand in vapor will be considered.
             */
            doLigandVaporElec = forceField.getBoolean(ForceFieldBoolean.LIGAND_VAPOR_ELEC, true);
            doNoLigandCondensedSCF = forceField.getBoolean(ForceFieldBoolean.NO_LIGAND_CONDENSED_SCF, true);

            /**
             * Flag to indicate application of an intermolecular softcore
             * potential.
             */
            intermolecularSoftcore = forceField.getBoolean(
                    ForceField.ForceFieldBoolean.INTERMOLECULAR_SOFTCORE, false);
        }

        String polar = forceField.getString(ForceFieldString.POLARIZATION, "MUTUAL");
        if (elecForm == ELEC_FORM.FIXED_CHARGE) {
            polar = "NONE";
        }
        boolean polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);
        if (polarizationTerm == false || polar.equalsIgnoreCase("NONE")) {
            polarization = Polarization.NONE;
        } else if (polar.equalsIgnoreCase("DIRECT")) {
            polarization = Polarization.DIRECT;
        } else {
            polarization = Polarization.MUTUAL;
        }

        cudaFFT = forceField.getBoolean(ForceField.ForceFieldBoolean.CUDAFFT, false);

        if (lambdaTerm) {
            shareddEdLambda = new SharedDouble();
            sharedd2EdLambda2 = new SharedDouble();
        } else {
            shareddEdLambda = null;
            sharedd2EdLambda2 = null;
            lambdaGrad = null;
            lambdaTorque = null;
            vaporCrystal = null;
            vaporLists = null;
            vaporPermanentSchedule = null;
            vaporEwaldSchedule = null;
            vacuumRanges = null;
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder("\n Electrostatics\n");
            sb.append(format("  Polarization:                        %8s\n", polarization.toString()));
            if (polarization == Polarization.MUTUAL) {
                sb.append(format("   SCF Convergence Criteria:          %8.3e\n", poleps));
                sb.append(format("   SCF Predictor:                      %8s\n", scfPredictor));
                sb.append(format("   SCF Algorithm:                      %8s\n", scfAlgorithm));
                if (scfAlgorithm == SCFAlgorithm.SOR) {
                    sb.append(format("   SOR Parameter:                      %8.3f\n", polsor));
                } else {
                    sb.append(format("   CG Preconditioner Cut-Off:          %8.3f\n", preconditionerCutoff));
                    sb.append(format("   CG Preconditioner Ewald Coefficient:%8.3f\n", preconditionerEwald));
                }
            }
            if (aewald > 0.0) {
                sb.append("  Particle-mesh Ewald\n");
                sb.append(format("   Ewald Coefficient:                  %8.3f\n", aewald));
                sb.append(format("   Particle Cut-Off:                   %8.3f (A)", off));
            } else {
                sb.append(format("   Electrostatics Cut-Off:             %8.3f (A)\n", off));
            }
            logger.info(sb.toString());
        }

        if (cudaFFT) {
            sectionThreads = 2;
            realSpaceThreads = parallelTeam.getThreadCount();
            reciprocalThreads = 1;
            sectionTeam = new ParallelTeam(sectionThreads);
            realSpaceTeam = parallelTeam;
            fftTeam = new ParallelTeam(reciprocalThreads);
        } else {
            boolean concurrent;
            int realThreads = 1;
            try {
                realThreads = forceField.getInteger(ForceField.ForceFieldInteger.PME_REAL_THREADS);
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
                /**
                 * If pme-real-threads is not defined, then do real and
                 * reciprocal space parts sequentially.
                 */
                sectionThreads = 1;
                realSpaceThreads = maxThreads;
                reciprocalThreads = maxThreads;
                sectionTeam = new ParallelTeam(sectionThreads);
                realSpaceTeam = parallelTeam;
                fftTeam = parallelTeam;
            }
        }

        realSpaceRanges = new Range[maxThreads];
        initializationRegion = new InitializationRegion(maxThreads);
        expandInducedDipolesRegion = new ExpandInducedDipolesRegion(maxThreads);
        initAtomArrays();

        /**
         * Note that we always pass on the unit cell crystal to ReciprocalSpace
         * instance even if the real space calculations require a
         * ReplicatesCrystal.
         */
        if (aewald > 0.0) {
            reciprocalSpace = new ReciprocalSpace(this, crystal.getUnitCell(), forceField,
                    atoms, aewald, fftTeam, parallelTeam);
            reciprocalEnergyRegion = new ReciprocalEnergyRegion(maxThreads);
        } else {
            reciprocalSpace = null;
            reciprocalEnergyRegion = null;
        }
        permanentFieldRegion = new PermanentFieldRegion(realSpaceTeam);
        inducedDipoleFieldRegion = new InducedDipoleFieldRegion(realSpaceTeam);
        directRegion = new DirectRegion(maxThreads);
        sorRegion = new SORRegion(maxThreads);
        realSpaceEnergyRegion = new RealSpaceEnergyRegion(maxThreads);
        reduceRegion = new ReduceRegion(maxThreads);
        realSpacePermTime = new long[maxThreads];
        realSpaceEnergyTime = new long[maxThreads];
        realSpaceSCFTime = new long[maxThreads];

        /**
         * Generalized Kirkwood currently requires aperiodic Ewald. The GK
         * reaction field is added to the intra-molecular to give a
         * self-consistent reaction field.
         */
        generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
        if (generalizedKirkwoodTerm) {
            generalizedKirkwood = new GeneralizedKirkwood(forceField, atoms, this, crystal, parallelTeam);
        } else {
            generalizedKirkwood = null;
        }

        if (lambdaTerm) {
            StringBuilder sb = new StringBuilder(" Lambda Parameters\n");
            sb.append(format(" Permanent Multipole Softcore Alpha:      %5.3f\n", permLambdaAlpha));
            sb.append(format(" Permanent Multipole Lambda Exponent:     %5.3f\n", permLambdaExponent));
            if (polarization != Polarization.NONE) {
                sb.append(format(" Polarization Lambda Exponent:            %5.3f\n", polLambdaExponent));
                sb.append(format(" Polarization Lambda Range:      %5.3f .. %5.3f\n",
                        polLambdaStart, polLambdaEnd));
                sb.append(format(" Condensed SCF Without Ligand:            %B\n", doNoLigandCondensedSCF));
            }
            sb.append(format(" Vapor Electrostatics:                    %B\n", doLigandVaporElec));
            logger.info(sb.toString());
        }
    }

    private void initAtomArrays() {
        if (localMultipole == null || localMultipole.length < nAtoms) {
            localMultipole = new double[nAtoms][10];
            frame = new MultipoleType.MultipoleFrameDefinition[nAtoms];
            axisAtom = new int[nAtoms][];
            cartMultipolePhi = new double[nAtoms][tensorCount];
            directDipole = new double[nAtoms][3];
            directDipoleCR = new double[nAtoms][3];
            cartesianDipolePhi = new double[nAtoms][tensorCount];
            cartesianDipolePhiCR = new double[nAtoms][tensorCount];
            ip11 = new int[nAtoms][];
            ip12 = new int[nAtoms][];
            ip13 = new int[nAtoms][];
            thole = new double[nAtoms];
            ipdamp = new double[nAtoms];
            polarizability = new double[nAtoms];
            realSpaceSchedule = new PairwiseSchedule(maxThreads, nAtoms, realSpaceRanges);
            if (scfAlgorithm == SCFAlgorithm.CG) {
                rsd = new double[3][nAtoms];
                rsdCR = new double[3][nAtoms];
                rsdPre = new double[3][nAtoms];
                rsdPreCR = new double[3][nAtoms];
                conj = new double[3][nAtoms];
                conjCR = new double[3][nAtoms];
                vec = new double[3][nAtoms];
                vecCR = new double[3][nAtoms];
            }
            if (scfPredictor != SCFPredictor.NONE) {
                if (lambdaTerm) {
                    predictorInducedDipole = new double[3][predictorOrder][nAtoms][3];
                    predictorInducedDipoleCR = new double[3][predictorOrder][nAtoms][3];
                } else {
                    predictorInducedDipole = new double[1][predictorOrder][nAtoms][3];
                    predictorInducedDipoleCR = new double[1][predictorOrder][nAtoms][3];
                }
            }
            /**
             * Initialize per-thread memory for collecting the gradient, torque,
             * field and chain-rule field.
             */
            grad = new double[maxThreads][3][nAtoms];
            torque = new double[maxThreads][3][nAtoms];
            field = new double[maxThreads][3][nAtoms];
            fieldCR = new double[maxThreads][3][nAtoms];
            if (lambdaTerm) {
                lambdaGrad = new double[maxThreads][3][nAtoms];
                lambdaTorque = new double[maxThreads][3][nAtoms];
            }
            isSoft = new boolean[nAtoms];
            use = new boolean[nAtoms];

            coordinates = new double[nSymm][3][nAtoms];
            globalMultipole = new double[nSymm][nAtoms][10];
            inducedDipole = new double[nSymm][nAtoms][3];
            inducedDipoleCR = new double[nSymm][nAtoms][3];
            /**
             * The size of reduced neighbor list depends on the size of the real
             * space cutoff.
             */
            realSpaceLists = new int[nSymm][nAtoms][];
            realSpaceCounts = new int[nSymm][nAtoms];
            preconditionerLists = new int[nSymm][nAtoms][preconditionerListSize];
            preconditionerCounts = new int[nSymm][nAtoms];
        }

        /**
         * Initialize the soft core lambda mask to false for all atoms.
         */
        fill(isSoft, false);
        /**
         * Initialize the use mask to true for all atoms.
         */
        fill(use, true);
        /**
         * Assign multipole parameters.
         */
        assignMultipoles();
        /**
         * Assign polarization groups.
         */
        assignPolarizationGroups();
        /**
         * Fill the thole, inverse polarization damping and polarizability
         * arrays.
         */
        for (Atom ai : atoms) {
            PolarizeType polarizeType = ai.getPolarizeType();
            int index = ai.xyzIndex - 1;
            thole[index] = polarizeType.thole;
            ipdamp[index] = polarizeType.pdamp;
            if (!(ipdamp[index] > 0.0)) {
                ipdamp[index] = Double.POSITIVE_INFINITY;
            } else {
                ipdamp[index] = 1.0 / ipdamp[index];
            }
            polarizability[index] = polarizeType.polarizability;
        }
    }

    /**
     * Pass in atoms that have been assigned electrostatics from a fixed
     * charge force field.
     *
     * @param atoms
     */
    public void setFixedCharges(Atom atoms[]) {
        for (Atom ai : atoms) {
            if (ai.getResolution() == Resolution.FIXEDCHARGE) {
                int index = ai.xyzIndex - 1;
                polarizability[index] = 0.0;
                localMultipole[index][t000] = ai.getMultipoleType().charge;
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

    /**
     * Initialize a boolean array of soft atoms and, if requested, ligand vapor
     * electrostatics.
     */
    private void initSoftCoreInit(boolean print) {
        if (initSoftCore) {
            return;
        }

        /**
         * Initialize a boolean array that marks soft atoms.
         */
        StringBuilder sb = new StringBuilder("\n Softcore Atoms:\n");
        int count = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            if (ai.applyLambda()) {
                isSoft[i] = true;
                if (print) {
                    sb.append(ai.toString()).append("\n");
                }
                count++;
            }
        }
        if (count > 0 && print) {
            logger.info(sb.toString());
        }

        /**
         * Initialize boundary conditions, an n^2 neighbor list and parallel
         * scheduling for ligand vapor electrostatics.
         */
        if (doLigandVaporElec) {
            double maxr = 10.0;
            for (int i = 0; i < nAtoms; i++) {
                Atom ai = atoms[i];
                if (ai.applyLambda()) {
                    /**
                     * Determine ligand size.
                     */
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
            vaporCrystal = new Crystal(3 * vacuumOff, 3 * vacuumOff, 3 * vacuumOff, 90.0, 90.0, 90.0, "P1");
            vaporCrystal.setAperiodic(true);
            NeighborList vacuumNeighborList = new NeighborList(null, vaporCrystal, atoms, vacuumOff, 2.0, parallelTeam);

            vacuumNeighborList.setIntermolecular(false, molecule);

            vaporLists = new int[1][nAtoms][];
            double coords[][] = new double[1][nAtoms * 3];
            for (int i = 0; i < nAtoms; i++) {
                coords[0][i * 3] = atoms[i].getX();
                coords[0][i * 3 + 1] = atoms[i].getY();
                coords[0][i * 3 + 2] = atoms[i].getZ();
            }
            vacuumNeighborList.buildList(coords, vaporLists, isSoft, true, true);

            /*
             sb = new StringBuilder();
             for (int i=0; i<nAtoms; i++) {
             int list[] = vaporLists[0][i];
             sb.append(String.format(" Atom %d:", i+1));
             for (int j=0; j<list.length; j++) {
             sb.append(String.format(" %d", list[j]+1));
             }
             sb.append("\n");
             }
             logger.info(sb.toString());
             */
            vaporPermanentSchedule = vacuumNeighborList.getPairwiseSchedule();
            vaporEwaldSchedule = vaporPermanentSchedule;
            vacuumRanges = new Range[maxThreads];
        } else {
            vaporCrystal = null;
            vaporLists = null;
            vaporPermanentSchedule = null;
            vaporEwaldSchedule = null;
            vacuumRanges = null;
        }

        /**
         * Set this flag to true to avoid re-initialization.
         */
        initSoftCore = true;
    }

    public void setResolution(Resolution resolution) {
        this.resolution = resolution;
    }

    public void setAtoms(Atom atoms[], int molecule[]) {
        if (lambdaTerm) {
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
        /**
         * Check if memory allocation is required.
         */
        int nSymmNew = crystal.spaceGroup.getNumberOfSymOps();
        if (nSymm < nSymmNew) {
            coordinates = new double[nSymmNew][3][nAtoms];
            realSpaceLists = new int[nSymmNew][nAtoms][];
            realSpaceCounts = new int[nSymmNew][nAtoms];
            preconditionerLists = new int[nSymmNew][nAtoms][preconditionerListSize];
            preconditionerCounts = new int[nSymmNew][nAtoms];
            globalMultipole = new double[nSymmNew][nAtoms][10];
            inducedDipole = new double[nSymmNew][nAtoms][3];
            inducedDipoleCR = new double[nSymmNew][nAtoms][3];
        }
        nSymm = nSymmNew;
        neighborLists = neighborList.getNeighborList();
        this.crystal = crystal;
        /**
         * Production NPT simulations will include reciprocal space
         * contributions, but just in case there is a check for a NP.
         */
        if (reciprocalSpace != null) {
            reciprocalSpace.setCrystal(crystal.getUnitCell());
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

        /**
         * Initialize energy variables.
         */
        permanentMultipoleEnergy = 0.0;
        polarizationEnergy = 0.0;
        generalizedKirkwoodEnergy = 0.0;
        /**
         * Initialize number of interactions.
         */
        interactions = 0;
        gkInteractions = 0;
        /**
         * Initialize timing variables.
         */
        for (int i = 0; i < maxThreads; i++) {
            realSpacePermTime[i] = 0;
            realSpaceEnergyTime[i] = 0;
            realSpaceSCFTime[i] = 0;
        }
        realSpacePermTotal = 0;
        realSpaceEnergyTotal = 0;
        realSpaceSCFTotal = 0;
        gkEnergyTotal = 0;

        if (reciprocalSpace != null) {
            reciprocalSpace.initTimings();
        }

        /**
         * Initialize Lambda variables.
         */
        if (lambdaTerm) {
            shareddEdLambda.set(0.0);
            sharedd2EdLambda2.set(0.0);
        }
        doPermanentRealSpace = true;
        permanentScale = 1.0;
        doPolarization = true;
        polarizationScale = 1.0;

        /**
         * Total permanent + polarization energy.
         */
        double energy;

        /**
         * Expand the coordinates and rotate multipoles into the global frame.
         */
        try {
            parallelTeam.execute(initializationRegion);
        } catch (Exception e) {
            String message = "Fatal exception expanding coordinates and rotating multipoles.\n";
            logger.log(Level.SEVERE, message, e);
        }

        if (!lambdaTerm) {
            lambdaMode = LambdaMode.OFF;
            energy = computeEnergy(print);
        } else {
            /**
             * Condensed phase with all atoms.
             */
            lambdaMode = LambdaMode.CONDENSED;
            energy = condensedEnergy();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Solvated energy: %20.8f", energy));
            }
            /**
             * Condensed phase SCF without ligand atoms.
             */
            lambdaMode = LambdaMode.CONDENSED_NO_LIGAND;
            double temp = energy;
            energy = condensedNoLigandSCF();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Step 2 energy:   %20.8f", energy - temp));
            }

            /**
             * Vapor ligand electrostatics.
             */
            if (doLigandVaporElec) {
                lambdaMode = LambdaMode.VAPOR;
                temp = energy;
                energy = vaporElec();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(String.format(" Vacuum energy:   %20.8f", energy - temp));
                }
            }
        }

        /**
         * Convert torques to gradients on multipole frame defining atoms. Add
         * to electrostatic gradient to the total XYZ gradient.
         */
        if (gradient || lambdaTerm) {
            try {
                parallelTeam.execute(reduceRegion);
            } catch (Exception e) {
                String message = "Exception calculating torques.";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Log some timings.
         */
        if (logger.isLoggable(Level.FINE)) {
            printRealSpaceTimings();
            if (aewald > 0.0) {
                reciprocalSpace.printTimings();
            }
        }

        return permanentMultipoleEnergy + polarizationEnergy;
    }

    private void printRealSpaceTimings() {

        double total = (realSpacePermTotal + realSpaceSCFTotal + realSpaceEnergyTotal) * toSeconds;

        logger.info(String.format("\n Real Space: %7.4f (sec)", total));
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
            int count = realSpaceEnergyRegion.realSpaceEnergyLoop[i].getCount();
            logger.info(String.format("    %3d   %7.4f %7.4f %7.4f %10d", i,
                    realSpacePermTime[i] * toSeconds, realSpaceSCFTime[i] * toSeconds,
                    realSpaceEnergyTime[i] * toSeconds, count));
            minPerm = min(realSpacePermTime[i], minPerm);
            maxPerm = max(realSpacePermTime[i], maxPerm);
            minSCF = min(realSpaceSCFTime[i], minSCF);
            maxSCF = max(realSpaceSCFTime[i], maxSCF);
            minEnergy = min(realSpaceEnergyTime[i], minEnergy);
            maxEnergy = max(realSpaceEnergyTime[i], maxEnergy);
            minCount = min(count, minCount);
            maxCount = max(count, maxCount);
        }
        logger.info(String.format(" Min      %7.4f %7.4f %7.4f %10d",
                minPerm * toSeconds, minSCF * toSeconds,
                minEnergy * toSeconds, minCount));
        logger.info(String.format(" Max      %7.4f %7.4f %7.4f %10d",
                maxPerm * toSeconds, maxSCF * toSeconds,
                maxEnergy * toSeconds, maxCount));
        logger.info(String.format(" Delta    %7.4f %7.4f %7.4f %10d",
                (maxPerm - minPerm) * toSeconds, (maxSCF - minSCF) * toSeconds,
                (maxEnergy - minEnergy) * toSeconds, (maxCount - minCount)));
        logger.info(String.format(" Actual   %7.4f %7.4f %7.4f %10d",
                realSpacePermTotal * toSeconds, realSpaceSCFTotal * toSeconds,
                realSpaceEnergyTotal * toSeconds, realSpaceEnergyRegion.getInteractions()));
    }

    /**
     * 1.) Total system under PBC. A.) Softcore real space for Ligand-Protein
     * and Ligand-Ligand. B.) Reciprocal space scaled by lambda. C.)
     * Polarization scaled by lambda.
     */
    private double condensedEnergy() {
        if (lambda < polLambdaStart) {
            /**
             * If the polarization has been completely decoupled, the
             * contribution of the complete system is zero.
             *
             * We can skip the SCF for part 1 for efficiency.
             */
            polarizationScale = 0.0;
            doPolarization = false;
        } else if (lambda <= polLambdaEnd) {
            polarizationScale = lPowPol;
            doPolarization = true;
        } else {
            polarizationScale = 1.0;
            doPolarization = true;
        }
        doPermanentRealSpace = true;
        permanentScale = lPowPerm;
        dEdLSign = 1.0;

        double energy = computeEnergy(false);

        return energy;
    }

    /**
     * 2.) Condensed phase system without the ligand. A.) No permanent real
     * space electrostatics needs to be calculated because this was handled
     * analytically in step 1.
     *
     * B.) Permanent reciprocal space scaled by (1 - lambda).
     *
     * C.) Polarization scaled by (1 - lambda).
     */
    private double condensedNoLigandSCF() {
        /**
         * Turn off the ligand.
         */
        boolean skip = true;
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].applyLambda()) {
                use[i] = false;
            } else {
                use[i] = true;
                skip = false;
            }
        }

        /**
         * Permanent real space is done for the condensed phase. Scale the
         * reciprocal space part.
         */
        doPermanentRealSpace = false;
        permanentScale = 1.0 - lPowPerm;
        dEdLSign = -1.0;

        /**
         * If we are past the end of the polarization lambda window, then only
         * the condensed phase is necessary.
         */
        if (lambda <= polLambdaEnd && doNoLigandCondensedSCF) {
            doPolarization = true;
            polarizationScale = 1.0 - lPowPol;
        } else {
            doPolarization = false;
            polarizationScale = 0.0;
        }

        /**
         * Turn off GK.
         */
        boolean gkBack = generalizedKirkwoodTerm;
        generalizedKirkwoodTerm = false;

        /*
         * If we are disappearing the entire system (ie. a small crystal) then
         * the energy of this step is 0 and we can skip it.
         */
        double energy;
        if (skip) {
            energy = permanentMultipoleEnergy + polarizationEnergy + generalizedKirkwoodEnergy;
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
     * 3.) Ligand in vapor
     *
     * A.) Real space with an Ewald coefficient of 0.0 (no reciprocal space).
     *
     * B.) Polarization scaled as in Step 2 by (1 - lambda).
     */
    private double vaporElec() {
        for (int i = 0; i < nAtoms; i++) {
            use[i] = atoms[i].applyLambda();
        }

        /**
         * Scale the permanent vacuum electrostatics. The softcore alpha is not
         * necessary (nothing in vacuum to collide with).
         */
        doPermanentRealSpace = true;
        permanentScale = 1.0 - lPowPerm;
        dEdLSign = -1.0;
        double lAlphaBack = lAlpha;
        double dlAlphaBack = dlAlpha;
        double d2lAlphaBack = d2lAlpha;
        lAlpha = 0.0;
        dlAlpha = 0.0;
        d2lAlpha = 0.0;

        /**
         * If we are past the end of the polarization lambda window, then only
         * the condensed phase is necessary.
         */
        if (lambda <= polLambdaEnd) {
            doPolarization = true;
            polarizationScale = 1.0 - lPowPol;
        } else {
            doPolarization = false;
            polarizationScale = 0.0;
        }

        /**
         * Save the current real space PME parameters.
         */
        double offBack = off;
        double aewaldBack = aewald;
        off = Double.MAX_VALUE;
        aewald = 0.0;
        setEwaldParameters(off, aewald);

        /**
         * Save the current parallelization schedule.
         */
        IntegerSchedule permanentScheduleBack = permanentSchedule;
        IntegerSchedule ewaldScheduleBack = realSpaceSchedule;
        Range rangesBack[] = realSpaceRanges;
        permanentSchedule = vaporPermanentSchedule;
        realSpaceSchedule = vaporEwaldSchedule;
        realSpaceRanges = vacuumRanges;

        /**
         * Use vacuum crystal / vacuum neighborLists.
         */
        Crystal crystalBack = crystal;
        int nSymmBack = nSymm;
        int listsBack[][][] = neighborLists;
        neighborLists = vaporLists;
        crystal = vaporCrystal;
        nSymm = 1;

        /**
         * Turn off GK if in use.
         */
        boolean gkBack = generalizedKirkwoodTerm;
        generalizedKirkwoodTerm = false;

        double energy = computeEnergy(false);

        /**
         * Revert to the saved parameters.
         */
        off = offBack;
        aewald = aewaldBack;
        setEwaldParameters(off, aewald);
        neighborLists = listsBack;
        crystal = crystalBack;
        nSymm = nSymmBack;
        permanentSchedule = permanentScheduleBack;
        realSpaceSchedule = ewaldScheduleBack;
        realSpaceRanges = rangesBack;
        lAlpha = lAlphaBack;
        dlAlpha = dlAlphaBack;
        d2lAlpha = d2lAlphaBack;
        generalizedKirkwoodTerm = gkBack;

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
        /**
         * Initialize the energy components to zero.
         */
        double eself = 0.0;
        double erecip = 0.0;
        double ereal = 0.0;
        double eselfi = 0.0;
        double erecipi = 0.0;
        double ereali = 0.0;

        /**
         * Find the permanent multipole potential, field, etc.
         */
        try {
            /**
             * Compute b-Splines and permanent density.
             */
            if (aewald > 0.0) {
                reciprocalSpace.computeBSplines();
                reciprocalSpace.splinePermanentMultipoles(globalMultipole, use);
            }

            /**
             * The real space contribution can be calculated at the same time
             * the reciprocal space convolution is being done.
             */
            sectionTeam.execute(permanentFieldRegion);

            /**
             * Collect the reciprocal space field.
             */
            if (aewald > 0.0) {
                reciprocalSpace.computePermanentPhi(cartMultipolePhi);
            }
        } catch (Exception e) {
            String message = "Fatal exception computing the permanent multipole field.\n";
            logger.log(Level.SEVERE, message, e);
        }

        /**
         * Compute Born radii if necessary.
         */
        if (generalizedKirkwoodTerm) {
            bornRadiiTotal -= System.nanoTime();
            generalizedKirkwood.setUse(use);
            generalizedKirkwood.computeBornRadii();
            bornRadiiTotal += System.nanoTime();
        }

        /**
         * Do the self-consistent field calculation.
         */
        if (polarization != Polarization.NONE && doPolarization) {
            selfConsistentField(logger.isLoggable(Level.FINE));
            if (aewald > 0.0) {
                if (gradient && polarization == Polarization.DIRECT) {
                    try {
                        reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
                        sectionTeam.execute(inducedDipoleFieldRegion);
                        reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
                    } catch (Exception ex) {
                        String message = "Fatal exception computing the induced reciprocal space field.\n";
                        logger.log(Level.SEVERE, message, ex);
                    }
                } else {
                    reciprocalSpace.cartToFracInducedDipoles(inducedDipole, inducedDipoleCR);
                }
            }
            if (scfPredictor != SCFPredictor.NONE) {
                saveMutualInducedDipoles();
            }
        }

        /**
         * Find the total real space energy. This includes the permanent
         * multipoles in their own real space potential and the interaction of
         * permanent multipoles with induced dipoles.
         *
         * Then compute the permanent and reciprocal space energy.
         */
        try {
            if (aewald > 0.0) {
                parallelTeam.execute(reciprocalEnergyRegion);
                interactions += nAtoms;
                eself = reciprocalEnergyRegion.getPermanentSelfEnergy();
                erecip = reciprocalEnergyRegion.getPermanentReciprocalEnergy();
                eselfi = reciprocalEnergyRegion.getInducedDipoleSelfEnergy();
                erecipi = reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
            }
            realSpaceEnergyTotal -= System.nanoTime();
            parallelTeam.execute(realSpaceEnergyRegion);
            realSpaceEnergyTotal += System.nanoTime();
            ereal = realSpaceEnergyRegion.getPermanentEnergy();
            ereali = realSpaceEnergyRegion.getPolarizationEnergy();
            interactions += realSpaceEnergyRegion.getInteractions();
        } catch (Exception e) {
            String message = "Exception computing the electrostatic energy.\n";
            logger.log(Level.SEVERE, message, e);
        }

        /**
         * Compute the generalized Kirkwood solvation free energy.
         */
        if (generalizedKirkwoodTerm) {
            gkEnergyTotal -= System.nanoTime();
            generalizedKirkwoodEnergy += generalizedKirkwood.solvationEnergy(gradient, print);
            gkInteractions += generalizedKirkwood.getInteractions();
            gkEnergyTotal += System.nanoTime();
        }

        /**
         * Collect energy terms.
         */
        permanentMultipoleEnergy += eself + erecip + ereal;
        polarizationEnergy += eselfi + erecipi + ereali;

        /**
         * Log some info.
         */
        if (logger.isLoggable(Level.FINE)) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n Multipole Self-Energy:   %16.8f\n", eself));
            sb.append(format(" Multipole Reciprocal:    %16.8f\n", erecip));
            sb.append(format(" Multipole Real Space:    %16.8f\n", ereal));
            sb.append(format(" Polarization Self-Energy:%16.8f\n", eselfi));
            sb.append(format(" Polarization Reciprocal: %16.8f\n", erecipi));
            sb.append(format(" Polarization Real Space: %16.8f\n", ereali));
            if (generalizedKirkwoodTerm) {
                sb.append(format(" Generalized Kirkwood:    %16.8f\n", generalizedKirkwoodEnergy));
            }
            logger.fine(sb.toString());
        }

        return permanentMultipoleEnergy + polarizationEnergy + generalizedKirkwoodEnergy;
    }

    /**
     * <p>
     * Getter for the field <code>interactions</code>.</p>
     *
     * @return a int.
     */
    public int getInteractions() {
        return interactions;
    }

    /**
     * <p>
     * getPermanentEnergy</p>
     *
     * @return a double.
     */
    public double getPermanentEnergy() {
        return permanentMultipoleEnergy;
    }

    /**
     * <p>
     * Getter for the field <code>polarizationEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    /**
     * <p>
     * getGKEnergy</p>
     *
     * @return a double.
     */
    public double getGKEnergy() {
        return generalizedKirkwoodEnergy;
    }

    /**
     * <p>
     * getGKInteractions</p>
     *
     * @return a int.
     */
    public int getGKInteractions() {
        return gkInteractions;
    }

    /**
     * <p>
     * getGradients</p>
     *
     * @param grad an array of double.
     */
    public void getGradients(double grad[][]) {
        double gx[] = this.grad[0][0];
        double gy[] = this.grad[0][1];
        double gz[] = this.grad[0][2];
        double x[] = grad[0];
        double y[] = grad[1];
        double z[] = grad[2];
        for (int i = 0; i < nAtoms; i++) {
            x[i] = gx[i];
            y[i] = gy[i];
            z[i] = gz[i];
        }
    }

    /**
     * <p>
     * Getter for the field <code>gradient</code>.</p>
     *
     * @return an array of double.
     */
    protected double[][][] getGradient() {
        return grad;
    }

    /**
     * <p>
     * Getter for the field <code>torque</code>.</p>
     *
     * @return an array of double.
     */
    protected double[][][] getTorque() {
        return torque;
    }

    protected double[][][] getLambdaGradient() {
        return lambdaGrad;
    }

    protected double[][][] getLambdaTorque() {
        return lambdaTorque;
    }

    /**
     * Apply the selected polarization model (NONE, Direct or Mutual).
     */
    private int selfConsistentField(boolean print) {
        if (polarization == Polarization.NONE) {
            return -1;
        }
        long startTime = System.nanoTime();

        /**
         * Compute the direct induced dipoles.
         */
        try {
            if (generalizedKirkwoodTerm) {
                gkEnergyTotal = -System.nanoTime();
                generalizedKirkwood.computePermanentGKField();
                gkEnergyTotal += System.nanoTime();
                logger.fine(String.format(" Computed GK permanent field %8.3f (sec)", gkEnergyTotal * 1.0e-9));
            }
            parallelTeam.execute(directRegion);
        } catch (Exception e) {
            String message = " Exception computing direct induced dipoles.";
            logger.log(Level.SEVERE, message, e);
        }

        /**
         * Return unless mutual polarization is selected.
         */
        if (polarization != Polarization.MUTUAL) {
            if (nSymm > 1) {
                try {
                    parallelTeam.execute(expandInducedDipolesRegion);
                } catch (Exception e) {
                    String message = " Exception expanding direct induced dipoles.";
                    logger.log(Level.SEVERE, message, e);
                }
            }
            return 0;
        }

        /**
         * Predict the current self-consistent induced dipoles using information
         * from previous steps.
         */
        if (scfPredictor != SCFPredictor.NONE) {
            switch (scfPredictor) {
                case ASPC:
                    aspcPredictor();
                    break;
                case LS:
                    leastSquaresPredictor();
                    break;
                case POLY:
                    polynomialPredictor();
                    break;
                case NONE:
                default:
                    break;
            }
        }

        /**
         * Expand the initial induced dipoles to P1 symmetry, if necessary.
         */
        if (nSymm > 1) {
            try {
                parallelTeam.execute(expandInducedDipolesRegion);
            } catch (Exception e) {
                String message = " Exception expanding initial induced dipoles.";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Converge the self-consistent field.
         */
        int iterations;
        switch (scfAlgorithm) {
            case SOR:
                iterations = scfBySOR(print, startTime);
                break;
            case CG:
            default:
                //iterations = scfByCG();
                iterations = scfByPCG(print, startTime);
                break;
        }

        /**
         * for (int i = 0; i < nAtoms; i++) { logger.info(format(" %d ID (%8.3f
         * %8.3f %8.3f) CR (%8.3f %8.3f %8.3f)", i, inducedDipole[0][i][0],
         * inducedDipole[0][i][1], inducedDipole[0][i][2],
         * inducedDipoleCR[0][i][0], inducedDipoleCR[0][i][1],
         * inducedDipoleCR[0][i][2])); }
         */
        return iterations;
    }

    /**
     * Converge the SCF using Successive Over-Relaxation (SOR).
     */
    private int scfBySOR(boolean print, long startTime) {
        long directTime = System.nanoTime() - startTime;
        /**
         * A request of 0 SCF cycles simplifies mutual polarization to direct
         * polarization.
         */
        StringBuilder sb = null;
        if (print) {
            sb = new StringBuilder(
                    "\n Self-Consistent Field\n"
                    + " Iter  RMS Change (Debye)  Time\n");
        }
        int completedSCFCycles = 0;
        int maxSCFCycles = 1000;
        double eps = 100.0;
        double previousEps;
        boolean done = false;
        while (!done) {
            long cycleTime = -System.nanoTime();
            try {
                if (aewald > 0.0) {
                    reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
                }
                sectionTeam.execute(inducedDipoleFieldRegion);
                if (aewald > 0.0) {
                    reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
                }

                if (generalizedKirkwoodTerm) {
                    /**
                     * GK field.
                     */
                    gkEnergyTotal = -System.nanoTime();
                    generalizedKirkwood.computeInducedGKField();
                    gkEnergyTotal += System.nanoTime();
                    logger.fine(String.format(" Computed GK induced field %8.3f (sec)", gkEnergyTotal * 1.0e-9));
                }
                parallelTeam.execute(sorRegion);
                if (nSymm > 1) {
                    parallelTeam.execute(expandInducedDipolesRegion);
                }
            } catch (Exception e) {
                String message = "Exception computing mutual induced dipoles.";
                logger.log(Level.SEVERE, message, e);
            }
            completedSCFCycles++;
            previousEps = eps;
            eps = sorRegion.getEps();
            eps = MultipoleType.DEBYE * sqrt(eps / (double) nAtoms);
            cycleTime += System.nanoTime();
            if (print) {
                sb.append(format(
                        " %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * toSeconds));
            }
            /**
             * If the RMS Debye change increases, fail the SCF process.
             */
            if (eps > previousEps) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Fatal SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
                throw new ArithmeticException(message);
            }
            /**
             * The SCF should converge well before the max iteration check.
             * Otherwise, fail the SCF process.
             */
            if (completedSCFCycles >= maxSCFCycles) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Maximum SCF iterations reached: (%d)\n", completedSCFCycles);
                throw new ArithmeticException(message);
            }
            /**
             * Check if the convergence criteria has been achieved.
             */
            if (eps < poleps) {
                done = true;
            }
        }
        if (print) {
            sb.append(format(" Direct:                  %7.4f\n",
                    toSeconds * directTime));
            startTime = System.nanoTime() - startTime;
            sb.append(format(" Total:                   %7.4f",
                    startTime * toSeconds));
            logger.info(sb.toString());
        }
        return completedSCFCycles;
    }

    public void destroy() throws Exception {
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
     * The Permanent Field Region should be executed by a ParallelTeam with
     * exactly 2 threads. The Real Space and Reciprocal Space Sections will be
     * run concurrently, each with the number of threads defined by their
     * respective ParallelTeam instances.
     */
    private class PermanentFieldRegion extends ParallelRegion {

        private PermanentRealSpaceFieldSection permanentRealSpaceFieldSection;
        private PermanentReciprocalSection permanentReciprocalSection;

        public PermanentFieldRegion(ParallelTeam pt) {
            permanentRealSpaceFieldSection = new PermanentRealSpaceFieldSection(pt);
            permanentReciprocalSection = new PermanentReciprocalSection();
        }

        @Override
        public void run() {
            try {
                execute(permanentRealSpaceFieldSection, permanentReciprocalSection);
            } catch (Exception e) {
                String message = "Fatal exception computing the permanent multipole field.\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Computes the Permanent Multipole Real Space Field.
         */
        private class PermanentRealSpaceFieldSection extends ParallelSection {

            private final PermanentRealSpaceFieldRegion permanentRealSpaceFieldRegion;
            private final ParallelTeam parallelTeam;

            public PermanentRealSpaceFieldSection(ParallelTeam pt) {
                this.parallelTeam = pt;
                int nt = pt.getThreadCount();
                permanentRealSpaceFieldRegion = new PermanentRealSpaceFieldRegion(nt);
            }

            @Override
            public void run() {
                try {
                    realSpacePermTotal -= System.nanoTime();
                    parallelTeam.execute(permanentRealSpaceFieldRegion);
                    realSpacePermTotal += System.nanoTime();
                } catch (Exception e) {
                    String message = "Fatal exception computing the real space field.\n";
                    logger.log(Level.SEVERE, message, e);
                }
            }
        }

        /**
         * Compute the permanent multipole reciprocal space contribution to the
         * electric potential, field, etc. using the number of threads specified
         * by the ParallelTeam used to construct the ReciprocalSpace instance.
         */
        private class PermanentReciprocalSection extends ParallelSection {

            @Override
            public void run() {
                if (aewald > 0.0) {
                    reciprocalSpace.permanentMultipoleConvolution();
                }
            }
        }

        private class PermanentRealSpaceFieldRegion extends ParallelRegion {

            private final InitializationLoop initializationLoop[];
            private final PermanentRealSpaceFieldLoop permanentRealSpaceFieldLoop[];
            private final SharedInteger sharedCount;
            private final int threadCount;

            public PermanentRealSpaceFieldRegion(int nt) {
                threadCount = nt;
                permanentRealSpaceFieldLoop = new PermanentRealSpaceFieldLoop[threadCount];
                initializationLoop = new InitializationLoop[threadCount];
                sharedCount = new SharedInteger();
            }

            @Override
            public void start() {
                sharedCount.set(0);
            }

            @Override
            public void run() {
                int threadIndex = getThreadIndex();
                if (initializationLoop[threadIndex] == null) {
                    initializationLoop[threadIndex] = new InitializationLoop();
                    permanentRealSpaceFieldLoop[threadIndex] = new PermanentRealSpaceFieldLoop();
                }
                try {
                    execute(0, nAtoms - 1, initializationLoop[threadIndex]);
                    execute(0, nAtoms - 1, permanentRealSpaceFieldLoop[threadIndex]);
                } catch (Exception e) {
                    String message = "Fatal exception computing the real space field in thread " + getThreadIndex() + "\n";
                    logger.log(Level.SEVERE, message, e);
                }
            }

            @Override
            public void finish() {
                /**
                 * Load balancing.
                 */
                int id = 0;
                int goal = sharedCount.get() / threadCount;
                int num = 0;
                int start = 0;
                for (int i = 0; i < nAtoms; i++) {
                    for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                        num += realSpaceCounts[iSymm][i];
                    }
                    if (num >= goal) {
                        /**
                         * Last thread gets the remaining atoms in its range.
                         */
                        if (id == threadCount - 1) {
                            realSpaceRanges[id] = new Range(start, nAtoms - 1);
                            break;
                        }

                        realSpaceRanges[id] = new Range(start, i);

                        // Reset the count.
                        num = 0;

                        // Next thread.
                        id++;

                        // Next range starts at i+1.
                        start = i + 1;

                        /**
                         * Out of atoms. Threads remaining get a null range.
                         */
                        if (start == nAtoms) {
                            for (int j = id; j < threadCount; j++) {
                                realSpaceRanges[j] = null;
                            }
                            break;
                        }
                    } else if (i == nAtoms - 1) {
                        /**
                         * Last atom without reaching goal for current thread.
                         */
                        realSpaceRanges[id] = new Range(start, nAtoms - 1);
                        for (int j = id + 1; j < threadCount; j++) {
                            realSpaceRanges[j] = null;
                        }
                    }
                }
            }

            private class InitializationLoop extends IntegerForLoop {

                @Override
                public IntegerSchedule schedule() {
                    return IntegerSchedule.fixed();
                }

                /**
                 * Initialize the field arrays.
                 */
                @Override
                public void start() {
                    int threadIndex = getThreadIndex();
                    realSpacePermTime[threadIndex] -= System.nanoTime();
                    double fX[] = field[threadIndex][0];
                    double fY[] = field[threadIndex][1];
                    double fZ[] = field[threadIndex][2];
                    double fXCR[] = fieldCR[threadIndex][0];
                    double fYCR[] = fieldCR[threadIndex][1];
                    double fZCR[] = fieldCR[threadIndex][2];
                    fill(fX, 0.0);
                    fill(fY, 0.0);
                    fill(fZ, 0.0);
                    fill(fXCR, 0.0);
                    fill(fYCR, 0.0);
                    fill(fZCR, 0.0);
                }

                @Override
                public void finish() {
                    int threadIndex = getThreadIndex();
                    realSpacePermTime[threadIndex] += System.nanoTime();
                }

                @Override
                public void run(int lb, int ub) {
                    /**
                     * Initialize the induced dipole arrays.
                     */
                    for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                        double ind0[][] = inducedDipole[0];
                        double indCR0[][] = inducedDipoleCR[0];
                        for (int i = lb; i <= ub; i++) {
                            double ind[] = ind0[i];
                            double indCR[] = indCR0[i];
                            ind[0] = 0.0;
                            ind[1] = 0.0;
                            ind[2] = 0.0;
                            indCR[0] = 0.0;
                            indCR[1] = 0.0;
                            indCR[2] = 0.0;
                        }
                    }
                }
            }

            private class PermanentRealSpaceFieldLoop extends IntegerForLoop {

                private final double dx_local[];
                private final double transOp[][];
                private double fX[], fY[], fZ[];
                private double fXCR[], fYCR[], fZCR[];
                private double mask_local[];
                private double maskp_local[];
                private int count;
                // Extra padding to avert cache interference.
                private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
                private long pad8, pad9, pada, padb, padc, padd, pade, padf;

                public PermanentRealSpaceFieldLoop() {
                    super();
                    dx_local = new double[3];
                    transOp = new double[3][3];
                }

                @Override
                public void start() {
                    int threadIndex = getThreadIndex();
                    realSpacePermTime[threadIndex] -= System.nanoTime();
                    count = 0;
                    fX = field[threadIndex][0];
                    fY = field[threadIndex][1];
                    fZ = field[threadIndex][2];
                    fXCR = fieldCR[threadIndex][0];
                    fYCR = fieldCR[threadIndex][1];
                    fZCR = fieldCR[threadIndex][2];
                    if (mask_local == null || mask_local.length < nAtoms) {
                        mask_local = new double[nAtoms];
                        maskp_local = new double[nAtoms];
                        fill(mask_local, 1.0);
                        fill(maskp_local, 1.0);
                    }
                }

                @Override
                public void finish() {
                    int threadIndex = getThreadIndex();
                    sharedCount.addAndGet(count);
                    realSpacePermTime[threadIndex] += System.nanoTime();
                }

                @Override
                public IntegerSchedule schedule() {
                    return permanentSchedule;
                }

                @Override
                public void run(int lb, int ub) {
                    int lists[][] = neighborLists[0];
                    int ewalds[][] = realSpaceLists[0];
                    int counts[] = realSpaceCounts[0];
                    int preLists[][] = preconditionerLists[0];
                    int preCounts[] = preconditionerCounts[0];
                    final double x[] = coordinates[0][0];
                    final double y[] = coordinates[0][1];
                    final double z[] = coordinates[0][2];
                    final double mpole[][] = globalMultipole[0];
                    /**
                     * Loop over atom chunk.
                     */
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        final int moleculei = molecule[i];
                        final double pdi = ipdamp[i];
                        final double pti = thole[i];
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final double globalMultipolei[] = mpole[i];
                        final double ci = globalMultipolei[0];
                        final double dix = globalMultipolei[t100];
                        final double diy = globalMultipolei[t010];
                        final double diz = globalMultipolei[t001];
                        final double qixx = globalMultipolei[t200] * oneThird;
                        final double qiyy = globalMultipolei[t020] * oneThird;
                        final double qizz = globalMultipolei[t002] * oneThird;
                        final double qixy = globalMultipolei[t110] * oneThird;
                        final double qixz = globalMultipolei[t101] * oneThird;
                        final double qiyz = globalMultipolei[t011] * oneThird;
                        /**
                         * Apply energy masking rules.
                         */
                        Atom ai = atoms[i];
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                for (int k : ip11[i]) {
                                    if (k == index) {
                                        maskp_local[index] = 0.5;
                                    }
                                }
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                maskp_local[index] = p13scale;
                            }
                        }
                        for (Bond bond : ai.getBonds()) {
                            int index = bond.get1_2(ai).xyzIndex - 1;
                            maskp_local[index] = p12scale;
                        }
                        /**
                         * Apply group based polarization masking rule.
                         */
                        for (int index : ip11[i]) {
                            mask_local[index] = d11scale;
                        }
                        /**
                         * Loop over the neighbor list.
                         */
                        final int list[] = lists[i];
                        int npair = list.length;
                        counts[i] = 0;
                        preCounts[i] = 0;
                        final int ewald[] = ewalds[i];
                        int preList[] = preLists[i];
                        for (int j = 0; j < npair; j++) {
                            int k = list[j];
                            if (!use[k]
                                    || // No intermolecular electrostatics
                                    (lambdaMode == LambdaMode.VAPOR && intermolecularSoftcore
                                    && moleculei != molecule[k])) {
                                continue;
                            }
                            final double xk = x[k];
                            final double yk = y[k];
                            final double zk = z[k];
                            dx_local[0] = xk - xi;
                            dx_local[1] = yk - yi;
                            dx_local[2] = zk - zi;
                            final double r2 = crystal.image(dx_local);
                            if (r2 <= off2) {
                                count++;
                                ewald[counts[i]++] = k;
                                final double xr = dx_local[0];
                                final double yr = dx_local[1];
                                final double zr = dx_local[2];
                                final double pdk = ipdamp[k];
                                final double ptk = thole[k];
                                final double globalMultipolek[] = mpole[k];
                                final double ck = globalMultipolek[t000];
                                final double dkx = globalMultipolek[t100];
                                final double dky = globalMultipolek[t010];
                                final double dkz = globalMultipolek[t001];
                                final double qkxx = globalMultipolek[t200] * oneThird;
                                final double qkyy = globalMultipolek[t020] * oneThird;
                                final double qkzz = globalMultipolek[t002] * oneThird;
                                final double qkxy = globalMultipolek[t110] * oneThird;
                                final double qkxz = globalMultipolek[t101] * oneThird;
                                final double qkyz = globalMultipolek[t011] * oneThird;
                                double r = sqrt(r2);
                                if (r < preconditionerCutoff) {
                                    if (preList.length <= preCounts[i]) {
                                        int len = preList.length;
                                        preLists[i] = copyOf(preList, len + 10);
                                        preList = preLists[i];
                                    }
                                    preList[preCounts[i]++] = k;
                                }
                                /**
                                 * Calculate the error function damping terms.
                                 */
                                final double ralpha = aewald * r;
                                final double exp2a = exp(-ralpha * ralpha);
                                final double rr1 = 1.0 / r;
                                final double rr2 = rr1 * rr1;
                                final double bn0 = erfc(ralpha) * rr1;
                                final double bn1 = (bn0 + an0 * exp2a) * rr2;
                                final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                                final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;
                                /**
                                 * Compute the error function scaled and
                                 * unscaled terms.
                                 */
                                double scale3 = 1.0;
                                double scale5 = 1.0;
                                double scale7 = 1.0;
                                double damp = pdi * pdk;
                                //if (damp != 0.0) {
                                final double pgamma = min(pti, ptk);
                                final double rdamp = r * damp;
                                damp = -pgamma * rdamp * rdamp * rdamp;
                                if (damp > -50.0) {
                                    double expdamp = exp(damp);
                                    scale3 = 1.0 - expdamp;
                                    scale5 = 1.0 - expdamp * (1.0 - damp);
                                    scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                }
                                //}
                                final double scale = mask_local[k];
                                final double scalep = maskp_local[k];
                                final double dsc3 = scale3 * scale;
                                final double dsc5 = scale5 * scale;
                                final double dsc7 = scale7 * scale;
                                final double psc3 = scale3 * scalep;
                                final double psc5 = scale5 * scalep;
                                final double psc7 = scale7 * scalep;
                                final double rr3 = rr1 * rr2;
                                final double rr5 = 3.0 * rr3 * rr2;
                                final double rr7 = 5.0 * rr5 * rr2;
                                final double drr3 = (1.0 - dsc3) * rr3;
                                final double drr5 = (1.0 - dsc5) * rr5;
                                final double drr7 = (1.0 - dsc7) * rr7;
                                final double prr3 = (1.0 - psc3) * rr3;
                                final double prr5 = (1.0 - psc5) * rr5;
                                final double prr7 = (1.0 - psc7) * rr7;
                                final double dir = dix * xr + diy * yr + diz * zr;
                                final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                                final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                                final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                                final double qir = (qix * xr + qiy * yr + qiz * zr) * 0.5;
                                final double bn123i = bn1 * ci + bn2 * dir + bn3 * qir;
                                final double fkmx = xr * bn123i - bn1 * dix - bn2 * qix;
                                final double fkmy = yr * bn123i - bn1 * diy - bn2 * qiy;
                                final double fkmz = zr * bn123i - bn1 * diz - bn2 * qiz;
                                final double ddr357i = drr3 * ci + drr5 * dir + drr7 * qir;
                                final double fkdx = xr * ddr357i - drr3 * dix - drr5 * qix;
                                final double fkdy = yr * ddr357i - drr3 * diy - drr5 * qiy;
                                final double fkdz = zr * ddr357i - drr3 * diz - drr5 * qiz;
                                fX[k] += (fkmx - fkdx);
                                fY[k] += (fkmy - fkdy);
                                fZ[k] += (fkmz - fkdz);
                                final double prr357i = prr3 * ci + prr5 * dir + prr7 * qir;
                                final double fkpx = xr * prr357i - prr3 * dix - prr5 * qix;
                                final double fkpy = yr * prr357i - prr3 * diy - prr5 * qiy;
                                final double fkpz = zr * prr357i - prr3 * diz - prr5 * qiz;
                                fXCR[k] += (fkmx - fkpx);
                                fYCR[k] += (fkmy - fkpy);
                                fZCR[k] += (fkmz - fkpz);
                                final double dkr = dkx * xr + dky * yr + dkz * zr;
                                final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                                final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                                final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                                final double qkr = (qkx * xr + qky * yr + qkz * zr) * 0.5;
                                final double bn123k = bn1 * ck - bn2 * dkr + bn3 * qkr;
                                final double fimx = -xr * bn123k - bn1 * dkx + bn2 * qkx;
                                final double fimy = -yr * bn123k - bn1 * dky + bn2 * qky;
                                final double fimz = -zr * bn123k - bn1 * dkz + bn2 * qkz;
                                final double drr357k = drr3 * ck - drr5 * dkr + drr7 * qkr;
                                final double fidx = -xr * drr357k - drr3 * dkx + drr5 * qkx;
                                final double fidy = -yr * drr357k - drr3 * dky + drr5 * qky;
                                final double fidz = -zr * drr357k - drr3 * dkz + drr5 * qkz;
                                fX[i] += (fimx - fidx);
                                fY[i] += (fimy - fidy);
                                fZ[i] += (fimz - fidz);
                                final double prr357k = prr3 * ck - prr5 * dkr + prr7 * qkr;
                                final double fipx = -xr * prr357k - prr3 * dkx + prr5 * qkx;
                                final double fipy = -yr * prr357k - prr3 * dky + prr5 * qky;
                                final double fipz = -zr * prr357k - prr3 * dkz + prr5 * qkz;
                                fXCR[i] += (fimx - fipx);
                                fYCR[i] += (fimy - fipy);
                                fZCR[i] += (fimz - fipz);
                            }
                        }
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                if (index < 0) {
                                    ak.print();
                                }
                                maskp_local[index] = 1.0;
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                maskp_local[index] = 1.0;
                            }
                        }
                        for (Bond bond : ai.getBonds()) {
                            int index = bond.get1_2(ai).xyzIndex - 1;
                            maskp_local[index] = 1.0;
                        }
                        for (int index : ip11[i]) {
                            mask_local[index] = 1.0;
                        }
                    }
                    /**
                     * Loop over symmetry mates.
                     */
                    for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                        SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                        crystal.getTransformationOperator(symOp, transOp);
                        lists = neighborLists[iSymm];
                        ewalds = realSpaceLists[iSymm];
                        counts = realSpaceCounts[iSymm];
                        preLists = preconditionerLists[iSymm];
                        preCounts = preconditionerCounts[iSymm];
                        double xs[] = coordinates[iSymm][0];
                        double ys[] = coordinates[iSymm][1];
                        double zs[] = coordinates[iSymm][2];
                        double mpoles[][] = globalMultipole[iSymm];
                        /**
                         * Loop over atoms in a chunk of the asymmetric unit.
                         */
                        for (int i = lb; i <= ub; i++) {
                            if (!use[i]) {
                                continue;
                            }
                            final double pdi = ipdamp[i];
                            final double pti = thole[i];
                            final double multipolei[] = mpole[i];
                            final double ci = multipolei[t000];
                            final double dix = multipolei[t100];
                            final double diy = multipolei[t010];
                            final double diz = multipolei[t001];
                            final double qixx = multipolei[t200] * oneThird;
                            final double qiyy = multipolei[t020] * oneThird;
                            final double qizz = multipolei[t002] * oneThird;
                            final double qixy = multipolei[t110] * oneThird;
                            final double qixz = multipolei[t101] * oneThird;
                            final double qiyz = multipolei[t011] * oneThird;
                            final double xi = x[i];
                            final double yi = y[i];
                            final double zi = z[i];
                            /**
                             * Loop over the neighbor list.
                             */
                            final int list[] = lists[i];
                            final int npair = list.length;
                            counts[i] = 0;
                            preCounts[i] = 0;
                            final int ewald[] = ewalds[i];
                            final int preList[] = preLists[i];
                            for (int j = 0; j < npair; j++) {
                                int k = list[j];
                                if (!use[k]) {
                                    continue;
                                }
                                final double xk = xs[k];
                                final double yk = ys[k];
                                final double zk = zs[k];
                                dx_local[0] = xk - xi;
                                dx_local[1] = yk - yi;
                                dx_local[2] = zk - zi;
                                final double r2 = crystal.image(dx_local);
                                if (r2 <= off2) {
                                    count++;
                                    ewald[counts[i]++] = k;
                                    double selfScale = 1.0;
                                    if (i == k) {
                                        selfScale = 0.5;
                                    }
                                    final double xr = dx_local[0];
                                    final double yr = dx_local[1];
                                    final double zr = dx_local[2];
                                    final double pdk = ipdamp[k];
                                    final double ptk = thole[k];
                                    final double multipolek[] = mpoles[k];
                                    final double ck = multipolek[t000];
                                    final double dkx = multipolek[t100];
                                    final double dky = multipolek[t010];
                                    final double dkz = multipolek[t001];
                                    final double qkxx = multipolek[t200] * oneThird;
                                    final double qkyy = multipolek[t020] * oneThird;
                                    final double qkzz = multipolek[t002] * oneThird;
                                    final double qkxy = multipolek[t110] * oneThird;
                                    final double qkxz = multipolek[t101] * oneThird;
                                    final double qkyz = multipolek[t011] * oneThird;
                                    final double r = sqrt(r2);
                                    if (r < preconditionerCutoff) {
                                        preList[preCounts[i]++] = k;
                                    }
                                    /**
                                     * Calculate the error function damping
                                     * terms.
                                     */
                                    final double ralpha = aewald * r;
                                    final double exp2a = exp(-ralpha * ralpha);
                                    final double rr1 = 1.0 / r;
                                    final double rr2 = rr1 * rr1;
                                    final double bn0 = erfc(ralpha) * rr1;
                                    final double bn1 = (bn0 + an0 * exp2a) * rr2;
                                    final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                                    final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;
                                    /**
                                     * Compute the error function scaled and
                                     * unscaled terms.
                                     */
                                    double scale3 = 1.0;
                                    double scale5 = 1.0;
                                    double scale7 = 1.0;
                                    double damp = pdi * pdk;
                                    //if (damp != 0.0) {
                                    final double pgamma = min(pti, ptk);
                                    final double rdamp = r * damp;
                                    damp = -pgamma * rdamp * rdamp * rdamp;
                                    if (damp > -50.0) {
                                        double expdamp = exp(damp);
                                        scale3 = 1.0 - expdamp;
                                        scale5 = 1.0 - expdamp * (1.0 - damp);
                                        scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                    }
                                    //}
                                    final double dsc3 = scale3;
                                    final double dsc5 = scale5;
                                    final double dsc7 = scale7;
                                    final double rr3 = rr1 * rr2;
                                    final double rr5 = 3.0 * rr3 * rr2;
                                    final double rr7 = 5.0 * rr5 * rr2;
                                    final double drr3 = (1.0 - dsc3) * rr3;
                                    final double drr5 = (1.0 - dsc5) * rr5;
                                    final double drr7 = (1.0 - dsc7) * rr7;

                                    final double dkr = dkx * xr + dky * yr + dkz * zr;
                                    final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                                    final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                                    final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                                    final double qkr = (qkx * xr + qky * yr + qkz * zr) * 0.5;
                                    final double bn123k = bn1 * ck - bn2 * dkr + bn3 * qkr;
                                    final double drr357k = drr3 * ck - drr5 * dkr + drr7 * qkr;
                                    final double fimx = -xr * bn123k - bn1 * dkx + bn2 * qkx;
                                    final double fimy = -yr * bn123k - bn1 * dky + bn2 * qky;
                                    final double fimz = -zr * bn123k - bn1 * dkz + bn2 * qkz;
                                    final double fidx = -xr * drr357k - drr3 * dkx + drr5 * qkx;
                                    final double fidy = -yr * drr357k - drr3 * dky + drr5 * qky;
                                    final double fidz = -zr * drr357k - drr3 * dkz + drr5 * qkz;

                                    final double dir = dix * xr + diy * yr + diz * zr;
                                    final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                                    final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                                    final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                                    final double qir = (qix * xr + qiy * yr + qiz * zr) * 0.5;
                                    final double bn123i = bn1 * ci + bn2 * dir + bn3 * qir;
                                    final double ddr357i = drr3 * ci + drr5 * dir + drr7 * qir;
                                    final double fkmx = xr * bn123i - bn1 * dix - bn2 * qix;
                                    final double fkmy = yr * bn123i - bn1 * diy - bn2 * qiy;
                                    final double fkmz = zr * bn123i - bn1 * diz - bn2 * qiz;
                                    final double fkdx = xr * ddr357i - drr3 * dix - drr5 * qix;
                                    final double fkdy = yr * ddr357i - drr3 * diy - drr5 * qiy;
                                    final double fkdz = zr * ddr357i - drr3 * diz - drr5 * qiz;

                                    final double fix = selfScale * (fimx - fidx);
                                    final double fiy = selfScale * (fimy - fidy);
                                    final double fiz = selfScale * (fimz - fidz);
                                    fX[i] += fix;
                                    fY[i] += fiy;
                                    fZ[i] += fiz;
                                    fXCR[i] += fix;
                                    fYCR[i] += fiy;
                                    fZCR[i] += fiz;
                                    final double xc = selfScale * (fkmx - fkdx);
                                    final double yc = selfScale * (fkmy - fkdy);
                                    final double zc = selfScale * (fkmz - fkdz);
                                    final double fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
                                    final double fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
                                    final double fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
                                    fX[k] += fkx;
                                    fY[k] += fky;
                                    fZ[k] += fkz;
                                    fXCR[k] += fkx;
                                    fYCR[k] += fky;
                                    fZCR[k] += fkz;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * The Induced Dipole Field Region should be executed by a ParallelTeam with
     * exactly 2 threads. The Real Space and Reciprocal Space Sections will be
     * run concurrently, each with the number of threads defined by their
     * respective ParallelTeam instances.
     */
    private class InducedDipoleFieldRegion extends ParallelRegion {

        private InducedDipoleRealSpaceFieldSection inducedRealSpaceFieldSection;
        private InducedDipoleReciprocalFieldSection inducedReciprocalFieldSection;

        public InducedDipoleFieldRegion(ParallelTeam pt) {
            inducedRealSpaceFieldSection = new InducedDipoleRealSpaceFieldSection(pt);
            inducedReciprocalFieldSection = new InducedDipoleReciprocalFieldSection();
        }

        @Override
        public void run() {
            try {
                if (aewald > 0.0) {
                    execute(inducedRealSpaceFieldSection, inducedReciprocalFieldSection);
                } else {
                    execute(inducedRealSpaceFieldSection);
                }
            } catch (Exception e) {
                String message = "Fatal exception computing the induced dipole field.\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class InducedDipoleRealSpaceFieldSection extends ParallelSection {

            private final InducedDipoleRealSpaceFieldRegion polarizationRealSpaceFieldRegion;
            private final ParallelTeam pt;

            public InducedDipoleRealSpaceFieldSection(ParallelTeam pt) {
                this.pt = pt;
                int nt = pt.getThreadCount();
                polarizationRealSpaceFieldRegion = new InducedDipoleRealSpaceFieldRegion(nt);
            }

            @Override
            public void run() {
                try {
                    realSpaceSCFTotal -= System.nanoTime();
                    pt.execute(polarizationRealSpaceFieldRegion);
                    realSpaceSCFTotal += System.nanoTime();
                } catch (Exception e) {
                    String message = "Fatal exception computing the real space field.\n";
                    logger.log(Level.SEVERE, message, e);
                }
            }
        }

        private class InducedDipoleReciprocalFieldSection extends ParallelSection {

            @Override
            public void run() {
                reciprocalSpace.inducedDipoleConvolution();
            }
        }

        private class InducedDipoleRealSpaceFieldRegion extends ParallelRegion {

            private final InducedRealSpaceFieldLoop inducedRealSpaceFieldLoop[];

            public InducedDipoleRealSpaceFieldRegion(int threadCount) {
                inducedRealSpaceFieldLoop = new InducedRealSpaceFieldLoop[threadCount];
            }

            @Override
            public void run() {
                int threadIndex = getThreadIndex();
                if (inducedRealSpaceFieldLoop[threadIndex] == null) {
                    inducedRealSpaceFieldLoop[threadIndex] = new InducedRealSpaceFieldLoop();
                }
                try {
                    execute(0, nAtoms - 1, inducedRealSpaceFieldLoop[threadIndex]);
                } catch (Exception e) {
                    String message = "Fatal exception computing the induced real space field in thread " + getThreadIndex() + "\n";
                    logger.log(Level.SEVERE, message, e);
                }
            }

            private class InducedRealSpaceFieldLoop extends IntegerForLoop {

                private double ind[][], indCR[][];
                private double x[], y[], z[];
                private double fX[], fY[], fZ[];
                private double fXCR[], fYCR[], fZCR[];

                public InducedRealSpaceFieldLoop() {
                }

                @Override
                public IntegerSchedule schedule() {
                    return realSpaceSchedule;
                }

                @Override
                public void start() {
                    int threadIndex = getThreadIndex();
                    realSpaceSCFTime[threadIndex] -= System.nanoTime();
                    fX = field[threadIndex][0];
                    fY = field[threadIndex][1];
                    fZ = field[threadIndex][2];
                    fXCR = fieldCR[threadIndex][0];
                    fYCR = fieldCR[threadIndex][1];
                    fZCR = fieldCR[threadIndex][2];
                    fill(fX, 0.0);
                    fill(fY, 0.0);
                    fill(fZ, 0.0);
                    fill(fXCR, 0.0);
                    fill(fYCR, 0.0);
                    fill(fZCR, 0.0);
                    x = coordinates[0][0];
                    y = coordinates[0][1];
                    z = coordinates[0][2];
                    ind = inducedDipole[0];
                    indCR = inducedDipoleCR[0];
                }

                @Override
                public void finish() {
                    int threadIndex = getThreadIndex();
                    realSpaceSCFTime[threadIndex] += System.nanoTime();
                }

                @Override
                public void run(int lb, int ub) {
                    final double dx[] = new double[3];
                    final double transOp[][] = new double[3][3];
                    /**
                     * Loop over a chunk of atoms.
                     */
                    int lists[][] = realSpaceLists[0];
                    int counts[] = realSpaceCounts[0];
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        final int moleculei = molecule[i];
                        double fx = 0.0;
                        double fy = 0.0;
                        double fz = 0.0;
                        double px = 0.0;
                        double py = 0.0;
                        double pz = 0.0;
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final double dipolei[] = ind[i];
                        final double uix = dipolei[0];
                        final double uiy = dipolei[1];
                        final double uiz = dipolei[2];
                        final double dipoleCRi[] = indCR[i];
                        final double pix = dipoleCRi[0];
                        final double piy = dipoleCRi[1];
                        final double piz = dipoleCRi[2];
                        final double pdi = ipdamp[i];
                        final double pti = thole[i];
                        /**
                         * Loop over the neighbor list.
                         */
                        final int list[] = lists[i];
                        final int npair = counts[i];
                        for (int j = 0; j < npair; j++) {
                            final int k = list[j];
                            if (!use[k]
                                    || // No intermolecular electrostatics
                                    (lambdaMode == LambdaMode.VAPOR && intermolecularSoftcore
                                    && moleculei != molecule[k])) {
                                continue;
                            }
                            final double pdk = ipdamp[k];
                            final double ptk = thole[k];
                            dx[0] = x[k] - xi;
                            dx[1] = y[k] - yi;
                            dx[2] = z[k] - zi;
                            final double r2 = crystal.image(dx);
                            /**
                             * Calculate the error function damping terms.
                             */
                            final double r = sqrt(r2);
                            final double rr1 = 1.0 / r;
                            final double rr2 = rr1 * rr1;
                            final double ralpha = aewald * r;
                            final double exp2a = exp(-ralpha * ralpha);
                            final double bn0 = erfc(ralpha) * rr1;
                            final double bn1 = (bn0 + an0 * exp2a) * rr2;
                            final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                            double scale3 = 1.0;
                            double scale5 = 1.0;
                            double damp = pdi * pdk;
                            //if (damp != 0.0) {
                            final double pgamma = min(pti, ptk);
                            final double rdamp = r * damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                final double expdamp = exp(damp);
                                scale3 = 1.0 - expdamp;
                                scale5 = 1.0 - expdamp * (1.0 - damp);
                            }
                            //}
                            double rr3 = rr1 * rr2;
                            double rr5 = 3.0 * rr3 * rr2;
                            rr3 *= (1.0 - scale3);
                            rr5 *= (1.0 - scale5);
                            final double xr = dx[0];
                            final double yr = dx[1];
                            final double zr = dx[2];
                            final double dipolek[] = ind[k];
                            final double ukx = dipolek[0];
                            final double uky = dipolek[1];
                            final double ukz = dipolek[2];
                            final double ukr = ukx * xr + uky * yr + ukz * zr;
                            final double bn2ukr = bn2 * ukr;
                            final double fimx = -bn1 * ukx + bn2ukr * xr;
                            final double fimy = -bn1 * uky + bn2ukr * yr;
                            final double fimz = -bn1 * ukz + bn2ukr * zr;
                            final double rr5ukr = rr5 * ukr;
                            final double fidx = -rr3 * ukx + rr5ukr * xr;
                            final double fidy = -rr3 * uky + rr5ukr * yr;
                            final double fidz = -rr3 * ukz + rr5ukr * zr;
                            fx += (fimx - fidx);
                            fy += (fimy - fidy);
                            fz += (fimz - fidz);
                            final double dipolepk[] = indCR[k];
                            final double pkx = dipolepk[0];
                            final double pky = dipolepk[1];
                            final double pkz = dipolepk[2];
                            final double pkr = pkx * xr + pky * yr + pkz * zr;
                            final double bn2pkr = bn2 * pkr;
                            final double pimx = -bn1 * pkx + bn2pkr * xr;
                            final double pimy = -bn1 * pky + bn2pkr * yr;
                            final double pimz = -bn1 * pkz + bn2pkr * zr;
                            final double rr5pkr = rr5 * pkr;
                            final double pidx = -rr3 * pkx + rr5pkr * xr;
                            final double pidy = -rr3 * pky + rr5pkr * yr;
                            final double pidz = -rr3 * pkz + rr5pkr * zr;
                            px += (pimx - pidx);
                            py += (pimy - pidy);
                            pz += (pimz - pidz);
                            final double uir = uix * xr + uiy * yr + uiz * zr;
                            final double bn2uir = bn2 * uir;
                            final double fkmx = -bn1 * uix + bn2uir * xr;
                            final double fkmy = -bn1 * uiy + bn2uir * yr;
                            final double fkmz = -bn1 * uiz + bn2uir * zr;
                            final double rr5uir = rr5 * uir;
                            final double fkdx = -rr3 * uix + rr5uir * xr;
                            final double fkdy = -rr3 * uiy + rr5uir * yr;
                            final double fkdz = -rr3 * uiz + rr5uir * zr;
                            fX[k] += (fkmx - fkdx);
                            fY[k] += (fkmy - fkdy);
                            fZ[k] += (fkmz - fkdz);
                            final double pir = pix * xr + piy * yr + piz * zr;
                            final double bn2pir = bn2 * pir;
                            final double pkmx = -bn1 * pix + bn2pir * xr;
                            final double pkmy = -bn1 * piy + bn2pir * yr;
                            final double pkmz = -bn1 * piz + bn2pir * zr;
                            final double rr5pir = rr5 * pir;
                            final double pkdx = -rr3 * pix + rr5pir * xr;
                            final double pkdy = -rr3 * piy + rr5pir * yr;
                            final double pkdz = -rr3 * piz + rr5pir * zr;
                            fXCR[k] += (pkmx - pkdx);
                            fYCR[k] += (pkmy - pkdy);
                            fZCR[k] += (pkmz - pkdz);
                        }
                        fX[i] += fx;
                        fY[i] += fy;
                        fZ[i] += fz;
                        fXCR[i] += px;
                        fYCR[i] += py;
                        fZCR[i] += pz;
                    }
                    /**
                     * Loop over symmetry mates.
                     */
                    for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                        SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                        crystal.getTransformationOperator(symOp, transOp);
                        lists = realSpaceLists[iSymm];
                        counts = realSpaceCounts[iSymm];
                        final double xs[] = coordinates[iSymm][0];
                        final double ys[] = coordinates[iSymm][1];
                        final double zs[] = coordinates[iSymm][2];
                        final double inds[][] = inducedDipole[iSymm];
                        final double indCRs[][] = inducedDipoleCR[iSymm];
                        /**
                         * Loop over a chunk of atoms.
                         */
                        for (int i = lb; i <= ub; i++) {
                            if (!use[i]) {
                                continue;
                            }
                            double fx = 0.0;
                            double fy = 0.0;
                            double fz = 0.0;
                            double px = 0.0;
                            double py = 0.0;
                            double pz = 0.0;
                            final double xi = x[i];
                            final double yi = y[i];
                            final double zi = z[i];
                            final double dipolei[] = ind[i];
                            final double uix = dipolei[0];
                            final double uiy = dipolei[1];
                            final double uiz = dipolei[2];
                            final double dipoleCRi[] = indCR[i];
                            final double pix = dipoleCRi[0];
                            final double piy = dipoleCRi[1];
                            final double piz = dipoleCRi[2];
                            final double pdi = ipdamp[i];
                            final double pti = thole[i];
                            /**
                             * Loop over the neighbor list.
                             */
                            final int list[] = lists[i];
                            final int npair = counts[i];
                            for (int j = 0; j < npair; j++) {
                                final int k = list[j];
                                if (!use[k]) {
                                    continue;
                                }
                                double selfScale = 1.0;
                                if (i == k) {
                                    selfScale = 0.5;
                                }
                                final double pdk = ipdamp[k];
                                final double ptk = thole[k];
                                dx[0] = xs[k] - xi;
                                dx[1] = ys[k] - yi;
                                dx[2] = zs[k] - zi;
                                final double r2 = crystal.image(dx);
                                /**
                                 * Calculate the error function damping terms.
                                 */
                                final double r = sqrt(r2);
                                final double rr1 = 1.0 / r;
                                final double rr2 = rr1 * rr1;
                                final double ralpha = aewald * r;
                                final double exp2a = exp(-ralpha * ralpha);
                                final double bn0 = erfc(ralpha) * rr1;
                                final double bn1 = (bn0 + an0 * exp2a) * rr2;
                                final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                                double scale3 = 1.0;
                                double scale5 = 1.0;
                                double damp = pdi * pdk;
                                //if (damp != 0.0) {
                                final double pgamma = min(pti, ptk);
                                final double rdamp = r * damp;
                                damp = -pgamma * rdamp * rdamp * rdamp;
                                if (damp > -50.0) {
                                    final double expdamp = exp(damp);
                                    scale3 = 1.0 - expdamp;
                                    scale5 = 1.0 - expdamp * (1.0 - damp);
                                }
                                //}
                                double rr3 = rr1 * rr2;
                                double rr5 = 3.0 * rr3 * rr2;
                                rr3 *= (1.0 - scale3);
                                rr5 *= (1.0 - scale5);
                                final double xr = dx[0];
                                final double yr = dx[1];
                                final double zr = dx[2];
                                final double dipolek[] = inds[k];
                                final double ukx = dipolek[0];
                                final double uky = dipolek[1];
                                final double ukz = dipolek[2];
                                final double dipolepk[] = indCRs[k];
                                final double pkx = dipolepk[0];
                                final double pky = dipolepk[1];
                                final double pkz = dipolepk[2];
                                final double ukr = ukx * xr + uky * yr + ukz * zr;
                                final double bn2ukr = bn2 * ukr;
                                final double fimx = -bn1 * ukx + bn2ukr * xr;
                                final double fimy = -bn1 * uky + bn2ukr * yr;
                                final double fimz = -bn1 * ukz + bn2ukr * zr;
                                final double rr5ukr = rr5 * ukr;
                                final double fidx = -rr3 * ukx + rr5ukr * xr;
                                final double fidy = -rr3 * uky + rr5ukr * yr;
                                final double fidz = -rr3 * ukz + rr5ukr * zr;
                                fx += selfScale * (fimx - fidx);
                                fy += selfScale * (fimy - fidy);
                                fz += selfScale * (fimz - fidz);
                                final double pkr = pkx * xr + pky * yr + pkz * zr;
                                final double bn2pkr = bn2 * pkr;
                                final double pimx = -bn1 * pkx + bn2pkr * xr;
                                final double pimy = -bn1 * pky + bn2pkr * yr;
                                final double pimz = -bn1 * pkz + bn2pkr * zr;
                                final double rr5pkr = rr5 * pkr;
                                final double pidx = -rr3 * pkx + rr5pkr * xr;
                                final double pidy = -rr3 * pky + rr5pkr * yr;
                                final double pidz = -rr3 * pkz + rr5pkr * zr;
                                px += selfScale * (pimx - pidx);
                                py += selfScale * (pimy - pidy);
                                pz += selfScale * (pimz - pidz);
                                final double uir = uix * xr + uiy * yr + uiz * zr;
                                final double bn2uir = bn2 * uir;
                                final double fkmx = -bn1 * uix + bn2uir * xr;
                                final double fkmy = -bn1 * uiy + bn2uir * yr;
                                final double fkmz = -bn1 * uiz + bn2uir * zr;
                                final double rr5uir = rr5 * uir;
                                final double fkdx = -rr3 * uix + rr5uir * xr;
                                final double fkdy = -rr3 * uiy + rr5uir * yr;
                                final double fkdz = -rr3 * uiz + rr5uir * zr;
                                double xc = selfScale * (fkmx - fkdx);
                                double yc = selfScale * (fkmy - fkdy);
                                double zc = selfScale * (fkmz - fkdz);
                                fX[k] += (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                                fY[k] += (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                                fZ[k] += (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                                final double pir = pix * xr + piy * yr + piz * zr;
                                final double bn2pir = bn2 * pir;
                                final double pkmx = -bn1 * pix + bn2pir * xr;
                                final double pkmy = -bn1 * piy + bn2pir * yr;
                                final double pkmz = -bn1 * piz + bn2pir * zr;
                                final double rr5pir = rr5 * pir;
                                final double pkdx = -rr3 * pix + rr5pir * xr;
                                final double pkdy = -rr3 * piy + rr5pir * yr;
                                final double pkdz = -rr3 * piz + rr5pir * zr;
                                xc = selfScale * (pkmx - pkdx);
                                yc = selfScale * (pkmy - pkdy);
                                zc = selfScale * (pkmz - pkdz);
                                fXCR[k] += (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                                fYCR[k] += (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                                fZCR[k] += (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                            }
                            fX[i] += fx;
                            fY[i] += fy;
                            fZ[i] += fz;
                            fXCR[i] += px;
                            fYCR[i] += py;
                            fZCR[i] += pz;
                        }
                    }
                }
            }
        }
    }

    private class DirectRegion extends ParallelRegion {

        private final DirectLoop directLoop[];

        public DirectRegion(int nt) {
            directLoop = new DirectLoop[nt];
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();
            if (directLoop[ti] == null) {
                directLoop[ti] = new DirectLoop();
            }
            try {
                execute(0, nAtoms - 1, directLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the direct induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class DirectLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                /**
                 * Reduce the direct field.
                 */
                for (int i = lb; i <= ub; i++) {
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double fxCR = 0.0;
                    double fyCR = 0.0;
                    double fzCR = 0.0;
                    for (int j = 1; j < maxThreads; j++) {
                        fx += field[j][0][i];
                        fy += field[j][1][i];
                        fz += field[j][2][i];
                        fxCR += fieldCR[j][0][i];
                        fyCR += fieldCR[j][1][i];
                        fzCR += fieldCR[j][2][i];
                    }
                    field[0][0][i] += fx;
                    field[0][1][i] += fy;
                    field[0][2][i] += fz;
                    fieldCR[0][0][i] += fxCR;
                    fieldCR[0][1][i] += fyCR;
                    fieldCR[0][2][i] += fzCR;
                }
                if (aewald > 0.0) {
                    /**
                     * Add the self and reciprocal space contributions.
                     */
                    for (int i = lb; i <= ub; i++) {
                        double mpolei[] = globalMultipole[0][i];
                        double phii[] = cartMultipolePhi[i];
                        double fx = aewald3 * mpolei[t100] - phii[t100];
                        double fy = aewald3 * mpolei[t010] - phii[t010];
                        double fz = aewald3 * mpolei[t001] - phii[t001];
                        field[0][0][i] += fx;
                        field[0][1][i] += fy;
                        field[0][2][i] += fz;
                        fieldCR[0][0][i] += fx;
                        fieldCR[0][1][i] += fy;
                        fieldCR[0][2][i] += fz;
                    }
                }
                if (generalizedKirkwoodTerm) {
                    /**
                     * Initialize the electric field to the direct field plus
                     * the permanent GK reaction field.
                     */
                    SharedDoubleArray gkField[] = generalizedKirkwood.sharedGKField;
                    for (int i = lb; i <= ub; i++) {
                        double fx = gkField[0].get(i);
                        double fy = gkField[1].get(i);
                        double fz = gkField[2].get(i);
                        field[0][0][i] += fx;
                        field[0][1][i] += fy;
                        field[0][2][i] += fz;
                        fieldCR[0][0][i] += fx;
                        fieldCR[0][1][i] += fy;
                        fieldCR[0][2][i] += fz;
                    }
                }
                /**
                 * Set the direct induced dipoles to the polarizability
                 * multiplied by the direct field.
                 */
                final double induced0[][] = inducedDipole[0];
                final double inducedCR0[][] = inducedDipoleCR[0];
                for (int i = lb; i <= ub; i++) {
                    final double polar = polarizability[i];
                    final double ind[] = induced0[i];
                    final double directi[] = directDipole[i];
                    ind[0] = polar * field[0][0][i];
                    ind[1] = polar * field[0][1][i];
                    ind[2] = polar * field[0][2][i];
                    directi[0] = ind[0];
                    directi[1] = ind[1];
                    directi[2] = ind[2];
                    final double indCR[] = inducedCR0[i];
                    final double directCRi[] = directDipoleCR[i];
                    indCR[0] = polar * fieldCR[0][0][i];
                    indCR[1] = polar * fieldCR[0][1][i];
                    indCR[2] = polar * fieldCR[0][2][i];
                    directCRi[0] = indCR[0];
                    directCRi[1] = indCR[1];
                    directCRi[2] = indCR[2];
                }
            }
        }
    }

    private class SORRegion extends ParallelRegion {

        private final SORLoop sorLoop[];
        private final SharedDouble sharedEps;
        private final SharedDouble sharedEpsCR;

        public SORRegion(int nt) {
            sorLoop = new SORLoop[nt];
            sharedEps = new SharedDouble();
            sharedEpsCR = new SharedDouble();
        }

        public double getEps() {
            double eps = sharedEps.get();
            double epsCR = sharedEpsCR.get();
            return max(eps, epsCR);
        }

        @Override
        public void start() {
            sharedEps.set(0.0);
            sharedEpsCR.set(0.0);
        }

        @Override
        public void run() throws Exception {
            try {
                int ti = getThreadIndex();
                if (sorLoop[ti] == null) {
                    sorLoop[ti] = new SORLoop();
                }
                execute(0, nAtoms - 1, sorLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class SORLoop extends IntegerForLoop {

            private double eps, epsCR;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eps = 0.0;
                epsCR = 0.0;
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                final double induced0[][] = inducedDipole[0];
                final double inducedCR0[][] = inducedDipoleCR[0];
                /**
                 * Reduce the real space field.
                 */
                for (int i = lb; i <= ub; i++) {
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double fxCR = 0.0;
                    double fyCR = 0.0;
                    double fzCR = 0.0;
                    for (int j = 1; j < maxThreads; j++) {
                        fx += field[j][0][i];
                        fy += field[j][1][i];
                        fz += field[j][2][i];
                        fxCR += fieldCR[j][0][i];
                        fyCR += fieldCR[j][1][i];
                        fzCR += fieldCR[j][2][i];
                    }
                    field[0][0][i] += fx;
                    field[0][1][i] += fy;
                    field[0][2][i] += fz;
                    fieldCR[0][0][i] += fxCR;
                    fieldCR[0][1][i] += fyCR;
                    fieldCR[0][2][i] += fzCR;
                }
                if (aewald > 0.0) {
                    /**
                     * Add the self and reciprocal space fields to the real
                     * space field.
                     */
                    for (int i = lb; i <= ub; i++) {
                        double dipolei[] = induced0[i];
                        double dipoleCRi[] = inducedCR0[i];
                        final double phii[] = cartesianDipolePhi[i];
                        final double phiCRi[] = cartesianDipolePhiCR[i];
                        double fx = aewald3 * dipolei[0] - phii[t100];
                        double fy = aewald3 * dipolei[1] - phii[t010];
                        double fz = aewald3 * dipolei[2] - phii[t001];
                        double fxCR = aewald3 * dipoleCRi[0] - phiCRi[t100];
                        double fyCR = aewald3 * dipoleCRi[1] - phiCRi[t010];
                        double fzCR = aewald3 * dipoleCRi[2] - phiCRi[t001];
                        field[0][0][i] += fx;
                        field[0][1][i] += fy;
                        field[0][2][i] += fz;
                        fieldCR[0][0][i] += fxCR;
                        fieldCR[0][1][i] += fyCR;
                        fieldCR[0][2][i] += fzCR;
                    }
                }
                if (generalizedKirkwoodTerm) {
                    SharedDoubleArray gkField[] = generalizedKirkwood.sharedGKField;
                    SharedDoubleArray gkFieldCR[] = generalizedKirkwood.sharedGKFieldCR;
                    /**
                     * Add the GK reaction field to the intramolecular field.
                     */
                    for (int i = lb; i <= ub; i++) {
                        field[0][0][i] += gkField[0].get(i);
                        field[0][1][i] += gkField[1].get(i);
                        field[0][2][i] += gkField[2].get(i);
                        fieldCR[0][0][i] += gkFieldCR[0].get(i);
                        fieldCR[0][1][i] += gkFieldCR[1].get(i);
                        fieldCR[0][2][i] += gkFieldCR[2].get(i);
                    }
                }

                /**
                 * Apply Successive Over-Relaxation (SOR).
                 */
                for (int i = lb; i <= ub; i++) {
                    final double ind[] = induced0[i];
                    final double indCR[] = inducedCR0[i];
                    final double direct[] = directDipole[i];
                    final double directCR[] = directDipoleCR[i];
                    final double polar = polarizability[i];
                    for (int j = 0; j < 3; j++) {
                        double previous = ind[j];
                        double mutual = polar * field[0][j][i];
                        ind[j] = direct[j] + mutual;
                        double delta = polsor * (ind[j] - previous);
                        ind[j] = previous + delta;
                        eps += delta * delta;
                        previous = indCR[j];
                        mutual = polar * fieldCR[0][j][i];
                        indCR[j] = directCR[j] + mutual;
                        delta = polsor * (indCR[j] - previous);
                        indCR[j] = previous + delta;
                        epsCR += delta * delta;
                    }
                }
            }

            @Override
            public void finish() {
                sharedEps.addAndGet(eps);
                sharedEpsCR.addAndGet(epsCR);
            }
        }
    }

    /**
     * The Real Space Energy Region class parallelizes evaluation of the real
     * space energy and gradient.
     */
    private class RealSpaceEnergyRegion extends ParallelRegion {

        private double permanentEnergy;
        private double polarizationEnergy;
        private final SharedInteger sharedInteractions;
        private final RealSpaceEnergyLoop realSpaceEnergyLoop[];

        public RealSpaceEnergyRegion(int nt) {
            sharedInteractions = new SharedInteger();
            realSpaceEnergyLoop = new RealSpaceEnergyLoop[nt];
        }

        public double getPermanentEnergy() {
            return permanentEnergy;
        }

        public double getPolarizationEnergy() {
            return polarizationEnergy;
        }

        public int getInteractions() {
            return sharedInteractions.get();
        }

        @Override
        public void start() {
            sharedInteractions.set(0);
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            if (realSpaceEnergyLoop[threadIndex] == null) {
                realSpaceEnergyLoop[threadIndex] = new RealSpaceEnergyLoop();
            }
            try {
                execute(0, nAtoms - 1, realSpaceEnergyLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing the real space energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        @Override
        public void finish() {
            permanentEnergy = 0.0;
            polarizationEnergy = 0.0;
            for (int i = 0; i < maxThreads; i++) {
                double e = realSpaceEnergyLoop[i].permanentEnergy;
                if (Double.isNaN(e)) {
                    logger.severe(String.format(" The permanent multipole energy of thread %d is %16.8f", i, e));
                }
                permanentEnergy += e;
                double ei = realSpaceEnergyLoop[i].inducedEnergy;
                if (Double.isNaN(ei)) {
                    logger.severe(String.format(" The polarization energy of thread %d is %16.8f", i, ei));
                }
                polarizationEnergy += ei;
            }
            permanentEnergy *= ELECTRIC;
            polarizationEnergy *= ELECTRIC;
        }

        /**
         * The Real Space Gradient Loop class contains methods and thread local
         * variables to parallelize the evaluation of the real space permanent
         * and polarization energies and gradients.
         */
        private class RealSpaceEnergyLoop extends IntegerForLoop {

            private double ci;
            private double dix, diy, diz;
            private double qixx, qiyy, qizz, qixy, qixz, qiyz;
            private double ck;
            private double dkx, dky, dkz;
            private double qkxx, qkyy, qkzz, qkxy, qkxz, qkyz;
            private double uix, uiy, uiz;
            private double pix, piy, piz;
            private double xr, yr, zr;
            private double ukx, uky, ukz;
            private double pkx, pky, pkz;
            private double bn0, bn1, bn2, bn3, bn4, bn5, bn6;
            private double r2, rr1, rr2, rr3, rr5, rr7, rr9, rr11, rr13;
            private double scale, scale3, scale5, scale7;
            private double scalep, scaled;
            private double ddsc3x, ddsc3y, ddsc3z;
            private double ddsc5x, ddsc5y, ddsc5z;
            private double ddsc7x, ddsc7y, ddsc7z;
            private double beta, l2;
            private boolean soft;
            private double selfScale;
            private double permanentEnergy;
            private double inducedEnergy;
            private double dUdL, d2UdL2;
            private int i, k, iSymm, count;
            private SymOp symOp;
            private double gX[], gY[], gZ[], tX[], tY[], tZ[];
            private double lgX[], lgY[], lgZ[], ltX[], ltY[], ltZ[];
            private double gxk_local[], gyk_local[], gzk_local[];
            private double txk_local[], tyk_local[], tzk_local[];
            private double lxk_local[], lyk_local[], lzk_local[];
            private double ltxk_local[], ltyk_local[], ltzk_local[];
            private double masking_local[];
            private double maskingp_local[];
            private double maskingd_local[];
            private final double dx_local[];
            private final double rot_local[][];
            private final double work[][];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public RealSpaceEnergyLoop() {
                super();
                dx_local = new double[3];
                work = new double[15][3];
                rot_local = new double[3][3];
            }

            private void init() {
                if (masking_local == null || masking_local.length < nAtoms) {
                    txk_local = new double[nAtoms];
                    tyk_local = new double[nAtoms];
                    tzk_local = new double[nAtoms];
                    gxk_local = new double[nAtoms];
                    gyk_local = new double[nAtoms];
                    gzk_local = new double[nAtoms];
                    lxk_local = new double[nAtoms];
                    lyk_local = new double[nAtoms];
                    lzk_local = new double[nAtoms];
                    ltxk_local = new double[nAtoms];
                    ltyk_local = new double[nAtoms];
                    ltzk_local = new double[nAtoms];
                    masking_local = new double[nAtoms];
                    maskingp_local = new double[nAtoms];
                    maskingd_local = new double[nAtoms];
                    fill(masking_local, 1.0);
                    fill(maskingp_local, 1.0);
                    fill(maskingd_local, 1.0);
                }
            }

            @Override
            public IntegerSchedule schedule() {
                return realSpaceSchedule;
            }

            @Override
            public void start() {
                init();
                int threadIndex = getThreadIndex();
                realSpaceEnergyTime[threadIndex] -= System.nanoTime();
                permanentEnergy = 0.0;
                inducedEnergy = 0.0;
                count = 0;
                gX = grad[threadIndex][0];
                gY = grad[threadIndex][1];
                gZ = grad[threadIndex][2];
                tX = torque[threadIndex][0];
                tY = torque[threadIndex][1];
                tZ = torque[threadIndex][2];
                if (lambdaTerm) {
                    dUdL = 0.0;
                    d2UdL2 = 0.0;
                    lgX = lambdaGrad[threadIndex][0];
                    lgY = lambdaGrad[threadIndex][1];
                    lgZ = lambdaGrad[threadIndex][2];
                    ltX = lambdaTorque[threadIndex][0];
                    ltY = lambdaTorque[threadIndex][1];
                    ltZ = lambdaTorque[threadIndex][2];
                }
            }

            @Override
            public void run(int lb, int ub) {
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (iSymm = 0; iSymm < nSymm; iSymm++) {
                    symOp = symOps.get(iSymm);
                    if (gradient) {
                        fill(gxk_local, 0.0);
                        fill(gyk_local, 0.0);
                        fill(gzk_local, 0.0);
                        fill(txk_local, 0.0);
                        fill(tyk_local, 0.0);
                        fill(tzk_local, 0.0);
                    }
                    if (lambdaTerm) {
                        fill(lxk_local, 0.0);
                        fill(lyk_local, 0.0);
                        fill(lzk_local, 0.0);
                        fill(ltxk_local, 0.0);
                        fill(ltyk_local, 0.0);
                        fill(ltzk_local, 0.0);
                    }
                    realSpaceChunk(lb, ub);
                    if (gradient) {
                        // Turn symmetry mate torques into gradients
                        torque(iSymm, txk_local, tyk_local, tzk_local,
                                gxk_local, gyk_local, gzk_local,
                                work[0], work[1], work[2], work[3], work[4],
                                work[5], work[6], work[7], work[8], work[9],
                                work[10], work[11], work[12], work[13], work[14]);
                        // Rotate symmetry mate gradients
                        if (iSymm != 0) {
                            crystal.applyTransSymRot(nAtoms,
                                    gxk_local, gyk_local, gzk_local,
                                    gxk_local, gyk_local, gzk_local,
                                    symOp, rot_local);
                        }
                        // Sum symmetry mate gradients into asymmetric unit gradients
                        for (int j = 0; j < nAtoms; j++) {
                            gX[j] += gxk_local[j];
                            gY[j] += gyk_local[j];
                            gZ[j] += gzk_local[j];
                        }
                    }
                    if (lambdaTerm) {
                        // Turn symmetry mate torques into gradients
                        torque(iSymm, ltxk_local, ltyk_local, ltzk_local,
                                lxk_local, lyk_local, lzk_local,
                                work[0], work[1], work[2], work[3], work[4],
                                work[5], work[6], work[7], work[8], work[9],
                                work[10], work[11], work[12], work[13], work[14]);
                        // Rotate symmetry mate gradients
                        if (iSymm != 0) {
                            crystal.applyTransSymRot(nAtoms, lxk_local, lyk_local, lzk_local,
                                    lxk_local, lyk_local, lzk_local, symOp, rot_local);
                        }
                        // Sum symmetry mate gradients into asymmetric unit gradients
                        for (int j = 0; j < nAtoms; j++) {
                            lgX[j] += lxk_local[j];
                            lgY[j] += lyk_local[j];
                            lgZ[j] += lzk_local[j];
                        }
                    }

                }
            }

            public int getCount() {
                return count;
            }

            @Override
            public void finish() {
                sharedInteractions.addAndGet(count);
                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(dUdL * ELECTRIC);
                    sharedd2EdLambda2.addAndGet(d2UdL2 * ELECTRIC);
                }
                realSpaceEnergyTime[getThreadIndex()] += System.nanoTime();
            }

            /**
             * Evaluate the real space permanent energy and polarization energy
             * for a chunk of atoms.
             *
             * @param lb The lower bound of the chunk.
             * @param ub The upper bound of the chunk.
             */
            private void realSpaceChunk(final int lb, final int ub) {
                final double x[] = coordinates[0][0];
                final double y[] = coordinates[0][1];
                final double z[] = coordinates[0][2];
                final double mpole[][] = globalMultipole[0];
                final double ind[][] = inducedDipole[0];
                final double indp[][] = inducedDipoleCR[0];
                final int lists[][] = realSpaceLists[iSymm];
                final double neighborX[] = coordinates[iSymm][0];
                final double neighborY[] = coordinates[iSymm][1];
                final double neighborZ[] = coordinates[iSymm][2];
                final double neighborMultipole[][] = globalMultipole[iSymm];
                final double neighborInducedDipole[][] = inducedDipole[iSymm];
                final double neighborInducedDipolep[][] = inducedDipoleCR[iSymm];
                for (i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final Atom ai = atoms[i];
                    final int moleculei = molecule[i];
                    if (iSymm == 0) {
                        for (Atom ak : ai.get1_5s()) {
                            masking_local[ak.xyzIndex - 1] = m15scale;
                        }
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = m14scale;
                                for (int j : ip11[i]) {
                                    if (j == index) {
                                        maskingp_local[index] = 0.5;
                                    }
                                }
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = m13scale;
                                maskingp_local[index] = p13scale;
                            }
                        }
                        for (Bond bond : ai.getBonds()) {
                            int index = bond.get1_2(ai).xyzIndex - 1;
                            masking_local[index] = m12scale;
                            maskingp_local[index] = p12scale;
                        }
                        for (int j : ip11[i]) {
                            maskingd_local[j] = d11scale;
                        }
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double globalMultipolei[] = mpole[i];
                    final double inducedDipolei[] = ind[i];
                    final double inducedDipolepi[] = indp[i];
                    ci = globalMultipolei[t000];
                    dix = globalMultipolei[t100];
                    diy = globalMultipolei[t010];
                    diz = globalMultipolei[t001];
                    qixx = globalMultipolei[t200] * oneThird;
                    qiyy = globalMultipolei[t020] * oneThird;
                    qizz = globalMultipolei[t002] * oneThird;
                    qixy = globalMultipolei[t110] * oneThird;
                    qixz = globalMultipolei[t101] * oneThird;
                    qiyz = globalMultipolei[t011] * oneThird;
                    uix = inducedDipolei[0];
                    uiy = inducedDipolei[1];
                    uiz = inducedDipolei[2];
                    pix = inducedDipolepi[0];
                    piy = inducedDipolepi[1];
                    piz = inducedDipolepi[2];
                    final boolean softi = isSoft[i];
                    final double pdi = ipdamp[i];
                    final double pti = thole[i];
                    final int list[] = lists[i];
                    final int npair = realSpaceCounts[iSymm][i];
                    for (int j = 0; j < npair; j++) {
                        k = list[j];
                        if (!use[k]
                                || // No intermolecular electrostatics
                                (lambdaMode == LambdaMode.VAPOR && intermolecularSoftcore
                                && moleculei != molecule[k])) {
                            continue;
                        }
                        selfScale = 1.0;
                        if (i == k) {
                            selfScale = 0.5;
                        }
                        beta = 0.0;
                        l2 = 1.0;
                        soft = (softi || isSoft[k]);
                        if (soft && doPermanentRealSpace) {
                            beta = lAlpha;
                            l2 = permanentScale;
                        }
                        final double xk = neighborX[k];
                        final double yk = neighborY[k];
                        final double zk = neighborZ[k];
                        dx_local[0] = xk - xi;
                        dx_local[1] = yk - yi;
                        dx_local[2] = zk - zi;
                        r2 = crystal.image(dx_local);
                        xr = dx_local[0];
                        yr = dx_local[1];
                        zr = dx_local[2];
                        final double globalMultipolek[] = neighborMultipole[k];
                        final double inducedDipolek[] = neighborInducedDipole[k];
                        final double inducedDipolepk[] = neighborInducedDipolep[k];
                        ck = globalMultipolek[t000];
                        dkx = globalMultipolek[t100];
                        dky = globalMultipolek[t010];
                        dkz = globalMultipolek[t001];
                        qkxx = globalMultipolek[t200] * oneThird;
                        qkyy = globalMultipolek[t020] * oneThird;
                        qkzz = globalMultipolek[t002] * oneThird;
                        qkxy = globalMultipolek[t110] * oneThird;
                        qkxz = globalMultipolek[t101] * oneThird;
                        qkyz = globalMultipolek[t011] * oneThird;
                        ukx = inducedDipolek[0];
                        uky = inducedDipolek[1];
                        ukz = inducedDipolek[2];
                        pkx = inducedDipolepk[0];
                        pky = inducedDipolepk[1];
                        pkz = inducedDipolepk[2];
                        final double pdk = ipdamp[k];
                        final double ptk = thole[k];
                        scale = masking_local[k];
                        scalep = maskingp_local[k];
                        scaled = maskingd_local[k];
                        scale3 = 1.0;
                        scale5 = 1.0;
                        scale7 = 1.0;
                        double r = sqrt(r2 + beta);
                        double ralpha = aewald * r;
                        double exp2a = exp(-ralpha * ralpha);
                        rr1 = 1.0 / r;
                        rr2 = rr1 * rr1;
                        bn0 = erfc(ralpha) * rr1;
                        bn1 = (bn0 + an0 * exp2a) * rr2;
                        bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                        bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;
                        bn4 = (7.0 * bn3 + an3 * exp2a) * rr2;
                        bn5 = (9.0 * bn4 + an4 * exp2a) * rr2;
                        bn6 = (11.0 * bn5 + an5 * exp2a) * rr2;
                        rr3 = rr1 * rr2;
                        rr5 = 3.0 * rr3 * rr2;
                        rr7 = 5.0 * rr5 * rr2;
                        rr9 = 7.0 * rr7 * rr2;
                        rr11 = 9.0 * rr9 * rr2;
                        rr13 = 11.0 * rr11 * rr2;
                        ddsc3x = 0.0;
                        ddsc3y = 0.0;
                        ddsc3z = 0.0;
                        ddsc5x = 0.0;
                        ddsc5y = 0.0;
                        ddsc5z = 0.0;
                        ddsc7x = 0.0;
                        ddsc7y = 0.0;
                        ddsc7z = 0.0;
                        double damp = pdi * pdk;
                        //if (damp != 0.0) {
                        double pgamma = min(pti, ptk);
                        double rdamp = r * damp;
                        damp = -pgamma * rdamp * rdamp * rdamp;
                        if (damp > -50.0) {
                            final double expdamp = exp(damp);
                            scale3 = 1.0 - expdamp;
                            scale5 = 1.0 - expdamp * (1.0 - damp);
                            scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                            final double temp3 = -3.0 * damp * expdamp * rr2;
                            final double temp5 = -damp;
                            final double temp7 = -0.2 - 0.6 * damp;
                            ddsc3x = temp3 * xr;
                            ddsc3y = temp3 * yr;
                            ddsc3z = temp3 * zr;
                            ddsc5x = temp5 * ddsc3x;
                            ddsc5y = temp5 * ddsc3y;
                            ddsc5z = temp5 * ddsc3z;
                            ddsc7x = temp7 * ddsc5x;
                            ddsc7y = temp7 * ddsc5y;
                            ddsc7z = temp7 * ddsc5z;
                        }
                        //}
                        if (doPermanentRealSpace) {
                            double ei = permanentPair();
                            //log(i,k,r,ei);
                            if (Double.isNaN(ei) || Double.isInfinite(ei)) {
                                logger.info(crystal.getUnitCell().toString());
                                logger.info(atoms[i].toString());
                                logger.info(atoms[k].toString());
                                logger.severe(String.format(" The permanent multipole energy between atoms %d and %d (%d) is %16.8f at %16.8f A.", i, k, iSymm, ei, r));
                            }
                            permanentEnergy += ei;
                            count++;
                        }
                        if (polarization != Polarization.NONE && doPolarization) {
                            /**
                             * Polarization does not use the softcore tensors.
                             */
                            if (soft && doPermanentRealSpace) {
                                scale3 = 1.0;
                                scale5 = 1.0;
                                scale7 = 1.0;
                                r = sqrt(r2);
                                ralpha = aewald * r;
                                exp2a = exp(-ralpha * ralpha);
                                rr1 = 1.0 / r;
                                rr2 = rr1 * rr1;
                                bn0 = erfc(ralpha) * rr1;
                                bn1 = (bn0 + an0 * exp2a) * rr2;
                                bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                                bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;
                                bn4 = (7.0 * bn3 + an3 * exp2a) * rr2;
                                bn5 = (9.0 * bn4 + an4 * exp2a) * rr2;
                                bn6 = (11.0 * bn5 + an5 * exp2a) * rr2;
                                rr3 = rr1 * rr2;
                                rr5 = 3.0 * rr3 * rr2;
                                rr7 = 5.0 * rr5 * rr2;
                                rr9 = 7.0 * rr7 * rr2;
                                rr11 = 9.0 * rr9 * rr2;
                                ddsc3x = 0.0;
                                ddsc3y = 0.0;
                                ddsc3z = 0.0;
                                ddsc5x = 0.0;
                                ddsc5y = 0.0;
                                ddsc5z = 0.0;
                                ddsc7x = 0.0;
                                ddsc7y = 0.0;
                                ddsc7z = 0.0;
                                damp = pdi * pdk;
                                //if (damp != 0.0) {
                                pgamma = min(pti, ptk);
                                rdamp = r * damp;
                                damp = -pgamma * rdamp * rdamp * rdamp;
                                if (damp > -50.0) {
                                    final double expdamp = exp(damp);
                                    scale3 = 1.0 - expdamp;
                                    scale5 = 1.0 - expdamp * (1.0 - damp);
                                    scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                    final double temp3 = -3.0 * damp * expdamp * rr2;
                                    final double temp5 = -damp;
                                    final double temp7 = -0.2 - 0.6 * damp;
                                    ddsc3x = temp3 * xr;
                                    ddsc3y = temp3 * yr;
                                    ddsc3z = temp3 * zr;
                                    ddsc5x = temp5 * ddsc3x;
                                    ddsc5y = temp5 * ddsc3y;
                                    ddsc5z = temp5 * ddsc3z;
                                    ddsc7x = temp7 * ddsc5x;
                                    ddsc7y = temp7 * ddsc5y;
                                    ddsc7z = temp7 * ddsc5z;
                                }
                                //}
                            }
                            double ei = polarizationPair();
                            if (Double.isNaN(ei) || Double.isInfinite(ei)) {
                                logger.info(crystal.getUnitCell().toString());
                                logger.info(atoms[i].toString());
                                logger.info(format(" with induced dipole: %8.3f %8.3f %8.3f", uix, uiy, uiz));
                                logger.info(atoms[k].toString());
                                logger.info(format(" with induced dipole: %8.3f %8.3f %8.3f", ukx, uky, ukz));
                                logger.severe(String.format(" The polarization energy due to atoms %d and %d (%d) is %10.6f at %10.6f A.", i + 1, k + 1, iSymm, ei, r));
                            }
                            inducedEnergy += ei;
                        }
                    }
                    if (iSymm == 0) {
                        for (Atom ak : ai.get1_5s()) {
                            int index = ak.xyzIndex - 1;
                            masking_local[index] = 1.0;
                        }
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = 1.0;
                                for (int j : ip11[i]) {
                                    if (j == index) {
                                        maskingp_local[index] = 1.0;
                                    }
                                }
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = 1.0;
                                maskingp_local[index] = 1.0;
                            }
                        }
                        for (Bond bond : ai.getBonds()) {
                            int index = bond.get1_2(ai).xyzIndex - 1;
                            masking_local[index] = 1.0;
                            maskingp_local[index] = 1.0;
                        }
                        for (int j : ip11[i]) {
                            maskingd_local[j] = 1.0;
                        }
                    }
                }
            }

            /**
             * Evaluate the real space permanent energy for a pair of multipole
             * sites.
             *
             * @return the permanent multipole energy.
             */
            private double permanentPair() {
                final double dixdkx = diy * dkz - diz * dky;
                final double dixdky = diz * dkx - dix * dkz;
                final double dixdkz = dix * dky - diy * dkx;
                final double dixrx = diy * zr - diz * yr;
                final double dixry = diz * xr - dix * zr;
                final double dixrz = dix * yr - diy * xr;
                final double dkxrx = dky * zr - dkz * yr;
                final double dkxry = dkz * xr - dkx * zr;
                final double dkxrz = dkx * yr - dky * xr;
                final double qirx = qixx * xr + qixy * yr + qixz * zr;
                final double qiry = qixy * xr + qiyy * yr + qiyz * zr;
                final double qirz = qixz * xr + qiyz * yr + qizz * zr;
                final double qkrx = qkxx * xr + qkxy * yr + qkxz * zr;
                final double qkry = qkxy * xr + qkyy * yr + qkyz * zr;
                final double qkrz = qkxz * xr + qkyz * yr + qkzz * zr;
                final double qiqkrx = qixx * qkrx + qixy * qkry + qixz * qkrz;
                final double qiqkry = qixy * qkrx + qiyy * qkry + qiyz * qkrz;
                final double qiqkrz = qixz * qkrx + qiyz * qkry + qizz * qkrz;
                final double qkqirx = qkxx * qirx + qkxy * qiry + qkxz * qirz;
                final double qkqiry = qkxy * qirx + qkyy * qiry + qkyz * qirz;
                final double qkqirz = qkxz * qirx + qkyz * qiry + qkzz * qirz;
                final double qixqkx = qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy - qiyz * qkyy - qizz * qkyz;
                final double qixqky = qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz - qixy * qkyz - qixz * qkzz;
                final double qixqkz = qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx - qiyy * qkxy - qiyz * qkxz;
                final double rxqirx = yr * qirz - zr * qiry;
                final double rxqiry = zr * qirx - xr * qirz;
                final double rxqirz = xr * qiry - yr * qirx;
                final double rxqkrx = yr * qkrz - zr * qkry;
                final double rxqkry = zr * qkrx - xr * qkrz;
                final double rxqkrz = xr * qkry - yr * qkrx;
                final double rxqikrx = yr * qiqkrz - zr * qiqkry;
                final double rxqikry = zr * qiqkrx - xr * qiqkrz;
                final double rxqikrz = xr * qiqkry - yr * qiqkrx;
                final double rxqkirx = yr * qkqirz - zr * qkqiry;
                final double rxqkiry = zr * qkqirx - xr * qkqirz;
                final double rxqkirz = xr * qkqiry - yr * qkqirx;
                final double qkrxqirx = qkry * qirz - qkrz * qiry;
                final double qkrxqiry = qkrz * qirx - qkrx * qirz;
                final double qkrxqirz = qkrx * qiry - qkry * qirx;
                final double qidkx = qixx * dkx + qixy * dky + qixz * dkz;
                final double qidky = qixy * dkx + qiyy * dky + qiyz * dkz;
                final double qidkz = qixz * dkx + qiyz * dky + qizz * dkz;
                final double qkdix = qkxx * dix + qkxy * diy + qkxz * diz;
                final double qkdiy = qkxy * dix + qkyy * diy + qkyz * diz;
                final double qkdiz = qkxz * dix + qkyz * diy + qkzz * diz;
                final double dixqkrx = diy * qkrz - diz * qkry;
                final double dixqkry = diz * qkrx - dix * qkrz;
                final double dixqkrz = dix * qkry - diy * qkrx;
                final double dkxqirx = dky * qirz - dkz * qiry;
                final double dkxqiry = dkz * qirx - dkx * qirz;
                final double dkxqirz = dkx * qiry - dky * qirx;
                final double rxqidkx = yr * qidkz - zr * qidky;
                final double rxqidky = zr * qidkx - xr * qidkz;
                final double rxqidkz = xr * qidky - yr * qidkx;
                final double rxqkdix = yr * qkdiz - zr * qkdiy;
                final double rxqkdiy = zr * qkdix - xr * qkdiz;
                final double rxqkdiz = xr * qkdiy - yr * qkdix;
                /**
                 * Calculate the scalar products for permanent multipoles.
                 */
                final double sc2 = dix * dkx + diy * dky + diz * dkz;
                final double sc3 = dix * xr + diy * yr + diz * zr;
                final double sc4 = dkx * xr + dky * yr + dkz * zr;
                final double sc5 = qirx * xr + qiry * yr + qirz * zr;
                final double sc6 = qkrx * xr + qkry * yr + qkrz * zr;
                final double sc7 = qirx * dkx + qiry * dky + qirz * dkz;
                final double sc8 = qkrx * dix + qkry * diy + qkrz * diz;
                final double sc9 = qirx * qkrx + qiry * qkry + qirz * qkrz;
                final double sc10 = 2.0 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) + qixx * qkxx + qiyy * qkyy + qizz * qkzz;
                /**
                 * Calculate the gl functions for permanent multipoles.
                 */
                final double gl0 = ci * ck;
                final double gl1 = ck * sc3 - ci * sc4;
                final double gl2 = ci * sc6 + ck * sc5 - sc3 * sc4;
                final double gl3 = sc3 * sc6 - sc4 * sc5;
                final double gl4 = sc5 * sc6;
                final double gl5 = -4.0 * sc9;
                final double gl6 = sc2;
                final double gl7 = 2.0 * (sc7 - sc8);
                final double gl8 = 2.0 * sc10;
                /**
                 * Compute the energy contributions for this interaction.
                 */
                final double scale1 = 1.0 - scale;
                final double ereal = gl0 * bn0 + (gl1 + gl6) * bn1 + (gl2 + gl7 + gl8) * bn2 + (gl3 + gl5) * bn3 + gl4 * bn4;
                final double efix = scale1 * (gl0 * rr1 + (gl1 + gl6) * rr3 + (gl2 + gl7 + gl8) * rr5 + (gl3 + gl5) * rr7 + gl4 * rr9);
                final double e = selfScale * l2 * (ereal - efix);
                if (gradient) {
                    final double gf1 = bn1 * gl0 + bn2 * (gl1 + gl6) + bn3 * (gl2 + gl7 + gl8) + bn4 * (gl3 + gl5) + bn5 * gl4;
                    final double gf2 = -ck * bn1 + sc4 * bn2 - sc6 * bn3;
                    final double gf3 = ci * bn1 + sc3 * bn2 + sc5 * bn3;
                    final double gf4 = 2.0 * bn2;
                    final double gf5 = 2.0 * (-ck * bn2 + sc4 * bn3 - sc6 * bn4);
                    final double gf6 = 2.0 * (-ci * bn2 - sc3 * bn3 - sc5 * bn4);
                    final double gf7 = 4.0 * bn3;
                    /*
                     * Get the permanent force with screening.
                     */
                    double ftm2x = gf1 * xr + gf2 * dix + gf3 * dkx + gf4 * (qkdix - qidkx) + gf5 * qirx + gf6 * qkrx + gf7 * (qiqkrx + qkqirx);
                    double ftm2y = gf1 * yr + gf2 * diy + gf3 * dky + gf4 * (qkdiy - qidky) + gf5 * qiry + gf6 * qkry + gf7 * (qiqkry + qkqiry);
                    double ftm2z = gf1 * zr + gf2 * diz + gf3 * dkz + gf4 * (qkdiz - qidkz) + gf5 * qirz + gf6 * qkrz + gf7 * (qiqkrz + qkqirz);
                    /*
                     * Get the permanent torque with screening.
                     */
                    double ttm2x = -bn1 * dixdkx + gf2 * dixrx + gf4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx) - gf5 * rxqirx - gf7 * (rxqikrx + qkrxqirx);
                    double ttm2y = -bn1 * dixdky + gf2 * dixry + gf4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky) - gf5 * rxqiry - gf7 * (rxqikry + qkrxqiry);
                    double ttm2z = -bn1 * dixdkz + gf2 * dixrz + gf4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz) - gf5 * rxqirz - gf7 * (rxqikrz + qkrxqirz);
                    double ttm3x = bn1 * dixdkx + gf3 * dkxrx - gf4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx) - gf6 * rxqkrx - gf7 * (rxqkirx - qkrxqirx);
                    double ttm3y = bn1 * dixdky + gf3 * dkxry - gf4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky) - gf6 * rxqkry - gf7 * (rxqkiry - qkrxqiry);
                    double ttm3z = bn1 * dixdkz + gf3 * dkxrz - gf4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz) - gf6 * rxqkrz - gf7 * (rxqkirz - qkrxqirz);
                    /**
                     * Handle the case where scaling is used.
                     */
                    if (scale1 != 0.0) {
                        final double gfr1 = rr3 * gl0 + rr5 * (gl1 + gl6) + rr7 * (gl2 + gl7 + gl8) + rr9 * (gl3 + gl5) + rr11 * gl4;
                        final double gfr2 = -ck * rr3 + sc4 * rr5 - sc6 * rr7;
                        final double gfr3 = ci * rr3 + sc3 * rr5 + sc5 * rr7;
                        final double gfr4 = 2.0 * rr5;
                        final double gfr5 = 2.0 * (-ck * rr5 + sc4 * rr7 - sc6 * rr9);
                        final double gfr6 = 2.0 * (-ci * rr5 - sc3 * rr7 - sc5 * rr9);
                        final double gfr7 = 4.0 * rr7;
                        /*
                         * Get the permanent force without screening.
                         */
                        final double ftm2rx = gfr1 * xr + gfr2 * dix + gfr3 * dkx + gfr4 * (qkdix - qidkx) + gfr5 * qirx + gfr6 * qkrx + gfr7 * (qiqkrx + qkqirx);
                        final double ftm2ry = gfr1 * yr + gfr2 * diy + gfr3 * dky + gfr4 * (qkdiy - qidky) + gfr5 * qiry + gfr6 * qkry + gfr7 * (qiqkry + qkqiry);
                        final double ftm2rz = gfr1 * zr + gfr2 * diz + gfr3 * dkz + gfr4 * (qkdiz - qidkz) + gfr5 * qirz + gfr6 * qkrz + gfr7 * (qiqkrz + qkqirz);
                        /*
                         * Get the permanent torque without screening.
                         */
                        final double ttm2rx = -rr3 * dixdkx + gfr2 * dixrx + gfr4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx) - gfr5 * rxqirx - gfr7 * (rxqikrx + qkrxqirx);
                        final double ttm2ry = -rr3 * dixdky + gfr2 * dixry + gfr4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky) - gfr5 * rxqiry - gfr7 * (rxqikry + qkrxqiry);
                        final double ttm2rz = -rr3 * dixdkz + gfr2 * dixrz + gfr4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz) - gfr5 * rxqirz - gfr7 * (rxqikrz + qkrxqirz);
                        final double ttm3rx = rr3 * dixdkx + gfr3 * dkxrx - gfr4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx) - gfr6 * rxqkrx - gfr7 * (rxqkirx - qkrxqirx);
                        final double ttm3ry = rr3 * dixdky + gfr3 * dkxry - gfr4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky) - gfr6 * rxqkry - gfr7 * (rxqkiry - qkrxqiry);
                        final double ttm3rz = rr3 * dixdkz + gfr3 * dkxrz - gfr4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz) - gfr6 * rxqkrz - gfr7 * (rxqkirz - qkrxqirz);
                        ftm2x -= scale1 * ftm2rx;
                        ftm2y -= scale1 * ftm2ry;
                        ftm2z -= scale1 * ftm2rz;
                        ttm2x -= scale1 * ttm2rx;
                        ttm2y -= scale1 * ttm2ry;
                        ttm2z -= scale1 * ttm2rz;
                        ttm3x -= scale1 * ttm3rx;
                        ttm3y -= scale1 * ttm3ry;
                        ttm3z -= scale1 * ttm3rz;
                    }
                    double prefactor = ELECTRIC * selfScale * l2;
                    gX[i] += prefactor * ftm2x;
                    gY[i] += prefactor * ftm2y;
                    gZ[i] += prefactor * ftm2z;
                    tX[i] += prefactor * ttm2x;
                    tY[i] += prefactor * ttm2y;
                    tZ[i] += prefactor * ttm2z;
                    gxk_local[k] -= prefactor * ftm2x;
                    gyk_local[k] -= prefactor * ftm2y;
                    gzk_local[k] -= prefactor * ftm2z;
                    txk_local[k] += prefactor * ttm3x;
                    tyk_local[k] += prefactor * ttm3y;
                    tzk_local[k] += prefactor * ttm3z;
                    /**
                     * This is dU/dL/dX for the first term of dU/dL: d[dlPow *
                     * ereal]/dx
                     */
                    if (lambdaTerm && soft) {
                        prefactor = ELECTRIC * selfScale * dEdLSign * dlPowPerm;
                        lgX[i] += prefactor * ftm2x;
                        lgY[i] += prefactor * ftm2y;
                        lgZ[i] += prefactor * ftm2z;
                        ltX[i] += prefactor * ttm2x;
                        ltY[i] += prefactor * ttm2y;
                        ltZ[i] += prefactor * ttm2z;
                        lxk_local[k] -= prefactor * ftm2x;
                        lyk_local[k] -= prefactor * ftm2y;
                        lzk_local[k] -= prefactor * ftm2z;
                        ltxk_local[k] += prefactor * ttm3x;
                        ltyk_local[k] += prefactor * ttm3y;
                        ltzk_local[k] += prefactor * ttm3z;
                    }
                }
                if (lambdaTerm && soft) {
                    double dRealdL = gl0 * bn1 + (gl1 + gl6) * bn2 + (gl2 + gl7 + gl8) * bn3 + (gl3 + gl5) * bn4 + gl4 * bn5;
                    double d2RealdL2 = gl0 * bn2 + (gl1 + gl6) * bn3 + (gl2 + gl7 + gl8) * bn4 + (gl3 + gl5) * bn5 + gl4 * bn6;

                    dUdL += selfScale * (dEdLSign * dlPowPerm * ereal + l2 * dlAlpha * dRealdL);
                    d2UdL2 += selfScale * (dEdLSign * (d2lPowPerm * ereal
                            + dlPowPerm * dlAlpha * dRealdL
                            + dlPowPerm * dlAlpha * dRealdL)
                            + l2 * d2lAlpha * dRealdL
                            + l2 * dlAlpha * dlAlpha * d2RealdL2);

                    double dFixdL = gl0 * rr3 + (gl1 + gl6) * rr5 + (gl2 + gl7 + gl8) * rr7 + (gl3 + gl5) * rr9 + gl4 * rr11;
                    double d2FixdL2 = gl0 * rr5 + (gl1 + gl6) * rr7 + (gl2 + gl7 + gl8) * rr9 + (gl3 + gl5) * rr11 + gl4 * rr13;
                    dFixdL *= scale1;
                    d2FixdL2 *= scale1;
                    dUdL -= selfScale * (dEdLSign * dlPowPerm * efix + l2 * dlAlpha * dFixdL);
                    d2UdL2 -= selfScale * (dEdLSign * (d2lPowPerm * efix
                            + dlPowPerm * dlAlpha * dFixdL
                            + dlPowPerm * dlAlpha * dFixdL)
                            + l2 * d2lAlpha * dFixdL
                            + l2 * dlAlpha * dlAlpha * d2FixdL2);
                    /**
                     * Collect terms for dU/dL/dX for the second term of dU/dL:
                     * d[fL2*dfL1dL*dRealdL]/dX
                     */
                    final double gf1 = bn2 * gl0 + bn3 * (gl1 + gl6)
                            + bn4 * (gl2 + gl7 + gl8)
                            + bn5 * (gl3 + gl5) + bn6 * gl4;
                    final double gf2 = -ck * bn2 + sc4 * bn3 - sc6 * bn4;
                    final double gf3 = ci * bn2 + sc3 * bn3 + sc5 * bn4;
                    final double gf4 = 2.0 * bn3;
                    final double gf5 = 2.0 * (-ck * bn3 + sc4 * bn4 - sc6 * bn5);
                    final double gf6 = 2.0 * (-ci * bn3 - sc3 * bn4 - sc5 * bn5);
                    final double gf7 = 4.0 * bn4;
                    /*
                     * Get the permanent force with screening.
                     */
                    double ftm2x = gf1 * xr + gf2 * dix + gf3 * dkx
                            + gf4 * (qkdix - qidkx) + gf5 * qirx
                            + gf6 * qkrx + gf7 * (qiqkrx + qkqirx);
                    double ftm2y = gf1 * yr + gf2 * diy + gf3 * dky
                            + gf4 * (qkdiy - qidky) + gf5 * qiry
                            + gf6 * qkry + gf7 * (qiqkry + qkqiry);
                    double ftm2z = gf1 * zr + gf2 * diz + gf3 * dkz
                            + gf4 * (qkdiz - qidkz) + gf5 * qirz
                            + gf6 * qkrz + gf7 * (qiqkrz + qkqirz);
                    /*
                     * Get the permanent torque with screening.
                     */
                    double ttm2x = -bn2 * dixdkx + gf2 * dixrx
                            + gf4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx)
                            - gf5 * rxqirx - gf7 * (rxqikrx + qkrxqirx);
                    double ttm2y = -bn2 * dixdky + gf2 * dixry
                            + gf4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky)
                            - gf5 * rxqiry - gf7 * (rxqikry + qkrxqiry);
                    double ttm2z = -bn2 * dixdkz + gf2 * dixrz
                            + gf4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz)
                            - gf5 * rxqirz - gf7 * (rxqikrz + qkrxqirz);
                    double ttm3x = bn2 * dixdkx + gf3 * dkxrx
                            - gf4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx)
                            - gf6 * rxqkrx - gf7 * (rxqkirx - qkrxqirx);
                    double ttm3y = bn2 * dixdky + gf3 * dkxry
                            - gf4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky)
                            - gf6 * rxqkry - gf7 * (rxqkiry - qkrxqiry);
                    double ttm3z = bn2 * dixdkz + gf3 * dkxrz
                            - gf4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz)
                            - gf6 * rxqkrz - gf7 * (rxqkirz - qkrxqirz);

                    /**
                     * Handle the case where scaling is used.
                     */
                    if (scale1 != 0.0) {
                        final double gfr1 = rr5 * gl0 + rr7 * (gl1 + gl6) + rr9 * (gl2 + gl7 + gl8) + rr11 * (gl3 + gl5) + rr13 * gl4;
                        final double gfr2 = -ck * rr5 + sc4 * rr7 - sc6 * rr9;
                        final double gfr3 = ci * rr5 + sc3 * rr7 + sc5 * rr9;
                        final double gfr4 = 2.0 * rr7;
                        final double gfr5 = 2.0 * (-ck * rr7 + sc4 * rr9 - sc6 * rr11);
                        final double gfr6 = 2.0 * (-ci * rr7 - sc3 * rr9 - sc5 * rr11);
                        final double gfr7 = 4.0 * rr9;

                        //Get the permanent force without screening.
                        final double ftm2rx = gfr1 * xr + gfr2 * dix + gfr3 * dkx + gfr4 * (qkdix - qidkx) + gfr5 * qirx + gfr6 * qkrx + gfr7 * (qiqkrx + qkqirx);
                        final double ftm2ry = gfr1 * yr + gfr2 * diy + gfr3 * dky + gfr4 * (qkdiy - qidky) + gfr5 * qiry + gfr6 * qkry + gfr7 * (qiqkry + qkqiry);
                        final double ftm2rz = gfr1 * zr + gfr2 * diz + gfr3 * dkz + gfr4 * (qkdiz - qidkz) + gfr5 * qirz + gfr6 * qkrz + gfr7 * (qiqkrz + qkqirz);

                        // Get the permanent torque without screening.
                        final double ttm2rx = -rr5 * dixdkx + gfr2 * dixrx + gfr4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx) - gfr5 * rxqirx - gfr7 * (rxqikrx + qkrxqirx);
                        final double ttm2ry = -rr5 * dixdky + gfr2 * dixry + gfr4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky) - gfr5 * rxqiry - gfr7 * (rxqikry + qkrxqiry);
                        final double ttm2rz = -rr5 * dixdkz + gfr2 * dixrz + gfr4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz) - gfr5 * rxqirz - gfr7 * (rxqikrz + qkrxqirz);
                        final double ttm3rx = rr5 * dixdkx + gfr3 * dkxrx - gfr4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx) - gfr6 * rxqkrx - gfr7 * (rxqkirx - qkrxqirx);
                        final double ttm3ry = rr5 * dixdky + gfr3 * dkxry - gfr4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky) - gfr6 * rxqkry - gfr7 * (rxqkiry - qkrxqiry);
                        final double ttm3rz = rr5 * dixdkz + gfr3 * dkxrz - gfr4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz) - gfr6 * rxqkrz - gfr7 * (rxqkirz - qkrxqirz);
                        ftm2x -= scale1 * ftm2rx;
                        ftm2y -= scale1 * ftm2ry;
                        ftm2z -= scale1 * ftm2rz;
                        ttm2x -= scale1 * ttm2rx;
                        ttm2y -= scale1 * ttm2ry;
                        ttm2z -= scale1 * ttm2rz;
                        ttm3x -= scale1 * ttm3rx;
                        ttm3y -= scale1 * ttm3ry;
                        ttm3z -= scale1 * ttm3rz;
                    }
                    /**
                     * Add in dU/dL/dX for the second term of dU/dL:
                     * d[lPow*dlAlpha*dRealdL]/dX
                     */
                    double prefactor = ELECTRIC * selfScale * l2 * dlAlpha;
                    lgX[i] += prefactor * ftm2x;
                    lgY[i] += prefactor * ftm2y;
                    lgZ[i] += prefactor * ftm2z;
                    ltX[i] += prefactor * ttm2x;
                    ltY[i] += prefactor * ttm2y;
                    ltZ[i] += prefactor * ttm2z;
                    lxk_local[k] -= prefactor * ftm2x;
                    lyk_local[k] -= prefactor * ftm2y;
                    lzk_local[k] -= prefactor * ftm2z;
                    ltxk_local[k] += prefactor * ttm3x;
                    ltyk_local[k] += prefactor * ttm3y;
                    ltzk_local[k] += prefactor * ttm3z;
                }
                return e;
            }

            /**
             * Evaluate the polarization energy for a pair of polarizable
             * multipole sites.
             *
             * @return the polarization energy.
             */
            private double polarizationPair() {
                final double dsc3 = 1.0 - scale3 * scaled;
                final double dsc5 = 1.0 - scale5 * scaled;
                final double dsc7 = 1.0 - scale7 * scaled;
                final double psc3 = 1.0 - scale3 * scalep;
                final double psc5 = 1.0 - scale5 * scalep;
                final double psc7 = 1.0 - scale7 * scalep;
                final double usc3 = 1.0 - scale3;
                final double usc5 = 1.0 - scale5;
                final double dixukx = diy * ukz - diz * uky;
                final double dixuky = diz * ukx - dix * ukz;
                final double dixukz = dix * uky - diy * ukx;
                final double dkxuix = dky * uiz - dkz * uiy;
                final double dkxuiy = dkz * uix - dkx * uiz;
                final double dkxuiz = dkx * uiy - dky * uix;
                final double dixukpx = diy * pkz - diz * pky;
                final double dixukpy = diz * pkx - dix * pkz;
                final double dixukpz = dix * pky - diy * pkx;
                final double dkxuipx = dky * piz - dkz * piy;
                final double dkxuipy = dkz * pix - dkx * piz;
                final double dkxuipz = dkx * piy - dky * pix;
                final double dixrx = diy * zr - diz * yr;
                final double dixry = diz * xr - dix * zr;
                final double dixrz = dix * yr - diy * xr;
                final double dkxrx = dky * zr - dkz * yr;
                final double dkxry = dkz * xr - dkx * zr;
                final double dkxrz = dkx * yr - dky * xr;
                final double qirx = qixx * xr + qixy * yr + qixz * zr;
                final double qiry = qixy * xr + qiyy * yr + qiyz * zr;
                final double qirz = qixz * xr + qiyz * yr + qizz * zr;
                final double qkrx = qkxx * xr + qkxy * yr + qkxz * zr;
                final double qkry = qkxy * xr + qkyy * yr + qkyz * zr;
                final double qkrz = qkxz * xr + qkyz * yr + qkzz * zr;
                final double rxqirx = yr * qirz - zr * qiry;
                final double rxqiry = zr * qirx - xr * qirz;
                final double rxqirz = xr * qiry - yr * qirx;
                final double rxqkrx = yr * qkrz - zr * qkry;
                final double rxqkry = zr * qkrx - xr * qkrz;
                final double rxqkrz = xr * qkry - yr * qkrx;
                final double qiukx = qixx * ukx + qixy * uky + qixz * ukz;
                final double qiuky = qixy * ukx + qiyy * uky + qiyz * ukz;
                final double qiukz = qixz * ukx + qiyz * uky + qizz * ukz;
                final double qkuix = qkxx * uix + qkxy * uiy + qkxz * uiz;
                final double qkuiy = qkxy * uix + qkyy * uiy + qkyz * uiz;
                final double qkuiz = qkxz * uix + qkyz * uiy + qkzz * uiz;
                final double qiukpx = qixx * pkx + qixy * pky + qixz * pkz;
                final double qiukpy = qixy * pkx + qiyy * pky + qiyz * pkz;
                final double qiukpz = qixz * pkx + qiyz * pky + qizz * pkz;
                final double qkuipx = qkxx * pix + qkxy * piy + qkxz * piz;
                final double qkuipy = qkxy * pix + qkyy * piy + qkyz * piz;
                final double qkuipz = qkxz * pix + qkyz * piy + qkzz * piz;
                final double uixqkrx = uiy * qkrz - uiz * qkry;
                final double uixqkry = uiz * qkrx - uix * qkrz;
                final double uixqkrz = uix * qkry - uiy * qkrx;
                final double ukxqirx = uky * qirz - ukz * qiry;
                final double ukxqiry = ukz * qirx - ukx * qirz;
                final double ukxqirz = ukx * qiry - uky * qirx;
                final double uixqkrpx = piy * qkrz - piz * qkry;
                final double uixqkrpy = piz * qkrx - pix * qkrz;
                final double uixqkrpz = pix * qkry - piy * qkrx;
                final double ukxqirpx = pky * qirz - pkz * qiry;
                final double ukxqirpy = pkz * qirx - pkx * qirz;
                final double ukxqirpz = pkx * qiry - pky * qirx;
                final double rxqiukx = yr * qiukz - zr * qiuky;
                final double rxqiuky = zr * qiukx - xr * qiukz;
                final double rxqiukz = xr * qiuky - yr * qiukx;
                final double rxqkuix = yr * qkuiz - zr * qkuiy;
                final double rxqkuiy = zr * qkuix - xr * qkuiz;
                final double rxqkuiz = xr * qkuiy - yr * qkuix;
                final double rxqiukpx = yr * qiukpz - zr * qiukpy;
                final double rxqiukpy = zr * qiukpx - xr * qiukpz;
                final double rxqiukpz = xr * qiukpy - yr * qiukpx;
                final double rxqkuipx = yr * qkuipz - zr * qkuipy;
                final double rxqkuipy = zr * qkuipx - xr * qkuipz;
                final double rxqkuipz = xr * qkuipy - yr * qkuipx;
                /**
                 * Calculate the scalar products for permanent multipoles.
                 */
                final double sc3 = dix * xr + diy * yr + diz * zr;
                final double sc4 = dkx * xr + dky * yr + dkz * zr;
                final double sc5 = qirx * xr + qiry * yr + qirz * zr;
                final double sc6 = qkrx * xr + qkry * yr + qkrz * zr;
                /**
                 * Calculate the scalar products for polarization components.
                 */
                final double sci1 = uix * dkx + uiy * dky + uiz * dkz + dix * ukx + diy * uky + diz * ukz;
                final double sci3 = uix * xr + uiy * yr + uiz * zr;
                final double sci4 = ukx * xr + uky * yr + ukz * zr;
                final double sci7 = qirx * ukx + qiry * uky + qirz * ukz;
                final double sci8 = qkrx * uix + qkry * uiy + qkrz * uiz;
                final double scip1 = pix * dkx + piy * dky + piz * dkz + dix * pkx + diy * pky + diz * pkz;
                final double scip2 = uix * pkx + uiy * pky + uiz * pkz + pix * ukx + piy * uky + piz * ukz;
                final double scip3 = pix * xr + piy * yr + piz * zr;
                final double scip4 = pkx * xr + pky * yr + pkz * zr;
                final double scip7 = qirx * pkx + qiry * pky + qirz * pkz;
                final double scip8 = qkrx * pix + qkry * piy + qkrz * piz;
                /**
                 * Calculate the gl functions for polarization components.
                 */
                final double gli1 = ck * sci3 - ci * sci4;
                final double gli2 = -sc3 * sci4 - sci3 * sc4;
                final double gli3 = sci3 * sc6 - sci4 * sc5;
                final double gli6 = sci1;
                final double gli7 = 2.0 * (sci7 - sci8);
                final double glip1 = ck * scip3 - ci * scip4;
                final double glip2 = -sc3 * scip4 - scip3 * sc4;
                final double glip3 = scip3 * sc6 - scip4 * sc5;
                final double glip6 = scip1;
                final double glip7 = 2.0 * (scip7 - scip8);
                /**
                 * Compute the energy contributions for this interaction.
                 */
                final double ereal = (gli1 + gli6) * bn1 + (gli2 + gli7) * bn2 + gli3 * bn3;
                final double efix = (gli1 + gli6) * rr3 * psc3 + (gli2 + gli7) * rr5 * psc5 + gli3 * rr7 * psc7;
                final double e = selfScale * 0.5 * (ereal - efix);
                if (!(gradient || lambdaTerm)) {
                    return polarizationScale * e;
                }
                boolean dorli = false;
                if (psc3 != 0.0 || dsc3 != 0.0 || usc3 != 0.0) {
                    dorli = true;
                }
                /*
                 * Get the induced force with screening.
                 */
                final double gfi1 = 0.5 * bn2 * (gli1 + glip1 + gli6 + glip6) + 0.5 * bn2 * scip2 + 0.5 * bn3 * (gli2 + glip2 + gli7 + glip7) - 0.5 * bn3 * (sci3 * scip4 + scip3 * sci4) + 0.5 * bn4 * (gli3 + glip3);
                final double gfi2 = -ck * bn1 + sc4 * bn2 - sc6 * bn3;
                final double gfi3 = ci * bn1 + sc3 * bn2 + sc5 * bn3;
                final double gfi4 = 2.0 * bn2;
                final double gfi5 = bn3 * (sci4 + scip4);
                final double gfi6 = -bn3 * (sci3 + scip3);
                double ftm2ix = gfi1 * xr + 0.5 * (gfi2 * (uix + pix) + bn2 * (sci4 * pix + scip4 * uix) + gfi3 * (ukx + pkx) + bn2 * (sci3 * pkx + scip3 * ukx) + (sci4 + scip4) * bn2 * dix + (sci3 + scip3) * bn2 * dkx + gfi4 * (qkuix + qkuipx - qiukx - qiukpx)) + gfi5 * qirx + gfi6 * qkrx;
                double ftm2iy = gfi1 * yr + 0.5 * (gfi2 * (uiy + piy) + bn2 * (sci4 * piy + scip4 * uiy) + gfi3 * (uky + pky) + bn2 * (sci3 * pky + scip3 * uky) + (sci4 + scip4) * bn2 * diy + (sci3 + scip3) * bn2 * dky + gfi4 * (qkuiy + qkuipy - qiuky - qiukpy)) + gfi5 * qiry + gfi6 * qkry;
                double ftm2iz = gfi1 * zr + 0.5 * (gfi2 * (uiz + piz) + bn2 * (sci4 * piz + scip4 * uiz) + gfi3 * (ukz + pkz) + bn2 * (sci3 * pkz + scip3 * ukz) + (sci4 + scip4) * bn2 * diz + (sci3 + scip3) * bn2 * dkz + gfi4 * (qkuiz + qkuipz - qiukz - qiukpz)) + gfi5 * qirz + gfi6 * qkrz;
                /*
                 * Get the induced torque with screening.
                 */
                final double gti2 = 0.5 * bn2 * (sci4 + scip4);
                final double gti3 = 0.5 * bn2 * (sci3 + scip3);
                final double gti4 = gfi4;
                final double gti5 = gfi5;
                final double gti6 = gfi6;
                double ttm2ix = -0.5 * bn1 * (dixukx + dixukpx) + gti2 * dixrx - gti5 * rxqirx + 0.5 * gti4 * (ukxqirx + rxqiukx + ukxqirpx + rxqiukpx);
                double ttm2iy = -0.5 * bn1 * (dixuky + dixukpy) + gti2 * dixry - gti5 * rxqiry + 0.5 * gti4 * (ukxqiry + rxqiuky + ukxqirpy + rxqiukpy);
                double ttm2iz = -0.5 * bn1 * (dixukz + dixukpz) + gti2 * dixrz - gti5 * rxqirz + 0.5 * gti4 * (ukxqirz + rxqiukz + ukxqirpz + rxqiukpz);
                double ttm3ix = -0.5 * bn1 * (dkxuix + dkxuipx) + gti3 * dkxrx - gti6 * rxqkrx - 0.5 * gti4 * (uixqkrx + rxqkuix + uixqkrpx + rxqkuipx);
                double ttm3iy = -0.5 * bn1 * (dkxuiy + dkxuipy) + gti3 * dkxry - gti6 * rxqkry - 0.5 * gti4 * (uixqkry + rxqkuiy + uixqkrpy + rxqkuipy);
                double ttm3iz = -0.5 * bn1 * (dkxuiz + dkxuipz) + gti3 * dkxrz - gti6 * rxqkrz - 0.5 * gti4 * (uixqkrz + rxqkuiz + uixqkrpz + rxqkuipz);
                double ftm2rix = 0.0;
                double ftm2riy = 0.0;
                double ftm2riz = 0.0;
                double ttm2rix = 0.0;
                double ttm2riy = 0.0;
                double ttm2riz = 0.0;
                double ttm3rix = 0.0;
                double ttm3riy = 0.0;
                double ttm3riz = 0.0;
                if (dorli) {
                    /*
                     * Get the induced force without screening.
                     */
                    final double gfri1 = 0.5 * rr5 * ((gli1 + gli6) * psc3 + (glip1 + glip6) * dsc3 + scip2 * usc3) + 0.5 * rr7 * ((gli7 + gli2) * psc5 + (glip7 + glip2) * dsc5 - (sci3 * scip4 + scip3 * sci4) * usc5) + 0.5 * rr9 * (gli3 * psc7 + glip3 * dsc7);
                    final double gfri4 = 2.0 * rr5;
                    final double gfri5 = rr7 * (sci4 * psc7 + scip4 * dsc7);
                    final double gfri6 = -rr7 * (sci3 * psc7 + scip3 * dsc7);
                    ftm2rix = gfri1 * xr + 0.5 * (-rr3 * ck * (uix * psc3 + pix * dsc3) + rr5 * sc4 * (uix * psc5 + pix * dsc5) - rr7 * sc6 * (uix * psc7 + pix * dsc7)) + (rr3 * ci * (ukx * psc3 + pkx * dsc3) + rr5 * sc3 * (ukx * psc5 + pkx * dsc5) + rr7 * sc5 * (ukx * psc7 + pkx * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * dix + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dkx + 0.5 * gfri4 * ((qkuix - qiukx) * psc5 + (qkuipx - qiukpx) * dsc5) + gfri5 * qirx + gfri6 * qkrx;
                    ftm2riy = gfri1 * yr + 0.5 * (-rr3 * ck * (uiy * psc3 + piy * dsc3) + rr5 * sc4 * (uiy * psc5 + piy * dsc5) - rr7 * sc6 * (uiy * psc7 + piy * dsc7)) + (rr3 * ci * (uky * psc3 + pky * dsc3) + rr5 * sc3 * (uky * psc5 + pky * dsc5) + rr7 * sc5 * (uky * psc7 + pky * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * diy + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dky + 0.5 * gfri4 * ((qkuiy - qiuky) * psc5 + (qkuipy - qiukpy) * dsc5) + gfri5 * qiry + gfri6 * qkry;
                    ftm2riz = gfri1 * zr + 0.5 * (-rr3 * ck * (uiz * psc3 + piz * dsc3) + rr5 * sc4 * (uiz * psc5 + piz * dsc5) - rr7 * sc6 * (uiz * psc7 + piz * dsc7)) + (rr3 * ci * (ukz * psc3 + pkz * dsc3) + rr5 * sc3 * (ukz * psc5 + pkz * dsc5) + rr7 * sc5 * (ukz * psc7 + pkz * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * diz + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dkz + 0.5 * gfri4 * ((qkuiz - qiukz) * psc5 + (qkuipz - qiukpz) * dsc5) + gfri5 * qirz + gfri6 * qkrz;
                    /*
                     * Get the induced torque without screening.
                     */
                    final double gtri2 = 0.5 * rr5 * (sci4 * psc5 + scip4 * dsc5);
                    final double gtri3 = 0.5 * rr5 * (sci3 * psc5 + scip3 * dsc5);
                    final double gtri4 = gfri4;
                    final double gtri5 = gfri5;
                    final double gtri6 = gfri6;
                    ttm2rix = -rr3 * (dixukx * psc3 + dixukpx * dsc3) * 0.5 + gtri2 * dixrx - gtri5 * rxqirx + gtri4 * ((ukxqirx + rxqiukx) * psc5 + (ukxqirpx + rxqiukpx) * dsc5) * 0.5;
                    ttm2riy = -rr3 * (dixuky * psc3 + dixukpy * dsc3) * 0.5 + gtri2 * dixry - gtri5 * rxqiry + gtri4 * ((ukxqiry + rxqiuky) * psc5 + (ukxqirpy + rxqiukpy) * dsc5) * 0.5;
                    ttm2riz = -rr3 * (dixukz * psc3 + dixukpz * dsc3) * 0.5 + gtri2 * dixrz - gtri5 * rxqirz + gtri4 * ((ukxqirz + rxqiukz) * psc5 + (ukxqirpz + rxqiukpz) * dsc5) * 0.5;
                    ttm3rix = -rr3 * (dkxuix * psc3 + dkxuipx * dsc3) * 0.5 + gtri3 * dkxrx - gtri6 * rxqkrx - gtri4 * ((uixqkrx + rxqkuix) * psc5 + (uixqkrpx + rxqkuipx) * dsc5) * 0.5;
                    ttm3riy = -rr3 * (dkxuiy * psc3 + dkxuipy * dsc3) * 0.5 + gtri3 * dkxry - gtri6 * rxqkry - gtri4 * ((uixqkry + rxqkuiy) * psc5 + (uixqkrpy + rxqkuipy) * dsc5) * 0.5;
                    ttm3riz = -rr3 * (dkxuiz * psc3 + dkxuipz * dsc3) * 0.5 + gtri3 * dkxrz - gtri6 * rxqkrz - gtri4 * ((uixqkrz + rxqkuiz) * psc5 + (uixqkrpz + rxqkuipz) * dsc5) * 0.5;
                }
                /*
                 * Account for partially excluded induced interactions.
                 */
                double temp3 = 0.5 * rr3 * ((gli1 + gli6) * scalep + (glip1 + glip6) * scaled);
                double temp5 = 0.5 * rr5 * ((gli2 + gli7) * scalep + (glip2 + glip7) * scaled);
                final double temp7 = 0.5 * rr7 * (gli3 * scalep + glip3 * scaled);
                final double fridmpx = temp3 * ddsc3x + temp5 * ddsc5x + temp7 * ddsc7x;
                final double fridmpy = temp3 * ddsc3y + temp5 * ddsc5y + temp7 * ddsc7y;
                final double fridmpz = temp3 * ddsc3z + temp5 * ddsc5z + temp7 * ddsc7z;
                /*
                 * Find some scaling terms for induced-induced force.
                 */
                temp3 = 0.5 * rr3 * scip2;
                temp5 = -0.5 * rr5 * (sci3 * scip4 + scip3 * sci4);
                final double findmpx = temp3 * ddsc3x + temp5 * ddsc5x;
                final double findmpy = temp3 * ddsc3y + temp5 * ddsc5y;
                final double findmpz = temp3 * ddsc3z + temp5 * ddsc5z;
                /*
                 * Modify the forces for partially excluded interactions.
                 */
                ftm2ix = ftm2ix - fridmpx - findmpx;
                ftm2iy = ftm2iy - fridmpy - findmpy;
                ftm2iz = ftm2iz - fridmpz - findmpz;
                /*
                 * Correction to convert mutual to direct polarization force.
                 */
                if (polarization == Polarization.DIRECT) {
                    final double gfd = 0.5 * (bn2 * scip2 - bn3 * (scip3 * sci4 + sci3 * scip4));
                    final double gfdr = 0.5 * (rr5 * scip2 * usc3 - rr7 * (scip3 * sci4 + sci3 * scip4) * usc5);
                    ftm2ix = ftm2ix - gfd * xr - 0.5 * bn2 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx);
                    ftm2iy = ftm2iy - gfd * yr - 0.5 * bn2 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky);
                    ftm2iz = ftm2iz - gfd * zr - 0.5 * bn2 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz);
                    final double fdirx = gfdr * xr + 0.5 * usc5 * rr5 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx);
                    final double fdiry = gfdr * yr + 0.5 * usc5 * rr5 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky);
                    final double fdirz = gfdr * zr + 0.5 * usc5 * rr5 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz);
                    ftm2ix = ftm2ix + fdirx + findmpx;
                    ftm2iy = ftm2iy + fdiry + findmpy;
                    ftm2iz = ftm2iz + fdirz + findmpz;
                }
                /**
                 * Handle the case where scaling is used.
                 */
                ftm2ix = ftm2ix - ftm2rix;
                ftm2iy = ftm2iy - ftm2riy;
                ftm2iz = ftm2iz - ftm2riz;
                ttm2ix = ttm2ix - ttm2rix;
                ttm2iy = ttm2iy - ttm2riy;
                ttm2iz = ttm2iz - ttm2riz;
                ttm3ix = ttm3ix - ttm3rix;
                ttm3iy = ttm3iy - ttm3riy;
                ttm3iz = ttm3iz - ttm3riz;
                double scalar = ELECTRIC * polarizationScale * selfScale;
                gX[i] += scalar * ftm2ix;
                gY[i] += scalar * ftm2iy;
                gZ[i] += scalar * ftm2iz;
                tX[i] += scalar * ttm2ix;
                tY[i] += scalar * ttm2iy;
                tZ[i] += scalar * ttm2iz;
                gxk_local[k] -= scalar * ftm2ix;
                gyk_local[k] -= scalar * ftm2iy;
                gzk_local[k] -= scalar * ftm2iz;
                txk_local[k] += scalar * ttm3ix;
                tyk_local[k] += scalar * ttm3iy;
                tzk_local[k] += scalar * ttm3iz;
                if (lambdaTerm) {
                    dUdL += dEdLSign * dlPowPol * e;
                    d2UdL2 += dEdLSign * d2lPowPol * e;
                    scalar = ELECTRIC * dEdLSign * dlPowPol * selfScale;
                    lgX[i] += scalar * ftm2ix;
                    lgY[i] += scalar * ftm2iy;
                    lgZ[i] += scalar * ftm2iz;
                    ltX[i] += scalar * ttm2ix;
                    ltY[i] += scalar * ttm2iy;
                    ltZ[i] += scalar * ttm2iz;
                    lxk_local[k] -= scalar * ftm2ix;
                    lyk_local[k] -= scalar * ftm2iy;
                    lzk_local[k] -= scalar * ftm2iz;
                    ltxk_local[k] += scalar * ttm3ix;
                    ltyk_local[k] += scalar * ttm3iy;
                    ltzk_local[k] += scalar * ttm3iz;
                }
                return polarizationScale * e;
            }
        }
    }

    private class ReciprocalEnergyRegion extends ParallelRegion {

        private final double aewald1 = -ELECTRIC * aewald / sqrtPi;
        private final double aewald2 = 2.0 * aewald * aewald;
        private final double aewald3 = -2.0 / 3.0 * ELECTRIC * aewald * aewald * aewald / sqrtPi;
        private final double aewald4 = -2.0 * aewald3;
        private final double twoThirds = 2.0 / 3.0;
        private double nfftX, nfftY, nfftZ;
        private double multipole[][];
        private double ind[][];
        private double indCR[][];
        private double fracMultipoles[][];
        private double fracInd[][];
        private double fracIndCR[][];
        private double fracMultipolePhi[][];
        private double fracInducedDipolePhi[][];
        private double fracInducedDipoleCRPhi[][];
        private double permanentSelfEnergy;
        private double permanentReciprocalEnergy;
        private final SharedDouble inducedDipoleSelfEnergy;
        private final SharedDouble inducedDipoleRecipEnergy;
        private final PermanentReciprocalEnergyLoop permanentReciprocalEnergyLoop[];
        private final InducedDipoleReciprocalEnergyLoop inducedDipoleReciprocalEnergyLoop[];

        public ReciprocalEnergyRegion(int nt) {
            permanentReciprocalEnergyLoop = new PermanentReciprocalEnergyLoop[nt];
            inducedDipoleReciprocalEnergyLoop = new InducedDipoleReciprocalEnergyLoop[nt];
            inducedDipoleSelfEnergy = new SharedDouble();
            inducedDipoleRecipEnergy = new SharedDouble();
        }

        public double getPermanentSelfEnergy() {
            return permanentSelfEnergy;
        }

        public double getPermanentReciprocalEnergy() {
            return permanentReciprocalEnergy;
        }

        public double getInducedDipoleSelfEnergy() {
            return inducedDipoleSelfEnergy.get();
        }

        public double getInducedDipoleReciprocalEnergy() {
            return inducedDipoleRecipEnergy.get();
        }

        @Override
        public void start() {
            multipole = globalMultipole[0];
            ind = inducedDipole[0];
            indCR = inducedDipoleCR[0];
            fracMultipoles = reciprocalSpace.getFracMultipoles();
            fracInd = reciprocalSpace.getFracInducedDipoles();
            fracIndCR = reciprocalSpace.getFracInducedDipolesCR();
            fracMultipolePhi = reciprocalSpace.getFracMultipolePhi();
            fracInducedDipolePhi = reciprocalSpace.getFracInducedDipolePhi();
            fracInducedDipoleCRPhi = reciprocalSpace.getFracInducedDipoleCRPhi();
            inducedDipoleSelfEnergy.set(0.0);
            inducedDipoleRecipEnergy.set(0.0);
            nfftX = reciprocalSpace.getXDim();
            nfftY = reciprocalSpace.getYDim();
            nfftZ = reciprocalSpace.getZDim();
        }

        @Override
        public void run() throws Exception {
            int threadIndex = getThreadIndex();
            if (permanentReciprocalEnergyLoop[threadIndex] == null) {
                permanentReciprocalEnergyLoop[threadIndex] = new PermanentReciprocalEnergyLoop();
                inducedDipoleReciprocalEnergyLoop[threadIndex] = new InducedDipoleReciprocalEnergyLoop();
            }
            try {
                execute(0, nAtoms - 1, permanentReciprocalEnergyLoop[threadIndex]);
                if (polarization != Polarization.NONE) {
                    execute(0, nAtoms - 1, inducedDipoleReciprocalEnergyLoop[threadIndex]);
                }
            } catch (Exception e) {
                String message = "Fatal exception computing the real space field in thread " + threadIndex + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        @Override
        public void finish() {
            /**
             * The permanent multipole self energy contributions are large
             * enough that rounding differences that result from threads
             * finishing in different orders removes deterministic behavior.
             */
            permanentSelfEnergy = 0.0;
            permanentReciprocalEnergy = 0.0;
            for (int i = 0; i < maxThreads; i++) {
                permanentSelfEnergy += permanentReciprocalEnergyLoop[i].eSelf;
                permanentReciprocalEnergy += permanentReciprocalEnergyLoop[i].eRecip;
            }
        }

        private class PermanentReciprocalEnergyLoop extends IntegerForLoop {

            private double gX[], gY[], gZ[], tX[], tY[], tZ[];
            private double lgX[], lgY[], lgZ[], ltX[], ltY[], ltZ[];
            protected double eSelf;
            protected double eRecip;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eSelf = 0.0;
                eRecip = 0.0;
                int ti = getThreadIndex();
                gX = grad[ti][0];
                gY = grad[ti][1];
                gZ = grad[ti][2];
                tX = torque[ti][0];
                tY = torque[ti][1];
                tZ = torque[ti][2];
                if (lambdaTerm) {
                    lgX = lambdaGrad[ti][0];
                    lgY = lambdaGrad[ti][1];
                    lgZ = lambdaGrad[ti][2];
                    ltX = lambdaTorque[ti][0];
                    ltY = lambdaTorque[ti][1];
                    ltZ = lambdaTorque[ti][2];
                }
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                /**
                 * Permanent multipole self energy and gradient.
                 */
                for (int i = lb; i <= ub; i++) {
                    if (use[i]) {
                        double in[] = globalMultipole[0][i];
                        double cii = in[t000] * in[t000];
                        double dii = in[t100] * in[t100] + in[t010] * in[t010] + in[t001] * in[t001];
                        double qii = in[t200] * in[t200] + in[t020] * in[t020] + in[t002] * in[t002]
                                + 2.0 * (in[t110] * in[t110] + in[t101] * in[t101] + in[t011] * in[t011]);
                        eSelf += aewald1 * (cii + aewald2 * (dii / 3.0 + 2.0 * aewald2 * qii / 45.0));
                    }
                }
                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(eSelf * dlPowPerm * dEdLSign);
                    sharedd2EdLambda2.addAndGet(eSelf * d2lPowPerm * dEdLSign);
                }
                /**
                 * Permanent multipole reciprocal space energy and gradient.
                 */
                final double recip[][] = crystal.getUnitCell().A;

                /*
                 if (getThreadIndex() == 0) {
                 logger.info(String.format(" %16.8f %16.8f %16.8f", recip[0][0], recip[0][1], recip[0][2]));
                 logger.info(String.format(" %16.8f %16.8f %16.8f", recip[1][0], recip[1][1], recip[1][2]));
                 logger.info(String.format(" %16.8f %16.8f %16.8f", recip[2][0], recip[2][1], recip[2][2]));
                 } */
                double dUdL = 0.0;
                double d2UdL2 = 0.0;
                for (int i = lb; i <= ub; i++) {
                    if (use[i]) {
                        final double phi[] = cartMultipolePhi[i];
                        final double mpole[] = multipole[i];
                        final double fmpole[] = fracMultipoles[i];

                        /*
                         if (i == 0) {
                         logger.info(String.format(" %16.8f %16.8f %16.8f", phi[0], phi[1], phi[2]));
                         logger.info(String.format(" %16.8f %16.8f %16.8f", mpole[0], mpole[1], mpole[2]));
                         logger.info(String.format(" %16.8f %16.8f %16.8f", fmpole[0], fmpole[1], fmpole[2]));
                         } */
                        double e = mpole[t000] * phi[t000] + mpole[t100] * phi[t100]
                                + mpole[t010] * phi[t010] + mpole[t001] * phi[t001]
                                + oneThird * (mpole[t200] * phi[t200]
                                + mpole[t020] * phi[t020]
                                + mpole[t002] * phi[t002]
                                + 2.0 * (mpole[t110] * phi[t110]
                                + mpole[t101] * phi[t101]
                                + mpole[t011] * phi[t011]));
                        eRecip += e;
                        if (gradient || lambdaTerm) {
                            final double fPhi[] = fracMultipolePhi[i];
                            double gx = fmpole[t000] * fPhi[t100] + fmpole[t100] * fPhi[t200] + fmpole[t010] * fPhi[t110]
                                    + fmpole[t001] * fPhi[t101]
                                    + fmpole[t200] * fPhi[t300] + fmpole[t020] * fPhi[t120]
                                    + fmpole[t002] * fPhi[t102] + fmpole[t110] * fPhi[t210]
                                    + fmpole[t101] * fPhi[t201] + fmpole[t011] * fPhi[t111];
                            double gy = fmpole[t000] * fPhi[t010] + fmpole[t100] * fPhi[t110] + fmpole[t010] * fPhi[t020]
                                    + fmpole[t001] * fPhi[t011] + fmpole[t200] * fPhi[t210] + fmpole[t020] * fPhi[t030]
                                    + fmpole[t002] * fPhi[t012] + fmpole[t110] * fPhi[t120] + fmpole[t101] * fPhi[t111]
                                    + fmpole[t011] * fPhi[t021];
                            double gz = fmpole[t000] * fPhi[t001] + fmpole[t100] * fPhi[t101] + fmpole[t010] * fPhi[t011]
                                    + fmpole[t001] * fPhi[t002] + fmpole[t200] * fPhi[t201] + fmpole[t020] * fPhi[t021]
                                    + fmpole[t002] * fPhi[t003] + fmpole[t110] * fPhi[t111] + fmpole[t101] * fPhi[t102]
                                    + fmpole[t011] * fPhi[t012];
                            gx *= nfftX;
                            gy *= nfftY;
                            gz *= nfftZ;
                            final double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
                            final double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
                            final double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
                            // Compute dipole torques
                            double tqx = -mpole[t010] * phi[t001] + mpole[t001] * phi[t010];
                            double tqy = -mpole[t001] * phi[t100] + mpole[t100] * phi[t001];
                            double tqz = -mpole[t100] * phi[t010] + mpole[t010] * phi[t100];
                            // Compute quadrupole torques
                            tqx -= twoThirds * (mpole[t110] * phi[t101] + mpole[t020] * phi[t011] + mpole[t011] * phi[t002]
                                    - mpole[t101] * phi[t110] - mpole[t011] * phi[t020] - mpole[t002] * phi[t011]);
                            tqy -= twoThirds * (mpole[t101] * phi[t200] + mpole[t011] * phi[t110] + mpole[t002] * phi[t101]
                                    - mpole[t200] * phi[t101] - mpole[t110] * phi[t011] - mpole[t101] * phi[t002]);
                            tqz -= twoThirds * (mpole[t200] * phi[t110] + mpole[t110] * phi[t020] + mpole[t101] * phi[t011]
                                    - mpole[t110] * phi[t200] - mpole[t020] * phi[t110] - mpole[t011] * phi[t101]);
                            if (gradient) {
                                gX[i] += permanentScale * ELECTRIC * dfx;
                                gY[i] += permanentScale * ELECTRIC * dfy;
                                gZ[i] += permanentScale * ELECTRIC * dfz;
                                tX[i] += permanentScale * ELECTRIC * tqx;
                                tY[i] += permanentScale * ELECTRIC * tqy;
                                tZ[i] += permanentScale * ELECTRIC * tqz;
                            }
                            if (lambdaTerm) {
                                dUdL += dEdLSign * dlPowPerm * e;
                                d2UdL2 += dEdLSign * d2lPowPerm * e;
                                lgX[i] += dEdLSign * dlPowPerm * ELECTRIC * dfx;
                                lgY[i] += dEdLSign * dlPowPerm * ELECTRIC * dfy;
                                lgZ[i] += dEdLSign * dlPowPerm * ELECTRIC * dfz;
                                ltX[i] += dEdLSign * dlPowPerm * ELECTRIC * tqx;
                                ltY[i] += dEdLSign * dlPowPerm * ELECTRIC * tqy;
                                ltZ[i] += dEdLSign * dlPowPerm * ELECTRIC * tqz;
                            }
                        }

                    }
                }

                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(0.5 * dUdL * ELECTRIC);
                    sharedd2EdLambda2.addAndGet(0.5 * d2UdL2 * ELECTRIC);
                }
            }

            @Override
            public void finish() {
                eSelf *= permanentScale;
                eRecip *= permanentScale * 0.5 * ELECTRIC;
            }
        }

        private class InducedDipoleReciprocalEnergyLoop extends IntegerForLoop {

            private double eSelf;
            private double eRecip;
            private double gX[], gY[], gZ[], tX[], tY[], tZ[];
            private double lgX[], lgY[], lgZ[], ltX[], ltY[], ltZ[];
            private final double sfPhi[] = new double[tensorCount];
            private final double sPhi[] = new double[tensorCount];

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eSelf = 0.0;
                eRecip = 0.0;
                int threadID = getThreadIndex();
                gX = grad[threadID][0];
                gY = grad[threadID][1];
                gZ = grad[threadID][2];
                tX = torque[threadID][0];
                tY = torque[threadID][1];
                tZ = torque[threadID][2];
                if (lambdaTerm) {
                    lgX = lambdaGrad[threadID][0];
                    lgY = lambdaGrad[threadID][1];
                    lgZ = lambdaGrad[threadID][2];
                    ltX = lambdaTorque[threadID][0];
                    ltY = lambdaTorque[threadID][1];
                    ltZ = lambdaTorque[threadID][2];
                }
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                /**
                 * Induced dipole self energy and gradient.
                 */
                for (int i = lb; i <= ub; i++) {
                    if (use[i]) {
                        final double indi[] = ind[i];
                        final double multipolei[] = multipole[i];
                        final double dix = multipolei[t100];
                        final double diy = multipolei[t010];
                        final double diz = multipolei[t001];
                        final double dii = indi[0] * dix + indi[1] * diy + indi[2] * diz;
                        eSelf += aewald3 * dii;
                    }
                }
                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(dEdLSign * dlPowPol * eSelf);
                    sharedd2EdLambda2.addAndGet(dEdLSign * d2lPowPol * eSelf);
                }
                if (gradient) {
                    for (int i = lb; i <= ub; i++) {
                        if (use[i]) {
                            final double indi[] = ind[i];
                            final double indpi[] = indCR[i];
                            final double multipolei[] = multipole[i];
                            final double dix = multipolei[t100];
                            final double diy = multipolei[t010];
                            final double diz = multipolei[t001];
                            final double uix = 0.5 * (indi[0] + indpi[0]);
                            final double uiy = 0.5 * (indi[1] + indpi[1]);
                            final double uiz = 0.5 * (indi[2] + indpi[2]);
                            final double tix = aewald4 * (diy * uiz - diz * uiy);
                            final double tiy = aewald4 * (diz * uix - dix * uiz);
                            final double tiz = aewald4 * (dix * uiy - diy * uix);
                            tX[i] += polarizationScale * tix;
                            tY[i] += polarizationScale * tiy;
                            tZ[i] += polarizationScale * tiz;
                            if (lambdaTerm) {
                                ltX[i] += dEdLSign * dlPowPol * tix;
                                ltY[i] += dEdLSign * dlPowPol * tiy;
                                ltZ[i] += dEdLSign * dlPowPol * tiz;
                            }
                        }
                    }
                }
                /**
                 * Induced dipole reciprocal space energy and gradient.
                 */
                for (int i = lb; i <= ub; i++) {
                    if (use[i]) {
                        final double fPhi[] = fracMultipolePhi[i];
                        final double findi[] = fracInd[i];
                        final double indx = findi[0];
                        final double indy = findi[1];
                        final double indz = findi[2];
                        eRecip += indx * fPhi[t100] + indy * fPhi[t010] + indz * fPhi[t001];
                        if (gradient) {
                            final double iPhi[] = cartesianDipolePhi[i];
                            final double iCRPhi[] = cartesianDipolePhiCR[i];
                            final double fiPhi[] = fracInducedDipolePhi[i];
                            final double fiCRPhi[] = fracInducedDipoleCRPhi[i];
                            final double mpolei[] = multipole[i];
                            final double fmpolei[] = fracMultipoles[i];
                            final double findCRi[] = fracIndCR[i];
                            final double inpx = findCRi[0];
                            final double inpy = findCRi[1];
                            final double inpz = findCRi[2];
                            final double insx = indx + inpx;
                            final double insy = indy + inpy;
                            final double insz = indz + inpz;
                            for (int t = 0; t < tensorCount; t++) {
                                sPhi[t] = 0.5 * (iPhi[t] + iCRPhi[t]);
                                sfPhi[t] = fiPhi[t] + fiCRPhi[t];
                            }
                            double gx = insx * fPhi[t200] + insy * fPhi[t110] + insz * fPhi[t101];
                            double gy = insx * fPhi[t110] + insy * fPhi[t020] + insz * fPhi[t011];
                            double gz = insx * fPhi[t101] + insy * fPhi[t011] + insz * fPhi[t002];
                            if (polarization == Polarization.MUTUAL) {
                                gx += indx * fiCRPhi[t200] + inpx * fiPhi[t200] + indy * fiCRPhi[t110] + inpy * fiPhi[t110] + indz * fiCRPhi[t101] + inpz * fiPhi[t101];
                                gy += indx * fiCRPhi[t110] + inpx * fiPhi[t110] + indy * fiCRPhi[t020] + inpy * fiPhi[t020] + indz * fiCRPhi[t011] + inpz * fiPhi[t011];
                                gz += indx * fiCRPhi[t101] + inpx * fiPhi[t101] + indy * fiCRPhi[t011] + inpy * fiPhi[t011] + indz * fiCRPhi[t002] + inpz * fiPhi[t002];
                            }
                            gx += fmpolei[t000] * sfPhi[t100] + fmpolei[t100] * sfPhi[t200] + fmpolei[t010] * sfPhi[t110] + fmpolei[t001] * sfPhi[t101] + fmpolei[t200] * sfPhi[t300] + fmpolei[t020] * sfPhi[t120] + fmpolei[t002] * sfPhi[t102] + fmpolei[t110] * sfPhi[t210] + fmpolei[t101] * sfPhi[t201] + fmpolei[t011] * sfPhi[t111];
                            gy += fmpolei[t000] * sfPhi[t010] + fmpolei[t100] * sfPhi[t110] + fmpolei[t010] * sfPhi[t020] + fmpolei[t001] * sfPhi[t011] + fmpolei[t200] * sfPhi[t210] + fmpolei[t020] * sfPhi[t030] + fmpolei[t002] * sfPhi[t012] + fmpolei[t110] * sfPhi[t120] + fmpolei[t101] * sfPhi[t111] + fmpolei[t011] * sfPhi[t021];
                            gz += fmpolei[t000] * sfPhi[t001] + fmpolei[t100] * sfPhi[t101] + fmpolei[t010] * sfPhi[t011] + fmpolei[t001] * sfPhi[t002] + fmpolei[t200] * sfPhi[t201] + fmpolei[t020] * sfPhi[t021] + fmpolei[t002] * sfPhi[t003] + fmpolei[t110] * sfPhi[t111] + fmpolei[t101] * sfPhi[t102] + fmpolei[t011] * sfPhi[t012];
                            gx *= nfftX;
                            gy *= nfftY;
                            gz *= nfftZ;
                            double recip[][] = crystal.getUnitCell().A;
                            double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
                            double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
                            double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
                            dfx *= 0.5 * ELECTRIC;
                            dfy *= 0.5 * ELECTRIC;
                            dfz *= 0.5 * ELECTRIC;
                            // Compute dipole torques
                            double tqx = -mpolei[t010] * sPhi[t001] + mpolei[t001] * sPhi[t010];
                            double tqy = -mpolei[t001] * sPhi[t100] + mpolei[t100] * sPhi[t001];
                            double tqz = -mpolei[t100] * sPhi[t010] + mpolei[t010] * sPhi[t100];
                            // Compute quadrupole torques
                            tqx -= twoThirds * (mpolei[t110] * sPhi[t101] + mpolei[t020] * sPhi[t011] + mpolei[t011] * sPhi[t002] - mpolei[t101] * sPhi[t110] - mpolei[t011] * sPhi[t020] - mpolei[t002] * sPhi[t011]);
                            tqy -= twoThirds * (mpolei[t101] * sPhi[t200] + mpolei[t011] * sPhi[t110] + mpolei[t002] * sPhi[t101] - mpolei[t200] * sPhi[t101] - mpolei[t110] * sPhi[t011] - mpolei[t101] * sPhi[t002]);
                            tqz -= twoThirds * (mpolei[t200] * sPhi[t110] + mpolei[t110] * sPhi[t020] + mpolei[t101] * sPhi[t011] - mpolei[t110] * sPhi[t200] - mpolei[t020] * sPhi[t110] - mpolei[t011] * sPhi[t101]);
                            tqx *= ELECTRIC;
                            tqy *= ELECTRIC;
                            tqz *= ELECTRIC;
                            gX[i] += polarizationScale * dfx;
                            gY[i] += polarizationScale * dfy;
                            gZ[i] += polarizationScale * dfz;
                            tX[i] += polarizationScale * tqx;
                            tY[i] += polarizationScale * tqy;
                            tZ[i] += polarizationScale * tqz;
                            if (lambdaTerm) {
                                lgX[i] += dEdLSign * dlPowPol * dfx;
                                lgY[i] += dEdLSign * dlPowPol * dfy;
                                lgZ[i] += dEdLSign * dlPowPol * dfz;
                                ltX[i] += dEdLSign * dlPowPol * tqx;
                                ltY[i] += dEdLSign * dlPowPol * tqy;
                                ltZ[i] += dEdLSign * dlPowPol * tqz;
                            }
                        }
                    }
                }
                eRecip *= 0.5 * ELECTRIC;
                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(dEdLSign * dlPowPol * eRecip);
                    sharedd2EdLambda2.addAndGet(dEdLSign * d2lPowPol * eRecip);
                }
            }

            @Override
            public void finish() {
                inducedDipoleSelfEnergy.addAndGet(polarizationScale * eSelf);
                inducedDipoleRecipEnergy.addAndGet(polarizationScale * eRecip);
            }
        }
    }

    private class InitializationRegion extends ParallelRegion {

        private final InitializationLoop initializationLoop[];
        private final RotateMultipolesLoop rotateMultipolesLoop[];

        public InitializationRegion(int maxThreads) {
            initializationLoop = new InitializationLoop[maxThreads];
            rotateMultipolesLoop = new RotateMultipolesLoop[maxThreads];
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            if (initializationLoop[threadIndex] == null) {
                initializationLoop[threadIndex] = new InitializationLoop();
                rotateMultipolesLoop[threadIndex] = new RotateMultipolesLoop();
            }
            try {
                execute(0, nAtoms - 1, initializationLoop[threadIndex]);
                execute(0, nAtoms - 1, rotateMultipolesLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception initializing coordinates in thread: " + threadIndex + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class InitializationLoop extends IntegerForLoop {

            private final double in[] = new double[3];
            private final double out[] = new double[3];
            private double x[];
            private double y[];
            private double z[];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                x = coordinates[0][0];
                y = coordinates[0][1];
                z = coordinates[0][2];
                int threadID = getThreadIndex();
                if (gradient) {
                    double gX[] = grad[threadID][0];
                    double gY[] = grad[threadID][1];
                    double gZ[] = grad[threadID][2];
                    double tX[] = torque[threadID][0];
                    double tY[] = torque[threadID][1];
                    double tZ[] = torque[threadID][2];
                    fill(gX, 0.0);
                    fill(gY, 0.0);
                    fill(gZ, 0.0);
                    fill(tX, 0.0);
                    fill(tY, 0.0);
                    fill(tZ, 0.0);
                }
                if (lambdaTerm) {
                    double lgX[] = lambdaGrad[threadID][0];
                    double lgY[] = lambdaGrad[threadID][1];
                    double lgZ[] = lambdaGrad[threadID][2];
                    double ltX[] = lambdaTorque[threadID][0];
                    double ltY[] = lambdaTorque[threadID][1];
                    double ltZ[] = lambdaTorque[threadID][2];
                    fill(lgX, 0.0);
                    fill(lgY, 0.0);
                    fill(lgZ, 0.0);
                    fill(ltX, 0.0);
                    fill(ltY, 0.0);
                    fill(ltZ, 0.0);
                }
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Initialize the local coordinate arrays.
                 */
                for (int i = lb; i <= ub; i++) {
                    Atom atom = atoms[i];
                    double xyz[] = atom.getXYZ();
                    x[i] = xyz[0];
                    y[i] = xyz[1];
                    z[i] = xyz[2];
                    use[i] = atom.isActive();

                    /**
                     * Real space Ewald is cutoff at ~7 A, compared to ~12 A for
                     * vdW, so the number of neighbors is much more compact. A
                     * specific list for real space Ewald is filled during
                     * computation of the permanent real space field that
                     * includes only evaluated interactions. Subsequent real
                     * space loops, especially the SCF, then do not spend time
                     * evaluating pairwise distances outside the cutoff.
                     */
                    int size = neighborLists[0][i].length;
                    if (vaporLists != null) {
                        size = max(size, vaporLists[0][i].length);
                    }
                    if (realSpaceLists[0][i] == null || realSpaceLists[0][i].length < size) {
                        realSpaceLists[0][i] = new int[size];
                    }
                }

                /**
                 * Expand coordinates.
                 */
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                    SymOp symOp = symOps.get(iSymm);
                    double xs[] = coordinates[iSymm][0];
                    double ys[] = coordinates[iSymm][1];
                    double zs[] = coordinates[iSymm][2];
                    for (int i = lb; i <= ub; i++) {
                        in[0] = x[i];
                        in[1] = y[i];
                        in[2] = z[i];
                        crystal.applySymOp(in, out, symOp);
                        xs[i] = out[0];
                        ys[i] = out[1];
                        zs[i] = out[2];
                        int size = neighborLists[iSymm][i].length;
                        if (realSpaceLists[iSymm][i] == null || realSpaceLists[iSymm][i].length < size) {
                            realSpaceLists[iSymm][i] = new int[size];
                        }
                    }
                }
            }
        }

        private class RotateMultipolesLoop extends IntegerForLoop {

            // Local variables
            private final double localOrigin[] = new double[3];
            private final double xAxis[] = new double[3];
            private final double yAxis[] = new double[3];
            private final double zAxis[] = new double[3];
            private final double rotmat[][] = new double[3][3];
            private final double tempDipole[] = new double[3];
            private final double tempQuadrupole[][] = new double[3][3];
            private final double dipole[] = new double[3];
            private final double quadrupole[][] = new double[3][3];
            private double chargeScale, dipoleScale, quadrupoleScale;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                chargeScale = 1.0;
                dipoleScale = 1.0;
                quadrupoleScale = 1.0;
                if (!useCharges) {
                    chargeScale = 0.0;
                }
                if (!useDipoles) {
                    dipoleScale = 0.0;
                }
                if (!useQuadrupoles) {
                    quadrupoleScale = 0.0;
                }
            }

            @Override
            public void run(int lb, int ub) {
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    final double x[] = coordinates[iSymm][0];
                    final double y[] = coordinates[iSymm][1];
                    final double z[] = coordinates[iSymm][2];
                    for (int ii = lb; ii <= ub; ii++) {
                        final double in[] = localMultipole[ii];
                        final double out[] = globalMultipole[iSymm][ii];
                        double softcoreScale = 1.0;
                        if (isSoft[ii] && noSoftcoreElectrostatics) {
                            softcoreScale = 0.0;
                        }
                        if (rotateMultipoles) {
                            localOrigin[0] = x[ii];
                            localOrigin[1] = y[ii];
                            localOrigin[2] = z[ii];
                            int referenceSites[] = axisAtom[ii];
                            for (int i = 0; i < 3; i++) {
                                zAxis[i] = 0.0;
                                xAxis[i] = 0.0;
                                dipole[i] = 0.0;
                                for (int j = 0; j < 3; j++) {
                                    quadrupole[i][j] = 0.0;
                                }
                            }
                            if (referenceSites == null || referenceSites.length < 2) {
                                out[t000] = in[0] * chargeScale * softcoreScale;
                                out[t100] = 0.0;
                                out[t010] = 0.0;
                                out[t001] = 0.0;
                                out[t200] = 0.0;
                                out[t020] = 0.0;
                                out[t002] = 0.0;
                                out[t110] = 0.0;
                                out[t101] = 0.0;
                                out[t011] = 0.0;
                                PolarizeType polarizeType = atoms[ii].getPolarizeType();
                                polarizability[ii] = polarizeType.polarizability * softcoreScale;
                                continue;
                            }
                            switch (frame[ii]) {
                                case BISECTOR:
                                    int index = referenceSites[0];
                                    zAxis[0] = x[index];
                                    zAxis[1] = y[index];
                                    zAxis[2] = z[index];
                                    index = referenceSites[1];
                                    xAxis[0] = x[index];
                                    xAxis[1] = y[index];
                                    xAxis[2] = z[index];
                                    diff(zAxis, localOrigin, zAxis);
                                    norm(zAxis, zAxis);
                                    diff(xAxis, localOrigin, xAxis);
                                    norm(xAxis, xAxis);
                                    sum(xAxis, zAxis, zAxis);
                                    norm(zAxis, zAxis);
                                    rotmat[0][2] = zAxis[0];
                                    rotmat[1][2] = zAxis[1];
                                    rotmat[2][2] = zAxis[2];
                                    double dot = dot(xAxis, zAxis);
                                    scalar(zAxis, dot, zAxis);
                                    diff(xAxis, zAxis, xAxis);
                                    norm(xAxis, xAxis);
                                    rotmat[0][0] = xAxis[0];
                                    rotmat[1][0] = xAxis[1];
                                    rotmat[2][0] = xAxis[2];
                                    break;
                                case ZTHENBISECTOR:
                                    index = referenceSites[0];
                                    zAxis[0] = x[index];
                                    zAxis[1] = y[index];
                                    zAxis[2] = z[index];
                                    index = referenceSites[1];
                                    xAxis[0] = x[index];
                                    xAxis[1] = y[index];
                                    xAxis[2] = z[index];
                                    index = referenceSites[2];
                                    yAxis[0] = x[index];
                                    yAxis[1] = y[index];
                                    yAxis[2] = z[index];
                                    diff(zAxis, localOrigin, zAxis);
                                    norm(zAxis, zAxis);
                                    rotmat[0][2] = zAxis[0];
                                    rotmat[1][2] = zAxis[1];
                                    rotmat[2][2] = zAxis[2];
                                    diff(xAxis, localOrigin, xAxis);
                                    norm(xAxis, xAxis);
                                    diff(yAxis, localOrigin, yAxis);
                                    norm(yAxis, yAxis);
                                    sum(xAxis, yAxis, xAxis);
                                    norm(xAxis, xAxis);
                                    dot = dot(xAxis, zAxis);
                                    scalar(zAxis, dot, zAxis);
                                    diff(xAxis, zAxis, xAxis);
                                    norm(xAxis, xAxis);
                                    rotmat[0][0] = xAxis[0];
                                    rotmat[1][0] = xAxis[1];
                                    rotmat[2][0] = xAxis[2];
                                    break;
                                case ZTHENX:
                                default:
                                    index = referenceSites[0];
                                    zAxis[0] = x[index];
                                    zAxis[1] = y[index];
                                    zAxis[2] = z[index];
                                    index = referenceSites[1];
                                    xAxis[0] = x[index];
                                    xAxis[1] = y[index];
                                    xAxis[2] = z[index];
                                    diff(zAxis, localOrigin, zAxis);
                                    norm(zAxis, zAxis);
                                    rotmat[0][2] = zAxis[0];
                                    rotmat[1][2] = zAxis[1];
                                    rotmat[2][2] = zAxis[2];
                                    diff(xAxis, localOrigin, xAxis);
                                    dot = dot(xAxis, zAxis);
                                    scalar(zAxis, dot, zAxis);
                                    diff(xAxis, zAxis, xAxis);
                                    norm(xAxis, xAxis);
                                    rotmat[0][0] = xAxis[0];
                                    rotmat[1][0] = xAxis[1];
                                    rotmat[2][0] = xAxis[2];
                            }
                            // Finally the Y elements.
                            rotmat[0][1] = rotmat[2][0] * rotmat[1][2] - rotmat[1][0] * rotmat[2][2];
                            rotmat[1][1] = rotmat[0][0] * rotmat[2][2] - rotmat[2][0] * rotmat[0][2];
                            rotmat[2][1] = rotmat[1][0] * rotmat[0][2] - rotmat[0][0] * rotmat[1][2];
                            // Do the rotation.
                            tempDipole[0] = in[t100];
                            tempDipole[1] = in[t010];
                            tempDipole[2] = in[t001];
                            tempQuadrupole[0][0] = in[t200];
                            tempQuadrupole[1][1] = in[t020];
                            tempQuadrupole[2][2] = in[t002];
                            tempQuadrupole[0][1] = in[t110];
                            tempQuadrupole[0][2] = in[t101];
                            tempQuadrupole[1][2] = in[t011];
                            tempQuadrupole[1][0] = in[t110];
                            tempQuadrupole[2][0] = in[t101];
                            tempQuadrupole[2][1] = in[t011];

                            // Check for chiral flipping.
                            if (frame[ii] == MultipoleType.MultipoleFrameDefinition.ZTHENX
                                    && referenceSites.length == 3) {
                                localOrigin[0] = x[ii];
                                localOrigin[1] = y[ii];
                                localOrigin[2] = z[ii];
                                int index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                xAxis[0] = x[index];
                                xAxis[1] = y[index];
                                xAxis[2] = z[index];
                                index = referenceSites[2];
                                yAxis[0] = x[index];
                                yAxis[1] = y[index];
                                yAxis[2] = z[index];
                                diff(localOrigin, yAxis, localOrigin);
                                diff(zAxis, yAxis, zAxis);
                                diff(xAxis, yAxis, xAxis);
                                double c1 = zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1];
                                double c2 = xAxis[1] * localOrigin[2] - xAxis[2] * localOrigin[1];
                                double c3 = localOrigin[1] * zAxis[2] - localOrigin[2] * zAxis[1];
                                double vol = localOrigin[0] * c1 + zAxis[0] * c2 + xAxis[0] * c3;
                                if (vol < 0.0) {
                                    tempDipole[1] = -tempDipole[1];
                                    tempQuadrupole[0][1] = -tempQuadrupole[0][1];
                                    tempQuadrupole[1][0] = -tempQuadrupole[1][0];
                                    tempQuadrupole[1][2] = -tempQuadrupole[1][2];
                                    tempQuadrupole[2][1] = -tempQuadrupole[2][1];
                                }
                            }
                            for (int i = 0; i < 3; i++) {
                                double[] rotmati = rotmat[i];
                                double[] quadrupolei = quadrupole[i];
                                for (int j = 0; j < 3; j++) {
                                    double[] rotmatj = rotmat[j];
                                    dipole[i] += rotmati[j] * tempDipole[j];
                                    if (j < i) {
                                        quadrupolei[j] = quadrupole[j][i];
                                    } else {
                                        for (int k = 0; k < 3; k++) {
                                            double[] localQuadrupolek = tempQuadrupole[k];
                                            quadrupolei[j] += rotmati[k]
                                                    * (rotmatj[0] * localQuadrupolek[0]
                                                    + rotmatj[1] * localQuadrupolek[1]
                                                    + rotmatj[2] * localQuadrupolek[2]);
                                        }
                                    }
                                }
                            }
                            out[t000] = in[0] * chargeScale * softcoreScale;
                            out[t100] = dipole[0] * dipoleScale * softcoreScale;
                            out[t010] = dipole[1] * dipoleScale * softcoreScale;
                            out[t001] = dipole[2] * dipoleScale * softcoreScale;
                            out[t200] = quadrupole[0][0] * quadrupoleScale * softcoreScale;
                            out[t020] = quadrupole[1][1] * quadrupoleScale * softcoreScale;
                            out[t002] = quadrupole[2][2] * quadrupoleScale * softcoreScale;
                            out[t110] = quadrupole[0][1] * quadrupoleScale * softcoreScale;
                            out[t101] = quadrupole[0][2] * quadrupoleScale * softcoreScale;
                            out[t011] = quadrupole[1][2] * quadrupoleScale * softcoreScale;
                        } else {
                            /**
                             * No multipole rotation for isolating torque vs.
                             * non-torque pieces of the multipole energy
                             * gradient.
                             */
                            out[t000] = in[t000] * chargeScale * softcoreScale;
                            out[t100] = in[t100] * dipoleScale * softcoreScale;
                            out[t010] = in[t010] * dipoleScale * softcoreScale;
                            out[t001] = in[t001] * dipoleScale * softcoreScale;
                            out[t200] = in[t200] * quadrupoleScale * softcoreScale;
                            out[t020] = in[t020] * quadrupoleScale * softcoreScale;
                            out[t002] = in[t002] * quadrupoleScale * softcoreScale;
                            out[t110] = in[t110] * quadrupoleScale * softcoreScale;
                            out[t101] = in[t101] * quadrupoleScale * softcoreScale;
                            out[t011] = in[t011] * quadrupoleScale * softcoreScale;
                        }
                        PolarizeType polarizeType = atoms[ii].getPolarizeType();
                        polarizability[ii] = polarizeType.polarizability * softcoreScale;
                    }
                }
            }
        }
    }

    private class ExpandInducedDipolesRegion extends ParallelRegion {

        private final ExpandInducedDipoleLoop expandInducedDipoleLoop[];

        public ExpandInducedDipolesRegion(int maxThreads) {
            expandInducedDipoleLoop = new ExpandInducedDipoleLoop[maxThreads];
            for (int i = 0; i < maxThreads; i++) {
                expandInducedDipoleLoop[i] = new ExpandInducedDipoleLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, expandInducedDipoleLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception expanding coordinates in thread: " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class ExpandInducedDipoleLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) {
                for (int s = 1; s < nSymm; s++) {
                    SymOp symOp = crystal.spaceGroup.symOps.get(s);
                    for (int ii = lb; ii <= ub; ii++) {
                        crystal.applySymRot(inducedDipole[0][ii], inducedDipole[s][ii], symOp);
                        crystal.applySymRot(inducedDipoleCR[0][ii], inducedDipoleCR[s][ii], symOp);
                    }
                }
            }
        }
    }

    private class ReduceRegion extends ParallelRegion {

        private final TorqueLoop torqueLoop[];
        private final ReduceLoop reduceLoop[];

        public ReduceRegion(int threadCount) {
            torqueLoop = new TorqueLoop[threadCount];
            reduceLoop = new ReduceLoop[threadCount];
        }

        @Override
        public void run() {
            try {
                int threadIndex = getThreadIndex();
                if (torqueLoop[threadIndex] == null) {
                    torqueLoop[threadIndex] = new TorqueLoop();
                    reduceLoop[threadIndex] = new ReduceLoop();
                }
                if (rotateMultipoles) {
                    execute(0, nAtoms - 1, torqueLoop[threadIndex]);
                }
                execute(0, nAtoms - 1, reduceLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing torque in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class TorqueLoop extends IntegerForLoop {

            private final double trq[] = new double[3];
            private final double u[] = new double[3];
            private final double v[] = new double[3];
            private final double w[] = new double[3];
            private final double r[] = new double[3];
            private final double s[] = new double[3];
            private final double uv[] = new double[3];
            private final double uw[] = new double[3];
            private final double vw[] = new double[3];
            private final double ur[] = new double[3];
            private final double us[] = new double[3];
            private final double vs[] = new double[3];
            private final double ws[] = new double[3];
            private final double t1[] = new double[3];
            private final double t2[] = new double[3];
            private final double localOrigin[] = new double[3];
            private double g[][];
            private double lg[][];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                int threadID = getThreadIndex();
                g = grad[threadID];
                if (lambdaTerm) {
                    lg = lambdaGrad[threadID];
                }
            }

            @Override
            public void run(int lb, int ub) {
                if (gradient) {
                    for (int i = lb; i <= ub; i++) {
                        torque(i, torque, g);
                    }
                }
                if (lambdaTerm) {
                    for (int i = lb; i <= ub; i++) {
                        torque(i, lambdaTorque, lg);
                    }
                }
            }

            public void torque(int i, double[][][] tq, double[][] gd) {
                final int ax[] = axisAtom[i];
                // Ions, for example, have no torque.
                if (ax == null || ax.length < 2) {
                    return;
                }
                final int ia = ax[0];
                final int ib = i;
                final int ic = ax[1];
                int id = 0;
                /**
                 * Reduce the torque for atom i.
                 */
                trq[0] = tq[0][0][i];
                trq[1] = tq[0][1][i];
                trq[2] = tq[0][2][i];
                for (int j = 1; j < maxThreads; j++) {
                    trq[0] += tq[j][0][i];
                    trq[1] += tq[j][1][i];
                    trq[2] += tq[j][2][i];
                }
                double x[] = coordinates[0][0];
                double y[] = coordinates[0][1];
                double z[] = coordinates[0][2];
                localOrigin[0] = x[ib];
                localOrigin[1] = y[ib];
                localOrigin[2] = z[ib];
                u[0] = x[ia];
                u[1] = y[ia];
                u[2] = z[ia];
                v[0] = x[ic];
                v[1] = y[ic];
                v[2] = z[ic];
                // Construct the three rotation axes for the local frame
                diff(u, localOrigin, u);
                diff(v, localOrigin, v);
                switch (frame[i]) {
                    default:
                    case ZTHENX:
                    case BISECTOR:
                        cross(u, v, w);
                        break;
                    case TRISECTOR:
                    case ZTHENBISECTOR:
                        id = ax[2];
                        w[0] = x[id];
                        w[1] = y[id];
                        w[2] = z[id];
                        diff(w, localOrigin, w);
                }

                double ru = r(u);
                double rv = r(v);
                double rw = r(w);
                scalar(u, 1.0 / ru, u);
                scalar(v, 1.0 / rv, v);
                scalar(w, 1.0 / rw, w);
                // Find the perpendicular and angle for each pair of axes.
                cross(v, u, uv);
                cross(w, u, uw);
                cross(w, v, vw);
                double ruv = r(uv);
                double ruw = r(uw);
                double rvw = r(vw);
                scalar(uv, 1.0 / ruv, uv);
                scalar(uw, 1.0 / ruw, uw);
                scalar(vw, 1.0 / rvw, vw);
                // Compute the sine of the angle between the rotation axes.
                double uvcos = dot(u, v);
                double uvsin = sqrt(1.0 - uvcos * uvcos);
                //double uwcos = dot(u, w);
                //double uwsin = sqrt(1.0 - uwcos * uwcos);
                //double vwcos = dot(v, w);
                //double vwsin = sqrt(1.0 - vwcos * vwcos);
                    /*
                 * Negative of dot product of torque with unit vectors gives
                 * result of infinitesimal rotation along these vectors.
                 */
                double dphidu = -(trq[0] * u[0] + trq[1] * u[1] + trq[2] * u[2]);
                double dphidv = -(trq[0] * v[0] + trq[1] * v[1] + trq[2] * v[2]);
                double dphidw = -(trq[0] * w[0] + trq[1] * w[1] + trq[2] * w[2]);
                switch (frame[i]) {
                    case ZTHENBISECTOR:
                        // Build some additional axes needed for the Z-then-Bisector method
                        sum(v, w, r);
                        cross(u, r, s);
                        double rr = r(r);
                        double rs = r(s);
                        scalar(r, 1.0 / rr, r);
                        scalar(s, 1.0 / rs, s);
                        // Find the perpendicular and angle for each pair of axes.
                        cross(r, u, ur);
                        cross(s, u, us);
                        cross(s, v, vs);
                        cross(s, w, ws);
                        double rur = r(ur);
                        double rus = r(us);
                        double rvs = r(vs);
                        double rws = r(ws);
                        scalar(ur, 1.0 / rur, ur);
                        scalar(us, 1.0 / rus, us);
                        scalar(vs, 1.0 / rvs, vs);
                        scalar(ws, 1.0 / rws, ws);
                        // Compute the sine of the angle between the rotation axes
                        double urcos = dot(u, r);
                        double ursin = sqrt(1.0 - urcos * urcos);
                        //double uscos = dot(u, s);
                        //double ussin = sqrt(1.0 - uscos * uscos);
                        double vscos = dot(v, s);
                        double vssin = sqrt(1.0 - vscos * vscos);
                        double wscos = dot(w, s);
                        double wssin = sqrt(1.0 - wscos * wscos);
                        // Compute the projection of v and w onto the ru-plane
                        scalar(s, -vscos, t1);
                        scalar(s, -wscos, t2);
                        sum(v, t1, t1);
                        sum(w, t2, t2);
                        double rt1 = r(t1);
                        double rt2 = r(t2);
                        scalar(t1, 1.0 / rt1, t1);
                        scalar(t2, 1.0 / rt2, t2);
                        double ut1cos = dot(u, t1);
                        double ut1sin = sqrt(1.0 - ut1cos * ut1cos);
                        double ut2cos = dot(u, t2);
                        double ut2sin = sqrt(1.0 - ut2cos * ut2cos);
                        double dphidr = -(trq[0] * r[0] + trq[1] * r[1] + trq[2] * r[2]);
                        double dphids = -(trq[0] * s[0] + trq[1] * s[1] + trq[2] * s[2]);
                        for (int j = 0; j < 3; j++) {
                            double du = ur[j] * dphidr / (ru * ursin) + us[j] * dphids / ru;
                            double dv = (vssin * s[j] - vscos * t1[j]) * dphidu / (rv * (ut1sin + ut2sin));
                            double dw = (wssin * s[j] - wscos * t2[j]) * dphidu / (rw * (ut1sin + ut2sin));
                            gd[j][ia] += du;
                            gd[j][ic] += dv;
                            gd[j][id] += dw;
                            gd[j][ib] -= (du + dv + dw);
                        }
                        break;
                    case ZTHENX:
                        for (int j = 0; j < 3; j++) {
                            double du = uv[j] * dphidv / (ru * uvsin) + uw[j] * dphidw / ru;
                            double dv = -uv[j] * dphidu / (rv * uvsin);
                            gd[j][ia] += du;
                            gd[j][ic] += dv;
                            gd[j][ib] -= (du + dv);
                        }
                        break;
                    case BISECTOR:
                        for (int j = 0; j < 3; j++) {
                            double du = uv[j] * dphidv / (ru * uvsin) + 0.5 * uw[j] * dphidw / ru;
                            double dv = -uv[j] * dphidu / (rv * uvsin) + 0.5 * vw[j] * dphidw / rv;
                            gd[j][ia] += du;
                            gd[j][ic] += dv;
                            gd[j][ib] -= (du + dv);
                        }
                        break;
                    default:
                        String message = "Fatal exception: Unknown frame definition: " + frame[i] + "\n";
                        logger.log(Level.SEVERE, message);
                }

            }
        }

        private class ReduceLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                if (gradient) {
                    double gx[] = grad[0][0];
                    double gy[] = grad[0][1];
                    double gz[] = grad[0][2];
                    for (int j = 1; j < maxThreads; j++) {
                        double tx[] = grad[j][0];
                        double ty[] = grad[j][1];
                        double tz[] = grad[j][2];
                        for (int i = lb; i <= ub; i++) {
                            gx[i] += tx[i];
                            gy[i] += ty[i];
                            gz[i] += tz[i];
                        }
                    }
                    for (int i = lb; i <= ub; i++) {
                        Atom ai = atoms[i];
                        ai.addToXYZGradient(gx[i], gy[i], gz[i]);
                    }
                }
                if (lambdaTerm) {
                    double lx[] = lambdaGrad[0][0];
                    double ly[] = lambdaGrad[0][1];
                    double lz[] = lambdaGrad[0][2];
                    for (int j = 1; j < maxThreads; j++) {
                        double tx[] = lambdaGrad[j][0];
                        double ty[] = lambdaGrad[j][1];
                        double tz[] = lambdaGrad[j][2];
                        for (int i = lb; i <= ub; i++) {
                            lx[i] += tx[i];
                            ly[i] += ty[i];
                            lz[i] += tz[i];
                        }
                    }
                }
            }
        }
    }

    /**
     * Determine the real space Ewald parameters and permanent multipole self
     * energy.
     *
     * @param off Real space cutoff.
     * @param aewald Ewald convergence parameter (0.0 turns off reciprocal
     * space).
     */
    private void setEwaldParameters(double off, double aewald) {
        off2 = off * off;
        alsq2 = 2.0 * aewald * aewald;
        if (aewald <= 0.0) {
            piEwald = Double.POSITIVE_INFINITY;
        } else {
            piEwald = 1.0 / (sqrtPi * aewald);
        }
        aewald3 = 4.0 / 3.0 * pow(aewald, 3.0) / sqrtPi;
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

    /**
     * A precision of 1.0e-8 results in an Ewald coefficient that ensures
     * continuity in the real space gradient, but at the cost of increased
     * amplitudes for high frequency reciprocal space structure factors.
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

    /**
     * <p>
     * ewaldCutoff</p>
     *
     * @param coeff a double.
     * @param maxCutoff a double.
     * @param eps a double.
     * @return a double.
     */
    public static double ewaldCutoff(double coeff, double maxCutoff, double eps) {
        /*
         * Set the tolerance value; use of 1.0d-8 requires strict convergence of
         * the real Space sum.
         */
        double ratio = erfc(coeff * maxCutoff) / maxCutoff;

        if (ratio > eps) {
            return maxCutoff;
        }

        /*
         * Use a binary search to refine the coefficient.
         */
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

    public double getEwaldCutoff() {
        return off;
    }

    /**
     * Given an array of atoms (with atom types), assign multipole types and
     * reference sites.
     *
     * @param atoms List
     * @param forceField ForceField
     */
    private void assignMultipoles() {
        if (forceField == null) {
            String message = "No force field is defined.\n";
            logger.log(Level.SEVERE, message);
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
            if (!assignMultipole(i)) {
                Atom atom = atoms[i];
                String message = "No multipole could be assigned to atom:\n"
                        + atom + "\nof type:\n" + atom.getAtomType();
                logger.log(Level.SEVERE, message);
            }
        }
        /**
         * Check for multipoles that were not assigned correctly.
         */
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
                sb.append("\n" + atoms[i].toString() + "\n");
                sb.append(format("%d", i + 1));
                for (int j = 0; j < 10; j++) {
                    sb.append(format(" %8.3f", localMultipole[i][j]));
                }
                sb.append("\n");
            }
        }
        if (sb.length() > 0) {
            String message = "Fatal exception: Error assigning multipoles. " + sb.toString();
            logger.log(Level.SEVERE, message);
            System.exit(-1);
        }
    }

    private boolean assignMultipole(int i) {
        Atom atom = atoms[i];
        AtomType atomType = atoms[i].getAtomType();
        if (atomType == null) {
            String message = " Multipoles can only be assigned to atoms that have been typed.";
            logger.severe(message);
            return false;
        }

        PolarizeType polarizeType = forceField.getPolarizeType(atomType.getKey());
        if (polarizeType != null) {
            atom.setPolarizeType(polarizeType);
        } else {
            String message = " No polarization type was found for " + atom.toString();
            logger.fine(message);
            double polarizability = 0.0;
            double thole = 0.0;
            int polarizationGroup[] = null;
            polarizeType = new PolarizeType(atomType.type,
                    polarizability, thole, polarizationGroup);
            forceField.addForceFieldType(polarizeType);
            atom.setPolarizeType(polarizeType);
        }

        String key;
        // No reference atoms.
        key = atomType.getKey() + " 0 0";
        MultipoleType multipoleType = forceField.getMultipoleType(key);
        if (multipoleType != null) {
            atom.setMultipoleType(multipoleType, null);
            localMultipole[i][t000] = multipoleType.charge;
            localMultipole[i][t100] = multipoleType.dipole[0];
            localMultipole[i][t010] = multipoleType.dipole[1];
            localMultipole[i][t001] = multipoleType.dipole[2];
            localMultipole[i][t200] = multipoleType.quadrupole[0][0];
            localMultipole[i][t020] = multipoleType.quadrupole[1][1];
            localMultipole[i][t002] = multipoleType.quadrupole[2][2];
            localMultipole[i][t110] = multipoleType.quadrupole[0][1];
            localMultipole[i][t101] = multipoleType.quadrupole[0][2];
            localMultipole[i][t011] = multipoleType.quadrupole[1][2];
            axisAtom[i] = null;
            frame[i] = multipoleType.frameDefinition;
            return true;
        }

        // No bonds.
        List<Bond> bonds = atom.getBonds();
        if (bonds == null || bonds.size() < 1) {
            String message = "Multipoles can only be assigned after bonded relationships are defined.\n";
            logger.severe(message);
        }

        // 1 reference atom.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            key = atomType.getKey() + " " + atom2.getAtomType().getKey() + " 0";
            multipoleType = multipoleType = forceField.getMultipoleType(key);
            if (multipoleType != null) {
                int multipoleReferenceAtoms[] = new int[1];
                multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                atom.setMultipoleType(multipoleType, null);
                localMultipole[i][0] = multipoleType.charge;
                localMultipole[i][1] = multipoleType.dipole[0];
                localMultipole[i][2] = multipoleType.dipole[1];
                localMultipole[i][3] = multipoleType.dipole[2];
                localMultipole[i][4] = multipoleType.quadrupole[0][0];
                localMultipole[i][5] = multipoleType.quadrupole[1][1];
                localMultipole[i][6] = multipoleType.quadrupole[2][2];
                localMultipole[i][7] = multipoleType.quadrupole[0][1];
                localMultipole[i][8] = multipoleType.quadrupole[0][2];
                localMultipole[i][9] = multipoleType.quadrupole[1][2];
                axisAtom[i] = multipoleReferenceAtoms;
                frame[i] = multipoleType.frameDefinition;
                return true;
            }
        }

        // 2 reference atoms.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                key = atomType.getKey() + " " + key2 + " " + key3;
                multipoleType = forceField.getMultipoleType(key);
                if (multipoleType != null) {
                    int multipoleReferenceAtoms[] = new int[2];
                    multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                    multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                    atom.setMultipoleType(multipoleType, null);
                    localMultipole[i][0] = multipoleType.charge;
                    localMultipole[i][1] = multipoleType.dipole[0];
                    localMultipole[i][2] = multipoleType.dipole[1];
                    localMultipole[i][3] = multipoleType.dipole[2];
                    localMultipole[i][4] = multipoleType.quadrupole[0][0];
                    localMultipole[i][5] = multipoleType.quadrupole[1][1];
                    localMultipole[i][6] = multipoleType.quadrupole[2][2];
                    localMultipole[i][7] = multipoleType.quadrupole[0][1];
                    localMultipole[i][8] = multipoleType.quadrupole[0][2];
                    localMultipole[i][9] = multipoleType.quadrupole[1][2];
                    axisAtom[i] = multipoleReferenceAtoms;
                    frame[i] = multipoleType.frameDefinition;
                    return true;
                }
            }
        }

        /**
         * 3 reference atoms.
         */
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                for (Bond b3 : bonds) {
                    if (b == b3 || b2 == b3) {
                        continue;
                    }
                    Atom atom4 = b3.get1_2(atom);
                    String key4 = atom4.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        int multipoleReferenceAtoms[] = new int[3];
                        multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                        multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                        multipoleReferenceAtoms[2] = atom4.xyzIndex - 1;
                        atom.setMultipoleType(multipoleType, null);
                        localMultipole[i][0] = multipoleType.charge;
                        localMultipole[i][1] = multipoleType.dipole[0];
                        localMultipole[i][2] = multipoleType.dipole[1];
                        localMultipole[i][3] = multipoleType.dipole[2];
                        localMultipole[i][4] = multipoleType.quadrupole[0][0];
                        localMultipole[i][5] = multipoleType.quadrupole[1][1];
                        localMultipole[i][6] = multipoleType.quadrupole[2][2];
                        localMultipole[i][7] = multipoleType.quadrupole[0][1];
                        localMultipole[i][8] = multipoleType.quadrupole[0][2];
                        localMultipole[i][9] = multipoleType.quadrupole[1][2];
                        axisAtom[i] = multipoleReferenceAtoms;
                        frame[i] = multipoleType.frameDefinition;
                        return true;
                    }
                }
                List<Angle> angles = atom.getAngles();
                for (Angle angle : angles) {
                    Atom atom4 = angle.get1_3(atom);
                    if (atom4 != null) {
                        String key4 = atom4.getAtomType().getKey();
                        key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                        multipoleType = forceField.getMultipoleType(key);
                        if (multipoleType != null) {
                            int multipoleReferenceAtoms[] = new int[3];
                            multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                            multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                            multipoleReferenceAtoms[2] = atom4.xyzIndex - 1;
                            atom.setMultipoleType(multipoleType, null);
                            localMultipole[i][0] = multipoleType.charge;
                            localMultipole[i][1] = multipoleType.dipole[0];
                            localMultipole[i][2] = multipoleType.dipole[1];
                            localMultipole[i][3] = multipoleType.dipole[2];
                            localMultipole[i][4] = multipoleType.quadrupole[0][0];
                            localMultipole[i][5] = multipoleType.quadrupole[1][1];
                            localMultipole[i][6] = multipoleType.quadrupole[2][2];
                            localMultipole[i][7] = multipoleType.quadrupole[0][1];
                            localMultipole[i][8] = multipoleType.quadrupole[0][2];
                            localMultipole[i][9] = multipoleType.quadrupole[1][2];
                            axisAtom[i] = multipoleReferenceAtoms;
                            frame[i] = multipoleType.frameDefinition;
                            return true;
                        }
                    }
                }
            }
        }

        /**
         * Revert to a 2 reference atom definition that may include a 1-3 site.
         * For example a hydrogen on water.
         */
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            List<Angle> angles = atom.getAngles();
            for (Angle angle : angles) {
                Atom atom3 = angle.get1_3(atom);
                if (atom3 != null) {
                    String key3 = atom3.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        int multipoleReferenceAtoms[] = new int[2];
                        multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                        multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                        atom.setMultipoleType(multipoleType, null);
                        localMultipole[i][0] = multipoleType.charge;
                        localMultipole[i][1] = multipoleType.dipole[0];
                        localMultipole[i][2] = multipoleType.dipole[1];
                        localMultipole[i][3] = multipoleType.dipole[2];
                        localMultipole[i][4] = multipoleType.quadrupole[0][0];
                        localMultipole[i][5] = multipoleType.quadrupole[1][1];
                        localMultipole[i][6] = multipoleType.quadrupole[2][2];
                        localMultipole[i][7] = multipoleType.quadrupole[0][1];
                        localMultipole[i][8] = multipoleType.quadrupole[0][2];
                        localMultipole[i][9] = multipoleType.quadrupole[1][2];
                        axisAtom[i] = multipoleReferenceAtoms;
                        frame[i] = multipoleType.frameDefinition;
                        return true;
                    }
                    for (Angle angle2 : angles) {
                        Atom atom4 = angle2.get1_3(atom);
                        if (atom4 != null && atom4 != atom3) {
                            String key4 = atom4.getAtomType().getKey();
                            key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                            multipoleType = forceField.getMultipoleType(key);
                            if (multipoleType != null) {
                                int multipoleReferenceAtoms[] = new int[3];
                                multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                                multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                                multipoleReferenceAtoms[2] = atom4.xyzIndex - 1;
                                atom.setMultipoleType(multipoleType, null);
                                localMultipole[i][0] = multipoleType.charge;
                                localMultipole[i][1] = multipoleType.dipole[0];
                                localMultipole[i][2] = multipoleType.dipole[1];
                                localMultipole[i][3] = multipoleType.dipole[2];
                                localMultipole[i][4] = multipoleType.quadrupole[0][0];
                                localMultipole[i][5] = multipoleType.quadrupole[1][1];
                                localMultipole[i][6] = multipoleType.quadrupole[2][2];
                                localMultipole[i][7] = multipoleType.quadrupole[0][1];
                                localMultipole[i][8] = multipoleType.quadrupole[0][2];
                                localMultipole[i][9] = multipoleType.quadrupole[1][2];
                                axisAtom[i] = multipoleReferenceAtoms;
                                frame[i] = multipoleType.frameDefinition;
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }

    private void assignPolarizationGroups() {
        /**
         * Find directly connected group members for each atom.
         */
        List<Integer> group = new ArrayList<>();
        List<Integer> polarizationGroup = new ArrayList<>();
        //int g11 = 0;
        for (Atom ai : atoms) {
            group.clear();
            polarizationGroup.clear();
            Integer index = ai.getXYZIndex() - 1;
            group.add(index);
            polarizationGroup.add(ai.getType());
            PolarizeType polarizeType = ai.getPolarizeType();
            if (polarizeType != null) {
                if (polarizeType.polarizationGroup != null) {
                    for (int i : polarizeType.polarizationGroup) {
                        if (!polarizationGroup.contains(i)) {
                            polarizationGroup.add(i);
                        }
                    }
                    growGroup(polarizationGroup, group, ai);
                    Collections.sort(group);
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
                //g11 += ip11[index].length;
                //System.out.println(format("%d %d", index + 1, g11));
            } else {
                String message = "The polarize keyword was not found for atom "
                        + (index + 1) + " with type " + ai.getType();
                logger.severe(message);
            }
        }
        /**
         * Find 1-2 group relationships.
         */
        int mask[] = new int[nAtoms];
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
                ArrayList<Bond> bonds = aj.getBonds();
                for (Bond b : bonds) {
                    Atom ak = b.get1_2(aj);
                    int k = ak.getXYZIndex() - 1;
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
            Collections.sort(list);
            ip12[i] = new int[list.size()];
            int j = 0;
            for (int k : list) {
                ip12[i][j++] = k;
            }
        }
        /**
         * Find 1-3 group relationships.
         */
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
            Collections.sort(list);
            int j = 0;
            for (int k : list) {
                ip13[i][j++] = k;
            }
        }
    }

    /**
     * A recursive method that checks all atoms bonded to the seed atom for
     * inclusion in the polarization group. The method is called on each newly
     * found group member.
     *
     * @param polarizationGroup Atom types that should be included in the group.
     * @param group XYZ indeces of current group members.
     * @param seed The bonds of the seed atom are queried for inclusion in the
     * group.
     */
    private void growGroup(List<Integer> polarizationGroup,
            List<Integer> group, Atom seed) {
        List<Bond> bonds = seed.getBonds();
        for (Bond bi : bonds) {
            Atom aj = bi.get1_2(seed);
            int tj = aj.getType();
            boolean added = false;
            for (int g : polarizationGroup) {
                if (g == tj) {
                    Integer index = aj.getXYZIndex() - 1;
                    if (!group.contains(index)) {
                        group.add(index);
                        added = true;
                        break;
                    }
                }
            }
            if (added) {
                PolarizeType polarizeType = aj.getPolarizeType();
                for (int i : polarizeType.polarizationGroup) {
                    if (!polarizationGroup.contains(i)) {
                        polarizationGroup.add(i);
                    }
                }
                growGroup(polarizationGroup, group, aj);
            }
        }
    }

    private void torque(int iSymm,
            double tx[], double ty[], double tz[],
            double gx[], double gy[], double gz[],
            double origin[], double[] u,
            double v[], double w[], double uv[], double uw[],
            double vw[], double ur[], double us[], double vs[],
            double ws[], double t1[], double t2[], double r[],
            double s[]) {
        for (int i = 0; i < nAtoms; i++) {
            final int ax[] = axisAtom[i];
            // Ions, for example, have no torque.
            if (ax == null || ax.length < 2) {
                continue;
            }
            final int ia = ax[0];
            final int ib = i;
            final int ic = ax[1];
            int id = 0;
            double x[] = coordinates[iSymm][0];
            double y[] = coordinates[iSymm][1];
            double z[] = coordinates[iSymm][2];
            origin[0] = x[ib];
            origin[1] = y[ib];
            origin[2] = z[ib];
            u[0] = x[ia];
            u[1] = y[ia];
            u[2] = z[ia];
            v[0] = x[ic];
            v[1] = y[ic];
            v[2] = z[ic];
            // Construct the three rotation axes for the local frame
            diff(u, origin, u);
            diff(v, origin, v);
            switch (frame[i]) {
                default:
                case ZTHENX:
                case BISECTOR:
                    cross(u, v, w);
                    break;
                case TRISECTOR:
                case ZTHENBISECTOR:
                    id = ax[2];
                    w[0] = x[id];
                    w[1] = y[id];
                    w[2] = z[id];
                    diff(w, origin, w);
            }

            double ru = r(u);
            double rv = r(v);
            double rw = r(w);
            scalar(u, 1.0 / ru, u);
            scalar(v, 1.0 / rv, v);
            scalar(w, 1.0 / rw, w);
            // Find the perpendicular and angle for each pair of axes.
            cross(v, u, uv);
            cross(w, u, uw);
            cross(w, v, vw);
            double ruv = r(uv);
            double ruw = r(uw);
            double rvw = r(vw);
            scalar(uv, 1.0 / ruv, uv);
            scalar(uw, 1.0 / ruw, uw);
            scalar(vw, 1.0 / rvw, vw);
            // Compute the sine of the angle between the rotation axes.
            double uvcos = dot(u, v);
            double uvsin = sqrt(1.0 - uvcos * uvcos);
            //double uwcos = dot(u, w);
            //double uwsin = sqrt(1.0 - uwcos * uwcos);
            //double vwcos = dot(v, w);
            //double vwsin = sqrt(1.0 - vwcos * vwcos);
                    /*
             * Negative of dot product of torque with unit vectors gives result
             * of infinitesimal rotation along these vectors.
             */
            double dphidu = -(tx[i] * u[0] + ty[i] * u[1] + tz[i] * u[2]);
            double dphidv = -(tx[i] * v[0] + ty[i] * v[1] + tz[i] * v[2]);
            double dphidw = -(tx[i] * w[0] + ty[i] * w[1] + tz[i] * w[2]);
            switch (frame[i]) {
                case ZTHENBISECTOR:
                    // Build some additional axes needed for the Z-then-Bisector method
                    sum(v, w, r);
                    cross(u, r, s);
                    double rr = r(r);
                    double rs = r(s);
                    scalar(r, 1.0 / rr, r);
                    scalar(s, 1.0 / rs, s);
                    // Find the perpendicular and angle for each pair of axes.
                    cross(r, u, ur);
                    cross(s, u, us);
                    cross(s, v, vs);
                    cross(s, w, ws);
                    double rur = r(ur);
                    double rus = r(us);
                    double rvs = r(vs);
                    double rws = r(ws);
                    scalar(ur, 1.0 / rur, ur);
                    scalar(us, 1.0 / rus, us);
                    scalar(vs, 1.0 / rvs, vs);
                    scalar(ws, 1.0 / rws, ws);
                    // Compute the sine of the angle between the rotation axes
                    double urcos = dot(u, r);
                    double ursin = sqrt(1.0 - urcos * urcos);
                    //double uscos = dot(u, s);
                    //double ussin = sqrt(1.0 - uscos * uscos);
                    double vscos = dot(v, s);
                    double vssin = sqrt(1.0 - vscos * vscos);
                    double wscos = dot(w, s);
                    double wssin = sqrt(1.0 - wscos * wscos);
                    // Compute the projection of v and w onto the ru-plane
                    scalar(s, -vscos, t1);
                    scalar(s, -wscos, t2);
                    sum(v, t1, t1);
                    sum(w, t2, t2);
                    double rt1 = r(t1);
                    double rt2 = r(t2);
                    scalar(t1, 1.0 / rt1, t1);
                    scalar(t2, 1.0 / rt2, t2);
                    double ut1cos = dot(u, t1);
                    double ut1sin = sqrt(1.0 - ut1cos * ut1cos);
                    double ut2cos = dot(u, t2);
                    double ut2sin = sqrt(1.0 - ut2cos * ut2cos);
                    double dphidr = -(tx[i] * r[0] + ty[i] * r[1] + tz[i] * r[2]);
                    double dphids = -(tx[i] * s[0] + ty[i] * s[1] + tz[i] * s[2]);
                    for (int j = 0; j < 3; j++) {
                        double du = ur[j] * dphidr / (ru * ursin) + us[j] * dphids / ru;
                        double dv = (vssin * s[j] - vscos * t1[j]) * dphidu / (rv * (ut1sin + ut2sin));
                        double dw = (wssin * s[j] - wscos * t2[j]) * dphidu / (rw * (ut1sin + ut2sin));
                        u[j] = du;
                        v[j] = dv;
                        w[j] = dw;
                        r[j] = -du - dv - dw;
                    }
                    gx[ia] += u[0];
                    gy[ia] += u[1];
                    gz[ia] += u[2];
                    gx[ic] += v[0];
                    gy[ic] += v[1];
                    gz[ic] += v[2];
                    gx[id] += w[0];
                    gy[id] += w[1];
                    gz[id] += w[2];
                    gx[ib] += r[0];
                    gy[ib] += r[1];
                    gz[ib] += r[2];
                    break;
                case ZTHENX:
                    for (int j = 0; j < 3; j++) {
                        double du = uv[j] * dphidv / (ru * uvsin) + uw[j] * dphidw / ru;
                        double dv = -uv[j] * dphidu / (rv * uvsin);
                        u[j] = du;
                        v[j] = dv;
                        w[j] = -du - dv;
                    }
                    gx[ia] += u[0];
                    gy[ia] += u[1];
                    gz[ia] += u[2];
                    gx[ic] += v[0];
                    gy[ic] += v[1];
                    gz[ic] += v[2];
                    gx[ib] += w[0];
                    gy[ib] += w[1];
                    gz[ib] += w[2];
                    break;
                case BISECTOR:
                    for (int j = 0; j < 3; j++) {
                        double du = uv[j] * dphidv / (ru * uvsin) + 0.5 * uw[j] * dphidw / ru;
                        double dv = -uv[j] * dphidu / (rv * uvsin) + 0.5 * vw[j] * dphidw / rv;
                        u[j] = du;
                        v[j] = dv;
                        w[j] = -du - dv;
                    }
                    gx[ia] += u[0];
                    gy[ia] += u[1];
                    gz[ia] += u[2];
                    gx[ic] += v[0];
                    gy[ic] += v[1];
                    gz[ic] += v[2];
                    gx[ib] += w[0];
                    gy[ib] += w[1];
                    gz[ib] += w[2];
                    break;
                default:
                    String message = "Fatal exception: Unknown frame definition: " + frame[i] + "\n";
                    logger.log(Level.SEVERE, message);
            }
        }
    }

    /**
     * Should softcore atoms include electrostatic moments and polarizability.
     *
     * @param noElec true means moments and polarizability set to zero.
     */
    public void setNoSoftCoreElectrostatics(boolean noElec) {
        noSoftcoreElectrostatics = noElec;
    }

    /**
     * {@inheritDoc}
     *
     * Set the electrostatic lambda scaling factor.
     */
    @Override
    public void setLambda(double lambda) {
        assert (lambda >= 0.0 && lambda <= 1.0);
        this.lambda = lambda;

        if (!initSoftCore) {
            initSoftCoreInit(true);
        }

        /**
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
         */
        lAlpha = permLambdaAlpha * (1.0 - lambda) * (1.0 - lambda);
        dlAlpha = permLambdaAlpha * (1.0 - lambda);
        d2lAlpha = -permLambdaAlpha;

        lPowPerm = pow(lambda, permLambdaExponent);
        dlPowPerm = permLambdaExponent * pow(lambda, permLambdaExponent - 1.0);
        d2lPowPerm = 0.0;
        if (permLambdaExponent >= 2.0) {
            d2lPowPerm = permLambdaExponent * (permLambdaExponent - 1.0) * pow(lambda, permLambdaExponent - 2.0);
        }

        /**
         * Polarization is turned on from polarizationLambdaStart ..
         * polarizationLambdaEnd.
         */
        lPowPol = 1.0;
        dlPowPol = 0.0;
        d2lPowPol = 0.0;
        if (lambda < polLambdaStart) {
            lPowPol = 0.0;
        } else if (lambda <= polLambdaEnd) {
            double polWindow = polLambdaEnd - polLambdaStart;
            double polLambdaScale = 1.0 / polWindow;
            polLambda = polLambdaScale * (lambda - polLambdaStart);
            lPowPol = pow(polLambda, polLambdaExponent);
            if (polLambdaExponent >= 1.0) {
                dlPowPol = polLambdaExponent * pow(polLambda, polLambdaExponent - 1.0);
                if (polLambdaExponent >= 2.0) {
                    d2lPowPol = polLambdaExponent * (polLambdaExponent - 1.0)
                            * pow(polLambda, polLambdaExponent - 2.0);
                }
            }
            /**
             * Add the chain rule term due to shrinking the lambda range for the
             * polarization energy.
             */
            dlPowPol *= polLambdaScale;
            d2lPowPol *= (polLambdaScale * polLambdaScale);
        }

        if (generalizedKirkwoodTerm) {
            generalizedKirkwood.setLambda(lambda);
        }

    }

    /**
     * {@inheritDoc}
     *
     * Get the current lambda scale value.
     */
    @Override
    public double getLambda() {
        return lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        double dEdL = shareddEdLambda.get();
        if (generalizedKirkwoodTerm) {
            dEdL += generalizedKirkwood.getdEdL();
        }
        return dEdL;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        double d2EdL2 = sharedd2EdLambda2.get();
        if (generalizedKirkwoodTerm) {
            d2EdL2 += generalizedKirkwood.getd2EdL2();
        }
        return d2EdL2;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient) {
        if (!lambdaTerm) {
            logger.warning(" The lambdaterm property is false.");
            return;
        }
        /**
         * Note that the Generalized Kirkwood contributions are already in the
         * lambdaGrad array.
         */
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            gradient[index++] += lambdaGrad[0][0][i];
            gradient[index++] += lambdaGrad[0][1][i];
            gradient[index++] += lambdaGrad[0][2][i];
        }
    }

    private void computeInduceDipoleField() {
        try {
            if (nSymm > 1) {
                parallelTeam.execute(expandInducedDipolesRegion);
            }

            if (aewald > 0.0) {
                reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
            }
            sectionTeam.execute(inducedDipoleFieldRegion);
            if (aewald > 0.0) {
                reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolePhiCR);
            }
            if (generalizedKirkwoodTerm) {
                /**
                 * GK field.
                 */
                gkEnergyTotal = -System.nanoTime();
                generalizedKirkwood.computeInducedGKField();
                gkEnergyTotal += System.nanoTime();
                logger.fine(String.format(" Computed GK induced field %8.3f (sec)", gkEnergyTotal * 1.0e-9));
            }
            parallelTeam.execute(pcgRegion);
        } catch (Exception e) {
            String message = "Exception computing induced dipole field.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private int scfByPCG(boolean print, long startTime) {
        long directTime = System.nanoTime() - startTime;
        /**
         * A request of 0 SCF cycles simplifies mutual polarization to direct
         * polarization.
         */
        StringBuilder sb = null;
        if (print) {
            sb = new StringBuilder(
                    "\n Self-Consistent Field\n"
                    + " Iter  RMS Change (Debye)  Time\n");
        }
        /**
         * Find the induced dipole field due to direct dipoles (or predicted
         * induced dipoles from previous steps).
         */
        computeInduceDipoleField();

        try {
            /**
             * Set initial conjugate gradient residual (a field).
             *
             * Store the current induced dipoles and load the residual induced
             * dipole
             */
            parallelTeam.execute(pcgInitRegion1);

            /**
             * Compute preconditioner.
             */
            if (nSymm > 1) {
                parallelTeam.execute(expandInducedDipolesRegion);
            }
            parallelTeam.execute(inducedDipolePreconditionerRegion);

            /**
             * Revert to the stored induce dipoles.
             *
             * Set initial conjugate vector (induced dipoles).
             */
            parallelTeam.execute(pcgInitRegion2);
        } catch (Exception e) {
            String message = "Exception initializing preconditioned CG.";
            logger.log(Level.SEVERE, message, e);
        }

        /**
         * Conjugate gradient iteration of the mutual induced dipoles.
         */
        int completedSCFCycles = 0;
        int maxSCFCycles = 1000;
        double eps = 100.0;
        double previousEps;
        boolean done = false;
        while (!done) {
            long cycleTime = -System.nanoTime();
            /**
             * Store a copy of the current induced dipoles, then set the induced
             * dipoles to the conjugate vector.
             */
            for (int i = 0; i < nAtoms; i++) {
                vec[0][i] = inducedDipole[0][i][0];
                vec[1][i] = inducedDipole[0][i][1];
                vec[2][i] = inducedDipole[0][i][2];
                inducedDipole[0][i][0] = conj[0][i];
                inducedDipole[0][i][1] = conj[1][i];
                inducedDipole[0][i][2] = conj[2][i];
                vecCR[0][i] = inducedDipoleCR[0][i][0];
                vecCR[1][i] = inducedDipoleCR[0][i][1];
                vecCR[2][i] = inducedDipoleCR[0][i][2];
                inducedDipoleCR[0][i][0] = conjCR[0][i];
                inducedDipoleCR[0][i][1] = conjCR[1][i];
                inducedDipoleCR[0][i][2] = conjCR[2][i];
            }
            /**
             * Find the induced dipole field.
             */
            computeInduceDipoleField();

            try {
                /**
                 * Revert the induced dipoles to the saved values, then save the
                 * new residual field.
                 *
                 * Compute dot product of the conjugate vector and new residual.
                 *
                 * Reduce the residual field, add to the induced dipoles based
                 * on the scaled conjugate vector and finally set the induced
                 * dipoles to the polarizability times the residual field.
                 */
                parallelTeam.execute(pcgIterRegion1);

                /**
                 * Compute preconditioner.
                 */
                if (nSymm > 1) {
                    parallelTeam.execute(expandInducedDipolesRegion);
                }
                parallelTeam.execute(inducedDipolePreconditionerRegion);

                /**
                 * Revert the induced dipoles to the saved values.
                 *
                 * Compute the dot product of the residual and preconditioner.
                 *
                 * Update the conjugate vector and sum the square of the
                 * residual field.
                 */
                pcgIterRegion2.sum = pcgIterRegion1.sumShared.get();
                pcgIterRegion2.sumCR = pcgIterRegion1.sumCRShared.get();
                parallelTeam.execute(pcgIterRegion2);

            } catch (Exception e) {
                String message = "Exception in first CG iteration region.";
                logger.log(Level.SEVERE, message, e);
            }

            previousEps = eps;
            // eps = max(eps, epsCR);
            eps = max(pcgIterRegion2.epsShared.get(), pcgIterRegion2.epsCRShared.get());
            completedSCFCycles++;
            eps = MultipoleType.DEBYE * sqrt(eps / (double) nAtoms);
            cycleTime += System.nanoTime();
            if (print) {
                sb.append(format(
                        " %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * toSeconds));
            }
            /**
             * If the RMS Debye change increases, fail the SCF process.
             */
            if (eps > previousEps) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Fatal SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
                throw new ArithmeticException(message);
            }
            /**
             * The SCF should converge well before the max iteration check.
             * Otherwise, fail the SCF process.
             */
            if (completedSCFCycles >= maxSCFCycles) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Maximum SCF iterations reached: (%d)\n", completedSCFCycles);
                throw new ArithmeticException(message);
            }
            /**
             * Check if the convergence criteria has been achieved.
             */
            if (eps < poleps) {
                done = true;
            }
        }
        if (print) {
            sb.append(format(" Direct:                  %7.4f\n",
                    toSeconds * directTime));
            startTime = System.nanoTime() - startTime;
            sb.append(format(" Total:                   %7.4f",
                    startTime * toSeconds));
            logger.info(sb.toString());
        }

        /**
         * Find the final induced dipole field.
         */
        computeInduceDipoleField();

        return completedSCFCycles;
    }

    /**
     * Evaluate the real space field due to induced dipoles using a short cutoff
     * (~3-4 A).
     */
    private class InducedDipolePreconditionerRegion extends ParallelRegion {

        private final InducedPreconditionerFieldLoop inducedPreconditionerFieldLoop[];
        private final ReduceLoop reduceLoop[];
        private double aewaldCopy;

        public InducedDipolePreconditionerRegion(int threadCount) {
            inducedPreconditionerFieldLoop = new InducedPreconditionerFieldLoop[threadCount];
            reduceLoop = new ReduceLoop[threadCount];
        }

        @Override
        public void start() {
            // Save a copy of the Ewald parameter.
            aewaldCopy = aewald;

            // Set the Ewald parameter to a value that optimizes the preconditioner.
            aewald = preconditionerEwald;
            setEwaldParameters(off, aewald);
        }

        @Override
        public void finish() {
            // Revert the Ewald parameter.
            aewald = aewaldCopy;
            setEwaldParameters(off, aewald);
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            if (inducedPreconditionerFieldLoop[threadIndex] == null) {
                inducedPreconditionerFieldLoop[threadIndex] = new InducedPreconditionerFieldLoop();
                reduceLoop[threadIndex] = new ReduceLoop();
            }
            try {
                execute(0, nAtoms - 1, inducedPreconditionerFieldLoop[threadIndex]);
                execute(0, nAtoms - 1, reduceLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing the induced real space field in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class InducedPreconditionerFieldLoop extends IntegerForLoop {

            private double x[], y[], z[];
            private double ind[][], indCR[][];
            private double fX[], fY[], fZ[];
            private double fXCR[], fYCR[], fZCR[];

            public InducedPreconditionerFieldLoop() {
            }

            @Override
            public IntegerSchedule schedule() {
                return realSpaceSchedule;
            }

            @Override
            public void start() {
                int threadIndex = getThreadIndex();
                realSpaceSCFTime[threadIndex] -= System.nanoTime();
                fX = field[threadIndex][0];
                fY = field[threadIndex][1];
                fZ = field[threadIndex][2];
                fXCR = fieldCR[threadIndex][0];
                fYCR = fieldCR[threadIndex][1];
                fZCR = fieldCR[threadIndex][2];
                fill(fX, 0.0);
                fill(fY, 0.0);
                fill(fZ, 0.0);
                fill(fXCR, 0.0);
                fill(fYCR, 0.0);
                fill(fZCR, 0.0);
                x = coordinates[0][0];
                y = coordinates[0][1];
                z = coordinates[0][2];
                ind = inducedDipole[0];
                indCR = inducedDipoleCR[0];
            }

            @Override
            public void finish() {
                int threadIndex = getThreadIndex();
                realSpaceSCFTime[threadIndex] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                final double dx[] = new double[3];
                final double transOp[][] = new double[3][3];
                /**
                 * Loop over a chunk of atoms.
                 */
                int lists[][] = preconditionerLists[0];
                int counts[] = preconditionerCounts[0];
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double px = 0.0;
                    double py = 0.0;
                    double pz = 0.0;
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double dipolei[] = ind[i];
                    final double uix = dipolei[0];
                    final double uiy = dipolei[1];
                    final double uiz = dipolei[2];
                    final double dipoleCRi[] = indCR[i];
                    final double pix = dipoleCRi[0];
                    final double piy = dipoleCRi[1];
                    final double piz = dipoleCRi[2];
                    final double pdi = ipdamp[i];
                    final double pti = thole[i];
                    /**
                     * Loop over the neighbor list.
                     */
                    final int list[] = lists[i];
                    final int npair = counts[i];
                    for (int j = 0; j < npair; j++) {
                        final int k = list[j];
                        if (!use[k]) {
                            continue;
                        }
                        final double pdk = ipdamp[k];
                        final double ptk = thole[k];
                        dx[0] = x[k] - xi;
                        dx[1] = y[k] - yi;
                        dx[2] = z[k] - zi;
                        final double r2 = crystal.image(dx);
                        /**
                         * Calculate the error function damping terms.
                         */
                        final double r = sqrt(r2);
                        final double rr1 = 1.0 / r;
                        final double rr2 = rr1 * rr1;
                        final double ralpha = aewald * r;
                        final double exp2a = exp(-ralpha * ralpha);
                        final double bn0 = erfc(ralpha) * rr1;
                        // final double exp2a = 1.0;
                        // final double bn0 = rr1;
                        final double bn1 = (bn0 + an0 * exp2a) * rr2;
                        final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                        double scale3 = 1.0;
                        double scale5 = 1.0;
                        double damp = pdi * pdk;
                        final double pgamma = min(pti, ptk);
                        final double rdamp = r * damp;
                        damp = -pgamma * rdamp * rdamp * rdamp;
                        if (damp > -50.0) {
                            final double expdamp = exp(damp);
                            scale3 = 1.0 - expdamp;
                            scale5 = 1.0 - expdamp * (1.0 - damp);
                        }
                        double rr3 = rr1 * rr2;
                        double rr5 = 3.0 * rr3 * rr2;
                        rr3 *= (1.0 - scale3);
                        rr5 *= (1.0 - scale5);
                        final double xr = dx[0];
                        final double yr = dx[1];
                        final double zr = dx[2];
                        final double dipolek[] = ind[k];
                        final double ukx = dipolek[0];
                        final double uky = dipolek[1];
                        final double ukz = dipolek[2];
                        final double ukr = ukx * xr + uky * yr + ukz * zr;
                        final double bn2ukr = bn2 * ukr;
                        final double fimx = -bn1 * ukx + bn2ukr * xr;
                        final double fimy = -bn1 * uky + bn2ukr * yr;
                        final double fimz = -bn1 * ukz + bn2ukr * zr;
                        final double rr5ukr = rr5 * ukr;
                        final double fidx = -rr3 * ukx + rr5ukr * xr;
                        final double fidy = -rr3 * uky + rr5ukr * yr;
                        final double fidz = -rr3 * ukz + rr5ukr * zr;
                        fx += (fimx - fidx);
                        fy += (fimy - fidy);
                        fz += (fimz - fidz);
                        final double dipolepk[] = indCR[k];
                        final double pkx = dipolepk[0];
                        final double pky = dipolepk[1];
                        final double pkz = dipolepk[2];
                        final double pkr = pkx * xr + pky * yr + pkz * zr;
                        final double bn2pkr = bn2 * pkr;
                        final double pimx = -bn1 * pkx + bn2pkr * xr;
                        final double pimy = -bn1 * pky + bn2pkr * yr;
                        final double pimz = -bn1 * pkz + bn2pkr * zr;
                        final double rr5pkr = rr5 * pkr;
                        final double pidx = -rr3 * pkx + rr5pkr * xr;
                        final double pidy = -rr3 * pky + rr5pkr * yr;
                        final double pidz = -rr3 * pkz + rr5pkr * zr;
                        px += (pimx - pidx);
                        py += (pimy - pidy);
                        pz += (pimz - pidz);
                        final double uir = uix * xr + uiy * yr + uiz * zr;
                        final double bn2uir = bn2 * uir;
                        final double fkmx = -bn1 * uix + bn2uir * xr;
                        final double fkmy = -bn1 * uiy + bn2uir * yr;
                        final double fkmz = -bn1 * uiz + bn2uir * zr;
                        final double rr5uir = rr5 * uir;
                        final double fkdx = -rr3 * uix + rr5uir * xr;
                        final double fkdy = -rr3 * uiy + rr5uir * yr;
                        final double fkdz = -rr3 * uiz + rr5uir * zr;
                        fX[k] += (fkmx - fkdx);
                        fY[k] += (fkmy - fkdy);
                        fZ[k] += (fkmz - fkdz);
                        final double pir = pix * xr + piy * yr + piz * zr;
                        final double bn2pir = bn2 * pir;
                        final double pkmx = -bn1 * pix + bn2pir * xr;
                        final double pkmy = -bn1 * piy + bn2pir * yr;
                        final double pkmz = -bn1 * piz + bn2pir * zr;
                        final double rr5pir = rr5 * pir;
                        final double pkdx = -rr3 * pix + rr5pir * xr;
                        final double pkdy = -rr3 * piy + rr5pir * yr;
                        final double pkdz = -rr3 * piz + rr5pir * zr;
                        fXCR[k] += (pkmx - pkdx);
                        fYCR[k] += (pkmy - pkdy);
                        fZCR[k] += (pkmz - pkdz);
                    }
                    fX[i] += fx;
                    fY[i] += fy;
                    fZ[i] += fz;
                    fXCR[i] += px;
                    fYCR[i] += py;
                    fZCR[i] += pz;
                }
                /**
                 * Loop over symmetry mates.
                 */
                for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                    SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                    crystal.getTransformationOperator(symOp, transOp);
                    lists = preconditionerLists[iSymm];
                    counts = preconditionerCounts[iSymm];
                    final double xs[] = coordinates[iSymm][0];
                    final double ys[] = coordinates[iSymm][1];
                    final double zs[] = coordinates[iSymm][2];
                    final double inds[][] = inducedDipole[iSymm];
                    final double indCRs[][] = inducedDipoleCR[iSymm];
                    /**
                     * Loop over a chunk of atoms.
                     */
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        double fx = 0.0;
                        double fy = 0.0;
                        double fz = 0.0;
                        double px = 0.0;
                        double py = 0.0;
                        double pz = 0.0;
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final double dipolei[] = ind[i];
                        final double uix = dipolei[0];
                        final double uiy = dipolei[1];
                        final double uiz = dipolei[2];
                        final double dipoleCRi[] = indCR[i];
                        final double pix = dipoleCRi[0];
                        final double piy = dipoleCRi[1];
                        final double piz = dipoleCRi[2];
                        final double pdi = ipdamp[i];
                        final double pti = thole[i];
                        /**
                         * Loop over the neighbor list.
                         */
                        final int list[] = lists[i];
                        final int npair = counts[i];
                        for (int j = 0; j < npair; j++) {
                            final int k = list[j];
                            if (!use[k]) {
                                continue;
                            }
                            double selfScale = 1.0;
                            if (i == k) {
                                selfScale = 0.5;
                            }
                            final double pdk = ipdamp[k];
                            final double ptk = thole[k];
                            dx[0] = xs[k] - xi;
                            dx[1] = ys[k] - yi;
                            dx[2] = zs[k] - zi;
                            final double r2 = crystal.image(dx);
                            /**
                             * Calculate the error function damping terms.
                             */
                            final double r = sqrt(r2);
                            final double rr1 = 1.0 / r;
                            final double rr2 = rr1 * rr1;
                            final double ralpha = aewald * r;
                            final double exp2a = exp(-ralpha * ralpha);
                            final double bn0 = erfc(ralpha) * rr1;
                            //final double exp2a = 1.0;
                            //final double bn0 = rr1;
                            final double bn1 = (bn0 + an0 * exp2a) * rr2;
                            final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                            double scale3 = 1.0;
                            double scale5 = 1.0;
                            double damp = pdi * pdk;
                            final double pgamma = min(pti, ptk);
                            final double rdamp = r * damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                final double expdamp = exp(damp);
                                scale3 = 1.0 - expdamp;
                                scale5 = 1.0 - expdamp * (1.0 - damp);
                            }
                            double rr3 = rr1 * rr2;
                            double rr5 = 3.0 * rr3 * rr2;
                            rr3 *= (1.0 - scale3);
                            rr5 *= (1.0 - scale5);
                            final double xr = dx[0];
                            final double yr = dx[1];
                            final double zr = dx[2];
                            final double dipolek[] = inds[k];
                            final double ukx = dipolek[0];
                            final double uky = dipolek[1];
                            final double ukz = dipolek[2];
                            final double dipolepk[] = indCRs[k];
                            final double pkx = dipolepk[0];
                            final double pky = dipolepk[1];
                            final double pkz = dipolepk[2];
                            final double ukr = ukx * xr + uky * yr + ukz * zr;
                            final double bn2ukr = bn2 * ukr;
                            final double fimx = -bn1 * ukx + bn2ukr * xr;
                            final double fimy = -bn1 * uky + bn2ukr * yr;
                            final double fimz = -bn1 * ukz + bn2ukr * zr;
                            final double rr5ukr = rr5 * ukr;
                            final double fidx = -rr3 * ukx + rr5ukr * xr;
                            final double fidy = -rr3 * uky + rr5ukr * yr;
                            final double fidz = -rr3 * ukz + rr5ukr * zr;
                            fx += selfScale * (fimx - fidx);
                            fy += selfScale * (fimy - fidy);
                            fz += selfScale * (fimz - fidz);
                            final double pkr = pkx * xr + pky * yr + pkz * zr;
                            final double bn2pkr = bn2 * pkr;
                            final double pimx = -bn1 * pkx + bn2pkr * xr;
                            final double pimy = -bn1 * pky + bn2pkr * yr;
                            final double pimz = -bn1 * pkz + bn2pkr * zr;
                            final double rr5pkr = rr5 * pkr;
                            final double pidx = -rr3 * pkx + rr5pkr * xr;
                            final double pidy = -rr3 * pky + rr5pkr * yr;
                            final double pidz = -rr3 * pkz + rr5pkr * zr;
                            px += selfScale * (pimx - pidx);
                            py += selfScale * (pimy - pidy);
                            pz += selfScale * (pimz - pidz);
                            final double uir = uix * xr + uiy * yr + uiz * zr;
                            final double bn2uir = bn2 * uir;
                            final double fkmx = -bn1 * uix + bn2uir * xr;
                            final double fkmy = -bn1 * uiy + bn2uir * yr;
                            final double fkmz = -bn1 * uiz + bn2uir * zr;
                            final double rr5uir = rr5 * uir;
                            final double fkdx = -rr3 * uix + rr5uir * xr;
                            final double fkdy = -rr3 * uiy + rr5uir * yr;
                            final double fkdz = -rr3 * uiz + rr5uir * zr;
                            double xc = selfScale * (fkmx - fkdx);
                            double yc = selfScale * (fkmy - fkdy);
                            double zc = selfScale * (fkmz - fkdz);
                            fX[k] += (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                            fY[k] += (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                            fZ[k] += (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                            final double pir = pix * xr + piy * yr + piz * zr;
                            final double bn2pir = bn2 * pir;
                            final double pkmx = -bn1 * pix + bn2pir * xr;
                            final double pkmy = -bn1 * piy + bn2pir * yr;
                            final double pkmz = -bn1 * piz + bn2pir * zr;
                            final double rr5pir = rr5 * pir;
                            final double pkdx = -rr3 * pix + rr5pir * xr;
                            final double pkdy = -rr3 * piy + rr5pir * yr;
                            final double pkdz = -rr3 * piz + rr5pir * zr;
                            xc = selfScale * (pkmx - pkdx);
                            yc = selfScale * (pkmy - pkdy);
                            zc = selfScale * (pkmz - pkdz);
                            fXCR[k] += (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                            fYCR[k] += (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                            fZCR[k] += (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                        }
                        fX[i] += fx;
                        fY[i] += fy;
                        fZ[i] += fz;
                        fXCR[i] += px;
                        fYCR[i] += py;
                        fZCR[i] += pz;
                    }
                }
            }
        }

        private class ReduceLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                /**
                 * Reduce the real space field.
                 */
                for (int i = lb; i <= ub; i++) {
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double fxCR = 0.0;
                    double fyCR = 0.0;
                    double fzCR = 0.0;
                    for (int j = 1; j < maxThreads; j++) {
                        fx += field[j][0][i];
                        fy += field[j][1][i];
                        fz += field[j][2][i];
                        fxCR += fieldCR[j][0][i];
                        fyCR += fieldCR[j][1][i];
                        fzCR += fieldCR[j][2][i];
                    }
                    field[0][0][i] += fx;
                    field[0][1][i] += fy;
                    field[0][2][i] += fz;
                    fieldCR[0][0][i] += fxCR;
                    fieldCR[0][1][i] += fyCR;
                    fieldCR[0][2][i] += fzCR;
                }
            }
        }
    }

    private class PCGInitRegion1 extends ParallelRegion {

        private final PCGInitLoop pcgLoop[];

        public PCGInitRegion1(int nt) {
            pcgLoop = new PCGInitLoop[nt];
        }

        @Override
        public void run() throws Exception {
            try {
                int ti = getThreadIndex();
                if (pcgLoop[ti] == null) {
                    pcgLoop[ti] = new PCGInitLoop();
                }
                execute(0, nAtoms - 1, pcgLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class PCGInitLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {

                for (int i = lb; i <= ub; i++) {
                    /**
                     * Set initial conjugate gradient residual (a field).
                     */
                    double ipolar = 1.0 / polarizability[i];
                    rsd[0][i] = (directDipole[i][0] - inducedDipole[0][i][0]) * ipolar + field[0][0][i];
                    rsd[1][i] = (directDipole[i][1] - inducedDipole[0][i][1]) * ipolar + field[0][1][i];
                    rsd[2][i] = (directDipole[i][2] - inducedDipole[0][i][2]) * ipolar + field[0][2][i];
                    rsdCR[0][i] = (directDipoleCR[i][0] - inducedDipoleCR[0][i][0]) * ipolar + fieldCR[0][0][i];
                    rsdCR[1][i] = (directDipoleCR[i][1] - inducedDipoleCR[0][i][1]) * ipolar + fieldCR[0][1][i];
                    rsdCR[2][i] = (directDipoleCR[i][2] - inducedDipoleCR[0][i][2]) * ipolar + fieldCR[0][2][i];
                    /**
                     * Store the current induced dipoles and load the residual
                     * induced dipole
                     */
                    double polar = polarizability[i];
                    vec[0][i] = inducedDipole[0][i][0];
                    vec[1][i] = inducedDipole[0][i][1];
                    vec[2][i] = inducedDipole[0][i][2];
                    vecCR[0][i] = inducedDipoleCR[0][i][0];
                    vecCR[1][i] = inducedDipoleCR[0][i][1];
                    vecCR[2][i] = inducedDipoleCR[0][i][2];
                    inducedDipole[0][i][0] = polar * rsd[0][i];
                    inducedDipole[0][i][1] = polar * rsd[1][i];
                    inducedDipole[0][i][2] = polar * rsd[2][i];
                    inducedDipoleCR[0][i][0] = polar * rsdCR[0][i];
                    inducedDipoleCR[0][i][1] = polar * rsdCR[1][i];
                    inducedDipoleCR[0][i][2] = polar * rsdCR[2][i];
                }
            }
        }
    }

    private class PCGInitRegion2 extends ParallelRegion {

        private final PCGInitLoop pcgLoop[];

        public PCGInitRegion2(int nt) {
            pcgLoop = new PCGInitLoop[nt];
        }

        @Override
        public void run() throws Exception {
            try {
                int ti = getThreadIndex();
                if (pcgLoop[ti] == null) {
                    pcgLoop[ti] = new PCGInitLoop();
                }
                execute(0, nAtoms - 1, pcgLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class PCGInitLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {

                for (int i = lb; i <= ub; i++) {
                    /**
                     * Revert to the stored induce dipoles.
                     */
                    inducedDipole[0][i][0] = vec[0][i];
                    inducedDipole[0][i][1] = vec[1][i];
                    inducedDipole[0][i][2] = vec[2][i];
                    inducedDipoleCR[0][i][0] = vecCR[0][i];
                    inducedDipoleCR[0][i][1] = vecCR[1][i];
                    inducedDipoleCR[0][i][2] = vecCR[2][i];

                    /**
                     * Set initial conjugate vector (induced dipoles).
                     */
                    double udiag = 2.0;
                    double polar = polarizability[i];
                    rsdPre[0][i] = polar * (field[0][0][i] + udiag * rsd[0][i]);
                    rsdPre[1][i] = polar * (field[0][1][i] + udiag * rsd[1][i]);
                    rsdPre[2][i] = polar * (field[0][2][i] + udiag * rsd[2][i]);
                    rsdPreCR[0][i] = polar * (fieldCR[0][0][i] + udiag * rsdCR[0][i]);
                    rsdPreCR[1][i] = polar * (fieldCR[0][1][i] + udiag * rsdCR[1][i]);
                    rsdPreCR[2][i] = polar * (fieldCR[0][2][i] + udiag * rsdCR[2][i]);
                    conj[0][i] = rsdPre[0][i];
                    conj[1][i] = rsdPre[1][i];
                    conj[2][i] = rsdPre[2][i];
                    conjCR[0][i] = rsdPreCR[0][i];
                    conjCR[1][i] = rsdPreCR[1][i];
                    conjCR[2][i] = rsdPreCR[2][i];
                }
            }
        }
    }

    private class PCGIterRegion1 extends ParallelRegion {

        private final PCGIterLoop1 iterLoop1[];
        private final PCGIterLoop2 iterLoop2[];
        private final SharedDouble dotShared;
        private final SharedDouble dotCRShared;
        private final SharedDouble sumShared;
        private final SharedDouble sumCRShared;

        public PCGIterRegion1(int nt) {
            iterLoop1 = new PCGIterLoop1[nt];
            iterLoop2 = new PCGIterLoop2[nt];
            dotShared = new SharedDouble();
            dotCRShared = new SharedDouble();
            sumShared = new SharedDouble();
            sumCRShared = new SharedDouble();
        }

        @Override
        public void start() {
            dotShared.set(0.0);
            dotCRShared.set(0.0);
            sumShared.set(0.0);
            sumCRShared.set(0.0);
        }

        @Override
        public void run() throws Exception {
            try {
                int ti = getThreadIndex();
                if (iterLoop1[ti] == null) {
                    iterLoop1[ti] = new PCGIterLoop1();
                    iterLoop2[ti] = new PCGIterLoop2();
                }
                execute(0, nAtoms - 1, iterLoop1[ti]);
                if (ti == 0) {
                    if (dotShared.get() != 0.0) {
                        dotShared.set(sumShared.get() / dotShared.get());
                    }
                    if (dotCRShared.get() != 0.0) {
                        dotCRShared.set(sumCRShared.get() / dotCRShared.get());
                    }
                }
                barrier();
                execute(0, nAtoms - 1, iterLoop2[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class PCGIterLoop1 extends IntegerForLoop {

            public double dot;
            public double dotCR;
            public double sum;
            public double sumCR;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                dot = 0.0;
                dotCR = 0.0;
                sum = 0.0;
                sumCR = 0.0;
            }

            @Override
            public void finish() {
                dotShared.addAndGet(dot);
                dotCRShared.addAndGet(dotCR);
                sumShared.addAndGet(sum);
                sumCRShared.addAndGet(sumCR);
            }

            @Override
            public void run(int lb, int ub) throws Exception {

                for (int i = lb; i <= ub; i++) {
                    double ipolar = 1.0 / polarizability[i];
                    inducedDipole[0][i][0] = vec[0][i];
                    inducedDipole[0][i][1] = vec[1][i];
                    inducedDipole[0][i][2] = vec[2][i];
                    vec[0][i] = conj[0][i] * ipolar - field[0][0][i];
                    vec[1][i] = conj[1][i] * ipolar - field[0][1][i];
                    vec[2][i] = conj[2][i] * ipolar - field[0][2][i];
                    inducedDipoleCR[0][i][0] = vecCR[0][i];
                    inducedDipoleCR[0][i][1] = vecCR[1][i];
                    inducedDipoleCR[0][i][2] = vecCR[2][i];
                    vecCR[0][i] = conjCR[0][i] * ipolar - fieldCR[0][0][i];
                    vecCR[1][i] = conjCR[1][i] * ipolar - fieldCR[0][1][i];
                    vecCR[2][i] = conjCR[2][i] * ipolar - fieldCR[0][2][i];

                    // Compute dot product of the conjugate vector and new residual.
                    dot += conj[0][i] * vec[0][i]
                            + conj[1][i] * vec[1][i]
                            + conj[2][i] * vec[2][i];
                    dotCR += conjCR[0][i] * vecCR[0][i]
                            + conjCR[1][i] * vecCR[1][i]
                            + conjCR[2][i] * vecCR[2][i];
                    // Compute dot product of the previous residual and preconditioner.
                    sum += rsd[0][i] * rsdPre[0][i]
                            + rsd[1][i] * rsdPre[1][i]
                            + rsd[2][i] * rsdPre[2][i];
                    sumCR += rsdCR[0][i] * rsdPreCR[0][i]
                            + rsdCR[1][i] * rsdPreCR[1][i]
                            + rsdCR[2][i] * rsdPreCR[2][i];
                }

            }
        }

        private class PCGIterLoop2 extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                double dot = dotShared.get();
                double dotCR = dotCRShared.get();
                for (int i = lb; i <= ub; i++) {
                    /**
                     * Reduce the residual field, add to the induced dipoles
                     * based on the scaled conjugate vector and finally set the
                     * induced dipoles to the polarizability times the residual
                     * field.
                     */
                    rsd[0][i] -= dot * vec[0][i];
                    rsd[1][i] -= dot * vec[1][i];
                    rsd[2][i] -= dot * vec[2][i];
                    rsdCR[0][i] -= dotCR * vecCR[0][i];
                    rsdCR[1][i] -= dotCR * vecCR[1][i];
                    rsdCR[2][i] -= dotCR * vecCR[2][i];
                    vec[0][i] = inducedDipole[0][i][0] + dot * conj[0][i];
                    vec[1][i] = inducedDipole[0][i][1] + dot * conj[1][i];
                    vec[2][i] = inducedDipole[0][i][2] + dot * conj[2][i];
                    vecCR[0][i] = inducedDipoleCR[0][i][0] + dotCR * conjCR[0][i];
                    vecCR[1][i] = inducedDipoleCR[0][i][1] + dotCR * conjCR[1][i];
                    vecCR[2][i] = inducedDipoleCR[0][i][2] + dotCR * conjCR[2][i];
                    double polar = polarizability[i];
                    inducedDipole[0][i][0] = polar * rsd[0][i];
                    inducedDipole[0][i][1] = polar * rsd[1][i];
                    inducedDipole[0][i][2] = polar * rsd[2][i];
                    inducedDipoleCR[0][i][0] = polar * rsdCR[0][i];
                    inducedDipoleCR[0][i][1] = polar * rsdCR[1][i];
                    inducedDipoleCR[0][i][2] = polar * rsdCR[2][i];
                }
            }
        }
    }

    private class PCGIterRegion2 extends ParallelRegion {

        private final PCGIterLoop1 iterLoop1[];
        private final PCGIterLoop2 iterLoop2[];
        private final SharedDouble dotShared;
        private final SharedDouble dotCRShared;
        private final SharedDouble epsShared;
        private final SharedDouble epsCRShared;
        public double sum;
        public double sumCR;

        public PCGIterRegion2(int nt) {
            iterLoop1 = new PCGIterLoop1[nt];
            iterLoop2 = new PCGIterLoop2[nt];
            dotShared = new SharedDouble();
            dotCRShared = new SharedDouble();
            epsShared = new SharedDouble();
            epsCRShared = new SharedDouble();
        }

        @Override
        public void start() {
            dotShared.set(0.0);
            dotCRShared.set(0.0);
            epsShared.set(0.0);
            epsCRShared.set(0.0);
            if (sum == 0.0) {
                sum = 1.0;
            }
            if (sumCR == 0.0) {
                sumCR = 1.0;
            }
        }

        @Override
        public void run() throws Exception {
            try {
                int ti = getThreadIndex();
                if (iterLoop1[ti] == null) {
                    iterLoop1[ti] = new PCGIterLoop1();
                    iterLoop2[ti] = new PCGIterLoop2();
                }
                execute(0, nAtoms - 1, iterLoop1[ti]);
                execute(0, nAtoms - 1, iterLoop2[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class PCGIterLoop1 extends IntegerForLoop {

            public double dot;
            public double dotCR;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                dot = 0.0;
                dotCR = 0.0;
            }

            @Override
            public void finish() {
                dotShared.addAndGet(dot / sum);
                dotCRShared.addAndGet(dotCR / sumCR);
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                double udiag = 2.0;
                for (int i = lb; i <= ub; i++) {
                    /**
                     * Revert the induced dipoles to the saved values.
                     */
                    inducedDipole[0][i][0] = vec[0][i];
                    inducedDipole[0][i][1] = vec[1][i];
                    inducedDipole[0][i][2] = vec[2][i];
                    inducedDipoleCR[0][i][0] = vecCR[0][i];
                    inducedDipoleCR[0][i][1] = vecCR[1][i];
                    inducedDipoleCR[0][i][2] = vecCR[2][i];

                    /**
                     * Compute the dot product of the residual and
                     * preconditioner.
                     */
                    double polar = polarizability[i];
                    rsdPre[0][i] = polar * (field[0][0][i] + udiag * rsd[0][i]);
                    rsdPre[1][i] = polar * (field[0][1][i] + udiag * rsd[1][i]);
                    rsdPre[2][i] = polar * (field[0][2][i] + udiag * rsd[2][i]);
                    rsdPreCR[0][i] = polar * (fieldCR[0][0][i] + udiag * rsdCR[0][i]);
                    rsdPreCR[1][i] = polar * (fieldCR[0][1][i] + udiag * rsdCR[1][i]);
                    rsdPreCR[2][i] = polar * (fieldCR[0][2][i] + udiag * rsdCR[2][i]);
                    dot += rsd[0][i] * rsdPre[0][i]
                            + rsd[1][i] * rsdPre[1][i]
                            + rsd[2][i] * rsdPre[2][i];
                    dotCR += rsdCR[0][i] * rsdPreCR[0][i]
                            + rsdCR[1][i] * rsdPreCR[1][i]
                            + rsdCR[2][i] * rsdPreCR[2][i];
                }
            }
        }

        private class PCGIterLoop2 extends IntegerForLoop {

            public double eps;
            public double epsCR;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eps = 0.0;
                epsCR = 0.0;
            }

            @Override
            public void finish() {
                epsShared.addAndGet(eps);
                epsCRShared.addAndGet(epsCR);
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                double dot = dotShared.get();
                double dotCR = dotCRShared.get();
                for (int i = lb; i <= ub; i++) {
                    /**
                     * Update the conjugate vector and sum the square of the
                     * residual field.
                     */
                    conj[0][i] = rsdPre[0][i] + dot * conj[0][i];
                    conj[1][i] = rsdPre[1][i] + dot * conj[1][i];
                    conj[2][i] = rsdPre[2][i] + dot * conj[2][i];
                    conjCR[0][i] = rsdPreCR[0][i] + dotCR * conjCR[0][i];
                    conjCR[1][i] = rsdPreCR[1][i] + dotCR * conjCR[1][i];
                    conjCR[2][i] = rsdPreCR[2][i] + dotCR * conjCR[2][i];
                    eps += rsd[0][i] * rsd[0][i]
                            + rsd[1][i] * rsd[1][i]
                            + rsd[2][i] * rsd[2][i];
                    epsCR += rsdCR[0][i] * rsdCR[0][i]
                            + rsdCR[1][i] * rsdCR[1][i]
                            + rsdCR[2][i] * rsdCR[2][i];
                }
            }
        }
    }

    private class PCGRegion extends ParallelRegion {

        private final PCGLoop pcgLoop[];

        public PCGRegion(int nt) {
            pcgLoop = new PCGLoop[nt];
        }

        @Override
        public void run() throws Exception {
            try {
                int ti = getThreadIndex();
                if (pcgLoop[ti] == null) {
                    pcgLoop[ti] = new PCGLoop();
                }
                execute(0, nAtoms - 1, pcgLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }

        }

        private class PCGLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                final double induced0[][] = inducedDipole[0];
                final double inducedCR0[][] = inducedDipoleCR[0];
                /**
                 * Reduce the real space field.
                 */
                for (int i = lb; i <= ub; i++) {
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double fxCR = 0.0;
                    double fyCR = 0.0;
                    double fzCR = 0.0;
                    for (int j = 1; j < maxThreads; j++) {
                        fx += field[j][0][i];
                        fy += field[j][1][i];
                        fz += field[j][2][i];
                        fxCR += fieldCR[j][0][i];
                        fyCR += fieldCR[j][1][i];
                        fzCR += fieldCR[j][2][i];
                    }
                    field[0][0][i] += fx;
                    field[0][1][i] += fy;
                    field[0][2][i] += fz;
                    fieldCR[0][0][i] += fxCR;
                    fieldCR[0][1][i] += fyCR;
                    fieldCR[0][2][i] += fzCR;
                }
                if (aewald > 0.0) {
                    /**
                     * Add the self and reciprocal space fields to the real
                     * space field.
                     */
                    for (int i = lb; i <= ub; i++) {
                        double dipolei[] = induced0[i];
                        double dipoleCRi[] = inducedCR0[i];
                        final double phii[] = cartesianDipolePhi[i];
                        final double phiCRi[] = cartesianDipolePhiCR[i];
                        double fx = aewald3 * dipolei[0] - phii[t100];
                        double fy = aewald3 * dipolei[1] - phii[t010];
                        double fz = aewald3 * dipolei[2] - phii[t001];
                        double fxCR = aewald3 * dipoleCRi[0] - phiCRi[t100];
                        double fyCR = aewald3 * dipoleCRi[1] - phiCRi[t010];
                        double fzCR = aewald3 * dipoleCRi[2] - phiCRi[t001];
                        field[0][0][i] += fx;
                        field[0][1][i] += fy;
                        field[0][2][i] += fz;
                        fieldCR[0][0][i] += fxCR;
                        fieldCR[0][1][i] += fyCR;
                        fieldCR[0][2][i] += fzCR;
                    }
                }
                if (generalizedKirkwoodTerm) {
                    SharedDoubleArray gkField[] = generalizedKirkwood.sharedGKField;
                    SharedDoubleArray gkFieldCR[] = generalizedKirkwood.sharedGKFieldCR;
                    /**
                     * Add the GK reaction field to the intramolecular field.
                     */
                    for (int i = lb; i <= ub; i++) {
                        field[0][0][i] += gkField[0].get(i);
                        field[0][1][i] += gkField[1].get(i);
                        field[0][2][i] += gkField[2].get(i);
                        fieldCR[0][0][i] += gkFieldCR[0].get(i);
                        fieldCR[0][1][i] += gkFieldCR[1].get(i);
                        fieldCR[0][2][i] += gkFieldCR[2].get(i);
                    }
                }
            }
        }
    }

    /**
     * Save the current converged mutual induced dipoles.
     */
    private void saveMutualInducedDipoles() {

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
                predictorInducedDipole[mode][predictorStartIndex][i][j]
                        = inducedDipole[0][i][j] - directDipole[i][j];
                predictorInducedDipoleCR[mode][predictorStartIndex][i][j]
                        = inducedDipoleCR[0][i][j] - directDipoleCR[i][j];
            }
        }
    }

    /**
     * The least-squares predictor with induced dipole information from 8-10
     * previous steps reduces the number SCF iterations by ~50%.
     */
    private void leastSquaresPredictor() {
        if (predictorCount < 2) {
            return;
        }
        try {
            /**
             * The Jacobian and target do not change during the LS optimization,
             * so it's most efficient to update them once before the
             * Least-Squares optimizer starts.
             */
            leastSquaresPredictor.updateJacobianAndTarget();
            int maxEvals = 100;
            fill(leastSquaresPredictor.initialSolution, 0.0);
            leastSquaresPredictor.initialSolution[0] = 1.0;
            PointVectorValuePair optimum
                    = leastSquaresOptimizer.optimize(maxEvals,
                            leastSquaresPredictor,
                            leastSquaresPredictor.calculateTarget(),
                            leastSquaresPredictor.weights,
                            leastSquaresPredictor.initialSolution);
            double[] optimalValues = optimum.getPoint();
            if (logger.isLoggable(Level.FINEST)) {
                logger.finest(String.format("\n LS RMS:            %10.6f", leastSquaresOptimizer.getRMS()));
                logger.finest(String.format(" LS Iterations:     %10d", leastSquaresOptimizer.getEvaluations()));
                logger.finest(String.format(" Jacobian Evals:    %10d", leastSquaresOptimizer.getJacobianEvaluations()));
                logger.finest(String.format(" Chi Square:        %10.6f", leastSquaresOptimizer.getChiSquare()));
                logger.finest(String.format(" LS Coefficients"));
                for (int i = 0; i < predictorOrder - 1; i++) {
                    logger.finest(String.format(" %2d  %10.6f", i + 1, optimalValues[i]));
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

            /**
             * Initialize a pointer into predictor induced dipole array.
             */
            int index = predictorStartIndex;
            /**
             * Apply the LS coefficients in order to provide an initial guess at
             * the converged induced dipoles.
             */
            for (int k = 0; k < predictorOrder - 1; k++) {
                /**
                 * Set the current coefficient.
                 */
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

    private class LeastSquaresPredictor
            implements DifferentiableMultivariateVectorFunction {

        double weights[];
        double target[];
        double values[];
        double jacobian[][];
        double initialSolution[];

        public LeastSquaresPredictor() {
            weights = new double[2 * nAtoms * 3];
            target = new double[2 * nAtoms * 3];
            values = new double[2 * nAtoms * 3];
            jacobian = new double[2 * nAtoms * 3][predictorOrder - 1];
            initialSolution = new double[predictorOrder - 1];
            fill(weights, 1.0);
            initialSolution[0] = 1.0;
        }

        public double[] calculateTarget() {
            return target;
        }

        public void updateJacobianAndTarget() {
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

        @Override
        public MultivariateMatrixFunction jacobian() {
            return multivariateMatrixFunction;
        }
        private MultivariateMatrixFunction multivariateMatrixFunction
                = new MultivariateMatrixFunction() {
                    @Override
                    public double[][] value(double[] point) {
                        return jacobian(point);
                    }
                };
    }

    /**
     * Always-stable predictor-corrector for the mutual induced dipoles.
     */
    private void aspcPredictor() {

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

        final double aspc[] = {22.0 / 7.0, -55.0 / 14.0, 55.0 / 21.0, -22.0 / 21.0, 5.0 / 21.0, -1.0 / 42.0};
        /**
         * Initialize a pointer into predictor induced dipole array.
         */
        int index = predictorStartIndex;
        /**
         * Expansion loop.
         */
        for (int k = 0; k < 6; k++) {
            /**
             * Set the current predictor coefficient.
             */
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

    private class PolynomialPredictor extends ParallelRegion {

        public PolynomialPredictor() {
        }

        @Override
        public void run() throws Exception {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        private class PolynomialPredictorLoop extends ParallelForLoop {

            public PolynomialPredictorLoop() {
            }
        }
    }

    /**
     * Polynomial predictor for the mutual induced dipoles.
     */
    private void polynomialPredictor() {

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

        /**
         * Check the number of previous induced dipole vectors available.
         */
        int n = predictorOrder;
        if (predictorCount < predictorOrder) {
            n = predictorCount;
        }
        /**
         * Initialize a pointer into predictor induced dipole array.
         */
        int index = predictorStartIndex;
        /**
         * Initialize the sign of the polynomial expansion.
         */
        double sign = -1.0;
        /**
         * Expansion loop.
         */
        for (int k = 0; k < n; k++) {
            /**
             * Set the current predictor sign and coefficient.
             */
            sign *= -1.0;
            double c = sign * VectorMath.binomial(n, k);
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
        logger.info(String.format("%s %6d-%s %6d-%s %10.4f  %10.4f",
                "ELEC", atoms[i].xyzIndex, atoms[i].getAtomType().name,
                atoms[k].xyzIndex, atoms[k].getAtomType().name, r, eij));
    }
    /**
     * Number of unique tensors for given order.
     */
    private static final int tensorCount = TensorRecursion.tensorCount(3);
    private static final double oneThird = 1.0 / 3.0;

}
