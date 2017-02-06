/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
import ffx.numerics.MultipoleTensor;
import ffx.numerics.MultipoleTensor.COORDINATES;
import ffx.numerics.MultipoleTensor.OPERATOR;
import ffx.numerics.VectorMath;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Torsion;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.nonbonded.ReciprocalSpace.FFTMethod;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.ForceField.ForceFieldType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.utils.EnergyException;

import static ffx.numerics.Erf.erfc;
import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.extended.ExtUtils.DebugHandler.DEBUG;
import static ffx.potential.extended.ExtUtils.DebugHandler.VERBOSE;
import static ffx.potential.extended.ExtUtils.DebugHandler.debugIntI;
import static ffx.potential.extended.ExtUtils.DebugHandler.debugIntK;
import static ffx.potential.extended.ExtUtils.formatArray;
import static ffx.potential.extended.ExtUtils.logf;
import static ffx.potential.parameters.MultipoleType.ELECTRIC;
import static ffx.potential.parameters.MultipoleType.checkMultipoleChirality;
import static ffx.potential.parameters.MultipoleType.getRotationMatrix;
import static ffx.potential.parameters.MultipoleType.rotateMultipole;
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
public class ParticleMeshEwaldQI extends ParticleMeshEwald implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(ParticleMeshEwald.class.getName());

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
     * Extended System variables
     */
    private boolean esvTerm = false;
    private int numESVs = 0;
    private boolean esvAtoms[] = null;  // [nAtoms]
    private double esvLambda[] = null;  // [nAtoms]
    private int atomMultistates[] = null;  // [for on-the-fly DualTop handling]
    private SharedDouble[] esvRealSpaceDeriv = null;     // [numESVs]
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
     * Lambda and Extended state variables.
     */
    private enum LambdaMode {

        OFF, CONDENSED, CONDENSED_NO_LIGAND, VAPOR
    };
    private LambdaMode lambdaMode = LambdaMode.OFF;
    /**
     * Current state.
     */
    private double lambda = 1.0;
    private ExtendedSystem esvSystem = null;
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
     * Specify inter-molecular softcore.
     */
    private boolean intermolecularSoftcore = false;
    /**
     * Specify intra-molecular softcore.
     */
    private boolean intramolecularSoftcore = false;
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
     * Flag for ligand atoms; treats both OSRW and ESV lambdas.
     */
    private boolean isSoft[];
    /**
     * Flag indicating if softcore variables have been initialized.
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
    private final boolean gpuFFT;
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

    private final boolean reciprocalSpaceTerm;
    private final ReciprocalSpace reciprocalSpace;
    private final ReciprocalEnergyRegion reciprocalEnergyRegion;
    private final RealSpaceEnergyRegionQI realSpaceEnergyRegionQI;
    private final ReduceRegion reduceRegion;
    private final GeneralizedKirkwood generalizedKirkwood;
    /**
     * Timing variables.
     */
    private final long realSpacePermTime[];
    private final long realSpaceEnergyTime[];
    private final long realSpaceSCFTime[];
    private long realSpacePermTotalQI, realSpaceEnergyTotalQI, realSpaceSCFTotalQI;
    private long bornRadiiTotal, gkEnergyTotal;
    private ELEC_FORM elecForm = ELEC_FORM.PAM;
    private static final double TO_SECONDS = 1.0e-9;
    private static final double TO_MS = 1.0e-6;
    /**
     * Debug flags.
     */
    // private final int DEBUG() = (import from ExtConstants).
    private final COORDINATES bufferCoords;

    /**
     * The sqrt of PI.
     */
    private static final double SQRT_PI = sqrt(Math.PI);

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
    public ParticleMeshEwaldQI(Atom atoms[], int molecule[], ForceField forceField,
            Crystal crystal, NeighborList neighborList, ELEC_FORM elecForm, ParallelTeam parallelTeam) {
        /* REM
            Used to require the dlAlphaMode == FACTORED.
            ie. dlAlpha /= -2.0, d2lAlpha /= -2.0;
        */
        bufferCoords = COORDINATES.QI;

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
        maxThreads = parallelTeam.getThreadCount() + 1;

        polsor = forceField.getDouble(ForceFieldDouble.POLAR_SOR, 0.70);
        poleps = forceField.getDouble(ForceFieldDouble.POLAR_EPS, 1e-5);
        if (elecForm == ELEC_FORM.PAM) {
            m12scale = forceField.getDouble(ForceFieldDouble.MPOLE_12_SCALE, 0.0);
            m13scale = forceField.getDouble(ForceFieldDouble.MPOLE_13_SCALE, 0.0);
            m14scale = forceField.getDouble(ForceFieldDouble.MPOLE_14_SCALE, 0.4);
            m15scale = forceField.getDouble(ForceFieldDouble.MPOLE_15_SCALE, 0.8);
        } else {
            double mpole14 = 0.5;
            String name = forceField.toString().toUpperCase();
            if (name.contains("AMBER")) {
                mpole14 = 1.0 / 1.2;
            }
            m12scale = forceField.getDouble(ForceFieldDouble.MPOLE_12_SCALE, 0.0);
            m13scale = forceField.getDouble(ForceFieldDouble.MPOLE_13_SCALE, 0.0);
            m14scale = forceField.getDouble(ForceFieldDouble.MPOLE_14_SCALE, mpole14);
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
            off = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 1000.0);
        }
        double ewaldPrecision = forceField.getDouble(ForceFieldDouble.EWALD_PRECISION, 1.0e-8);
        aewald = forceField.getDouble(ForceFieldDouble.EWALD_ALPHA, ewaldCoefficient(off, ewaldPrecision));
        setEwaldParameters(off, aewald);

        reciprocalSpaceTerm = forceField.getBoolean(ForceFieldBoolean.RECIPTERM, true);

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

        if (lambdaTerm || esvTerm) {
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
            intramolecularSoftcore = forceField.getBoolean(
                    ForceField.ForceFieldBoolean.INTRAMOLECULAR_SOFTCORE, false);
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

        String temp = forceField.getString(ForceField.ForceFieldString.FFT_METHOD, "PJ");
        FFTMethod method;
        try {
            method = ReciprocalSpace.FFTMethod.valueOf(temp.toUpperCase().trim());
        } catch (Exception e) {
            method = ReciprocalSpace.FFTMethod.PJ;
        }
        gpuFFT = method != FFTMethod.PJ;

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

        StringBuilder config = new StringBuilder();
        config.append(format("\n Quasi-Internal PME Settings\n"));
        config.append(format("  Debug,Verbose,dbgI&K: %5b %5b %5d %5d\n", 
                DEBUG(), VERBOSE(), debugIntI().orElse(-1), debugIntK().orElse(-1)));
        config.append(format("  Buffer Coords:        %5s\n",
                bufferCoords.toString()));
        config.append(format("  Chrg,Dipl,Quad:       %5b %5b %5b\n",
                useCharges, useDipoles, useQuadrupoles));
        logger.info(config.toString());

        if (gpuFFT) {
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
        if (aewald > 0.0 && reciprocalSpaceTerm) {
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
        realSpaceEnergyRegionQI = new RealSpaceEnergyRegionQI(maxThreads);
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
        if (esvTerm) {
            StringBuilder sb = new StringBuilder(" ESV Parameters\n");
            sb.append(format(" Permanent Multipole Softcore Alpha:      %5.3f\n", permLambdaAlpha));
            sb.append(format(" Permanent Multipole lambda exponent constrained to unity.\n"));
            if (polarization != Polarization.NONE) {
                throw new UnsupportedOperationException();
            }
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
                if (lambdaTerm || esvTerm) {
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
            esvAtoms = new boolean[nAtoms];  // Needed regardless of esvTerm state.
            fill(esvAtoms, false);
            if (esvTerm) {
                numESVs = esvSystem.n();
                esvRealSpaceDeriv = new SharedDouble[numESVs];
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
     * Pass in atoms that have been assigned electrostatics from a fixed charge
     * force field.
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
    private void initSoftCore(boolean rebuild, boolean print) {
        if (initSoftCore && !rebuild) {
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
        
        if (esvTerm) {
            sb = new StringBuilder("\n ESV-PME Softcore:\n");
            for (int i = 0; i < nAtoms; i++) {
                // Only add softcores due to ESVs; don't interfere with existing.
                if (esvSystem.isExtAll(i)) {
                    isSoft[i] = true;
                    sb.append(atoms[i].toString()).append("\n");
                }
            }
            if (print) {
                logger.info(sb.toString());
            }
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
            vaporPermanentSchedule = vacuumNeighborList.getPairwiseSchedule();
            vaporEwaldSchedule = vaporPermanentSchedule;
            vacuumRanges = new Range[maxThreads];
            neighborList.setDisableUpdates(forceField.getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false));
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

    @Override
    public void setAtoms(Atom atoms[], int molecule[]) {
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
    @Override
    public double energy(boolean gradient, boolean print) {
        this.gradient = gradient;

        /**
         * Initialize energy, interaction, and timing variables.
         */
        permanentMultipoleEnergy = 0.0;
        polarizationEnergy = 0.0;
        generalizedKirkwoodEnergy = 0.0;
        interactions = 0;
        gkInteractions = 0;
        for (int i = 0; i < maxThreads; i++) {
            realSpacePermTime[i] = 0;
            realSpaceEnergyTime[i] = 0;
            realSpaceSCFTime[i] = 0;
        }
        realSpacePermTotalQI = 0;
        realSpaceEnergyTotalQI = 0;
        realSpaceSCFTotalQI = 0;
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
        if (esvTerm) {
            fill(esvRealSpaceDeriv, 0.0);
        }
        doPermanentRealSpace = true;
        permanentScale = 1.0;
        doPolarization = true;
        polarizationScale = 1.0;

        /**
         * Expand the coordinates and rotate multipoles into the global frame.
         */
        try {
            parallelTeam.execute(initializationRegion);
        } catch (Exception e) {
            String message = "Fatal exception expanding coordinates and rotating multipoles.\n";
            logger.log(Level.SEVERE, message, e);
        }

        double energyLog;
        if (!lambdaTerm && !esvTerm) {
            lambdaMode = LambdaMode.OFF;
            computeEnergy(print);
        } else {
            /**
             * Condensed phase with all atoms.
             */
            lambdaMode = LambdaMode.CONDENSED;
            energyLog = condensedEnergy();
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Solvated energy: %20.8f", energyLog));
            }

            /**
             * Condensed phase SCF without ligand atoms.
             */
            if (doNoLigandCondensedSCF) {
                lambdaMode = LambdaMode.CONDENSED_NO_LIGAND;
                double previous = energyLog;
                energyLog = condensedNoLigandSCF();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(String.format(" Step 2 energy:   %20.8f", energyLog - previous));
                }
            }

            /**
             * Vapor ligand electrostatics.
             */
            if (doLigandVaporElec) {
                lambdaMode = LambdaMode.VAPOR;
                double previous = energyLog;
                energyLog = vaporElec();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(String.format(" Vacuum energy:   %20.8f", energyLog - previous));
                }
            }
        }

        /**
         * Convert torques to gradients on multipole frame defining atoms. Add
         * to electrostatic gradient to the total XYZ gradient.
         */
        if (gradient || lambdaTerm || esvTerm) {
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
            if (aewald > 0.0 && reciprocalSpaceTerm) {
                reciprocalSpace.printTimings();
            }
        }

        if (polarizationEnergy != 0.0) {
            logger.warning(format("QI doesn't yet support polarization.\n"
                    + "    Non-zero polarization energy: %g", polarizationEnergy));
        }
        return permanentMultipoleEnergy + polarizationEnergy;
    }

    private void printRealSpaceTimings() {
        double total = (realSpacePermTotalQI + realSpaceSCFTotalQI + realSpaceEnergyTotalQI) * TO_SECONDS;

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
            int count = realSpaceEnergyRegionQI.realSpaceEnergyLoops[i].getCount();
            logger.info(String.format("    %3d   %7.4f %7.4f %7.4f %10d", i,
                    realSpacePermTime[i] * TO_SECONDS, realSpaceSCFTime[i] * TO_SECONDS,
                    realSpaceEnergyTime[i] * TO_SECONDS, count));
            minPerm = min(realSpacePermTime[i], minPerm);
            maxPerm = max(realSpacePermTime[i], maxPerm);
            minSCF = min(realSpaceSCFTime[i], minSCF);
            maxSCF = max(realSpaceSCFTime[i], maxSCF);
            minEnergy = min(realSpaceEnergyTime[i], minEnergy);
            maxEnergy = max(realSpaceEnergyTime[i], maxEnergy);
            minCount = min(count, minCount);
            maxCount = max(count, maxCount);
        }
        int inter = realSpaceEnergyRegionQI.getInteractions();
        logger.info(String.format(" Min      %7.4f %7.4f %7.4f %10d",
                minPerm * TO_SECONDS, minSCF * TO_SECONDS,
                minEnergy * TO_SECONDS, minCount));
        logger.info(String.format(" Max      %7.4f %7.4f %7.4f %10d",
                maxPerm * TO_SECONDS, maxSCF * TO_SECONDS,
                maxEnergy * TO_SECONDS, maxCount));
        logger.info(String.format(" Delta    %7.4f %7.4f %7.4f %10d",
                (maxPerm - minPerm) * TO_SECONDS, (maxSCF - minSCF) * TO_SECONDS,
                (maxEnergy - minEnergy) * TO_SECONDS, (maxCount - minCount)));
        logger.info(String.format(" Actual   %7.4f %7.4f %7.4f %10d",
                realSpacePermTotalQI * TO_SECONDS, realSpaceSCFTotalQI * TO_SECONDS,
                realSpaceEnergyTotalQI * TO_SECONDS, inter));
    }

    /**
     * 1.) Total system under PBC.
     *      A.) Softcore real space for Ligand-Protein and Ligand-Ligand.
     *      B.) Reciprocal space scaled by lambda.
     *      C.) Polarization scaled by lambda.
     */
    private double condensedEnergy() {
        if (!esvTerm) {
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
        } else {    // ESVs present
            double largestLambda = 0.0;
            for (ExtendedVariable esv : esvSystem) {
                largestLambda = Math.max(esv.getLambda(), largestLambda);
            }
            if (largestLambda < polLambdaStart) {
                // Skip polarization if(f) all lambdas are below its starting point.
                doPolarization = false;
                polarizationScale = 0.0;
            } else if (largestLambda <= polLambdaEnd) {
                doPolarization = true;
                // Scale the polarization values to account for the squished path.
                polarizationScale = lPowPol;    // TODO calculate @ inner loop
            } else {
                doPolarization = true;
                polarizationScale = 1.0;
            }            
        }
        doPermanentRealSpace = true;
        permanentScale = lPowPerm;
        dEdLSign = 1.0;

        double energy = computeEnergy(false);

        return energy;
    }

    /**
     * 2.) Condensed phase system without the ligand.
     *      A.) No permanent real space electrostatics needs to be calculated because
     *          this was handled analytically in step 1.
     *      B.) Permanent reciprocal space scaled by (1 - lambda).
     *      C.) Polarization scaled by (1 - lambda).
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
     *      A.) Real space with an Ewald coefficient of 0.0 (no reciprocal space).
     *      B.) Polarization scaled as in Step 2 by (1 - lambda).
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
            if (reciprocalSpaceTerm && aewald > 0.0) {
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
            if (reciprocalSpaceTerm && aewald > 0.0) {
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
            if (reciprocalSpaceTerm && aewald > 0.0) {
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
            if (reciprocalSpaceTerm && aewald > 0.0) {
                parallelTeam.execute(reciprocalEnergyRegion);
                interactions += nAtoms;
                eself = reciprocalEnergyRegion.getPermanentSelfEnergy();
                erecip = reciprocalEnergyRegion.getPermanentReciprocalEnergy();
                eselfi = reciprocalEnergyRegion.getInducedDipoleSelfEnergy();
                erecipi = reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
            }

            realSpaceEnergyTotalQI = -System.nanoTime();
            parallelTeam.execute(realSpaceEnergyRegionQI);
            realSpaceEnergyTotalQI += System.nanoTime();
            ereal = realSpaceEnergyRegionQI.getPermanentEnergy();
            ereali = realSpaceEnergyRegionQI.getPolarizationEnergy();
            interactions += realSpaceEnergyRegionQI.getInteractions();

            if (DEBUG() && (lambdaMode == LambdaMode.OFF || lambdaMode == LambdaMode.CONDENSED)) {
                logger.info(format(" (perm,pol,time):  qi (%12.6f  %12.6f) %8.3f ms",
                        ereal, ereali, realSpaceEnergyTotalQI * TO_MS));
            }
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


    @Override
    public int getInteractions() {
        return interactions;
    }

    @Override
    public double getPermanentEnergy() {
        return permanentMultipoleEnergy;
    }

    @Override
    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    @Override
    public double getGKEnergy() {
        return generalizedKirkwoodEnergy;
    }

    @Override
    public double getCavitationEnergy(boolean throwError) {
        return generalizedKirkwood.getCavitationEnergy(throwError);
    }

    @Override
    public double getDispersionEnergy(boolean throwError) {
        return generalizedKirkwood.getDispersionEnergy(throwError);
    }

    public double getCavitationEnergy() {
        return generalizedKirkwood.getCavitationEnergy(false);
    }

    public double getDispersionEnergy() {
        return generalizedKirkwood.getDispersionEnergy(false);
    }

    @Override
    public int getGKInteractions() {
        return gkInteractions;
    }

    public void getGradients(double grad[][]) {
        if (grad == null) {
            grad = new double[3][nAtoms];
        }
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

    @Override
    protected double[][][] getGradient() {
        return grad;
    }

    @Override
    protected double[][][] getTorque() {
        return torque;
    }

    @Override
    protected double[][][] getLambdaGradient() {
        return lambdaGrad;
    }

    @Override
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

        if (System.getProperty("printInducedDipoles") != null) {
            StringBuilder sb = new StringBuilder();
            sb.append("     Atom                                         Induced Dipole \n");
            sb.append("    ======                                       ================\n");
            for (int i = 0; i < nAtoms; i++) {
                sb.append(format("%-47s: (%+8.6f %+8.6f %+8.6f)\n", atoms[i],
                        inducedDipole[0][i][0], inducedDipole[0][i][1], inducedDipole[0][i][2]));
            }
            logger.info(sb.toString());
        }
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
                if (reciprocalSpaceTerm && aewald > 0.0) {
                    reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
                }
                sectionTeam.execute(inducedDipoleFieldRegion);
                if (reciprocalSpaceTerm && aewald > 0.0) {
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
                        " %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * TO_SECONDS));
            }
            /**
             * If the RMS Debye change increases, fail the SCF process.
             */
            if (eps > previousEps) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Fatal SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
                throw new EnergyException(message, false);
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
                throw new EnergyException(message, false);
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
                    TO_SECONDS * directTime));
            startTime = System.nanoTime() - startTime;
            sb.append(format(" Total:                   %7.4f",
                    startTime * TO_SECONDS));
            logger.info(sb.toString());
        }
        return completedSCFCycles;
    }

    @Override
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
                    realSpacePermTotalQI -= System.nanoTime();
                    parallelTeam.execute(permanentRealSpaceFieldRegion);
                    realSpacePermTotalQI += System.nanoTime();
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
                if (reciprocalSpaceTerm && aewald > 0.0) {
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
                            if (!use[k]) {
                                continue;
                            }
                            boolean sameMolecule = (moleculei == molecule[k]);
                            if (lambdaMode == LambdaMode.VAPOR) {
                                if ((intermolecularSoftcore && !sameMolecule)
                                        || (intramolecularSoftcore && sameMolecule)) {
                                    continue;
                                }
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
                                final double pgamma = min(pti, ptk);
                                final double rdamp = r * damp;
                                damp = -pgamma * rdamp * rdamp * rdamp;
                                if (damp > -50.0) {
                                    double expdamp = exp(damp);
                                    scale3 = 1.0 - expdamp;
                                    scale5 = 1.0 - expdamp * (1.0 - damp);
                                    scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                }
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
                if (reciprocalSpaceTerm && aewald > 0.0) {
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
                    realSpaceSCFTotalQI -= System.nanoTime();
                    pt.execute(polarizationRealSpaceFieldRegion);
                    realSpaceSCFTotalQI += System.nanoTime();
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
                            if (!use[k]) {
                                continue;
                            }
                            boolean sameMolecule = (moleculei == molecule[k]);
                            if (lambdaMode == LambdaMode.VAPOR) {
                                if ((intermolecularSoftcore && !sameMolecule)
                                        || (intramolecularSoftcore && sameMolecule)) {
                                    continue;
                                }
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
    private class RealSpaceEnergyRegionQI extends ParallelRegion {

        private double permanentEnergy;
        private double polarizationEnergy;
        private double mutualScale = (polarization == Polarization.DIRECT || polarization == Polarization.NONE)
                ? 0.0 : 1.0;
        private final int numThreads;

        private final SharedInteger sharedInteractions;
        private final RealSpaceEnergyLoopQI realSpaceEnergyLoops[];
        private final MultipoleTensor[] tensors;

        public RealSpaceEnergyRegionQI(int nt) {
            numThreads = nt;
            sharedInteractions = new SharedInteger();
            realSpaceEnergyLoops = new RealSpaceEnergyLoopQI[nt];
            tensors = new MultipoleTensor[nt];
            for (int thread = 0; thread < nt; thread++) {
                realSpaceEnergyLoops[thread] = new RealSpaceEnergyLoopQI();
                tensors[thread] = new MultipoleTensor(
                        OPERATOR.SCREENED_COULOMB, COORDINATES.QI, 5, aewald);
            }
            esvRealSpaceDeriv = new SharedDouble[numESVs];
            for (int i = 0; i < numESVs; i++) {
                esvRealSpaceDeriv[i] = new SharedDouble(0.0);
            }
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
            // [threadID][X/Y/Z][atomID]
            for (int thread = 0; thread < maxThreads; thread++) {
                for (int i = 0; i < 3; i++) {
                    fill(grad[thread][i], 0.0);
                    fill(torque[thread][i], 0.0);
                    fill(field[thread][i], 0.0);
                    fill(fieldCR[thread][i], 0.0);
                    if (lambdaTerm) {
                        fill(lambdaGrad[thread][i], 0.0);
                        fill(lambdaTorque[thread][i], 0.0);
                    }
                }
            }
            for (int i = 0; i < numESVs; i++) {
                esvRealSpaceDeriv[i] = new SharedDouble(0.0);
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            realSpaceEnergyLoops[threadIndex] = new RealSpaceEnergyLoopQI();
            try {
                execute(0, nAtoms - 1, realSpaceEnergyLoops[threadIndex]);
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
                double e = realSpaceEnergyLoops[i].permanentEnergy;
                if (Double.isNaN(e)) {
                    //logger.severe(String.format(" The permanent multipole energy of thread %d is %16.8f", i, e));
                    throw new EnergyException(String.format(" The permanent multipole energy of thread %d is %16.8f", i, e), true);
                }
                permanentEnergy += e;
                double ei = realSpaceEnergyLoops[i].inducedEnergy;
                if (Double.isNaN(ei)) {
                    //logger.severe(String.format(" The polarization energy of thread %d is %16.8f", i, ei));
                    throw new EnergyException(String.format(" The polarization energy of thread %d is %16.8f", i, ei), true);
                }
                polarizationEnergy += ei;
            }
            permanentEnergy *= ELECTRIC;
            polarizationEnergy *= ELECTRIC;
            
            /*
            if (DEBUG() && lambdaTerm && lambdaMode == LambdaMode.CONDENSED) {
                double[][][] compsReduced = new double[nAtoms][nAtoms][nComps];
                double[] termSum = new double[]{0.0, 0.0};
                double[] termSum2 = new double[]{0.0, 0.0, 0.0, 0.0, 0.0};
                double termSum2Sum = 0.0;
                for (int i = 0; i < nAtoms; i++) {
                    for (int k = 0; k < nAtoms; k++) {
                        for (int comp = 0; comp < nComps; comp++) {
                            compsReduced[i][k][comp] = 0.0;
                            for (int thread = 0; thread < maxThreads; thread++) {
                                compsReduced[i][k][comp] += compQI[thread][i][k][comp];
                            }
                        }
                        termSum[0] += compsReduced[i][k][1];
                        termSum[1] += compsReduced[i][k][2];
                        for (int c2 = 0; c2 < 5; c2++) {
                            termSum2[c2] += comp2QI[i][k][c2].get();
                            termSum2Sum += comp2QI[i][k][c2].get();
                        }
                        double[] vals = new double[]{termSum[0], termSum[1], compsReduced[i][k][1], compsReduced[i][k][2]};
                        if (DEBUG() > 1) {
                            logger.info(format("    Creating termSums: i,j,sums,comps: %d,%d,%s", i, k, formatArray(vals)));
                        }
                    }
                }
                if (DEBUG() > 0 || (lambdaMode == LambdaMode.CONDENSED || lambdaMode == LambdaMode.OFF)) {
                    for (int i = 0; i < nAtoms; i++) {
                        for (int k = 0; k < nAtoms; k++) {
                            for (int comp = 0; comp < nComps; comp++) {
                                if (compsReduced[i][k][comp] != compQIshared[i][k][comp].get()) {
                                    logger.info(format("COMP MISMATCH (%d,%d,%d): %g, %g",
                                            i, k, comp,
                                            compsReduced[i][k][comp], compQIshared[i][k][comp].get()));
                                }
                            }
                        }
                    }
                }
                logf("QiS TOTALS: dE/dL = %g + %g = %g --> %g + %g = %g    (shdEdLqi %g)",
                        termSum[0], termSum[1], termSum[0] + termSum[1],
                        termSum[0] * ELECTRIC, termSum[1] * ELECTRIC, termSum[0] * ELECTRIC + termSum[1] * ELECTRIC,
                        shareddEdLambdaQI.get()));
                logf("Qi2 TOTALS: d2E/dL2 = %g + %g + %g + %g + %g = %g",
                        termSum2[0], termSum2[1], termSum2[2], termSum2[3], termSum2[4], termSum2Sum));
            }
            */
        }

        /**
         * The Real Space Gradient Loop class contains methods and thread local
         * variables to parallelize the evaluation of the real space permanent
         * and polarization energies and gradients.
         */
        private class RealSpaceEnergyLoopQI extends IntegerForLoop {

            private double r2O, r2B;
            private double scale, scalep, scaled;
            private double lBufferDistance, l2;
            private boolean soft;
            private double selfScale;
            private double permanentEnergy;
            private double inducedEnergy;
            private double dUdL, d2UdL2;
            private int i, k, iSymm, count;
            private SymOp symOp;
            // Store contributions to the gradient.
            private double gX[], gY[], gZ[], tX[], tY[], tZ[];
            private double gxk_local[], gyk_local[], gzk_local[];
            private double txk_local[], tyk_local[], tzk_local[];
            // Store contributions to dE/dX/dL
            private double lgX[], lgY[], lgZ[], ltX[], ltY[], ltZ[];
            private double lxk_local[], lyk_local[], lzk_local[];
            private double ltxk_local[], ltyk_local[], ltzk_local[];
            // Store contributions to dE/dX/dLDH
            private double ldhgX[][], ldhgY[][], ldhgZ[][], ldhtX[][], ldhtY[][], ldhtZ[][];
            // Masking rules
            private double masking_local[];
            private double maskingp_local[];
            private double maskingd_local[];
            private final double dx_local[];
            private final double rot_local[][];
            private final double work[][];

            // Force and torque contributions for a single interaction.
            private MultipoleTensor tensor;
            private final double[] Fi, Ti, Tk;
            private final double[] permFi, permTi, permTk;
            private final double[] polFi, polTi, polTk;
            private final double[] FiC, TiC, TkC;
            private final double[] FiT, TiT, TkT;
            private final double[] energy;

            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public RealSpaceEnergyLoopQI() {
                super();
                dx_local = new double[3];
                work = new double[15][3];
                rot_local = new double[3][3];
                Fi = new double[3];
                Ti = new double[3];
                Tk = new double[3];
                permFi = new double[3];
                permTi = new double[3];
                permTk = new double[3];
                polFi = new double[3];
                polTi = new double[3];
                polTk = new double[3];
                FiC = new double[3];
                TiC = new double[3];
                TkC = new double[3];
                FiT = new double[3];
                TiT = new double[3];
                TkT = new double[3];
                energy = new double[2];
            }

            private void init() {
                if (masking_local == null || masking_local.length < nAtoms) {
                    txk_local = new double[nAtoms];
                    tyk_local = new double[nAtoms];
                    tzk_local = new double[nAtoms];
                    gxk_local = new double[nAtoms];
                    gyk_local = new double[nAtoms];
                    gzk_local = new double[nAtoms];
                    if (lambdaTerm) {
                        lxk_local = new double[nAtoms];
                        lyk_local = new double[nAtoms];
                        lzk_local = new double[nAtoms];
                        ltxk_local = new double[nAtoms];
                        ltyk_local = new double[nAtoms];
                        ltzk_local = new double[nAtoms];
                    }

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
                tensor = tensors[threadIndex];
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
                    // Do all the work.
                    realSpaceChunk(lb, ub);
                    // Collect results.
                    if (gradient) {
                        // Turn symmetry mate torques into gradients
                        if (rotateMultipoles) {
                            torque(iSymm, txk_local, tyk_local, tzk_local,
                                    gxk_local, gyk_local, gzk_local,
                                    work[0], work[1], work[2], work[3], work[4],
                                    work[5], work[6], work[7], work[8], work[9],
                                    work[10], work[11], work[12], work[13], work[14]);
                        }
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
                        if (rotateMultipoles) {
                            // Turn symmetry mate torques into gradients
                            torque(iSymm, ltxk_local, ltyk_local, ltzk_local,
                                    lxk_local, lyk_local, lzk_local,
                                    work[0], work[1], work[2], work[3], work[4],
                                    work[5], work[6], work[7], work[8], work[9],
                                    work[10], work[11], work[12], work[13], work[14]);
                        }
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
                if (esvTerm) {
                    for (int i = 0; i < numESVs; i++) {
                        // REM do this in the inner loop since dUdEsv obsoleted
//                        esvRealSpaceDeriv[i].addAndGet(dUdEsvLocal[i] * ELECTRIC);
                        esvRealSpaceDeriv[i].addAndGet(0.0 * ELECTRIC); // intermediate dUdL * ELECTRIC
                    }
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
                final double mpole[][] = globalMultipole[0];    // [nSymm][nAtoms][10]
                final double ind[][] = inducedDipole[0];        // [nsymm][nAtoms][3]
                final double indCR[][] = inducedDipoleCR[0];
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
                    /**
                     * Set masking scale factors.
                     */
                    if (iSymm == 0) {
                        applyScaleFactors(ai);
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double globalMultipolei[] = mpole[i];
                    final double inducedDipolei[] = ind[i];
                    final double inducedDipolepi[] = indCR[i];
                    final boolean softi = isSoft[i];        // includes ESV softs
                    final boolean esvi = esvAtoms[i];
                    final double pdi = ipdamp[i];
                    final double pti = thole[i];
                    final int list[] = lists[i];
                    final int npair = realSpaceCounts[iSymm][i];
                    for (int j = 0; j < npair; j++) {
                        k = list[j];
                        if (!use[k]) {
                            continue;
                        }
                        final boolean softk = isSoft[k];    // includes ESV softs
                        final boolean esvk = esvAtoms[k];
                        boolean sameMolecule = (moleculei == molecule[k]);
                        if (lambdaMode == LambdaMode.VAPOR) {
                            if ((intermolecularSoftcore && !sameMolecule)
                                    || (intramolecularSoftcore && sameMolecule)) {
                                continue;
                            }
                        }
                        selfScale = 1.0;
                        if (i == k) {
                            selfScale = 0.5;
                        }
                        lBufferDistance = 0.0;
                        l2 = 1.0;
                        soft = (softi || softk);
                        if (soft && doPermanentRealSpace) {
                            lBufferDistance = lAlpha;
                            l2 = permanentScale;
                        }
                        if (esvTerm && (esvi || esvk)) {
                            double esvLambdaProduct = esvLambda[i] * esvLambda[k] * lambda;
//                            initSoftCore = true;    // only needed on system expansion (and destruction)
                            setLambda(esvLambdaProduct);
                            
                            /*  EXAMPLE for system calls; double-check eqs.
                                final int idxi = esvSystem.atomEsvId(i);
                                final int idxk = esvSystem.atomEsvId(k);
                                if (esvi) {
                                    final double dlpdli = esvLambda[k] * lambda;
                                    final double dEsvPartI = dedlp * dlpdli;
                                    esvRealSpaceDeriv[idxi].addAndGet(dEsvPartI);
                                }
                            */                            
                            
                            lBufferDistance = lAlpha;
                            l2 = permanentScale;
                        }
                        final double xk = neighborX[k];
                        final double yk = neighborY[k];
                        final double zk = neighborZ[k];
                        dx_local[0] = xk - xi;
                        dx_local[1] = yk - yi;
                        dx_local[2] = zk - zi;
                        r2B = crystal.image(dx_local);

                        final double globalMultipolek[] = neighborMultipole[k];
                        final double inducedDipolek[] = neighborInducedDipole[k];
                        final double inducedDipolepk[] = neighborInducedDipolep[k];
                        final double pdk = ipdamp[k];   // == 1/polarizability^6
                        final double ptk = thole[k];
                        scale = masking_local[k];
                        scalep = maskingp_local[k];
                        scaled = maskingd_local[k];
                        double damp = min(pti, ptk);
                        double aiak = pdi * pdk;
                        if (doPermanentRealSpace) {     // TODO prefer pairPermPol once available
                            permanentEnergy += pairPerm(dx_local, globalMultipolei, globalMultipolek);
                            count++;
                        } else if (doPermanentRealSpace && doPolarization && polarization != Polarization.NONE) {
                            logf("Skipping unfinished QI polarization loop.");
//                            double eTotal = pairPermPol(dx_local, globalMultipolei, globalMultipolek,
//                                    inducedDipolei, inducedDipolek, inducedDipolepi, inducedDipolepk,
//                                    damp, aiak, energy);
//                            permanentEnergy += energy[0];
//                            inducedEnergy += energy[1];
//                            count++;
                        } else {
                            logf("Skipping unfinished QI induction loop.");
//                            inducedEnergy += pairPol(dx_local, globalMultipolei, globalMultipolek,
//                                    inducedDipolei, inducedDipolek, inducedDipolepi, inducedDipolepk,
//                                    damp, aiak);
                        }
                    }
                    /**
                     * Reset masking scale factors.
                     */
                    if (iSymm == 0) {
                        resetScaleFactors(ai);
                    }
                }
            }

            /*
            private double pairPermPol(double[] r, double[] Qi, double[] Qk,
                    double[] ui, double[] uk, double[] uiCR, double[] ukCR,
                    double damp, double aiak, double energy[]) {
                /**
                 * Compute screened real space interactions.
                tensor.setR(r, lBufferDistance);
                // Add buffer.

                tensor.setMultipolesQI(Qi, Qk);
                tensor.setDipolesQI(ui, uiCR, uk, ukCR);
                tensor.setOperator(OPERATOR.SCREENED_COULOMB);

                // Order 6 QI
                tensor.order5QI();

                double ePermScreened = tensor.multipoleEnergyQI(permFi, permTi, permTk);
                double ePolScreened = tensor.polarizationEnergyQI(1.0, 1.0, mutualScale, polFi, polTi, polTk);

                /**
                 * Subtract away masked Coulomb interactions included in PME.
                double scale1 = 1.0 - scale;
                double scaled1 = 1.0 - scaled;
                double scalep1 = 1.0 - scalep;
                double ePermCoulomb = 0.0;
                double ePolCoulomb = 0.0;
                if (scale1 != 0.0 || scaled1 != 0.0 || scalep1 != 0.0) {
                    tensor.setOperator(OPERATOR.COULOMB);
                    tensor.order5QI();
                    if (scale1 != 0.0) {
                        ePermCoulomb = scale1 * tensor.multipoleEnergyQI(FiC, TiC, TkC);
                        permFi[0] -= scale1 * FiC[0];
                        permFi[1] -= scale1 * FiC[1];
                        permFi[2] -= scale1 * FiC[2];
                        permTi[0] -= scale1 * TiC[0];
                        permTi[1] -= scale1 * TiC[1];
                        permTi[2] -= scale1 * TiC[2];
                        permTk[0] -= scale1 * TkC[0];
                        permTk[1] -= scale1 * TkC[1];
                        permTk[2] -= scale1 * TkC[2];
                    }
                    if (scaled1 != 0.0 || scalep1 != 0.0) {
                        ePolCoulomb += tensor.polarizationEnergyQI(
                                scaled1, scalep1, 0.0, FiC, TiC, TkC);
                        polFi[0] -= FiC[0];
                        polFi[1] -= FiC[1];
                        polFi[2] -= FiC[2];
                        polTi[0] -= TiC[0];
                        polTi[1] -= TiC[1];
                        polTi[2] -= TiC[2];
                        polTk[0] -= TkC[0];
                        polTk[1] -= TkC[1];
                        polTk[2] -= TkC[2];
                    }
                }

                /**
                 * Account for Thole Damping.
                double eThole = 0.0;
                tensor.setTholeDamping(damp, aiak);
                boolean applyThole = tensor.applyDamping();
                if (applyThole) {
                    tensor.setOperator(OPERATOR.THOLE_FIELD);
                    tensor.order4QI();
                    tensor.setDipolesQI(ui, uiCR, uk, ukCR);
                    eThole = tensor.polarizationEnergyQI(scaled, scalep, mutualScale, FiT, TiT, TkT);
                    polFi[0] -= FiT[0];
                    polFi[1] -= FiT[1];
                    polFi[2] -= FiT[2];
                    polTi[0] -= TiT[0];
                    polTi[1] -= TiT[1];
                    polTi[2] -= TiT[2];
                    polTk[0] -= TkT[0];
                    polTk[1] -= TkT[1];
                    polTk[2] -= TkT[2];
                }

                final double ePerm = selfScale * l2 * (ePermScreened - ePermCoulomb);
                energy[0] = ePerm;

                if (gradient) {
                    double prefactor = ELECTRIC * selfScale * l2;
                    gX[i] += prefactor * permFi[0];
                    gY[i] += prefactor * permFi[1];
                    gZ[i] += prefactor * permFi[2];
                    tX[i] += prefactor * permTi[0];
                    tY[i] += prefactor * permTi[1];
                    tZ[i] += prefactor * permTi[2];
                    gxk_local[k] -= prefactor * permFi[0];
                    gyk_local[k] -= prefactor * permFi[1];
                    gzk_local[k] -= prefactor * permFi[2];
                    txk_local[k] += prefactor * permTk[0];
                    tyk_local[k] += prefactor * permTk[1];
                    tzk_local[k] += prefactor * permTk[2];
                    /**
                     * This is dU/dL/dX for the first term of dU/dL: d[dlPow *
                     * ereal]/dx
                    if (lambdaTerm && soft) {
                        prefactor = ELECTRIC * selfScale * dEdLSign * dlPowPerm;
                        lgX[i] += prefactor * permFi[0];
                        lgY[i] += prefactor * permFi[1];
                        lgZ[i] += prefactor * permFi[2];
                        ltX[i] += prefactor * permTi[0];
                        ltY[i] += prefactor * permTi[1];
                        ltZ[i] += prefactor * permTi[2];
                        lxk_local[k] -= prefactor * permFi[0];
                        lyk_local[k] -= prefactor * permFi[1];
                        lzk_local[k] -= prefactor * permFi[2];
                        ltxk_local[k] += prefactor * permTk[0];
                        ltyk_local[k] += prefactor * permTk[1];
                        ltzk_local[k] += prefactor * permTk[2];
                    }
                }

                final double e = selfScale * 0.5 * (ePolScreened - ePolCoulomb - eThole);
                if (!(gradient || lambdaTerm || esvTerm)) {
                    double ePol = polarizationScale * e;
                    energy[1] = ePol;
                    return ePerm + ePol;
                }

                double scalar = ELECTRIC * polarizationScale * selfScale;
                gX[i] += scalar * polFi[0];
                gY[i] += scalar * polFi[1];
                gZ[i] += scalar * polFi[2];
                tX[i] += scalar * polTi[0];
                tY[i] += scalar * polTi[1];
                tZ[i] += scalar * polTi[2];
                gxk_local[k] -= scalar * polFi[0];
                gyk_local[k] -= scalar * polFi[1];
                gzk_local[k] -= scalar * polFi[2];
                txk_local[k] += scalar * polTk[0];
                tyk_local[k] += scalar * polTk[1];
                tzk_local[k] += scalar * polTk[2];
                if (lambdaTerm) {
                    dUdL += dEdLSign * dlPowPol * e;
                    d2UdL2 += dEdLSign * d2lPowPol * e;
                    scalar = ELECTRIC * dEdLSign * dlPowPol * selfScale;
                    lgX[i] += scalar * polFi[0];
                    lgY[i] += scalar * polFi[1];
                    lgZ[i] += scalar * polFi[2];
                    ltX[i] += scalar * polTi[0];
                    ltY[i] += scalar * polTi[1];
                    ltZ[i] += scalar * polTi[2];
                    lxk_local[k] -= scalar * polFi[0];
                    lyk_local[k] -= scalar * polFi[1];
                    lzk_local[k] -= scalar * polFi[2];
                    ltxk_local[k] += scalar * polTk[0];
                    ltyk_local[k] += scalar * polTk[1];
                    ltzk_local[k] += scalar * polTk[2];
                }

                double ePol = polarizationScale * e;
                energy[1] = ePol;

                return ePerm + ePol;
            }
             */
            private double pairPerm(double[] r, double[] Qio, double[] Qko) {
                double[] Qi = new double[10], Qk = new double[10];
                for (int i = 0; i < 10; i++) {
                    Qi[i] = 0.0;
                    Qk[i] = 0.0;
                }
                for (int i = 0; i < 10; i++) {
                    if ((i == 0 && useCharges)
                            || (i >= 1 && i < 4 && useDipoles)
                            || (i >= 4 && useQuadrupoles)) {
                        Qi[i] = Qio[i];
                        Qk[i] = Qko[i];
                    }
                }

                /* REM
                assert selfScale == 1.0
                assert soft     -> l2 == permanentScale
                assert !soft    -> l2 == 1.0
                */

                /**
                 * Set MultipoleTensor distance; handle lambda buffering.
                 */
                r2O = crystal.image(dx_local);
                if (soft && (lambdaTerm || esvTerm)) {
                    
                    tensor.setR(dx_local, lBufferDistance);
                } else {    // hard
                    logf("No softcore buffering > i,k: %d %d", i, k);
                    tensor.setR(dx_local);
                }
                r2B = tensor.getR()[4];
                final double rO = sqrt(r2O);
                final double rB = sqrt(r2B);

                /**
                 * Compute screened real space interactions.
                 */
                double ePerm, dPermdL;
                double scale1 = 1.0 - scale;
                if (aewald == 0.0 || scale == 1.0) {
                    tensor.setOperator(OPERATOR.COULOMB);
                } else {
                    tensor.setOperator(OPERATOR.SCREENED_COULOMB);
                }
                tensor.setR(dx_local, lBufferDistance);
                tensor.setMultipolesQI(Qi, Qk);
                tensor.order6QI();
                ePerm = tensor.multipoleEnergyQI(permFi, permTi, permTk);
                dPermdL = tensor.getdEdF();
                logf(" scrn > ePerm,dPerm,d2Perm:   %g %g", ePerm, dPermdL);

                if (scale != 1.0) {
                    tensor.setOperator(OPERATOR.COULOMB);
                    tensor.setR(dx_local, lBufferDistance);
                    tensor.order6QI();
                    ePerm -= scale1 * tensor.multipoleEnergyQI(FiC, TiC, TkC);
                    dPermdL -= scale1 * tensor.getdEdF();
                    logf(" -coul=res > ePerm,dPerm,d2Perm: %g %g", ePerm, dPermdL);

                    permFi[0] -= scale1 * FiC[0];
                    permFi[1] -= scale1 * FiC[1];
                    permFi[2] -= scale1 * FiC[2];
                    permTi[0] -= scale1 * TiC[0];
                    permTi[1] -= scale1 * TiC[1];
                    permTi[2] -= scale1 * TiC[2];
                    permTk[0] -= scale1 * TkC[0];
                    permTk[1] -= scale1 * TkC[1];
                    permTk[2] -= scale1 * TkC[2];
                }

                if (DEBUG()) {
                    if (selfScale != 1.0 || (permanentScale != 1.0 && permanentScale != lambda)) {
                        logger.severe(format("Non-unity selfScale: %g %g", selfScale, permanentScale));
                    }
                    if (lPowPerm != lambda || permanentScale != lPowPerm
                            || (!soft && l2 != 1.0) || (soft && l2 != permanentScale)) {
                        logger.severe(format("Inconsistency > l2,lPowPerm,lambda: %g %g %g %g", lambda, l2, lPowPerm, permanentScale));
                    }
                }

                final double e = selfScale * l2 * ePerm;

                if (!gradient && !lambdaTerm && !soft) {
                    return e;
                }

                double scalar;
                if (gradient) {
                    scalar = ELECTRIC * selfScale * l2;
                    gX[i] += scalar * permFi[0];
                    gY[i] += scalar * permFi[1];
                    gZ[i] += scalar * permFi[2];
                    tX[i] += scalar * permTi[0];
                    tY[i] += scalar * permTi[1];
                    tZ[i] += scalar * permTi[2];
                    gxk_local[k] -= scalar * permFi[0];
                    gyk_local[k] -= scalar * permFi[1];
                    gzk_local[k] -= scalar * permFi[2];
                    txk_local[k] += scalar * permTk[0];
                    tyk_local[k] += scalar * permTk[1];
                    tzk_local[k] += scalar * permTk[2];
                }

                if (lambdaTerm && soft) {
                    /**
                     * This is dU/dL/dX for the first term of dU/dL: d[dlPow *
                     * ereal]/dx
                     *
                     * But... MT returns as either d?/d[sqrt(r^2+a(1-L))] <--
                     * bufferCoords.QI or as d?/d[z+a(1-L)] <--
                     * bufferCoords.GLOBAL
                     */
                    scalar = ELECTRIC * selfScale * dEdLSign * dlPowPerm;
                    lgX[i] += scalar * permFi[0];
                    lgY[i] += scalar * permFi[1];
                    lgZ[i] += scalar * permFi[2];
                    ltX[i] += scalar * permTi[0];
                    ltY[i] += scalar * permTi[1];
                    ltZ[i] += scalar * permTi[2];
                    lxk_local[k] -= scalar * permFi[0];
                    lyk_local[k] -= scalar * permFi[1];
                    lzk_local[k] -= scalar * permFi[2];
                    ltxk_local[k] += scalar * permTk[0];
                    ltyk_local[k] += scalar * permTk[1];
                    ltzk_local[k] += scalar * permTk[2];

                    double S = System.getProperty("dedlSign0") != null
                            ? dEdLSign * lPowPerm : lPowPerm;       // 1.0
                    double dSdL = System.getProperty("dedlSign1") != null
                            ? dEdLSign * dlPowPerm : dlPowPerm;     // 1.0
                    double d2SdL2 = System.getProperty("dedlSign2") != null
                            ? dEdLSign * d2lPowPerm : d2lPowPerm;   // 0.0

                    /* [FMODE]
                     * Old Factoring Method f = sqrt(r^2 + lAlpha) df/dL =
                     * -alpha * (1.0 - lambda) / f g = 1 / sqrt(r^2 + lAlpha)
                     * dg/dL = alpha * (1.0 - lambda) / (r^2 + lAlpha)^(3/2)
                     * define dlAlpha = alpha * 1.0 - lambda) then df/dL =
                     * -dlAlpha / f and dg/dL = dlAlpha * g^3
                     * 
                     * These two working option sets reflect that first derivatives can
                     * be taken from either the Global or QI frame regardless of
                     * which was used for multipole interaction, provided the
                     * appropriate Jacobian (== trivial, it's dL
                    
                        case 3:     // works for 1st deriv; requires dlAlphaMode == FACTORED, mt-Rmode == INDEPENDENT
                            double totalDist = sqrt(crystal.image(dx_local) + lBufferDistance);
                            F = lAlpha;
                            dFdL = dlAlpha / rB;
                            // second deriv solution to dedz*Q == dedl:
                            // Q -> (\[Alpha] (B + R^2 - 3 \[Alpha] (1 - \[Lambda])^2))/(B - 2 R^2)
//                            d2FdL2 = a*(r2orig - 2*B) / (B - 2*r2orig);
                            d2FdL2 = a * (r2O - 2 * B) / (B - 2 * r2O);
                            logf(" Fmode%d > F:       %g\n"
                                    + "   dlAlpha,rB,dFdL: %g / %g (%g) = %g\n"
                                    + "   r2O,B,d2FdL2:    %g ... %g = %g",
                                    Fmode, F, dlAlpha, rB, totalDist, dFdL,
                                    r2O, B, d2FdL2));
                            break;
                        case 4:     // works for 1st deriv; requires dlAlphaMode == FACTORED, mt-Rmode == INDEPENDENT
                            F = lAlpha;
                            dFdL = -dlAlpha / rO;
                            // second deriv solution to dedz*Q == dedl:
                            // Q -> (\[Alpha] (B + R^2 - 3 \[Alpha] (1 - \[Lambda])^2))/(B - 2 R^2)
//                            d2FdL2 = a*(r2orig - 2*B) / (B - 2*r2orig);
                            d2FdL2 = a * (r2O - 2 * B) / (B - 2 * r2O);
                            logf(" Fmode%d > F: %g\n"
                                    + "   dlAlpha,rO,dFdL: -%g / %g = %g\n"
                                    + "   r2O,B,d2FdL2:     %g ... %g = %g",
                                    Fmode, F, dlAlpha, rO, dFdL,
                                    r2O, B, d2FdL2));
                            break;
                    */  // [/FMODE]

                    final double a = permLambdaAlpha, B = lBufferDistance;
                    final double F = lAlpha;
                    final double dFdL = dlAlpha / rB;
//                  (Fmode) -> final double dFdL = -dlAlpha / rO;
                    final double P = ePerm;
                    final double dPdF = dPermdL;

                    // d(SP)dL = (dSdL*P)+(S*dPdL) = (dSdL*P)+(S*dPdZ*dZdL)
                    final double termA = dEdLSign * (dSdL * P);
                    final double termB = (S * dPdF * dFdL);
                    final double thisInteraction = selfScale * (termA + termB);
                    final double[] components = new double[]{thisInteraction, selfScale * termA, selfScale * termB,
                        selfScale * (termA + termB), dSdL, P, S, dPdF, dFdL};
                    logf("i,k;termA,termB,sum:  %8d   %8d     %8g   %8g   %8g\n"
                            + "             comps:  (%8g * %8g) + (%8g * %8g * %8g) = %g\n",
                            i, k, termA, termB, termA + termB,
                            dSdL, P, S, dPdF, dFdL, thisInteraction);
                    if (!Double.isFinite(thisInteraction)) {
                        logger.warning(format("NaN output from PME.\n"
                                + "    comps thread %d: %s", getThreadIndex(), formatArray(components)));
                    }
                    dUdL += thisInteraction;

                    /* REFERENCE: original derivation
                        double S = dEdLSign * dlPowPerm, dSdL = dEdLSign, d2SdL2 = 0.0;
                        double P = ePerm, dPdL = dPermdL, d2PdL2 = d2PermdL2;
                        double F = lAlpha, dFdL = dlAlpha, d2FdL2 = d2lAlpha;
                        double dPdF = (dFdL != 0.0) ? dPdL / dFdL : 0.0;
                        double d2PdF2 = (d2FdL2 != 0.0) ? d2PdL2 / d2FdL2 : 0.0;
                        dUdL += selfScale * ((dSdL * P) + (S * dPdF * dFdL));
                        d2UdL2 += selfScale * ((d2SdL2 * P) + (dSdL * S * dPdF * dFdL)
                                + ((dSdL * dPdF) + (S * d2PdF2)) * dFdL
                                + (S * dPdF * d2FdL2)); */
                    /* REFERENCE: from ParticleMeshEwaldCart
                        final double e = selfScale * l2 * (ereal - efix);
                        dUdL += selfScale * (dEdLSign * dlPowPerm * ereal + l2 * dlAlpha * dRealdL);
                        d2UdL2 += selfScale * (dEdLSign * (d2lPowPerm * ereal
                                + dlPowPerm * dlAlpha * dRealdL
                                + dlPowPerm * dlAlpha * dRealdL)
                                + l2 * d2lAlpha * dRealdL
                                + l2 * dlAlpha * dlAlpha * d2RealdL2);  */
                    
                    /**
                     * Add in dU/dL/dX for the second term of dU/dL:
                     * d[lPow*dlAlpha*dRealdL]/dX
                     */
                    // No additional call to MT; use 6th order tensor instead.
                    scalar = ELECTRIC * selfScale * l2 * dlAlpha;
//                    if (bufferCoords == COORDINATES.QI) {
//                    switch dxlMode:
//                      1: scalar /= (rO * rB * rB * rB);
//                      2: scalar /= (rB * rB * rB);
//                      3: scalar /= rB;
//                    }
                    lgX[i] += scalar * permFi[0];
                    lgY[i] += scalar * permFi[1];
                    lgZ[i] += scalar * permFi[2];
                    ltX[i] += scalar * permTi[0];
                    ltY[i] += scalar * permTi[1];
                    ltZ[i] += scalar * permTi[2];
                    lxk_local[k] -= scalar * permFi[0];
                    lyk_local[k] -= scalar * permFi[1];
                    lzk_local[k] -= scalar * permFi[2];
                    ltxk_local[k] += scalar * permTk[0];
                    ltyk_local[k] += scalar * permTk[1];
                    ltzk_local[k] += scalar * permTk[2];
                }

                return e;
            }

            private double pairPerm_globalMT(double[] r, double[] Qi, double[] Qk) {
                MultipoleTensor tensor = new MultipoleTensor(
                        OPERATOR.SCREENED_COULOMB, COORDINATES.GLOBAL, 6, aewald);
                double[] dummy1 = new double[3], dummy2 = new double[3], dummy3 = new double[3];
                double[] permFi = dummy1, permTi = dummy2, permTk = dummy3;
                double[] gX = new double[nAtoms], gY = new double[nAtoms], gZ = new double[nAtoms];
                double[] tX = new double[nAtoms], tY = new double[nAtoms], tZ = new double[nAtoms];
                double[] gxk_local = new double[nAtoms], gyk_local = new double[nAtoms], gzk_local = new double[nAtoms];
                double[] txk_local = new double[nAtoms], tyk_local = new double[nAtoms], tzk_local = new double[nAtoms];
                double[] lgX = new double[nAtoms], lgY = new double[nAtoms], lgZ = new double[nAtoms];
                double[] ltX = new double[nAtoms], ltY = new double[nAtoms], ltZ = new double[nAtoms];
                double[] lxk_local = new double[nAtoms], lyk_local = new double[nAtoms], lzk_local = new double[nAtoms];
                double[] ltxk_local = new double[nAtoms], ltyk_local = new double[nAtoms], ltzk_local = new double[nAtoms];
                double dUdL = 0.0, d2UdL2 = 0.0;

                double dScreendL = 0.0, d2ScreendL2 = 0.0;
                double dCouldL = 0.0, d2CouldL2 = 0.0;
                double dPermdL = 0.0, d2PermdL2 = 0.0;

                /**
                 * Compute screened real space interactions.
                 */
                tensor.setR(r);
                tensor.setMultipolesQI(Qi, Qk);
                tensor.setOperator(OPERATOR.SCREENED_COULOMB);
                tensor.order6();

                double ePermScreened = tensor.multipoleEnergyQI(permFi, permTi, permTk);
                dScreendL = tensor.getdEdF();
                dPermdL = dScreendL;
                d2PermdL2 = d2ScreendL2;

                /**
                 * Subtract away masked Coulomb interactions included in PME.
                 */
                double scale1 = 1.0 - scale;
                double ePermCoulomb = 0.0;
                if (scale1 != 0.0) {
                    tensor.setOperator(OPERATOR.COULOMB);
                    ePermCoulomb = tensor.multipoleEnergyQI(FiC, TiC, TkC);
                    dCouldL = tensor.getdEdF();
                    dPermdL -= dCouldL;
                    d2PermdL2 -= d2CouldL2;

                    permFi[0] -= scale1 * FiC[0];
                    permFi[1] -= scale1 * FiC[1];
                    permFi[2] -= scale1 * FiC[2];
                    permTi[0] -= scale1 * TiC[0];
                    permTi[1] -= scale1 * TiC[1];
                    permTi[2] -= scale1 * TiC[2];
                    permTk[0] -= scale1 * TkC[0];
                    permTk[1] -= scale1 * TkC[1];
                    permTk[2] -= scale1 * TkC[2];
                }

                final double ePerm = selfScale * l2 * (ePermScreened - scale1 * ePermCoulomb);

                if (!gradient) {
                    return ePerm;
                }

                double scalar = ELECTRIC * selfScale * l2;
                gX[i] += scalar * permFi[0];
                gY[i] += scalar * permFi[1];
                gZ[i] += scalar * permFi[2];
                tX[i] += scalar * permTi[0];
                tY[i] += scalar * permTi[1];
                tZ[i] += scalar * permTi[2];
                gxk_local[k] -= scalar * permFi[0];
                gyk_local[k] -= scalar * permFi[1];
                gzk_local[k] -= scalar * permFi[2];
                txk_local[k] += scalar * permTk[0];
                tyk_local[k] += scalar * permTk[1];
                tzk_local[k] += scalar * permTk[2];

                double pref1 = 0.0, pref2 = 0.0;
                if (lambdaTerm) {
                    if (System.getProperty("pme-S-qi") != null) {
                        double S = selfScale * l2, dSdL = selfScale, d2SdL2 = 0.0;
                        double P = ePermScreened - scale1 * ePermCoulomb, dPdL = dPermdL, d2PdL2 = d2PermdL2;
                        double F = lAlpha, dFdL = dlAlpha, d2FdL2 = d2lAlpha;
                        double dPdF = (dFdL != 0.0) ? dPdL / dFdL : 0.0;
                        double d2PdF2 = (d2FdL2 != 0.0) ? d2PdL2 / d2FdL2 : 0.0;
                        dUdL += (dSdL * l2 * P) + (S * dPdF * dFdL);
                        d2UdL2 += (d2SdL2 * P) + (dSdL * S * dPdF * dFdL)
                                + ((dSdL * dPdF) + (S * d2PdF2)) * dFdL
                                + (S * dPdF * d2FdL2);
                    } else {
                        double dEdL = dPermdL;
                        double d2EdL2 = d2PermdL2;
                        dUdL += selfScale * (dEdLSign * dlPowPerm * ePerm + l2 * dlAlpha * dEdL);
                        d2UdL2 += selfScale * (dEdLSign * (d2lPowPerm * ePerm
                                + dlPowPerm * dlAlpha * dEdL
                                + dlPowPerm * dlAlpha * dEdL)
                                + l2 * d2lAlpha * dEdL
                                + l2 * dlAlpha * dlAlpha * d2EdL2);
                    }

                    /**
                     * This is dU/dL/dX for the first term of dU/dL: d[dlPow *
                     * ereal]/dx
                     */
                    scalar = ELECTRIC * selfScale * dEdLSign * dlPowPerm;
                    pref1 = scalar;
                    lgX[i] += scalar * permFi[0];
                    lgY[i] += scalar * permFi[1];
                    lgZ[i] += scalar * permFi[2];
                    ltX[i] += scalar * permTi[0];
                    ltY[i] += scalar * permTi[1];
                    ltZ[i] += scalar * permTi[2];
                    lxk_local[k] -= scalar * permFi[0];
                    lyk_local[k] -= scalar * permFi[1];
                    lzk_local[k] -= scalar * permFi[2];
                    ltxk_local[k] += scalar * permTk[0];
                    ltyk_local[k] += scalar * permTk[1];
                    ltzk_local[k] += scalar * permTk[2];
                    /**
                     * Add in dU/dL/dX for the second term of dU/dL:
                     * d[lPow*dlAlpha*dRealdL]/dX
                     */
                    // No additional call to MT; use 6th order tensor instead.
                    scalar = ELECTRIC * selfScale * l2 * dlAlpha;
                    pref2 = scalar;
                    lgX[i] += scalar * permFi[0];
                    lgY[i] += scalar * permFi[1];
                    lgZ[i] += scalar * permFi[2];
                    ltX[i] += scalar * permTi[0];
                    ltY[i] += scalar * permTi[1];
                    ltZ[i] += scalar * permTi[2];
                    lxk_local[k] -= scalar * permFi[0];
                    lyk_local[k] -= scalar * permFi[1];
                    lzk_local[k] -= scalar * permFi[2];
                    ltxk_local[k] += scalar * permTk[0];
                    ltyk_local[k] += scalar * permTk[1];
                    ltzk_local[k] += scalar * permTk[2];
                }

                return ePerm;
            }

//            private double pairPol(double[] r, double[] Qi, double[] Qk,
//                    double[] ui, double[] uk, double[] uiCR, double[] ukCR,
//                    double damp, double aiak) {
//
//                /**
//                 * Compute screened real space interactions.
//                 */
//                tensor.setR(r);
//                tensor.setMultipolesQI(Qi, Qk);
//                tensor.setDipolesQI(ui, uiCR, uk, ukCR);
//                tensor.setOperator(OPERATOR.SCREENED_COULOMB);
//                tensor.order5QI();
//
//                double mutualScale = 1.0;
//                if (polarization == Polarization.DIRECT) {
//                    mutualScale = 0.0;
//                }
//
//                double ePolScreened = tensor.polarizationEnergyQI(
//                        1.0, 1.0, mutualScale, polFi, polTi, polTk);
//
//                /**
//                 * Subtract away masked Coulomb interactions included in PME.
//                 */
//                double scaled1 = 1.0 - scaled;
//                double scalep1 = 1.0 - scalep;
//                double ePolCoulomb = 0.0;
//                if (scaled1 != 0.0 || scalep1 != 0.0) {
//                    tensor.setOperator(OPERATOR.COULOMB);
//                    tensor.order5QI();
//                    ePolCoulomb += tensor.polarizationEnergyQI(
//                            scaled1, scalep1, 0.0, FiC, TiC, TkC);
//                    polFi[0] -= FiC[0];
//                    polFi[1] -= FiC[1];
//                    polFi[2] -= FiC[2];
//                    polTi[0] -= TiC[0];
//                    polTi[1] -= TiC[1];
//                    polTi[2] -= TiC[2];
//                    polTk[0] -= TkC[0];
//                    polTk[1] -= TkC[1];
//                    polTk[2] -= TkC[2];
//                }
//
//                /**
//                 * Subtract away Thole Damped interactions included in PME.
//                 */
//                double eThole = 0.0;
//                tensor.setTholeDamping(damp, aiak);
//                boolean applyThole = tensor.applyDamping();
//                if (applyThole) {
//                    tensor.setOperator(OPERATOR.THOLE_FIELD);
//                    tensor.order4QI();
//                    tensor.setDipolesQI(ui, uiCR, uk, ukCR);
//                    eThole = tensor.polarizationEnergyQI(scaled, scalep, mutualScale, FiT, TiT, TkT);
//                    polFi[0] -= FiT[0];
//                    polFi[1] -= FiT[1];
//                    polFi[2] -= FiT[2];
//                    polTi[0] -= TiT[0];
//                    polTi[1] -= TiT[1];
//                    polTi[2] -= TiT[2];
//                    polTk[0] -= TkT[0];
//                    polTk[1] -= TkT[1];
//                    polTk[2] -= TkT[2];
//                }
//
//                final double e = selfScale * 0.5 * (ePolScreened - ePolCoulomb - eThole);
//                if (!(gradient || lambdaTerm || esvTerm)) {
//                    return polarizationScale * e;
//                }
//
//                double scalar = ELECTRIC * polarizationScale * selfScale;
//                gX[i] += scalar * polFi[0];
//                gY[i] += scalar * polFi[1];
//                gZ[i] += scalar * polFi[2];
//                tX[i] += scalar * polTi[0];
//                tY[i] += scalar * polTi[1];
//                tZ[i] += scalar * polTi[2];
//                gxk_local[k] -= scalar * polFi[0];
//                gyk_local[k] -= scalar * polFi[1];
//                gzk_local[k] -= scalar * polFi[2];
//                txk_local[k] += scalar * polTk[0];
//                tyk_local[k] += scalar * polTk[1];
//                tzk_local[k] += scalar * polTk[2];
//                if (lambdaTerm) {
//                    dUdL += dEdLSign * dlPowPol * e;
//                    d2UdL2 += dEdLSign * d2lPowPol * e;
//                    scalar = ELECTRIC * dEdLSign * dlPowPol * selfScale;
//                    lgX[i] += scalar * polFi[0];
//                    lgY[i] += scalar * polFi[1];
//                    lgZ[i] += scalar * polFi[2];
//                    ltX[i] += scalar * polTi[0];
//                    ltY[i] += scalar * polTi[1];
//                    ltZ[i] += scalar * polTi[2];
//                    lxk_local[k] -= scalar * polFi[0];
//                    lyk_local[k] -= scalar * polFi[1];
//                    lzk_local[k] -= scalar * polFi[2];
//                    ltxk_local[k] += scalar * polTk[0];
//                    ltyk_local[k] += scalar * polTk[1];
//                    ltzk_local[k] += scalar * polTk[2];
//                }
//                return polarizationScale * e;
//            }
            private void applyScaleFactors(Atom ai) {
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

            private void resetScaleFactors(Atom ai) {
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

            private void logPolarizationError(double ei, int i, int k, double indI[], double indK[]) {
                double r2 = 0.0;
                logger.info(crystal.getUnitCell().toString());
                logger.info(atoms[i].toString());
                logger.info(format(" with induced dipole: %8.3f %8.3f %8.3f",
                        indI[0], indI[1], indI[2]));
                logger.info(atoms[k].toString());
                logger.info(format(" with induced dipole: %8.3f %8.3f %8.3f",
                        indK[0], indK[1], indK[2]));
                String message = String.format(" %s\n "
                        + "%s\n with induced dipole: %8.3f %8.3f %8.3f\n "
                        + "%s\n with induced dipole: %8.3f %8.3f %8.3f\n"
                        + " The pol. energy for atoms %d and %d (%d) is %10.6f at %10.6f A.",
                        crystal.getUnitCell(), atoms[i], indI[0], indI[1], indI[2],
                        atoms[k], indK[0], indK[1], indK[2], i + 1, k + 1, iSymm, ei, sqrt(r2));
                throw new EnergyException(message, true);
            }

            private void logPermanentError(double ei, int i, int k) {
                double r2 = 0.0;
                logger.info(crystal.getUnitCell().toString());
                logger.info(atoms[i].toString());
                logger.info(atoms[k].toString());
                String message = String.format(" %s\n %s\n %s\n The permanent "
                        + "multipole energy between atoms %d and %d (%d) is "
                        + "%16.8f at %16.8f A.", crystal.getUnitCell(), atoms[i],
                        atoms[k], i, k, iSymm, ei, sqrt(r2));
                throw new EnergyException(message, true);
            }

        }
    }

    private class ReciprocalEnergyRegion extends ParallelRegion {

        private final double aewald1 = -ELECTRIC * aewald / SQRT_PI;
        private final double aewald2 = 2.0 * aewald * aewald;
        private final double aewald3 = -2.0 / 3.0 * ELECTRIC * aewald * aewald * aewald / SQRT_PI;
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
                        if (esvTerm && esvAtoms[i]) {
                            int esvi = esvSystem.exthEsvId(i);
                            // TODO verify that this is an add when eSelf above is += ...
                            esvRealSpaceDeriv[esvi].addAndGet(eSelf * dlPowPerm * dEdLSign);    // TODO
                        }
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

                double dUdL = 0.0;
                double d2UdL2 = 0.0;
                for (int i = lb; i <= ub; i++) {
                    if (use[i]) {
                        final boolean esvi = esvAtoms[i];
                        final double phi[] = cartMultipolePhi[i];
                        final double mpole[] = multipole[i];
                        final double fmpole[] = fracMultipoles[i];

                        double e = mpole[t000] * phi[t000] + mpole[t100] * phi[t100]
                                + mpole[t010] * phi[t010] + mpole[t001] * phi[t001]
                                + oneThird * (mpole[t200] * phi[t200]
                                + mpole[t020] * phi[t020]
                                + mpole[t002] * phi[t002]
                                + 2.0 * (mpole[t110] * phi[t110]
                                + mpole[t101] * phi[t101]
                                + mpole[t011] * phi[t011]));
                        eRecip += e;
                        if (gradient || lambdaTerm || esvTerm) {
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
                            if (esvTerm && esvi) {
                                // TODO multiply esv chain term
                                final int idxi = esvSystem.exthEsvId(i);
                                esvRealSpaceDeriv[idxi].addAndGet(dEdLSign * dlPowPerm * e);    // TODO check
                            }
                        }
                    }
                }

                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(0.5 * dUdL * ELECTRIC);
                    sharedd2EdLambda2.addAndGet(0.5 * d2UdL2 * ELECTRIC);
                }
//                if (esvTerm) {    // TODO REMOVE (Add as we go?)
//                    for (int i = 0; i < numESVs; i++) {
//                        esvRealSpaceDeriv[i].addAndGet(0.5 * dUdEsv[i] * ELECTRIC);      // TODO missing permanentScale?
//                    }
//                }
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
            private double ldhgX[][], ldhgY[][], ldhgZ[][], ldhtX[][], ldhtY[][], ldhtZ[][];
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
                if (esvTerm) {
                    for (int i = 0; i < nAtoms; i++) {
                        final int idxi = esvSystem.exthEsvId(i);    // TODO determine eSelf is += and this is add?
                        esvRealSpaceDeriv[idxi].addAndGet(dEdLSign * dlPowPol * eSelf);
                    }
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
                if (esvTerm) {
                    for (int i = 0; i < nAtoms; i++) {
                        final int idxi = esvSystem.exthEsvId(i);
                        esvRealSpaceDeriv[idxi].addAndGet(dEdLSign * dlPowPol * eRecip);
                    }
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
                    x[i] = atom.getX();
                    y[i] = atom.getY();
                    z[i] = atom.getZ();
                    use[i] = atom.getUse();

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
            private final double localDipole[] = new double[3];
            private final double localQuadrupole[][] = new double[3][3];
            private final double frameCoords[][] = new double[4][3];
            private final double xAxis[] = new double[3];
            private final double yAxis[] = new double[3];
            private final double zAxis[] = new double[3];
            private final double rotmat[][] = new double[3][3];
            private final double globalDipole[] = new double[3];
            private final double globalQuadrupole[][] = new double[3][3];
            private double chargeScale, dipoleScale, traceScale;
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
                traceScale = 1.0;
                if (!useCharges) {
                    chargeScale = 0.0;
                }
                if (!useDipoles) {
                    dipoleScale = 0.0;
                }
                if (!useQuadrupoles) {
                    traceScale = 0.0;
                }
            }

            @Override
            public void run(int lb, int ub) {
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    final double x[] = coordinates[iSymm][0];
                    final double y[] = coordinates[iSymm][1];
                    final double z[] = coordinates[iSymm][2];
                    for (int ii = lb; ii <= ub; ii++) {
                        Atom atom = atoms[ii];
                        final double in[] = localMultipole[ii];
                        final double out[] = globalMultipole[iSymm][ii];
                        double elecScale = 1.0;
                        if (!atom.getElectrostatics()) {
                            elecScale = 0.0;
                        }
                        if (rotateMultipoles) {
                            localOrigin[0] = x[ii];
                            localOrigin[1] = y[ii];
                            localOrigin[2] = z[ii];
                            int referenceSites[] = axisAtom[ii];
                            int nSites = 0;
                            if (referenceSites != null) {
                                nSites = referenceSites.length;
                            }
                            for (int i = 0; i < nSites; i++) {
                                int index = referenceSites[i];
                                frameCoords[i][0] = x[index];
                                frameCoords[i][1] = y[index];
                                frameCoords[i][2] = z[index];
                            }
                            for (int i = 0; i < 3; i++) {
                                zAxis[i] = 0.0;
                                xAxis[i] = 0.0;
                                yAxis[i] = 0.0;
                                globalDipole[i] = 0.0;
                                for (int j = 0; j < 3; j++) {
                                    globalQuadrupole[i][j] = 0.0;
                                }
                            }
                            if (nSites < 1) {
                                out[t000] = in[0] * chargeScale * elecScale;
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
                                polarizability[ii] = polarizeType.polarizability * elecScale;
                                continue;
                            }
                            localDipole[0] = in[t100];
                            localDipole[1] = in[t010];
                            localDipole[2] = in[t001];
                            localQuadrupole[0][0] = in[t200];
                            localQuadrupole[1][1] = in[t020];
                            localQuadrupole[2][2] = in[t002];
                            localQuadrupole[0][1] = in[t110];
                            localQuadrupole[0][2] = in[t101];
                            localQuadrupole[1][2] = in[t011];
                            localQuadrupole[1][0] = in[t110];
                            localQuadrupole[2][0] = in[t101];
                            localQuadrupole[2][1] = in[t011];
                            // Check for chiral flipping.
                            if (frame[ii] == MultipoleType.MultipoleFrameDefinition.ZTHENX
                                    && referenceSites.length == 3) {
                                checkMultipoleChirality(frame[ii],
                                        localOrigin, frameCoords, localDipole, localQuadrupole);
                            }
                            // Do the rotation.
                            getRotationMatrix(frame[ii], localOrigin, frameCoords, rotmat);
                            rotateMultipole(rotmat, localDipole, localQuadrupole,
                                    globalDipole, globalQuadrupole);
                            out[t000] = in[0] * chargeScale * elecScale;
                            out[t100] = globalDipole[0] * dipoleScale * elecScale;
                            out[t010] = globalDipole[1] * dipoleScale * elecScale;
                            out[t001] = globalDipole[2] * dipoleScale * elecScale;
                            out[t200] = globalQuadrupole[0][0] * traceScale * elecScale;
                            out[t020] = globalQuadrupole[1][1] * traceScale * elecScale;
                            out[t002] = globalQuadrupole[2][2] * traceScale * elecScale;
                            out[t110] = globalQuadrupole[0][1] * traceScale * elecScale;
                            out[t101] = globalQuadrupole[0][2] * traceScale * elecScale;
                            out[t011] = globalQuadrupole[1][2] * traceScale * elecScale;
                        } else {
                            /**
                             * Do not perform multipole rotation, which helps to
                             * isolate torque vs. non-torque pieces of the
                             * multipole energy gradient.
                             */
                            out[t000] = in[t000] * chargeScale * elecScale;
                            out[t100] = in[t100] * dipoleScale * elecScale;
                            out[t010] = in[t010] * dipoleScale * elecScale;
                            out[t001] = in[t001] * dipoleScale * elecScale;
                            out[t200] = in[t200] * traceScale * elecScale;
                            out[t020] = in[t020] * traceScale * elecScale;
                            out[t002] = in[t002] * traceScale * elecScale;
                            out[t110] = in[t110] * traceScale * elecScale;
                            out[t101] = in[t101] * traceScale * elecScale;
                            out[t011] = in[t011] * traceScale * elecScale;
                        }
                        PolarizeType polarizeType = atoms[ii].getPolarizeType();
                        polarizability[ii] = polarizeType.polarizability * elecScale;
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
            private double ldhg[][][];
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
                //double uwcos = dotK(u, w);
                //double uwsin = sqrt(1.0 - uwcos * uwcos);
                //double vwcos = dotK(v, w);
                //double vwsin = sqrt(1.0 - vwcos * vwcos);
                /*
                 * Negative of dotK product of torque with unit vectors gives
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
                        //double uscos = dotK(u, s);
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
            if (!MultipoleType.assignMultipole(atoms[i], forceField,
                    localMultipole[i], i, axisAtom, frame)) {
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
            //double uwcos = dotK(u, w);
            //double uwsin = sqrt(1.0 - uwcos * uwcos);
            //double vwcos = dotK(v, w);
            //double vwsin = sqrt(1.0 - vwcos * vwcos);
            /*
             * Negative of dotK product of torque with unit vectors gives result
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
                    //double uscos = dotK(u, s);
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
     * {@inheritDoc}
     *
     * Set the electrostatic lambda scaling factor.
     */
    @Override
    public void setLambda(double lambda) {
        assert (lambda >= 0.0 && lambda <= 1.0);
        this.lambda = lambda;

        // Also do this if ESVs have been changed.
        if (!initSoftCore) {
            initSoftCore(false, true);
        }

        /**
         * f = sqrt(r^2 + lAlpha) df/dL = -alpha * (1.0 - lambda) / f g = 1 /
         * sqrt(r^2 + lAlpha) dg/dL = alpha * (1.0 - lambda) / (r^2 +
         * lAlpha)^(3/2) define dlAlpha = alpha * 1.0 - lambda) then df/dL =
         * -dlAlpha / f and dg/dL = dlAlpha * g^3
         */
        lAlpha = permLambdaAlpha * (1.0 - lambda) * (1.0 - lambda);
        dlAlpha = -2.0 * permLambdaAlpha * (1.0 - lambda);
        d2lAlpha = 2.0 * permLambdaAlpha;

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

        if (esvTerm) {
            updateEsvLambda();
        }
    }

    /**
     * Attach system with extended variable such as titrations.
     */
    public void attachExtendedSystem(ExtendedSystem system) {
        esvTerm = true;
        esvSystem = system;
        numESVs = system.n();
        initAtomArrays();
        updateEsvLambda();
        
        /**
         * Enforce ESV-handling requirements slash best practices.
         */
        if (esvTerm) {
            if (permLambdaAlpha != 2.0 || permLambdaExponent != 1.0 || doLigandVaporElec || doNoLigandCondensedSCF
                    || !intermolecularSoftcore || !intramolecularSoftcore) {
                logf(" (EsvSys) Nonstandard PME-ESV configuration: %.2f %.2f %b %b %b %b\n"
                     + "          Enforcing failsafe defaults:       %.2f %.2f %b %b %b %b",
                        permLambdaAlpha, permLambdaExponent, doLigandVaporElec, doNoLigandCondensedSCF,
                        intermolecularSoftcore, intramolecularSoftcore, 2.0, 1.0, false, false, true, true);
            }
            permLambdaAlpha = 2.0;
            permLambdaExponent = 1.0;
            doLigandVaporElec = false;
            doNoLigandCondensedSCF = false;
            intermolecularSoftcore = true;
            intramolecularSoftcore = true;
        }
    }

    public void detachExtendedSystem() {
        fill(esvAtoms, false);
        esvTerm = false;
        esvSystem = null;
        esvLambda = null;
        esvRealSpaceDeriv = null;
        numESVs = 0;
        initSoftCore(true, false);
    }
    
    /**
     * Precalculate PME lambda terms; must be called when either OSRW or ESV lambdas are propagated.
     * TODO Assess cost of storing {lAlpha, lPowPol, their derivs} in the inner loops vs precomputing.
     * TODO Note! If we have atom -> esv *PAIR* (potentially) lookups, then it isn't n^2...
     */
    public void updateEsvLambda() {
        if (!esvTerm) {
            return; // ESV removal in detach().
        }
        numESVs = esvSystem.n();
        esvRealSpaceDeriv = new SharedDouble[numESVs];
        
        /**
         * Force rebuild of the softcore lists.
         */
        initSoftCore(true, true);
        
        
        for (int i = 0; i < nAtoms; i++) {
            final Atom ai = atoms[i];
            final double li = esvSystem.exthLambda(i);
            
            final double L = (lambdaTerm) ? lambda*li : li;
            
            // Set permanent electrostatics scaling.
            if (permLambdaExponent != 1.0) {
                logger.severe("Attempted to use ESV with non-unity lambda exponent.");
            }
            // We assume (as in VdW), that permLambdaExponent is equal to unity.
            // TODO HERE
            
            double polWindow = polLambdaEnd - polLambdaStart;
            double polLambdaScale = 1.0 / polWindow;
            
            lPowPerm = L;          // L^permLambdaExponent
            dlPowPerm = 1.0;       // (1-permExponent) * L^0;
            d2lPowPerm = 0.0;      // 0*...

            lAlpha = permLambdaAlpha * (1.0 - L) * (1.0 - L);       
            dlAlpha = -2.0 * permLambdaAlpha * (1.0 - L);
            d2lAlpha = 2.0 * permLambdaAlpha;
            
            // Set polarization scaling.
            lPowPol = 1.0;
            dlPowPol = 0.0;
            d2lPowPol = 0.0;
            double polLambda = 0.0;
            if (L < polLambdaStart) {
                lPowPol = 0.0;
            } else if (L <= polLambdaEnd) {
                polLambda = polLambdaScale * (L - polLambdaStart);
                lPowPol = pow(polLambda, polLambdaExponent);
                if (polLambdaExponent >= 1.0) {
                    dlPowPol = polLambdaExponent * pow(polLambda, polLambdaExponent - 1.0);
                    if (polLambdaExponent >= 2.0) {
                        d2lPowPol = polLambdaExponent * (polLambdaExponent - 1.0)
                                * pow(polLambda, polLambdaExponent - 2.0);
                    }
                }
                // Chain rule d/t shrunken (ie non-[0,1]) range.
                dlPowPol *= polarizationScale;
                d2lPowPol *= (polarizationScale * polarizationScale);
            }
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
        if (shareddEdLambda == null || !lambdaTerm) {
            logger.warning("Tried to get null/off lambda derivative.");
            return 0.0;
        }
        double dEdL = shareddEdLambda.get();
        if (generalizedKirkwoodTerm) {
            dEdL += generalizedKirkwood.getdEdL();
        }
        return dEdL;
    }

    @Override
    public double[] getdEdEsv() {
        if (!esvTerm) {
            throw new IllegalStateException();
        }
        StringBuilder sb = new StringBuilder();
        double[] dEdEsv = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            if (doPermanentRealSpace) {
                sb.append(format(" RealSpaceDeriv%d:  %16.8f", i,
                        dEdEsv[i] = esvRealSpaceDeriv[i].get()));
            }
            if (reciprocalSpaceTerm) {
                sb.append(format(" RealSpaceDeriv%d:  %16.8f", i,
                        dEdEsv[i] = esvRealSpaceDeriv[i].get()));
                // TODO dEdEsv[i] = esvReciprocalDeriv[i].get();
            }
            if (polarization != Polarization.NONE) {
                // TODO dEdEsv[i] = esvPolDeriv[i].get();
            }
            if (generalizedKirkwoodTerm) {
                // TODO add GK-ESV  // dEdEsv[i] = esvGkDeriv[i].get();
            }
        }
        return dEdEsv;
    }
    
    public double getdEdEsv(int esvID) {
        return getdEdEsv()[esvID];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        if (sharedd2EdLambda2 == null || !lambdaTerm) {
            logger.warning("Tried to get null/off lambda (second) derivative.");
            return 0.0;
        }
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
        if (lambdaGrad == null || !lambdaTerm) {
            return;
        }
        /**
         * Note that the Generalized Kirkwood contributions are already in the
         * lambdaGrad array.
         */
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                gradient[index++] += lambdaGrad[0][0][i];
                gradient[index++] += lambdaGrad[0][1][i];
                gradient[index++] += lambdaGrad[0][2][i];
            }
        }
    }

    private void computeInduceDipoleField() {
        try {
            if (nSymm > 1) {
                parallelTeam.execute(expandInducedDipolesRegion);
            }

            if (reciprocalSpaceTerm && aewald > 0.0) {
                reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipoleCR, use);
            }
            sectionTeam.execute(inducedDipoleFieldRegion);
            if (reciprocalSpaceTerm && aewald > 0.0) {
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
                        " %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * TO_SECONDS));
            }
            /**
             * If the RMS Debye change increases, fail the SCF process.
             */
            if (eps > previousEps) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Fatal SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
                throw new EnergyException(message, false);
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
                throw new EnergyException(message, false);
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
                    TO_SECONDS * directTime));
            startTime = System.nanoTime() - startTime;
            sb.append(format(" Total:                   %7.4f",
                    startTime * TO_SECONDS));
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
                    double ipolar;
                    if (polarizability[i] > 0) {
                        ipolar = 1.0 / polarizability[i];
                        rsd[0][i] = (directDipole[i][0] - inducedDipole[0][i][0]) * ipolar + field[0][0][i];
                        rsd[1][i] = (directDipole[i][1] - inducedDipole[0][i][1]) * ipolar + field[0][1][i];
                        rsd[2][i] = (directDipole[i][2] - inducedDipole[0][i][2]) * ipolar + field[0][2][i];
                        rsdCR[0][i] = (directDipoleCR[i][0] - inducedDipoleCR[0][i][0]) * ipolar + fieldCR[0][0][i];
                        rsdCR[1][i] = (directDipoleCR[i][1] - inducedDipoleCR[0][i][1]) * ipolar + fieldCR[0][1][i];
                        rsdCR[2][i] = (directDipoleCR[i][2] - inducedDipoleCR[0][i][2]) * ipolar + fieldCR[0][2][i];
                    } else {
                        rsd[0][i] = 0.0;
                        rsd[1][i] = 0.0;
                        rsd[2][i] = 0.0;
                        rsdCR[0][i] = 0.0;
                        rsdCR[1][i] = 0.0;
                        rsdCR[2][i] = 0.0;
                    }
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
                    if (polarizability[i] > 0) {
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
                    } else {
                        inducedDipole[0][i][0] = 0.0;
                        inducedDipole[0][i][1] = 0.0;
                        inducedDipole[0][i][2] = 0.0;
                        vec[0][i] = 0.0;
                        vec[1][i] = 0.0;
                        vec[2][i] = 0.0;
                        inducedDipoleCR[0][i][0] = 0.0;
                        inducedDipoleCR[0][i][1] = 0.0;
                        inducedDipoleCR[0][i][2] = 0.0;
                        vecCR[0][i] = 0.0;
                        vecCR[1][i] = 0.0;
                        vecCR[2][i] = 0.0;
                    }

                    // Compute dotK product of the conjugate vector and new residual.
                    dot += conj[0][i] * vec[0][i]
                            + conj[1][i] * vec[1][i]
                            + conj[2][i] * vec[2][i];
                    dotCR += conjCR[0][i] * vecCR[0][i]
                            + conjCR[1][i] * vecCR[1][i]
                            + conjCR[2][i] * vecCR[2][i];
                    // Compute dotK product of the previous residual and preconditioner.
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

    private void log(String type, int i, int k, double r, double eij) {
        logger.info(String.format("%s %6d-%s %6d-%s %10.4f  %10.4f",
                type, atoms[i].xyzIndex, atoms[i].getAtomType().name,
                atoms[k].xyzIndex, atoms[k].getAtomType().name, r, eij));
    }

    /**s
     * Number of unique tensors for given order.
     */
    private static final int tensorCount = MultipoleTensor.tensorCount(3);
    private static final double oneThird = 1.0 / 3.0;

}
