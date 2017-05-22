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

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
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
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Torsion;
import ffx.potential.extended.ExtUtils;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedSystem.ExtendedSystemConfig;
import ffx.potential.nonbonded.ReciprocalSpace.FFTMethod;
import ffx.potential.nonbonded.ScfPredictor.PredictorMode;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.utils.EnergyException;

import static ffx.numerics.Erf.erfc;
import static ffx.numerics.MultipoleTensor.OPERATOR.COULOMB;
import static ffx.numerics.MultipoleTensor.OPERATOR.SCREENED_COULOMB;
import static ffx.numerics.MultipoleTensor.OPERATOR.THOLE_FIELD;
import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.parameters.MultipoleType.ELECTRIC;
import static ffx.potential.parameters.MultipoleType.checkMultipoleChirality;
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
import static ffx.potential.parameters.MultipoleType.zeroD;
import static ffx.potential.parameters.MultipoleType.zeroM;

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

    private double permanentMultipoleEnergy;
    private double polarizationEnergy;
    private double generalizedKirkwoodEnergy;

    private double permanentRealSpaceEnergy;
    private double permanentSelfEnergy;
    private double permanentReciprocalEnergy;
    private double inducedRealSpaceEnergy;
    private double inducedSelfEnergy;
    private double inducedReciprocalEnergy;

    /**
     * Flag to indicate use of generalized Kirkwood.
     */
    private boolean generalizedKirkwoodTerm;
    /**
     * Unit cell and spacegroup information.
     */
    private Crystal crystal;
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
     * *****************************************
     * Extended System Variables 1. Scaling only, rather than softcoring,
     * pending consideration of scaled-soft interaction. 2. Interpolated
     * multipoles preloaded in updateEsvLambda(). 3. Polarizability affected
     * only for disappearing atoms (ie titrating H+).
     */
    private boolean esvTerm = false;
    private ExtendedSystem esvSystem;
    private ExtendedSystemConfig esvConfig = ExtendedSystem.DefaultConfig;
    private int numESVs = 0;
    /**
     * EsvID index to gradient arrays. [atom]
     */
    private Integer[] esvIndex;
    /**
     * Denotes if the atomic multipole is scaled by an ESV. [atom]
     */
    private boolean[] esvAtomsScaled;
    /**
     * Denotes if the atomic polarizability is scaled by an ESV. Note that only
     * titrating hydrogens have scaled polarizabilities. [atom]
     */
    private boolean[] esvAtomsScaledAlpha;
    /**
     * Per ESV lambda. [atom]
     */
    private double[] esvLambda;
    /**
     * ESV multipole derivatives; equal to Protonated value - Deprotonated
     * value. [sym][atom][10]
     */
    private double[][][] esvMultipoleDot;
    /**
     * Atomic polarizabilities of the AMOEBA model.
     */
    private double[] unscaledPolarizability;
    /**
     * Direct induced dipoles without polarizability scaling. [atom][3]
     */
    private double[][] unscaledDirectDipole;
    private double[][] unscaledDirectDipoleCR;
    /**
     * Induced dipoles without polarizability scaling. [sym][atom][3]
     */
    private double[][][] unscaledInducedDipole;
    private double[][][] unscaledInducedDipoleCR;
    /**
     * Induced dipole field without polarizability scaling. [sym][3][atom]
     */
    private double[][][] unscaledField;
    private double[][][] unscaledFieldCR;
    /**
     * Cartesian permanent multipole phi calculated from dotted multipoles.
     */
    private double cartMultipoleDotPhi[][];
    /**
     * Cartesian Dipole Phi without polarizability scaling. [3][atom]
     */
    private double[][] unscaledCartDipolePhi;
    private double[][] unscaledCartDipolePhiCR;
    /**
     * Shared ESV derivative component arrays. [numESVs]
     */
    private SharedDouble[] esvPermRealDeriv_shared;
    private SharedDouble[] esvPermSelfDeriv_shared;
    private SharedDouble[] esvPermRecipDeriv_shared;
    private SharedDouble[] esvInducedRealDeriv_shared;
    private SharedDouble[] esvInducedSelfDeriv_shared;
    private SharedDouble[] esvInducedRecipDeriv_shared;
    private final boolean printInducedDipoles = prop("pme.printInducedDipoles", false);
    private final boolean printDecomposition = prop("pme.decompose", false);

    /**
     * *****************************************
     */
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
     * Reference to the force field being used.
     */
    private final ForceField forceField;
    /**
     * Optionally predict induced dipoles prior to the SCF calculation.
     */
    private ScfPredictor scfPredictor = null;
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
    private LambdaMode lambdaMode = LambdaMode.OFF;
    /**
     * If lambdaTerm is true, some ligand atom interactions with the environment
     * are being turned on/off.
     */
    private boolean lambdaTerm;
    private double lambda = 1.0;
    /**
     * The polarization Lambda value goes from 0.0 .. 1.0 as the global lambda
     * value varies between polarizationLambdaStart .. 1.0.
     */
    private double polLambda = 1.0;
    /**
     * Constant α in: r' = sqrt(r^2 + α*(1 - L)^2)
     */
    private double polLambdaAlpha = 1.0;
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
     * Start turning on permanent electrostatics later in the Lambda path.
     */
    private double permLambdaStart = 0.5;
    private double permLambdaEnd = 1.0;
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
     * sc1 = lAlpha = α*(1 - L)^2 where L = (lambda - start) * (1 / (start -
     * end))
     */
    private double sc1 = 0.0;
    private double dsc1dL = 0.0;
    private double d2sc1dL2 = 0.0;
    private double dEdLSign = 1.0;
    /**
     * sc2 = lPowPerm = L^permanentLambdaExponent where L = (lambda - start) *
     * (1 / (start - end))
     */
    private double sc2 = 1.0;
    private double dsc2dL = 0.0;
    private double d2sc2dL2 = 0.0;
    private boolean realSpaceTerm;
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
    private int molecule[];

    /**
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand). 2.) Uenv = The polarization energy of the system without the
     * ligand. 3.) Uligand = The polarization energy of the ligand by itself.
     * 4.) Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     *
     * Set polarizationScale to L for part 1. Set polarizationScale to (1-L) for
     * parts 2 & 3. This is OSRW_lambda only; ESV_lambda scales by overwriting
     * polarizabilities in InitializationRegion.
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
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand). 2.) Uenv = The polarization energy of the system without the
     * ligand. 3.) Uligand = The polarization energy of the ligand by itself.
     * 4.) Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
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
    /**
     * By way of keeping up with the Tinker codebase; an as-yet not fully
     * explained scale factor of one-half applied only to 1-4 atoms that share a
     * polarization group.
     */
    private final double p14scale;
    private final double p15scale;
    private final double intra14Scale;

    private double ipdamp[];
    private double thole[];
    private double polarizability[];

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
    private double gradPerm[][][];
    private double gradInduced[][][];
    /**
     * Torque array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] torque;
    private double[][][] torquePerm;
    private double[][][] torqueInduced;
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

    private boolean reciprocalSpaceTerm;
    private final ReciprocalSpace reciprocalSpace;
    private final ReciprocalEnergyRegion reciprocalEnergyRegion;
    private final RealSpaceEnergyRegion realSpaceEnergyRegion;
    private final ReduceRegion reduceRegion;
    private final GeneralizedKirkwood generalizedKirkwood;
    /**
     * Timing variables.
     */
    private final long realSpacePermTimes[];
    private final long realSpaceTimes[];
    private final long realSpaceScfTimes[];
    private long realSpacePermTimeTotal, realSpaceTimeTotal, realSpaceScfTimeTotal;
    private long bornRadiiTotal, gkEnergyTotal;
    private ELEC_FORM elecForm = ELEC_FORM.PAM;
    private static final double TO_SECONDS = 1.0e-9;
    private static final double TO_MS = 1.0e-6;
    private static final boolean skipRealRegion = prop("pme.skipReal", false);
    private static final boolean skipRecipRegion = prop("pme.skipRecip", false);

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
        intra14Scale = forceField.getDouble(ForceFieldDouble.POLAR_14_INTRA, 0.5);
        d11scale = forceField.getDouble(ForceFieldDouble.DIRECT_11_SCALE, 0.0);
        p12scale = forceField.getDouble(ForceFieldDouble.POLAR_12_SCALE, 0.0);
        p13scale = forceField.getDouble(ForceFieldDouble.POLAR_13_SCALE, 0.0);
        p14scale = forceField.getDouble(ForceFieldDouble.POLAR_14_SCALE, 1.0);
        p15scale = forceField.getDouble(ForceFieldDouble.POLAR_15_SCALE, 1.0);
        useCharges = forceField.getBoolean(ForceFieldBoolean.USE_CHARGES, true);
        useDipoles = forceField.getBoolean(ForceFieldBoolean.USE_DIPOLES, true);
        useQuadrupoles = forceField.getBoolean(ForceFieldBoolean.USE_QUADRUPOLES, true);
        rotateMultipoles = forceField.getBoolean(ForceFieldBoolean.ROTATE_MULTIPOLES, true);
//        lambdaTerm = forceField.getBoolean(ForceFieldBoolean.LAMBDATERM, false);
        lambdaTerm = false; // TODO enable

        if (!crystal.aperiodic()) {
            off = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        } else {
            off = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 1000.0);
        }
        double ewaldPrecision = forceField.getDouble(ForceFieldDouble.EWALD_PRECISION, 1.0e-8);
        aewald = forceField.getDouble(ForceFieldDouble.EWALD_ALPHA, ewaldCoefficient(off, ewaldPrecision));
        setEwaldParameters(off, aewald);

        reciprocalSpaceTerm = forceField.getBoolean(ForceFieldBoolean.RECIPTERM, true);

        /**
         * Instantiate the requested SCF predictor; default is a 6th-order least
         * squares method.
         */
        PredictorMode predictorMode = ExtUtils.prop(PredictorMode.class, "scf-predictor", PredictorMode.NONE);
        int predictorOrder = ExtUtils.prop("scf-predictor-order", 6);
        if (predictorMode != PredictorMode.NONE) {
            scfPredictor = new ScfPredictor(predictorMode, predictorOrder, forceField);
        }

        String algorithm = forceField.getString(ForceFieldString.SCF_ALGORITHM, "CG");
        try {
            algorithm = algorithm.replaceAll("-", "_").toUpperCase();
            scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
        } catch (Exception e) {
            scfAlgorithm = SCFAlgorithm.CG;
        }

        generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
        if (generalizedKirkwoodTerm && scfAlgorithm == SCFAlgorithm.CG) {
            scfAlgorithm = SCFAlgorithm.SOR;
            logger.warning("Preconditioner does not yet support GK; setting scf-algorithm=SOR instead.");
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
             * A PERMANENT_LAMBDA_EXPONENT of 2 gives a non-zero d2U/dL2 at the
             * beginning of the permanent schedule. Choosing a power of 3 or
             * greater ensures a smooth dU/dL and d2U/dL2 over the schedule.
             */
            permLambdaExponent = forceField.getDouble(ForceFieldDouble.PERMANENT_LAMBDA_EXPONENT, 2.0);
            if (permLambdaExponent < 2.0) {
                permLambdaExponent = 2.0;
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
            permLambdaStart = forceField.getDouble(ForceFieldDouble.PERMANENT_LAMBDA_START, 0.5);
            permLambdaEnd = forceField.getDouble(ForceFieldDouble.PERMANENT_LAMBDA_END, 1.0);

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
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n Electrostatics       %25s\n", getClass().getSimpleName()));
            sb.append(format("  Polarization:                        %8s\n", polarization.toString()));
            if (polarization == Polarization.MUTUAL) {
                sb.append(format("   SCF Convergence Criteria:          %8.3e\n", poleps));
                if (scfPredictor != null) {
                    sb.append(format("   SCF Predictor:                      %8s\n", scfPredictor.toString()));
                }
                sb.append(format("   SCF Algorithm:                      %8s\n", scfAlgorithm));
                if (scfAlgorithm == SCFAlgorithm.SOR) {
                    sb.append(format("   SOR Parameter:                      %8.3f\n", polsor));
                } else {
                    sb.append(format("   CG Preconditioner Cut-Off:          %8.3f\n", preconditionerCutoff));
                    sb.append(format("   CG Preconditioner Ewald Coefficient:%8.3f\n", preconditionerEwald));
                }
            }
            if (aewald > 0.0) {
                sb.append(format("  Particle-mesh Ewald                  %8s\n", getClass().getSimpleName()));
                sb.append(format("   Ewald Coefficient:                  %8.3f\n", aewald));
                sb.append(format("   Particle Cut-Off:                   %8.3f (A)", off));
            } else {
                sb.append(format("   Electrostatics Cut-Off:             %8.3f (A)\n", off));
            }
            logger.info(sb.toString());
        }

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
        realSpaceEnergyRegion = new RealSpaceEnergyRegion(maxThreads);
        reduceRegion = new ReduceRegion(maxThreads);
        realSpaceTimes = new long[maxThreads];
        realSpacePermTimes = new long[maxThreads];
        realSpaceScfTimes = new long[maxThreads];

        /**
         * Generalized Kirkwood currently requires aperiodic Ewald. The GK
         * reaction field is added to the intra-molecular to give a
         * self-consistent reaction field.
         */
        if (generalizedKirkwoodTerm) {
            if (esvTerm) {
                throw new UnsupportedOperationException();
            }
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
        if (localMultipole == null || localMultipole.length < nAtoms
                || lambdaTerm || esvTerm) {
            localMultipole = new double[nAtoms][10];
            frame = new MultipoleType.MultipoleFrameDefinition[nAtoms];
            axisAtom = new int[nAtoms][];
            cartMultipolePhi = new double[nAtoms][tensorCount];
            if (esvTerm) {
                cartMultipoleDotPhi = new double[nAtoms][tensorCount];
            }
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

            /**
             * Initialize per-thread memory for collecting the gradient, torque,
             * field and chain-rule field.
             */
            grad = new double[maxThreads][3][nAtoms];
            gradPerm = new double[maxThreads][3][nAtoms];
            gradInduced = new double[maxThreads][3][nAtoms];
            torque = new double[maxThreads][3][nAtoms];
            torquePerm = new double[maxThreads][3][nAtoms];
            torqueInduced = new double[maxThreads][3][nAtoms];
            field = new double[maxThreads][3][nAtoms];
            fieldCR = new double[maxThreads][3][nAtoms];
            if (lambdaTerm) {
                lambdaGrad = new double[maxThreads][3][nAtoms];
                lambdaTorque = new double[maxThreads][3][nAtoms];
            }
            isSoft = new boolean[nAtoms];
            fill(isSoft, false);
            use = new boolean[nAtoms];
            fill(use, true);

            coordinates = new double[nSymm][3][nAtoms];
            globalMultipole = new double[nSymm][nAtoms][10];
            inducedDipole = new double[nSymm][nAtoms][3];
            inducedDipoleCR = new double[nSymm][nAtoms][3];
            if (scfPredictor != null) {
                if (esvTerm) {
                    throw new UnsupportedOperationException();
                }
                scfPredictor.setInducedDipoleReferences(inducedDipole, inducedDipoleCR, lambdaTerm);
            }

            /* ESV flag array, initialized regardless of esvTerm. */
            esvAtomsScaled = new boolean[nAtoms];	// True for other ESV residue atoms.
            esvAtomsScaledAlpha = new boolean[nAtoms];
            fill(esvAtomsScaled, false);
            fill(esvAtomsScaledAlpha, false);

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
         * Assign multipole parameters and polarization groups.
         */
        for (int i = 0; i < nAtoms; i++) {
            MultipoleType.multipoleTypeFactory(atoms[i], forceField);
            localMultipole[i] = atoms[i].getMultipoleType().getMultipole();
            axisAtom[i] = atoms[i].getAxisAtomIndices();
            frame[i] = atoms[i].getMultipoleType().frameDefinition;
        }
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

        if (esvTerm) {
            updateEsvLambda();
        }
    }

    /**
     * Pass in atoms that have been assigned electrostatics from a fixed charge
     * force field.
     *
     * @param atoms
     */
    @Override
    public void setFixedCharges(Atom atoms[]) {
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

    /**
     * Initialize a boolean array of soft atoms and, if requested, ligand vapor
     * electrostatics.
     */
    private void initSoftCore(boolean print) {
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
            vaporPermanentSchedule = vacuumNeighborList.getPairwiseSchedule();
            vaporEwaldSchedule = vaporPermanentSchedule;
            vacuumRanges = new Range[maxThreads];
            vacuumNeighborList.setDisableUpdates(forceField.getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false));
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
        if (esvTerm) {
            updateEsvLambda();
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
        permanentRealSpaceEnergy = 0.0;
        permanentSelfEnergy = 0.0;
        permanentReciprocalEnergy = 0.0;
        polarizationEnergy = 0.0;
        inducedRealSpaceEnergy = 0.0;
        inducedSelfEnergy = 0.0;
        inducedReciprocalEnergy = 0.0;
        generalizedKirkwoodEnergy = 0.0;
        interactions = 0;
        gkInteractions = 0;
        for (int i = 0; i < maxThreads; i++) {
            realSpacePermTimes[i] = 0;
            realSpaceTimes[i] = 0;
            realSpaceScfTimes[i] = 0;
        }
        realSpacePermTimeTotal = 0;
        realSpaceTimeTotal = 0;
        realSpaceScfTimeTotal = 0;
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
        realSpaceTerm = true;
        permanentScale = 1.0;
        doPolarization = true;
        polarizationScale = 1.0;

        /**
         * Expand the coordinates and rotate multipoles into the global frame.
         */
        try {
            parallelTeam.execute(initializationRegion);
        } catch (RuntimeException ex) {
            logger.warning("Runtime exception expanding coordinates and rotating multipoles.");
            throw ex;
        } catch (Exception ex) {
            logger.severe("Fatal exception expanding coordinates and rotating multipoles.");
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
            if (logger.isLoggable(Level.FINER)) {
                logger.fine(String.format(" Solvated energy: %20.8f", energyLog));
            }

            /**
             * Condensed phase SCF without ligand atoms.
             */
            if (doNoLigandCondensedSCF) {
                if (esvTerm) {
                    throw new UnsupportedOperationException();
                }
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
                if (esvTerm) {
                    throw new UnsupportedOperationException();
                }
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
        if (logger.isLoggable(Level.FINER)) {
            printRealSpaceTimings();
            if (aewald > 0.0 && reciprocalSpaceTerm) {
                reciprocalSpace.printTimings();
            }
        }

        return permanentMultipoleEnergy + polarizationEnergy;
    }

    private void printRealSpaceTimings() {
        double total = (realSpacePermTimeTotal + realSpaceScfTimeTotal + realSpaceTimeTotal) * TO_SECONDS;

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
            int count = realSpaceEnergyRegion.realSpaceEnergyLoops[i].getCount();
            logger.info(String.format("    %3d   %7.4f %7.4f %7.4f %10d", i,
                    realSpacePermTimes[i] * TO_SECONDS, realSpaceScfTimes[i] * TO_SECONDS,
                    realSpaceTimes[i] * TO_SECONDS, count));
            minPerm = min(realSpacePermTimes[i], minPerm);
            maxPerm = max(realSpacePermTimes[i], maxPerm);
            minSCF = min(realSpaceScfTimes[i], minSCF);
            maxSCF = max(realSpaceScfTimes[i], maxSCF);
            minEnergy = min(realSpaceTimes[i], minEnergy);
            maxEnergy = max(realSpaceTimes[i], maxEnergy);
            minCount = min(count, minCount);
            maxCount = max(count, maxCount);
        }
        int inter = realSpaceEnergyRegion.getInteractions();
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
                realSpacePermTimeTotal * TO_SECONDS, realSpaceScfTimeTotal * TO_SECONDS,
                realSpaceTimeTotal * TO_SECONDS, inter));
    }

    /**
     * 1.) Total system under PBC. A.) Softcore real space for Ligand-Protein
     * and Ligand-Ligand. B.) Reciprocal space scaled by lambda. C.)
     * Polarization scaled by lambda.
     */
    private double condensedEnergy() {
        if (esvTerm) {
            doPolarization = esvConfig.polarization;
            polarizationScale = 1.0;
        } else if (lambdaTerm) {
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
        }

        realSpaceTerm = true;
        permanentScale = sc2;
        dEdLSign = 1.0;

        double energy = computeEnergy(false);

        return energy;
    }

    /**
     * 2.) Condensed phase system without the ligand. A.) No permanent real
     * space electrostatics needs to be calculated because this was handled
     * analytically in step 1. B.) Permanent reciprocal space scaled by (1 -
     * lambda). C.) Polarization scaled by (1 - lambda).
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
        realSpaceTerm = false;
        permanentScale = 1.0 - sc2;
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
     * 3.) Ligand in vapor A.) Real space with an Ewald coefficient of 0.0 (no
     * reciprocal space). B.) Polarization scaled as in Step 2 by (1 - lambda).
     */
    private double vaporElec() {
        for (int i = 0; i < nAtoms; i++) {
            use[i] = atoms[i].applyLambda();
        }

        /**
         * Scale the permanent vacuum electrostatics. The softcore alpha is not
         * necessary (nothing in vacuum to collide with).
         */
        realSpaceTerm = true;
        permanentScale = 1.0 - sc2;
        dEdLSign = -1.0;
        double lAlphaBack = sc1;
        double dlAlphaBack = dsc1dL;
        double d2lAlphaBack = d2sc1dL2;
        sc1 = 0.0;
        dsc1dL = 0.0;
        d2sc1dL2 = 0.0;

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
        sc1 = lAlphaBack;
        dsc1dL = dlAlphaBack;
        d2sc1dL2 = d2lAlphaBack;
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
                if (esvTerm) {
                    reciprocalSpace.computePermanentDotPhi(cartMultipoleDotPhi);
                }
            }
        } catch (RuntimeException ex) {
            logger.warning("Runtime exception computing the permanent multipole field.");
            throw ex;
        } catch (Exception ex) {
            String msg = "Fatal exception computing the permanent multipole field.";
            logger.log(Level.SEVERE, msg, ex);
        }

        /**
         * Compute Born radii if necessary.
         */
        if (generalizedKirkwoodTerm) {
            if (esvTerm) {
                throw new UnsupportedOperationException();
            }
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
                        if (esvTerm) {
                            reciprocalSpace.splineInducedDipoles(unscaledInducedDipole, unscaledInducedDipoleCR, use);
                            sectionTeam.execute(inducedDipoleFieldRegion);
                            reciprocalSpace.computeInducedPhi(unscaledCartDipolePhi, unscaledCartDipolePhiCR);
                        }
                    } catch (RuntimeException ex) {
                        logger.warning("Runtime exception computing the induced reciprocal space field.");
                        throw ex;
                    } catch (Exception ex) {
                        String msg = "Fatal exception computing the induced reciprocal space field.";
                        logger.log(Level.SEVERE, msg, ex);
                    }
                } else {
                    reciprocalSpace.cartToFracInducedDipoles(inducedDipole, inducedDipoleCR);
                    if (esvTerm) {
                        reciprocalSpace.cartToFracUnscaledDipoles(unscaledInducedDipole, unscaledInducedDipoleCR);
                    }
                }
            }
            if (scfPredictor != null) {
                scfPredictor.saveMutualInducedDipoles(
                        inducedDipole, inducedDipoleCR, directDipole, directDipoleCR);
            }
        }

        /**
         * Find the total real space energy. This includes the permanent
         * multipoles in their own real space potential and the interaction of
         * permanent multipoles with induced dipoles.
         *
         * Then compute the permanent and reciprocal space energy.
         */
        if (reciprocalSpaceTerm && aewald > 0.0) {
            try {
                parallelTeam.execute(reciprocalEnergyRegion);
            } catch (RuntimeException ex) {
                logger.warning("Exception computing the reciprocal space energy.");
                throw ex;
            } catch (Exception ex) {
                logger.log(Level.SEVERE, "Exception computing the reciprocal space energy.", ex);
            }
            interactions += nAtoms;
            eself = reciprocalEnergyRegion.getPermanentSelfEnergy();
            erecip = reciprocalEnergyRegion.getPermanentReciprocalEnergy();
            eselfi = reciprocalEnergyRegion.getInducedDipoleSelfEnergy();
            erecipi = reciprocalEnergyRegion.getInducedDipoleReciprocalEnergy();
        }

        if (realSpaceTerm) {
            realSpaceTimeTotal = -System.nanoTime();
            try {
                parallelTeam.execute(realSpaceEnergyRegion);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the real space energy.");
                throw ex;
                } catch (Exception ex) {
                    String msg = "Fatal exception computing the real space energy";
                    logger.log(Level.SEVERE, msg, ex);
                }
            ereal = realSpaceEnergyRegion.getPermanentEnergy();
            ereali = realSpaceEnergyRegion.getPolarizationEnergy();
            interactions += realSpaceEnergyRegion.getInteractions();
            realSpaceTimeTotal += System.nanoTime();
        }

        /**
         * Compute the generalized Kirkwood solvation free energy.
         */
        if (generalizedKirkwoodTerm) {
            gkEnergyTotal -= System.nanoTime();
            generalizedKirkwoodEnergy += generalizedKirkwood.solvationEnergy(gradient, print);
            gkInteractions += generalizedKirkwood.getInteractions();
            if (esvTerm) {
                // TODO: Get GK derivative.
            }
            gkEnergyTotal += System.nanoTime();
        }

        /**
         * Collect energy terms.
         */
        permanentRealSpaceEnergy = ereal;
        inducedRealSpaceEnergy = ereali;
        permanentSelfEnergy = eself;
        inducedSelfEnergy = eselfi;
        permanentReciprocalEnergy = erecip;
        inducedReciprocalEnergy = erecipi;
        permanentMultipoleEnergy += eself + erecip + ereal;
        polarizationEnergy += eselfi + erecipi + ereali;

        /**
         * Log some info.
         */
        if (logger.isLoggable(Level.FINE) || printDecomposition) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n   Multipole Self-Energy:   %16.8f\n", permanentSelfEnergy));
            sb.append(format("   Multipole Reciprocal:    %16.8f\n", permanentReciprocalEnergy));
            sb.append(format("   Multipole Real Space:    %16.8f\n", permanentRealSpaceEnergy));
            sb.append(format("   Polarization Self-Energy:%16.8f\n", inducedSelfEnergy));
            sb.append(format("   Polarization Reciprocal: %16.8f\n", inducedReciprocalEnergy));
            sb.append(format("   Polarization Real Space: %16.8f\n", inducedRealSpaceEnergy));
            if (generalizedKirkwoodTerm) {
                sb.append(format("   Generalized Kirkwood:    %16.8f\n", generalizedKirkwoodEnergy));
            }
            logger.info(sb.toString());
        }

        return permanentMultipoleEnergy + polarizationEnergy + generalizedKirkwoodEnergy;
    }

    @Override
    public int getInteractions() {
        return interactions;
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
        } catch (RuntimeException ex) {
            logger.warning("Runtime exception computing direct induced dipoles.");
            throw ex;
        } catch (Exception ex) {
            String msg = "Fatal exception computing the direct induced dipoles";
            logger.log(Level.SEVERE, msg, ex);
        }

        /**
         * Return unless mutual polarization is selected.
         */
        if (polarization != Polarization.MUTUAL) {
            if (nSymm > 1) {
                try {
                    parallelTeam.execute(expandInducedDipolesRegion);
                } catch (RuntimeException ex) {
                    logger.warning("Exception expanding direct induced dipoles.");
                    throw ex;
                } catch (Exception ex) {
                    logger.log(Level.SEVERE, "Exception expanding direct induced dipoles.", ex);
                }
            }
            return 0;
        }

        /**
         * Predict the current self-consistent induced dipoles using information
         * from previous steps.
         */
        if (scfPredictor != null) {
            if (esvTerm) {
                throw new UnsupportedOperationException();
            }
            scfPredictor.run(lambdaMode);
        }

        /**
         * Expand the initial induced dipoles to P1 symmetry, if necessary.
         */
        if (nSymm > 1) {
            try {
                parallelTeam.execute(expandInducedDipolesRegion);
            } catch (RuntimeException ex) {
                logger.warning("Exception expanding initial induced dipoles.");
                throw ex;
            } catch (Exception ex) {
                logger.log(Level.SEVERE, "Exception expanding initial induced dipoles.", ex);
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

                if (esvTerm && reciprocalSpaceTerm && aewald > 0.0) {
                    reciprocalSpace.splineInducedDipoles(unscaledInducedDipole, unscaledInducedDipoleCR, use);
                    sectionTeam.execute(inducedDipoleFieldRegion);
                    reciprocalSpace.computeInducedPhi(unscaledCartDipolePhi, unscaledCartDipolePhiCR);
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
            } catch (RuntimeException ex) {
                logger.warning("Exception computing mutual induced dipoles.");
                throw ex;
            } catch (Exception ex) {
                logger.log(Level.SEVERE, "Exception computing mutual induced dipoles.", ex);
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
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the permanent multipole field.");
                throw ex;
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
                    realSpacePermTimeTotal -= System.nanoTime();
                    parallelTeam.execute(permanentRealSpaceFieldRegion);
                    realSpacePermTimeTotal += System.nanoTime();
                } catch (RuntimeException ex) {
                    logger.warning("Runtime exception computing the real space field.");
                    throw ex;
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
                } catch (RuntimeException ex) {
                    logger.warning("Runtime exception computing the real space field in thread " + getThreadIndex());
                    throw ex;
                } catch (Exception ex) {
                    String msg = "Fatal exception computing the real space field in thread " + getThreadIndex();
                    logger.log(Level.SEVERE, msg, ex);
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
                    realSpacePermTimes[threadIndex] -= System.nanoTime();
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
                    realSpacePermTimes[threadIndex] += System.nanoTime();
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
                private double fXuns[], fYuns[], fZuns[];
                private double fXCRuns[], fYCRuns[], fZCRuns[];
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
                    realSpacePermTimes[threadIndex] -= System.nanoTime();
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
                    if (esvTerm) {
                        fXuns = unscaledField[threadIndex][0];
                        fYuns = unscaledField[threadIndex][1];
                        fZuns = unscaledField[threadIndex][2];
                        fXCRuns = unscaledFieldCR[threadIndex][0];
                        fYCRuns = unscaledFieldCR[threadIndex][1];
                        fZCRuns = unscaledFieldCR[threadIndex][2];
                    }
                }

                @Override
                public void finish() {
                    int threadIndex = getThreadIndex();
                    sharedCount.addAndGet(count);
                    realSpacePermTimes[threadIndex] += System.nanoTime();
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
                        final Atom ai = atoms[i];
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.getIndex() - 1;
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
                                int index = ak.getIndex() - 1;
                                maskp_local[index] = p13scale;
                            }
                        }
                        for (Bond bond : ai.getBonds()) {
                            int index = bond.get1_2(ai).getIndex() - 1;
                            maskp_local[index] = p12scale;
                        }
                        /**
                         * Apply group based polarization masking rule.
                         */
                        for (int index : ip11[i]) {
                            mask_local[index] = d11scale;
                        }
//                        applyMaskingRules(true, i, null, mask_local, maskp_local);
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
                                    if (esvTerm) {
                                        throw new IllegalStateException();
                                    }
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
                                final double scale = (!esvTerm || esvConfig.allowMaskPolarD) ? mask_local[k] : 1.0;
                                final double scalep = (!esvTerm || esvConfig.allowMaskPolarP) ? maskp_local[k] : 1.0;
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
                                int index = ak.getIndex() - 1;
                                if (index < 0 || index >= nAtoms) {
                                    ak.print();
                                }
                                maskp_local[index] = 1.0;
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.getIndex() - 1;
                                maskp_local[index] = 1.0;
                            }
                        }
                        for (Bond bond : ai.getBonds()) {
                            int index = bond.get1_2(ai).getIndex() - 1;
                            maskp_local[index] = 1.0;
                        }
                        for (int index : ip11[i]) {
                            mask_local[index] = 1.0;
                        }
//                        applyMaskingRules(false, i, null, mask_local, maskp_local);
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
                if (reciprocalSpaceTerm && aewald > 0.0 && (!esvTerm || esvConfig.recipFieldEffects)) {
                    execute(inducedRealSpaceFieldSection, inducedReciprocalFieldSection);
                } else {
                    execute(inducedRealSpaceFieldSection);
                }
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the induced dipole field.");
                throw ex;
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
                    realSpaceScfTimeTotal -= System.nanoTime();
                    pt.execute(polarizationRealSpaceFieldRegion);
                    realSpaceScfTimeTotal += System.nanoTime();
                } catch (RuntimeException ex) {
                    logger.warning("Runtime exception computing the real space SCF in thread " + getThreadIndex());
                    throw ex;
                } catch (Exception ex) {
                    String msg = "Fatal exception computing the real space SCF in thread " + getThreadIndex();
                    logger.log(Level.SEVERE, msg, ex);
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
                } catch (RuntimeException ex) {
                    logger.warning("Runtime exception computing the induced real space field in thread " + getThreadIndex());
                    throw ex;
                } catch (Exception ex) {
                    String msg = "Fatal exception computing the induced real space field in thread " + getThreadIndex();
                    logger.log(Level.SEVERE, msg, ex);
                }
            }

            private class InducedRealSpaceFieldLoop extends IntegerForLoop {

                private double ind[][], indCR[][];
                private double uns[][], unsCR[][];
                private double x[], y[], z[];
                private double fX[], fY[], fZ[];
                private double fXCR[], fYCR[], fZCR[];
                private double[] fXuns, fYuns, fZuns;
                private double[] fXCRuns, fYCRuns, fZCRuns;

                public InducedRealSpaceFieldLoop() {
                }

                @Override
                public IntegerSchedule schedule() {
                    return realSpaceSchedule;
                }

                @Override
                public void start() {
                    int threadIndex = getThreadIndex();
                    realSpaceScfTimes[threadIndex] -= System.nanoTime();
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
                    if (esvTerm) {
                        uns = unscaledInducedDipole[0];
                        unsCR = unscaledInducedDipoleCR[0];
                        fXuns = unscaledField[threadIndex][0];
                        fYuns = unscaledField[threadIndex][1];
                        fZuns = unscaledField[threadIndex][2];
                        fXCRuns = unscaledFieldCR[threadIndex][0];
                        fYCRuns = unscaledFieldCR[threadIndex][1];
                        fZCRuns = unscaledFieldCR[threadIndex][2];
                        fill(fXuns, 0.0);
                        fill(fYuns, 0.0);
                        fill(fZuns, 0.0);
                        fill(fXCRuns, 0.0);
                        fill(fYCRuns, 0.0);
                        fill(fZCRuns, 0.0);
                    }
                }

                @Override
                public void finish() {
                    int threadIndex = getThreadIndex();
                    realSpaceScfTimes[threadIndex] += System.nanoTime();
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
                        double fxu = 0.0, fyu = 0.0, fzu = 0.0;
                        double pxu = 0.0, pyu = 0.0, pzu = 0.0;
                        final double uixu, uiyu, uizu;
                        final double pixu, piyu, pizu;
                        if (esvTerm) {
                            final double[] unsi = uns[i];
                            uixu = unsi[0];
                            uiyu = unsi[1];
                            uizu = unsi[2];
                            final double[] unsiCR = unsCR[i];
                            pixu = unsiCR[0];
                            piyu = unsiCR[1];
                            pizu = unsiCR[2];
                        } else {
                            uixu = uix;
                            uiyu = uiy;
                            uizu = uiz;
                            pixu = pix;
                            piyu = piy;
                            pizu = piz;
                        }
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
                                if (esvTerm) {
                                    throw new UnsupportedOperationException();
                                }
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
                            {	// scaled
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
                            }	// scaled
                            if (esvTerm) {	// unscaled
                                final double unscaledk[] = uns[k];
                                final double ukxu = unscaledk[0];
                                final double ukyu = unscaledk[1];
                                final double ukzu = unscaledk[2];
                                final double ukru = ukxu * xr + ukyu * yr + ukzu * zr;
                                final double bn2ukr = bn2 * ukru;
                                final double fimx = -bn1 * ukxu + bn2ukr * xr;
                                final double fimy = -bn1 * ukyu + bn2ukr * yr;
                                final double fimz = -bn1 * ukzu + bn2ukr * zr;
                                final double rr5ukr = rr5 * ukru;
                                final double fidx = -rr3 * ukxu + rr5ukr * xr;
                                final double fidy = -rr3 * ukyu + rr5ukr * yr;
                                final double fidz = -rr3 * ukzu + rr5ukr * zr;
                                fxu += (fimx - fidx);
                                fyu += (fimy - fidy);
                                fzu += (fimz - fidz);
                                final double unscaledPk[] = unsCR[k];
                                final double pkxu = unscaledPk[0];
                                final double pkyu = unscaledPk[1];
                                final double pkzu = unscaledPk[2];
                                final double pkru = pkxu * xr + pkyu * yr + pkzu * zr;
                                final double bn2pkr = bn2 * pkru;
                                final double pimx = -bn1 * pkxu + bn2pkr * xr;
                                final double pimy = -bn1 * pkyu + bn2pkr * yr;
                                final double pimz = -bn1 * pkzu + bn2pkr * zr;
                                final double rr5pkr = rr5 * pkru;
                                final double pidx = -rr3 * pkxu + rr5pkr * xr;
                                final double pidy = -rr3 * pkyu + rr5pkr * yr;
                                final double pidz = -rr3 * pkzu + rr5pkr * zr;
                                pxu += (pimx - pidx);
                                pyu += (pimy - pidy);
                                pzu += (pimz - pidz);
                                final double uixuxr = uixu * xr;
                                final double uiyuyr = uiyu * yr;
                                final double uizuzr = uizu * zr;
                                final double uiru = uixuxr + uiyuyr + uizuzr;
                                final double bn2uir = bn2 * uiru;
                                final double fkmx = -bn1 * uixuxr + bn2uir;
                                final double fkmy = -bn1 * uiyuyr + bn2uir;
                                final double fkmz = -bn1 * uizuzr + bn2uir;
                                final double rr5uir = rr5 * uiru;
                                final double fkdx = -rr3 * uixuxr + rr5uir;
                                final double fkdy = -rr3 * uiyuyr + rr5uir;
                                final double fkdz = -rr3 * uizuzr + rr5uir;
                                fXuns[k] += (fkmx - fkdx);
                                fYuns[k] += (fkmy - fkdy);
                                fZuns[k] += (fkmz - fkdz);
                                final double pixuxr = pixu * xr;
                                final double piyuyr = piyu * yr;
                                final double pizuzr = pizu * zr;
                                final double piru = pixuxr + piyuyr + pizuzr;
                                final double bn2pir = bn2 * piru;
                                final double pkmx = -bn1 * pixuxr + bn2pir;
                                final double pkmy = -bn1 * piyuyr + bn2pir;
                                final double pkmz = -bn1 * pizuzr + bn2pir;
                                final double rr5pir = rr5 * piru;
                                final double pkdx = -rr3 * pixuxr + rr5pir;
                                final double pkdy = -rr3 * piyuyr + rr5pir;
                                final double pkdz = -rr3 * pizuzr + rr5pir;
                                fXCRuns[k] += (pkmx - pkdx);
                                fYCRuns[k] += (pkmy - pkdy);
                                fZCRuns[k] += (pkmz - pkdz);
                            }
                        }
                        fX[i] += fx;
                        fY[i] += fy;
                        fZ[i] += fz;
                        fXCR[i] += px;
                        fYCR[i] += py;
                        fZCR[i] += pz;
                        if (esvTerm) {
                            fXuns[i] += fxu;
                            fYuns[i] += fyu;
                            fZuns[i] += fzu;
                            fXCRuns[i] += pxu;
                            fYCRuns[i] += pyu;
                            fZCRuns[i] += pzu;
                        }
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
                        final double unsym[][] = (esvTerm) ? unscaledInducedDipole[iSymm] : null;
                        final double unsymCR[][] = (esvTerm) ? unscaledInducedDipoleCR[iSymm] : null;
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
                            double fxu = 0.0;
                            double fyu = 0.0;
                            double fzu = 0.0;
                            double pxu = 0.0;
                            double pyu = 0.0;
                            double pzu = 0.0;
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
                            final double uixu, uiyu, uizu;
                            final double pixu, piyu, pizu;
                            if (esvTerm) {
                                final double[] unscaledi = uns[i];
                                uixu = unscaledi[0];
                                uiyu = unscaledi[1];
                                uizu = unscaledi[2];
                                final double[] unscaledCRi = unsCR[i];
                                pixu = unscaledCRi[0];
                                piyu = unscaledCRi[1];
                                pizu = unscaledCRi[2];
                            } else {
                                uixu = uix;
                                uiyu = uiy;
                                uizu = uiz;
                                pixu = pix;
                                piyu = piy;
                                pizu = piz;
                            }
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
                                {	// scaled
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
                                }	// scaled
                                if (esvTerm) {	// unscaled
                                    final double[] unsymk = unsym[k];
                                    final double ukxu = unsymk[0];
                                    final double ukyu = unsymk[1];
                                    final double ukzu = unsymk[2];
                                    final double[] unsymCRk = unsymCR[k];
                                    final double pkxu = unsymCRk[0];
                                    final double pkyu = unsymCRk[1];
                                    final double pkzu = unsymCRk[2];
                                    final double ukru = ukxu * xr + ukyu * yr + ukzu * zr;
                                    final double bn2ukr = bn2 * ukru;
                                    final double fimx = -bn1 * ukxu + bn2ukr * xr;
                                    final double fimy = -bn1 * ukyu + bn2ukr * yr;
                                    final double fimz = -bn1 * ukzu + bn2ukr * zr;
                                    final double rr5ukr = rr5 * ukru;
                                    final double fidx = -rr3 * ukxu + rr5ukr * xr;
                                    final double fidy = -rr3 * ukyu + rr5ukr * yr;
                                    final double fidz = -rr3 * ukzu + rr5ukr * zr;
                                    fxu += selfScale * (fimx - fidx);
                                    fyu += selfScale * (fimy - fidy);
                                    fzu += selfScale * (fimz - fidz);
                                    final double pkru = pkxu * xr + pkyu * yr + pkzu * zr;
                                    final double bn2pkr = bn2 * pkru;
                                    final double pimx = -bn1 * pkxu + bn2pkr * xr;
                                    final double pimy = -bn1 * pkyu + bn2pkr * yr;
                                    final double pimz = -bn1 * pkzu + bn2pkr * zr;
                                    final double rr5pkr = rr5 * pkru;
                                    final double pidx = -rr3 * pkxu + rr5pkr * xr;
                                    final double pidy = -rr3 * pkyu + rr5pkr * yr;
                                    final double pidz = -rr3 * pkzu + rr5pkr * zr;
                                    pxu += selfScale * (pimx - pidx);
                                    pyu += selfScale * (pimy - pidy);
                                    pzu += selfScale * (pimz - pidz);
                                    final double uir = uixu * xr + uiyu * yr + uizu * zr;
                                    final double bn2uir = bn2 * uir;
                                    final double fkmx = -bn1 * uixu + bn2uir * xr;
                                    final double fkmy = -bn1 * uiyu + bn2uir * yr;
                                    final double fkmz = -bn1 * uizu + bn2uir * zr;
                                    final double rr5uir = rr5 * uir;
                                    final double fkdx = -rr3 * uixu + rr5uir * xr;
                                    final double fkdy = -rr3 * uiyu + rr5uir * yr;
                                    final double fkdz = -rr3 * uizu + rr5uir * zr;
                                    final double xcu = selfScale * (fkmx - fkdx);
                                    final double ycu = selfScale * (fkmy - fkdy);
                                    final double zcu = selfScale * (fkmz - fkdz);
                                    fXuns[k] += (xcu * transOp[0][0] + ycu * transOp[1][0] + zcu * transOp[2][0]);
                                    fYuns[k] += (xcu * transOp[0][1] + ycu * transOp[1][1] + zcu * transOp[2][1]);
                                    fZuns[k] += (xcu * transOp[0][2] + ycu * transOp[1][2] + zcu * transOp[2][2]);
                                    final double pir = pixu * xr + piyu * yr + pizu * zr;
                                    final double bn2pir = bn2 * pir;
                                    final double pkmx = -bn1 * pixu + bn2pir * xr;
                                    final double pkmy = -bn1 * piyu + bn2pir * yr;
                                    final double pkmz = -bn1 * pizu + bn2pir * zr;
                                    final double rr5pir = rr5 * pir;
                                    final double pkdx = -rr3 * pixu + rr5pir * xr;
                                    final double pkdy = -rr3 * piyu + rr5pir * yr;
                                    final double pkdz = -rr3 * pizu + rr5pir * zr;
                                    final double xcp = selfScale * (pkmx - pkdx);
                                    final double ycp = selfScale * (pkmy - pkdy);
                                    final double zcp = selfScale * (pkmz - pkdz);
                                    fXCRuns[k] += (xcp * transOp[0][0] + ycp * transOp[1][0] + zcp * transOp[2][0]);
                                    fYCRuns[k] += (xcp * transOp[0][1] + ycp * transOp[1][1] + zcp * transOp[2][1]);
                                    fZCRuns[k] += (xcp * transOp[0][2] + ycp * transOp[1][2] + zcp * transOp[2][2]);
                                }	// unscaled
                            }
                            fX[i] += fx;
                            fY[i] += fy;
                            fZ[i] += fz;
                            fXCR[i] += px;
                            fYCR[i] += py;
                            fZCR[i] += pz;
                            if (esvTerm) {
                                fXuns[i] += fxu;
                                fYuns[i] += fyu;
                                fZuns[i] += fzu;
                                fXCRuns[i] += pxu;
                                fYCRuns[i] += pyu;
                                fZCRuns[i] += pzu;
                            }
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
        public void run() {
            int ti = getThreadIndex();
            if (directLoop[ti] == null) {
                directLoop[ti] = new DirectLoop();
            }
            try {
                execute(0, nAtoms - 1, directLoop[ti]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the direct induced dipoles in thread " + getThreadIndex());
                throw ex;
            } catch (Exception ex) {
                String msg = "Fatal exception computing the direct induced dipoles in thread " + getThreadIndex();
                logger.log(Level.SEVERE, msg, ex);
            }
        }

        @Override
        public void finish() {
            // Print the scaled permanent multipoles, reduced field, and direct dipoles for debug atoms.
            if (esvTerm && esvConfig.verbose) {
                for (int k = 0; k < nAtoms; k++) {
                    if ((esvAtomsScaled[k] || esvAtomsScaledAlpha[k])
                            && (esvConfig.permanentWhitelist(k) || !esvConfig.useWhitelists)) {
                        double[] mpole = globalMultipole[0][k];
                        double[] mdot = esvMultipoleDot[0][k];
                        logger.info(format("(%d) mpole_sym0: %.6g; %.6g %.6g %.6g; %.6g %.6g %.6g %.6g %.6g %.6g",
                                k, mpole[0], mpole[1], mpole[2], mpole[3],
                                mpole[4], mpole[5], mpole[6], mpole[7], mpole[8], mpole[9]));
                        logger.info(format("(%d) esv_sym0:   %.6g; %.6g %.6g %.6g; %.6g %.6g %.6g %.6g %.6g %.6g",
                                k, mdot[0], mdot[1], mdot[2], mdot[3],
                                mdot[4], mdot[5], mdot[6], mdot[7], mdot[8], mdot[9]));
                    }
                    if ((esvAtomsScaled[k] || esvAtomsScaledAlpha[k])
                            && (esvConfig.inducedWhitelist(k) || !esvConfig.useWhitelists)) {
                        logger.info(format("(%d) polar,field,direct: %.4g * (%.4g %.4g %.4g) = (%.4g %.4g %.4g)",
                                k, polarizability[k],
                                field[0][0][k], field[0][1][k], field[0][2][k],
                                directDipole[k][0], directDipole[k][1], directDipole[k][2]));
                        logger.info(format("(%d) polUns,fUns,dirUns: %.4g * (%.4g %.4g %.4g) = (%.4g %.4g %.4g)",
                                k, unscaledPolarizability[k],
                                unscaledField[0][0][k], unscaledField[0][1][k], unscaledField[0][2][k],
                                unscaledDirectDipole[k][0], unscaledDirectDipole[k][1], unscaledDirectDipole[k][2]));
                    }
                }
            }
        }

        private class DirectLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int lb, int ub) {
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
                    if (esvTerm) {
                        double fxu = 0.0;
                        double fyu = 0.0;
                        double fzu = 0.0;
                        double fxCRu = 0.0;
                        double fyCRu = 0.0;
                        double fzCRu = 0.0;
                        for (int j = 1; j < maxThreads; j++) {
                            fxu += unscaledField[j][0][i];
                            fyu += unscaledField[j][1][i];
                            fzu += unscaledField[j][2][i];
                            fxCRu += unscaledFieldCR[j][0][i];
                            fyCRu += unscaledFieldCR[j][1][i];
                            fzCRu += unscaledFieldCR[j][2][i];
                        }
                        unscaledField[0][0][i] += fxu;
                        unscaledField[0][1][i] += fyu;
                        unscaledField[0][2][i] += fzu;
                        unscaledFieldCR[0][0][i] += fxCRu;
                        unscaledFieldCR[0][1][i] += fyCRu;
                        unscaledFieldCR[0][2][i] += fzCRu;
                    }
                }
                if (aewald > 0.0 && (!esvTerm || esvConfig.recipFieldEffects)) {
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
                        if (esvTerm) {
                            unscaledField[0][0][i] += fx;
                            unscaledField[0][1][i] += fy;
                            unscaledField[0][2][i] += fz;
                            unscaledFieldCR[0][0][i] += fx;
                            unscaledFieldCR[0][1][i] += fy;
                            unscaledFieldCR[0][2][i] += fz;
                        }
                    }
                }
                if (generalizedKirkwoodTerm) {
                    if (esvTerm) {
                        throw new UnsupportedOperationException();
                    }
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
                    if (esvTerm) {
                        // Zero out dipoles not on the induced whitelist.
                        if (esvConfig.useWhitelists && !esvConfig.inducedWhitelist(i)) {
                            for (int axis = 0; axis < 3; axis++) {
                                directi[axis] = 0.0;
                                directCRi[axis] = 0.0;
                                ind[axis] = 0.0;
                                indCR[axis] = 0.0;
                            }
                            continue;
                        }
                        // Unscaled direct dipoles. Permanent field
                        final double polarUns = unscaledPolarizability[i];
                        final double[] unsDiri = unscaledDirectDipole[i];
                        final double[] unsDirCRi = unscaledDirectDipoleCR[i];
                        final double[] unsIndi = unscaledInducedDipole[0][i];
                        final double[] unsIndCRi = unscaledInducedDipoleCR[0][i];
                        unsDiri[0] = polarUns * field[0][0][i];
                        unsDiri[1] = polarUns * field[0][1][i];
                        unsDiri[2] = polarUns * field[0][2][i];
                        unsDirCRi[0] = polarUns * fieldCR[0][0][i];
                        unsDirCRi[1] = polarUns * fieldCR[0][1][i];
                        unsDirCRi[2] = polarUns * fieldCR[0][2][i];
                        unsIndi[0] = unsDiri[0];
                        unsIndi[1] = unsDiri[1];
                        unsIndi[2] = unsDiri[2];
                        unsIndCRi[0] = unsDirCRi[0];
                        unsIndCRi[1] = unsDirCRi[1];
                        unsIndCRi[2] = unsDirCRi[2];
                    }
                }
            }
        }
    }

    private class SORRegion extends ParallelRegion {

        private final SORLoop sorLoop[];
        private final SharedDouble sharedEps;
        private final SharedDouble sharedEpsCR;
        private SharedDouble sharedEpsUns;
        private SharedDouble sharedEpsUnsCR;

        public SORRegion(int nt) {
            sorLoop = new SORLoop[nt];
            sharedEps = new SharedDouble();
            sharedEpsCR = new SharedDouble();
            sharedEpsUns = (esvTerm) ? new SharedDouble() : null;
            sharedEpsUnsCR = (esvTerm) ? new SharedDouble() : null;
        }

        public double getEps() {
            double eps = sharedEps.get();
            double epsCR = sharedEpsCR.get();
            if (esvTerm) {
                double epsUns = sharedEpsUns.get();
                double epsUnsCR = sharedEpsUnsCR.get();
                return max(max(eps, epsCR), max(epsUns, epsUnsCR));
            } else {
                return max(eps, epsCR);
            }
        }

        @Override
        public void start() {
            sharedEps.set(0.0);
            sharedEpsCR.set(0.0);
            if (esvTerm) {
                if (sharedEpsUns == null) {
                    sharedEpsUns = new SharedDouble();
                    sharedEpsUnsCR = new SharedDouble();
                }
                sharedEpsUns.set(0.0);
                sharedEpsUnsCR.set(0.0);
            }
        }

        @Override
        public void run() {
            try {
                int ti = getThreadIndex();
                if (sorLoop[ti] == null) {
                    sorLoop[ti] = new SORLoop();
                }
                execute(0, nAtoms - 1, sorLoop[ti]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the mutual induced dipoles in thread " + getThreadIndex());
                throw ex;
            } catch (Exception ex) {
                String msg = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex();
                logger.log(Level.SEVERE, msg, ex);
            }

        }

        private class SORLoop extends IntegerForLoop {

            private double eps, epsCR;
            private double epsUns, epsUnsCR;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eps = 0.0;
                epsCR = 0.0;
                epsUns = 0.0;
                epsUnsCR = 0.0;
            }

            @Override
            public void run(int lb, int ub) {
                final double induced0[][] = inducedDipole[0];
                final double inducedCR0[][] = inducedDipoleCR[0];
                final double unscaled0[][] = (esvTerm) ? unscaledInducedDipole[0] : null;
                final double unscaledCR0[][] = (esvTerm) ? unscaledInducedDipoleCR[0] : null;
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
                    double fxu = 0.0;
                    double fyu = 0.0;
                    double fzu = 0.0;
                    double fxCRu = 0.0;
                    double fyCRu = 0.0;
                    double fzCRu = 0.0;
                    for (int j = 1; j < maxThreads; j++) {
                        fx += field[j][0][i];
                        fy += field[j][1][i];
                        fz += field[j][2][i];
                        fxCR += fieldCR[j][0][i];
                        fyCR += fieldCR[j][1][i];
                        fzCR += fieldCR[j][2][i];
                        if (esvTerm) {
                            fxu += unscaledField[j][0][i];
                            fyu += unscaledField[j][1][i];
                            fzu += unscaledField[j][2][i];
                            fxCRu += unscaledFieldCR[j][0][i];
                            fyCRu += unscaledFieldCR[j][1][i];
                            fzCRu += unscaledFieldCR[j][2][i];
                        }
                    }
                    field[0][0][i] += fx;
                    field[0][1][i] += fy;
                    field[0][2][i] += fz;
                    fieldCR[0][0][i] += fxCR;
                    fieldCR[0][1][i] += fyCR;
                    fieldCR[0][2][i] += fzCR;
                    if (esvTerm) {
                        unscaledField[0][0][i] += fxu;
                        unscaledField[0][1][i] += fyu;
                        unscaledField[0][2][i] += fzu;
                        unscaledFieldCR[0][0][i] += fxCRu;
                        unscaledFieldCR[0][1][i] += fyCRu;
                        unscaledFieldCR[0][2][i] += fzCRu;
                    }
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
                        if (esvTerm) {
                            final double dipoleiUns[] = unscaled0[i];
                            final double dipoleiUnsCR[] = unscaledCR0[i];
                            final double phiiUns[] = unscaledCartDipolePhi[i];
                            final double phiiUnsCR[] = unscaledCartDipolePhiCR[i];
                            unscaledField[0][0][i] += (aewald3 * dipoleiUns[0] - phiiUns[t100]);
                            unscaledField[0][1][i] += (aewald3 * dipoleiUns[1] - phiiUns[t010]);
                            unscaledField[0][2][i] += (aewald3 * dipoleiUns[2] - phiiUns[t001]);
                            unscaledFieldCR[0][0][i] += (aewald3 * dipoleiUnsCR[0] - phiiUnsCR[t100]);
                            unscaledFieldCR[0][1][i] += (aewald3 * dipoleiUnsCR[1] - phiiUnsCR[t010]);
                            unscaledFieldCR[0][2][i] += (aewald3 * dipoleiUnsCR[2] - phiiUnsCR[t001]);
                        }
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
                        fieldCR[0][0][i] += gkFieldCR[0].get(i);
                        field[0][1][i] += gkField[1].get(i);
                        fieldCR[0][1][i] += gkFieldCR[1].get(i);
                        field[0][2][i] += gkField[2].get(i);
                        fieldCR[0][2][i] += gkFieldCR[2].get(i);
                    }
                    if (esvTerm) {
                        throw new UnsupportedOperationException();	// TODO: Handle unscaled field.
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
                        if (esvTerm) {
                            final double uns[] = unscaled0[i];
                            final double unsCR[] = unscaledCR0[i];
                            final double directUns[] = unscaledDirectDipole[i];
                            final double directUnsCR[] = unscaledDirectDipoleCR[i];
                            final double polarUns = unscaledPolarizability[i];
                            double prevUns = uns[j];
                            double mutualUns = polarUns * unscaledField[0][j][i];
                            uns[j] = directUns[j] + mutualUns;
                            double deltaUns = polsor * (uns[j] - prevUns);
                            uns[j] = prevUns + deltaUns;
                            epsUns += deltaUns * deltaUns;
                            prevUns = unsCR[j];
                            mutualUns = polarUns * unscaledFieldCR[0][j][i];
                            unsCR[j] = directUnsCR[j] + mutualUns;
                            deltaUns = polsor * (unsCR[j] - prevUns);
                            unsCR[j] = prevUns + deltaUns;
                            epsUnsCR += deltaUns * deltaUns;
                        }
                    }
                }
            }

            @Override
            public void finish() {
                sharedEps.addAndGet(eps);
                sharedEpsCR.addAndGet(epsCR);
                if (esvTerm) {
                    sharedEpsUns.addAndGet(epsUns);
                    sharedEpsUnsCR.addAndGet(epsUnsCR);
                }
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
        //private final double mutualScale = (polarization == Polarization.DIRECT || polarization == Polarization.NONE) ? 0.0 : 1.0;
        private final int numThreads;

        private final SharedInteger sharedInteractions;
        private final RealSpaceEnergyLoop realSpaceEnergyLoops[];

        public RealSpaceEnergyRegion(int nt) {
            numThreads = nt;
            sharedInteractions = new SharedInteger();
            realSpaceEnergyLoops = new RealSpaceEnergyLoop[nt];
            // One-time, single-threaded instantiation of the loop and tensor arrays.
            for (int thread = 0; thread < nt; thread++) {
                realSpaceEnergyLoops[thread] = new RealSpaceEnergyLoop();
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
            // Zero out grad,torque, and field arrays.
            // [threadID][X/Y/Z][atomID]
            for (int thread = 0; thread < maxThreads; thread++) {
                for (int i = 0; i < 3; i++) {
                    fill(grad[thread][i], 0.0);
                    fill(torque[thread][i], 0.0);
//                    fill(gradPerm[thread][i], 0.0);
//                    fill(gradInduced[thread][i], 0.0);
//                    fill(torquePerm[thread][i], 0.0);
//                    fill(torqueInduced[thread][i], 0.0);
                    fill(field[thread][i], 0.0);
                    fill(fieldCR[thread][i], 0.0);
                    if (lambdaTerm) {
                        fill(lambdaGrad[thread][i], 0.0);
                        fill(lambdaTorque[thread][i], 0.0);
                    }
                    if (esvTerm) {
                        fill(unscaledField[thread][i], 0.0);
                        fill(unscaledFieldCR[thread][i], 0.0);
                    }
                }
            }
            if (esvTerm) {
                for (int i = 0; i < numESVs; i++) {
                    esvPermRealDeriv_shared[i].getAndSet(0.0);
                    esvInducedRealDeriv_shared[i].getAndSet(0.0);
                }
            }
        }

        @Override
        public void run() {
            if (skipRealRegion) {
                logger.warning("Skipping PME real region.");
                return;
            }
            int threadIndex = getThreadIndex();
            realSpaceEnergyLoops[threadIndex] = new RealSpaceEnergyLoop();
            try {
                execute(0, nAtoms - 1, realSpaceEnergyLoops[threadIndex]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the real space energy in thread " + getThreadIndex());
                throw ex;
            } catch (Exception ex) {
                String msg = "Fatal exception computing the real space energy in thread " + getThreadIndex();
                logger.log(Level.SEVERE, msg, ex);
            }
        }

        @Override
        public void finish() {
            permanentEnergy = 0.0;
            polarizationEnergy = 0.0;
            for (int i = 0; i < maxThreads; i++) {
                double e = realSpaceEnergyLoops[i].permanentEnergy;
                if (Double.isNaN(e)) {
                    logger.warning(format(" The permanent energy of thread %d is %16.8f", i, e));
                    throw new EnergyException(format(" The permanent multipole energy of thread %d is %16.8f", i, e), true);
                }
                permanentEnergy += e;
                double ei = realSpaceEnergyLoops[i].inducedEnergy;
                if (Double.isNaN(ei)) {
                    logger.warning(format(" The induced energy of thread %d is %16.8f", i, ei));
                    throw new EnergyException(format(" The polarization energy of thread %d is %16.8f", i, ei), true);
                }
                polarizationEnergy += ei;
            }
        }
        
        private void rotateIntoAxisArrays(double[] source, int i,
                double[] destx, double[] desty, double[] destz, double prefactor) {
            destx[i] += prefactor * source[0];
            desty[i] += prefactor * source[1];
            destz[i] += prefactor * source[2];
        }
        
        private void transposeForcesToGradByAxis(
                double[] Fi, double[] Ti, double[] Tk,
                int i, int k, double prefactor,
                double[] gX, double[] gY, double[] gZ,
                double[] tX, double[] tY, double[] tZ,
                double[] gxk_local, double[] gyk_local, double[] gzk_local,
                double[] txk_local, double[] tyk_local, double[] tzk_local) {
            rotateIntoAxisArrays(Fi, i, gX, gY, gZ, prefactor);
            rotateIntoAxisArrays(Ti, i, tX, tY, tZ, prefactor);
            rotateIntoAxisArrays(Fi, k, gxk_local, gyk_local, gzk_local, -prefactor);
            rotateIntoAxisArrays(Tk, k, txk_local, tyk_local, tzk_local, prefactor);
        }
        
        /**
         * The Real Space Gradient Loop class contains methods and thread local
         * variables to parallelize the evaluation of the real space permanent
         * and polarization energies and gradients.
         */
        private class RealSpaceEnergyLoop extends IntegerForLoop {

            private boolean soft;
            private double permanentEnergy;
            private double inducedEnergy;
            private double dUdL, d2UdL2;
            private int iSymm, count;
            private SymOp symOp;
            // Store contributions to the gradient.
            private double[] gXm, gYm, gZm, tXm, tYm, tZm;
            private double[] gxmk_local, gymk_local, gzmk_local;
            private double[] txmk_local, tymk_local, tzmk_local;
            private double[] gXi, gYi, gZi, tXi, tYi, tZi;
            private double[] gxik_local, gyik_local, gzik_local;
            private double[] txik_local, tyik_local, tzik_local;
            // Store contributions to dE/dX/dL
            private double[] lgX, lgY, lgZ, ltX, ltY, ltZ;
            private double[] lxk_local, lyk_local, lzk_local;
            private double[] ltxk_local, ltyk_local, ltzk_local;
            private double[] esvPermRealDeriv_local;
            private double[] esvInducedRealDeriv_local;
            private double[] masking_local, maskingd_local, maskingp_local;
            private MultipoleTensor scrnTensor
                    = new MultipoleTensor(SCREENED_COULOMB, COORDINATES.QI, 5, aewald);
            private MultipoleTensor coulTensor
                    = new MultipoleTensor(COULOMB, COORDINATES.QI, 5, aewald);
            private MultipoleTensor tholeTensor
                    = new MultipoleTensor(THOLE_FIELD, COORDINATES.QI, 4, aewald);
            private MultipoleTensor scrnTensorPolar
                    = new MultipoleTensor(SCREENED_COULOMB, COORDINATES.QI, 5, aewald);
            private MultipoleTensor coulTensorPolar
                    = new MultipoleTensor(COULOMB, COORDINATES.QI, 5, aewald);
            private MultipoleTensor tholeTensorPolar
                    = new MultipoleTensor(THOLE_FIELD, COORDINATES.QI, 4, aewald);

            // Force and torque contributions for a single interaction.
            private final double[] permFi = new double[3], permTi = new double[3], permTk = new double[3];
            private final double[] polFi = new double[3], polTi = new double[3], polTk = new double[3];
            private final double[] FiC = new double[3], TiC = new double[3], TkC = new double[3];
            private final double[] FiT = new double[3], TiT = new double[3], TkT = new double[3];
            private final double[] dx_local = new double[3];
            private final double[][] rot_local = new double[3][3];
            private final double[][] work = new double[15][3];

            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            private void init() {
                if (masking_local == null || masking_local.length < nAtoms) {
                    txmk_local = new double[nAtoms];
                    tymk_local = new double[nAtoms];
                    tzmk_local = new double[nAtoms];
                    gxmk_local = new double[nAtoms];
                    gymk_local = new double[nAtoms];
                    gzmk_local = new double[nAtoms];
                    txik_local = new double[nAtoms];
                    tyik_local = new double[nAtoms];
                    tzik_local = new double[nAtoms];
                    gxik_local = new double[nAtoms];
                    gyik_local = new double[nAtoms];
                    gzik_local = new double[nAtoms];
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
                if (esvTerm) {
                    if (esvPermRealDeriv_local == null || esvPermRealDeriv_local.length < numESVs) {
                        esvPermRealDeriv_local = new double[numESVs];
                        esvInducedRealDeriv_local = new double[numESVs];
                    }
                    fill(esvPermRealDeriv_local, 0.0);
                    fill(esvInducedRealDeriv_local, 0.0);
                }
            }

            @Override
            public IntegerSchedule schedule() {
                return realSpaceSchedule;
            }

            @Override
            public void start() {
                // Set local references to master grad, torque, and field arrays.
                init();
                int threadIndex = getThreadIndex();
                realSpaceTimes[threadIndex] -= System.nanoTime();
                permanentEnergy = 0.0;
                inducedEnergy = 0.0;
                count = 0;
                gXm = gradPerm[threadIndex][0];
                gYm = gradPerm[threadIndex][1];
                gZm = gradPerm[threadIndex][2];
                gXi = gradInduced[threadIndex][0];
                gYi = gradInduced[threadIndex][1];
                gZi = gradInduced[threadIndex][2];
                tXm = torquePerm[threadIndex][0];
                tYm = torquePerm[threadIndex][1];
                tZm = torquePerm[threadIndex][2];
                tXi = torqueInduced[threadIndex][0];
                tYi = torqueInduced[threadIndex][1];
                tZi = torqueInduced[threadIndex][2];
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
                        fill(gxmk_local, 0.0);
                        fill(gymk_local, 0.0);
                        fill(gzmk_local, 0.0);
                        fill(txmk_local, 0.0);
                        fill(tymk_local, 0.0);
                        fill(tzmk_local, 0.0);
                        fill(gxik_local, 0.0);
                        fill(gyik_local, 0.0);
                        fill(gzik_local, 0.0);
                        fill(txik_local, 0.0);
                        fill(tyik_local, 0.0);
                        fill(tzik_local, 0.0);
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
                    realSpaceChunk(lb, ub, iSymm);
                    // Collect results.
                    if (gradient) {
                        // Turn symmetry mate torques into gradients
                        if (rotateMultipoles) {
                            torque(iSymm, txmk_local, tymk_local, tzmk_local,
                                    gxmk_local, gymk_local, gzmk_local,
                                    work[0], work[1], work[2], work[3], work[4],
                                    work[5], work[6], work[7], work[8], work[9],
                                    work[10], work[11], work[12], work[13], work[14]);
                            torque(iSymm, txik_local, tyik_local, tzik_local,
                                    gxik_local, gyik_local, gzik_local,
                                    work[0], work[1], work[2], work[3], work[4],
                                    work[5], work[6], work[7], work[8], work[9],
                                    work[10], work[11], work[12], work[13], work[14]);
                        }
                        // Rotate symmetry mate gradients
                        if (iSymm != 0) {
                            crystal.applyTransSymRot(nAtoms,
                                    gxmk_local, gymk_local, gzmk_local,
                                    gxmk_local, gymk_local, gzmk_local,
                                    symOp, rot_local);
                            crystal.applyTransSymRot(nAtoms,
                                    gxik_local, gyik_local, gzik_local,
                                    gxik_local, gyik_local, gzik_local,
                                    symOp, rot_local);
                        }
                        // Sum symmetry mate gradients into asymmetric unit gradients
                        for (int j = 0; j < nAtoms; j++) {
                            gXm[j] += gxmk_local[j];
                            gYm[j] += gymk_local[j];
                            gZm[j] += gzmk_local[j];
                            gXi[j] += gxik_local[j];
                            gYi[j] += gyik_local[j];
                            gZi[j] += gzik_local[j];
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
                    shareddEdLambda.addAndGet(dUdL);
                    sharedd2EdLambda2.addAndGet(d2UdL2);
                }
                if (esvTerm) {
                    /* Every-time, parallel reduction to shared ESV deriv. */
                    for (int esv = 0; esv < numESVs; esv++) {
                        esvPermRealDeriv_shared[esv].addAndGet(esvPermRealDeriv_local[esv]);
                        esvInducedRealDeriv_shared[esv].addAndGet(esvInducedRealDeriv_local[esv]);
                    }
                }
                realSpaceTimes[getThreadIndex()] += System.nanoTime();
            }

            /**
             * Evaluate the real space permanent energy and polarization energy
             * for a chunk of atoms.
             *
             * @param lb The lower bound of the chunk.
             * @param ub The upper bound of the chunk.
             */
            private void realSpaceChunk(final int lb, final int ub, final int iSymm) {
                final double x[] = coordinates[0][0];
                final double y[] = coordinates[0][1];
                final double z[] = coordinates[0][2];
                final int lists[][] = realSpaceLists[iSymm];
                final double neighborX[] = coordinates[iSymm][0];
                final double neighborY[] = coordinates[iSymm][1];
                final double neighborZ[] = coordinates[iSymm][2];
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    if (iSymm == 0) {
                        applyScaleFactors(i);
                        //applyMaskingRules(true, i, masking_local, maskingd_local, maskingp_local);
                    }
                    final double pdi = ipdamp[i];
                    final double pti = thole[i];
                    final int list[] = lists[i];
                    final int npair = realSpaceCounts[iSymm][i];
                    for (int j = 0; j < npair; j++) {
                        final int k = list[j];
                        if (!use[k]) {
                            continue;
                        }
                        final double xk = neighborX[k];
                        final double yk = neighborY[k];
                        final double zk = neighborZ[k];
                        boolean sameMolecule = (molecule[i] == molecule[k]);
                        if (lambdaMode == LambdaMode.VAPOR
                                && ((intermolecularSoftcore && !sameMolecule) || (intramolecularSoftcore && sameMolecule))) {
                            if (esvTerm) {
                                throw new UnsupportedOperationException();
                            }
                            continue;
                        }
                        soft = (isSoft[i] || isSoft[k]);
                        final LambdaFactors lfPerm, lfPol;
                        if (lambdaTerm && soft) {
                            if (true) throw new UnsupportedOperationException();
                            // Lambda Factors: Permanent
                            final double permLambdaProduct = (esvTerm)
                                    ? lambda * esvLambda[i] * esvLambda[k]
                                    : lambda;
                            final double sc1 = permLambdaAlpha * (1.0 - permLambdaProduct) * (1.0 - permLambdaProduct);
                            final double dsc1dL = -2.0 * permLambdaAlpha * (1.0 - permLambdaProduct);
                            final double d2sc1dL2 = 2.0 * permLambdaAlpha;
                            final double sc2 = pow(permLambdaProduct, permLambdaExponent);
                            final double dsc2dL = permLambdaExponent * pow(permLambdaProduct, permLambdaExponent - 1);
                            final double d2sc2dL2 = (permLambdaExponent < 2.0) ? 0.0
                                    : permLambdaExponent * (permLambdaExponent - 1) * pow(permLambdaProduct, permLambdaExponent - 2);
                            lfPerm = new LambdaFactors(sc1, dsc1dL, d2sc1dL2, sc2, dsc2dL, d2sc2dL2);

                            // Lambda Factors: Polarization
                            final double polLambdaProduct = (esvTerm)
                                    ? lambda * esvLambda[i] * esvLambda[k]
                                    : lambda;
                            final double sc1p = polLambdaAlpha * (1.0 - polLambdaProduct) * (1.0 - polLambdaProduct);
                            final double dsc1dLp = -2.0 * polLambdaAlpha * (1.0 - polLambdaProduct);
                            final double d2sc1dL2p = -2.0 * polLambdaAlpha;
                            final double sc2p = pow(polLambdaProduct, polLambdaExponent);
                            final double dsc2dLp = polLambdaExponent * pow(polLambdaProduct, polLambdaExponent - 1);
                            final double d2sc2dL2p = (polLambdaExponent < 2.0) ? 0.0
                                    : polLambdaExponent * (polLambdaExponent - 1) * pow(polLambdaProduct, polLambdaExponent - 2);
                            lfPol = new LambdaFactors(sc1p, dsc1dLp, d2sc1dL2p, sc2p, dsc2dLp, d2sc2dL2p);
                        } else {
                            // Not soft; use defaults (0.0, 0.0, 0.0, 1.0, 0.0, 0.0).
                            lfPerm = LambdaFactors.Defaults;
                            lfPol = LambdaFactors.Defaults;
                        }
                        dx_local[0] = xk - xi;
                        dx_local[1] = yk - yi;
                        dx_local[2] = zk - zi;
                        crystal.image(dx_local);

                        final double pdk = ipdamp[k];   // == 1/polarizability^6
                        final double ptk = thole[k];
                        double pgamma = min(pti, ptk);
                        double aiak = pdi * pdk;

                        final double scalePerm = masking_local[k];
                        final double scalePolGroup = maskingd_local[k];
                        final double scalePolEnergy = maskingp_local[k];
                        final boolean applyScalePerm = scalePerm != 1.0;
                        final boolean applyScalePol = (scalePolGroup != 1.0 || scalePolEnergy != 1.0);
                        final boolean damped = MultipoleTensor.checkDampingCriterion(dx_local, pgamma, aiak);

                        permanentEnergy += interactPermanent(i, k, iSymm, dx_local, lfPerm, applyScalePerm, scalePerm);
                        inducedEnergy += interactPolarization(i, k, iSymm, dx_local, lfPol, applyScalePol, damped, scalePolGroup, scalePolEnergy , pgamma, aiak);
                        count++;
                    }
                    if (iSymm == 0) {
                        resetScaleFactors(i);
                        //applyMaskingRules(false, i, masking_local, maskingd_local, maskingp_local);
                    }
                }	// loop over i
            }

            private double interactPolarization(
                    final int i, final int k, final int iSymm, final double[] dx_local,
                    LambdaFactors lf, final boolean applyScalePol, final boolean damped,
                    final double scaleD, final double scaleP,
                    final double pgamma, final double aiak) {
                final double[] Qi = globalMultipole[0][i];
                final double[] Qk = globalMultipole[iSymm][k];
                final double[] ui = inducedDipole[0][i];
                final double[] vi = inducedDipoleCR[0][i];
                final double[] uk = inducedDipole[iSymm][k];
                final double[] vk = inducedDipoleCR[iSymm][k];
                final double selfScale = (i == k) ? 0.5 : 1.0;
                final double mutualScale = (polarization == Polarization.MUTUAL) ? 1.0 : 0.0;
                final double prefactor = (0.5 * ELECTRIC * polarizationScale * selfScale);

                scrnTensorPolar.generateTensor(dx_local, Qi, Qk, ui, vi, uk, vk);
                double energy = scrnTensorPolar.polarizationEnergy(1.0, 1.0, mutualScale, polFi, polTi, polTk);
                /**
                 * Subtract away masked Coulomb interactions included in PME.
                 */
                if (applyScalePol) {
                    coulTensorPolar.generateTensor(dx_local, Qi, Qk, ui, vi, uk, vk);
                    energy -= coulTensorPolar.polarizationEnergy(1.0 - scaleD, 1.0 - scaleP, 0.0, FiC, TiC, TkC);
                    for (int axis = 0; axis < 3; axis++) {
                        polFi[axis] -= FiC[axis];
                        polTi[axis] -= TiC[axis];
                        polTk[axis] -= TkC[axis];
                    }
                }

                if (damped) {
                    tholeTensorPolar.generateTensor(dx_local, Qi, Qk, ui, vi, uk, vk, pgamma, aiak);
                    energy -= tholeTensorPolar.polarizationEnergy(scaleD, scaleP, mutualScale, FiT, TiT, TkT);
                    for (int axis = 0; axis < 3; axis++) {
                        polFi[axis] -= FiT[axis];
                        polTi[axis] -= TiT[axis];
                        polTk[axis] -= TkT[axis];
                    }
                }
                energy *= prefactor;

                final double gPrefactor = ELECTRIC * polarizationScale * selfScale;
                transposeForcesToGradByAxis(
                        polFi, polTi, polTk, i, k, gPrefactor,
                        gXi, gYi, gZi, tXi, tYi, tZi,
                        gxik_local, gyik_local, gzik_local,
                        txik_local, tyik_local, tzik_local);
                if (lambdaTerm) {
                    final double lPrefactor = ELECTRIC * dEdLSign * lf.dsc2dL * selfScale;
                    transposeForcesToGradByAxis(
                            polFi, polTi, polTk, i, k, lPrefactor,
                            lgX, lgY, lgZ, ltX, ltY, ltZ,
                            lxk_local, lyk_local, lzk_local,
                            ltxk_local, ltyk_local, ltzk_local);

                    final double dLp_dL = (esvTerm) ? esvLambda[i] * esvLambda[k] : 1.0;
                    dUdL += dEdLSign * dlPowPol * energy * dLp_dL;
                    d2UdL2 += dEdLSign * d2lPowPol * energy * dLp_dL;
                }

                // Compute ExtendedSystem lambda derivative(s).
                if (esvTerm) {
                    final double crScale = scaleD, regScale = scaleP;
                    // Set common derivative components: dotted multipoles, dipole scaling, prefactor.
                    final double[] Qidot = esvMultipoleDot[0][i];
                    final double[] Qkdot = esvMultipoleDot[iSymm][k];
                    // Collect influence on dipoles due to permanent scaling of atom i.
                    if (esvAtomsScaled[i]) {
                        scrnTensorPolar.generateTensor(dx_local, Qidot, zeroM, zeroD, zeroD, uk, vk);
                        double dUdLi = scrnTensor.indkFieldsDotI(1.0, 1.0);
                        if (applyScalePol) {
                            coulTensorPolar.generateTensor(dx_local, Qidot, zeroM, zeroD, zeroD, uk, vk);
                            dUdLi -= coulTensor.indkFieldsDotI(1.0 - crScale, 1.0 - regScale);
                        }
                        if (damped) {
                            tholeTensorPolar.setTholeDamping(pgamma, aiak);
                            tholeTensorPolar.generateTensor(dx_local, Qidot, zeroM, zeroD, zeroD, uk, vk, pgamma, aiak);
                            dUdLi -= tholeTensor.indkFieldsDotI(crScale, regScale);
                        }
                        esvInducedRealDeriv_local[esvIndex[i]] += dUdLi * prefactor;
                    }
                    // Collect influence on dipoles due to permanent scaling of atom k.
                    if (esvAtomsScaled[k]) {
                        scrnTensorPolar.generateTensor(dx_local, zeroM, Qkdot, ui, vi, zeroD, zeroD);
                        double dUdLk = scrnTensorPolar.indiFieldsDotK(1.0, 1.0);
                        if (applyScalePol) {
                            coulTensorPolar.generateTensor(dx_local, zeroM, Qkdot, ui, vi, zeroD, zeroD);
                            dUdLk -= coulTensorPolar.indiFieldsDotK(1.0 - crScale, 1.0 - regScale);
                        }
                        if (damped) {
                            tholeTensorPolar.setTholeDamping(pgamma, aiak);
                            tholeTensorPolar.generateTensor(dx_local, zeroM, Qkdot, ui, vi, zeroD, zeroD, pgamma, aiak);
                            dUdLk -= tholeTensorPolar.indiFieldsDotK(crScale, regScale);
                        }
                        esvInducedRealDeriv_local[esvIndex[k]] += dUdLk * prefactor;
                    }

                    /* Collect influence on mutual polarization due to alpha (and, implicitly, mu) scaling. */
                    final boolean titrHi = esvAtomsScaledAlpha[i];
                    final boolean titrHk = esvAtomsScaledAlpha[k];
                    if (esvConfig.scaleAlpha && (titrHi || titrHk)) {
                        final double polarI = polarizability[i];
                        final double polarK = polarizability[k];
                        final double polarUnsI = unscaledPolarizability[i];
                        final double polarUnsK = unscaledPolarizability[k];
                        final double[] ui_uns = unscaledInducedDipole[0][i];
                        final double[] vi_uns = unscaledInducedDipoleCR[0][i];
                        final double[] uk_uns = unscaledInducedDipole[iSymm][k];
                        final double[] vk_uns = unscaledInducedDipoleCR[iSymm][k];

                        double ui_unsDotvk = 0.0, uk_unsDotvi = 0.0;
                        if (!applyScalePol) {
                            coulTensorPolar.generateTensor(dx_local, Qi, Qk,
                                    ui_uns, vi, uk_uns, vk);
                            ui_unsDotvk = coulTensorPolar.uiDotvk();
                            uk_unsDotvi = coulTensorPolar.ukDotvi();
                        } else {
                            scrnTensorPolar.generateTensor(dx_local, Qi, Qk,
                                    ui_uns, vi, uk_uns, vk);
                            ui_unsDotvk = scrnTensorPolar.uiDotvk();
                            uk_unsDotvi = scrnTensorPolar.ukDotvi();
                            coulTensorPolar.generateTensor(dx_local, Qi, Qk,
                                    ui_uns, vi, uk_uns, vk);
                            ui_unsDotvk -= (1.0 - scaleP) * coulTensorPolar.uiDotvk();
                            uk_unsDotvi -= (1.0 - scaleP) * coulTensorPolar.ukDotvi();
                        }
                        if (damped) {
                            tholeTensorPolar.setTholeDamping(pgamma, aiak);
                            tholeTensorPolar.generateTensor(dx_local, Qi, Qk,
                                    ui_uns, vi, uk_uns, vk, pgamma, aiak);
                            ui_unsDotvk -= scaleP * tholeTensorPolar.uiDotvk();
                            uk_unsDotvi -= scaleP * tholeTensorPolar.ukDotvi();
                        }
                        
                        if (titrHi) {
                            final double alphaLambdaChainI = -1.0 / (polarUnsI * esvLambda[i] * esvLambda[i]);
                            esvInducedRealDeriv_local[esvIndex[i]]
                                    += ui_unsDotvk * prefactor * alphaLambdaChainI;
                        }
                        if (titrHk) {
                            final double alphaLambdaChainK = -1.0 / (polarUnsK * esvLambda[k] * esvLambda[k]);
                            esvInducedRealDeriv_local[esvIndex[k]]
                                    += uk_unsDotvi * prefactor * alphaLambdaChainK;
                        }
                    }
                }
                return energy;
            }

            private double interactPermanent(
                    final int i, final int k, final int iSymm, final double[] dx_local,
                    LambdaFactors lf, final boolean applyScalePerm, final double scalePerm) {
                final double[] Qi = globalMultipole[0][i];
                final double[] Qk = globalMultipole[iSymm][k];
                final double selfScale = (i == k) ? 0.5 : 1.0;

                /* Compute screened real space interactions.    */
                final double scale1 = 1.0 - scalePerm;
                final double prefactor = ELECTRIC * permanentScale * lf.sc2 * selfScale;
                double dz = 0.0;

                scrnTensor.generateTensor(dx_local, lf.sc1, Qi, Qk);
                double energy = scrnTensor.multipoleEnergy(permFi, permTi, permTk);
                if (soft) {
                    dz = prefactor * scrnTensor.getdEdZbuff() * lf.dsc1dL;
                }
                if (applyScalePerm) {
                    coulTensor.generateTensor(dx_local, lf.sc1, Qi, Qk);
                    energy -= scale1 * coulTensor.multipoleEnergy(FiC, TiC, TkC);
                    for (int axis = 0; axis < 3; axis++) {
                        permFi[axis] -= scale1*FiC[axis];
                        permTi[axis] -= scale1*TiC[axis];
                        permTk[axis] -= scale1*TkC[axis];
                    }
                    if (soft) {
                        dz -= scale1 * coulTensor.getdEdZbuff() * lf.dsc1dL;
                    }
                }
                energy *= prefactor;

                if (gradient) {
                    final double gPrefactor = ELECTRIC * permanentScale * selfScale;
                    transposeForcesToGradByAxis(
                            permFi, permTi, permTk, i, k, gPrefactor,
                            gXm, gYm, gZm, tXm, tYm, tZm,
                            gxmk_local, gymk_local, gzmk_local,
                            txmk_local, tymk_local, tzmk_local);
                }

                /**
                 * If both atoms are scaled by the same ESV, their contributions
                 * are additive. Scaled Mdot multipoles are pre-weighted by the
                 * switch derivative.
                 */
                if (esvTerm) {
                    if (esvAtomsScaled[i]) {
                        final double[] Qidot = esvMultipoleDot[0][i];
                        scrnTensor.generateTensor(dx_local, lf.sc1, Qidot, Qk);
                        double dUdLi = scrnTensor.mpoleKFieldDotIqi();
                        if (applyScalePerm) {
                            coulTensor.generateTensor(dx_local, lf.sc1, Qidot, Qk);
                            dUdLi -= scale1 * coulTensor.mpoleKFieldDotIqi();
                        }
                        esvPermRealDeriv_local[esvIndex[i]] += dUdLi * prefactor;
                    }
                    if (esvAtomsScaled[k]) {
                        final double[] Qkdot = esvMultipoleDot[iSymm][k];
                        scrnTensor.generateTensor(dx_local, lf.sc1, Qi, Qkdot);
                        double dUdLk = scrnTensor.mpoleIFieldDotKqi();
                        if (applyScalePerm) {
                            coulTensor.generateTensor(dx_local, lf.sc1, Qi, Qkdot);
                            dUdLk -= scale1 * coulTensor.mpoleIFieldDotKqi();
                        }
                        esvPermRealDeriv_local[esvIndex[k]] += dUdLk * prefactor;
                    }
                }

                if (lambdaTerm && soft) {
                    final double F = lf.sc1;            // Buffer, alpha*(1-Lp)^2
                    final double dFdL = lf.dsc1dL;      // dBuffer_dLp
                    final double S = lf.sc2;            // Lp
                    final double dSdL = lf.dsc2dL;      // dLpdL
                    double P = energy;					// Energy
                    double dPdF = dz;					// dEnergy/dBuffer

                    final double termA = dSdL * P;
                    final double termB = S * dPdF * dFdL;
                    final double dU_dLambdaProduct = selfScale * (termA + termB);
                    final double dLambdaProduct_dL = (esvTerm) ? esvLambda[i] * esvLambda[k] : 1.0;
                    dUdL += dU_dLambdaProduct * dLambdaProduct_dL;

                    /* dU/dL/dX, of first term: d[dlPow * ereal]/dx             */
                    final double lbdScale1 = scalePerm * lf.dsc2dL;
                    transposeForcesToGradByAxis(
                            permFi, permTi, permTk, i, k, lbdScale1,
                            lgX, lgY, lgZ, ltX, ltY, ltZ,
                            lxk_local, lyk_local, lzk_local,
                            ltxk_local, ltyk_local, ltzk_local);

                    /* dU/dL/dX, of second term: d[lPow*dlAlpha*dRealdL]/dX     */
                    // No additional call to MT; use 6th order tensor instead.
                    final double lbdScale2 = scalePerm * lf.sc2 * lf.dsc1dL;
                    transposeForcesToGradByAxis(
                            permFi, permTi, permTk, i, k, lbdScale2,
                            lgX, lgY, lgZ, ltX, ltY, ltZ,
                            lxk_local, lyk_local, lzk_local,
                            ltxk_local, ltyk_local, ltzk_local);
                }
                return energy;
            }

            private void applyScaleFactors(int i) {
                final Atom ai = atoms[i];
                for (Atom ak : ai.get1_5s()) {
                    masking_local[ak.getIndex() - 1] = m15scale;
                }
                for (Torsion torsion : ai.getTorsions()) {
                    Atom ak = torsion.get1_4(ai);
                    if (ak != null) {
                        int index = ak.getIndex() - 1;
                        masking_local[index] = m14scale;
                        for (int j : ip11[i]) {
                            if (j == index) {
                                maskingp_local[index] = intra14Scale;
                            }
                        }
                    }
                }
                for (Angle angle : ai.getAngles()) {
                    Atom ak = angle.get1_3(ai);
                    if (ak != null) {
                        int index = ak.getIndex() - 1;
                        masking_local[index] = m13scale;
                        maskingp_local[index] = p13scale;
                    }
                }
                for (Bond bond : ai.getBonds()) {
                    int index = bond.get1_2(ai).getIndex() - 1;
                    masking_local[index] = m12scale;
                    maskingp_local[index] = p12scale;
                }
                for (int j : ip11[i]) {
                    maskingd_local[j] = d11scale;
                }
            }

            private void resetScaleFactors(int i) {
                final Atom ai = atoms[i];
                for (Atom ak : ai.get1_5s()) {
                    int index = ak.getIndex() - 1;
                    masking_local[index] = 1.0;
                }
                for (Torsion torsion : ai.getTorsions()) {
                    Atom ak = torsion.get1_4(ai);
                    if (ak != null) {
                        int index = ak.getIndex() - 1;
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
                        int index = ak.getIndex() - 1;
                        masking_local[index] = 1.0;
                        maskingp_local[index] = 1.0;
                    }
                }
                for (Bond bond : ai.getBonds()) {
                    int index = bond.get1_2(ai).getIndex() - 1;
                    masking_local[index] = 1.0;
                    maskingp_local[index] = 1.0;
                }
                for (int j : ip11[i]) {
                    maskingd_local[j] = 1.0;
                }
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
        private PermanentReciprocalEnergyLoop permanentReciprocalEnergyLoop[];
        private InducedDipoleReciprocalEnergyLoop inducedDipoleReciprocalEnergyLoop[];

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
            // Every-time, single-threaded reset of esv deriv.
            if (esvTerm) {
                for (int i = 0; i < numESVs; i++) {
                    esvPermSelfDeriv_shared[i].getAndSet(0.0);
                    esvPermRecipDeriv_shared[i].getAndSet(0.0);
                    esvInducedSelfDeriv_shared[i].getAndSet(0.0);
                    esvInducedRecipDeriv_shared[i].getAndSet(0.0);
                }
            }
        }

        @Override
        public void run() {
            if (skipRecipRegion || (esvTerm && !esvConfig.allowRecipRegion)) {
                logger.warning("Skipping PME reciprocal region.");
                return;
            }
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
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the reciprocal space energy in thread " + getThreadIndex());
                throw ex;
            } catch (Exception ex) {
                String msg = "Fatal exception computing the reciprocal space energy in thread " + getThreadIndex();
                logger.log(Level.SEVERE, msg, ex);
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
            for (int thread = 0; thread < maxThreads; thread++) {
                if (permanentReciprocalEnergyLoop[thread] != null) {
                    permanentSelfEnergy += permanentReciprocalEnergyLoop[thread].eSelf;
                    permanentReciprocalEnergy += permanentReciprocalEnergyLoop[thread].eRecip;
                }
            }
        }

        private class PermanentReciprocalEnergyLoop extends IntegerForLoop {

            private double gX[], gY[], gZ[], tX[], tY[], tZ[];
            private double lgX[], lgY[], lgZ[], ltX[], ltY[], ltZ[];
            protected double eSelf;
            protected double eRecip;
            private double[] esvPermSelfDeriv_local;
            private double[] esvPermRecipDeriv_local;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eSelf = 0.0;
                eRecip = 0.0;
                int ti = getThreadIndex();
                gX = gradPerm[ti][0];
                gY = gradPerm[ti][1];
                gZ = gradPerm[ti][2];
                tX = torquePerm[ti][0];
                tY = torquePerm[ti][1];
                tZ = torquePerm[ti][2];
                if (lambdaTerm) {
                    lgX = lambdaGrad[ti][0];
                    lgY = lambdaGrad[ti][1];
                    lgZ = lambdaGrad[ti][2];
                    ltX = lambdaTorque[ti][0];
                    ltY = lambdaTorque[ti][1];
                    ltZ = lambdaTorque[ti][2];
                }
                if (esvTerm) {
                    if (esvPermSelfDeriv_local == null || esvPermSelfDeriv_local.length < numESVs) {
                        esvPermSelfDeriv_local = new double[numESVs];
                        esvPermRecipDeriv_local = new double[numESVs];
                    }
                    fill(esvPermSelfDeriv_local, 0.0);
                    fill(esvPermRecipDeriv_local, 0.0);
                }
            }

            @Override
            public void run(int lb, int ub) {
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
                        if (esvAtomsScaled[i]) {
                            final double indot[] = esvMultipoleDot[0][i];
                            final double ciidot = indot[t000] * in[t000];
                            final double diidot = indot[t100] * in[t100] + indot[t010] * in[t010] + indot[t001] * in[t001];
                            final double qiidot = indot[t200] * in[t200] + indot[t020] * in[t020] + indot[t002] * in[t002]
                                    + 2.0 * (indot[t110] * in[t110] + indot[t101] * in[t101] + indot[t011] * in[t011]);
                            final double eSelfDot = 2.0 * aewald1 * (ciidot + aewald2 * (diidot / 3.0 + 2.0 * aewald2 * qiidot / 45.0));
                            esvPermSelfDeriv_local[esvIndex[i]] += eSelfDot;
                        }
                    }
                }
                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(eSelf * dsc2dL * dEdLSign);
                    sharedd2EdLambda2.addAndGet(eSelf * d2sc2dL2 * dEdLSign);
                }
                /**
                 * Permanent multipole reciprocal space energy and gradient.
                 */
                final double recip[][] = crystal.getUnitCell().A;

                double dUdL = 0.0;
                double d2UdL2 = 0.0;
                for (int i = lb; i <= ub; i++) {
                    if (use[i]) {
                        final double phi[] = cartMultipolePhi[i];
                        final double mpole[] = multipole[i];

                        double e = mpole[t000] * phi[t000] + mpole[t100] * phi[t100]
                                + mpole[t010] * phi[t010] + mpole[t001] * phi[t001]
                                + oneThird
                                * (mpole[t200] * phi[t200]
                                + mpole[t020] * phi[t020]
                                + mpole[t002] * phi[t002]
                                + 2.0 * (mpole[t110] * phi[t110]
                                + mpole[t101] * phi[t101]
                                + mpole[t011] * phi[t011]));
                        eRecip += e;

                        if (esvAtomsScaled[i]) {
                            final double mpoleDot[] = esvMultipoleDot[0][i];
                            final double edot = mpoleDot[t000] * phi[t000]
                                    + mpoleDot[t100] * phi[t100]
                                    + mpoleDot[t010] * phi[t010]
                                    + mpoleDot[t001] * phi[t001]
                                    + oneThird * (mpoleDot[t200] * phi[t200]
                                    + mpoleDot[t020] * phi[t020]
                                    + mpoleDot[t002] * phi[t002]
                                    + 2.0 * (mpoleDot[t110] * phi[t110]
                                    + mpoleDot[t101] * phi[t101]
                                    + mpoleDot[t011] * phi[t011]));
                            esvPermRecipDeriv_local[esvIndex[i]] += 2.0 * edot;
                        }

                        if (gradient || lambdaTerm) {
                            final double fmpole[] = fracMultipoles[i];
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
                                dUdL += dEdLSign * dsc2dL * e;
                                d2UdL2 += dEdLSign * d2sc2dL2 * e;
                                lgX[i] += dEdLSign * dsc2dL * ELECTRIC * dfx;
                                lgY[i] += dEdLSign * dsc2dL * ELECTRIC * dfy;
                                lgZ[i] += dEdLSign * dsc2dL * ELECTRIC * dfz;
                                ltX[i] += dEdLSign * dsc2dL * ELECTRIC * tqx;
                                ltY[i] += dEdLSign * dsc2dL * ELECTRIC * tqy;
                                ltZ[i] += dEdLSign * dsc2dL * ELECTRIC * tqz;
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
                if (esvTerm) {
                    // Every-time, parallel reduction to shared esv deriv.
                    for (int i = 0; i < numESVs; i++) {
                        esvPermSelfDeriv_shared[i].addAndGet(permanentScale * esvPermSelfDeriv_local[i]);
                        esvPermRecipDeriv_shared[i].addAndGet(permanentScale * 0.5 * ELECTRIC * esvPermRecipDeriv_local[i]);
                    }
                }
            }
        }

        private class InducedDipoleReciprocalEnergyLoop extends IntegerForLoop {

            private double eSelf;
            private double eRecip;
            private double gX[], gY[], gZ[], tX[], tY[], tZ[];
            private double lgX[], lgY[], lgZ[], ltX[], ltY[], ltZ[];
            private final double sfPhi[] = new double[tensorCount];
            private final double sPhi[] = new double[tensorCount];
            private double[] esvInducedSelfDeriv_local;
            private double[] esvInducedRecipDeriv_local;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                eSelf = 0.0;
                eRecip = 0.0;
                int threadID = getThreadIndex();
                gX = gradInduced[threadID][0];
                gY = gradInduced[threadID][1];
                gZ = gradInduced[threadID][2];
                tX = torqueInduced[threadID][0];
                tY = torqueInduced[threadID][1];
                tZ = torqueInduced[threadID][2];
                if (lambdaTerm) {
                    lgX = lambdaGrad[threadID][0];
                    lgY = lambdaGrad[threadID][1];
                    lgZ = lambdaGrad[threadID][2];
                    ltX = lambdaTorque[threadID][0];
                    ltY = lambdaTorque[threadID][1];
                    ltZ = lambdaTorque[threadID][2];
                }
                if (esvTerm) {
                    if (esvInducedSelfDeriv_local == null || esvInducedSelfDeriv_local.length < numESVs) {
                        esvInducedSelfDeriv_local = new double[numESVs];
                        esvInducedRecipDeriv_local = new double[numESVs];
                    }
                    fill(esvInducedSelfDeriv_local, 0.0);
                    fill(esvInducedRecipDeriv_local, 0.0);
                }
            }

            @Override
            public void run(int lb, int ub) {
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
                        if (esvAtomsScaled[i]) {
                            double[] mDot = esvMultipoleDot[0][i];
                            double[] unsi = unscaledInducedDipole[0][i];
                            final double indiDotMdot
                                    = (indi[0] + unsi[0]) * mDot[t100]
                                    + (indi[1] + unsi[1]) * mDot[t010]
                                    + (indi[2] + unsi[2]) * mDot[t001];
                            final double dIndSelfdLi = aewald3 * indiDotMdot;
                            esvInducedSelfDeriv_local[esvIndex[i]] += dIndSelfdLi;
                        }
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
                        if (esvAtomsScaled[i]) {
                            final double[] fracMDotPhi = reciprocalSpace.getFracMultipoleDotPhi()[i];
                            final double[] indUns = unscaledInducedDipole[0][i];
                            final double indUnsx = indUns[0];
                            final double indUnsy = indUns[1];
                            final double indUnsz = indUns[2];
                            final double dRecipdL = 2.0 * (indx * fracMDotPhi[t100] + indUnsx * fPhi[t100]
                                    + indy * fracMDotPhi[t010] + indUnsy * fPhi[t010]
                                    + indz * fracMDotPhi[t001] + indUnsz * fPhi[t001]);
                            /*
								In reciprocal space, induced dipoles interact with their own permanent multipoles,
							via a self-energy term, albeit only those of symmetry mate copies.
								In addition, induced dipoles contribute to the electric field used in the self-consistent
							iterator. The lambda derivative of a titrating hydrogen's induced dipole site with
							lambda-scaled polarizability thus consists of four terms. They are the result of a product
							rule for both self and reciprocal energy components, each of which takes the form
							of a dot product between dipole and field).
								The induced dipole self energy is obtained by a dot product
                             */
                            esvInducedRecipDeriv_local[esvIndex[i]] += dRecipdL;
                        }
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
                if (lambdaTerm) {
                    shareddEdLambda.addAndGet(dEdLSign * dlPowPol * eRecip);
                    sharedd2EdLambda2.addAndGet(dEdLSign * d2lPowPol * eRecip);
                }
            }

            @Override
            public void finish() {
                inducedDipoleSelfEnergy.addAndGet(eSelf * polarizationScale);
                inducedDipoleRecipEnergy.addAndGet(eRecip * polarizationScale * 0.5 * ELECTRIC);
                if (esvTerm) {
                    for (int i = 0; i < numESVs; i++) {
                        esvInducedSelfDeriv_shared[i].addAndGet(esvInducedSelfDeriv_local[i] * polarizationScale);
                        esvInducedRecipDeriv_shared[i].addAndGet(esvInducedRecipDeriv_local[i] * polarizationScale * 0.5 * ELECTRIC);
                    }
                }
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
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception initializing coordinates in thread: " + threadIndex);
                throw ex;
            } catch (Exception ex) {
                String msg = "Fatal exception initializing coordinates in thread: " + threadIndex;
                logger.log(Level.SEVERE, msg, ex);
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
                    for (int axis = 0; axis < 3; axis++) {
                        fill(grad[threadID][axis], 0.0);
                        fill(gradPerm[threadID][axis], 0.0);
                        fill(torquePerm[threadID][axis], 0.0);
                        fill(gradInduced[threadID][axis], 0.0);
                        fill(torqueInduced[threadID][axis], 0.0);
                    }
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

                    for (int j = 0; j < maxThreads; j++) {
                        grad[j][0][i] = 0.0;
                        grad[j][1][i] = 0.0;
                        grad[j][2][i] = 0.0;
                        torque[j][0][i] = 0.0;
                        torque[j][1][i] = 0.0;
                        torque[j][2][i] = 0.0;
                    }

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
                if (esvTerm) {
                    List<Double> scales = esvConfig.cdqScales;
                    chargeScale = scales.get(0);
                    dipoleScale = scales.get(1);
                    traceScale = scales.get(2);
                }
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
                        /* For shared ESV atoms, pipe in the interpolated multipole instead. */
                        final double in[] = (esvTerm && esvAtomsScaled[ii])
                                ? atom.getEsvMultipole().getMultipole()
                                : localMultipole[ii];
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
                            localQuadrupole[1][0] = in[t110];
                            localQuadrupole[1][2] = in[t011];
                            localQuadrupole[2][0] = in[t101];
                            localQuadrupole[2][1] = in[t011];
                            // Check for chiral flipping.
                            boolean needsChiralInversion = false;
                            if (frame[ii] == MultipoleType.MultipoleFrameDefinition.ZTHENX
                                    && referenceSites.length == 3) {
                                needsChiralInversion = checkMultipoleChirality(in, frame[ii], localOrigin, frameCoords);
                                if (needsChiralInversion) {
                                    localDipole[1] = -localDipole[1];
                                    localQuadrupole[0][1] = -localQuadrupole[0][1];
                                    localQuadrupole[1][0] = -localQuadrupole[1][0];
                                    localQuadrupole[1][2] = -localQuadrupole[1][2];
                                    localQuadrupole[2][1] = -localQuadrupole[2][1];
                                }
                            }
                            // Do the rotation.
                            double[][] rotmat = MultipoleType.getRotationMatrix(frame[ii], localOrigin, frameCoords);
                            MultipoleType.rotateMultipole(rotmat, localDipole, localQuadrupole,
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
                            /* For ESV atoms, also rotate and scale the Mdot multipole. */
                            if (esvTerm && esvAtomsScaled[ii]) {
                                final double mdotCharge = atom.getEsvMultipoleDot().getCharge();
                                final double[] mdotDipole = atom.getEsvMultipoleDot().getDipole();
                                final double[][] mdotQuad = atom.getEsvMultipoleDot().getQuadrupole();
                                if (needsChiralInversion) {
                                    mdotDipole[1] = -mdotDipole[1];
                                    mdotQuad[0][1] = -mdotQuad[0][1];
                                    mdotQuad[1][0] = -mdotQuad[1][0];
                                    mdotQuad[1][2] = -mdotQuad[1][2];
                                    mdotQuad[2][1] = -mdotQuad[2][1];
                                }
                                double[] rotDipole = new double[3];
                                double[][] rotQuad = new double[3][3];
                                MultipoleType.rotateMultipole(rotmat, mdotDipole, mdotQuad, rotDipole, rotQuad);
                                esvMultipoleDot[iSymm][ii][t000] = mdotCharge * chargeScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t100] = rotDipole[0] * dipoleScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t010] = rotDipole[1] * dipoleScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t001] = rotDipole[2] * dipoleScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t200] = rotQuad[0][0] * traceScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t020] = rotQuad[1][1] * traceScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t002] = rotQuad[2][2] * traceScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t110] = rotQuad[0][1] * traceScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t101] = rotQuad[0][2] * traceScale * elecScale;
                                esvMultipoleDot[iSymm][ii][t011] = rotQuad[1][2] * traceScale * elecScale;
                            }
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
                        double polar = polarizeType.polarizability * elecScale;
                        if (esvTerm && esvAtomsScaled[ii]) {
                            final double lii = esvLambda[ii];
                            polarizability[ii] = lii * polar;
                            unscaledPolarizability[ii] = polar;
                        } else {
                            polarizability[ii] = polar;
                        }
                    }
                }
            }
        }
    }

    private class ExpandInducedDipolesRegion extends ParallelRegion {

        private final ExpandInducedDipoleLoop expandInducedDipoleLoop[];

        public ExpandInducedDipolesRegion(int maxThreads) {
            expandInducedDipoleLoop = new ExpandInducedDipoleLoop[maxThreads];
            for (int thread = 0; thread < maxThreads; thread++) {
                expandInducedDipoleLoop[thread] = new ExpandInducedDipoleLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, expandInducedDipoleLoop[getThreadIndex()]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception expanding coordinates in thread: " + getThreadIndex());
                throw ex;
            } catch (Exception e) {
                String msg = "Fatal exception expanding coordinates in thread: " + getThreadIndex();
                logger.log(Level.SEVERE, msg, e);
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
                        if (esvTerm) {
                            crystal.applySymRot(unscaledInducedDipole[0][ii], unscaledInducedDipole[s][ii], symOp);
                            crystal.applySymRot(unscaledInducedDipoleCR[0][ii], unscaledInducedDipoleCR[s][ii], symOp);
                        }
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
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing torque in thread " + getThreadIndex());
                throw ex;
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
            private double gm[][];
            private double gi[][];
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
                gm = gradPerm[threadID];
                gi = gradInduced[threadID];
                if (lambdaTerm) {
                    lg = lambdaGrad[threadID];
                }
            }

            @Override
            public void run(int lb, int ub) {
                if (gradient) {
                    for (int i = lb; i <= ub; i++) {
                        torque(i, torquePerm, gm);
                        torque(i, torqueInduced, gi);
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
            public void run(int lb, int ub) {
                if (gradient) {
                    // Sum permanent and induced into combined array.
                    for (int thread = 0; thread < maxThreads; thread++) {
                        for (int atom = lb; atom <= ub; atom++) {
                            for (int axis = 0; axis < 3; axis++) {
                                grad[thread][axis][atom] =
                                        gradPerm[thread][axis][atom] + gradInduced[thread][axis][atom];
                            }
                        }
                    }
                    // Reduce combined array.
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
                    // Add gradient to Atom objects.
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

    @Override
    public double getEwaldCutoff() {
        return off;
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
            Integer index = ai.getIndex() - 1;
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
                logger.log(Level.SEVERE, message, new IllegalStateException());
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
                    Integer index = aj.getIndex() - 1;
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
            initSoftCore(true);
        }

        sc1 = 0.0;
        dsc1dL = 0.0;
        d2sc1dL2 = 0.0;
        sc2 = 1.0;
        dsc2dL = 0.0;
        d2sc2dL2 = 0.0;
        if (lambda < permLambdaStart) {
            sc2 = 0.0;
        } else if (lambda <= permLambdaEnd) {
            double permLambdaScale = 1.0 / (permLambdaStart - permLambdaEnd);
            double permLambda = (lambda - permLambdaStart) * permLambdaScale;
            sc1 = permLambdaAlpha * (1.0 - permLambda) * (1.0 - permLambda);
            dsc1dL = -2.0 * permLambdaAlpha * (1.0 - permLambda);
            d2sc1dL2 = 2.0 * permLambdaAlpha;

            sc2 = pow(permLambda, permLambdaExponent);
            dsc2dL = permLambdaExponent * pow(permLambda, permLambdaExponent - 1.0);
            d2sc2dL2 = 0.0;
            if (permLambdaExponent >= 2.0) {
                d2sc2dL2 = permLambdaExponent * (permLambdaExponent - 1.0) * pow(permLambda, permLambdaExponent - 2.0);
            }
            // Treat chain rule for the activation window.
            dsc1dL *= permLambdaScale;
            d2sc1dL2 *= (permLambdaScale * permLambdaScale);
            dsc2dL *= permLambdaScale;
            d2sc2dL2 *= (permLambdaScale * permLambdaScale);
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
        if (lambdaTerm) {
            logger.warning("Extended system attachment resulted in lambdaTerm being disabled.");
            lambdaTerm = false;
        }
        // Set object handles.
        esvTerm = true;
        esvSystem = system;
        numESVs = esvSystem.size();
        esvConfig = esvSystem.getConfig();
        // Update atoms and reinitialize arrays.
        /* Only include ExtH atoms in the PME arrays.
            The background heavy atoms have their multipoles interpolated into
            the foreground outside this class.              */
        setAtoms(esvSystem.getExtendedAtoms(), esvSystem.getExtendedMolecule());
        // Allocate shared derivative storage.
        esvPermRealDeriv_shared = new SharedDouble[numESVs];
        esvPermSelfDeriv_shared = new SharedDouble[numESVs];
        esvPermRecipDeriv_shared = new SharedDouble[numESVs];
        esvInducedRealDeriv_shared = new SharedDouble[numESVs];
        esvInducedSelfDeriv_shared = new SharedDouble[numESVs];
        esvInducedRecipDeriv_shared = new SharedDouble[numESVs];
        for (int i = 0; i < numESVs; i++) {
            esvPermRealDeriv_shared[i] = new SharedDouble(0.0);
            esvPermSelfDeriv_shared[i] = new SharedDouble(0.0);
            esvPermRecipDeriv_shared[i] = new SharedDouble(0.0);
            esvInducedRealDeriv_shared[i] = new SharedDouble(0.0);
            esvInducedSelfDeriv_shared[i] = new SharedDouble(0.0);
            esvInducedRecipDeriv_shared[i] = new SharedDouble(0.0);
        }
        updateEsvLambda();

        logger.info(String.format(" Attached extended system (%d variables) to PME.", numESVs));
        logger.info(String.format("   RecipTerm, ScaleAlpha,Polar:  %6b, %6s, %6s",
                reciprocalSpaceTerm, esvConfig.scaleAlpha, polarization.toString()));
        logger.info(String.format("   PermAlpha,Exp; PolAlpha,Exp:  %6.2f, %6.2f, %6.2f, %6.2f",
                permLambdaAlpha, permLambdaExponent, polLambdaAlpha, polLambdaExponent));

        /**
         * Enforce ESV-handling requirements slash best practices.
         */
        if (doLigandVaporElec || doNoLigandCondensedSCF) {
            logger.warning(format(" (EsvSys) Nonstandard PME-ESV configuration: %b %b\n"
                    + "          Enforcing vaporElec,noLigCondSCF:  %b %b",
                    doLigandVaporElec, doNoLigandCondensedSCF, false, false));
            doLigandVaporElec = false;
            doNoLigandCondensedSCF = false;
        }
    }

    public void detachExtendedSystem() {
        fill(esvAtomsScaled, false);
        esvTerm = false;
        esvSystem = null;
        esvPermRealDeriv_shared = null;
        esvPermSelfDeriv_shared = null;
        esvPermRecipDeriv_shared = null;
        esvInducedRealDeriv_shared = null;
        esvInducedSelfDeriv_shared = null;
        esvInducedRecipDeriv_shared = null;
        numESVs = 0;
    }

    /**
     * Precalculate ESV factors subsequent to lambda propagation.
     */
    public void updateEsvLambda() {
        if (!esvTerm) {
            return;
        }
        // Query ExtendedSystem to create local preloads of all lambda quantities.
        numESVs = esvSystem.size();
        if (esvLambda == null || esvLambda.length < nAtoms) {
            esvAtomsScaled = new boolean[nAtoms];
            esvAtomsScaledAlpha = new boolean[nAtoms];
            esvLambda = new double[nAtoms];
            esvIndex = new Integer[nAtoms];
            fill(esvAtomsScaled, false);
            fill(esvAtomsScaledAlpha, false);
            fill(esvLambda, 1.0);
            fill(esvIndex, null);
        }
        /* Preload components for permanent electrostatics. */
        for (int i = 0; i < nAtoms; i++) {
            esvAtomsScaled[i] = esvSystem.isExtended(i);
            esvAtomsScaledAlpha[i] = esvSystem.isAlphaScaled(i);
            esvLambda[i] = esvSystem.getLambda(i);
            esvIndex[i] = (esvSystem.getEsvIndex(i) != null) ? esvSystem.getEsvIndex(i) : null;
        }

        // For atoms with both a foreground and background multipole, preload interpolated multipole and derivative.
        // Initialize unscaled variables.
        esvMultipoleDot = new double[nSymm][nAtoms][10];
        unscaledPolarizability = new double[nAtoms];
        unscaledField = new double[maxThreads][3][nAtoms];
        unscaledFieldCR = new double[maxThreads][3][nAtoms];
        unscaledDirectDipole = new double[nAtoms][3];
        unscaledDirectDipoleCR = new double[nAtoms][3];
        unscaledInducedDipole = new double[nSymm][nAtoms][3];
        unscaledInducedDipoleCR = new double[nSymm][nAtoms][3];
        unscaledCartDipolePhi = new double[nAtoms][tensorCount];
        unscaledCartDipolePhiCR = new double[nAtoms][tensorCount];
        for (int i = 0; i < nAtoms; i++) {
            final Atom ai = atoms[i];
            if (ai.getPolarizeType() == null) {
                logger.warning("Null polarize type during ESV init.");
                continue;
            }
            unscaledPolarizability[i] = ai.getUnscaledPolarizability();
            polarizability[i] = (esvConfig.scaleAlpha)
                    ? ai.getScaledPolarizability()
                    : ai.getUnscaledPolarizability();
        }
    }

    /**
     * This enables PermanentRealSpaceFieldLoop and RealSpaceEnergyLoop to share
     * a code path for masking rule treatment. Defaults: d11scale = 0.0;
     * (scaleD) p12scale = p13scale = 0.0;	(scaleP) m12scale = m13scale = 0.0;
     * (scale) m14scale = 0.4 amoeba, 1/1.2 amber, 0.5 else; m15scale = 0.8
     * amoeba, 1 else.
     */
    private void applyMaskingRules(boolean set, int i, double[] maskPerm, double[] maskPolD, double[] maskPolP) {
        final Atom ai = atoms[i];
        for (Atom ak : ai.get1_5s()) {
            if (maskPerm != null) {
                maskPerm[ak.getArrayIndex()] = (set) ? m15scale : 1.0;
            }
        }
        for (Torsion torsion : ai.getTorsions()) {
            Atom ak = torsion.get1_4(ai);
            if (ak != null) {
                final int index = ak.getArrayIndex();
                if (maskPerm != null) {
                    maskPerm[index] = (set) ? m14scale : 1.0;
                }
                for (int member : ip11[i]) {
                    if (member == index) {
                        maskPolP[index] = (set) ? intra14Scale : 1.0;
                    }
                }
            }
        }
        for (Angle angle : ai.getAngles()) {
            Atom ak = angle.get1_3(ai);
            if (ak != null) {
                final int index = ak.getArrayIndex();
                if (maskPerm != null) {
                    maskPerm[index] = (set) ? m13scale : 1.0;
                }
                maskPolP[index] = (set) ? p13scale : 1.0;
            }
        }
        for (Bond bond : ai.getBonds()) {
            Atom ak = bond.get1_2(ai);
            final int index = ak.getArrayIndex();
            if (maskPerm != null) {
                maskPerm[index] = (set) ? m12scale : 1.0;
            }
            maskPolP[index] = (set) ? p12scale : 1.0;
        }
        for (int index : ip11[i]) {
            maskPolD[index] = (set) ? d11scale : 1.0;
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

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

    public double getEsvDerivative(int id) {
        return getEsvDeriv_Permanent(id) + getEsvDeriv_Induced(id);
    }

    public double getEsvDeriv_Permanent(int id) {
        return esvPermRealDeriv_shared[id].get()
                + esvPermRecipDeriv_shared[id].get()
                + esvPermSelfDeriv_shared[id].get();
    }

    public double getEsvDeriv_PermReal(int id) {
        return esvPermRealDeriv_shared[id].get();
    }

    public double getEsvDeriv_PermSelf(int id) {
        return esvPermSelfDeriv_shared[id].get();
    }

    public double getEsvDeriv_PermRecip(int id) {
        return esvPermRecipDeriv_shared[id].get();
    }

    public double getEsvDeriv_Induced(int id) {
        return esvInducedRealDeriv_shared[id].get()
                + esvInducedRecipDeriv_shared[id].get()
                + esvInducedSelfDeriv_shared[id].get();
    }

    public double getEsvDeriv_IndReal(int id) {
        return esvInducedRealDeriv_shared[id].get();
    }

    public double getEsvDeriv_IndSelf(int id) {
        return esvInducedSelfDeriv_shared[id].get();
    }

    public double getEsvDeriv_IndRecip(int id) {
        return esvInducedRecipDeriv_shared[id].get();
    }

    public double getEsvDeriv_GK(int id) {
        throw new UnsupportedOperationException();
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
                "ELEC", atoms[i].getIndex(), atoms[i].getAtomType().name,
                atoms[k].getIndex(), atoms[k].getAtomType().name, r, eij));
    }

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
            if (esvTerm) {
                if (reciprocalSpaceTerm && aewald > 0.0) {
                    reciprocalSpace.splineInducedDipoles(unscaledInducedDipole, unscaledInducedDipoleCR, use);
                }
                sectionTeam.execute(inducedDipoleFieldRegion);
                if (reciprocalSpaceTerm && aewald > 0.0) {
                    reciprocalSpace.computeInducedPhi(unscaledCartDipolePhi, unscaledCartDipolePhiCR);
                }
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
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the induced real space field in thread " + getThreadIndex());
                throw ex;
            } catch (Exception ex) {
                String msg = "Fatal exception computing the induced real space field in thread " + getThreadIndex();
                logger.log(Level.SEVERE, msg, ex);
            }
        }

        private class InducedPreconditionerFieldLoop extends IntegerForLoop {

            private double x[], y[], z[];
            private double ind[][], indCR[][];
            private double indUns[][], indUnsCR[][];
            private double fX[], fY[], fZ[];
            private double fXCR[], fYCR[], fZCR[];
            private double fXuns[], fYuns[], fZuns[];
            private double fXCRuns[], fYCRuns[], fZCRuns[];

            public InducedPreconditionerFieldLoop() {
            }

            @Override
            public IntegerSchedule schedule() {
                return realSpaceSchedule;
            }

            @Override
            public void start() {
                int threadIndex = getThreadIndex();
                realSpaceScfTimes[threadIndex] -= System.nanoTime();
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
                if (esvTerm) {
                    indUns = unscaledInducedDipole[0];
                    indUnsCR = unscaledInducedDipoleCR[0];
                    fXuns = unscaledField[threadIndex][0];
                    fYuns = unscaledField[threadIndex][1];
                    fZuns = unscaledField[threadIndex][2];
                    fXCRuns = unscaledFieldCR[threadIndex][0];
                    fYCRuns = unscaledFieldCR[threadIndex][1];
                    fZCRuns = unscaledFieldCR[threadIndex][2];
                    fill(fXuns, 0.0);
                    fill(fYuns, 0.0);
                    fill(fZuns, 0.0);
                    fill(fXCRuns, 0.0);
                    fill(fYCRuns, 0.0);
                    fill(fZCRuns, 0.0);
                }
            }

            @Override
            public void finish() {
                int threadIndex = getThreadIndex();
                realSpaceScfTimes[threadIndex] += System.nanoTime();
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
                    double fxu = 0.0;
                    double fyu = 0.0;
                    double fzu = 0.0;
                    double pxu = 0.0;
                    double pyu = 0.0;
                    double pzu = 0.0;
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
                        {
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
                        if (esvTerm) {	// unscaled
                            final double ukxUns = indUns[k][0];
                            final double ukyUns = indUns[k][1];
                            final double ukzUns = indUns[k][2];
                            final double ukr = ukxUns * xr + ukyUns * yr + ukzUns * zr;
                            final double bn2ukr = bn2 * ukr;
                            final double fimx = -bn1 * ukxUns + bn2ukr * xr;
                            final double fimy = -bn1 * ukyUns + bn2ukr * yr;
                            final double fimz = -bn1 * ukzUns + bn2ukr * zr;
                            final double rr5ukr = rr5 * ukr;
                            final double fidx = -rr3 * ukxUns + rr5ukr * xr;
                            final double fidy = -rr3 * ukyUns + rr5ukr * yr;
                            final double fidz = -rr3 * ukzUns + rr5ukr * zr;
                            fxu += (fimx - fidx);
                            fyu += (fimy - fidy);
                            fzu += (fimz - fidz);
                            final double pkxUns = indUnsCR[k][0];
                            final double pkyUns = indUnsCR[k][1];
                            final double pkzUns = indUnsCR[k][2];
                            final double pkr = pkxUns * xr + pkyUns * yr + pkzUns * zr;
                            final double bn2pkr = bn2 * pkr;
                            final double pimx = -bn1 * pkxUns + bn2pkr * xr;
                            final double pimy = -bn1 * pkyUns + bn2pkr * yr;
                            final double pimz = -bn1 * pkzUns + bn2pkr * zr;
                            final double rr5pkr = rr5 * pkr;
                            final double pidx = -rr3 * pkxUns + rr5pkr * xr;
                            final double pidy = -rr3 * pkyUns + rr5pkr * yr;
                            final double pidz = -rr3 * pkzUns + rr5pkr * zr;
                            pxu += (pimx - pidx);
                            pyu += (pimy - pidy);
                            pzu += (pimz - pidz);
                            final double uir = uix * xr + uiy * yr + uiz * zr;
                            final double bn2uir = bn2 * uir;
                            final double fkmx = -bn1 * uix + bn2uir * xr;
                            final double fkmy = -bn1 * uiy + bn2uir * yr;
                            final double fkmz = -bn1 * uiz + bn2uir * zr;
                            final double rr5uir = rr5 * uir;
                            final double fkdx = -rr3 * uix + rr5uir * xr;
                            final double fkdy = -rr3 * uiy + rr5uir * yr;
                            final double fkdz = -rr3 * uiz + rr5uir * zr;
                            fXuns[k] += (fkmx - fkdx);
                            fYuns[k] += (fkmy - fkdy);
                            fZuns[k] += (fkmz - fkdz);
                            final double pir = pix * xr + piy * yr + piz * zr;
                            final double bn2pir = bn2 * pir;
                            final double pkmx = -bn1 * pix + bn2pir * xr;
                            final double pkmy = -bn1 * piy + bn2pir * yr;
                            final double pkmz = -bn1 * piz + bn2pir * zr;
                            final double rr5pir = rr5 * pir;
                            final double pkdx = -rr3 * pix + rr5pir * xr;
                            final double pkdy = -rr3 * piy + rr5pir * yr;
                            final double pkdz = -rr3 * piz + rr5pir * zr;
                            fXCRuns[k] += (pkmx - pkdx);
                            fYCRuns[k] += (pkmy - pkdy);
                            fZCRuns[k] += (pkmz - pkdz);
                        }	// unscaled
                    }
                    fX[i] += fx;
                    fY[i] += fy;
                    fZ[i] += fz;
                    fXCR[i] += px;
                    fYCR[i] += py;
                    fZCR[i] += pz;
                    if (esvTerm) {
                        fXuns[i] += fxu;
                        fYuns[i] += fyu;
                        fZuns[i] += fzu;
                        fXCRuns[i] += pxu;
                        fYCRuns[i] += pyu;
                        fZCRuns[i] += pzu;
                    }
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
                    final double indsUns[][] = (esvTerm) ? unscaledInducedDipole[iSymm] : null;
                    final double indsUnsCR[][] = (esvTerm) ? unscaledInducedDipoleCR[iSymm] : null;
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
                        double fxu = 0.0;
                        double fyu = 0.0;
                        double fzu = 0.0;
                        double pxu = 0.0;
                        double pyu = 0.0;
                        double pzu = 0.0;
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

                        final double uixUns, uiyUns, uizUns;
                        final double pixUns, piyUns, pizUns;
                        if (esvTerm) {
                            uixUns = indUns[i][0];
                            uiyUns = indUns[i][1];
                            uizUns = indUns[i][2];
                            pixUns = indUnsCR[i][0];
                            piyUns = indUnsCR[i][1];
                            pizUns = indUnsCR[i][2];
                        } else {
                            uixUns = uix;
                            uiyUns = uiy;
                            uizUns = uiz;
                            pixUns = pix;
                            piyUns = piy;
                            pizUns = piz;
                        }
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
                            {
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
                            if (esvTerm) {	// unscaled
                                final double ukxUns = indsUns[k][0];
                                final double ukyUns = indsUns[k][1];
                                final double ukzUns = indsUns[k][2];
                                final double pkxUns = indsUnsCR[k][0];
                                final double pkyUns = indsUnsCR[k][1];
                                final double pkzUns = indsUnsCR[k][2];
                                final double ukr = ukxUns * xr + ukyUns * yr + ukzUns * zr;
                                final double bn2ukr = bn2 * ukr;
                                final double fimx = -bn1 * ukxUns + bn2ukr * xr;
                                final double fimy = -bn1 * ukyUns + bn2ukr * yr;
                                final double fimz = -bn1 * ukzUns + bn2ukr * zr;
                                final double rr5ukr = rr5 * ukr;
                                final double fidx = -rr3 * ukxUns + rr5ukr * xr;
                                final double fidy = -rr3 * ukyUns + rr5ukr * yr;
                                final double fidz = -rr3 * ukzUns + rr5ukr * zr;
                                fxu += selfScale * (fimx - fidx);
                                fyu += selfScale * (fimy - fidy);
                                fzu += selfScale * (fimz - fidz);
                                final double pkr = pkxUns * xr + pkyUns * yr + pkzUns * zr;
                                final double bn2pkr = bn2 * pkr;
                                final double pimx = -bn1 * pkxUns + bn2pkr * xr;
                                final double pimy = -bn1 * pkyUns + bn2pkr * yr;
                                final double pimz = -bn1 * pkzUns + bn2pkr * zr;
                                final double rr5pkr = rr5 * pkr;
                                final double pidx = -rr3 * pkxUns + rr5pkr * xr;
                                final double pidy = -rr3 * pkyUns + rr5pkr * yr;
                                final double pidz = -rr3 * pkzUns + rr5pkr * zr;
                                pxu += selfScale * (pimx - pidx);
                                pyu += selfScale * (pimy - pidy);
                                pzu += selfScale * (pimz - pidz);
                                final double uir = uixUns * xr + uiyUns * yr + uizUns * zr;
                                final double bn2uir = bn2 * uir;
                                final double fkmx = -bn1 * uixUns + bn2uir * xr;
                                final double fkmy = -bn1 * uiyUns + bn2uir * yr;
                                final double fkmz = -bn1 * uizUns + bn2uir * zr;
                                final double rr5uir = rr5 * uir;
                                final double fkdx = -rr3 * uixUns + rr5uir * xr;
                                final double fkdy = -rr3 * uiyUns + rr5uir * yr;
                                final double fkdz = -rr3 * uizUns + rr5uir * zr;
                                final double xcf = selfScale * (fkmx - fkdx);
                                final double ycf = selfScale * (fkmy - fkdy);
                                final double zcf = selfScale * (fkmz - fkdz);
                                fX[k] += (xcf * transOp[0][0] + ycf * transOp[1][0] + zcf * transOp[2][0]);
                                fY[k] += (xcf * transOp[0][1] + ycf * transOp[1][1] + zcf * transOp[2][1]);
                                fZ[k] += (xcf * transOp[0][2] + ycf * transOp[1][2] + zcf * transOp[2][2]);
                                final double pir = pixUns * xr + piyUns * yr + pizUns * zr;
                                final double bn2pir = bn2 * pir;
                                final double pkmx = -bn1 * pixUns + bn2pir * xr;
                                final double pkmy = -bn1 * piyUns + bn2pir * yr;
                                final double pkmz = -bn1 * pizUns + bn2pir * zr;
                                final double rr5pir = rr5 * pir;
                                final double pkdx = -rr3 * pixUns + rr5pir * xr;
                                final double pkdy = -rr3 * piyUns + rr5pir * yr;
                                final double pkdz = -rr3 * pizUns + rr5pir * zr;
                                final double xcp = selfScale * (pkmx - pkdx);
                                final double ycp = selfScale * (pkmy - pkdy);
                                final double zcp = selfScale * (pkmz - pkdz);
                                fXCR[k] += (xcp * transOp[0][0] + ycp * transOp[1][0] + zcp * transOp[2][0]);
                                fYCR[k] += (xcp * transOp[0][1] + ycp * transOp[1][1] + zcp * transOp[2][1]);
                                fZCR[k] += (xcp * transOp[0][2] + ycp * transOp[1][2] + zcp * transOp[2][2]);
                            }	// unscaled
                        }
                        fX[i] += fx;
                        fY[i] += fy;
                        fZ[i] += fz;
                        fXCR[i] += px;
                        fYCR[i] += py;
                        fZCR[i] += pz;
                        if (esvTerm) {
                            fXuns[i] += fxu;
                            fYuns[i] += fyu;
                            fZuns[i] += fzu;
                            fXCRuns[i] += pxu;
                            fYCRuns[i] += pyu;
                            fZCRuns[i] += pzu;
                        }
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
            public void run(int lb, int ub) {
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
        public void run() {
            try {
                int ti = getThreadIndex();
                if (pcgLoop[ti] == null) {
                    pcgLoop[ti] = new PCGInitLoop();
                }
                execute(0, nAtoms - 1, pcgLoop[ti]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the mutual induced dipoles in thread " + getThreadIndex());
                throw ex;
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
            public void run(int lb, int ub) {

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
        public void run() {
            try {
                int ti = getThreadIndex();
                if (pcgLoop[ti] == null) {
                    pcgLoop[ti] = new PCGInitLoop();
                }
                execute(0, nAtoms - 1, pcgLoop[ti]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the mutual induced dipoles in thread " + getThreadIndex());
                throw ex;
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
            public void run(int lb, int ub) {

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
        public void run() {
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
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the mutual induced dipoles in thread " + getThreadIndex());
                throw ex;
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
            public void run(int lb, int ub) {
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
            public void run(int lb, int ub) {
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
        public void run() {
            try {
                int ti = getThreadIndex();
                if (iterLoop1[ti] == null) {
                    iterLoop1[ti] = new PCGIterLoop1();
                    iterLoop2[ti] = new PCGIterLoop2();
                }
                execute(0, nAtoms - 1, iterLoop1[ti]);
                execute(0, nAtoms - 1, iterLoop2[ti]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the mutual induced dipoles in thread " + getThreadIndex());
                throw ex;
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
            public void run(int lb, int ub) {
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
            public void run(int lb, int ub) {
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
        public void run() {
            try {
                int ti = getThreadIndex();
                if (pcgLoop[ti] == null) {
                    pcgLoop[ti] = new PCGLoop();
                }
                execute(0, nAtoms - 1, pcgLoop[ti]);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception computing the mutual induced dipoles in thread " + getThreadIndex());
                throw ex;
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
            public void run(int lb, int ub) {
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
                    if (esvTerm) {
                        fx = 0.0;
                        fy = 0.0;
                        fz = 0.0;
                        fxCR = 0.0;
                        fyCR = 0.0;
                        fzCR = 0.0;
                        for (int j = 1; j < maxThreads; j++) {
                            fx += unscaledField[j][0][i];
                            fy += unscaledField[j][1][i];
                            fz += unscaledField[j][2][i];
                            fxCR += unscaledFieldCR[j][0][i];
                            fyCR += unscaledFieldCR[j][1][i];
                            fzCR += unscaledFieldCR[j][2][i];
                        }
                        unscaledField[0][0][i] += fx;
                        unscaledField[0][1][i] += fy;
                        unscaledField[0][2][i] += fz;
                        unscaledFieldCR[0][0][i] += fxCR;
                        unscaledFieldCR[0][1][i] += fyCR;
                        unscaledFieldCR[0][2][i] += fzCR;
                    }
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
                        if (esvTerm) {
                            unscaledField[0][0][i] += fx;
                            unscaledField[0][1][i] += fy;
                            unscaledField[0][2][i] += fz;
                            unscaledFieldCR[0][0][i] += fxCR;
                            unscaledFieldCR[0][1][i] += fyCR;
                            unscaledFieldCR[0][2][i] += fzCR;
                        }
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
                        // TODO: Add unscaled field handling to GK.
                    }
                }
            }
        }
    }

    @Override
    public double getTotalMultipoleEnergy() {
        return permanentMultipoleEnergy + polarizationEnergy;
    }

    @Override
    public double getPermanentEnergy() {
        return permanentMultipoleEnergy;
    }

    @Override
    public double getPermRealEnergy() {
        return permanentRealSpaceEnergy;
    }

    @Override
    public double getPermSelfEnergy() {
        return permanentSelfEnergy;
    }

    @Override
    public double getPermRecipEnergy() {
        return permanentReciprocalEnergy;
    }

    @Override
    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    @Override
    public double getIndRealEnergy() {
        return inducedRealSpaceEnergy;
    }

    @Override
    public double getIndSelfEnergy() {
        return inducedSelfEnergy;
    }

    @Override
    public double getIndRecipEnergy() {
        return inducedReciprocalEnergy;
    }

    @Override
    public double getGKEnergy() {
        return generalizedKirkwoodEnergy;
    }

    @Override
    public GeneralizedKirkwood getGK() {
        return generalizedKirkwood;
    }

    /**
     * Number of unique tensors for given order.
     */
    private static final int tensorCount = MultipoleTensor.tensorCount(3);
    /**
     * Taylor expansion multiplier for quadrupole interactions.
     */
    private double oneThird = 1.0 / 3.0;
    /**
     ******************************
     * Access methods for OpenMM.
     ******************************
     */
    @Override
    public double[][][] getCoordinates() {
        return coordinates;
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
    @Override
    public Polarization getPolarizationType() {
        return polarization;
    }
    @Override
    public int[][] getAxisAtoms() {
        return axisAtom;
    }
    @Override
    public double getScale14() {
        return m14scale;
    }
    @Override
    public double getEwaldCoefficient() {
        return aewald;
    }
    @Override
    public ReciprocalSpace getReciprocalSpace() {
        return reciprocalSpace;
    }
    @Override
    public ELEC_FORM getElecForm() {
        return elecForm;
    }

    @Override
    public String getName() {
        return "Quasi-internal";
    }

}
