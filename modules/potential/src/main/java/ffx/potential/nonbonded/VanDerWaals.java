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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.AdderDoubleArray;
import ffx.numerics.AtomicDoubleArray;
import ffx.numerics.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.MultiDoubleArray;
import ffx.numerics.PJDoubleArray;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Torsion;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.VDWType;

import static ffx.numerics.AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI;
import static ffx.potential.nonbonded.VanDerWaalsForm.EPS;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADMIN;
import static ffx.potential.parameters.ForceField.ForceFieldString.ARRAY_REDUCTION;
import static ffx.potential.parameters.ForceField.toEnumForm;

/**
 * The van der Waals class computes van der Waals interaction in parallel using
 * a {@link NeighborList} for any {@link Crystal}. The repulsive power (e.g.
 * 12), attractive power (e.g. 6) and buffering (e.g. for the AMOEBA
 * buffered-14-7) can all be specified such that both Lennard-Jones and AMOEBA
 * are supported.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class VanDerWaals implements MaskingInterface,
        LambdaInterface {

    private static final Logger logger = Logger.getLogger(VanDerWaals.class.getName());

    /**
     * This field specifies resolution for multi-scale modeling.
     */
    private Resolution resolution = null;

    /**
     * Boundary conditions and crystal symmetry.
     */
    private Crystal crystal;
    /**
     * An array of all atoms in the system.
     */
    private Atom atoms[];
    private Atom previousAtoms[];
    /**
     * Specification of the molecular index for each atom.
     */
    private int molecule[];
    private int previousMolecule[];
    /**
     * The Force Field that defines the van der Waals interactions.
     */
    private ForceField forceField;
    /**
     * An array of whether each atom in the system should be used in the
     * calculations.
     */
    private boolean use[] = null;
    /**
     * A local convenience variable equal to atoms.length.
     */
    private int nAtoms;
    /**
     * A local convenience variable equal to the number of crystal symmetry
     * operators.
     */
    private int nSymm;
    /**
     * *************************************************************************
     * Lambda variables.
     */
    private boolean gradient;
    private boolean lambdaTerm;
    private boolean esvTerm;
    private boolean isSoft[];

    /**
     * There are 2 softCore arrays of length nAtoms.
     *
     * The first is used for atoms in the outer loop that are hard. This mask
     * equals: false for inner loop hard atoms true for inner loop soft atoms
     *
     * The second is used for atoms in the outer loop that are soft. This mask
     * equals: true for inner loop hard atoms false for inner loop soft atoms
     */
    private boolean softCore[][];
    private boolean softCoreInit;
    private static final byte HARD = 0;
    private static final byte SOFT = 1;

    /**
     * Turn on inter-molecular softcore interactions using molecular index.
     */
    private boolean intermolecularSoftcore = false;
    /**
     * Turn on intra-molecular softcore interactions using molecular index.
     */
    private boolean intramolecularSoftcore = false;
    /**
     * Current value of the lambda state variable.
     */
    private double lambda = 1.0;
    /**
     * Exponent on lambda.
     */
    private double vdwLambdaExponent = 1.0;
    /**
     * Offset in Angstroms.
     */
    private double vdwLambdaAlpha = 0.05;
    /**
     * Polymorphic inner class to set sc1,sc2,dsc1,etc only when necessary.
     * [nThreads]
     */
    private LambdaFactors[] lambdaFactors = null;
    private double sc1 = 0.0;       // alpha * (1 - lambdaProduct)^2
    private double sc2 = 1.0;       // lambdaProduct
    private double dsc1dL = 0.0;
    private double dsc2dL = 0.0;
    private double d2sc1dL2 = 0.0;
    private double d2sc2dL2 = 0.0;
    /**
     * Generalized extended system variables.
     */
    private ExtendedSystem esvSystem;
    private int numESVs = 0;
    private double esvLambda[];
    private double esvLambdaSwitch[];
    private double esvSwitchDeriv[];
    private boolean esvAtoms[];
    private int atomEsvID[];
    /**
     * TODO: To enable multi-dimensional lambda variables. Preload this with the
     * effective (combined) lambda for each atom state. [nAtoms][nStates]
     */
    private double esvStateLambda[][];
    /**
     * TODO: To enable multi-dimensional lambda variables. Preload this with the
     * effective (combined) radEps for each atom state. [nAtoms][nStates]
     */
    private double esvStateRadEps[][];

    /**
     * *************************************************************************
     * Coordinate arrays.
     */
    /**
     * A local copy of atomic coordinates, including reductions on the hydrogen
     * atoms.
     */
    private double coordinates[];
    /**
     * Reduced coordinates of size: [nSymm][nAtoms * 3]
     */
    private double reduced[][];
    private double reducedXYZ[];
    /**
     * Neighbor lists for each atom. Size: [nSymm][nAtoms][nNeighbors]
     */
    private int[][][] neighborLists;
    private static final byte XX = 0;
    private static final byte YY = 1;
    private static final byte ZZ = 2;
    /**
     * *************************************************************************
     * Force field parameters and constants for the Buffered-14-7 potential.
     */
    /**
     * A local reference to the atom class of each atom in the system.
     */
    private int atomClass[];
    /**
     * Hydrogen atom vdW sites are located toward their heavy atom relative to
     * their nucleus. This is a look-up that gives the heavy atom index for each
     * hydrogen.
     */
    private int reductionIndex[];
    private int bondMask[][];
    private int angleMask[][];
    private int torsionMask[][];
    /**
     * Each hydrogen vdW site is located a fraction of the way from the heavy
     * atom nucleus to the hydrogen nucleus (~0.9).
     */
    private double reductionValue[];
    private double longRangeCorrection;
    private final boolean doLongRangeCorrection;
    /**
     * *************************************************************************
     * Parallel variables.
     */
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final IntegerSchedule pairwiseSchedule;
    private final SharedInteger sharedInteractions;
    private final SharedDouble sharedEnergy;
    private final SharedDouble shareddEdL;
    private final SharedDouble sharedd2EdL2;
    private SharedDouble[] esvDeriv;

    private AtomicDoubleArrayImpl atomicDoubleArrayImpl = MULTI;
    /**
     * X-component of the Cartesian coordinate gradient.
     */
    private AtomicDoubleArray gradX;
    /**
     * Y-component of the Cartesian coordinate gradient.
     */
    private AtomicDoubleArray gradY;
    /**
     * Z-component of the Cartesian coordinate gradient.
     */
    private AtomicDoubleArray gradZ;
    /**
     * X-component of the lambda derivative of the Cartesian coordinate
     * gradient.
     */
    private AtomicDoubleArray lambdaGradX;
    /**
     * Y-component of the lambda derivative of the Cartesian coordinate
     * gradient.
     */
    private AtomicDoubleArray lambdaGradY;
    /**
     * Z-component of the lambda derivative of the Cartesian coordinate
     * gradient.
     */
    private AtomicDoubleArray lambdaGradZ;

    /**
     * The neighbor-list includes 1-2 and 1-3 interactions, which are masked out
     * in the van der Waals energy code. The AMOEBA force field includes 1-4
     * interactions fully.
     */
    private NeighborList neighborList;
    private final VanDerWaalsRegion vanDerWaalsRegion;
    private boolean neighborListOnly = true;
    /**
     * Timing variables.
     */
    private boolean print = false;
    private final long initializationTime[];
    private final long vdwTime[];
    private final long reductionTime[];
    private long initializationTotal, vdwTotal, reductionTotal;
    private final VanDerWaalsForm vdwForm;
    private final NonbondedCutoff nonbondedCutoff;
    private final MultiplicativeSwitch multiplicativeSwitch;

    /**
     * The VanDerWaals class constructor.
     *
     * @param atoms the Atom array to do van Der Waals calculations on.
     * @param molecule the molecule number for each atom.
     * @param crystal The boundary conditions.
     * @param forceField the ForceField parameters to apply.
     * @param parallelTeam The parallel environment.
     *
     * @since 1.0
     */
    public VanDerWaals(Atom atoms[], int molecule[], Crystal crystal, ForceField forceField,
            ParallelTeam parallelTeam) {
        this.atoms = atoms;
        this.molecule = molecule;
        this.crystal = crystal;
        this.parallelTeam = parallelTeam;
        this.forceField = forceField;
        nAtoms = atoms.length;
        nSymm = crystal.spaceGroup.getNumberOfSymOps();

        vdwForm = new VanDerWaalsForm(forceField);

        /**
         * Lambda parameters.
         */
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);
        if (lambdaTerm) {
            shareddEdL = new SharedDouble();
            sharedd2EdL2 = new SharedDouble();
            vdwLambdaAlpha = forceField.getDouble(ForceFieldDouble.VDW_LAMBDA_ALPHA, 0.05);
            vdwLambdaExponent = forceField.getDouble(ForceFieldDouble.VDW_LAMBDA_EXPONENT, 1.0);
            if (vdwLambdaAlpha < 0.0) {
                vdwLambdaAlpha = 0.05;
            }
            if (vdwLambdaExponent < 1.0) {
                vdwLambdaExponent = 1.0;
            }
            intermolecularSoftcore = forceField.getBoolean(
                    ForceField.ForceFieldBoolean.INTERMOLECULAR_SOFTCORE, false);
            intramolecularSoftcore = forceField.getBoolean(
                    ForceField.ForceFieldBoolean.INTRAMOLECULAR_SOFTCORE, false);
        } else {
            shareddEdL = null;
            sharedd2EdL2 = null;
        }

        /**
         * Parallel constructs.
         */
        threadCount = parallelTeam.getThreadCount();
        sharedInteractions = new SharedInteger();
        sharedEnergy = new SharedDouble();
        doLongRangeCorrection = forceField.getBoolean(ForceField.ForceFieldBoolean.VDWLRTERM, false);
        vanDerWaalsRegion = new VanDerWaalsRegion();
        initializationTime = new long[threadCount];
        vdwTime = new long[threadCount];
        reductionTime = new long[threadCount];

        /**
         * Define how force arrays will be accumulated.
         */
        atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
        String value = forceField.getString(ARRAY_REDUCTION, "MULTI");
        try {
            atomicDoubleArrayImpl = AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
        } catch (Exception e) {
            logger.info(format(" Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value, atomicDoubleArrayImpl));
        }
        logger.info(format("  Using %s arrays.", atomicDoubleArrayImpl.toString()));

        /**
         * Allocate coordinate arrays and set up reduction indices and values.
         */
        initAtomArrays();

        /**
         * Set up the cutoff and polynomial switch.
         */
        double buff = 2.0;
        double vdwcut;
        if (!crystal.aperiodic()) {
            vdwcut = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        } else {
            vdwcut = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, crystal.a / 2.0 - (buff + 1.0));
            // If aperiodic, set the vdW cutoff to cover everything.
        }

        // Ensure van der Waals cutoff is at least as large as Ewald cutoff.
        double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        if (ewaldOff > vdwcut) {
            vdwcut = ewaldOff;
            logger.info(" The van der Waals cutoff must be at least as large as the Ewald cutoff.");
            logger.info(String.format(" The van der Waals cutoff has been set to %f", ewaldOff));
        }

        /**
         * Define the multiplicative switch, which sets vdW energy to zero at
         * the cutoff distance using a window that begin at 90% of the cutoff
         * distance.
         */
        double vdwtaper = 0.9 * vdwcut;
        double cut = vdwtaper;
        double off = vdwcut;
        nonbondedCutoff = new NonbondedCutoff(off, cut, buff);
        multiplicativeSwitch = new MultiplicativeSwitch(off, cut);

        /**
         * Parallel neighbor list builder.
         */
        neighborList = new NeighborList(null, this.crystal, atoms, off, buff, parallelTeam);
        pairwiseSchedule = neighborList.getPairwiseSchedule();
        neighborLists = new int[nSymm][][];

        /**
         * Reduce and expand the coordinates of the asymmetric unit. Then build
         * the first neighborlist.
         */
        buildNeighborList(atoms);

        logger.info("  Van der Waals");
        logger.info(format("   Switch Start:                         %6.3f (A)", cut));
        logger.info(format("   Cut-Off:                              %6.3f (A)", off));
        //logger.info(format(" Long-Range Correction:                   %B", doLongRangeCorrection));

        if (lambdaTerm) {
            logger.info("  Lambda,ESV Parameters");
            logger.info(format("   Softcore Alpha:                        %5.3f", vdwLambdaAlpha));
            logger.info(format("   Lambda Exponent:                       %5.3f\n", vdwLambdaExponent));
        }
    }

    /**
     * Allocate coordinate arrays and set up reduction indices and values.
     */
    private void initAtomArrays() {
        if (esvTerm) {
            atoms = esvSystem.getAtomsExtH();
            nAtoms = atoms.length;
        }
        if (atomClass == null || nAtoms > atomClass.length
                || lambdaTerm || esvTerm) {
            atomClass = new int[nAtoms];
            coordinates = new double[nAtoms * 3];
            reduced = new double[nSymm][nAtoms * 3];
            reducedXYZ = reduced[0];
            reductionIndex = new int[nAtoms];
            reductionValue = new double[nAtoms];
            bondMask = new int[nAtoms][];
            angleMask = new int[nAtoms][];
            if (vdwForm.vdwType == VanDerWaalsForm.VDW_TYPE.LENNARD_JONES) {
                torsionMask = new int[nAtoms][];
            } else {
                torsionMask = null;
            }
            use = new boolean[nAtoms];
            isSoft = new boolean[nAtoms];
            softCore = new boolean[2][nAtoms];
            lambdaGradX = null;
            lambdaGradY = null;
            lambdaGradZ = null;

            switch (atomicDoubleArrayImpl) {
                case MULTI:
                    gradX = new MultiDoubleArray(threadCount, nAtoms);
                    gradY = new MultiDoubleArray(threadCount, nAtoms);
                    gradZ = new MultiDoubleArray(threadCount, nAtoms);
                    if (lambdaTerm) {
                        lambdaGradX = new MultiDoubleArray(threadCount, nAtoms);
                        lambdaGradY = new MultiDoubleArray(threadCount, nAtoms);
                        lambdaGradZ = new MultiDoubleArray(threadCount, nAtoms);
                    }
                    break;
                case PJ:
                    gradX = new PJDoubleArray(threadCount, nAtoms);
                    gradY = new PJDoubleArray(threadCount, nAtoms);
                    gradZ = new PJDoubleArray(threadCount, nAtoms);
                    if (lambdaTerm) {
                        lambdaGradX = new PJDoubleArray(threadCount, nAtoms);
                        lambdaGradY = new PJDoubleArray(threadCount, nAtoms);
                        lambdaGradZ = new PJDoubleArray(threadCount, nAtoms);
                    }
                    break;
                case ADDER:
                default:
                    gradX = new AdderDoubleArray(threadCount, nAtoms);
                    gradY = new AdderDoubleArray(threadCount, nAtoms);
                    gradZ = new AdderDoubleArray(threadCount, nAtoms);
                    if (lambdaTerm) {
                        lambdaGradX = new AdderDoubleArray(threadCount, nAtoms);
                        lambdaGradY = new AdderDoubleArray(threadCount, nAtoms);
                        lambdaGradZ = new AdderDoubleArray(threadCount, nAtoms);
                    }
                    break;
            }
        }

        /**
         * Initialize all atoms to be used in the energy.
         */
        fill(use, true);
        fill(isSoft, false);
        fill(softCore[HARD], false);
        fill(softCore[SOFT], false);
        softCoreInit = false;

        esvAtoms = new boolean[nAtoms]; // Needs initialized regardless of esvTerm.
        esvLambda = new double[nAtoms];
        atomEsvID = new int[nAtoms];
        fill(esvAtoms, false);
        fill(esvLambda, 1.0);
        fill(atomEsvID, -1);
        if (esvTerm) {
            updateEsvLambda();
        }

        lambdaFactors = new LambdaFactors[threadCount];
        for (int i = 0; i < threadCount; i++) {
            if (esvTerm) {
                lambdaFactors[i] = new LambdaFactorsESV();
            } else if (lambdaTerm) {
                lambdaFactors[i] = new LambdaFactorsOSRW();
            } else {
                lambdaFactors[i] = new LambdaFactors();
            }
        }

        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            assert (i == ai.xyzIndex - 1);
            double xyz[] = ai.getXYZ(null);
            int i3 = i * 3;
            coordinates[i3 + XX] = xyz[XX];
            coordinates[i3 + YY] = xyz[YY];
            coordinates[i3 + ZZ] = xyz[ZZ];
            AtomType atomType = ai.getAtomType();
            if (atomType == null) {
                logger.severe(ai.toString());
                continue;   // Severe no longer guarantees program crash.
            }
            String vdwIndex = forceField.getString(ForceField.ForceFieldString.VDWINDEX, "Class");
            if (vdwIndex.equalsIgnoreCase("Type")) {
                atomClass[i] = atomType.type;
            } else {
                atomClass[i] = atomType.atomClass;
            }
            VDWType type = forceField.getVDWType(Integer.toString(atomClass[i]));
            if (type == null) {
                logger.info(" No vdW type for atom class " + atomClass[i]);
                logger.severe(" No vdW type for atom " + ai.toString());
                return;
            }
            ai.setVDWType(type);
            ArrayList<Bond> bonds = ai.getBonds();
            int numBonds = bonds.size();
            if (type.reductionFactor > 0.0 && numBonds == 1) {
                Bond bond = bonds.get(0);
                Atom heavyAtom = bond.get1_2(ai);
                // Atom indexes start at 1
                reductionIndex[i] = heavyAtom.xyzIndex - 1;
                reductionValue[i] = type.reductionFactor;
            } else {
                reductionIndex[i] = i;
                reductionValue[i] = 0.0;
            }
            bondMask[i] = new int[numBonds];
            for (int j = 0; j < numBonds; j++) {
                Bond bond = bonds.get(j);
                bondMask[i][j] = bond.get1_2(ai).xyzIndex - 1;
            }
            ArrayList<Angle> angles = ai.getAngles();
            int numAngles = 0;
            for (Angle angle : angles) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    numAngles++;
                }
            }
            angleMask[i] = new int[numAngles];
            int j = 0;
            for (Angle angle : angles) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    angleMask[i][j++] = ak.xyzIndex - 1;
                }
            }
            if (vdwForm.scale14 != 1.0) {
                ArrayList<Torsion> torsions = ai.getTorsions();
                int numTorsions = 0;
                for (Torsion torsion : torsions) {
                    Atom ak = torsion.get1_4(ai);
                    if (ak != null) {
                        numTorsions++;
                    }
                }
                torsionMask[i] = new int[numTorsions];
                j = 0;
                for (Torsion torsion : torsions) {
                    Atom ak = torsion.get1_4(ai);
                    if (ak != null) {
                        torsionMask[i][j++] = ak.xyzIndex - 1;
                    }
                }
            }
        }
    }

    public void setResolution(Resolution resolution) {
        this.resolution = resolution;
    }

    public final void buildNeighborList(Atom[] atoms) {
        neighborList.setAtoms(atoms);
        if (esvTerm) {  // TODO: Move ESV neighborlist construction into the parallel team.
            neighborList.buildList(reduced, neighborLists, null, neighborListOnly, true);
        } else {
            neighborListOnly = true;
            try {
                parallelTeam.execute(vanDerWaalsRegion);
            } catch (Exception e) {
                String message = " Fatal exception expanding coordinates.\n";
                logger.log(Level.SEVERE, message, e);
            }
            neighborListOnly = false;
        }
    }

    public void setAtoms(Atom atoms[], int molecule[]) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.molecule = molecule;

        if (nAtoms != molecule.length) {
            logger.warning("Atom and molecule arrays are of different lengths.");
            throw new IllegalArgumentException();
        }
        initAtomArrays();
        buildNeighborList(atoms);
    }

    /**
     * Allow sharing the of the VanDerWaals NeighborList with ParticleMeshEwald.
     *
     * @return The NeighborList.
     */
    public NeighborList getNeighborList() {
        return neighborList;
    }

    /**
     * <p>
     * Getter for the field <code>pairwiseSchedule</code>.</p>
     *
     * @return a {@link edu.rit.pj.IntegerSchedule} object.
     */
    public IntegerSchedule getPairwiseSchedule() {
        return pairwiseSchedule;
    }

    private double getLongRangeCorrection() {
        if (true) {     // Current implementation of VdW-LR is not to be trusted.
            throw new UnsupportedOperationException();
        }
        if (esvTerm) {  // Need to treat esvLambda chain terms below before you can do this.
            throw new UnsupportedOperationException();
        }

        int maxClass = vdwForm.maxClass;

        /**
         * Long range correction.
         */
        int radCount[] = new int[maxClass + 1];
        int softRadCount[] = new int[maxClass + 1];
        for (int i = 0; i < maxClass; i++) {
            radCount[i] = 0;
            softRadCount[i] = 0;
        }

        for (int i = 0; i < nAtoms; i++) {
            radCount[atomClass[i]]++;
            if (isSoft[i]) {
                softRadCount[atomClass[i]]++;
            }
        }

        /**
         * Integrate to maxR = 60 Angstroms or ~20 sigma. Integration step size
         * of delR to be 0.01 Angstroms.
         */
        double maxR = 60.0;
        int n = (int) (100.0 * (maxR - nonbondedCutoff.cut));
        double delR = (maxR - nonbondedCutoff.cut) / n;
        double total = 0.0;
        /*
         * logger.info(String.format(" Long range correction integral: Steps %d,
         * Step Size %8.3f, Window %8.3f-%8.3f", n, delR, cut, cut + delR * n));
         */
        /**
         * Loop over vdW types.
         */
        for (int i = 0; i < maxClass; i++) {
            for (int j = 0; j < maxClass; j++) {
                int j2 = j * 2;
                double irv = vdwForm.radEps[i][j2 + vdwForm.RADMIN];
                double ev = vdwForm.radEps[i][j2 + vdwForm.EPS];
                double sume = 0.0;
                for (int k = 0; k <= n; k++) {
                    double r = nonbondedCutoff.cut + k * delR;
                    double r2 = r * r;
                    final double rho = r * irv;
                    final double rho3 = rho * rho * rho;
                    final double rho7 = rho3 * rho3 * rho;
                    final double rhod = rho + vdwForm.delta;
                    final double rhod3 = rhod * rhod * rhod;
                    final double rhod7 = rhod3 * rhod3 * rhod;
                    final double t1 = vdwForm.t1n / rhod7;
                    final double t2 = vdwForm.gamma1 / rho7 + vdwForm.gamma;
                    final double eij = ev * t1 * (t2 - 2.0);
                    /**
                     * Apply one minus the multiplicative switch if the
                     * interaction distance is less than the end of the
                     * switching window.
                     */
                    double taper = 1.0;
                    if (r2 < nonbondedCutoff.off2) {
                        double r3 = r * r2;
                        double r4 = r2 * r2;
                        double r5 = r2 * r3;
                        taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                        taper = 1.0 - taper;
                    }
                    double jacobian = 4.0 * PI * r2;
                    double e = jacobian * eij * taper;
                    if (k != 0 && k != n) {
                        sume += e;
                    } else {
                        sume += 0.5 * e;
                    }
                }
                double trapezoid = delR * sume;
                // Normal correction
                total += radCount[i] * radCount[j] * trapezoid;
                // Correct for softCore vdW that are being turned off.
                // TODO add accounting for esvLambda softcoring
                if (lambda < 1.0) {
                    total -= (softRadCount[i] * radCount[j]
                            + (radCount[i] - softRadCount[i]) * softRadCount[j])
                            * (1.0 - lambda) * trapezoid;
                }
            }
        }
        /**
         * Note the factor of 2 to avoid double counting.
         */
        return total / 2;
    }

    /**
     * Get the total van der Waals potential energy.
     *
     * @return The energy.
     * @since 1.0
     */
    public double getEnergy() {
        return sharedEnergy.get();
    }

    /**
     * Get the number of interacting pairs.
     *
     * @return The interaction count.
     * @since 1.0
     */
    public int getInteractions() {
        return sharedInteractions.get();
    }

    /**
     * Get the buffer size.
     *
     * @return The buffer.
     * @since 1.0
     */
    public double getBuffer() {
        return nonbondedCutoff.buff;
    }

    /**
     * The energy routine may be called repeatedly.
     *
     * @param gradient If true, gradients with respect to atomic coordinates are
     * computed.
     * @param print If true, there is verbose printing.
     * @return The energy.
     * @since 1.0
     */
    public double energy(boolean gradient, boolean print) {
        this.gradient = gradient;
        this.print = print;

        try {
            parallelTeam.execute(vanDerWaalsRegion);
        } catch (Exception e) {
            String message = " Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
        }

        return sharedEnergy.get();
    }

    /**
     * {@inheritDoc}
     *
     * Apply masking rules for 1-2 and 1-3 interactions.
     */
    @Override
    public void applyMask(final double mask[], final int i) {
        if (vdwForm.scale14 != 1.0) {
            final int[] torsionMaski = torsionMask[i];
            final int n14 = torsionMaski.length;
            for (int m = 0; m < n14; m++) {
                mask[torsionMaski[m]] = vdwForm.scale14;
            }
        }
        final int[] angleMaski = angleMask[i];
        final int n13 = angleMaski.length;
        for (int m = 0; m < n13; m++) {
            mask[angleMaski[m]] = vdwForm.scale13;
        }
        final int[] bondMaski = bondMask[i];
        final int n12 = bondMaski.length;
        for (int m = 0; m < n12; m++) {
            mask[bondMaski[m]] = vdwForm.scale12;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Remove the masking rules for 1-2 and 1-3 interactions.
     */
    @Override
    public void removeMask(final double mask[], final int i) {
        if (vdwForm.scale14 != 1.0) {
            final int[] torsionMaski = torsionMask[i];
            final int n14 = torsionMaski.length;
            for (int m = 0; m < n14; m++) {
                mask[torsionMaski[m]] = 1.0;
            }
        }
        final int[] angleMaski = angleMask[i];
        final int n13 = angleMaski.length;
        for (int m = 0; m < n13; m++) {
            mask[angleMaski[m]] = 1.0;
        }
        final int[] bondMaski = bondMask[i];
        final int n12 = bondMaski.length;
        for (int m = 0; m < n12; m++) {
            mask[bondMaski[m]] = 1.0;
        }
    }

    /**
     * <p>
     * Getter for the field <code>neighborLists</code>.</p>
     *
     * @return an array of int.
     */
    public int[][][] getNeighborLists() {
        return neighborLists;
    }

    /**
     * Log the van der Waals interaction.
     *
     * @param i Atom i.
     * @param k Atom j.
     * @param r The distance rij.
     * @param eij The interaction energy.
     * @since 1.0
     */
    private void log(int i, int k, double r, double eij) {
        final Atom ai = atoms[i];
        final Atom ak = atoms[k];
        int classi = ai.getAtomType().atomClass;
        int classk = ak.getAtomType().atomClass;
        double combined = 1.0 / vdwForm.radEps[classi][classk * 2 + VanDerWaalsForm.RADMIN];
        logger.info(format("%s %6d-%s %6d-%s %10.4f  %10.4f  %10.4f",
                "VDW", atoms[i].xyzIndex, atoms[i].getAtomType().name,
                atoms[k].xyzIndex, atoms[k].getAtomType().name,
                combined, r, eij));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        assert (lambda >= 0.0 && lambda <= 1.0);
        this.lambda = lambda;
        sc1 = vdwLambdaAlpha * (1.0 - lambda) * (1.0 - lambda);
        dsc1dL = -2.0 * vdwLambdaAlpha * (1.0 - lambda);
        d2sc1dL2 = 2.0 * vdwLambdaAlpha;
        sc2 = pow(lambda, vdwLambdaExponent);
        dsc2dL = vdwLambdaExponent * pow(lambda, vdwLambdaExponent - 1.0);
        if (vdwLambdaExponent >= 2.0) {
            d2sc2dL2 = vdwLambdaExponent * (vdwLambdaExponent - 1.0) * pow(lambda, vdwLambdaExponent - 2.0);
        } else {
            d2sc2dL2 = 0.0;
        }

        /**
         * If LambdaFactors are in OSRW mode, update them now.
         */
        if (!esvTerm) {
            for (LambdaFactors lf : lambdaFactors) {
                lf.setFactors();
            }
        }

        initSoftCore(false);

        // Redo the long range correction.
        if (doLongRangeCorrection) {
            longRangeCorrection = getLongRangeCorrection();
            logger.info(String.format(" Long-range vdW correction %12.8f (kcal/mole).",
                    longRangeCorrection));
        } else {
            longRangeCorrection = 0.0;
        }
    }

    /**
     * VdW version should get only the ExtH version of ExtendedSystem lists.
     * This is equivalent to the mola atom array on loading the fully-protonated
     * system. Note: we assume that heavy-atom radii do not differ between
     * protonation states; this is violated only by {Cys-SG, Asp-OD[12], and
     * Glu-OD[12]} (plus the bugged/missing Am'13_Tyr-OH).
     */
    public void updateEsvLambda() {
        if (!esvTerm) {
            logger.warning("Improper method call: updateEsvLambda().");
            return;
        }
        numESVs = esvSystem.n();
        if (esvLambdaSwitch == null || esvLambdaSwitch.length < nAtoms) {
            esvLambdaSwitch = new double[nAtoms];
            esvSwitchDeriv = new double[nAtoms];
            atomEsvID = new int[nAtoms];
            fill(esvLambdaSwitch, 1.0);
            fill(esvSwitchDeriv, 0.0);
            fill(atomEsvID, -1);
        }
        for (int i = 0; i < nAtoms; i++) {
            if (esvSystem.isExtH(i)) {
                esvAtoms[i] = true;
                esvLambda[i] = esvSystem.exthLambda(i);
                esvLambdaSwitch[i] = esvSystem.exthEsv(i).getLambdaSwitch();
                esvSwitchDeriv[i] = esvSystem.exthEsv(i).getSwitchDeriv();
                atomEsvID[i] = esvSystem.exthEsvId(i);
            }
        }
        if (esvDeriv == null || esvDeriv.length < numESVs) {
            esvDeriv = new SharedDouble[numESVs];
            for (int i = 0; i < numESVs; i++) {
                esvDeriv[i] = new SharedDouble(0.0);
            }
        }
        initSoftCore(true);
        // Call to long-range correction here, when it's trustworthy.
    }

    /**
     * The trick: The setFactors(i,k) method is called every time through the
     * inner VdW loop, avoiding an "if (esv)" branch statement. A plain OSRW run
     * will have an object of type LambdaFactorsOSRW instead, which contains an
     * empty version of setFactors(i,k). The OSRW version sets new factors only
     * on lambda updates, in setLambda().
     */
    public class LambdaFactors {

        protected double sc1 = 0.0;
        protected double dsc1dL = 0.0;
        protected double d2sc1dL2 = 0.0;
        protected double sc2 = 1.0;
        protected double dsc2dL = 0.0;
        protected double d2sc2dL2 = 0.0;

        /**
         * Overriden by the OSRW version which updates only during setLambda().
         */
        public void setFactors() {
        }

        /**
         * Overriden by the ESV version which updates with every softcore
         * interaction.
         */
        public void setFactors(int i, int k) {
        }
    }

    public class LambdaFactorsOSRW extends LambdaFactors {

        @Override
        public void setFactors() {
            sc1 = VanDerWaals.this.sc1;
            dsc1dL = VanDerWaals.this.dsc1dL;
            d2sc1dL2 = VanDerWaals.this.d2sc1dL2;
            sc2 = VanDerWaals.this.sc2;
            dsc2dL = VanDerWaals.this.dsc2dL;
            d2sc2dL2 = VanDerWaals.this.d2sc2dL2;
        }
    }

    public class LambdaFactorsESV extends LambdaFactors {

        @Override
        public void setFactors(int i, int k) {
            final double esvLambdaProduct = esvLambda[i] * esvLambda[k] * lambda;
            sc1 = vdwLambdaAlpha * (1.0 - esvLambdaProduct) * (1.0 - esvLambdaProduct);
            dsc1dL = -2.0 * vdwLambdaAlpha * (1.0 - esvLambdaProduct);
            d2sc1dL2 = 2.0 * vdwLambdaAlpha;
            sc2 = esvLambdaProduct;
            dsc2dL = 1.0;
            d2sc2dL2 = 0.0;
        }
    }

    private void initSoftCore(boolean rebuild) {
        /**
         * Initialize the softcore atom masks.
         */
        if (!softCoreInit || rebuild) {
            for (int i = 0; i < nAtoms; i++) {
                isSoft[i] = atoms[i].applyLambda();
                if (esvTerm && esvSystem.isExtH(i)) {
                    isSoft[i] = true;
                }
                if (isSoft[i]) {
                    // Outer loop atom hard, inner loop atom soft.
                    softCore[HARD][i] = true;
                    // Both soft: full intramolecular ligand interactions.
                    softCore[SOFT][i] = false;
                } else {
                    // Both hard: full interaction between atoms.
                    softCore[HARD][i] = false;
                    // Outer loop atom soft, inner loop atom hard.
                    softCore[SOFT][i] = true;
                }
            }
            softCoreInit = true;
        }
    }

    public void attachExtendedSystem(ExtendedSystem system) {
        if (system == null) {
            logger.severe("Tried to attach null extended system.");
        }
        esvTerm = true;
        esvSystem = system;
        numESVs = esvSystem.n();

        // Launch shared lambda/esvLambda initializers if missed (ie. !lambdaTerm) in constructor.
        vdwLambdaAlpha = forceField.getDouble(ForceFieldDouble.VDW_LAMBDA_ALPHA, 0.05);
        vdwLambdaExponent = forceField.getDouble(ForceFieldDouble.VDW_LAMBDA_EXPONENT, 1.0);
        if (vdwLambdaExponent != 1.0) {
            logger.warning(format("ESVs are compatible only with a vdwLambdaExponent of unity!"
                    + " (found %g, resetting to 1.0)", vdwLambdaExponent));
            vdwLambdaExponent = 1.0;
        }
        if (vdwLambdaAlpha < 0.0) {
            vdwLambdaAlpha = 0.05;
        }
        if (vdwLambdaExponent < 1.0) {
            vdwLambdaExponent = 1.0;
        }
        intermolecularSoftcore = forceField.getBoolean(
                ForceField.ForceFieldBoolean.INTERMOLECULAR_SOFTCORE, false);
        intramolecularSoftcore = forceField.getBoolean(
                ForceField.ForceFieldBoolean.INTRAMOLECULAR_SOFTCORE, false);

        previousAtoms = atoms;
        previousMolecule = molecule;
        Atom[] atomsExt = esvSystem.getAtomsExtH();
        int[] moleculeExt = esvSystem.getMoleculeExtH();
        setAtoms(atomsExt, moleculeExt);
        updateEsvLambda();
    }

    public void detachExtendedSystem() {
        setAtoms(previousAtoms, molecule);
        fill(esvAtoms, false);
        fill(esvLambda, 1.0);
        esvTerm = false;
        esvSystem = null;
        esvDeriv = null;
        numESVs = 0;
        initSoftCore(true); // To remove entries from isSoft[] that were due only to ESVs.
    }

    public void setIntermolecularSoftcore(boolean intermolecularSoftcore) {
        if (!(lambdaTerm || esvTerm)) {
            logger.warning("Illegal softcoring.");
            throw new IllegalArgumentException();
        }
        this.intermolecularSoftcore = intermolecularSoftcore;
    }

    public void setIntramolecularSoftcore(boolean intramolecularSoftcore) {
        if (!(lambdaTerm || esvTerm)) {
            logger.warning("Illegal softcoring.");
            throw new IllegalArgumentException();
        }
        this.intramolecularSoftcore = intramolecularSoftcore;
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
     */
    @Override
    public double getdEdL() {
        if (shareddEdL == null || !lambdaTerm) {
            return 0.0;
        }
        return shareddEdL.get();
    }

    public double[] getdEdEsv() {
        if (!esvTerm || esvSystem == null) {
            logger.warning("Suspicious call to non-existent ESV derivative.");
        }
        double[] dEdEsv = new double[numESVs];
        for (int i = 0; i < numESVs; i++) {
            dEdEsv[i] = esvDeriv[i].get();
        }
        return dEdEsv;
    }

    public double getdEdEsv(int esvID) {
        return esvDeriv[esvID].get();
    }

    /**
     * {@inheritDoc}
     *
     * @param lambdaGradient the lambda Gradient array (dU/dL/dX).
     */
    @Override
    public void getdEdXdL(double[] lambdaGradient) {
        if (lambdaGradX == null || !lambdaTerm) {
            return;
        }
        int index = 0;
        // double lgx[] = lambdaGradX[0];
        // double lgy[] = lambdaGradY[0];
        // double lgz[] = lambdaGradZ[0];
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                lambdaGradient[index++] += lambdaGradX.get(i);
                lambdaGradient[index++] += lambdaGradY.get(i);
                lambdaGradient[index++] += lambdaGradZ.get(i);
                // lambdaGradient[index++] += lgx[i];
                // lambdaGradient[index++] += lgy[i];
                // lambdaGradient[index++] += lgz[i];
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        if (sharedd2EdL2 == null || !lambdaTerm) {
            return 0.0;
        }
        return sharedd2EdL2.get();
    }

    /**
     * If the crystal being passed in is not equal to the current crystal, then
     * some van der Waals data structures may need to updated. If
     * <code>nSymm</code> has changed, update arrays dimensioned by nSymm.
     * Finally, rebuild the neighbor-lists.
     *
     * @param crystal The new crystal instance defining the symmetry and
     * boundary conditions.
     */
    public void setCrystal(Crystal crystal) {
        this.crystal = crystal;
        int newNSymm = crystal.spaceGroup.getNumberOfSymOps();
        if (nSymm != newNSymm) {
            nSymm = newNSymm;
            /**
             * Allocate memory if necessary.
             */
            if (reduced == null || reduced.length < nSymm) {
                reduced = new double[nSymm][nAtoms * 3];
                reducedXYZ = reduced[0];
                neighborLists = new int[nSymm][][];
            }
        }
        neighborList.setCrystal(crystal);
        neighborListOnly = true;
        try {
            print = false;
            parallelTeam.execute(vanDerWaalsRegion);
        } catch (Exception e) {
            String message = " Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    public void destroy() throws Exception {
        if (neighborList != null) {
            neighborList.destroy();
        }
    }

    /**
     * Test if both atoms match the set Resolution (or true when unset).
     */
    private boolean include(Atom atom1, Atom atom2) {
        return ((resolution == null)
                || (atom1.getResolution() == resolution && atom2.getResolution() == resolution));
    }

    private class VanDerWaalsRegion extends ParallelRegion {

        private final InitializationLoop initializationLoop[];
        private final ExpandLoop expandLoop[];
        private final VanDerWaalsLoop vanDerWaalsLoop[];
        private final ReductionLoop reductionLoop[];

        public VanDerWaalsRegion() {
            initializationLoop = new InitializationLoop[threadCount];
            expandLoop = new ExpandLoop[threadCount];
            vanDerWaalsLoop = new VanDerWaalsLoop[threadCount];
            reductionLoop = new ReductionLoop[threadCount];
        }

        /**
         * {@inheritDoc}
         *
         * This is method should not be called; it is invoked by Parallel Java.
         *
         * @since 1.0
         */
        @Override
        public void start() throws IOException {
            /**
             * Initialize the shared variables.
             */
            if (doLongRangeCorrection) {
                sharedEnergy.set(longRangeCorrection);
            } else {
                sharedEnergy.set(0.0);
            }
            sharedInteractions.set(0);
            if (lambdaTerm) {
                shareddEdL.set(0.0);
                sharedd2EdL2.set(0.0);
            }
            if (esvTerm) {
                for (int i = 0; i < numESVs; i++) {
                    esvDeriv[i].set(0.0);
                }
                lambdaFactors = new LambdaFactorsESV[threadCount];
                for (int i = 0; i < threadCount; i++) {
                    lambdaFactors[i] = new LambdaFactorsESV();
                }
            }

            gradX.alloc(nAtoms);
            gradY.alloc(nAtoms);
            gradZ.alloc(nAtoms);
            if (lambdaTerm) {
                lambdaGradX.alloc(nAtoms);
                lambdaGradY.alloc(nAtoms);
                lambdaGradZ.alloc(nAtoms);
            }
        }

        @Override
        public void finish() throws IOException {
            neighborListOnly = false;
        }

        @Override
        public void run() throws Exception {
            int threadIndex = getThreadIndex();

            /**
             * Locally initialize the Loops to help with NUMA?
             */
            if (initializationLoop[threadIndex] == null) {
                initializationLoop[threadIndex] = new InitializationLoop();
                expandLoop[threadIndex] = new ExpandLoop();
                vanDerWaalsLoop[threadIndex] = new VanDerWaalsLoop();
                reductionLoop[threadIndex] = new ReductionLoop();
            }

            /**
             * Initialize and expand coordinates.
             */
            try {
                if (threadIndex == 0) {
                    initializationTotal = -System.nanoTime();
                }
                execute(0, nAtoms - 1, initializationLoop[threadIndex]);
                execute(0, nAtoms - 1, expandLoop[threadIndex]);
                if (threadIndex == 0) {
                    initializationTotal += System.nanoTime();
                }
            } catch (Exception e) {
                String message = "Fatal exception expanding coordinates in thread: " + threadIndex + "\n";
                logger.log(Level.SEVERE, message, e);
            }

            /**
             * Build the neighbor-list (if necessary) using reduced coordinates.
             */
            if (threadIndex == 0) {
                neighborList.buildList(reduced, neighborLists, null, neighborListOnly, false);
            }
            barrier();

            if (neighborListOnly) {
                return;
            }

            /**
             * Compute van der Waals energy and gradient.
             */
            try {
                if (threadIndex == 0) {
                    vdwTotal = -System.nanoTime();
                }
                execute(0, nAtoms - 1, vanDerWaalsLoop[threadIndex]);
                if (threadIndex == 0) {
                    vdwTotal += System.nanoTime();
                }
            } catch (Exception e) {
                String message = "Fatal exception evaluating van der Waals energy in thread: " + threadIndex + "\n";
                logger.log(Level.SEVERE, message, e);
            }

            /**
             * Reduce derivatives.
             */
            if (gradient || lambdaTerm) {
                try {
                    if (threadIndex == 0) {
                        reductionTotal = -System.nanoTime();
                    }
                    execute(0, nAtoms - 1, reductionLoop[threadIndex]);
                    if (threadIndex == 0) {
                        reductionTotal += System.nanoTime();
                    }
                } catch (Exception e) {
                    String message = "Fatal exception reducing van der Waals gradient in thread: " + threadIndex + "\n";
                    logger.log(Level.SEVERE, message, e);
                }
            }

            /**
             * Log timings.
             */
            if (threadIndex == 0 && logger.isLoggable(Level.FINE)) {
                double total = (initializationTotal + vdwTotal + reductionTotal) * 1e-9;
                logger.fine(format("\n van der Waals: %7.4f (sec)", total));
                logger.fine(" Thread    Init    Energy  Reduce  Total     Counts");
                long initMax = 0;
                long vdwMax = 0;
                long reductionMax = 0;
                long initMin = Long.MAX_VALUE;
                long vdwMin = Long.MAX_VALUE;
                long reductionMin = Long.MAX_VALUE;
                int countMin = Integer.MAX_VALUE;
                int countMax = 0;
                for (int i = 0; i < threadCount; i++) {
                    int count = vanDerWaalsLoop[i].getCount();
                    long totalTime = initializationTime[i] + vdwTime[i] + reductionTime[i];
                    logger.fine(format("    %3d   %7.4f %7.4f %7.4f %7.4f %10d",
                            i, initializationTime[i] * 1e-9, vdwTime[i] * 1e-9,
                            reductionTime[i] * 1e-9, totalTime * 1e-9, count));
                    initMax = max(initializationTime[i], initMax);
                    vdwMax = max(vdwTime[i], vdwMax);
                    reductionMax = max(reductionTime[i], reductionMax);
                    countMax = max(countMax, count);
                    initMin = min(initializationTime[i], initMin);
                    vdwMin = min(vdwTime[i], vdwMin);
                    reductionMin = min(reductionTime[i], reductionMin);
                    countMin = min(countMin, count);
                }
                long totalMin = initMin + vdwMin + reductionMin;
                long totalMax = initMax + vdwMax + reductionMax;
                long totalActual = initializationTotal + vdwTotal + reductionTotal;
                logger.fine(format(" Min      %7.4f %7.4f %7.4f %7.4f %10d",
                        initMin * 1e-9, vdwMin * 1e-9,
                        reductionMin * 1e-9, totalMin * 1e-9, countMin));
                logger.fine(format(" Max      %7.4f %7.4f %7.4f %7.4f %10d",
                        initMax * 1e-9, vdwMax * 1e-9,
                        reductionMax * 1e-9, totalMax * 1e-9, countMax));
                logger.fine(format(" Delta    %7.4f %7.4f %7.4f %7.4f %10d",
                        (initMax - initMin) * 1e-9, (vdwMax - vdwMin) * 1e-9,
                        (reductionMax - reductionMin) * 1e-9, (totalMax - totalMin) * 1e-9,
                        (countMax - countMin)));
                logger.fine(format(" Actual   %7.4f %7.4f %7.4f %7.4f %10d\n",
                        initializationTotal * 1e-9, vdwTotal * 1e-9,
                        reductionTotal * 1e-9, totalActual * 1e-9, sharedInteractions.get()));
            }
        }

        /**
         * Update the local coordinate array and initialize reduction variables.
         */
        private class InitializationLoop extends IntegerForLoop {

            private int threadID;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                threadID = getThreadIndex();
                initializationTime[threadID] = -System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb, i3 = 3 * lb; i <= ub; i++, i3 += 3) {
                    Atom atom = atoms[i];
                    coordinates[i3 + XX] = atom.getX();
                    coordinates[i3 + YY] = atom.getY();
                    coordinates[i3 + ZZ] = atom.getZ();
                    use[i] = atom.getUse();
                }
                if (gradient) {
                    gradX.reset(threadID, lb, ub);
                    gradY.reset(threadID, lb, ub);
                    gradZ.reset(threadID, lb, ub);
                }
                if (lambdaTerm) {
                    lambdaGradX.reset(threadID, lb, ub);
                    lambdaGradY.reset(threadID, lb, ub);
                    lambdaGradZ.reset(threadID, lb, ub);
                }
            }
        }

        private class ExpandLoop extends IntegerForLoop {

            private int threadID;
            private final double in[] = new double[3];
            private final double out[] = new double[3];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                threadID = getThreadIndex();
            }

            @Override
            public void finish() {
                initializationTime[threadID] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Set up the local coordinate array for the asymmetric unit,
                 * applying reduction factors to the hydrogen atoms.
                 */
                for (int i = lb; i <= ub; i++) {
                    int i3 = i * 3;
                    int iX = i3 + XX;
                    int iY = i3 + YY;
                    int iZ = i3 + ZZ;
                    double x = coordinates[iX];
                    double y = coordinates[iY];
                    double z = coordinates[iZ];
                    int redIndex = reductionIndex[i];
                    if (redIndex >= 0) {
                        int r3 = redIndex * 3;
                        double rx = coordinates[r3++];
                        double ry = coordinates[r3++];
                        double rz = coordinates[r3];
                        double a = reductionValue[i];
                        reducedXYZ[iX] = a * (x - rx) + rx;
                        reducedXYZ[iY] = a * (y - ry) + ry;
                        reducedXYZ[iZ] = a * (z - rz) + rz;
                        double[] rxyz = {reducedXYZ[iX], reducedXYZ[iY], reducedXYZ[iZ]};
                        atoms[i].setRedXYZ(rxyz);
                    } else {
                        reducedXYZ[iX] = x;
                        reducedXYZ[iY] = y;
                        reducedXYZ[iZ] = z;
                    }
                }

                List<SymOp> symOps = crystal.spaceGroup.symOps;

                if (symOps.size() != nSymm) {
                    logger.info(String.format(" Programming Error: nSymm %d != symOps.size %d", nSymm, symOps.size()));
                    logger.log(Level.INFO, " Replicates\n{0}", crystal.toString());
                    logger.log(Level.INFO, " Unit Cell\n{0}", crystal.getUnitCell().toString());
                }

                double sp2 = crystal.getSpecialPositionCutoff();
                sp2 *= sp2;
                for (int iSymOp = 1; iSymOp < nSymm; iSymOp++) {
                    SymOp symOp = symOps.get(iSymOp);
                    double xyz[] = reduced[iSymOp];
                    for (int i = lb; i <= ub; i++) {
                        int i3 = i * 3;
                        int iX = i3 + XX;
                        int iY = i3 + YY;
                        int iZ = i3 + ZZ;
                        in[0] = reducedXYZ[iX];
                        in[1] = reducedXYZ[iY];
                        in[2] = reducedXYZ[iZ];
                        crystal.applySymOp(in, out, symOp);
                        xyz[iX] = out[0];
                        xyz[iY] = out[1];
                        xyz[iZ] = out[2];

                        /**
                         * Check if the atom is at a special position.
                         */
                        double dx = in[0] - out[0];
                        double dy = in[1] - out[1];
                        double dz = in[2] - out[2];
                        double r2 = dx * dx + dy * dy + dz * dz;
                        if (r2 < sp2) {
                            logger.log(Level.WARNING, " Atom may be at a special position: {0}", atoms[i].toString());
                        }
                    }
                }
            }
        }

        /**
         * The van der Waals loop class contains methods and thread local
         * variables used to evaluate the van der Waals energy and gradients
         * with respect to atomic coordinates.
         *
         * @author Michael J. Schnieders
         * @since 1.0
         */
        private class VanDerWaalsLoop extends IntegerForLoop {

            private int count;
            private double energy;
            private int threadID;
            private double dEdL;
            private double d2EdL2;
            private double mask[];
            private final double dx_local[];
            private final double transOp[][];
            private LambdaFactors lambdaFactorsLocal;

            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public VanDerWaalsLoop() {
                super();
                dx_local = new double[3];
                transOp = new double[3][3];
            }

            public int getCount() {
                return count;
            }

            @Override
            public IntegerSchedule schedule() {
                return pairwiseSchedule;
            }

            @Override
            public void start() {
                threadID = getThreadIndex();
                vdwTime[threadID] = -System.nanoTime();
                energy = 0.0;
                count = 0;
                if (lambdaTerm) {
                    dEdL = 0.0;
                    d2EdL2 = 0.0;
                }
                lambdaFactorsLocal = lambdaFactors[threadID];
                if (lambdaFactorsLocal == null) {
                    System.exit(1);
                }
                if (mask == null || mask.length < nAtoms) {
                    mask = new double[nAtoms];
                    fill(mask, 1.0);
                }
            }

            @Override
            public void finish() {
                /**
                 * Reduce the energy, interaction count and gradients from this
                 * thread into the shared variables.
                 */
                sharedEnergy.addAndGet(energy);
                sharedInteractions.addAndGet(count);
                if (lambdaTerm) {
                    shareddEdL.addAndGet(dEdL);
                    sharedd2EdL2.addAndGet(d2EdL2);
                }
                vdwTime[threadID] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                double e = 0.0;
                double xyzS[] = reduced[0];
                int list[][] = neighborLists[0];        // neighborLists array: [nSymm][nAtoms][nNeighbors]
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    Atom atomi = atoms[i];
                    final boolean esvi = esvAtoms[i];
                    final int idxi = atomEsvID[i];
                    int i3 = i * 3;
                    final double xi = reducedXYZ[i3++];
                    final double yi = reducedXYZ[i3++];
                    final double zi = reducedXYZ[i3];
                    final int redi = reductionIndex[i];
                    final double redv = reductionValue[i];
                    final double rediv = 1.0 - redv;
                    final int classi = atomClass[i];
                    final double radEpsi[] = vdwForm.radEps[classi];
                    final int moleculei = molecule[i];
                    double gxi = 0.0;
                    double gyi = 0.0;
                    double gzi = 0.0;
                    double gxredi = 0.0;
                    double gyredi = 0.0;
                    double gzredi = 0.0;
                    double lxi = 0.0;
                    double lyi = 0.0;
                    double lzi = 0.0;
                    double lxredi = 0.0;
                    double lyredi = 0.0;
                    double lzredi = 0.0;
                    double localEsvDerivI = 0.0;
                    applyMask(mask, i);
                    // Default is that the outer loop atom is hard.
                    boolean softCorei[] = softCore[HARD];
                    if (isSoft[i]) {
                        softCorei = softCore[SOFT];
                    }
                    /**
                     * Loop over the neighbor list.
                     */
                    final int neighbors[] = list[i];
                    final int npair = neighbors.length;

                    for (int j = 0; j < npair; j++) {
                        final int k = neighbors[j];
                        Atom atomk = atoms[k];
                        if (!use[k] || !include(atomi, atomk)) {
                            continue;
                        }
                        final boolean esvk = esvAtoms[k];
                        final int idxk = atomEsvID[k];
                        // Hide these global variable names for thread safety.
                        final double sc1, dsc1dL, d2sc1dL2;
                        final double sc2, dsc2dL, d2sc2dL2;
                        int k3 = k * 3;
                        final double xk = xyzS[k3++];
                        final double yk = xyzS[k3++];
                        final double zk = xyzS[k3];
                        dx_local[0] = xi - xk;
                        dx_local[1] = yi - yk;
                        dx_local[2] = zi - zk;
                        final double r2 = crystal.image(dx_local);
                        int a2 = atomClass[k] * 2;
                        final double irv = radEpsi[a2 + RADMIN];
                        if (r2 <= nonbondedCutoff.off2 && mask[k] > 0 && irv > 0) {
                            final double r = sqrt(r2);
                            boolean sameMolecule = (moleculei == molecule[k]);
                            boolean soft = softCorei[k]
                                    || (intermolecularSoftcore && !sameMolecule)
                                    || (intramolecularSoftcore && sameMolecule)
                                    || esvi || esvk;
                            /**
                             * The setFactors(i,k) method is empty unless ESVs
                             * are present. If OSRW lambda present,
                             * lambdaFactors will already have been updated
                             * during setLambda().
                             */
                            if (soft) {
                                lambdaFactorsLocal.setFactors(i, k);
                                sc1 = lambdaFactorsLocal.sc1;
                                dsc1dL = lambdaFactorsLocal.dsc1dL;
                                d2sc1dL2 = lambdaFactorsLocal.d2sc1dL2;
                                sc2 = lambdaFactorsLocal.sc2;
                                dsc2dL = lambdaFactorsLocal.dsc2dL;
                                d2sc2dL2 = lambdaFactorsLocal.d2sc2dL2;
                            } else {
                                sc1 = 0.0;
                                dsc1dL = 0.0;
                                d2sc1dL2 = 0.0;
                                sc2 = 1.0;
                                dsc2dL = 0.0;
                                d2sc2dL2 = 0.0;
                            }
                            final double alpha = sc1;
                            final double lambda5 = sc2;
                            /**
                             * Calculate van der Waals interaction energy.
                             * Notation of Schnieders et al. The structure,
                             * thermodynamics, and solubility of organic
                             * crystals from simulation with a polarizable force
                             * field. J. Chem. Theory Comput. 8, 1721–1736
                             * (2012).
                             */
                            final double ev = mask[k] * radEpsi[a2 + EPS];
                            final double eps_lambda = ev * lambda5;
                            final double rho = r * irv;
                            final double rhoDisp1 = vdwForm.rhoDisp1(rho);
                            final double rhoDisp = rhoDisp1 * rho;
                            final double rhoDelta1 = vdwForm.rhoDelta1(rho + vdwForm.delta);
                            final double rhoDelta = rhoDelta1 * (rho + vdwForm.delta);
                            final double alphaRhoDelta = alpha + rhoDelta;
                            final double alphaRhoDispGamma = alpha + rhoDisp + vdwForm.gamma;
                            final double t1d = 1.0 / alphaRhoDelta;
                            final double t2d = 1.0 / alphaRhoDispGamma;
                            final double t1 = vdwForm.t1n * t1d;
                            final double t2a = vdwForm.gamma1 * t2d;
                            final double t2 = t2a - 2.0;
                            double eik = eps_lambda * t1 * t2;
                            /**
                             * Apply a multiplicative switch if the interaction
                             * distance is greater than the beginning of the
                             * taper.
                             */
                            double taper = 1.0;
                            double dtaper = 0.0;
                            if (r2 > nonbondedCutoff.cut2) {
                                final double r3 = r2 * r;
                                final double r4 = r2 * r2;
                                final double r5 = r2 * r3;
                                taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                                dtaper = multiplicativeSwitch.dtaper(r, r2, r3, r4);
                            }
                            eik *= taper;
                            final double eik_preswitch = eik;
                            if (esvi || esvk) {
                                eik *= esvLambdaSwitch[i] * esvLambdaSwitch[k];
                            }
                            e += eik;
                            count++;
                            if (!gradient && !soft) {
                                continue;
                            }
                            final int redk = reductionIndex[k];
                            final double red = reductionValue[k];
                            final double redkv = 1.0 - red;
                            final double dt1d_dr = vdwForm.repDispPower * rhoDelta1 * irv;
                            final double dt2d_dr = vdwForm.dispersivePower * rhoDisp1 * irv;
                            final double dt1_dr = t1 * dt1d_dr * t1d;
                            final double dt2_dr = t2a * dt2d_dr * t2d;
                            final double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);
                            final double ir = 1.0 / r;
                            final double drdx = dx_local[0] * ir;
                            final double drdy = dx_local[1] * ir;
                            final double drdz = dx_local[2] * ir;
                            if (gradient) {
                                final double dswitch = (eik * dtaper + dedr * taper);
                                final double dedx = dswitch * drdx;
                                final double dedy = dswitch * drdy;
                                final double dedz = dswitch * drdz;
                                gxi += dedx * redv;
                                gyi += dedy * redv;
                                gzi += dedz * redv;
                                gxredi += dedx * rediv;
                                gyredi += dedy * rediv;
                                gzredi += dedz * rediv;
                                gradX.sub(threadID, k, red * dedx);
                                gradY.sub(threadID, k, red * dedy);
                                gradZ.sub(threadID, k, red * dedz);
                                gradX.sub(threadID, redk, redkv * dedx);
                                gradY.sub(threadID, redk, redkv * dedy);
                                gradZ.sub(threadID, redk, redkv * dedz);
                            }
                            if (gradient && soft) {
                                final double dt1 = -t1 * t1d * dsc1dL;
                                final double dt2 = -t2a * t2d * dsc1dL;
                                final double f1 = dsc2dL * t1 * t2;
                                final double f2 = sc2 * dt1 * t2;
                                final double f3 = sc2 * t1 * dt2;
                                final double dedl = ev * (f1 + f2 + f3);
                                dEdL += dedl * taper;
                                final double t1d2 = -dsc1dL * t1d * t1d;
                                final double t2d2 = -dsc1dL * t2d * t2d;
                                final double d2t1 = -dt1 * t1d * dsc1dL - t1 * t1d * d2sc1dL2 - t1 * t1d2 * dsc1dL;
                                final double d2t2 = -dt2 * t2d * dsc1dL - t2a * t2d * d2sc1dL2 - t2a * t2d2 * dsc1dL;
                                final double df1 = d2sc2dL2 * t1 * t2 + dsc2dL * dt1 * t2 + dsc2dL * t1 * dt2;
                                final double df2 = dsc2dL * dt1 * t2 + sc2 * d2t1 * t2 + sc2 * dt1 * dt2;
                                final double df3 = dsc2dL * t1 * dt2 + sc2 * dt1 * dt2 + sc2 * t1 * d2t2;
                                final double de2dl2 = ev * (df1 + df2 + df3);
                                d2EdL2 += de2dl2 * taper;
                                final double t11 = -dsc2dL * t2 * dt1_dr;
                                final double t21 = -dsc2dL * t1 * dt2_dr;
                                final double t13 = 2.0 * sc2 * t2 * dt1_dr * dsc1dL * t1d;
                                final double t23 = 2.0 * sc2 * t1 * dt2_dr * dsc1dL * t2d;
                                final double t12 = -sc2 * dt2 * dt1_dr;
                                final double t22 = -sc2 * dt1 * dt2_dr;
                                final double dedldr = ev * (t11 + t12 + t13 + t21 + t22 + t23);
                                final double dswitch = dedl * dtaper + dedldr * taper;
                                final double dedldx = dswitch * drdx;
                                final double dedldy = dswitch * drdy;
                                final double dedldz = dswitch * drdz;
                                lxi += dedldx * redv;
                                lyi += dedldy * redv;
                                lzi += dedldz * redv;
                                lxredi += dedldx * rediv;
                                lyredi += dedldy * rediv;
                                lzredi += dedldz * rediv;
                                if (lambdaTerm) {
                                    lambdaGradX.sub(threadID, k, red * dedldx);
                                    lambdaGradY.sub(threadID, k, red * dedldy);
                                    lambdaGradZ.sub(threadID, k, red * dedldz);
                                    lambdaGradX.sub(threadID, redk, redkv * dedldx);
                                    lambdaGradY.sub(threadID, redk, redkv * dedldy);
                                    lambdaGradZ.sub(threadID, redk, redkv * dedldz);
                                }
                                if (esvi || esvk) {
                                    // Assign this gradient to attached ESVs.
                                    final double dedlp = dedl * taper;
                                    if (esvi) {
                                        final double dlpdli = esvLambda[k] * lambda;
                                        final double dEsvI = dedlp * dlpdli;
                                        // d[S*E] = S'E + E'S
                                        final double dSwEsvI
                                                = esvSwitchDeriv[i] * esvLambdaSwitch[k] * eik_preswitch
                                                + dEsvI * esvLambdaSwitch[i] * esvLambdaSwitch[k];
                                        localEsvDerivI += dSwEsvI;
                                    }
                                    if (esvk) {
                                        final double dlpdlk = esvLambda[i] * lambda;
                                        final double dEsvK = dedlp * dlpdlk;
                                        // d[S*E] = S'E + E'S
                                        final double dSwEsvK
                                                = esvLambdaSwitch[i] * esvSwitchDeriv[k] * eik_preswitch
                                                + dEsvK * esvLambdaSwitch[i] * esvLambdaSwitch[k];
                                        esvDeriv[idxk].addAndGet(dSwEsvK);
                                    }
                                }
                            }
                        }
                    }
                    if (gradient) {
                        gradX.add(threadID, i, gxi);
                        gradY.add(threadID, i, gyi);
                        gradZ.add(threadID, i, gzi);
                        gradX.add(threadID, redi, gxredi);
                        gradY.add(threadID, redi, gyredi);
                        gradZ.add(threadID, redi, gzredi);
                        if (lambdaTerm) {
                            lambdaGradX.add(threadID, i, lxi);
                            lambdaGradY.add(threadID, i, lyi);
                            lambdaGradZ.add(threadID, i, lzi);
                            lambdaGradX.add(threadID, redi, lxredi);
                            lambdaGradY.add(threadID, redi, lyredi);
                            lambdaGradZ.add(threadID, redi, lzredi);
                        }
                        if (esvi) {
                            esvDeriv[idxi].addAndGet(localEsvDerivI);
                        }
                    }
                    removeMask(mask, i);
                }
                energy += e;
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (int iSymOp = 1; iSymOp < nSymm; iSymOp++) {
                    e = 0.0;
                    SymOp symOp = symOps.get(iSymOp);
                    /**
                     * Compute the total transformation operator: R = ToCart *
                     * Rot * ToFrac.
                     */
                    crystal.getTransformationOperator(symOp, transOp);
                    xyzS = reduced[iSymOp];
                    list = neighborLists[iSymOp];
                    for (int i = lb; i <= ub; i++) {
                        int i3 = i * 3;
                        if (!use[i]) {
                            continue;
                        }
                        Atom atomi = atoms[i];
                        final boolean esvi = esvAtoms[i];
                        final int idxi = atomEsvID[i];
                        final double xi = reducedXYZ[i3++];
                        final double yi = reducedXYZ[i3++];
                        final double zi = reducedXYZ[i3];
                        final int redi = reductionIndex[i];
                        final double redv = reductionValue[i];
                        final double rediv = 1.0 - redv;
                        final int classi = atomClass[i];
                        final double radEpsi[] = vdwForm.radEps[classi];
                        double gxi = 0.0;
                        double gyi = 0.0;
                        double gzi = 0.0;
                        double gxredi = 0.0;
                        double gyredi = 0.0;
                        double gzredi = 0.0;
                        double lxi = 0.0;
                        double lyi = 0.0;
                        double lzi = 0.0;
                        double lxredi = 0.0;
                        double lyredi = 0.0;
                        double lzredi = 0.0;
                        double localEsvDerivI = 0.0;
                        // Default is that the outer loop atom is hard.
                        boolean softCorei[] = softCore[HARD];
                        if (isSoft[i]) {
                            softCorei = softCore[SOFT];
                        }
                        /**
                         * Loop over the neighbor list.
                         */
                        final int neighbors[] = list[i];
                        final int npair = neighbors.length;
                        for (int j = 0; j < npair; j++) {
                            final int k = neighbors[j];
                            Atom atomk = atoms[k];
                            if (!use[k] || !include(atomi, atomk)) {
                                continue;
                            }
                            final boolean esvk = esvAtoms[k];
                            final int idxk = atomEsvID[k];
                            // Hide these global variable names for thread safety.
                            final double sc1, dsc1dL, d2sc1dL2;
                            final double sc2, dsc2dL, d2sc2dL2;
                            int k3 = k * 3;
                            final double xk = xyzS[k3++];
                            final double yk = xyzS[k3++];
                            final double zk = xyzS[k3];
                            dx_local[0] = xi - xk;
                            dx_local[1] = yi - yk;
                            dx_local[2] = zi - zk;
                            final double r2 = crystal.image(dx_local);
                            int a2 = atomClass[k] * 2;
                            final double irv = radEpsi[a2 + RADMIN];
                            if (r2 <= nonbondedCutoff.off2 && irv > 0) {
                                final double selfScale = (i == k) ? 0.5 : 1.0;
                                final double r = sqrt(r2);
                                boolean soft = isSoft[i] || softCorei[k]
                                        || esvi || esvk;
                                if (soft) {
                                    lambdaFactorsLocal.setFactors(i, k);
                                    sc1 = lambdaFactorsLocal.sc1;
                                    dsc1dL = lambdaFactorsLocal.dsc1dL;
                                    d2sc1dL2 = lambdaFactorsLocal.d2sc1dL2;
                                    sc2 = lambdaFactorsLocal.sc2;
                                    dsc2dL = lambdaFactorsLocal.dsc2dL;
                                    d2sc2dL2 = lambdaFactorsLocal.d2sc2dL2;
                                } else {
                                    sc1 = 0.0;
                                    dsc1dL = 0.0;
                                    d2sc1dL2 = 0.0;
                                    sc2 = 1.0;
                                    dsc2dL = 0.0;
                                    d2sc2dL2 = 0.0;
                                }
                                final double alpha = sc1;
                                final double lambda5 = sc2;
                                final double ev = radEpsi[a2 + EPS];
                                final double eps_lambda = ev * lambda5;
                                final double rho = r * irv;
                                // final double rhoDisp1 = pow(rho, vdwForm.dispersivePower1);
                                final double rhoDisp1 = vdwForm.rhoDisp1(rho);
                                final double rhoDisp = rhoDisp1 * rho;
                                // final double rhoDelta1 = pow(rho + vdwForm.delta, vdwForm.repDispPower1);
                                final double rhoDelta1 = vdwForm.rhoDelta1(rho + vdwForm.delta);
                                final double rhoDelta = rhoDelta1 * (rho + vdwForm.delta);
                                final double alphaRhoDelta = alpha + rhoDelta;
                                final double alphaRhoDispGamma = alpha + rhoDisp + vdwForm.gamma;
                                final double t1d = 1.0 / alphaRhoDelta;
                                final double t2d = 1.0 / alphaRhoDispGamma;
                                final double t1 = vdwForm.t1n * t1d;
                                final double t2a = vdwForm.gamma1 * t2d;
                                final double t2 = t2a - 2.0;
                                double eik = eps_lambda * t1 * t2;
                                /**
                                 * Apply a multiplicative switch if the
                                 * interaction distance is greater than the
                                 * beginning of the taper.
                                 */
                                double taper = 1.0;
                                double dtaper = 0.0;
                                if (r2 > nonbondedCutoff.cut2) {
                                    final double r3 = r2 * r;
                                    final double r4 = r2 * r2;
                                    final double r5 = r2 * r3;
                                    taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                                    dtaper = multiplicativeSwitch.dtaper(r, r2, r3, r4);
                                }
                                final double eik_preswitch = eik;
                                if (esvi || esvk) {
                                    eik *= esvLambdaSwitch[i] * esvLambdaSwitch[k];
                                }
                                e += selfScale * eik * taper;
                                count++;
                                if (!gradient && !soft) {
                                    continue;
                                }
                                final int redk = reductionIndex[k];
                                final double red = reductionValue[k];
                                final double redkv = 1.0 - red;
                                final double dt1d_dr = vdwForm.repDispPower * rhoDelta1 * irv;
                                final double dt2d_dr = vdwForm.dispersivePower * rhoDisp1 * irv;
                                final double dt1_dr = t1 * dt1d_dr * t1d;
                                final double dt2_dr = t2a * dt2d_dr * t2d;
                                double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);
                                final double ir = 1.0 / r;
                                double drdx = dx_local[0] * ir;
                                double drdy = dx_local[1] * ir;
                                double drdz = dx_local[2] * ir;
                                dedr = (eik * dtaper + dedr * taper);
                                if (gradient) {
                                    double dedx = selfScale * dedr * drdx;
                                    double dedy = selfScale * dedr * drdy;
                                    double dedz = selfScale * dedr * drdz;
                                    gxi += dedx * redv;
                                    gyi += dedy * redv;
                                    gzi += dedz * redv;
                                    gxredi += dedx * rediv;
                                    gyredi += dedy * rediv;
                                    gzredi += dedz * rediv;
                                    /**
                                     * Apply the transpose of the transformation
                                     * operator.
                                     */
                                    final double dedxk = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
                                    final double dedyk = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
                                    final double dedzk = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];
                                    gradX.sub(threadID, k, red * dedxk);
                                    gradY.sub(threadID, k, red * dedyk);
                                    gradZ.sub(threadID, k, red * dedzk);
                                    gradX.sub(threadID, redk, redkv * dedxk);
                                    gradY.sub(threadID, redk, redkv * dedyk);
                                    gradZ.sub(threadID, redk, redkv * dedzk);
                                }
                                if (gradient && soft) {
                                    double dt1 = -t1 * t1d * dsc1dL;
                                    double dt2 = -t2a * t2d * dsc1dL;
                                    double f1 = dsc2dL * t1 * t2;
                                    double f2 = sc2 * dt1 * t2;
                                    double f3 = sc2 * t1 * dt2;
                                    final double dedl = ev * (f1 + f2 + f3);
                                    dEdL += selfScale * dedl * taper;
                                    double t1d2 = -dsc1dL * t1d * t1d;
                                    double t2d2 = -dsc1dL * t2d * t2d;
                                    double d2t1 = -dt1 * t1d * dsc1dL - t1 * t1d * d2sc1dL2 - t1 * t1d2 * dsc1dL;
                                    double d2t2 = -dt2 * t2d * dsc1dL - t2a * t2d * d2sc1dL2 - t2a * t2d2 * dsc1dL;
                                    double df1 = d2sc2dL2 * t1 * t2 + dsc2dL * dt1 * t2 + dsc2dL * t1 * dt2;
                                    double df2 = dsc2dL * dt1 * t2 + sc2 * d2t1 * t2 + sc2 * dt1 * dt2;
                                    double df3 = dsc2dL * t1 * dt2 + sc2 * dt1 * dt2 + sc2 * t1 * d2t2;
                                    double de2dl2 = ev * (df1 + df2 + df3);
                                    d2EdL2 += selfScale * de2dl2 * taper;
                                    double t11 = -dsc2dL * t2 * dt1_dr;
                                    double t12 = -sc2 * dt2 * dt1_dr;
                                    double t13 = 2.0 * sc2 * t2 * dt1_dr * dsc1dL * t1d;
                                    double t21 = -dsc2dL * t1 * dt2_dr;
                                    double t22 = -sc2 * dt1 * dt2_dr;
                                    double t23 = 2.0 * sc2 * t1 * dt2_dr * dsc1dL * t2d;
                                    double dedldr = ev * (t11 + t12 + t13 + t21 + t22 + t23);
                                    dedldr = dedl * dtaper + dedldr * taper;
                                    double dedldx = selfScale * dedldr * drdx;
                                    double dedldy = selfScale * dedldr * drdy;
                                    double dedldz = selfScale * dedldr * drdz;
                                    lxi += dedldx * redv;
                                    lyi += dedldy * redv;
                                    lzi += dedldz * redv;
                                    lxredi += dedldx * rediv;
                                    lyredi += dedldy * rediv;
                                    lzredi += dedldz * rediv;
                                    /**
                                     * Apply the transpose of the transformation
                                     * operator.
                                     */
                                    final double dedldxk = dedldx * transOp[0][0] + dedldy * transOp[1][0] + dedldz * transOp[2][0];
                                    final double dedldyk = dedldx * transOp[0][1] + dedldy * transOp[1][1] + dedldz * transOp[2][1];
                                    final double dedldzk = dedldx * transOp[0][2] + dedldy * transOp[1][2] + dedldz * transOp[2][2];
                                    lambdaGradX.sub(threadID, k, red * dedldxk);
                                    lambdaGradY.sub(threadID, k, red * dedldyk);
                                    lambdaGradZ.sub(threadID, k, red * dedldzk);
                                    lambdaGradX.sub(threadID, redk, redkv * dedldxk);
                                    lambdaGradY.sub(threadID, redk, redkv * dedldyk);
                                    lambdaGradZ.sub(threadID, redk, redkv * dedldzk);
                                    if (esvi || esvk) {
                                        // Assign this gradient to attached ESVs.
                                        final double dedlp = selfScale * dedl * taper;
                                        if (esvi) {
                                            final double dlpdli = esvLambda[k] * lambda;
                                            final double dEsvI = dedlp * dlpdli;
                                            // d[S*E] = S'E + E'S
                                            final double dSwEsvI
                                                    = esvSwitchDeriv[i] * esvLambdaSwitch[k] * eik_preswitch
                                                    + dEsvI * esvLambdaSwitch[i] * esvLambdaSwitch[k];
                                            localEsvDerivI += dSwEsvI;
                                        }
                                        if (esvk) {
                                            final double dlpdlk = esvLambda[i] * lambda;
                                            final double dEsvK = dedlp * dlpdlk;
                                            // d[S*E] = S'E + E'S
                                            final double dSwEsvK
                                                    = esvLambdaSwitch[i] * esvSwitchDeriv[k] * eik_preswitch
                                                    + dEsvK * esvLambdaSwitch[i] * esvLambdaSwitch[k];
                                            esvDeriv[idxk].addAndGet(dSwEsvK);
                                        }
                                    }
                                }
                            }
                        }
                        if (gradient) {
                            gradX.add(threadID, i, gxi);
                            gradY.add(threadID, i, gyi);
                            gradZ.add(threadID, i, gzi);
                            gradX.add(threadID, redi, gxredi);
                            gradY.add(threadID, redi, gyredi);
                            gradZ.add(threadID, redi, gzredi);
                            if (lambdaTerm) {
                                lambdaGradX.add(threadID, i, lxi);
                                lambdaGradY.add(threadID, i, lyi);
                                lambdaGradZ.add(threadID, i, lzi);
                                lambdaGradX.add(threadID, redi, lxredi);
                                lambdaGradY.add(threadID, redi, lyredi);
                                lambdaGradZ.add(threadID, redi, lzredi);
                            }
                            if (esvi) {
                                esvDeriv[idxi].addAndGet(localEsvDerivI);
                            }
                        }
                    }
                    energy += e;
                }
            }
        }

        /**
         * Reduce van der Waals gradient.
         */
        private class ReductionLoop extends IntegerForLoop {

            int threadID;

            @Override
            public void start() {
                threadID = getThreadIndex();
                reductionTime[threadID] = -System.nanoTime();
            }

            @Override
            public void finish() {
                reductionTime[threadID] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                if (gradient) {
                    gradX.reduce(lb, ub);
                    gradY.reduce(lb, ub);
                    gradZ.reduce(lb, ub);
                    for (int i = lb; i <= ub; i++) {
                        Atom ai = atoms[i];
                        ai.addToXYZGradient(gradX.get(i), gradY.get(i), gradZ.get(i));
                    }
                }
                if (lambdaTerm) {
                    lambdaGradX.reduce(lb, ub);
                    lambdaGradY.reduce(lb, ub);
                    lambdaGradZ.reduce(lb, ub);
                }
            }
        }
    }

}
