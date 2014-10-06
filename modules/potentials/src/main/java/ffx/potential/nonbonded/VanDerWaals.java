/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
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
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Torsion;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.VDWType;

/**
 * The van der Waals class computes the buffered 14-7 van der Waals interaction
 * used by the AMOEBA force field in parallel using a {@link NeighborList} for
 * any {@link Crystal}.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class VanDerWaals implements MaskingInterface,
        LambdaInterface {

    private static final Logger logger = Logger.getLogger(VanDerWaals.class.getName());

    public enum VDW_FORM {

        BUFFERED_14_7, LENNARD_JONES_6_12
    };

    /**
     * MolecularAssembly
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * Crystal parameters.
     */
    private Crystal crystal;
    /**
     * An array of all atoms in the system.
     */
    private Atom[] atoms;
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
    private boolean lambdaTerm;
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
    private int molecule[];
    private boolean intermolecularSoftcore = false;
    private double lambda = 1.0;
    private double vdwLambdaExponent = 1.0;
    private double vdwLambdaAlpha = 0.05;
    private double sc1 = 0.0;
    private double sc2 = 1.0;
    private double dsc1dL = 0.0;
    private double dsc2dL = 0.0;
    private double d2sc1dL2 = 0.0;
    private double d2sc2dL2 = 0.0;
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
    private final double radEps[][];
    private double longRangeCorrection;
    private final boolean doLongRangeCorrection;
    private int maxClass;
    private static final byte RADMIN = 0;
    private static final byte EPS = 1;
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
    private boolean gradient;
    /**
     * X-component of the Cartesian coordinate gradient. Size:
     * [threadCount][nAtoms]
     */
    private double gradX[][];
    /**
     * Y-component of the Cartesian coordinate gradient. Size:
     * [threadCount][nAtoms]
     */
    private double gradY[][];
    /**
     * Z-component of the Cartesian coordinate gradient. Size:
     * [threadCount][nAtoms]
     */
    private double gradZ[][];
    /**
     * X-component of the lambda derivative of the Cartesian coordinate
     * gradient. Size: [threadCount][nAtoms]
     */
    private double lambdaGradX[][];
    /**
     * Y-component of the lambda derivative of the Cartesian coordinate
     * gradient. Size: [threadCount][nAtoms]
     */
    private double lambdaGradY[][];
    /**
     * Z-component of the lambda derivative of the Cartesian coordinate
     * gradient. Size: [threadCount][nAtoms]
     */
    private double lambdaGradZ[][];
    /**
     * The neighbor-list includes 1-2 and 1-3 interactions, which are masked out
     * in the van der Waals energy code. The AMOEBA force field includes 1-4
     * interactions fully.
     */
    private final NeighborList neighborList;
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

    private final MultiplicativeSwitch multiplicativeSwitch;

    private VDW_FORM vdwForm = VDW_FORM.BUFFERED_14_7;

    /**
     * The VanDerWaals class constructor.
     *
     * @param molecularAssembly The MolecularAssembly to compute the van der
     * Waals energy of.
     * @param crystal The boundary conditions.
     * @param parallelTeam The parallel environment.
     * @param vdwForm
     * @since 1.0
     */
    public VanDerWaals(MolecularAssembly molecularAssembly, Crystal crystal,
            ParallelTeam parallelTeam, VDW_FORM vdwForm) {
        this.molecularAssembly = molecularAssembly;
        this.atoms = molecularAssembly.getAtomArray();
        this.crystal = crystal;
        this.parallelTeam = parallelTeam;
        this.vdwForm = vdwForm;

        nAtoms = atoms.length;
        nSymm = crystal.spaceGroup.getNumberOfSymOps();

        /**
         * Configure van der Waals well shape parameters.
         */
        switch (vdwForm) {
            case LENNARD_JONES_6_12:
                repulsivePower = 12;
                dispersivePower = 6;
                this.delta = 0.0;
                this.gamma = 0.0;
                break;
            case BUFFERED_14_7:
            default:
                repulsivePower = 14;
                dispersivePower = 7;
                this.delta = 0.07;
                this.gamma = 0.12;
                break;
        }

        repDispPower = repulsivePower - dispersivePower;
        dispersivePower1 = dispersivePower - 1;
        repDispPower1 = repDispPower - 1;
        delta1 = 1.0 + delta;
        t1n = pow(delta1, dispersivePower);
        gamma1 = 1.0 + gamma;

        ForceField forceField = molecularAssembly.getForceField();
        Map<String, VDWType> map = forceField.getVDWTypes();
        TreeMap<String, VDWType> vdwTypes = new TreeMap<>(map);
        maxClass = 0;
        for (VDWType vdwType : vdwTypes.values()) {
            if (vdwType.atomClass > maxClass) {
                maxClass = vdwType.atomClass;
            }
        }
        radEps = new double[maxClass + 1][2 * (maxClass + 1)];
        /**
         * Atom Class numbering starts at 1.
         */
        double twosix = 1.122462048309372981;
        for (VDWType vdwi : vdwTypes.values()) {
            int i = vdwi.atomClass;
            double ri = 0.5 * vdwi.radius;
            double ri2 = ri * ri;
            double ri3 = ri * ri2;
            double e1 = vdwi.wellDepth;
            double se1 = sqrt(e1);
            for (VDWType vdwj : vdwTypes.tailMap(vdwi.getKey()).values()) {
                int j = vdwj.atomClass;
                double rj = 0.5 * vdwj.radius;
                double rj2 = rj * rj;
                double rj3 = rj * rj2;
                double e2 = vdwj.wellDepth;
                double se2 = sqrt(e2);

                double radmin;
                double eps;

                switch (vdwForm) {
                    case LENNARD_JONES_6_12:
                        /**
                         * Geometric
                         */
                        radmin = 2.0 * sqrt(ri * twosix) * sqrt(rj * twosix);
                        /**
                         * Geometric
                         */
                        eps = se1 * se2;
                        break;
                    default:
                    case BUFFERED_14_7:
                        /**
                         * Cubic-mean.
                         */
                        radmin = 2.0 * (ri3 + rj3) / (ri2 + rj2);
                        /**
                         * HHG
                         */
                        eps = 4.0 * (e1 * e2) / ((se1 + se2) * (se1 + se2));
                        break;
                }

                radEps[i][j * 2 + RADMIN] = 1.0 / radmin;
                radEps[i][j * 2 + EPS] = eps;
                radEps[j][i * 2 + RADMIN] = 1.0 / radmin;
                radEps[j][i * 2 + EPS] = eps;
            }
        }

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
         * Allocate coordinate arrays and set up reduction indices and values.
         */
        initAtomArrays();

        /**
         * Set up the cutoff and polynomial switch.
         */
        double vdwcut;
        if (!crystal.aperiodic()) {
            vdwcut = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        } else {
            vdwcut = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, crystal.a / 2.0 - 3.0);
        }
        double vdwtaper = 0.9 * vdwcut;
        cut = vdwtaper;
        off = vdwcut;
        buff = 2.0;
        cut2 = cut * cut;
        off2 = off * off;
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
        try {
            parallelTeam.execute(vanDerWaalsRegion);
        } catch (Exception e) {
            String message = " Fatal exception executing van Der Waals Region.\n";
            logger.log(Level.SEVERE, message, e);
        }

        logger.info("  Van der Waals");
        logger.info(format("   Switch Start:                         %6.3f (A)", cut));
        logger.info(format("   Cut-Off:                              %6.3f (A)", off));
        //logger.info(format(" Long-Range Correction:                   %B", doLongRangeCorrection));

        if (lambdaTerm) {
            logger.info("  Lambda Parameters");
            logger.info(format("   Softcore Alpha:                        %5.3f", vdwLambdaAlpha));
            logger.info(format("   Lambda Exponent:                       %5.3f", vdwLambdaExponent));
        }
    }

    /**
     * Allocate coordinate arrays and set up reduction indices and values.
     */
    private void initAtomArrays() {
        if (atomClass == null || nAtoms > atomClass.length) {
            atomClass = new int[nAtoms];
            coordinates = new double[nAtoms * 3];
            reduced = new double[nSymm][nAtoms * 3];
            reducedXYZ = reduced[0];
            reductionIndex = new int[nAtoms];
            reductionValue = new double[nAtoms];
            bondMask = new int[nAtoms][];
            angleMask = new int[nAtoms][];
            if (vdwForm == VDW_FORM.LENNARD_JONES_6_12) {
                torsionMask = new int[nAtoms][];
            } else {
                torsionMask = null;
            }
            use = new boolean[nAtoms];
            isSoft = new boolean[nAtoms];
            softCore = new boolean[2][nAtoms];
            gradX = new double[threadCount][];
            gradY = new double[threadCount][];
            gradZ = new double[threadCount][];
            if (lambdaTerm) {
                lambdaGradX = new double[threadCount][];
                lambdaGradY = new double[threadCount][];
                lambdaGradZ = new double[threadCount][];
            } else {
                lambdaGradX = null;
                lambdaGradY = null;
                lambdaGradZ = null;
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
        molecule = molecularAssembly.getMoleculeNumbers();

        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            assert (i == ai.xyzIndex - 1);
            double xyz[] = ai.getXYZ();
            int i3 = i * 3;
            coordinates[i3 + XX] = xyz[XX];
            coordinates[i3 + YY] = xyz[YY];
            coordinates[i3 + ZZ] = xyz[ZZ];
            AtomType type = ai.getAtomType();
            if (type == null) {
                logger.severe(ai.toString());
                continue;
            }
            atomClass[i] = type.atomClass;
            ForceField forceField = molecularAssembly.getForceField();
            VDWType vdwType = forceField.getVDWType(Integer.toString(atomClass[i]));
            ai.setVDWType(vdwType);
            ArrayList<Bond> bonds = ai.getBonds();
            int numBonds = bonds.size();
            if (vdwType.reductionFactor > 0.0 && numBonds == 1) {
                Bond bond = bonds.get(0);
                Atom heavyAtom = bond.get1_2(ai);
                // Atom indexes start at 1
                reductionIndex[i] = heavyAtom.xyzIndex - 1;
                reductionValue[i] = vdwType.reductionFactor;
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
            if (vdwForm == VDW_FORM.LENNARD_JONES_6_12) {
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

    public void setAtoms(Atom atoms[]) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        initAtomArrays();

        /**
         * Rebuild the NeighborList.
         */
        neighborList.setAtoms(atoms);
        neighborListOnly = true;
        print = false;
        try {
            parallelTeam.execute(vanDerWaalsRegion);
        } catch (Exception e) {
            String message = " Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
        }
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
        int n = (int) (100.0 * (maxR - cut));
        double delR = (maxR - cut) / n;
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
                double irv = radEps[i][j2 + RADMIN];
                double ev = radEps[i][j2 + EPS];
                double sume = 0.0;
                for (int k = 0; k <= n; k++) {
                    double r = cut + k * delR;
                    double r2 = r * r;
                    final double rho = r * irv;
                    final double rho3 = rho * rho * rho;
                    final double rho7 = rho3 * rho3 * rho;
                    final double rhod = rho + delta;
                    final double rhod3 = rhod * rhod * rhod;
                    final double rhod7 = rhod3 * rhod3 * rhod;
                    final double t1 = t1n / rhod7;
                    final double t2 = gamma1 / rho7 + gamma;
                    final double eij = ev * t1 * (t2 - 2.0);
                    /**
                     * Apply one minus the multiplicative switch if the
                     * interaction distance is less than the end of the
                     * switching window.
                     */
                    double taper = 1.0;
                    if (r2 < off2) {
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
        return this.buff;
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
        if (vdwForm == VDW_FORM.LENNARD_JONES_6_12) {
            final int[] torsionMaski = torsionMask[i];
            final int n14 = torsionMaski.length;
            for (int m = 0; m < n14; m++) {
                mask[torsionMaski[m]] = 0.5;
            }
        }
        final int[] angleMaski = angleMask[i];
        final int n13 = angleMaski.length;
        for (int m = 0; m < n13; m++) {
            mask[angleMaski[m]] = 0;
        }
        final int[] bondMaski = bondMask[i];
        final int n12 = bondMaski.length;
        for (int m = 0; m < n12; m++) {
            mask[bondMaski[m]] = 0;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Remove the masking rules for 1-2 and 1-3 interactions.
     */
    @Override
    public void removeMask(final double mask[], final int i) {
        if (vdwForm == VDW_FORM.LENNARD_JONES_6_12) {
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
        int classi = atoms[i].getAtomType().atomClass;
        int classk = atoms[k].getAtomType().atomClass;
        logger.info(String.format("%s %6d-%s %6d-%s %10.4f  %10.4f  %10.4f",
                "VDW", atoms[i].xyzIndex, atoms[i].getAtomType().name,
                atoms[k].xyzIndex, atoms[k].getAtomType().name,
                1.0 / radEps[classi][classk * 2 + RADMIN], r, eij));
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
         * Initialize the softcore atom masks.
         */
        if (!softCoreInit) {
            for (int i = 0; i < nAtoms; i++) {
                isSoft[i] = atoms[i].applyLambda();
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

        // Redo the long range correction.
        if (doLongRangeCorrection) {
            longRangeCorrection = getLongRangeCorrection();
            logger.info(String.format(" Long-range vdW correction %12.8f (kcal/mole).",
                    longRangeCorrection));
        } else {
            longRangeCorrection = 0.0;
        }
    }

    public void setIntermolecularSoftcore(boolean intermolecularSoftcore) {
        this.intermolecularSoftcore = intermolecularSoftcore;
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
        return shareddEdL.get();
    }

    /**
     * {@inheritDoc}
     *
     * @param lambdaGradient
     */
    @Override
    public void getdEdXdL(double[] lambdaGradient) {
        int index = 0;
        double lgx[] = lambdaGradX[0];
        double lgy[] = lambdaGradY[0];
        double lgz[] = lambdaGradZ[0];
        for (int i = 0; i < nAtoms; i++) {
            lambdaGradient[index++] += lgx[i];
            lambdaGradient[index++] += lgy[i];
            lambdaGradient[index++] += lgz[i];
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
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
        public void start() {
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
        }

        @Override
        public void finish() {
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
                boolean forceRebuild = false;
                if (neighborListOnly) {
                    forceRebuild = true;
                }
                neighborList.buildList(reduced, neighborLists, null, forceRebuild, print);
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

                logger.info(String.format("\n van der Waals: %7.4f (sec)", total));
                logger.info(" Thread    Init    Energy  Reduce  Total     Counts");
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
                    logger.info(String.format("    %3d   %7.4f %7.4f %7.4f %7.4f %10d",
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
                logger.info(String.format(" Min      %7.4f %7.4f %7.4f %7.4f %10d",
                        initMin * 1e-9, vdwMin * 1e-9,
                        reductionMin * 1e-9, totalMin * 1e-9, countMin));
                logger.info(String.format(" Max      %7.4f %7.4f %7.4f %7.4f %10d",
                        initMax * 1e-9, vdwMax * 1e-9,
                        reductionMax * 1e-9, totalMax * 1e-9, countMax));
                logger.info(String.format(" Delta    %7.4f %7.4f %7.4f %7.4f %10d",
                        (initMax - initMin) * 1e-9, (vdwMax - vdwMin) * 1e-9,
                        (reductionMax - reductionMin) * 1e-9, (totalMax - totalMin) * 1e-9,
                        (countMax - countMin)));
                logger.info(String.format(" Actual   %7.4f %7.4f %7.4f %7.4f %10d\n",
                        initializationTotal * 1e-9, vdwTotal * 1e-9,
                        reductionTotal * 1e-9, totalActual * 1e-9, sharedInteractions.get()));
            }
        }

        /**
         * Update the local coordinate array and initialize reduction variables.
         */
        private class InitializationLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                initializationTime[getThreadIndex()] = -System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {

                for (int i = lb, i3 = 3 * lb; i <= ub; i++, i3 += 3) {
                    Atom atom = atoms[i];
                    final double xyz[] = atom.getXYZ();
                    coordinates[i3 + XX] = xyz[XX];
                    coordinates[i3 + YY] = xyz[YY];
                    coordinates[i3 + ZZ] = xyz[ZZ];
                    use[i] = atom.isActive();
                }

                int threadIndex = getThreadIndex();

                if (gradient) {
                    if (gradX[threadIndex] == null || gradX[threadIndex].length < nAtoms) {
                        gradX[threadIndex] = new double[nAtoms];
                        gradY[threadIndex] = new double[nAtoms];
                        gradZ[threadIndex] = new double[nAtoms];
                    } else {
                        fill(gradX[threadIndex], 0.0);
                        fill(gradY[threadIndex], 0.0);
                        fill(gradZ[threadIndex], 0.0);
                    }
                }

                if (lambdaTerm) {
                    if (lambdaGradX[threadIndex] == null || lambdaGradX[threadIndex].length < nAtoms) {
                        lambdaGradX[threadIndex] = new double[nAtoms];
                        lambdaGradY[threadIndex] = new double[nAtoms];
                        lambdaGradZ[threadIndex] = new double[nAtoms];
                    } else {
                        fill(lambdaGradX[threadIndex], 0.0);
                        fill(lambdaGradY[threadIndex], 0.0);
                        fill(lambdaGradZ[threadIndex], 0.0);
                    }
                }
            }
        }

        private class ExpandLoop extends IntegerForLoop {

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
            public void finish() {
                initializationTime[getThreadIndex()] += System.nanoTime();
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
            private double gxi_local[];
            private double gyi_local[];
            private double gzi_local[];
            private double dEdL;
            private double d2EdL2;
            private double lxi_local[];
            private double lyi_local[];
            private double lzi_local[];
            private final double dx_local[];
            private final double transOp[][];
            private final double mask[];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public VanDerWaalsLoop() {
                super();
                mask = new double[nAtoms];
                dx_local = new double[3];
                transOp = new double[3][3];
                fill(mask, 1.0);
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
                int threadId = getThreadIndex();
                gxi_local = gradX[threadId];
                gyi_local = gradY[threadId];
                gzi_local = gradZ[threadId];
                energy = 0.0;
                count = 0;
                if (lambdaTerm) {
                    dEdL = 0.0;
                    d2EdL2 = 0.0;
                    lxi_local = lambdaGradX[threadId];
                    lyi_local = lambdaGradY[threadId];
                    lzi_local = lambdaGradZ[threadId];
                } else {
                    lxi_local = null;
                    lyi_local = null;
                    lzi_local = null;
                }
                vdwTime[threadId] = -System.nanoTime();
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
                vdwTime[getThreadIndex()] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                double e = 0.0;
                double xyzS[] = reduced[0];
                int list[][] = neighborLists[0];
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    int i3 = i * 3;
                    final double xi = reducedXYZ[i3++];
                    final double yi = reducedXYZ[i3++];
                    final double zi = reducedXYZ[i3];
                    final int redi = reductionIndex[i];
                    final double redv = reductionValue[i];
                    final double rediv = 1.0 - redv;
                    final int classi = atomClass[i];
                    final double radEpsi[] = radEps[classi];
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
                        if (!use[k]) {
                            continue;
                        }
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
                        if (r2 <= off2 && mask[k] > 0 && irv > 0) {
                            final double r = sqrt(r2);
                            double alpha = 0.0;
                            double lambda5 = 1.0;
                            boolean soft = softCorei[k] || (intermolecularSoftcore && (moleculei != molecule[k]));
                            if (soft) {
                                alpha = sc1;
                                lambda5 = sc2;
                            }
                            final double ev = mask[k] * radEpsi[a2 + EPS];
                            final double eps_lambda = ev * lambda5;
                            final double rho = r * irv;
                            final double rhoDisp1 = pow(rho, dispersivePower1);
                            final double rhoDisp = rhoDisp1 * rho;
                            final double rhoDelta1 = pow(rho + delta, repDispPower1);
                            final double rhoDelta = rhoDelta1 * (rho + delta);
                            final double alphaRhoDelta = alpha + rhoDelta;
                            final double alphaRhoDispGamma = alpha + rhoDisp + gamma;
                            final double t1d = 1.0 / alphaRhoDelta;
                            final double t2d = 1.0 / alphaRhoDispGamma;
                            final double t1 = t1n * t1d;
                            final double t2a = gamma1 * t2d;
                            final double t2 = t2a - 2.0;
                            final double eij = eps_lambda * t1 * t2;
                            /**
                             * Apply a multiplicative switch if the interaction
                             * distance is greater than the beginning of the
                             * taper.
                             */
                            double taper = 1.0;
                            double dtaper = 0.0;
                            if (r2 > cut2) {
                                final double r3 = r2 * r;
                                final double r4 = r2 * r2;
                                final double r5 = r2 * r3;
                                taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                                dtaper = multiplicativeSwitch.dtaper(r, r2, r3, r4);
                            }
                            e += eij * taper;
                            count++;
                            if (!(gradient || (lambdaTerm && soft))) {
                                continue;
                            }
                            final int redk = reductionIndex[k];
                            final double red = reductionValue[k];
                            final double redkv = 1.0 - red;
                            final double dt1d_dr = repDispPower * rhoDelta1 * irv;
                            final double dt2d_dr = dispersivePower * rhoDisp1 * irv;
                            final double dt1_dr = t1 * dt1d_dr * t1d;
                            final double dt2_dr = t2a * dt2d_dr * t2d;
                            final double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);
                            final double ir = 1.0 / r;
                            final double drdx = dx_local[0] * ir;
                            final double drdy = dx_local[1] * ir;
                            final double drdz = dx_local[2] * ir;
                            if (gradient) {
                                final double dswitch = (eij * dtaper + dedr * taper);
                                final double dedx = dswitch * drdx;
                                final double dedy = dswitch * drdy;
                                final double dedz = dswitch * drdz;
                                gxi += dedx * redv;
                                gyi += dedy * redv;
                                gzi += dedz * redv;
                                gxredi += dedx * rediv;
                                gyredi += dedy * rediv;
                                gzredi += dedz * rediv;
                                gxi_local[k] -= red * dedx;
                                gyi_local[k] -= red * dedy;
                                gzi_local[k] -= red * dedz;
                                gxi_local[redk] -= redkv * dedx;
                                gyi_local[redk] -= redkv * dedy;
                                gzi_local[redk] -= redkv * dedz;
                            }
                            if (lambdaTerm && soft) {
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
                                lxi_local[k] -= red * dedldx;
                                lyi_local[k] -= red * dedldy;
                                lzi_local[k] -= red * dedldz;
                                lxi_local[redk] -= redkv * dedldx;
                                lyi_local[redk] -= redkv * dedldy;
                                lzi_local[redk] -= redkv * dedldz;
                            }
                        }
                    }
                    if (gradient) {
                        gxi_local[i] += gxi;
                        gyi_local[i] += gyi;
                        gzi_local[i] += gzi;
                        gxi_local[redi] += gxredi;
                        gyi_local[redi] += gyredi;
                        gzi_local[redi] += gzredi;
                    }
                    if (lambdaTerm) {
                        lxi_local[i] += lxi;
                        lyi_local[i] += lyi;
                        lzi_local[i] += lzi;
                        lxi_local[redi] += lxredi;
                        lyi_local[redi] += lyredi;
                        lzi_local[redi] += lzredi;
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
                        final double xi = reducedXYZ[i3++];
                        final double yi = reducedXYZ[i3++];
                        final double zi = reducedXYZ[i3];
                        final int redi = reductionIndex[i];
                        final double redv = reductionValue[i];
                        final double rediv = 1.0 - redv;
                        final int classi = atomClass[i];
                        final double radEpsi[] = radEps[classi];
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
                            if (r2 <= off2 && irv > 0) {
                                double selfScale = 1.0;
                                if (i == k) {
                                    selfScale = 0.5;
                                }
                                final double r = sqrt(r2);
                                double alpha = 0.0;
                                double lambda5 = 1.0;
                                boolean soft = (isSoft[i] || softCorei[k]);
                                if (soft) {
                                    alpha = sc1;
                                    lambda5 = sc2;
                                }
                                final double ev = radEpsi[a2 + EPS];
                                final double eps_lambda = ev * lambda5;
                                final double rho = r * irv;
                                final double rhoDisp1 = pow(rho, dispersivePower1);
                                final double rhoDisp = rhoDisp1 * rho;
                                final double rhoDelta1 = pow(rho + delta, repDispPower1);
                                final double rhoDelta = rhoDelta1 * (rho + delta);
                                final double alphaRhoDelta = alpha + rhoDelta;
                                final double alphaRhoDispGamma = alpha + rhoDisp + gamma;
                                final double t1d = 1.0 / alphaRhoDelta;
                                final double t2d = 1.0 / alphaRhoDispGamma;
                                final double t1 = t1n * t1d;
                                final double t2a = gamma1 * t2d;
                                final double t2 = t2a - 2.0;
                                double eij = eps_lambda * t1 * t2;
                                /**
                                 * Apply a multiplicative switch if the
                                 * interaction distance is greater than the
                                 * beginning of the taper.
                                 */
                                double taper = 1.0;
                                double dtaper = 0.0;
                                if (r2 > cut2) {
                                    final double r3 = r2 * r;
                                    final double r4 = r2 * r2;
                                    final double r5 = r2 * r3;
                                    taper = multiplicativeSwitch.taper(r, r2, r3, r4, r5);
                                    dtaper = multiplicativeSwitch.dtaper(r, r2, r3, r4);
                                }
                                e += selfScale * eij * taper;
                                count++;
                                if (!(gradient || (lambdaTerm && soft))) {
                                    continue;
                                }
                                final int redk = reductionIndex[k];
                                final double red = reductionValue[k];
                                final double redkv = 1.0 - red;
                                final double dt1d_dr = repDispPower * rhoDelta1 * irv;
                                final double dt2d_dr = dispersivePower * rhoDisp1 * irv;
                                final double dt1_dr = t1 * dt1d_dr * t1d;
                                final double dt2_dr = t2a * dt2d_dr * t2d;
                                double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);
                                final double ir = 1.0 / r;
                                double drdx = dx_local[0] * ir;
                                double drdy = dx_local[1] * ir;
                                double drdz = dx_local[2] * ir;
                                dedr = (eij * dtaper + dedr * taper);
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
                                    gxi_local[k] -= red * dedxk;
                                    gyi_local[k] -= red * dedyk;
                                    gzi_local[k] -= red * dedzk;
                                    gxi_local[redk] -= redkv * dedxk;
                                    gyi_local[redk] -= redkv * dedyk;
                                    gzi_local[redk] -= redkv * dedzk;
                                }
                                if (lambdaTerm && soft) {
                                    double dt1 = -t1 * t1d * dsc1dL;
                                    double dt2 = -t2a * t2d * dsc1dL;
                                    double f1 = dsc2dL * t1 * t2;
                                    double f2 = sc2 * dt1 * t2;
                                    double f3 = sc2 * t1 * dt2;
                                    double dedl = ev * (f1 + f2 + f3);
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
                                    lxi_local[k] -= red * dedldxk;
                                    lyi_local[k] -= red * dedldyk;
                                    lzi_local[k] -= red * dedldzk;
                                    lxi_local[redk] -= redkv * dedldxk;
                                    lyi_local[redk] -= redkv * dedldyk;
                                    lzi_local[redk] -= redkv * dedldzk;
                                }
                            }
                        }
                        if (gradient) {
                            gxi_local[i] += gxi;
                            gyi_local[i] += gyi;
                            gzi_local[i] += gzi;
                            gxi_local[redi] += gxredi;
                            gyi_local[redi] += gyredi;
                            gzi_local[redi] += gzredi;
                        }
                        if (lambdaTerm) {
                            lxi_local[i] += lxi;
                            lyi_local[i] += lyi;
                            lzi_local[i] += lzi;
                            lxi_local[redi] += lxredi;
                            lyi_local[redi] += lyredi;
                            lzi_local[redi] += lzredi;
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

            @Override
            public void start() {
                reductionTime[getThreadIndex()] = -System.nanoTime();
            }

            @Override
            public void finish() {
                reductionTime[getThreadIndex()] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                if (gradient) {
                    double gx[] = gradX[0];
                    double gy[] = gradY[0];
                    double gz[] = gradZ[0];
                    for (int t = 1; t < threadCount; t++) {
                        double gxt[] = gradX[t];
                        double gyt[] = gradY[t];
                        double gzt[] = gradZ[t];
                        for (int i = lb; i <= ub; i++) {
                            gx[i] += gxt[i];
                            gy[i] += gyt[i];
                            gz[i] += gzt[i];
                        }
                    }

                    for (int i = lb; i <= ub; i++) {
                        Atom ai = atoms[i];
                        ai.addToXYZGradient(gradX[0][i], gradY[0][i], gradZ[0][i]);
                    }
                }
                if (lambdaTerm) {
                    double lx[] = lambdaGradX[0];
                    double ly[] = lambdaGradY[0];
                    double lz[] = lambdaGradZ[0];
                    for (int t = 1; t < threadCount; t++) {
                        double lxt[] = lambdaGradX[t];
                        double lyt[] = lambdaGradY[t];
                        double lzt[] = lambdaGradZ[t];
                        for (int i = lb; i <= ub; i++) {
                            lx[i] += lxt[i];
                            ly[i] += lyt[i];
                            lz[i] += lzt[i];
                        }
                    }
                }
            }
        }
    }
    /**
     * *************************************************************************
     * Cutoff and switching constants.
     */
    /**
     * At the distance "cut", a multiplicative switch begins to be applied.
     */
    private final double cut;
    /**
     * The distance cut squared.
     */
    private final double cut2;
    /**
     * All vdW interactions are 0 at the distance "off".
     */
    private final double off;
    /**
     * The distance off squared.
     */
    private final double off2;
    private final double buff;

    /**
     * *************************************************************************
     * Buffered-14-7 constants.
     */
    /**
     * First constant suggested by Halgren for the Buffered-14-7 potential.
     */
    private final double gamma;
    private final double gamma1;
    /**
     * Second constant suggested by Halgren for the Buffered-14-7 potential.
     */
    private final double delta;
    private final double delta1;
    private final double t1n;
    private final int repulsivePower;
    private final int dispersivePower;
    private final int dispersivePower1;
    private final int repDispPower;
    private final int repDispPower1;

    /**
     * The 6 coefficients of the multiplicative polynomial switch are unique
     * given the distances "off" and "cut". They are found by solving a system
     * of 6 equations, which define the boundary conditions of the switch.
     * <br>
     * f(cut) = 1
     * <br>
     * f'(cut) = f"(cut) = 0
     * <br>
     * f(off) = f'(off) = f"(off) = 0
     */
    private class MultiplicativeSwitch {

        protected final double off;
        protected final double off2;
        protected final double cut;
        protected final double cut2;

        protected final double c0;
        protected final double c1;
        protected final double c2;
        protected final double c3;
        protected final double c4;
        protected final double c5;
        protected final double twoC2;
        protected final double threeC3;
        protected final double fourC4;
        protected final double fiveC5;

        public MultiplicativeSwitch(double off, double cut) {
            this.off = off;
            this.off2 = off * off;
            this.cut = cut;
            this.cut2 = cut * cut;

            double denom = pow(off - cut, 5.0);
            c0 = off * off2 * (off2 - 5.0 * off * cut + 10.0 * cut2) / denom;
            c1 = -30.0 * off2 * cut2 / denom;
            c2 = 30.0 * (off2 * cut + off * cut2) / denom;
            c3 = -10.0 * (off2 + 4.0 * off * cut + cut2) / denom;
            c4 = 15.0 * (off + cut) / denom;
            c5 = -6.0 / denom;
            twoC2 = 2.0 * c2;
            threeC3 = 3.0 * c3;
            fourC4 = 4.0 * c4;
            fiveC5 = 5.0 * c5;
        }

        public double taper(double r, double r2, double r3, double r4, double r5) {
            return c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
        }

        public double dtaper(double r, double r2, double r3, double r4) {
            return fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
        }

    }

}
