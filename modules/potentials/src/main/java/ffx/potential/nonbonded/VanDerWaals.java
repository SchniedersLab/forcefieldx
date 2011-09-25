/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.LambdaInterface;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.VDWType;

/**
 * The van der Waals class computes the buffered 14-7 van der Waals interaction
 * used by the AMOEBA force field in parallel using a {@link NeighborList} for any
 * {@link Crystal}.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class VanDerWaals extends ParallelRegion implements MaskingInterface,
        LambdaInterface {

    private static final Logger logger = Logger.getLogger(VanDerWaals.class.getName());
    /**
     * Crystal parameters.
     */
    private final Crystal crystal;
    /**
     * An array of all atoms in the system.
     */
    private final Atom[] atoms;
    /**
     * A local convenience variable equal to atoms.length.
     */
    private final int nAtoms;
    private final int nSymm;
    /***************************************************************************
     * Lambda variables.
     */
    private boolean lambdaTerm;
    private boolean gradient;
    private final boolean isSoft[];
    /**
     * There are 2 softCore arrays of length nAtoms.
     *
     * The first is used for atoms in the outer loop that are hard.
     * This mask equals:
     * false    for inner loop hard atoms
     * true     for inner loop soft atoms
     *
     * The second is used for atoms in the outer loop that are soft.
     * This mask equals:
     * true     for inner loop hard atoms
     * false    for inner loop soft atoms
     */
    private final boolean softCore[][];
    private static final int HARD = 0;
    private static final int SOFT = 1;
    private double lambda = 1.0;
    private double vdwLambdaExponent = 1.0;
    private double vdwLambdaAlpha = 0.05;
    private double sc1 = 0.0;
    private double sc2 = 1.0;
    private double dsc1dL = 0.0;
    private double dsc2dL = 0.0;
    private double d2sc1dL2 = 0.0;
    private double d2sc2dL2 = 0.0;
    /***************************************************************************
     * Coordinate arrays.
     */
    /**
     * A local copy of atomic coordinates, including reductions on the hydrogen
     * atoms.
     */
    private final double coordinates[];
    private final double reduced[][];
    private final double reducedXYZ[];
    private final int[][][] neighborLists;
    private static final int XX = 0;
    private static final int YY = 1;
    private static final int ZZ = 2;
    /***************************************************************************
     * Force field parameters and constants for the Buffered-14-7 potential.
     */
    /**
     * A local reference to the atom class of each atom in the system.
     */
    private final int atomClass[];
    /**
     * Hydrogen atom vdW sites are located toward their heavy atom relative to
     * their nucleus. This is a look-up that gives the heavy atom index for each
     * hydrogen.
     */
    private final int reductionIndex[];
    private final int bondMask[][];
    private final int angleMask[][];
    /**
     * Each hydrogen vdW site is located a fraction of the way from the heavy
     * atom nucleus to the hydrogen nucleus (~0.9).
     */
    private final double reductionValue[];
    private final double radEps[][];
    private double longRangeCorrection;
    private final boolean doLongRangeCorrection;
    private int maxClass;
    private static final int RADMIN = 0;
    private static final int EPS = 1;
    /***************************************************************************
     * Parallel variables.
     */
    /**
     * The Parallel Team.
     */
    private final ParallelTeam parallelTeam;
    private final IntegerSchedule pairwiseSchedule;
    private final int threadCount;
    private final SharedInteger sharedInteractions;
    private final SharedDouble sharedEnergy;
    private final SharedDouble shareddEdL;
    private final SharedDouble sharedd2EdL2;
    private final SharedDoubleArray shareddEdLdX[];
    private final double gradX[][];
    private final double gradY[][];
    private final double gradZ[][];
    private final double lambdaGradX[][];
    private final double lambdaGradY[][];
    private final double lambdaGradZ[][];
    /**
     * The neighbor-list includes 1-2 or 1-3 interactions, but the interactions
     * are masked out. The AMOEBA force field includes 1-4 interactions fully.
     */
    private final InitializationRegion initializationRegion;
    private final NeighborList neighborListBuilder;
    private final VanDerWaalsLoop vanDerWaalsLoop[];
    private final ReductionRegion reductionRegion;

    /**
     * The VanDerWaals class constructor.
     *
     * @param forceField The ForceField instance contains {@link VDWType} parameters.
     * @param atoms An ordered (by xyzIndex) atom array.
     * @param crystal A valid Crystal is required.
     * @param parallelTeam The parallel environment.
     * @since 1.0
     */
    public VanDerWaals(ForceField forceField, Atom[] atoms,
            Crystal crystal, ParallelTeam parallelTeam) {
        this.atoms = atoms;
        this.parallelTeam = parallelTeam;
        this.crystal = crystal;
        nAtoms = atoms.length;
        nSymm = this.crystal.spaceGroup.getNumberOfSymOps();

        /**
         * Set up the Buffered-14-7 parameters.
         */
        Map<String, VDWType> map = forceField.getVDWTypes();
        TreeMap<String, VDWType> vdwTypes = new TreeMap<String, VDWType>(map);
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
                /**
                 * Cubic-mean.
                 */
                double radmin = 2.0 * (ri3 + rj3) / (ri2 + rj2);
                /**
                 * HHG
                 */
                double eps = 4.0 * (e1 * e2) / ((se1 + se2) * (se1 + se2));
                radEps[i][j * 2 + RADMIN] = radmin;
                radEps[i][j * 2 + EPS] = eps;
                radEps[j][i * 2 + RADMIN] = radmin;
                radEps[j][i * 2 + EPS] = eps;
            }
        }

        /**
         * Allocate coordinate arrays and set up reduction indices and values.
         */
        coordinates = new double[nAtoms * 3];
        reduced = new double[nSymm][nAtoms * 3];
        reducedXYZ = reduced[0];
        atomClass = new int[nAtoms];
        reductionIndex = new int[nAtoms];
        reductionValue = new double[nAtoms];
        bondMask = new int[nAtoms][];
        angleMask = new int[nAtoms][];
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            double xyz[] = ai.getXYZ();
            int i3 = i * 3;
            coordinates[i3 + XX] = xyz[XX];
            coordinates[i3 + YY] = xyz[YY];
            coordinates[i3 + ZZ] = xyz[ZZ];
            atomClass[i] = ai.getAtomType().atomClass;
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
                Bond b = bonds.get(j);
                bondMask[i][j] = b.get1_2(ai).xyzIndex - 1;
            }
            ArrayList<Angle> angles = ai.getAngles();
            int numAngles = 0;
            for (Angle a : angles) {
                Atom ak = a.get1_3(ai);
                if (ak != null) {
                    numAngles++;
                }
            }
            angleMask[i] = new int[numAngles];
            int j = 0;
            for (Angle a : angles) {
                Atom ak = a.get1_3(ai);
                if (ak != null) {
                    angleMask[i][j++] = ak.xyzIndex - 1;
                }
            }
        }

        /**
         * Set up the cutoff and polynomial switch.
         */
        double vdwcut = 0.0;
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

        /**
         * Initialize the soft core lambda masks.
         */
        isSoft = new boolean[nAtoms];
        softCore = new boolean[2][nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            isSoft[i] = false;
            softCore[HARD][i] = false;
            softCore[SOFT][i] = false;
        }
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);
        if (lambdaTerm) {
            vdwLambdaAlpha = forceField.getDouble(ForceFieldDouble.VDW_LAMBDA_ALPHA, 0.05);
            vdwLambdaExponent = forceField.getDouble(ForceFieldDouble.VDW_LAMBDA_EXPONENT, 1.0);
            if (vdwLambdaAlpha < 0.0) {
                vdwLambdaAlpha = 0.05;
            }
            if (vdwLambdaExponent < 1.0) {
                vdwLambdaExponent = 1.0;
            }
        }
        /**
         * Parallel constructs.
         */
        threadCount = parallelTeam.getThreadCount();
        if (lambdaTerm) {
            shareddEdL = new SharedDouble();
            sharedd2EdL2 = new SharedDouble();
            shareddEdLdX = new SharedDoubleArray[3];
            shareddEdLdX[0] = new SharedDoubleArray(nAtoms);
            shareddEdLdX[1] = new SharedDoubleArray(nAtoms);
            shareddEdLdX[2] = new SharedDoubleArray(nAtoms);
            lambdaGradX = new double[threadCount][nAtoms];
            lambdaGradY = new double[threadCount][nAtoms];
            lambdaGradZ = new double[threadCount][nAtoms];
        } else {
            shareddEdL = null;
            sharedd2EdL2 = null;
            shareddEdLdX = null;
            lambdaGradX = null;
            lambdaGradY = null;
            lambdaGradZ = null;
        }
        doLongRangeCorrection = forceField.getBoolean(ForceField.ForceFieldBoolean.VDWLRTERM, false);

        sharedInteractions = new SharedInteger();
        sharedEnergy = new SharedDouble();
        gradX = new double[threadCount][nAtoms];
        gradY = new double[threadCount][nAtoms];
        gradZ = new double[threadCount][nAtoms];

        initializationRegion = new InitializationRegion();
        reductionRegion = new ReductionRegion();
        vanDerWaalsLoop = new VanDerWaalsLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            vanDerWaalsLoop[i] = new VanDerWaalsLoop(i);
        }

        /**
         * Parallel neighbor list builder.
         */
        neighborListBuilder = new NeighborList(null, this.crystal, atoms, off, buff, parallelTeam);
        neighborLists = new int[nSymm][][];

        /**
         * Reduce and expand the coordinates of the asymmetric unit.
         */
        try {
            parallelTeam.execute(initializationRegion);
        } catch (Exception e) {
            String message = "Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
        }

        /**
         * Build the neighbor-list using the reduced coordinates.
         */
        neighborListBuilder.buildList(reduced, neighborLists, null, true, true);
        pairwiseSchedule = neighborListBuilder.getPairwiseSchedule();

        logger.info(" Van der Waals");
        logger.info(format(" Switch Start:                            %5.2f (A)", cut));
        logger.info(format(" Cut-Off:                                 %5.2f (A)", off));
        //logger.info(format(" Long-Range Correction:                   %B", doLongRangeCorrection));        

        if (lambdaTerm) {
            logger.info("\n Van der Waals Lambda Parameters");
            logger.info(String.format(" Softcore alpha:  %5.3f", vdwLambdaAlpha));
            logger.info(String.format(" Lambda exponent: %5.3f", vdwLambdaExponent));
        }
    }

    /**
     * <p>Getter for the field <code>pairwiseSchedule</code>.</p>
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
         * Integrate to maxR = 60 Angstroms or ~20 sigma.
         * Integration step size of delR to be 0.01 Angstroms.
         */
        double maxR = 60.0;
        int n = (int) (100.0 * (maxR - cut));
        double delR = (maxR - cut) / n;
        double total = 0.0;
        /*
        logger.info(String.format(" Long range correction integral: Steps %d, Step Size %8.3f, Window %8.3f-%8.3f",
        n, delR, cut, cut + delR * n));
         */
        /**
         * Loop over vdW types.
         */
        for (int i = 0; i < maxClass; i++) {
            for (int j = 0; j < maxClass; j++) {
                int j2 = j * 2;
                double rv = radEps[i][j2 + RADMIN];
                double ev = radEps[i][j2 + EPS];
                double sume = 0.0;
                for (int k = 0; k <= n; k++) {
                    double r = cut + k * delR;
                    double r2 = r * r;
                    final double rho = r / rv;
                    final double rho3 = rho * rho * rho;
                    final double rho7 = rho3 * rho3 * rho;
                    final double rhod = rho + ZERO_07;
                    final double rhod3 = rhod * rhod * rhod;
                    final double rhod7 = rhod3 * rhod3 * rhod;
                    final double t1 = t1n / rhod7;
                    final double t2 = t2n / rho7 + ZERO_12;
                    final double eij = ev * t1 * (t2 - 2.0);
                    /**
                     * Apply one minus the multiplicative switch if the interaction
                     * distance is less than the end of the switching window.
                     */
                    double taper = 1.0;
                    if (r2 < off2) {
                        double r3 = r * r2;
                        double r4 = r2 * r2;
                        double r5 = r2 * r3;
                        taper = c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
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
     *  computed.
     * @param print If true, there is verbose printing.
     * @return The energy.
     * @since 1.0
     */
    public double energy(boolean gradient, boolean print) {
        this.gradient = gradient;

        /**
         * Reduced and expand the coordinates of the asymmetric unit.
         */
        try {
            parallelTeam.execute(initializationRegion);
        } catch (Exception e) {
            String message = " Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
        }

        /**
         * Build the neighbor-list (if necessary) using reduced coordinates.
         */
        neighborListBuilder.buildList(reduced, neighborLists, null, false, false);

        /**
         * Calculate the van der Waals energy.
         */
        try {
            parallelTeam.execute(this);
            /**
             * Reduce the gradient array.
             */
            if (gradient) {
                parallelTeam.execute(reductionRegion);
            }

        } catch (Exception e) {
            String message = " Fatal exception evaluating van der Waals energy.\n";
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
        final int[] bondMaski = bondMask[i];
        final int n12 = bondMaski.length;
        for (int m = 0; m < n12; m++) {
            mask[bondMaski[m]] = 0.0;
        }
        final int[] angleMaski = angleMask[i];
        final int n13 = angleMaski.length;
        for (int m = 0; m < n13; m++) {
            mask[angleMaski[m]] = 0.0;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Remove the masking rules for 1-2 and 1-3 interactions.
     */
    @Override
    public void removeMask(final double mask[], final int i) {
        final int[] bondMaski = bondMask[i];
        final int n12 = bondMaski.length;
        for (int m = 0; m < n12; m++) {
            mask[bondMaski[m]] = 1.0;
        }
        final int[] angleMaski = angleMask[i];
        final int n13 = angleMaski.length;
        for (int m = 0; m < n13; m++) {
            mask[angleMaski[m]] = 1.0;
        }
    }

    /**
     * <p>Getter for the field <code>neighborLists</code>.</p>
     *
     * @return an array of int.
     */
    public int[][][] getNeighborLists() {
        return neighborLists;
    }

    /**
     * {@inheritDoc}
     *
     * This is method should not be called; it is invoked by Parallel Java.
     * @since 0.l
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

    /**
     * {@inheritDoc}
     *
     * This is method should not be called; it is invoked by Parallel Java.
     * @since 0.l
     */
    @Override
    public void run() {
        try {
            execute(0, nAtoms - 1, vanDerWaalsLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception evaluating van der Waals energy in thread: " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
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
                radEps[classi][classk * 2 + RADMIN] / ZERO_07, r, eij));
    }

    /** {@inheritDoc} */
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
        for (int i = 0; i < nAtoms; i++) {
            isSoft[i] = atoms[i].applyLambda();
            if (isSoft[i]) {
                // Outer loop atom hard, inner loop atom soft.
                softCore[HARD][i] = true;
                // Both soft - full interaction.
                softCore[SOFT][i] = false;
            } else {
                // Both hard - full interaction.
                softCore[HARD][i] = false;
                // Outer loop atom soft, inner loop atom hard.
                softCore[SOFT][i] = true;
            }
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

    /** {@inheritDoc} */
    @Override
    public double getLambda() {
        return lambda;
    }

    /** {@inheritDoc} */
    @Override
    public double getdEdL() {
        return shareddEdL.get();
    }

    /** {@inheritDoc} */
    @Override
    public void getdEdXdL(double[] gradient) {
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            gradient[index++] += shareddEdLdX[0].get(i);
            gradient[index++] += shareddEdLdX[1].get(i);
            gradient[index++] += shareddEdLdX[2].get(i);
        }
    }

    /** {@inheritDoc} */
    @Override
    public double getd2EdL2() {
        return sharedd2EdL2.get();
    }

    /**
     * 1.) Initialize the local coordinate array.
     * 2.) Initialize reduction variables.
     * 3.) Apply hydrogen reductions and then expand the coordinates to P1.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class InitializationRegion extends ParallelRegion {

        private final InitializationLoop initializationLoop[];
        private final ExpandLoop expandLoop[];

        public InitializationRegion() {
            initializationLoop = new InitializationLoop[threadCount];
            expandLoop = new ExpandLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                initializationLoop[i] = new InitializationLoop();
                expandLoop[i] = new ExpandLoop();
            }
        }

        @Override
        public void run() {
            try {
                int ti = getThreadIndex();
                execute(0, nAtoms - 1, initializationLoop[ti]);
                execute(0, nAtoms - 1, expandLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception expanding coordinates in thread: " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
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
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    final double xyz[] = atoms[i].getXYZ();
                    int i3 = i * 3;
                    coordinates[i3++] = xyz[XX];
                    coordinates[i3++] = xyz[YY];
                    coordinates[i3] = xyz[ZZ];
                }

                if (gradient) {
                    for (int j = 0; j < threadCount; j++) {
                        for (int i = lb; i <= ub; i++) {
                            gradX[j][i] = 0.0;
                            gradY[j][i] = 0.0;
                            gradZ[j][i] = 0.0;
                        }
                    }
                }
                if (lambdaTerm) {
                    for (int j = 0; j < threadCount; j++) {
                        for (int i = lb; i <= ub; i++) {
                            lambdaGradX[j][i] = 0.0;
                            lambdaGradY[j][i] = 0.0;
                            lambdaGradZ[j][i] = 0.0;
                        }
                    }
                    for (int i = lb; i <= ub; i++) {
                        shareddEdLdX[0].set(i, 0.0);
                        shareddEdLdX[1].set(i, 0.0);
                        shareddEdLdX[2].set(i, 0.0);
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
                    /*
                    if (i == 0) {
                    in[0] = reducedXYZ[iX];
                    in[1] = reducedXYZ[iY];
                    in[2] = reducedXYZ[iZ];
                    crystal.toFractionalCoordinates(in, out);
                    atoms[i].print();
                    logger.info(format("%5d %d CART %10.8f %10.8f %10.8f", i, 0, in[0], in[1], in[2]));
                    logger.info(format("%5d %d FRAC %10.8f %10.8f %10.8f", i, 0, out[0], out[1], out[2]));
                    } */
                }

                List<SymOp> symOps = crystal.spaceGroup.symOps;
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
                            logger.severe(" Atom %d is at a special position: " + atoms[i].toString());
                        }
                    }
                }
            }
        }
    }

    /**
     * The van der Waals loop class contains methods and thread local variables
     * used to evaluate the van der Waals energy and gradients with respect to
     * atomic coordinates.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class VanDerWaalsLoop extends IntegerForLoop {

        private int count;
        private double energy;
        private long computeTime;
        private final double gxi_local[];
        private final double gyi_local[];
        private final double gzi_local[];
        private double dEdL;
        private double d2EdL2;
        private final double lxi_local[];
        private final double lyi_local[];
        private final double lzi_local[];
        private final double dx_local[];
        private final double rotmat[][];
        private final double mask[];
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        public VanDerWaalsLoop(int threadId) {
            super();
            gxi_local = gradX[threadId];
            gyi_local = gradY[threadId];
            gzi_local = gradZ[threadId];
            if (lambdaTerm) {
                lxi_local = lambdaGradX[threadId];
                lyi_local = lambdaGradY[threadId];
                lzi_local = lambdaGradZ[threadId];
            } else {
                lxi_local = null;
                lyi_local = null;
                lzi_local = null;
            }
            mask = new double[nAtoms];
            dx_local = new double[3];
            rotmat = new double[3][3];
            for (int i = 0; i < nAtoms; i++) {
                mask[i] = 1.0;
            }
        }

        public long getComputeTime() {
            return computeTime;
        }

        @Override
        public IntegerSchedule schedule() {
            return pairwiseSchedule;
        }

        @Override
        public void start() {
            energy = 0.0;
            count = 0;
            if (lambdaTerm) {
                dEdL = 0.0;
                d2EdL2 = 0.0;
            }
            computeTime = 0;
        }

        @Override
        public void run(int lb, int ub) {
            long startTime = System.nanoTime();
            /**
             * Loop over symmetry operators.
             */
            List<SymOp> symOps = crystal.spaceGroup.symOps;
            for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                double e = 0.0;
                SymOp symOp = symOps.get(iSymOp);
                double xyzS[] = reduced[iSymOp];
                int list[][] = neighborLists[iSymOp];
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
                    if (iSymOp == 0) {
                        applyMask(mask, i);
                    }
                    // Default is that the outer loop atom is hard.
                    boolean softCorei[] = softCore[HARD];
                    boolean softSymmetry = false;
                    if (isSoft[i]) {
                        softCorei = softCore[SOFT];
                        /**
                         * All interactions between a soft atom and a symmetry
                         * mate atom are turned off.
                         */
                        if (iSymOp > 0) {
                            softSymmetry = true;
                        }
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
                        if (r2 <= off2 && mask[k] > 0.0) {
                            // This will only happen for iSymm > 0.
                            double selfScale = 1.0;
                            if (i == k) {
                                selfScale = 0.5;
                            }
                            final double r = sqrt(r2);
                            final double r3 = r2 * r;
                            final double r4 = r2 * r2;
                            int a2 = atomClass[k] * 2;
                            double alpha = 0.0;
                            double lambda5 = 1.0;
                            boolean soft = (softCorei[k] || softSymmetry);
                            if (soft) {
                                alpha = sc1;
                                lambda5 = sc2;
                            }
                            final double rv = radEpsi[a2 + RADMIN];
                            final double ev = radEpsi[a2 + EPS];
                            final double eps_lambda = ev * lambda5;
                            final double rho = r / rv;
                            final double rho3 = rho * rho * rho;
                            final double rho6 = rho3 * rho3;
                            final double rho7 = rho6 * rho;
                            final double rho_07 = rho + ZERO_07;
                            final double rho_07_pow3 = rho_07 * rho_07 * rho_07;
                            final double rho_07_pow7 = rho_07_pow3 * rho_07_pow3 * rho_07;
                            final double a_rho_07_pow7 = alpha + rho_07_pow7;
                            final double a_rho7_ZERO_12 = alpha + rho7 + ZERO_12;
                            final double t1d = 1.0 / a_rho_07_pow7;
                            final double t2d = 1.0 / a_rho7_ZERO_12;
                            final double t1 = t1n * t1d;
                            final double t2a = t2n * t2d;
                            final double t2 = t2a - 2.0;
                            double eij = eps_lambda * t1 * t2;
                            /**
                             * Apply a multiplicative switch if the interaction
                             * distance is greater than the beginning of the
                             * taper.
                             */
                            double taper = 1.0;
                            double dtaper = 0.0;
                            if (r2 > cut2) {
                                final double r5 = r2 * r3;
                                taper = c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
                                dtaper = fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
                            }
                            e += selfScale * eij * taper;
                            count++;
                            if (!(gradient || (lambdaTerm && soft))) {
                                continue;
                            }
                            final int redk = reductionIndex[k];
                            final double red = reductionValue[k];
                            final double redkv = 1.0 - red;
                            final double rho_07_pow6 = rho_07_pow3 * rho_07_pow3;
                            final double drho_dr = 1.0 / rv;
                            final double dt1d_dr = 7.0 * rho_07_pow6 * drho_dr;
                            final double dt2d_dr = 7.0 * rho6 * drho_dr;
                            final double dt1_dr = t1 * dt1d_dr * t1d;
                            final double dt2_dr = t2a * dt2d_dr * t2d;
                            double dedr = -eps_lambda * (dt1_dr * t2 + t1 * dt2_dr);
                            double ir = 1.0 / r;
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
                                dx_local[0] = dedx;
                                dx_local[1] = dedy;
                                dx_local[2] = dedz;
                                crystal.applyTransSymRot(dx_local, dx_local, symOp, rotmat);
                                dedx = dx_local[0];
                                dedy = dx_local[1];
                                dedz = dx_local[2];
                                gxi_local[k] -= red * dedx;
                                gyi_local[k] -= red * dedy;
                                gzi_local[k] -= red * dedz;
                                gxi_local[redk] -= redkv * dedx;
                                gyi_local[redk] -= redkv * dedy;
                                gzi_local[redk] -= redkv * dedz;
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
                                double d2t1 = -dt1 * t1d * dsc1dL
                                        - t1 * t1d * d2sc1dL2
                                        - t1 * t1d2 * dsc1dL;

                                double d2t2 = -dt2 * t2d * dsc1dL
                                        - t2a * t2d * d2sc1dL2
                                        - t2a * t2d2 * dsc1dL;

                                double df1 = d2sc2dL2 * t1 * t2
                                        + dsc2dL * dt1 * t2
                                        + dsc2dL * t1 * dt2;
                                double df2 = dsc2dL * dt1 * t2
                                        + sc2 * d2t1 * t2
                                        + sc2 * dt1 * dt2;
                                double df3 = dsc2dL * t1 * dt2
                                        + sc2 * dt1 * dt2
                                        + sc2 * t1 * d2t2;
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
                                dx_local[0] = dedldx;
                                dx_local[1] = dedldy;
                                dx_local[2] = dedldz;
                                crystal.applyTransSymRot(dx_local, dx_local, symOp, rotmat);
                                dedldx = dx_local[0];
                                dedldy = dx_local[1];
                                dedldz = dx_local[2];
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
                    if (iSymOp == 0) {
                        removeMask(mask, i);
                    }
                }
                energy += e;
            }
            computeTime += System.nanoTime() - startTime;
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

            /* 
            logger.info(String.format(" Thread %d computed %10d interactions in %8.3f sec.", 
            getThreadIndex(), count, computeTime * toSeconds)); */

        }
    }

    /**
     * Reduce the van der Waals gradient.
     */
    private class ReductionRegion extends ParallelRegion {

        private final ReductionLoop reductionLoop[];

        public ReductionRegion() {
            reductionLoop = new ReductionLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                reductionLoop[i] = new ReductionLoop();
            }
        }

        @Override
        public void run() {
            try {
                int ti = getThreadIndex();
                execute(0, nAtoms - 1, reductionLoop[ti]);
            } catch (Exception e) {
                String message = "Fatal exception reducing van der Waals gradient in thread: " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Reduce van der Waals gradient.
         */
        private class ReductionLoop extends IntegerForLoop {

            @Override
            public void run(int lb, int ub) {
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

                if (lambdaTerm) {
                    for (int i = lb; i <= ub; i++) {
                        double lx = 0.0;
                        double ly = 0.0;
                        double lz = 0.0;
                        for (int t = 0; t < threadCount; t++) {
                            lx += lambdaGradX[t][i];
                            ly += lambdaGradY[t][i];
                            lz += lambdaGradZ[t][i];
                        }
                        shareddEdLdX[0].set(i, lx);
                        shareddEdLdX[1].set(i, ly);
                        shareddEdLdX[2].set(i, lz);
                    }
                }
            }
        }
    }
    /***************************************************************************
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
     * The 6 coefficients of the multiplicative polynomial switch are unique
     * given the distances "off" and "cut". They are found by solving a system
     * of 6 equations, which define the boundary conditions of the switch.
     * f(cut) = 1 f'(cut) = f"(cut) = 0 f(off) = f'(off) = f"(off) = 0
     */
    private final double c0;
    private final double c1;
    private final double c2;
    private final double c3;
    private final double c4;
    private final double c5;
    private final double twoC2;
    private final double threeC3;
    private final double fourC4;
    private final double fiveC5;
    /***************************************************************************
     * Buffered-14-7 constants.
     */
    /**
     * First constant suggested by Halgren for the Buffered-14-7
     * potential.
     */
    private static final double ZERO_12 = 0.12;
    private static final double t2n = 1.12;
    /**
     * Second constant suggested by Halgren for the Buffered-14-7
     * potential.
     */
    private static final double ZERO_07 = 0.07;
    private static final double ONE_07 = 1.07;
    private static final double t1n = pow(ONE_07, 7.0);
}
