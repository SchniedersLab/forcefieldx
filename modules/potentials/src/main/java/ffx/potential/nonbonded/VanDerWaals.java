/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.TreeMap;
import java.util.Vector;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
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
 */
public class VanDerWaals extends ParallelRegion implements MaskingInterface {

    private static final Logger logger = Logger.getLogger(VanDerWaals.class.getName());
    /**
     * An array of all atoms in the system.
     */
    private final Atom[] atoms;
    /**
     * A local convenience variable equal to atoms.length.
     */
    private final int nAtoms;
    /**
     * Crystal parameters.
     */
    private final Crystal crystal;
    private final int nSymm;
    private boolean gradient;
    private boolean print;
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
    private final int XX = 0;
    private final int YY = 1;
    private final int ZZ = 2;
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
    /**
     * Each hydrogen vdW site is located a fraction of the way from the heavy
     * atom nucleus to the hydrogen nucleus (~0.9).
     */
    private final double reductionValue[];
    private final double longRangeCorrection;
    private final double radEps[][];
    private static final int RADMIN = 0;
    private static final int RADMIN7 = 1;
    private static final int EPS = 2;
    private final int bondMask[][];
    private final int angleMask[][];
    /***************************************************************************
     * Parallel variables.
     */
    /**
     * The Parallel Team.
     */
    private final ParallelTeam parallelTeam;
    private final SharedInteger interactions;
    private final SharedDouble vdwEnergy;
    private final SharedDoubleArray grad[];
    private final int threadCount;
    private long overheadTime;
    /**
     * The neighbor-list includes 1-2 or 1-3 interactions, but the interactions
     * are masked out. The AMOEBA force field includes 1-4 interactions fully.
     */
    private final NeighborList neighborList;
    private final ExpandRegion expandRegion;
    private final VanDerWaalsLoop vanDerWaalsLoop[];

    /**
     * The VanDerWaals class constructor.
     *
     * @param forceField The ForceField instance contains {@link VDWType} parameters.
     * @param unOrderedAtoms The Atom array can be unordered.
     * @param crystal A valid Crystal is required.
     * @param parallelTeam The parallel environment.
     * @since 1.0
     */
    public VanDerWaals(ForceField forceField, Atom[] unOrderedAtoms,
            Crystal crystal, ParallelTeam parallelTeam) {
        this.atoms = new Atom[unOrderedAtoms.length];
        this.crystal = crystal;
        this.parallelTeam = parallelTeam;
        for (Atom ai : unOrderedAtoms) {
            atoms[ai.xyzIndex - 1] = ai;
        }
        nAtoms = atoms.length;
        nSymm = crystal.spaceGroup.symOps.size();
        // Set up the Buffered-14-7 parameters.
        TreeMap<String, VDWType> vdwTypes = forceField.getVDWTypes();
        int maxClass = 0;
        for (VDWType vdwType : vdwTypes.values()) {
            if (vdwType.atomClass > maxClass) {
                maxClass = vdwType.atomClass;
            }
        }
        radEps = new double[maxClass + 1][3 * (maxClass + 1)];
        // Atom Class numbering starts at 1
        for (VDWType vdwi : vdwTypes.values()) {
            int i = vdwi.atomClass;
            double ri = 0.5 * vdwi.radius;
            double ri2 = ri * ri;
            double ri3 = ri * ri2;
            double e1 = vdwi.wellDepth;
            double se1 = sqrt(e1);
            for (VDWType vdwj : vdwTypes.tailMap(vdwi.key).values()) {
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
                double radmin7 = pow(radmin, 7.0);
                radmin = dhal * radmin;
                /**
                 * HHG
                 */
                double eps = 4.0 * (e1 * e2) / ((se1 + se2) * (se1 + se2));
                radEps[i][j * 3 + RADMIN] = radmin;
                radEps[i][j * 3 + RADMIN7] = radmin7;
                radEps[i][j * 3 + EPS] = eps;
                radEps[j][i * 3 + RADMIN] = radmin;
                radEps[j][i * 3 + RADMIN7] = radmin7;
                radEps[j][i * 3 + EPS] = eps;
            }
        }
        /**
         * Set up the cutoff and polynomial switch.
         */
        double vdwcut = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        
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
         * Long range correction.
         */
        int radCount[] = new int[maxClass + 1];
        for (int i = 0; i < maxClass; i++) {
            radCount[i] = 0;
        }
        for (int i = 0; i < nAtoms; i++) {
            radCount[atomClass[i]]++;
        }
        /**
         * Choose maxR = 60 Angstroms or ~20 sigma.
         * Choose delR to be 0.01 Angstroms.
         */
        double maxR = 60.0;
        double delR = 0.01;

        longRangeCorrection = 0.0;



        // Create the neighbor list.
        neighborLists = new int[nSymm][][];
        // Parallel constructs.
        threadCount = parallelTeam.getThreadCount();
        vanDerWaalsLoop = new VanDerWaalsLoop[threadCount];
        expandRegion = new ExpandRegion();
        for (int i = 0; i < threadCount; i++) {
            vanDerWaalsLoop[i] = new VanDerWaalsLoop();
        }
        grad = new SharedDoubleArray[3];
        grad[0] = new SharedDoubleArray(nAtoms);
        grad[1] = new SharedDoubleArray(nAtoms);
        grad[2] = new SharedDoubleArray(nAtoms);
        interactions = new SharedInteger();
        vdwEnergy = new SharedDouble();
        // Parallel neighbor list builder.
        neighborList = new NeighborList(null, crystal, atoms, off, buff, parallelTeam);

        /**
         * Reduced and expand the coordinates of the asymmetric unit.
         */
        try {
            parallelTeam.execute(expandRegion);
        } catch (Exception e) {
            String message = "Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
        /**
         * Build the neighbor-list using the reduced coordinates.
         */
        neighborList.buildList(reduced, neighborLists, true, true);
    }

    /**
     * Get the total van der Waals potential energy.
     * @return The energy.
     * @since 1.0
     */
    public double getEnergy() {
        return vdwEnergy.get();
    }

    /**
     * Get the number of interacting pairs.
     * @return The interaction count.
     * @since 1.0
     */
    public int getInteractions() {
        return interactions.get();
    }

    /**
     * Get the buffer size.
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
        this.print = print;

        /**
         * Update the local coordinate array.
         */
        for (int i = 0; i < nAtoms; i++) {
            final double xyz[] = atoms[i].getXYZ();
            int i3 = i * 3;
            coordinates[i3++] = xyz[XX];
            coordinates[i3++] = xyz[YY];
            coordinates[i3] = xyz[ZZ];
        }

        /**
         * Reduced and expand the coordinates of the asymmetric unit.
         */
        try {
            parallelTeam.execute(expandRegion);
        } catch (Exception e) {
            String message = "Fatal exception expanding coordinates.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }

        /**
         * Build the neighbor-list (if necessary) using reduced coordinates.
         */
        neighborList.buildList(reduced, neighborLists, false, false);
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = "Fatal exception evaluating van der Waals energy.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
        return vdwEnergy.get();
    }

    /**
     * Apply masking rules for 1-2 and 1-3 interactions.
     *
     * @param mask The masking array.
     * @param i The atom whose mask should be applied.
     */
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
     * Remove the masking rules for 1-2 and 1-3 interactions.
     *
     * @param mask The masking array.
     * @param i The atom whose mask should be removed.
     */
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

    public int[][][] getNeighborLists() {
        return neighborLists;
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 0.l
     */
    @Override
    public void start() {
        overheadTime = System.nanoTime();
        /**
         * Initialize the shared variables.
         */
        vdwEnergy.set(0.0);
        interactions.set(0);
        if (gradient) {
            for (int i = 0; i < nAtoms; i++) {
                grad[0].set(i, 0.0);
                grad[1].set(i, 0.0);
                grad[2].set(i, 0.0);
            }
        }
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 0.l
     */
    @Override
    public void run() {
        try {
            execute(0, nAtoms - 1, vanDerWaalsLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception evaluating van der Waals energy in thread: " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
    }

    /**
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 0.l
     */
    @Override
    public void finish() {
        /**
         * Accumulate the gradients.
         */
        if (gradient) {
            for (int i = 0; i < nAtoms; i++) {
                Atom ai = atoms[i];
                ai.addToGradient(grad[0].get(i), grad[1].get(i), grad[2].get(i));
            }
        }
        long computeTime = 0;
        for (int i = 0; i < threadCount; i++) {
            computeTime += vanDerWaalsLoop[i].getComputeTime();
        }
        overheadTime = System.nanoTime() - overheadTime;
        overheadTime = overheadTime - computeTime / threadCount;
        double compute = (double) computeTime / threadCount * toSeconds;
        double overhead = (double) overheadTime * toSeconds;
        double efficiency = compute / (compute + overhead) * 100;
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
                radEps[classi][classk * 3 + RADMIN] / dhal, r, eij));
    }

    /**
     * Apply hydrogen reductions and then expand the coordinates to P1.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class ExpandRegion extends ParallelRegion {

        private final ExpandLoop expandLoop[];

        public ExpandRegion() {
            super();
            expandLoop = new ExpandLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                expandLoop[i] = new ExpandLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, expandLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception expanding coordinates in thread: " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
                System.exit(-1);
            }
        }

        private class ExpandLoop extends IntegerForLoop {

            private final double in[] = new double[3];
            private final double out[] = new double[3];
            private final IntegerSchedule schedule = IntegerSchedule.fixed();

            @Override
            public IntegerSchedule schedule() {
                return schedule;
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
                Vector<SymOp> symOps = crystal.spaceGroup.symOps;
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
        private final double dx_local[];
        private final double dx2_local[];
        private final double mask[];
        private final IntegerSchedule schedule = IntegerSchedule.fixed();
        // 128 bytes of extra padding to avert cache interference.
        private long p0, p1, p2, p3, p4, p5, p6, p7;
        private long p8, p9, pa, pb, pc, pd, pe, pf;

        public VanDerWaalsLoop() {
            super();
            gxi_local = new double[nAtoms];
            gyi_local = new double[nAtoms];
            gzi_local = new double[nAtoms];
            mask = new double[nAtoms];
            dx_local = new double[3];
            dx2_local = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                mask[i] = 1.0;
            }
        }

        public long getComputeTime() {
            return computeTime;
        }

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            energy = 0.0;
            count = 0;
            if (gradient) {
                for (int i = 0; i < nAtoms; i++) {
                    gxi_local[i] = 0.0;
                    gyi_local[i] = 0.0;
                    gzi_local[i] = 0.0;
                }
            }
            computeTime = 0;
        }

        @Override
        public void run(int lb, int ub) {
            long startTime = System.nanoTime();
            /**
             * Loop over symmetry operators. Interactions between atoms in the
             * asymmetric unit and those in symmetry mates count 1/2 strength.
             */
            Vector<SymOp> symOps = crystal.spaceGroup.symOps;
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
                    if (iSymOp == 0) {
                        applyMask(mask, i);
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
                            final double r = sqrt(r2);
                            final double r3 = r2 * r;
                            final double r4 = r2 * r2;
                            final double r7 = r4 * r3;
                            int a3 = atomClass[k] * 3;
                            final double rv = radEpsi[a3++];
                            final double rv7 = radEpsi[a3++];
                            final double ev = radEpsi[a3];
                            final double rho = r7 + ghal * rv7;
                            final double tau = dhal_plus_one / (r + rv);
                            final double tau3 = tau * tau * tau;
                            final double tau7 = tau3 * tau3 * tau;
                            final double rv7_rho = rv7 / rho;
                            final double evtau7 = ev * tau7;
                            final double eij = evtau7 * rv7 * (ghal_plus_one * rv7_rho - 2.0);
                            /**
                             * Apply a multiplicative switch if the interaction
                             * distance is greater than the beginning of the
                             * taper.
                             */
                            double taper = 1.0;
                            if (r2 > cut2) {
                                final double r5 = r2 * r3;
                                taper = c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
                            }
                            e += eij * taper;
                            count++;
                            if (gradient) {
                                final double r6 = r3 * r3;
                                final double dtau = tau / dhal_plus_one;
                                final double gtau = evtau7 * r6 * ghal_plus_one * rv7_rho * rv7_rho;
                                double de = -7.0 * (dtau * eij + gtau);
                                if (r2 > cut2) {
                                    final double dtaper = fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
                                    de = eij * dtaper + de * taper;
                                }
                                de /= r;
                                final double dedx = de * dx_local[0];
                                final double dedy = de * dx_local[1];
                                final double dedz = de * dx_local[2];
                                gxi += dedx * redv;
                                gyi += dedy * redv;
                                gzi += dedz * redv;
                                gxredi += dedx * rediv;
                                gyredi += dedy * rediv;
                                gzredi += dedz * rediv;
                                final int redk = reductionIndex[k];
                                final double red = reductionValue[k];
                                final double redkv = 1.0 - red;
                                if (iSymOp > 0) {
                                    dx_local[0] = dedx * red;
                                    dx_local[1] = dedy * red;
                                    dx_local[2] = dedz * red;
                                    dx2_local[0] = dedx * redkv;
                                    dx2_local[1] = dedy * redkv;
                                    dx2_local[2] = dedz * redkv;
                                    crystal.applySymRot(dx_local, dx_local,
                                            symOp);
                                    crystal.applySymRot(dx2_local, dx2_local,
                                            symOp);
                                    gxi_local[k] -= 0.5 * dx_local[0];
                                    gyi_local[k] -= 0.5 * dx_local[1];
                                    gzi_local[k] -= 0.5 * dx_local[2];
                                    gxi_local[redk] -= 0.5 * dx2_local[0];
                                    gyi_local[redk] -= 0.5 * dx2_local[1];
                                    gzi_local[redk] -= 0.5 * dx2_local[2];
                                } else {
                                    gxi_local[k] -= dedx * red;
                                    gyi_local[k] -= dedy * red;
                                    gzi_local[k] -= dedz * red;
                                    gxi_local[redk] -= dedx * redkv;
                                    gyi_local[redk] -= dedy * redkv;
                                    gzi_local[redk] -= dedz * redkv;
                                }
                            }
                        }
                    }
                    if (gradient) {
                        if (iSymOp > 0) {
                            gxi_local[i] += 0.5 * gxi;
                            gyi_local[i] += 0.5 * gyi;
                            gzi_local[i] += 0.5 * gzi;
                            gxi_local[redi] += 0.5 * gxredi;
                            gyi_local[redi] += 0.5 * gyredi;
                            gzi_local[redi] += 0.5 * gzredi;
                        } else {
                            gxi_local[i] += gxi;
                            gyi_local[i] += gyi;
                            gzi_local[i] += gzi;
                            gxi_local[redi] += gxredi;
                            gyi_local[redi] += gyredi;
                            gzi_local[redi] += gzredi;
                        }
                    }
                    if (iSymOp == 0) {
                        removeMask(mask, i);
                    }
                }
                if (iSymOp > 0) {
                    energy += e * 0.5;
                } else {
                    energy += e;
                }
            }
            computeTime += System.nanoTime() - startTime;
        }

        @Override
        public void finish() {
            /**
             * Reduce the energy, interaction count and gradients from this
             * thread into the shared variables.
             */
            vdwEnergy.addAndGet(energy);
            interactions.addAndGet(count);
            if (gradient) {
                grad[0].reduce(gxi_local, DoubleOp.SUM);
                grad[1].reduce(gyi_local, DoubleOp.SUM);
                grad[2].reduce(gzi_local, DoubleOp.SUM);
            }
        }
    }
    private final double toSeconds = 0.000000001;
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
     * The constant ghal was suggested by Halgren for the Buffered-14-7
     * potential.
     */
    private static final double ghal = 0.12;
    /**
     * The constant ghal plus 1.0.
     */
    private static final double ghal_plus_one = 1.12;
    /**
     * The constant dhal was suggested by Halgren for the Buffered-14-7
     * potential.
     */
    private static final double dhal = 0.07;
    /**
     * The constant dhal plue 1.0.
     */
    private static final double dhal_plus_one = 1.07;
}
