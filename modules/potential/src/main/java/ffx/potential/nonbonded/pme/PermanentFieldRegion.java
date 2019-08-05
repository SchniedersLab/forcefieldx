//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.potential.nonbonded.pme;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.util.Range;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Torsion;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.LambdaMode;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.PCGVariables;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.ScaleFactors;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.parameters.ForceField;
import static ffx.numerics.special.Erf.erfc;
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

/**
 * The Permanent Field Region should be executed by a ParallelTeam with
 * exactly 2 threads. The Real Space and Reciprocal Space Sections will be
 * run concurrently, each with the number of threads defined by their
 * respective ParallelTeam instances.
 */
public class PermanentFieldRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PermanentFieldRegion.class.getName());

    private PermanentRealSpaceFieldSection permanentRealSpaceFieldSection;
    private PermanentReciprocalSection permanentReciprocalSection;

    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    /**
     * Unit cell and spacegroup information.
     */
    private Crystal crystal;
    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    private double[][][] coordinates;
    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    private double[][][] globalMultipole;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public double[][][] inducedDipole;
    public double[][][] inducedDipoleCR;
    /**
     * Neighbor lists, including atoms beyond the real space cutoff.
     * [nsymm][nAtoms][nAllNeighbors]
     */
    private int[][][] neighborLists;
    /**
     * Neighbor lists, without atoms beyond the real space cutoff.
     * [nSymm][nAtoms][nIncludedNeighbors]
     */
    private int[][][] realSpaceLists;
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
    /**
     * When computing the polarization energy at Lambda there are 3 pieces.
     * <p>
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand).
     * <p>
     * 2.) Uenv = The polarization energy of the system without the ligand.
     * <p>
     * 3.) Uligand = The polarization energy of the ligand by itself.
     * <p>
     * Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     * <p>
     * Set the "use" array to true for all atoms for part 1. Set the "use" array
     * to true for all atoms except the ligand for part 2. Set the "use" array
     * to true only for the ligand atoms for part 3.
     * <p>
     * The "use" array can also be employed to turn off atoms for computing the
     * electrostatic energy of sub-structures.
     */
    private boolean[] use;
    /**
     * Molecule number for each atom.
     */
    private int[] molecule;
    private double[] ipdamp;
    private double[] thole;
    /**
     * Polarization groups.
     */
    protected int[][] ip11;
    /**
     * The current LambdaMode of this PME instance (or OFF for no lambda dependence).
     */
    private LambdaMode lambdaMode = LambdaMode.OFF;

    private IntegerSchedule permanentSchedule;
    private boolean reciprocalSpaceTerm;
    private double off2;
    private double preconditionerCutoff;
    private double an0, an1, an2;
    private double aewald;
    /**
     * Reciprocal space instance.
     */
    private ReciprocalSpace reciprocalSpace;
    /**
     * Number of neighboring atoms within the real space cutoff. [nSymm][nAtoms]
     */
    private int[][] realSpaceCounts;
    private Range[] realSpaceRanges;
    /**
     * Field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] field;
    /**
     * Chain rule field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] fieldCR;
    /**
     * Timing variables.
     */
    private long realSpacePermTotal;
    private long[] realSpacePermTime;

    private ScaleFactors scaleFactors;

    /**
     * Specify inter-molecular softcore.
     */
    private final boolean intermolecularSoftcore;
    /**
     * Specify intra-molecular softcore.
     */
    private final boolean intramolecularSoftcore;
    /**
     * Constant applied to multipole interactions.
     */
    private static final double oneThird = 1.0 / 3.0;

    public PermanentFieldRegion(ParallelTeam pt, ForceField forceField, boolean lambdaTerm) {
        permanentRealSpaceFieldSection = new PermanentRealSpaceFieldSection(pt);
        permanentReciprocalSection = new PermanentReciprocalSection();

        // Flag to indicate application of an intermolecular softcore potential.
        if (lambdaTerm) {
            intermolecularSoftcore = forceField.getBoolean(
                    ForceField.ForceFieldBoolean.INTERMOLECULAR_SOFTCORE, false);
            intramolecularSoftcore = forceField.getBoolean(
                    ForceField.ForceFieldBoolean.INTRAMOLECULAR_SOFTCORE, false);
        } else {
            intermolecularSoftcore = false;
            intramolecularSoftcore = false;
        }

    }

    public void init(Atom[] atoms, Crystal crystal, double[][][] coordinates, double[][][] globalMultipole,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR, int[][][] neighborLists,
                     int[][][] realSpaceLists, ScaleFactors scaleFactors,
                     boolean[] use, int[] molecule, double[] ipdamp, double[] thole, int[][] ip11,
                     LambdaMode lambdaMode, IntegerSchedule permanentSchedule, boolean reciprocalSpaceTerm,
                     EwaldParameters ewaldParameters, ReciprocalSpace reciprocalSpace, PCGVariables pcgVariables,
                     int[][] realSpaceCounts, Range[] realSpaceRanges, double[][][] field, double[][][] fieldCR,
                     long realSpacePermTotal, long[] realSpacePermTime) {
        this.atoms = atoms;
        this.crystal = crystal;
        this.coordinates = coordinates;
        this.globalMultipole = globalMultipole;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.neighborLists = neighborLists;
        this.realSpaceLists = realSpaceLists;
        this.scaleFactors = scaleFactors;
        this.use = use;
        this.molecule = molecule;
        this.ipdamp = ipdamp;
        this.thole = thole;
        this.ip11 = ip11;
        this.lambdaMode = lambdaMode;
        this.permanentSchedule = permanentSchedule;
        this.reciprocalSpaceTerm = reciprocalSpaceTerm;
        if (pcgVariables != null) {
            this.preconditionerCutoff = pcgVariables.preconditionerCutoff;
            this.preconditionerLists = pcgVariables.preconditionerLists;
            this.preconditionerCounts = pcgVariables.preconditionerCounts;
        }
        this.aewald = ewaldParameters.aewald;
        this.an0 = ewaldParameters.an0;
        this.an1 = ewaldParameters.an1;
        this.an2 = ewaldParameters.an2;
        this.off2 = ewaldParameters.off2;
        this.reciprocalSpace = reciprocalSpace;
        this.realSpaceCounts = realSpaceCounts;
        this.realSpaceRanges = realSpaceRanges;
        this.field = field;
        this.fieldCR = fieldCR;
        this.realSpacePermTotal = realSpacePermTotal;
        this.realSpacePermTime = realSpacePermTime;
    }

    public long getRealSpacePermTotal() {
        return realSpacePermTotal;
    }

    @Override
    public void run() {
        try {
            execute(permanentRealSpaceFieldSection, permanentReciprocalSection);
        } catch (RuntimeException e) {
            String message = "Runtime exception computing the permanent multipole field.\n";
            logger.log(Level.WARNING, message, e);
            throw e;
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

        PermanentRealSpaceFieldSection(ParallelTeam pt) {
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
            } catch (RuntimeException e) {
                String message = "Fatal exception computing the real space field.\n";
                logger.log(Level.WARNING, message, e);
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

        private final InitializationLoop[] initializationLoop;
        private final PermanentRealSpaceFieldLoop[] permanentRealSpaceFieldLoop;
        private final SharedInteger sharedCount;
        private final int threadCount;

        PermanentRealSpaceFieldRegion(int nt) {
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
                int nAtoms = atoms.length;
                execute(0, nAtoms - 1, initializationLoop[threadIndex]);
                execute(0, nAtoms - 1, permanentRealSpaceFieldLoop[threadIndex]);
            } catch (RuntimeException e) {
                String message = "Runtime exception computing the real space field.\n";
                logger.log(Level.SEVERE, message, e);
            } catch (Exception e) {
                String message = "Fatal exception computing the real space field in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        @Override
        public void finish() {
            if (realSpaceRanges == null) {
                logger.severe(" RealSpaceRange array is null");
            }

            int nAtoms = atoms.length;

            // Load balancing.
            int id = 0;
            int goal = sharedCount.get() / threadCount;
            int num = 0;
            int start = 0;

            for (int i = 0; i < nAtoms; i++) {
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                int nSymm = symOps.size();
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    num += realSpaceCounts[iSymm][i];
                }
                if (num >= goal) {
                    // Last thread gets the remaining atoms in its range.
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

                    // Out of atoms. Threads remaining get a null range.
                    if (start == nAtoms) {
                        for (int j = id; j < threadCount; j++) {
                            realSpaceRanges[j] = null;
                        }
                        break;
                    }
                } else if (i == nAtoms - 1) {

                    // Last atom without reaching goal for current thread.
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
                double[] fX = field[threadIndex][0];
                double[] fY = field[threadIndex][1];
                double[] fZ = field[threadIndex][2];
                double[] fXCR = fieldCR[threadIndex][0];
                double[] fYCR = fieldCR[threadIndex][1];
                double[] fZCR = fieldCR[threadIndex][2];
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
                // Initialize the induced dipole arrays.
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                int nSymm = symOps.size();
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    double[][] ind0 = inducedDipole[0];
                    double[][] indCR0 = inducedDipoleCR[0];
                    for (int i = lb; i <= ub; i++) {
                        double[] ind = ind0[i];
                        double[] indCR = indCR0[i];
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

            private final double[] dx_local;
            private final double[][] transOp;
            private double[] fX, fY, fZ;
            private double[] fXCR, fYCR, fZCR;
            private double[] mask_local;
            private double[] maskp_local;
            private int count;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            PermanentRealSpaceFieldLoop() {
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
                int nAtoms = atoms.length;
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
                int[][] lists = neighborLists[0];
                int[][] ewalds = realSpaceLists[0];
                int[] counts = realSpaceCounts[0];
                int[][] preLists = preconditionerLists[0];
                int[] preCounts = preconditionerCounts[0];
                final double[] x = coordinates[0][0];
                final double[] y = coordinates[0][1];
                final double[] z = coordinates[0][2];
                final double[][] mpole = globalMultipole[0];
                // Loop over atom chunk.
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
                    final double[] globalMultipolei = mpole[i];
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

                    // Apply energy masking rules.
                    Atom ai = atoms[i];
                    for (Torsion torsion : ai.getTorsions()) {
                        Atom ak = torsion.get1_4(ai);
                        if (ak != null) {
                            int index = ak.getIndex() - 1;
                            maskp_local[index] = scaleFactors.p14scale;
                            for (int k : ip11[i]) {
                                if (k == index) {
                                    maskp_local[index] = scaleFactors.intra14Scale * scaleFactors.p14scale;
                                }
                            }
                        }
                    }
                    for (Angle angle : ai.getAngles()) {
                        Atom ak = angle.get1_3(ai);
                        if (ak != null) {
                            int index = ak.getIndex() - 1;
                            maskp_local[index] = scaleFactors.p13scale;
                        }
                    }
                    for (Bond bond : ai.getBonds()) {
                        int index = bond.get1_2(ai).getIndex() - 1;
                        maskp_local[index] = scaleFactors.p12scale;
                    }

                    // Apply group based polarization masking rule.
                    for (int index : ip11[i]) {
                        mask_local[index] = scaleFactors.d11scale;
                    }

                    // Loop over the neighbor list.
                    final int[] list = lists[i];
                    counts[i] = 0;
                    preCounts[i] = 0;
                    final int[] ewald = ewalds[i];
                    int[] preList = preLists[i];

                    for (int k : list) {
                        if (!use[k]) {
                            continue;
                        }
                        if (lambdaMode == ParticleMeshEwald.LambdaMode.VAPOR) {
                            boolean sameMolecule = (moleculei == molecule[k]);
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
                            final double[] globalMultipolek = mpole[k];
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

                            // Calculate the error function damping terms.
                            final double ralpha = aewald * r;
                            final double exp2a = exp(-ralpha * ralpha);
                            final double rr1 = 1.0 / r;
                            final double rr2 = rr1 * rr1;
                            final double bn0 = erfc(ralpha) * rr1;
                            final double bn1 = (bn0 + an0 * exp2a) * rr2;
                            final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                            final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;

                            // Compute the error function scaled and unscaled terms.
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

                    for (Atom ak : ai.get1_5s()) {
                        maskp_local[ak.getIndex() - 1] = 1.0;
                    }
                    for (Torsion torsion : ai.getTorsions()) {
                        Atom ak = torsion.get1_4(ai);
                        if (ak != null) {
                            int index = ak.getIndex() - 1;
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
                }

                // Loop over symmetry mates.
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                int nSymm = symOps.size();
                for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                    SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                    crystal.getTransformationOperator(symOp, transOp);
                    lists = neighborLists[iSymm];
                    ewalds = realSpaceLists[iSymm];
                    counts = realSpaceCounts[iSymm];
                    preLists = preconditionerLists[iSymm];
                    preCounts = preconditionerCounts[iSymm];
                    double[] xs = coordinates[iSymm][0];
                    double[] ys = coordinates[iSymm][1];
                    double[] zs = coordinates[iSymm][2];
                    double[][] mpoles = globalMultipole[iSymm];

                    // Loop over atoms in a chunk of the asymmetric unit.
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        final double pdi = ipdamp[i];
                        final double pti = thole[i];
                        final double[] multipolei = mpole[i];
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

                        // Loop over the neighbor list.
                        final int[] list = lists[i];
                        final int npair = list.length;
                        counts[i] = 0;
                        preCounts[i] = 0;
                        final int[] ewald = ewalds[i];
                        final int[] preList = preLists[i];
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
                                final double[] multipolek = mpoles[k];
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

                                // Calculate the error function damping terms.
                                final double ralpha = aewald * r;
                                final double exp2a = exp(-ralpha * ralpha);
                                final double rr1 = 1.0 / r;
                                final double rr2 = rr1 * rr1;
                                final double bn0 = erfc(ralpha) * rr1;
                                final double bn1 = (bn0 + an0 * exp2a) * rr2;
                                final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                                final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;

                                // Compute the error function scaled and unscaled terms.
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
