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
import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Torsion;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.ELEC_FORM;
import ffx.potential.nonbonded.ParticleMeshEwald.LambdaMode;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.nonbonded.ParticleMeshEwaldCart;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.utils.EnergyException;
import static ffx.numerics.math.VectorMath.cross;
import static ffx.numerics.math.VectorMath.diff;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.r;
import static ffx.numerics.math.VectorMath.scalar;
import static ffx.numerics.math.VectorMath.sum;
import static ffx.numerics.special.Erf.erfc;
import static ffx.potential.nonbonded.ParticleMeshEwald.DEFAULT_ELECTRIC;
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
 * The Real Space Energy Region class parallelizes evaluation of the real
 * space energy and gradient.
 */
public class RealSpaceEnergyRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(RealSpaceEnergyRegion.class.getName());

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
     * Multipole frame definition.
     */
    private MultipoleFrameDefinition[] frame;
    /**
     * Multipole frame defining atoms.
     */
    private int[][] axisAtom;
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
    /**
     * Polarization groups.
     */
    protected int[][] ip11;
    /**
     * Flag for ligand atoms.
     */
    private boolean[] isSoft;
    private double[] ipdamp;
    private double[] thole;

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
     * Pairwise schedule for load balancing.
     */
    private IntegerSchedule realSpaceSchedule;
    private long[] realSpaceEnergyTime;

    /**
     * Gradient array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] grad;
    /**
     * Torque array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] torque;
    /**
     * Partial derivative of the gradient with respect to Lambda.
     * [threadID][X/Y/Z][atomID]
     */
    private double[][][] lambdaGrad;
    /**
     * Partial derivative of the torque with respect to Lambda.
     * [threadID][X/Y/Z][atomID]
     */
    private double[][][] lambdaTorque;
    /**
     * Partial derivative with respect to Lambda.
     */
    private SharedDouble shareddEdLambda;
    /**
     * Second partial derivative with respect to Lambda.
     */
    private SharedDouble sharedd2EdLambda2;
    /**
     * If true, compute coordinate gradient.
     */
    private boolean gradient;
    /**
     * If lambdaTerm is true, some ligand atom interactions with the environment
     * are being turned on/off.
     */
    private boolean lambdaTerm;

    /**
     * The current LambdaMode of this PME instance (or OFF for no lambda dependence).
     */
    private LambdaMode lambdaMode = ParticleMeshEwald.LambdaMode.OFF;
    private Polarization polarization;
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
     * When computing the polarization energy at L there are 3 pieces.
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
     * Set polarizationScale to L for part 1. Set polarizationScale to (1-L) for
     * parts 2 & 3.
     */
    private double polarizationScale = 1.0;

    // *************************************************************************
    // Mutable Particle Mesh Ewald constants.
    private double aewald;
    private double an0;
    private double an1;
    private double an2;
    private double an3;
    private double an4;
    private double an5;

    private final int maxThreads;
    private final double electric;
    /**
     * Constant applied to multipole interactions.
     */
    private static final double oneThird = 1.0 / 3.0;
    private double permanentEnergy;
    private double polarizationEnergy;
    private final SharedInteger sharedInteractions;
    public final RealSpaceEnergyLoop[] realSpaceEnergyLoop;

    private ParticleMeshEwaldCart.ScaleFactors scaleFactors;

    /**
     * Specify inter-molecular softcore.
     */
    private final boolean intermolecularSoftcore;
    /**
     * Specify intra-molecular softcore.
     */
    private final boolean intramolecularSoftcore;

    public RealSpaceEnergyRegion(int nt, ForceField forceField, ELEC_FORM elecForm, boolean lambdaTerm) {
        maxThreads = nt;
        electric = forceField.getDouble(ForceField.ForceFieldDouble.ELECTRIC, DEFAULT_ELECTRIC);
        sharedInteractions = new SharedInteger();
        realSpaceEnergyLoop = new RealSpaceEnergyLoop[nt];

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

    public void init(Atom[] atoms, Crystal crystal, double[][][] coordinates, MultipoleFrameDefinition[] frame,
                     int[][] axisAtom, double[][][] globalMultipole, double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     boolean[] use, int[] molecule, int[][] ip11, boolean[] isSoft, double[] ipdamp, double[] thole,
                     int[][][] realSpaceLists, int[][] realSpaceCounts, IntegerSchedule realSpaceSchedule, long[] realSpaceEnergyTime,
                     double[][][] grad, double[][][] torque, double[][][] lambdaGrad, double[][][] lambdaTorque,
                     SharedDouble shareddEdLambda, SharedDouble sharedd2EdLambda2, boolean gradient, boolean lambdaTerm,
                     LambdaMode lambdaMode, Polarization polarization,
                     ParticleMeshEwaldCart.EwaldParameters ewaldParameters,
                     ParticleMeshEwaldCart.ScaleFactors scaleFactors,
                     ParticleMeshEwaldCart.AlchemicalFactors alchemicalFactors) {
        this.atoms = atoms;
        this.crystal = crystal;
        this.coordinates = coordinates;
        this.frame = frame;
        this.axisAtom = axisAtom;
        this.globalMultipole = globalMultipole;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.use = use;
        this.molecule = molecule;
        this.ip11 = ip11;
        this.isSoft = isSoft;
        this.ipdamp = ipdamp;
        this.thole = thole;
        this.realSpaceLists = realSpaceLists;
        this.realSpaceCounts = realSpaceCounts;
        this.realSpaceSchedule = realSpaceSchedule;
        this.realSpaceEnergyTime = realSpaceEnergyTime;
        this.grad = grad;
        this.torque = torque;
        this.lambdaGrad = lambdaGrad;
        this.lambdaTorque = lambdaTorque;
        this.shareddEdLambda = shareddEdLambda;
        this.sharedd2EdLambda2 = sharedd2EdLambda2;
        this.gradient = gradient;
        this.lambdaTerm = lambdaTerm;
        this.lambdaMode = lambdaMode;
        this.polarization = polarization;
        this.lAlpha = alchemicalFactors.lAlpha;
        this.dlAlpha = alchemicalFactors.dlAlpha;
        this.d2lAlpha = alchemicalFactors.d2lAlpha;
        this.dEdLSign = alchemicalFactors.dEdLSign;
        this.dlPowPerm = alchemicalFactors.dlPowPerm;
        this.d2lPowPerm = alchemicalFactors.d2lPowPerm;
        this.doPermanentRealSpace = alchemicalFactors.doPermanentRealSpace;
        this.permanentScale = alchemicalFactors.permanentScale;
        this.lPowPol = alchemicalFactors.lPowPol;
        this.dlPowPol = alchemicalFactors.dlPowPol;
        this.d2lPowPol = alchemicalFactors.d2lPowPol;
        this.doPolarization = alchemicalFactors.doPolarization;
        this.polarizationScale = alchemicalFactors.polarizationScale;
        this.aewald = ewaldParameters.aewald;
        this.an0 = ewaldParameters.an0;
        this.an1 = ewaldParameters.an1;
        this.an2 = ewaldParameters.an2;
        this.an3 = ewaldParameters.an3;
        this.an4 = ewaldParameters.an4;
        this.an5 = ewaldParameters.an5;
        this.scaleFactors = scaleFactors;
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
            int nAtoms = atoms.length;
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
                throw new EnergyException(format(" The permanent multipole energy of thread %d is %16.8f", i, e), false);
            }
            permanentEnergy += e;
            double ei = realSpaceEnergyLoop[i].inducedEnergy;
            if (Double.isNaN(ei)) {
                throw new EnergyException(format(" The polarization energyof thread %d is %16.8f", i, ei), false);
            }
            polarizationEnergy += ei;
        }
        permanentEnergy *= electric;
        polarizationEnergy *= electric;
    }

    public int getCount(int i) {
        return realSpaceEnergyLoop[i].getCount();
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
        private double[] gX, gY, gZ, tX, tY, tZ;
        private double[] lgX, lgY, lgZ, ltX, ltY, ltZ;
        private double[] gxk_local, gyk_local, gzk_local;
        private double[] txk_local, tyk_local, tzk_local;
        private double[] lxk_local, lyk_local, lzk_local;
        private double[] ltxk_local, ltyk_local, ltzk_local;
        private double[] masking_local;
        private double[] maskingp_local;
        private double[] maskingd_local;
        private final double[] dx_local;
        private final double[][] rot_local;
        private final double[][] work;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        RealSpaceEnergyLoop() {
            super();
            dx_local = new double[3];
            work = new double[15][3];
            rot_local = new double[3][3];
        }

        private void init() {
            int nAtoms = atoms.length;
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
            int nSymm = symOps.size();
            int nAtoms = atoms.length;
            for (iSymm = 0; iSymm < nSymm; iSymm++) {
                SymOp symOp = symOps.get(iSymm);
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
                shareddEdLambda.addAndGet(dUdL * electric);
                sharedd2EdLambda2.addAndGet(d2UdL2 * electric);
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
            final double[] x = coordinates[0][0];
            final double[] y = coordinates[0][1];
            final double[] z = coordinates[0][2];
            final double[][] mpole = globalMultipole[0];
            final double[][] ind = inducedDipole[0];
            final double[][] indp = inducedDipoleCR[0];
            final int[][] lists = realSpaceLists[iSymm];
            final double[] neighborX = coordinates[iSymm][0];
            final double[] neighborY = coordinates[iSymm][1];
            final double[] neighborZ = coordinates[iSymm][2];
            final double[][] neighborMultipole = globalMultipole[iSymm];
            final double[][] neighborInducedDipole = inducedDipole[iSymm];
            final double[][] neighborInducedDipolep = inducedDipoleCR[iSymm];
            for (i = lb; i <= ub; i++) {
                if (!use[i]) {
                    continue;
                }
                final Atom ai = atoms[i];
                final int moleculei = molecule[i];
                if (iSymm == 0) {
                    for (Atom ak : ai.get1_5s()) {
                        masking_local[ak.getIndex() - 1] = scaleFactors.m15scale;
                    }
                    for (Torsion torsion : ai.getTorsions()) {
                        Atom ak = torsion.get1_4(ai);
                        if (ak != null) {
                            int index = ak.getIndex() - 1;
                            masking_local[index] = scaleFactors.m14scale;
                            maskingp_local[index] = scaleFactors.p14scale;
                            for (int j : ip11[i]) {
                                if (j == index) {
                                    maskingp_local[index] = scaleFactors.intra14Scale * scaleFactors.p14scale;
                                }
                            }
                        }
                    }
                    for (Angle angle : ai.getAngles()) {
                        Atom ak = angle.get1_3(ai);
                        if (ak != null) {
                            int index = ak.getIndex() - 1;
                            maskingp_local[index] = scaleFactors.p13scale;
                            masking_local[index] = scaleFactors.m13scale;
                        }
                    }
                    for (Bond bond : ai.getBonds()) {
                        int index = bond.get1_2(ai).getIndex() - 1;
                        maskingp_local[index] = scaleFactors.p12scale;
                        masking_local[index] = scaleFactors.m12scale;
                    }
                    for (int j : ip11[i]) {
                        maskingd_local[j] = scaleFactors.d11scale;
                    }
                }
                final double xi = x[i];
                final double yi = y[i];
                final double zi = z[i];
                final double[] globalMultipolei = mpole[i];
                final double[] inducedDipolei = ind[i];
                final double[] inducedDipolepi = indp[i];
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
                final int[] list = lists[i];
                final int npair = realSpaceCounts[iSymm][i];
                for (int j = 0; j < npair; j++) {
                    k = list[j];
                    if (!use[k]) {
                        continue;
                    }
                    boolean sameMolecule = (moleculei == molecule[k]);
                    if (lambdaMode == ParticleMeshEwald.LambdaMode.VAPOR) {
                        if ((intermolecularSoftcore && !sameMolecule)
                                || (intramolecularSoftcore && sameMolecule)) {
                            continue;
                        }
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
                    final double[] globalMultipolek = neighborMultipole[k];
                    final double[] inducedDipolek = neighborInducedDipole[k];
                    final double[] inducedDipolepk = neighborInducedDipolep[k];
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
                    if (doPermanentRealSpace) {
                        double ei = permanentPair();
                        //log(i,k,r,ei);
                        if (Double.isNaN(ei) || Double.isInfinite(ei)) {
                            String message = format(" %s\n %s\n %s\n "
                                            + "The permanent multipole energy between "
                                            + "atoms %d and %d (%d) is %16.8f at "
                                            + "%16.8f A.", crystal.getUnitCell().toString(),
                                    atoms[i].toString(), atoms[k].toString(), i, k, iSymm, ei, r);
                            throw new EnergyException(message, false);
                        }
                        permanentEnergy += ei;
                        count++;
                    }
                    if (polarization != ParticleMeshEwald.Polarization.NONE && doPolarization) {
                        // Polarization does not use the softcore tensors.
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
                        }
                        double ei = polarizationPair();
                        if (Double.isNaN(ei) || Double.isInfinite(ei)) {
                            String message = format(" %s\n"
                                            + " %s\n with induced dipole: %8.3f %8.3f %8.3f\n"
                                            + " %s\n with induced dipole: %8.3f %8.3f %8.3f\n"
                                            + " The polarization energy due to atoms "
                                            + "%d and %d (%d) is %10.6f at %10.6f A.",
                                    crystal.getUnitCell(), atoms[i], uix, uiy, uiz,
                                    atoms[k], ukx, uky, ukz, i + 1, k + 1, iSymm, ei, r);
                            throw new EnergyException(message, false);
                        }
                        inducedEnergy += ei;
                    }
                }
                if (iSymm == 0) {
                    for (Atom ak : ai.get1_5s()) {
                        int index = ak.getIndex() - 1;
                        masking_local[index] = 1.0;
                        maskingp_local[index] = 1.0;
                    }
                    for (Torsion torsion : ai.getTorsions()) {
                        Atom ak = torsion.get1_4(ai);
                        if (ak != null) {
                            int index = ak.getIndex() - 1;
                            masking_local[index] = 1.0;
                            maskingp_local[index] = 1.0;
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

            // Calculate the scalar products for permanent multipoles.
            final double sc2 = dix * dkx + diy * dky + diz * dkz;
            final double sc3 = dix * xr + diy * yr + diz * zr;
            final double sc4 = dkx * xr + dky * yr + dkz * zr;
            final double sc5 = qirx * xr + qiry * yr + qirz * zr;
            final double sc6 = qkrx * xr + qkry * yr + qkrz * zr;
            final double sc7 = qirx * dkx + qiry * dky + qirz * dkz;
            final double sc8 = qkrx * dix + qkry * diy + qkrz * diz;
            final double sc9 = qirx * qkrx + qiry * qkry + qirz * qkrz;
            final double sc10 = 2.0 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) + qixx * qkxx + qiyy * qkyy + qizz * qkzz;

            // Calculate the gl functions for permanent multipoles.
            final double gl0 = ci * ck;
            final double gl1 = ck * sc3 - ci * sc4;
            final double gl2 = ci * sc6 + ck * sc5 - sc3 * sc4;
            final double gl3 = sc3 * sc6 - sc4 * sc5;
            final double gl4 = sc5 * sc6;
            final double gl5 = -4.0 * sc9;
            final double gl6 = sc2;
            final double gl7 = 2.0 * (sc7 - sc8);
            final double gl8 = 2.0 * sc10;

            // Compute the energy contributions for this interaction.
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

                // Get the permanent force with screening.
                double ftm2x = gf1 * xr + gf2 * dix + gf3 * dkx + gf4 * (qkdix - qidkx) + gf5 * qirx + gf6 * qkrx + gf7 * (qiqkrx + qkqirx);
                double ftm2y = gf1 * yr + gf2 * diy + gf3 * dky + gf4 * (qkdiy - qidky) + gf5 * qiry + gf6 * qkry + gf7 * (qiqkry + qkqiry);
                double ftm2z = gf1 * zr + gf2 * diz + gf3 * dkz + gf4 * (qkdiz - qidkz) + gf5 * qirz + gf6 * qkrz + gf7 * (qiqkrz + qkqirz);

                // Get the permanent torque with screening.
                double ttm2x = -bn1 * dixdkx + gf2 * dixrx + gf4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx) - gf5 * rxqirx - gf7 * (rxqikrx + qkrxqirx);
                double ttm2y = -bn1 * dixdky + gf2 * dixry + gf4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky) - gf5 * rxqiry - gf7 * (rxqikry + qkrxqiry);
                double ttm2z = -bn1 * dixdkz + gf2 * dixrz + gf4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz) - gf5 * rxqirz - gf7 * (rxqikrz + qkrxqirz);
                double ttm3x = bn1 * dixdkx + gf3 * dkxrx - gf4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx) - gf6 * rxqkrx - gf7 * (rxqkirx - qkrxqirx);
                double ttm3y = bn1 * dixdky + gf3 * dkxry - gf4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky) - gf6 * rxqkry - gf7 * (rxqkiry - qkrxqiry);
                double ttm3z = bn1 * dixdkz + gf3 * dkxrz - gf4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz) - gf6 * rxqkrz - gf7 * (rxqkirz - qkrxqirz);

                // Handle the case where scaling is used.
                if (scale1 != 0.0) {
                    final double gfr1 = rr3 * gl0 + rr5 * (gl1 + gl6) + rr7 * (gl2 + gl7 + gl8) + rr9 * (gl3 + gl5) + rr11 * gl4;
                    final double gfr2 = -ck * rr3 + sc4 * rr5 - sc6 * rr7;
                    final double gfr3 = ci * rr3 + sc3 * rr5 + sc5 * rr7;
                    final double gfr4 = 2.0 * rr5;
                    final double gfr5 = 2.0 * (-ck * rr5 + sc4 * rr7 - sc6 * rr9);
                    final double gfr6 = 2.0 * (-ci * rr5 - sc3 * rr7 - sc5 * rr9);
                    final double gfr7 = 4.0 * rr7;

                    // Get the permanent force without screening.
                    final double ftm2rx = gfr1 * xr + gfr2 * dix + gfr3 * dkx + gfr4 * (qkdix - qidkx) + gfr5 * qirx + gfr6 * qkrx + gfr7 * (qiqkrx + qkqirx);
                    final double ftm2ry = gfr1 * yr + gfr2 * diy + gfr3 * dky + gfr4 * (qkdiy - qidky) + gfr5 * qiry + gfr6 * qkry + gfr7 * (qiqkry + qkqiry);
                    final double ftm2rz = gfr1 * zr + gfr2 * diz + gfr3 * dkz + gfr4 * (qkdiz - qidkz) + gfr5 * qirz + gfr6 * qkrz + gfr7 * (qiqkrz + qkqirz);

                    // Get the permanent torque without screening.
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
                double prefactor = electric * selfScale * l2;
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

                // This is dU/dL/dX for the first term of dU/dL: d[dlPow * ereal]/dx
                if (lambdaTerm && soft) {
                    prefactor = electric * selfScale * dEdLSign * dlPowPerm;
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

                // Collect terms for dU/dL/dX for the second term of dU/dL: d[fL2*dfL1dL*dRealdL]/dX
                final double gf1 = bn2 * gl0 + bn3 * (gl1 + gl6)
                        + bn4 * (gl2 + gl7 + gl8)
                        + bn5 * (gl3 + gl5) + bn6 * gl4;
                final double gf2 = -ck * bn2 + sc4 * bn3 - sc6 * bn4;
                final double gf3 = ci * bn2 + sc3 * bn3 + sc5 * bn4;
                final double gf4 = 2.0 * bn3;
                final double gf5 = 2.0 * (-ck * bn3 + sc4 * bn4 - sc6 * bn5);
                final double gf6 = 2.0 * (-ci * bn3 - sc3 * bn4 - sc5 * bn5);
                final double gf7 = 4.0 * bn4;

                // Get the permanent force with screening.
                double ftm2x = gf1 * xr + gf2 * dix + gf3 * dkx
                        + gf4 * (qkdix - qidkx) + gf5 * qirx
                        + gf6 * qkrx + gf7 * (qiqkrx + qkqirx);
                double ftm2y = gf1 * yr + gf2 * diy + gf3 * dky
                        + gf4 * (qkdiy - qidky) + gf5 * qiry
                        + gf6 * qkry + gf7 * (qiqkry + qkqiry);
                double ftm2z = gf1 * zr + gf2 * diz + gf3 * dkz
                        + gf4 * (qkdiz - qidkz) + gf5 * qirz
                        + gf6 * qkrz + gf7 * (qiqkrz + qkqirz);

                // Get the permanent torque with screening.
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

                // Handle the case where scaling is used.
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

                // Add in dU/dL/dX for the second term of dU/dL: d[lPow*dlAlpha*dRealdL]/dX
                double prefactor = electric * selfScale * l2 * dlAlpha;
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
            final double usr5 = bn2 - usc5 * rr5;

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

            // Calculate the scalar products for permanent multipoles.
            final double sc3 = dix * xr + diy * yr + diz * zr;
            final double sc4 = dkx * xr + dky * yr + dkz * zr;
            final double sc5 = qirx * xr + qiry * yr + qirz * zr;
            final double sc6 = qkrx * xr + qkry * yr + qkrz * zr;

            // Calculate the scalar products for polarization components.
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

            // Calculate the gl functions for polarization components.
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

            // Compute the energy contributions for this interaction.
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

            // Get the induced force with screening.
            final double gfi1 = 0.5 * bn2 * (gli1 + glip1 + gli6 + glip6) + 0.5 * bn2 * scip2 + 0.5 * bn3 * (gli2 + glip2 + gli7 + glip7) - 0.5 * bn3 * (sci3 * scip4 + scip3 * sci4) + 0.5 * bn4 * (gli3 + glip3);
            final double gfi2 = -ck * bn1 + sc4 * bn2 - sc6 * bn3;
            final double gfi3 = ci * bn1 + sc3 * bn2 + sc5 * bn3;
            final double gfi4 = 2.0 * bn2;
            final double gfi5 = bn3 * (sci4 + scip4);
            final double gfi6 = -bn3 * (sci3 + scip3);
            double ftm2ix = gfi1 * xr + 0.5 * (gfi2 * (uix + pix) + bn2 * (sci4 * pix + scip4 * uix) + gfi3 * (ukx + pkx) + bn2 * (sci3 * pkx + scip3 * ukx) + (sci4 + scip4) * bn2 * dix + (sci3 + scip3) * bn2 * dkx + gfi4 * (qkuix + qkuipx - qiukx - qiukpx)) + gfi5 * qirx + gfi6 * qkrx;
            double ftm2iy = gfi1 * yr + 0.5 * (gfi2 * (uiy + piy) + bn2 * (sci4 * piy + scip4 * uiy) + gfi3 * (uky + pky) + bn2 * (sci3 * pky + scip3 * uky) + (sci4 + scip4) * bn2 * diy + (sci3 + scip3) * bn2 * dky + gfi4 * (qkuiy + qkuipy - qiuky - qiukpy)) + gfi5 * qiry + gfi6 * qkry;
            double ftm2iz = gfi1 * zr + 0.5 * (gfi2 * (uiz + piz) + bn2 * (sci4 * piz + scip4 * uiz) + gfi3 * (ukz + pkz) + bn2 * (sci3 * pkz + scip3 * ukz) + (sci4 + scip4) * bn2 * diz + (sci3 + scip3) * bn2 * dkz + gfi4 * (qkuiz + qkuipz - qiukz - qiukpz)) + gfi5 * qirz + gfi6 * qkrz;

            // Get the induced torque with screening.
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
                // Get the induced force without screening.
                final double gfri1 = 0.5 * rr5 * ((gli1 + gli6) * psc3 + (glip1 + glip6) * dsc3 + scip2 * usc3) + 0.5 * rr7 * ((gli7 + gli2) * psc5 + (glip7 + glip2) * dsc5 - (sci3 * scip4 + scip3 * sci4) * usc5) + 0.5 * rr9 * (gli3 * psc7 + glip3 * dsc7);
                final double gfri4 = 2.0 * rr5;
                final double gfri5 = rr7 * (sci4 * psc7 + scip4 * dsc7);
                final double gfri6 = -rr7 * (sci3 * psc7 + scip3 * dsc7);
                ftm2rix = gfri1 * xr + 0.5 * (-rr3 * ck * (uix * psc3 + pix * dsc3) + rr5 * sc4 * (uix * psc5 + pix * dsc5) - rr7 * sc6 * (uix * psc7 + pix * dsc7)) + (rr3 * ci * (ukx * psc3 + pkx * dsc3) + rr5 * sc3 * (ukx * psc5 + pkx * dsc5) + rr7 * sc5 * (ukx * psc7 + pkx * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * dix + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dkx + 0.5 * gfri4 * ((qkuix - qiukx) * psc5 + (qkuipx - qiukpx) * dsc5) + gfri5 * qirx + gfri6 * qkrx;
                ftm2riy = gfri1 * yr + 0.5 * (-rr3 * ck * (uiy * psc3 + piy * dsc3) + rr5 * sc4 * (uiy * psc5 + piy * dsc5) - rr7 * sc6 * (uiy * psc7 + piy * dsc7)) + (rr3 * ci * (uky * psc3 + pky * dsc3) + rr5 * sc3 * (uky * psc5 + pky * dsc5) + rr7 * sc5 * (uky * psc7 + pky * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * diy + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dky + 0.5 * gfri4 * ((qkuiy - qiuky) * psc5 + (qkuipy - qiukpy) * dsc5) + gfri5 * qiry + gfri6 * qkry;
                ftm2riz = gfri1 * zr + 0.5 * (-rr3 * ck * (uiz * psc3 + piz * dsc3) + rr5 * sc4 * (uiz * psc5 + piz * dsc5) - rr7 * sc6 * (uiz * psc7 + piz * dsc7)) + (rr3 * ci * (ukz * psc3 + pkz * dsc3) + rr5 * sc3 * (ukz * psc5 + pkz * dsc5) + rr7 * sc5 * (ukz * psc7 + pkz * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * diz + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dkz + 0.5 * gfri4 * ((qkuiz - qiukz) * psc5 + (qkuipz - qiukpz) * dsc5) + gfri5 * qirz + gfri6 * qkrz;

                // Get the induced torque without screening.
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

            // Account for partially excluded induced interactions.
            double temp3 = 0.5 * rr3 * ((gli1 + gli6) * scalep + (glip1 + glip6) * scaled);
            double temp5 = 0.5 * rr5 * ((gli2 + gli7) * scalep + (glip2 + glip7) * scaled);
            final double temp7 = 0.5 * rr7 * (gli3 * scalep + glip3 * scaled);
            final double fridmpx = temp3 * ddsc3x + temp5 * ddsc5x + temp7 * ddsc7x;
            final double fridmpy = temp3 * ddsc3y + temp5 * ddsc5y + temp7 * ddsc7y;
            final double fridmpz = temp3 * ddsc3z + temp5 * ddsc5z + temp7 * ddsc7z;

            // Find some scaling terms for induced-induced force.
            temp3 = 0.5 * rr3 * scip2;
            temp5 = -0.5 * rr5 * (sci3 * scip4 + scip3 * sci4);
            final double findmpx = temp3 * ddsc3x + temp5 * ddsc5x;
            final double findmpy = temp3 * ddsc3y + temp5 * ddsc5y;
            final double findmpz = temp3 * ddsc3z + temp5 * ddsc5z;

            // Modify the forces for partially excluded interactions.
            ftm2ix = ftm2ix - fridmpx - findmpx;
            ftm2iy = ftm2iy - fridmpy - findmpy;
            ftm2iz = ftm2iz - fridmpz - findmpz;

            // Correction to convert mutual to direct polarization force.
            if (polarization == ParticleMeshEwald.Polarization.DIRECT) {
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

            // Correction for OPT induced dipoles.
//            if (scfAlgorithm == ParticleMeshEwald.SCFAlgorithm.EPT) {
//                double eptx = 0.0;
//                double epty = 0.0;
//                double eptz = 0.0;
//                for (int jj = 0; jj < optOrder; jj++) {
//                    double optix = optDipole[jj][i][0];
//                    double optiy = optDipole[jj][i][1];
//                    double optiz = optDipole[jj][i][2];
//                    double optixp = optDipoleCR[jj][i][0];
//                    double optiyp = optDipoleCR[jj][i][1];
//                    double optizp = optDipoleCR[jj][i][2];
//                    double uirm = optix * xr + optiy * yr + optiz * zr;
//                    for (int mm = 0; mm < optOrder - jj; mm++) {
//                        double optkx = optDipole[mm][k][0];
//                        double optky = optDipole[mm][k][1];
//                        double optkz = optDipole[mm][k][2];
//                        double optkxp = optDipoleCR[mm][k][0];
//                        double optkyp = optDipoleCR[mm][k][1];
//                        double optkzp = optDipoleCR[mm][k][2];
//                        double ukrm = optkx * xr + optky * yr + optkz * zr;
//                        double term1 = bn2 - usc3 * rr5;
//                        double term2 = bn3 - usc5 * rr7;
//                        double term3 = usr5 + term1;
//                        double term4 = rr3;
//                        double term5 = -xr * term3 + ddsc3x * term4;
//                        double term6 = -(bn2 - usc5 * rr5) + xr * xr * term2 - rr5 * xr * ddsc5x;
//                        double tixx = optix * term5 + uirm * term6;
//                        double tkxx = optkx * term5 + ukrm * term6;
//                        term5 = -yr * term3 + ddsc3y * term4;
//                        term6 = -usr5 + yr * yr * term2 - rr5 * yr * ddsc5y;
//                        double tiyy = optiy * term5 + uirm * term6;
//                        double tkyy = optky * term5 + ukrm * term6;
//                        term5 = -zr * term3 + ddsc3z * term4;
//                        term6 = -usr5 + zr * zr * term2 - rr5 * zr * ddsc5z;
//                        double tizz = optiz * term5 + uirm * term6;
//                        double tkzz = optkz * term5 + ukrm * term6;
//                        term4 = -usr5 * yr;
//                        term5 = -xr * term1 + rr3 * ddsc3x;
//                        term6 = xr * yr * term2 - rr5 * yr * ddsc5x;
//                        double tixy = optix * term4 + optiy * term5 + uirm * term6;
//                        double tkxy = optkx * term4 + optky * term5 + ukrm * term6;
//                        term4 = -usr5 * zr;
//                        term6 = xr * zr * term2 - rr5 * zr * ddsc5x;
//                        double tixz = optix * term4 + optiz * term5 + uirm * term6;
//                        double tkxz = optkx * term4 + optkz * term5 + ukrm * term6;
//                        term5 = -yr * term1 + rr3 * ddsc3y;
//                        term6 = yr * zr * term2 - rr5 * zr * ddsc5y;
//                        double tiyz = optiy * term4 + optiz * term5 + uirm * term6;
//                        double tkyz = optkz * term4 + optkz * term5 + ukrm * term6;
//                        double depx = tixx * optkxp + tkxx * optixp
//                                + tixy * optkyp + tkxy * optiyp
//                                + tixz * optkzp + tkxz * optizp;
//                        double depy = tixy * optkxp + tkxy * optixp
//                                + tiyy * optkyp + tkyy * optiyp
//                                + tiyz * optkzp + tkyz * optizp;
//                        double depz = tixz * optkxp + tkxz * optixp
//                                + tiyz * optkyp + tkyz * optiyp
//                                + tizz * optkzp + tkzz * optizp;
//                        double optCoefSum = optRegion.optCoefficientsSum[jj + mm + 1];
//                        double fx = optCoefSum * depx;
//                        double fy = optCoefSum * depy;
//                        double fz = optCoefSum * depz;
//                        eptx += fx;
//                        epty += fy;
//                        eptz += fz;
//                    }
//                }
//                if (i == 0 && k == 1) {
//                    logger.info(format(" %s %16.8f %16.8f %16.8f", "ept force", eptx, epty, eptz));
//                }
//                ftm2ix += eptx;
//                ftm2iy += epty;
//                ftm2iz += eptz;
//            }

            // Handle the case where scaling is used.
            ftm2ix = ftm2ix - ftm2rix;
            ftm2iy = ftm2iy - ftm2riy;
            ftm2iz = ftm2iz - ftm2riz;
            ttm2ix = ttm2ix - ttm2rix;
            ttm2iy = ttm2iy - ttm2riy;
            ttm2iz = ttm2iz - ttm2riz;
            ttm3ix = ttm3ix - ttm3rix;
            ttm3iy = ttm3iy - ttm3riy;
            ttm3iz = ttm3iz - ttm3riz;

            double scalar = electric * polarizationScale * selfScale;
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
                scalar = electric * dEdLSign * dlPowPol * selfScale;
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

    private void torque(int iSymm, double[] tx, double[] ty, double[] tz, double[] gx, double[] gy, double[] gz,
                        double[] origin, double u[],
                        double[] v, double[] w, double[] uv, double[] uw,
                        double[] vw, double[] ur, double[] us, double[] vs,
                        double[] ws, double[] t1, double[] t2, double[] r, double[] s) {
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            final int[] ax = axisAtom[i];
            // Ions, for example, have no torque.
            if (ax == null || ax.length < 2) {
                continue;
            }
            final int ia = ax[0];
            final int ib = i;
            final int ic = ax[1];
            int id = 0;
            double[] x = coordinates[iSymm][0];
            double[] y = coordinates[iSymm][1];
            double[] z = coordinates[iSymm][2];
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
}