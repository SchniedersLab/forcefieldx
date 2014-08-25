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

import java.nio.DoubleBuffer;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.numerics.TensorRecursion;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Complex3DCuda;
import ffx.numerics.fft.Complex3DParallel;
import ffx.numerics.fft.Real3DParallel;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;

import static ffx.crystal.Crystal.mod;
import static ffx.numerics.UniformBSpline.bSpline;
import static ffx.numerics.UniformBSpline.bSplineDerivatives;
import static ffx.numerics.fft.Complex3D.iComplex3D;
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
 * The Reciprocal Space class computes the reciprocal space contribution to
 * {@link ParticleMeshEwald} for the AMOEBA force field.
 *
 * <ol>
 * <li>
 * Assignment of polarizable multipole charge density to the 3D grid, via
 * b-Splines, is parallelized using a spatial decomposition. </li>
 * <li>
 * The convolution depends on methods of the {@link Real3DParallel} and
 * {@link Complex3DParallel} classes.
 * </li>
 * <li> Finally, the electric potential and its gradients are collected, in
 * parallel, off the grid using b-Splines.
 * </li>
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class ReciprocalSpace {

    private static final Logger logger = Logger.getLogger(ReciprocalSpace.class.getName());
    private final ParticleMeshEwald particleMeshEwald;
    private final ForceField forceField;
    /**
     * The unit cell for the simulation. A ReplicatesCrystal will give
     * equivalent results up to the limits of double precision, but is more
     * expensive.
     */
    private Crystal crystal;
    /**
     * Array of atoms in the asymmetric unit.
     */
    private Atom atoms[];
    /**
     * Number of atoms in the asymmetric unit.
     */
    private int nAtoms;
    /**
     * Number of unit cell symmetry operators.
     */
    private final int nSymm;
    /**
     * Number of threads.
     */
    private final int threadCount;
    /**
     * The b-Spline order to use for discretization to/from the reciprocal grid.
     */
    private final int bSplineOrder;
    /**
     * Three derivatives of the potential are needed for AMOEBA. Specifically,
     * the field gradient is used to compute the energy of the quadrupole moment
     * and the 2nd field gradient (3rd derivative of the reciprocal potential)
     * the energy gradient.
     */
    private final int derivOrder = 3;
    /**
     * The X-dimension of the FFT grid.
     */
    private int fftX;
    /**
     * The Y-dimension of the FFT grid.
     */
    private int fftY;
    /**
     * The Z-dimension of the FFT grid.
     */
    private int fftZ;
    /**
     * Number of doubles needed to compute a complex to complex 3D FFT (fftX *
     * fftY * fftZ * 2).
     */
    private int fftSpace;
    /**
     * Reciprocal space grid. [fftSpace]
     */
    private double splineGrid[];
    /**
     * Wraps the splineGrid.
     */
    private DoubleBuffer splineBuffer;
    /**
     * Ewald convergence parameters.
     */
    private final double aEwald;
    /**
     * Reference to atomic coordinate array.
     */
    private double coordinates[][][];
    /**
     * Fractional multipole array.
     */
    private double fracMultipole[][][];
    /**
     * Fractional induced dipole array.
     */
    private double fracInducedDipole[][][];
    /**
     * Fractional induced dipole chain rule array.
     */
    private double fracInducedDipoleCR[][][];
    /**
     * Fractional multipole potential and its derivatives.
     */
    private double fracMultipolePhi[][];
    /**
     * Fractional induced dipole potential and its derivatives.
     */
    private double fracInducedDipolePhi[][];
    /**
     * Fractional induced dipole chain rule potential and its derivatives.
     */
    private double fracInducedDipolePhiCR[][];
    /**
     * Parallel team instance.
     */
    private final ParallelTeam parallelTeam;
    private final ParallelTeam fftTeam;
    private final BSplineRegion bSplineRegion;

    private SpatialDensityRegion spatialDensityRegion;
    private final SpatialPermanentLoop spatialPermanentLoops[];
    private final SpatialInducedLoop spatialInducedLoops[];

    private SliceRegion sliceRegion;
    private final SlicePermanentLoop slicePermanentLoops[];
    private final SliceInducedLoop sliceInducedLoops[];
    /**
     * Number of atoms for a given symmetry operator that a given thread is
     * responsible for applying to the FFT grid. gridAtomCount[nSymm][nThread]
     */
    private final int gridAtomCount[][];
    /**
     * Atom indices for a given symmetry operator that a given thread is
     * responsible for applying to the FFT grid.
     * gridAtomCount[nSymm][nThread][nAtoms]
     */
    private final int gridAtomList[][][];

    private final PermanentPhiRegion permanentPhiRegion;

    private final InducedPhiRegion polarizationPhiRegion;
    private Complex3DParallel complexFFT3D;
    private final IntegerSchedule recipSchedule;
    /**
     * Timing variables.
     */
    private final long bSplineTime[];
    private final long splinePermanentTime[];
    private final long permanentPhiTime[];
    private final long splineInducedTime[];
    private final int splineCount[];
    private final long inducedPhiTime[];
    private long bSplineTotal, splinePermanentTotal, splineInducedTotal;
    private long permanentPhiTotal, inducedPhiTotal, convTotal;
    /**
     * CUDA convolution variables.
     */
    private final boolean cudaFFT;
    private Thread cudaThread;
    private Complex3DCuda cudaFFT3D;
    private GridMethod gridMethod = GridMethod.SPATIAL;

    enum GridMethod {

        SPATIAL, SLICE
    };

    /**
     * Reciprocal Space PME contribution.
     *
     * @param particleMeshEwald a
     * {@link ffx.potential.nonbonded.ParticleMeshEwald} object.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     * @param aewald the Ewald parameter.
     * @param fftTeam a {@link edu.rit.pj.ParallelTeam} object.
     * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
     */
    public ReciprocalSpace(ParticleMeshEwald particleMeshEwald,
            Crystal crystal, ForceField forceField, Atom atoms[],
            double aewald, ParallelTeam fftTeam, ParallelTeam parallelTeam) {

        this.particleMeshEwald = particleMeshEwald;
        this.crystal = crystal.getUnitCell();
        this.forceField = forceField;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.aEwald = aewald;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;

        coordinates = particleMeshEwald.coordinates;
        threadCount = parallelTeam.getThreadCount();
        nSymm = this.crystal.spaceGroup.getNumberOfSymOps();

        /**
         * Construct the parallel convolution object.
         */
        boolean available;
        String recipStrategy = null;
        try {
            recipStrategy = forceField.getString(ForceField.ForceFieldString.RECIP_SCHEDULE);
            IntegerSchedule.parse(recipStrategy);
            available = true;
        } catch (Exception e) {
            available = false;
        }
        if (available) {
            recipSchedule = IntegerSchedule.parse(recipStrategy);
            logger.log(Level.INFO, "  Convolution schedule {0}", recipStrategy);
        } else {
            recipSchedule = IntegerSchedule.fixed();
        }
        cudaFFT = forceField.getBoolean(ForceField.ForceFieldBoolean.CUDAFFT, false);
        bSplineOrder = forceField.getInteger(ForceFieldInteger.PME_ORDER, 5);

        /**
         * Initialize convolution objects that may be re-allocated during NPT
         * simulations.
         */
        double density = initConvolution();

        /**
         * Log PME reciprocal space parameters.
         */
        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("   B-Spline Order:                     %8d\n", bSplineOrder));
            sb.append(format("   Mesh Density:                       %8.3f\n", density));
            sb.append(format("   Mesh Dimensions:               (%3d,%3d,%3d)", fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        initAtomArrays();

        /**
         * Construct the parallel BSplineRegion, DensityLoops and Phi objects.
         */
        bSplineRegion = new BSplineRegion();
        switch (gridMethod) {
            case SPATIAL:
                spatialPermanentLoops = new SpatialPermanentLoop[threadCount];
                spatialInducedLoops = new SpatialInducedLoop[threadCount];
                for (int i = 0; i < threadCount; i++) {
                    spatialPermanentLoops[i] = new SpatialPermanentLoop(spatialDensityRegion, bSplineRegion);
                    spatialInducedLoops[i] = new SpatialInducedLoop(spatialDensityRegion, bSplineRegion);
                }
                slicePermanentLoops = null;
                sliceInducedLoops = null;
                gridAtomCount = null;
                gridAtomList = null;
                break;
            case SLICE:
            default:
                slicePermanentLoops = new SlicePermanentLoop[threadCount];
                sliceInducedLoops = new SliceInducedLoop[threadCount];
                for (int i = 0; i < threadCount; i++) {
                    slicePermanentLoops[i] = new SlicePermanentLoop(sliceRegion, bSplineRegion);
                    sliceInducedLoops[i] = new SliceInducedLoop(sliceRegion, bSplineRegion);
                }
                gridAtomCount = new int[nSymm][threadCount];
                gridAtomList = new int[nSymm][threadCount][nAtoms];
                spatialPermanentLoops = null;
                spatialInducedLoops = null;
        }
        permanentPhiRegion = new PermanentPhiRegion(bSplineRegion);
        polarizationPhiRegion = new InducedPhiRegion(bSplineRegion);

        /**
         * Initialize timing variables.
         */
        bSplineTime = new long[threadCount];
        splinePermanentTime = new long[threadCount];
        splineInducedTime = new long[threadCount];
        splineCount = new int[threadCount];
        permanentPhiTime = new long[threadCount];
        inducedPhiTime = new long[threadCount];
    }

    public void setAtoms(Atom atoms[]) {
        this.atoms = atoms;
        nAtoms = atoms.length;
        initAtomArrays();
    }

    private void initAtomArrays() {
        if (fracMultipolePhi == null || fracMultipolePhi.length < nAtoms) {
            /**
             * Allocate memory for fractional multipoles and induced dipoles.
             */
            fracMultipole = new double[nSymm][nAtoms][10];
            fracInducedDipole = new double[nSymm][nAtoms][3];
            fracInducedDipoleCR = new double[nSymm][nAtoms][3];
            /**
             * Allocate memory for fractional permanent Phi and induced Phi.
             */
            fracMultipolePhi = new double[nAtoms][tensorCount];
            fracInducedDipolePhi = new double[nAtoms][tensorCount];
            fracInducedDipolePhiCR = new double[nAtoms][tensorCount];
        }
    }

    public void setCrystal(Crystal crystal) {
        /**
         * Check if the number of symmetry operators has changed.
         */
        this.crystal = crystal.getUnitCell();
        if (nSymm != this.crystal.spaceGroup.getNumberOfSymOps()) {
            logger.info(this.crystal.toString());
            logger.info(crystal.toString());
            logger.severe(" The reciprocal space class does not currently allow changes in the number of symmetry operators.");
        }
        this.coordinates = particleMeshEwald.coordinates;
        initConvolution();
    }

    private double initConvolution() {

        /**
         * Store the current reciprocal space grid dimensions.
         */
        int fftXCurrent = fftX;
        int fftYCurrent = fftY;
        int fftZCurrent = fftZ;

        double density = forceField.getDouble(ForceFieldDouble.PME_MESH_DENSITY, 1.2);

        int nX = forceField.getInteger(ForceFieldInteger.PME_GRID_X, -1);
        if (nX < 2) {
            nX = (int) Math.floor(crystal.a * density) + 1;
            if (nX % 2 != 0) {
                nX += 1;
            }
            while (!Complex.preferredDimension(nX)) {
                nX += 2;
            }
        }
        int nY = forceField.getInteger(ForceFieldInteger.PME_GRID_Y, -1);
        if (nY < 2) {
            nY = (int) Math.floor(crystal.b * density) + 1;
            if (nY % 2 != 0) {
                nY += 1;
            }
            while (!Complex.preferredDimension(nY)) {
                nY += 2;
            }
        }
        int nZ = forceField.getInteger(ForceFieldInteger.PME_GRID_Z, -1);
        if (nZ < 2) {
            nZ = (int) Math.floor(crystal.c * density) + 1;
            if (nZ % 2 != 0) {
                nZ += 1;
            }
            while (!Complex.preferredDimension(nZ)) {
                nZ += 2;
            }
        }

        fftX = nX;
        fftY = nY;
        fftZ = nZ;

        /**
         * Populate the matrix that fractionalizes multipoles.
         */
        transformMultipoleMatrix();
        /**
         * Populate the matrix that convert fractional potential components into
         * orthogonal Cartesian coordinates.
         */
        transformFieldMatrix();
        /**
         * Compute the Cartesian to fractional matrix.
         */
        for (int i = 0; i < 3; i++) {
            a[0][i] = fftX * crystal.A[i][0];
            a[1][i] = fftY * crystal.A[i][1];
            a[2][i] = fftZ * crystal.A[i][2];
        }

        fftSpace = fftX * fftY * fftZ * 2;
        boolean dimChanged = fftX != fftXCurrent || fftY != fftYCurrent || fftZ != fftZCurrent;

        if (!cudaFFT) {
            if (complexFFT3D == null || dimChanged) {
                complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam, recipSchedule);
                if (splineGrid == null || splineGrid.length < fftSpace) {
                    splineGrid = new double[fftSpace];
                }
                splineBuffer = DoubleBuffer.wrap(splineGrid);
            }
            complexFFT3D.setRecip(generalizedInfluenceFunction());
            cudaFFT3D = null;
            cudaThread = null;
        } else {
            if (cudaFFT3D == null || dimChanged) {
                if (cudaFFT3D != null) {
                    cudaFFT3D.free();
                }
                cudaFFT3D = new Complex3DCuda(fftX, fftY, fftZ);
                cudaThread = new Thread(cudaFFT3D);
                cudaThread.setPriority(Thread.MAX_PRIORITY);
                cudaThread.start();
                splineBuffer = cudaFFT3D.getDoubleBuffer();
            }
            cudaFFT3D.setRecip(generalizedInfluenceFunction());
            complexFFT3D = null;
        }

        switch (gridMethod) {
            case SPATIAL:
                if (spatialDensityRegion == null || dimChanged) {
                    spatialDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ, splineGrid, bSplineOrder, nSymm,
                            10, threadCount, crystal, atoms, coordinates);
                    if (cudaFFT) {
                        spatialDensityRegion.setGridBuffer(splineBuffer);
                    }
                } else {
                    spatialDensityRegion.setCrystal(crystal, fftX, fftY, fftZ);
                    spatialDensityRegion.coordinates = coordinates;
                }
                break;
            case SLICE:
            default:
                if (sliceRegion == null || dimChanged) {
                    sliceRegion = new SliceRegion(fftX, fftY, fftZ, splineGrid, bSplineOrder, nSymm,
                            threadCount, crystal, atoms, coordinates);
                    if (cudaFFT) {
                        sliceRegion.setGridBuffer(splineBuffer);
                    }
                } else {
                    sliceRegion.setCrystal(crystal, fftX, fftY, fftZ);
                    sliceRegion.coordinates = coordinates;
                }
        }
        return density;
    }

    public void printTimings() {
        if (logger.isLoggable(Level.FINE)) {
            if (complexFFT3D != null) {
                double total = (bSplineTotal + convTotal + splinePermanentTotal + permanentPhiTotal
                        + splineInducedTotal + inducedPhiTotal) * toSeconds;

                logger.fine(String.format("\n Reciprocal Space: %7.4f (sec)", total));
                long convTime[] = complexFFT3D.getTimings();
                logger.fine("                           Direct Field    SCF Field");
                logger.fine(" Thread  B-Spline  3DConv  Spline  Phi     Spline  Phi      Count");

                long minBSpline = Long.MAX_VALUE;
                long maxBSpline = 0;
                long minConv = Long.MAX_VALUE;
                long maxConv = 0;
                long minBSPerm = Long.MAX_VALUE;
                long maxBSPerm = 0;
                long minPhiPerm = Long.MAX_VALUE;
                long maxPhiPerm = 0;
                long minBSInduced = Long.MAX_VALUE;
                long maxBSInduced = 0;
                long minPhiInduced = Long.MAX_VALUE;
                long maxPhiInduced = 0;
                int minCount = Integer.MAX_VALUE;
                int maxCount = 0;

                for (int i = 0; i < threadCount; i++) {
                    logger.fine(String.format("    %3d   %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f  %6d", i,
                            bSplineTime[i] * toSeconds, convTime[i] * toSeconds, splinePermanentTime[i] * toSeconds,
                            permanentPhiTime[i] * toSeconds, splineInducedTime[i] * toSeconds,
                            inducedPhiTime[i] * toSeconds, splineCount[i]));
                    minBSpline = min(bSplineTime[i],minBSpline);
                    maxBSpline = max(bSplineTime[i],maxBSpline);
                    minConv = min(convTime[i],minConv);
                    maxConv = max(convTime[i],maxConv);
                    minBSPerm = min(splinePermanentTime[i],minBSPerm);
                    maxBSPerm = max(splinePermanentTime[i],maxBSPerm);
                    minPhiPerm = min(permanentPhiTime[i],minPhiPerm);
                    maxPhiPerm = max(permanentPhiTime[i],maxPhiPerm);
                    minBSInduced = min(splineInducedTime[i],minBSInduced);
                    maxBSInduced = max(splineInducedTime[i],maxBSInduced);
                    minPhiInduced = min(inducedPhiTime[i],minPhiInduced);
                    maxPhiInduced = max(inducedPhiTime[i],maxPhiInduced);
                    minCount = min(splineCount[i],minCount);
                    maxCount = max(splineCount[i],maxCount);
                }
                logger.fine(String.format(" Min      %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f  %6d",
                        minBSpline * toSeconds, minConv * toSeconds, minBSPerm * toSeconds,
                        minPhiPerm * toSeconds, minBSInduced * toSeconds,
                        minPhiInduced * toSeconds, minCount));
                logger.fine(String.format(" Max      %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f  %6d",
                        maxBSpline * toSeconds, maxConv * toSeconds, maxBSPerm * toSeconds,
                        maxPhiPerm * toSeconds, maxBSInduced * toSeconds,
                        maxPhiInduced * toSeconds, maxCount));
                logger.fine(String.format(" Delta    %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f  %6d",
                        (maxBSpline-minBSpline) * toSeconds, (maxConv-minConv) * toSeconds,
                        (maxBSPerm-minBSPerm) * toSeconds, (maxPhiPerm-minPhiPerm) * toSeconds,
                        (maxBSInduced-minBSInduced) * toSeconds, (maxPhiInduced-minPhiInduced) * toSeconds,
                        (maxCount-minCount)));
                logger.fine(String.format(" Actual   %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f  %6d",
                        bSplineTotal * toSeconds, convTotal * toSeconds, splinePermanentTotal * toSeconds,
                        permanentPhiTotal * toSeconds, splineInducedTotal * toSeconds,
                        inducedPhiTotal * toSeconds, nAtoms * nSymm));
            }
        }
    }

    public void initTimings() {
        /**
         * Reset total timings.
         */
        bSplineTotal = 0;
        splinePermanentTotal = 0;
        splineInducedTotal = 0;
        permanentPhiTotal = 0;
        inducedPhiTotal = 0;
        convTotal = 0;

        /**
         * Reset timing arrays.
         */
        for (int i = 0; i < threadCount; i++) {
            bSplineTime[i] = 0;
            splinePermanentTime[i] = 0;
            splineInducedTime[i] = 0;
            permanentPhiTime[i] = 0;
            inducedPhiTime[i] = 0;
            splineCount[i] = 0;
        }

        /**
         * Reset 3D convolution timing.
         */
        if (complexFFT3D != null) {
            complexFFT3D.initTiming();
        }
    }

    /**
     * <p>
     * computeBSplines</p>
     */
    public void computeBSplines() {
        try {
            bSplineTotal -= System.nanoTime();
            parallelTeam.execute(bSplineRegion);
            bSplineTotal += System.nanoTime();
        } catch (Exception e) {
            String message = " Fatal exception evaluating b-Splines.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Use b-Splines to place the permanent multipoles onto the FFT grid for the
     * atoms in use.
     *
     * @param globalMultipoles an array of double.
     * @param use an array of boolean.
     */
    public void splinePermanentMultipoles(double globalMultipoles[][][], boolean use[]) {
        splinePermanentTotal -= System.nanoTime();
        switch (gridMethod) {
            case SPATIAL:
                spatialDensityRegion.setCrystal(crystal.getUnitCell(), fftX, fftY, fftZ);
                spatialDensityRegion.assignAtomsToCells();
                spatialDensityRegion.setDensityLoop(spatialPermanentLoops);
                for (int i = 0; i < threadCount; i++) {
                    spatialPermanentLoops[i].setPermanent(globalMultipoles);
                    spatialPermanentLoops[i].setUse(use);
                    spatialPermanentLoops[i].setRegion(spatialDensityRegion);
                }
                try {
                    parallelTeam.execute(spatialDensityRegion);
                } catch (Exception e) {
                    String message = " Fatal exception evaluating permanent multipole density.";
                    logger.log(Level.SEVERE, message, e);
                }
                break;
            case SLICE:
            default:
                sliceRegion.setCrystal(crystal.getUnitCell(), fftX, fftY, fftZ);
                sliceRegion.setDensityLoop(slicePermanentLoops);
                for (int i = 0; i < threadCount; i++) {
                    slicePermanentLoops[i].setPermanent(globalMultipoles);
                    slicePermanentLoops[i].setUse(use);
                }
                try {
                    parallelTeam.execute(sliceRegion);
                } catch (Exception e) {
                    String message = " Fatal exception evaluating permanent multipole density.";
                    logger.log(Level.SEVERE, message, e);
                }
                break;

        }
        splinePermanentTotal += System.nanoTime();
    }

    /**
     * <p>
     * permanentMultipoleConvolution</p>
     */
    public void permanentMultipoleConvolution() {
        convTotal -= System.nanoTime();
        try {
            if (cudaFFT) {
                cudaFFT3D.convolution(splineGrid);
            } else {
                complexFFT3D.convolution(splineGrid);
            }
        } catch (Exception e) {
            String message = " Fatal exception evaluating permanent convolution.";
            logger.log(Level.SEVERE, message, e);
        }
        convTotal += System.nanoTime();
    }

    /**
     * Compute the potential Phi and its derivatives for all atoms.
     *
     * @param cartPermanentPhi an array of double.
     */
    public void computePermanentPhi(double cartPermanentPhi[][]) {
        permanentPhiTotal -= System.nanoTime();
        try {
            permanentPhiRegion.setCartPermanentPhi(cartPermanentPhi);
            parallelTeam.execute(permanentPhiRegion);
        } catch (Exception e) {
            String message = " Fatal exception evaluating permanent reciprocal space potential.";
            logger.log(Level.SEVERE, message, e);
        }
        permanentPhiTotal += System.nanoTime();
    }

    /**
     * Place the induced dipoles onto the FFT grid for the atoms in use.
     *
     * @param inducedDipole Induced dipoles.
     * @param inducedDipoleCR Chain rule term for induced dipole gradient.
     * @param use The atoms in use.
     */
    public void splineInducedDipoles(double inducedDipole[][][],
            double inducedDipoleCR[][][],
            boolean use[]) {
        splineInducedTotal -= System.nanoTime();
        switch (gridMethod) {
            case SPATIAL:
                spatialDensityRegion.setDensityLoop(spatialInducedLoops);
                for (int i = 0; i < threadCount; i++) {
                    spatialInducedLoops[i].setInducedDipoles(inducedDipole, inducedDipoleCR);
                    spatialInducedLoops[i].setUse(use);
                    spatialInducedLoops[i].setRegion(spatialDensityRegion);
                }
                try {
                    parallelTeam.execute(spatialDensityRegion);
                } catch (Exception e) {
                    String message = " Fatal exception evaluating induced density.\n";
                    logger.log(Level.SEVERE, message, e);
                }
                break;
            case SLICE:
            default:
                sliceRegion.setDensityLoop(sliceInducedLoops);
                for (int i = 0; i < threadCount; i++) {
                    sliceInducedLoops[i].setInducedDipoles(inducedDipole, inducedDipoleCR);
                    sliceInducedLoops[i].setUse(use);
                }
                try {
                    parallelTeam.execute(sliceRegion);
                } catch (Exception e) {
                    String message = " Fatal exception evaluating permanent multipole density.";
                    logger.log(Level.SEVERE, message, e);
                }
                break;

        }
        splineInducedTotal += System.nanoTime();
    }

    /**
     * Note that the Java function "signum" and the FORTRAN version have
     * different definitions for an argument of zero.
     * <p>
     * JAVA: Math.signum(0.0) == 0.0
     * <p>
     * FORTRAN: signum(0.0) .eq. 1.0
     */
    public void inducedDipoleConvolution() {
        convTotal -= System.nanoTime();
        try {
            if (cudaFFT) {
                cudaFFT3D.convolution(splineGrid);
            } else {
                complexFFT3D.convolution(splineGrid);
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating induced convolution.";
            logger.log(Level.SEVERE, message, e);
        }
        convTotal += System.nanoTime();
    }

    /**
     * <p>
     * computeInducedPhi</p>
     *
     * @param cartInducedDipolePhi an array of double.
     * @param cartInducedDipoleCRPhi an array of double.
     */
    public void computeInducedPhi(double cartInducedDipolePhi[][],
            double cartInducedDipoleCRPhi[][]) {
        inducedPhiTotal -= System.nanoTime();
        try {
            polarizationPhiRegion.setCartInducedDipolePhi(cartInducedDipolePhi,
                    cartInducedDipoleCRPhi);
            parallelTeam.execute(polarizationPhiRegion);
        } catch (Exception e) {
            String message = "Fatal exception evaluating induced reciprocal space potential.";
            logger.log(Level.SEVERE, message, e);
        }
        inducedPhiTotal += System.nanoTime();
    }

    /**
     * <p>
     * cartToFracInducedDipoles</p>
     *
     * @param inducedDipole an array of double.
     * @param inducedDipoleCR an array of double.
     */
    public void cartToFracInducedDipoles(double inducedDipole[][][],
            double inducedDipoleCR[][][]) {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            for (int i = 0; i < nAtoms; i++) {
                double in[] = inducedDipole[iSymm][i];
                double out[] = fracInducedDipole[iSymm][i];
                out[0] = a[0][0] * in[0] + a[0][1] * in[1] + a[0][2] * in[2];
                out[1] = a[1][0] * in[0] + a[1][1] * in[1] + a[1][2] * in[2];
                out[2] = a[2][0] * in[0] + a[2][1] * in[1] + a[2][2] * in[2];
                in = inducedDipoleCR[iSymm][i];
                out = fracInducedDipoleCR[iSymm][i];
                out[0] = a[0][0] * in[0] + a[0][1] * in[1] + a[0][2] * in[2];
                out[1] = a[1][0] * in[0] + a[1][1] * in[1] + a[1][2] * in[2];
                out[2] = a[2][0] * in[0] + a[2][1] * in[1] + a[2][2] * in[2];
            }
        }
    }

    /**
     * <p>
     * Getter for the field <code>fracMultipolePhi</code>.</p>
     *
     * @return an array of double.
     */
    public double[][] getFracMultipolePhi() {
        return fracMultipolePhi;
    }

    /**
     * <p>
     * getFracMultipoles</p>
     *
     * @return an array of double.
     */
    public double[][] getFracMultipoles() {
        return fracMultipole[0];
    }

    /**
     * <p>
     * Getter for the field <code>fracInducedDipolePhi</code>.</p>
     *
     * @return an array of double.
     */
    public double[][] getFracInducedDipolePhi() {
        return fracInducedDipolePhi;
    }

    /**
     * <p>
     * getFracInducedDipoles</p>
     *
     * @return an array of double.
     */
    public double[][] getFracInducedDipoles() {
        return fracInducedDipole[0];
    }

    /**
     * <p>
     * getFracInducedDipoleCRPhi</p>
     *
     * @return an array of double.
     */
    public double[][] getFracInducedDipoleCRPhi() {
        return fracInducedDipolePhiCR;
    }

    /**
     * <p>
     * getFracInducedDipolesCR</p>
     *
     * @return an array of double.
     */
    public double[][] getFracInducedDipolesCR() {
        return fracInducedDipoleCR[0];
    }

    /**
     * <p>
     * getXDim</p>
     *
     * @return a double.
     */
    public double getXDim() {
        return fftX;
    }

    /**
     * <p>
     * getYDim</p>
     *
     * @return a double.
     */
    public double getYDim() {
        return fftY;
    }

    /**
     * <p>
     * getZDim</p>
     *
     * @return a double.
     */
    public double getZDim() {
        return fftZ;
    }

    /**
     * The BSplineRegion computes the b-Splines used spline multipoles onto the
     * reciprocal space grid and spline the reciprocal potential from the
     * reciprocal space grid. An instance of BSplineRegion automatically updates
     * itself to be consistent with the current
     */
    private class BSplineRegion extends ParallelRegion {

        private double r00;
        private double r01;
        private double r02;
        private double r10;
        private double r11;
        private double r12;
        private double r20;
        private double r21;
        private double r22;
        private final BSplineFillLoop bSplineFillLoop[];
        public double splineX[][][][];
        public double splineY[][][][];
        public double splineZ[][][][];
        public int initGrid[][][];

        public BSplineRegion() {
            bSplineFillLoop = new BSplineFillLoop[threadCount];
        }

        @Override
        public void start() {
            r00 = crystal.A[0][0];
            r01 = crystal.A[0][1];
            r02 = crystal.A[0][2];
            r10 = crystal.A[1][0];
            r11 = crystal.A[1][1];
            r12 = crystal.A[1][2];
            r20 = crystal.A[2][0];
            r21 = crystal.A[2][1];
            r22 = crystal.A[2][2];
            if (splineX == null || splineX[0].length < nAtoms) {
                initGrid = new int[nSymm][nAtoms][];
                splineX = new double[nSymm][nAtoms][][];
                splineY = new double[nSymm][nAtoms][][];
                splineZ = new double[nSymm][nAtoms][][];
            }
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            if (bSplineFillLoop[threadIndex] == null) {
                bSplineFillLoop[threadIndex] = new BSplineFillLoop();
            }
            /**
             * Currently this condition would indicate a programming bug, since
             * the space group is not allowed to change and the ReciprocalSpace
             * class should be operating on a unit cell with a fixed number of
             * space group symmetry operators (and not a replicated unit cell).
             */
            if (splineX.length < nSymm) {
                logger.severe(" Programming Error: the number of reciprocal space symmetry operators changed.");
            }
            try {
                execute(0, nAtoms - 1, bSplineFillLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class BSplineFillLoop extends IntegerForLoop {

            private final double bSplineWork[][];
            private final IntegerSchedule schedule = IntegerSchedule.fixed();
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public BSplineFillLoop() {
                bSplineWork = new double[bSplineOrder][bSplineOrder];
            }

            @Override
            public IntegerSchedule schedule() {
                return schedule;
            }

            @Override
            public void start() {
                bSplineTime[getThreadIndex()] -= System.nanoTime();
            }

            @Override
            public void finish() {
                bSplineTime[getThreadIndex()] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                    final double x[] = coordinates[iSymOp][0];
                    final double y[] = coordinates[iSymOp][1];
                    final double z[] = coordinates[iSymOp][2];
                    final int initgridi[][] = initGrid[iSymOp];
                    final double splineXi[][][] = splineX[iSymOp];
                    final double splineYi[][][] = splineY[iSymOp];
                    final double splineZi[][][] = splineZ[iSymOp];
                    for (int i = lb; i <= ub; i++) {
                        if (initgridi[i] == null) {
                            initgridi[i] = new int[3];
                            splineXi[i] = new double[bSplineOrder][derivOrder + 1];
                            splineYi[i] = new double[bSplineOrder][derivOrder + 1];
                            splineZi[i] = new double[bSplineOrder][derivOrder + 1];
                        }
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final int[] grd = initgridi[i];
                        final double wx = xi * r00 + yi * r10 + zi * r20;
                        final double ux = wx - round(wx) + 0.5;
                        final double frx = fftX * ux;
                        final int ifrx = (int) frx;
                        final double bx = frx - ifrx;
                        grd[0] = ifrx - bSplineOrder;
                        bSplineDerivatives(bx, bSplineOrder, derivOrder, splineXi[i],
                                bSplineWork);
                        final double wy = xi * r01 + yi * r11 + zi * r21;
                        final double uy = wy - round(wy) + 0.5;
                        final double fry = fftY * uy;
                        final int ifry = (int) fry;
                        final double by = fry - ifry;
                        grd[1] = ifry - bSplineOrder;
                        bSplineDerivatives(by, bSplineOrder, derivOrder, splineYi[i],
                                bSplineWork);
                        final double wz = xi * r02 + yi * r12 + zi * r22;
                        final double uz = wz - round(wz) + 0.5;
                        final double frz = fftZ * uz;
                        final int ifrz = (int) frz;
                        final double bz = frz - ifrz;
                        grd[2] = ifrz - bSplineOrder;
                        bSplineDerivatives(bz, bSplineOrder, derivOrder, splineZi[i],
                                bSplineWork);
                    }
                }
            }
        }
    }

    private class SpatialPermanentLoop extends SpatialDensityLoop {

        private double globalMultipoles[][][] = null;
        private boolean use[] = null;
        private int threadIndex;
        private final BSplineRegion bSplines;

        public SpatialPermanentLoop(SpatialDensityRegion region, BSplineRegion splines) {
            super(region, region.nSymm, region.actualCount);
            this.bSplines = splines;
        }

        public void setPermanent(double globalMultipoles[][][]) {
            this.globalMultipoles = globalMultipoles;
        }

        private void setUse(boolean use[]) {
            this.use = use;
        }

        @Override
        public void start() {
            threadIndex = getThreadIndex();
            splinePermanentTime[threadIndex] -= System.nanoTime();
        }

        @Override
        public void finish() {
            splinePermanentTime[threadIndex] += System.nanoTime();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            splineCount[threadIndex]++;
            /**
             * Convert Cartesian multipoles in the global frame to fractional
             * multipoles.
             */
            final double gm[] = globalMultipoles[iSymm][n];
            final double fm[] = fracMultipole[iSymm][n];
            /**
             * Charge
             */
            fm[0] = gm[0];
            /**
             * Dipole
             */
            for (int j = 1; j < 4; j++) {
                fm[j] = 0.0;
                for (int k = 1; k < 4; k++) {
                    fm[j] = fm[j] + transformMultipoleMatrix[j][k] * gm[k];
                }
            }
            /**
             * Quadrupole
             */
            for (int j = 4; j < 10; j++) {
                fm[j] = 0.0;
                for (int k = 4; k < 7; k++) {
                    fm[j] = fm[j] + transformMultipoleMatrix[j][k] * gm[k];
                }
                for (int k = 7; k < 10; k++) {
                    fm[j] = fm[j] + transformMultipoleMatrix[j][k] * 2.0 * gm[k];
                }
                /**
                 * Fractional quadrupole components are pre-multiplied by a
                 * factor of 1/3 that arises in their potential.
                 */
                fm[j] = fm[j] / 3.0;
            }

            /**
             * Some atoms are not used during Lambda dynamics.
             */
            if (use != null && !use[n]) {
                return;
            }

            final double[][] splx = bSplines.splineX[iSymm][n];
            final double[][] sply = bSplines.splineY[iSymm][n];
            final double[][] splz = bSplines.splineZ[iSymm][n];
            final int igrd0 = bSplines.initGrid[iSymm][n][0];
            final int jgrd0 = bSplines.initGrid[iSymm][n][1];
            int k0 = bSplines.initGrid[iSymm][n][2];
            final double c = fm[t000];
            final double dx = fm[t100];
            final double dy = fm[t010];
            final double dz = fm[t001];
            final double qxx = fm[t200];
            final double qyy = fm[t020];
            final double qzz = fm[t002];
            final double qxy = fm[t110];
            final double qxz = fm[t101];
            final double qyz = fm[t011];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final double splzi[] = splz[ith3];
                final double v0 = splzi[0];
                final double v1 = splzi[1];
                final double v2 = splzi[2];
                final double c0 = c * v0;
                final double dx0 = dx * v0;
                final double dy0 = dy * v0;
                final double dz1 = dz * v1;
                final double qxx0 = qxx * v0;
                final double qyy0 = qyy * v0;
                final double qzz2 = qzz * v2;
                final double qxy0 = qxy * v0;
                final double qxz1 = qxz * v1;
                final double qyz1 = qyz * v1;
                final int k = mod(++k0, fftZ);
                int j0 = jgrd0;
                for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                    final double splyi[] = sply[ith2];
                    final double u0 = splyi[0];
                    final double u1 = splyi[1];
                    final double u2 = splyi[2];
                    final double term0 = (c0 + dz1 + qzz2) * u0 + (dy0 + qyz1) * u1 + qyy0 * u2;
                    final double term1 = (dx0 + qxz1) * u0 + qxy0 * u1;
                    final double term2 = qxx0 * u0;
                    final int j = mod(++j0, fftY);
                    int i0 = igrd0;
                    for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                        final int i = mod(++i0, fftX);
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        final double splxi[] = splx[ith1];
                        final double add = splxi[0] * term0 + splxi[1] * term1 + splxi[2] * term2;
                        final double current = splineBuffer.get(ii);
                        splineBuffer.put(ii, current + add);
                        /*
                         if (n == 0) {
                         logger.info(String.format(" %d %16.8f", ii, current + add));
                         } */
                        //splineGrid[ii] += add;
                    }
                }
            }
        }
    }

    private class SpatialInducedLoop extends SpatialDensityLoop {

        private double inducedDipole[][][] = null;
        private double inducedDipoleCR[][][] = null;
        private boolean use[] = null;
        private double a00, a01, a02;
        private double a10, a11, a12;
        private double a20, a21, a22;
        private final double m[][] = new double[3][3];
        private final BSplineRegion bSplineRegion;

        public SpatialInducedLoop(SpatialDensityRegion region, BSplineRegion splines) {
            super(region, region.nSymm, region.actualCount);
            this.bSplineRegion = splines;
        }

        public void setInducedDipoles(double inducedDipole[][][],
                double inducedDipoleCR[][][]) {
            this.inducedDipole = inducedDipole;
            this.inducedDipoleCR = inducedDipoleCR;
        }

        public void setUse(boolean use[]) {
            this.use = use;
        }

        @Override
        public void start() {
            splineInducedTime[getThreadIndex()] -= System.nanoTime();
            for (int i = 0; i < 3; i++) {
                m[0][i] = fftX * crystal.A[i][0];
                m[1][i] = fftY * crystal.A[i][1];
                m[2][i] = fftZ * crystal.A[i][2];
            }
            a00 = m[0][0];
            a01 = m[0][1];
            a02 = m[0][2];
            a10 = m[1][0];
            a11 = m[1][1];
            a12 = m[1][2];
            a20 = m[2][0];
            a21 = m[2][1];
            a22 = m[2][2];
        }

        @Override
        public void finish() {
            splineInducedTime[getThreadIndex()] += System.nanoTime();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            /**
             * Convert Cartesian induced dipole to fractional induced dipole.
             */
            double ind[] = inducedDipole[iSymm][n];
            double dx = ind[0];
            double dy = ind[1];
            double dz = ind[2];
            final double ux = a00 * dx + a01 * dy + a02 * dz;
            final double uy = a10 * dx + a11 * dy + a12 * dz;
            final double uz = a20 * dx + a21 * dy + a22 * dz;
            final double find[] = fracInducedDipole[iSymm][n];
            find[0] = ux;
            find[1] = uy;
            find[2] = uz;
            /**
             * Convert Cartesian induced dipole CR term to fractional induced
             * dipole CR term.
             */
            double indCR[] = inducedDipoleCR[iSymm][n];
            dx = indCR[0];
            dy = indCR[1];
            dz = indCR[2];
            final double px = a00 * dx + a01 * dy + a02 * dz;
            final double py = a10 * dx + a11 * dy + a12 * dz;
            final double pz = a20 * dx + a21 * dy + a22 * dz;
            final double findCR[] = fracInducedDipoleCR[iSymm][n];
            findCR[0] = px;
            findCR[1] = py;
            findCR[2] = pz;
            if (use != null && !use[n]) {
                return;
            }
            final double[][] splx = bSplineRegion.splineX[iSymm][n];
            final double[][] sply = bSplineRegion.splineY[iSymm][n];
            final double[][] splz = bSplineRegion.splineZ[iSymm][n];
            final int igrd0 = bSplineRegion.initGrid[iSymm][n][0];
            final int jgrd0 = bSplineRegion.initGrid[iSymm][n][1];
            int k0 = bSplineRegion.initGrid[iSymm][n][2];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final double splzi[] = splz[ith3];
                final double v0 = splzi[0];
                final double v1 = splzi[1];
                final double dx0 = ux * v0;
                final double dy0 = uy * v0;
                final double dz1 = uz * v1;
                final double px0 = px * v0;
                final double py0 = py * v0;
                final double pz1 = pz * v1;
                final int k = mod(++k0, fftZ);
                int j0 = jgrd0;
                for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                    final double splyi[] = sply[ith2];
                    final double u0 = splyi[0];
                    final double u1 = splyi[1];
                    final double term0 = dz1 * u0 + dy0 * u1;
                    final double term1 = dx0 * u0;
                    final double termp0 = pz1 * u0 + py0 * u1;
                    final double termp1 = px0 * u0;
                    final int j = mod(++j0, fftY);
                    int i0 = igrd0;
                    for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                        final int i = mod(++i0, fftX);
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        final double splxi[] = splx[ith1];
                        final double add = splxi[0] * term0 + splxi[1] * term1;
                        final double addi = splxi[0] * termp0 + splxi[1] * termp1;
                        final double current = splineBuffer.get(ii);
                        final double currenti = splineBuffer.get(ii + 1);
                        splineBuffer.put(ii, current + add);
                        splineBuffer.put(ii + 1, currenti + addi);
                        //splineGrid[ii] += add;
                        //splineGrid[ii + 1] += addi;
                    }
                }
            }
        }
    }

    private class SlicePermanentLoop extends SliceLoop {

        private double globalMultipoles[][][] = null;
        private final double[] fracMPole = new double[10];
        private boolean use[] = null;
        private int threadIndex;
        private final BSplineRegion bSplines;

        public SlicePermanentLoop(SliceRegion region, BSplineRegion splines) {
            super(region.nAtoms, region.nSymm, region);
            this.bSplines = splines;
        }

        public void setPermanent(double globalMultipoles[][][]) {
            this.globalMultipoles = globalMultipoles;
        }

        private void setUse(boolean use[]) {
            this.use = use;
        }

        @Override
        public void start() {
            threadIndex = getThreadIndex();
            splinePermanentTime[threadIndex] -= System.nanoTime();
            for (int i = 0; i < nSymm; i++) {
                gridAtomCount[i][threadIndex] = 0;
            }
        }

        @Override
        public void finish() {
            splinePermanentTime[threadIndex] += System.nanoTime();
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            boolean atomContributes = false;
            int k0 = bSplines.initGrid[iSymm][iAtom][2];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final int k = mod(++k0, fftZ);
                if (lb <= k && k <= ub) {
                    atomContributes = true;
                    break;
                }
            }
            if (!atomContributes) {
                return;
            }

            splineCount[threadIndex]++;

            // Add atom n to the list for the current symmOp and thread.
            int index = gridAtomCount[iSymm][threadIndex];
            gridAtomList[iSymm][threadIndex][index] = iAtom;
            gridAtomCount[iSymm][threadIndex]++;

            /**
             * Convert Cartesian multipoles in the global frame to fractional
             * multipoles.
             */
            final double gm[] = globalMultipoles[iSymm][iAtom];
            final double fm[] = fracMPole;
            /**
             * Charge
             */
            fm[0] = gm[0];
            /**
             * Dipole
             */
            for (int j = 1; j < 4; j++) {
                fm[j] = 0.0;
                for (int k = 1; k < 4; k++) {
                    fm[j] = fm[j] + transformMultipoleMatrix[j][k] * gm[k];
                }
            }
            /**
             * Quadrupole
             */
            for (int j = 4; j < 10; j++) {
                fm[j] = 0.0;
                for (int k = 4; k < 7; k++) {
                    fm[j] = fm[j] + transformMultipoleMatrix[j][k] * gm[k];
                }
                for (int k = 7; k < 10; k++) {
                    fm[j] = fm[j] + transformMultipoleMatrix[j][k] * 2.0 * gm[k];
                }
                /**
                 * Fractional quadrupole components are pre-multiplied by a
                 * factor of 1/3 that arises in their potential.
                 */
                fm[j] = fm[j] / 3.0;
            }

            System.arraycopy(fm, 0, fracMultipole[iSymm][iAtom], 0, 10);

            /**
             * Some atoms are not used during Lambda dynamics.
             */
            if (use != null && !use[iAtom]) {
                return;
            }

            final double[][] splx = bSplines.splineX[iSymm][iAtom];
            final double[][] sply = bSplines.splineY[iSymm][iAtom];
            final double[][] splz = bSplines.splineZ[iSymm][iAtom];
            final int igrd0 = bSplines.initGrid[iSymm][iAtom][0];
            final int jgrd0 = bSplines.initGrid[iSymm][iAtom][1];
            k0 = bSplines.initGrid[iSymm][iAtom][2];
            final double c = fm[t000];
            final double dx = fm[t100];
            final double dy = fm[t010];
            final double dz = fm[t001];
            final double qxx = fm[t200];
            final double qyy = fm[t020];
            final double qzz = fm[t002];
            final double qxy = fm[t110];
            final double qxz = fm[t101];
            final double qyz = fm[t011];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final double splzi[] = splz[ith3];
                final double v0 = splzi[0];
                final double v1 = splzi[1];
                final double v2 = splzi[2];
                final double c0 = c * v0;
                final double dx0 = dx * v0;
                final double dy0 = dy * v0;
                final double dz1 = dz * v1;
                final double qxx0 = qxx * v0;
                final double qyy0 = qyy * v0;
                final double qzz2 = qzz * v2;
                final double qxy0 = qxy * v0;
                final double qxz1 = qxz * v1;
                final double qyz1 = qyz * v1;
                final int k = mod(++k0, fftZ);
                if (k < lb || k > ub) {
                    continue;
                }
                int j0 = jgrd0;
                for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                    final double splyi[] = sply[ith2];
                    final double u0 = splyi[0];
                    final double u1 = splyi[1];
                    final double u2 = splyi[2];
                    // Pieces of a multipole
                    final double term0 = (c0 + dz1 + qzz2) * u0 + (dy0 + qyz1) * u1 + qyy0 * u2;
                    final double term1 = (dx0 + qxz1) * u0 + qxy0 * u1;
                    final double term2 = qxx0 * u0;
                    final int j = mod(++j0, fftY);
                    int i0 = igrd0;
                    for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                        final int i = mod(++i0, fftX);
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        final double splxi[] = splx[ith1];
                        final double add = splxi[0] * term0 + splxi[1] * term1 + splxi[2] * term2;
                        final double current = splineBuffer.get(ii);
                        splineBuffer.put(ii, current + add);
                        /*
                         if (n == 0) {
                         logger.info(String.format(" %d %16.8f", ii, current + add));
                         } */
                        //splineGrid[ii] += add;
                    }
                }
            }
        }
    }

    private class SliceInducedLoop extends SliceLoop {

        private double inducedDipole[][][] = null;
        private double inducedDipoleCR[][][] = null;
        private final double m[][] = new double[3][3];
        private boolean use[] = null;
        private final BSplineRegion bSplineRegion;
        private double a00, a01, a02;
        private double a10, a11, a12;
        private double a20, a21, a22;

        public SliceInducedLoop(SliceRegion region, BSplineRegion splines) {
            super(region.nAtoms, region.nSymm, region);
            this.bSplineRegion = splines;
        }

        public void setInducedDipoles(double inducedDipole[][][],
                double inducedDipoleCR[][][]) {
            this.inducedDipole = inducedDipole;
            this.inducedDipoleCR = inducedDipoleCR;
        }

        public void setUse(boolean use[]) {
            this.use = use;
        }

        @Override
        public void start() {
            int threadIndex = getThreadIndex();
            splineInducedTime[threadIndex] -= System.nanoTime();
            for (int i = 0; i < 3; i++) {
                m[0][i] = fftX * crystal.A[i][0];
                m[1][i] = fftY * crystal.A[i][1];
                m[2][i] = fftZ * crystal.A[i][2];
            }
            a00 = m[0][0];
            a01 = m[0][1];
            a02 = m[0][2];
            a10 = m[1][0];
            a11 = m[1][1];
            a12 = m[1][2];
            a20 = m[2][0];
            a21 = m[2][1];
            a22 = m[2][2];
        }

        @Override
        public void finish() {
            int threadIndex = getThreadIndex();
            splineInducedTime[threadIndex] += System.nanoTime();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            int ti = getThreadIndex();
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                int list[] = gridAtomList[iSymm][ti];
                int n = gridAtomCount[iSymm][ti];
                for (int i = 0; i < n; i++) {
                    int iAtom = list[i];
                    gridDensity(iSymm, iAtom, lb, ub);
                }
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            /**
             * Convert Cartesian induced dipole to fractional induced dipole.
             */
            double ind[] = inducedDipole[iSymm][iAtom];
            double dx = ind[0];
            double dy = ind[1];
            double dz = ind[2];
            final double ux = a00 * dx + a01 * dy + a02 * dz;
            final double uy = a10 * dx + a11 * dy + a12 * dz;
            final double uz = a20 * dx + a21 * dy + a22 * dz;
            final double find[] = fracInducedDipole[iSymm][iAtom];
            find[0] = ux;
            find[1] = uy;
            find[2] = uz;
            /**
             * Convert Cartesian induced dipole CR term to fractional induced
             * dipole CR term.
             */
            double indCR[] = inducedDipoleCR[iSymm][iAtom];
            dx = indCR[0];
            dy = indCR[1];
            dz = indCR[2];
            final double px = a00 * dx + a01 * dy + a02 * dz;
            final double py = a10 * dx + a11 * dy + a12 * dz;
            final double pz = a20 * dx + a21 * dy + a22 * dz;
            final double findCR[] = fracInducedDipoleCR[iSymm][iAtom];
            findCR[0] = px;
            findCR[1] = py;
            findCR[2] = pz;
            if (use != null && !use[iAtom]) {
                return;
            }
            final double[][] splx = bSplineRegion.splineX[iSymm][iAtom];
            final double[][] sply = bSplineRegion.splineY[iSymm][iAtom];
            final double[][] splz = bSplineRegion.splineZ[iSymm][iAtom];
            final int igrd0 = bSplineRegion.initGrid[iSymm][iAtom][0];
            final int jgrd0 = bSplineRegion.initGrid[iSymm][iAtom][1];
            int k0 = bSplineRegion.initGrid[iSymm][iAtom][2];
            for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                final double splzi[] = splz[ith3];
                final double v0 = splzi[0];
                final double v1 = splzi[1];
                final double dx0 = ux * v0;
                final double dy0 = uy * v0;
                final double dz1 = uz * v1;
                final double px0 = px * v0;
                final double py0 = py * v0;
                final double pz1 = pz * v1;
                final int k = mod(++k0, fftZ);
                if (k < lb || k > ub) {
                    continue;
                }
                int j0 = jgrd0;
                for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                    final double splyi[] = sply[ith2];
                    final double u0 = splyi[0];
                    final double u1 = splyi[1];
                    final double term0 = dz1 * u0 + dy0 * u1;
                    final double term1 = dx0 * u0;
                    final double termp0 = pz1 * u0 + py0 * u1;
                    final double termp1 = px0 * u0;
                    final int j = mod(++j0, fftY);
                    int i0 = igrd0;
                    for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                        final int i = mod(++i0, fftX);
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        final double splxi[] = splx[ith1];
                        final double add = splxi[0] * term0 + splxi[1] * term1;
                        final double addi = splxi[0] * termp0 + splxi[1] * termp1;
                        final double current = splineBuffer.get(ii);
                        final double currenti = splineBuffer.get(ii + 1);
                        splineBuffer.put(ii, current + add);
                        splineBuffer.put(ii + 1, currenti + addi);
                        //splineGrid[ii] += add;
                        //splineGrid[ii + 1] += addi;
                    }
                }
            }
        }
    }

    private class PermanentPhiRegion extends ParallelRegion {

        private final FractionalPhiLoop fractionalPhiLoop[];
        private final BSplineRegion bSplineRegion;
        private double cartPermanentPhi[][];

        public PermanentPhiRegion(BSplineRegion bSplineRegion) {
            this.bSplineRegion = bSplineRegion;
            fractionalPhiLoop = new FractionalPhiLoop[threadCount];
        }

        public void setCartPermanentPhi(double cartPermanentPhi[][]) {
            this.cartPermanentPhi = cartPermanentPhi;
        }

        @Override
        public void run() {

            int threadIndex = getThreadIndex();

            if (fractionalPhiLoop[threadIndex] == null) {
                fractionalPhiLoop[threadIndex] = new FractionalPhiLoop();
            }

            try {
                execute(0, nAtoms - 1, fractionalPhiLoop[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class FractionalPhiLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                int threadIndex = getThreadIndex();
                permanentPhiTime[threadIndex] -= System.nanoTime();
            }

            @Override
            public void finish() {
                int threadIndex = getThreadIndex();
                permanentPhiTime[threadIndex] += System.nanoTime();
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    final double[][] splx = bSplineRegion.splineX[0][n];
                    final double[][] sply = bSplineRegion.splineY[0][n];
                    final double[][] splz = bSplineRegion.splineZ[0][n];
                    final int igrd[] = bSplineRegion.initGrid[0][n];
                    final int igrd0 = igrd[0];
                    final int jgrd0 = igrd[1];
                    int k0 = igrd[2];
                    double tuv000 = 0.0;
                    double tuv100 = 0.0;
                    double tuv010 = 0.0;
                    double tuv001 = 0.0;
                    double tuv200 = 0.0;
                    double tuv020 = 0.0;
                    double tuv002 = 0.0;
                    double tuv110 = 0.0;
                    double tuv101 = 0.0;
                    double tuv011 = 0.0;
                    double tuv300 = 0.0;
                    double tuv030 = 0.0;
                    double tuv003 = 0.0;
                    double tuv210 = 0.0;
                    double tuv201 = 0.0;
                    double tuv120 = 0.0;
                    double tuv021 = 0.0;
                    double tuv102 = 0.0;
                    double tuv012 = 0.0;
                    double tuv111 = 0.0;
                    for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                        final int k = mod(++k0, fftZ);
                        int j0 = jgrd0;
                        double tu00 = 0.0;
                        double tu10 = 0.0;
                        double tu01 = 0.0;
                        double tu20 = 0.0;
                        double tu11 = 0.0;
                        double tu02 = 0.0;
                        double tu30 = 0.0;
                        double tu21 = 0.0;
                        double tu12 = 0.0;
                        double tu03 = 0.0;
                        for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                            final int j = mod(++j0, fftY);
                            int i0 = igrd0;
                            double t0 = 0.0;
                            double t1 = 0.0;
                            double t2 = 0.0;
                            double t3 = 0.0;
                            for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                                final int i = mod(++i0, fftX);
                                final int ii = iComplex3D(i, j, k, fftX, fftY);
                                //final double tq = splineGrid[ii];
                                final double tq = splineBuffer.get(ii);
                                final double splxi[] = splx[ith1];
                                t0 += tq * splxi[0];
                                t1 += tq * splxi[1];
                                t2 += tq * splxi[2];
                                t3 += tq * splxi[3];
                            }
                            final double splyi[] = sply[ith2];
                            final double u0 = splyi[0];
                            final double u1 = splyi[1];
                            final double u2 = splyi[2];
                            final double u3 = splyi[3];
                            tu00 += t0 * u0;
                            tu10 += t1 * u0;
                            tu01 += t0 * u1;
                            tu20 += t2 * u0;
                            tu11 += t1 * u1;
                            tu02 += t0 * u2;
                            tu30 += t3 * u0;
                            tu21 += t2 * u1;
                            tu12 += t1 * u2;
                            tu03 += t0 * u3;
                        }
                        final double splzi[] = splz[ith3];
                        final double v0 = splzi[0];
                        final double v1 = splzi[1];
                        final double v2 = splzi[2];
                        final double v3 = splzi[3];
                        tuv000 += tu00 * v0;
                        tuv100 += tu10 * v0;
                        tuv010 += tu01 * v0;
                        tuv001 += tu00 * v1;
                        tuv200 += tu20 * v0;
                        tuv020 += tu02 * v0;
                        tuv002 += tu00 * v2;
                        tuv110 += tu11 * v0;
                        tuv101 += tu10 * v1;
                        tuv011 += tu01 * v1;
                        tuv300 += tu30 * v0;
                        tuv030 += tu03 * v0;
                        tuv003 += tu00 * v3;
                        tuv210 += tu21 * v0;
                        tuv201 += tu20 * v1;
                        tuv120 += tu12 * v0;
                        tuv021 += tu02 * v1;
                        tuv102 += tu10 * v2;
                        tuv012 += tu01 * v2;
                        tuv111 += tu11 * v1;
                    }
                    double out[] = fracMultipolePhi[n];
                    out[t000] = tuv000;
                    out[t100] = tuv100;
                    out[t010] = tuv010;
                    out[t001] = tuv001;
                    out[t200] = tuv200;
                    out[t020] = tuv020;
                    out[t002] = tuv002;
                    out[t110] = tuv110;
                    out[t101] = tuv101;
                    out[t011] = tuv011;
                    out[t300] = tuv300;
                    out[t030] = tuv030;
                    out[t003] = tuv003;
                    out[t210] = tuv210;
                    out[t201] = tuv201;
                    out[t120] = tuv120;
                    out[t021] = tuv021;
                    out[t102] = tuv102;
                    out[t012] = tuv012;
                    out[t111] = tuv111;
                    double in[] = out;
                    out = cartPermanentPhi[n];
                    out[0] = transformFieldMatrix[0][0] * in[0];
                    for (int j = 1; j < 4; j++) {
                        out[j] = 0.0;
                        for (int k = 1; k < 4; k++) {
                            out[j] += transformFieldMatrix[j][k] * in[k];
                        }
                    }
                    for (int j = 4; j < 10; j++) {
                        out[j] = 0.0;
                        for (int k = 4; k < 10; k++) {
                            out[j] += transformFieldMatrix[j][k] * in[k];
                        }
                    }
                }
            }
        }
    }

    private class InducedPhiRegion extends ParallelRegion {

        private final BSplineRegion bSplineRegion;
        private final PolarizationPhiInducedLoop polarizationPhiInducedLoops[];
        private double cartInducedDipolePhi[][];
        private double cartInducedDipoleCRPhi[][];

        public InducedPhiRegion(BSplineRegion bSplineRegion) {
            this.bSplineRegion = bSplineRegion;
            polarizationPhiInducedLoops = new PolarizationPhiInducedLoop[threadCount];
        }

        public void setCartInducedDipolePhi(double cartInducedDipolePhi[][],
                double cartInducedDipoleCRPhi[][]) {
            this.cartInducedDipolePhi = cartInducedDipolePhi;
            this.cartInducedDipoleCRPhi = cartInducedDipoleCRPhi;
        }

        @Override
        public void run() {

            int threadIndex = getThreadIndex();

            if (polarizationPhiInducedLoops[threadIndex] == null) {
                polarizationPhiInducedLoops[threadIndex] = new PolarizationPhiInducedLoop();
            }

            try {
                execute(0, nAtoms - 1, polarizationPhiInducedLoops[threadIndex]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }

        }

        private class PolarizationPhiInducedLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void start() {
                int threadIndex = getThreadIndex();
                inducedPhiTime[threadIndex] -= System.nanoTime();
            }

            @Override
            public void finish() {
                int threadIndex = getThreadIndex();
                inducedPhiTime[threadIndex] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                for (int n = lb; n <= ub; n++) {
                    final double[][] splx = bSplineRegion.splineX[0][n];
                    final double[][] sply = bSplineRegion.splineY[0][n];
                    final double[][] splz = bSplineRegion.splineZ[0][n];
                    final int igrd[] = bSplineRegion.initGrid[0][n];
                    final int igrd0 = igrd[0];
                    final int jgrd0 = igrd[1];
                    int k0 = igrd[2];
                    double tuv000 = 0.0;
                    double tuv100 = 0.0;
                    double tuv010 = 0.0;
                    double tuv001 = 0.0;
                    double tuv200 = 0.0;
                    double tuv020 = 0.0;
                    double tuv002 = 0.0;
                    double tuv110 = 0.0;
                    double tuv101 = 0.0;
                    double tuv011 = 0.0;
                    double tuv300 = 0.0;
                    double tuv030 = 0.0;
                    double tuv003 = 0.0;
                    double tuv210 = 0.0;
                    double tuv201 = 0.0;
                    double tuv120 = 0.0;
                    double tuv021 = 0.0;
                    double tuv102 = 0.0;
                    double tuv012 = 0.0;
                    double tuv111 = 0.0;
                    double tuv000p = 0.0;
                    double tuv100p = 0.0;
                    double tuv010p = 0.0;
                    double tuv001p = 0.0;
                    double tuv200p = 0.0;
                    double tuv020p = 0.0;
                    double tuv002p = 0.0;
                    double tuv110p = 0.0;
                    double tuv101p = 0.0;
                    double tuv011p = 0.0;
                    double tuv300p = 0.0;
                    double tuv030p = 0.0;
                    double tuv003p = 0.0;
                    double tuv210p = 0.0;
                    double tuv201p = 0.0;
                    double tuv120p = 0.0;
                    double tuv021p = 0.0;
                    double tuv102p = 0.0;
                    double tuv012p = 0.0;
                    double tuv111p = 0.0;
                    for (int ith3 = 0; ith3 < bSplineOrder; ith3++) {
                        final int k = mod(++k0, fftZ);
                        int j0 = jgrd0;
                        double tu00 = 0.0;
                        double tu10 = 0.0;
                        double tu01 = 0.0;
                        double tu20 = 0.0;
                        double tu11 = 0.0;
                        double tu02 = 0.0;
                        double tu30 = 0.0;
                        double tu21 = 0.0;
                        double tu12 = 0.0;
                        double tu03 = 0.0;
                        double tu00p = 0.0;
                        double tu10p = 0.0;
                        double tu01p = 0.0;
                        double tu20p = 0.0;
                        double tu11p = 0.0;
                        double tu02p = 0.0;
                        double tu30p = 0.0;
                        double tu21p = 0.0;
                        double tu12p = 0.0;
                        double tu03p = 0.0;
                        for (int ith2 = 0; ith2 < bSplineOrder; ith2++) {
                            final int j = mod(++j0, fftY);
                            int i0 = igrd0;
                            double t0 = 0.0;
                            double t1 = 0.0;
                            double t2 = 0.0;
                            double t3 = 0.0;
                            double t0p = 0.0;
                            double t1p = 0.0;
                            double t2p = 0.0;
                            double t3p = 0.0;
                            for (int ith1 = 0; ith1 < bSplineOrder; ith1++) {
                                final int i = mod(++i0, fftX);
                                final int ii = iComplex3D(i, j, k, fftX, fftY);
                                //final double tq = splineGrid[ii];
                                //final double tp = splineGrid[ii + 1];
                                final double tq = splineBuffer.get(ii);
                                final double tp = splineBuffer.get(ii + 1);
                                final double splxi[] = splx[ith1];
                                t0 += tq * splxi[0];
                                t1 += tq * splxi[1];
                                t2 += tq * splxi[2];
                                t3 += tq * splxi[3];
                                t0p += tp * splxi[0];
                                t1p += tp * splxi[1];
                                t2p += tp * splxi[2];
                                t3p += tp * splxi[3];
                            }
                            final double splyi[] = sply[ith2];
                            final double u0 = splyi[0];
                            final double u1 = splyi[1];
                            final double u2 = splyi[2];
                            final double u3 = splyi[3];
                            tu00 += t0 * u0;
                            tu10 += t1 * u0;
                            tu01 += t0 * u1;
                            tu20 += t2 * u0;
                            tu11 += t1 * u1;
                            tu02 += t0 * u2;
                            tu30 += t3 * u0;
                            tu21 += t2 * u1;
                            tu12 += t1 * u2;
                            tu03 += t0 * u3;
                            tu00p += t0p * u0;
                            tu10p += t1p * u0;
                            tu01p += t0p * u1;
                            tu20p += t2p * u0;
                            tu11p += t1p * u1;
                            tu02p += t0p * u2;
                            tu30p += t3p * u0;
                            tu21p += t2p * u1;
                            tu12p += t1p * u2;
                            tu03p += t0p * u3;
                        }
                        final double splzi[] = splz[ith3];
                        final double v0 = splzi[0];
                        final double v1 = splzi[1];
                        final double v2 = splzi[2];
                        final double v3 = splzi[3];
                        tuv000 += tu00 * v0;
                        tuv100 += tu10 * v0;
                        tuv010 += tu01 * v0;
                        tuv001 += tu00 * v1;
                        tuv200 += tu20 * v0;
                        tuv020 += tu02 * v0;
                        tuv002 += tu00 * v2;
                        tuv110 += tu11 * v0;
                        tuv101 += tu10 * v1;
                        tuv011 += tu01 * v1;
                        tuv300 += tu30 * v0;
                        tuv030 += tu03 * v0;
                        tuv003 += tu00 * v3;
                        tuv210 += tu21 * v0;
                        tuv201 += tu20 * v1;
                        tuv120 += tu12 * v0;
                        tuv021 += tu02 * v1;
                        tuv102 += tu10 * v2;
                        tuv012 += tu01 * v2;
                        tuv111 += tu11 * v1;
                        tuv000p += tu00p * v0;
                        tuv100p += tu10p * v0;
                        tuv010p += tu01p * v0;
                        tuv001p += tu00p * v1;
                        tuv200p += tu20p * v0;
                        tuv020p += tu02p * v0;
                        tuv002p += tu00p * v2;
                        tuv110p += tu11p * v0;
                        tuv101p += tu10p * v1;
                        tuv011p += tu01p * v1;
                        tuv300p += tu30p * v0;
                        tuv030p += tu03p * v0;
                        tuv003p += tu00p * v3;
                        tuv210p += tu21p * v0;
                        tuv201p += tu20p * v1;
                        tuv120p += tu12p * v0;
                        tuv021p += tu02p * v1;
                        tuv102p += tu10p * v2;
                        tuv012p += tu01p * v2;
                        tuv111p += tu11p * v1;
                    }
                    double out[] = fracInducedDipolePhi[n];
                    out[t000] = tuv000;
                    out[t100] = tuv100;
                    out[t010] = tuv010;
                    out[t001] = tuv001;
                    out[t200] = tuv200;
                    out[t020] = tuv020;
                    out[t002] = tuv002;
                    out[t110] = tuv110;
                    out[t101] = tuv101;
                    out[t011] = tuv011;
                    out[t300] = tuv300;
                    out[t030] = tuv030;
                    out[t003] = tuv003;
                    out[t210] = tuv210;
                    out[t201] = tuv201;
                    out[t120] = tuv120;
                    out[t021] = tuv021;
                    out[t102] = tuv102;
                    out[t012] = tuv012;
                    out[t111] = tuv111;
                    double in[] = out;
                    out = cartInducedDipolePhi[n];
                    out[0] = transformFieldMatrix[0][0] * in[0];
                    for (int j = 1; j < 4; j++) {
                        out[j] = 0.0;
                        for (int k = 1; k < 4; k++) {
                            out[j] += transformFieldMatrix[j][k] * in[k];
                        }
                    }
                    for (int j = 4; j < 10; j++) {
                        out[j] = 0.0;
                        for (int k = 4; k < 10; k++) {
                            out[j] += transformFieldMatrix[j][k] * in[k];
                        }
                    }

                    out = fracInducedDipolePhiCR[n];
                    out[t000] = tuv000p;
                    out[t100] = tuv100p;
                    out[t010] = tuv010p;
                    out[t001] = tuv001p;
                    out[t200] = tuv200p;
                    out[t020] = tuv020p;
                    out[t002] = tuv002p;
                    out[t110] = tuv110p;
                    out[t101] = tuv101p;
                    out[t011] = tuv011p;
                    out[t300] = tuv300p;
                    out[t030] = tuv030p;
                    out[t003] = tuv003p;
                    out[t210] = tuv210p;
                    out[t201] = tuv201p;
                    out[t120] = tuv120p;
                    out[t021] = tuv021p;
                    out[t102] = tuv102p;
                    out[t012] = tuv012p;
                    out[t111] = tuv111p;
                    in = out;
                    out = cartInducedDipoleCRPhi[n];
                    out[0] = transformFieldMatrix[0][0] * in[0];
                    for (int j = 1; j < 4; j++) {
                        out[j] = 0.0;
                        for (int k = 1; k < 4; k++) {
                            out[j] += transformFieldMatrix[j][k] * in[k];
                        }
                    }
                    for (int j = 4; j < 10; j++) {
                        out[j] = 0.0;
                        for (int k = 4; k < 10; k++) {
                            out[j] += transformFieldMatrix[j][k] * in[k];
                        }
                    }
                }
            }
        }
    }

    private double[] generalizedInfluenceFunction() {

        double influenceFunction[] = new double[fftSpace / 2];

        double bsModX[] = new double[fftX];
        double bsModY[] = new double[fftY];
        double bsModZ[] = new double[fftZ];
        int maxfft = max(max(fftX, fftY), fftZ);
        double bsArray[] = new double[maxfft];
        double c[] = new double[bSplineOrder];

        bSpline(0.0, bSplineOrder, c);
        for (int i = 1; i < bSplineOrder + 1; i++) {
            bsArray[i] = c[i - 1];
        }
        discreteFTMod(bsModX, bsArray, fftX, bSplineOrder);
        discreteFTMod(bsModY, bsArray, fftY, bSplineOrder);
        discreteFTMod(bsModZ, bsArray, fftZ, bSplineOrder);

        double r00 = crystal.A[0][0];
        double r01 = crystal.A[0][1];
        double r02 = crystal.A[0][2];
        double r10 = crystal.A[1][0];
        double r11 = crystal.A[1][1];
        double r12 = crystal.A[1][2];
        double r20 = crystal.A[2][0];
        double r21 = crystal.A[2][1];
        double r22 = crystal.A[2][2];
        int ntot = fftX * fftY * fftZ;
        double piTerm = (PI / aEwald) * (PI / aEwald);
        double volTerm = PI * crystal.volume;
        int nfXY = fftX * fftY;
        int nX_2 = (fftX + 1) / 2;
        int nY_2 = (fftY + 1) / 2;
        int nZ_2 = (fftZ + 1) / 2;

        for (int i = 0; i < ntot - 1; i++) {
            int kZ = (i + 1) / nfXY;
            int j = i - kZ * nfXY + 1;
            int kY = j / fftX;
            int kX = j - kY * fftX;
            int h = kX;
            int k = kY;
            int l = kZ;
            if (kX >= nX_2) {
                h -= fftX;
            }
            if (kY >= nY_2) {
                k -= fftY;
            }
            if (kZ >= nZ_2) {
                l -= fftZ;
            }
            double sX = r00 * h + r01 * k + r02 * l;
            double sY = r10 * h + r11 * k + r12 * l;
            double sZ = r20 * h + r21 * k + r22 * l;
            double sSquared = sX * sX + sY * sY + sZ * sZ;
            double term = -piTerm * sSquared;
            double expterm = 0.0;
            if (term > -50.0) {
                double denom = sSquared * volTerm * bsModX[kX] * bsModY[kY] * bsModZ[kZ];
                expterm = exp(term) / denom;
                if (crystal.aperiodic()) {
                    expterm *= (1.0 - cos(PI * crystal.a * sqrt(sSquared)));
                }
            }
            int ii = iComplex3D(kX, kY, kZ, fftX, fftY) / 2;
            influenceFunction[ii] = expterm;
        }
        /**
         * Account for the zeroth grid point for a periodic system.
         */
        influenceFunction[0] = 0.0;
        if (crystal.aperiodic()) {
            influenceFunction[0] = 0.5 * PI / crystal.a;
        }

        return influenceFunction;
    }

    private void transformMultipoleMatrix() {
        double a[][] = new double[3][3];
        for (int i = 0; i < 3; i++) {
            a[0][i] = fftX * crystal.A[i][0];
            a[1][i] = fftY * crystal.A[i][1];
            a[2][i] = fftZ * crystal.A[i][2];
        }

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                transformMultipoleMatrix[i][j] = 0.0;
            }
        }

        // Charge
        transformMultipoleMatrix[0][0] = 1.0;
        // Dipole
        for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
                transformMultipoleMatrix[i][j] = a[i - 1][j - 1];
            }
        }
        // Quadrupole
        for (int i1 = 0; i1 < 3; i1++) {
            int k = qi1[i1];
            for (int i2 = 0; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                transformMultipoleMatrix[i1 + 4][i2 + 4] = a[k][i] * a[k][j];
            }
        }
        for (int i1 = 3; i1 < 6; i1++) {
            int k = qi1[i1];
            int l = qi2[i1];
            for (int i2 = 0; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                transformMultipoleMatrix[i1 + 4][i2 + 4] = a[k][i] * a[l][j] + a[k][j] * a[l][i];
            }
        }
    }

    private void transformFieldMatrix() {
        double a[][] = new double[3][3];

        for (int i = 0; i < 3; i++) {
            a[i][0] = fftX * crystal.A[i][0];
            a[i][1] = fftY * crystal.A[i][1];
            a[i][2] = fftZ * crystal.A[i][2];
        }
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                transformFieldMatrix[i][j] = 0.0;
            }
        }
        transformFieldMatrix[0][0] = 1.0;
        for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
                transformFieldMatrix[i][j] = a[i - 1][j - 1];
            }
        }
        for (int i1 = 0; i1 < 3; i1++) {
            int k = qi1[i1];
            for (int i2 = 0; i2 < 3; i2++) {
                int i = qi1[i2];
                transformFieldMatrix[i1 + 4][i2 + 4] = a[k][i] * a[k][i];
            }
            for (int i2 = 3; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                transformFieldMatrix[i1 + 4][i2 + 4] = 2.0 * a[k][i] * a[k][j];
            }
        }
        for (int i1 = 3; i1 < 6; i1++) {
            int k = qi1[i1];
            int n = qi2[i1];
            for (int i2 = 0; i2 < 3; i2++) {
                int i = qi1[i2];
                transformFieldMatrix[i1 + 4][i2 + 4] = a[k][i] * a[n][i];
            }
            for (int i2 = 3; i2 < 6; i2++) {
                int i = qi1[i2];
                int j = qi2[i2];
                transformFieldMatrix[i1 + 4][i2 + 4] = a[k][i] * a[n][j] + a[n][i] * a[k][j];
            }
        }
    }

    /**
     * Computes the modulus of the discrete Fourier Transform of "bsarray" and
     * stores it in "bsmod".
     *
     * @param bsmod
     * @param bsarray
     * @param nfft
     * @param order
     */
    private static void discreteFTMod(double bsmod[], double bsarray[],
            int nfft, int order) {
        /**
         * Get the modulus of the discrete Fourier fft.
         */
        double factor = 2.0 * PI / nfft;
        for (int i = 0; i < nfft; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int j = 0; j < nfft; j++) {
                double arg = factor * (i * j);
                sum1 = sum1 + bsarray[j] * cos(arg);
                sum2 = sum2 + bsarray[j] * sin(arg);
            }
            bsmod[i] = sum1 * sum1 + sum2 * sum2;
        }
        /**
         * Fix for exponential Euler spline interpolation failure.
         */
        double eps = 1.0e-7;
        if (bsmod[0] < eps) {
            bsmod[0] = 0.5 * bsmod[1];
        }
        for (int i = 1; i < nfft - 1; i++) {
            if (bsmod[i] < eps) {
                bsmod[i] = 0.5 * (bsmod[i - 1] + bsmod[i + 1]);
            }
        }
        if (bsmod[nfft - 1] < eps) {
            bsmod[nfft - 1] = 0.5 * bsmod[nfft - 2];
        }
        /**
         * Compute and apply the optimal zeta coefficient.
         */
        int jcut = 50;
        int order2 = 2 * order;
        for (int i = 0; i < nfft; i++) {
            int k = i;
            double zeta;
            if (i > nfft / 2) {
                k = k - nfft;
            }
            if (k == 0) {
                zeta = 1.0;
            } else {
                double sum1 = 1.0;
                double sum2 = 1.0;
                factor = PI * k / nfft;
                for (int j = 0; j < jcut; j++) {
                    double arg = factor / (factor + PI * (j + 1));
                    sum1 = sum1 + pow(arg, order);
                    sum2 = sum2 + pow(arg, order2);
                }
                for (int j = 0; j < jcut; j++) {
                    double arg = factor / (factor - PI * (j + 1));
                    sum1 = sum1 + pow(arg, order);
                    sum2 = sum2 + pow(arg, order2);
                }
                zeta = sum2 / sum1;
            }
            bsmod[i] = bsmod[i] * zeta * zeta;
        }
    }
    /**
     * First lookup index to pack a 2D tensor into a 1D array.
     */
    private static final int qi1[] = {0, 1, 2, 0, 0, 1};
    /**
     * Second lookup index to pack a 2D tensor into a 1D array.
     */
    private static final int qi2[] = {0, 1, 2, 1, 2, 2};
    private final double transformFieldMatrix[][] = new double[10][10];
    private final double transformMultipoleMatrix[][] = new double[10][10];
    private final double a[][] = new double[3][3];
    private static final int tensorCount = TensorRecursion.tensorCount(3);
    private static final double toSeconds = 1.0e-9;
}
