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
package ffx.xray;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static org.apache.commons.math.util.FastMath.exp;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SymOp;
import ffx.numerics.ComplexNumber;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Complex3DParallel;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.SliceLoop;
import ffx.potential.nonbonded.SliceRegion;
import ffx.potential.nonbonded.SpatialDensityLoop;
import ffx.potential.nonbonded.SpatialDensityRegion;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.numerics.fft.Complex3D.iComplex3D;
import static ffx.xray.CrystalReciprocalSpace.SolventModel.BINARY;
import static ffx.xray.CrystalReciprocalSpace.SolventModel.GAUSSIAN;
import static ffx.xray.CrystalReciprocalSpace.SolventModel.POLYNOMIAL;

/**
 * Structure factor calculation (including bulk solvent structure factors)
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 *
 * @author Timothy D. Fenn
 */
public class CrystalReciprocalSpace {

    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());

    /**
     * The possible solvent model methods
     */
    public static interface SolventModel {

        /**
         * do not model solvent scattering
         */
        public static final int NONE = 1;
        /**
         * the classical binary (0, 1) model
         */
        public static final int BINARY = 2;
        /**
         * a Gaussian switch model
         */
        public static final int GAUSSIAN = 3;
        /**
         * a cubic polynomial switch model (default)
         */
        public static final int POLYNOMIAL = 4;
    }

    public enum GridMethod {

        SPATIAL, SLICE
    }

    private static double toSeconds = 0.000000001;
    private boolean solvent = false;
    private final boolean neutron;
    private boolean useThreeGaussians = true;
    protected boolean lambdaTerm = false;
    protected double lambda = 1.0;
    // not final for purposes of finite differences
    private double coordinates[][][];
    private final double bAdd;
    private double aRad;
    private double weight = 1.0;
    protected double solventA;
    protected double solventB;
    private final double fftScale;
    protected final double densityGrid[];
    protected final double solventGrid[];
    private int aRadGrid;
    protected int solventModel;
    private final int nSymm;
    private final int bulkNSymm;
    private final int nAtoms;
    private final int fftX, fftY, fftZ;
    private final int halfFFTX, halfFFTY, halfFFTZ;
    private final int complexFFT3DSpace;
    private final int threadCount;
    private final ParallelTeam parallelTeam;
    private final Atom atoms[];
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionList;
    private final FormFactor atomFormFactors[][];
    private final FormFactor solventFormFactors[][];
    private final GridMethod gridMethod;

    private final SpatialDensityRegion atomicDensityRegion;
    private final AtomicDensityLoop atomicDensityLoops[];
    private final SpatialDensityRegion solventDensityRegion;
    private final SolventDensityLoop solventDensityLoops[];
    private final SolventDensityLoop bulkSolventDensityLoops[];

    private final SliceRegion atomicSliceRegion;
    private final AtomicSliceLoop atomicSliceLoops[];
    private final SliceRegion solventSliceRegion;
    private final SolventSliceLoop solventSliceLoops[];
    private final SolventSliceLoop bulkSolventSliceLoops[];

    private final InitRegion initRegion;
    private final ExtractRegion extractRegion;
    private final AtomicScaleRegion atomicScaleRegion;
    private final AtomicGradientRegion atomicGradientRegion;

    private final BabinetRegion babinetRegion;
    private final SolventScaleRegion solventScaleRegion;
    private final SolventGradientRegion solventGradientRegion;

    private final ParallelTeam fftTeam;
    private final Complex3DParallel complexFFT3D;

    /**
     * Crystal Reciprocal Space constructor, assumes this is not a bulk solvent
     * mask and is not a neutron data set
     *
     * @param reflectionList the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     */
    public CrystalReciprocalSpace(ReflectionList reflectionList,
            Atom atoms[], ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this(reflectionList, atoms, fftTeam, parallelTeam, false, false,
                SolventModel.POLYNOMIAL);
    }

    /**
     * Crystal Reciprocal Space constructor, assumes this is not a neutron data
     * set and implements a polynomial bulk solvent mask if needed
     *
     * @param reflectionList the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventMask true if this is a bulk solvent mask
     */
    public CrystalReciprocalSpace(ReflectionList reflectionList,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventMask) {
        this(reflectionList, atoms, fftTeam, parallelTeam, solventMask, false,
                SolventModel.POLYNOMIAL);
    }

    /**
     * Crystal Reciprocal Space constructor, assumes a polynomial bulk solvent
     * mask if needed
     *
     * @param reflectionList the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventMask true if this is a bulk solvent mask
     * @param neutron true if this is a neutron structure
     */
    public CrystalReciprocalSpace(ReflectionList reflectionList,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventMask, boolean neutron) {
        this(reflectionList, atoms, fftTeam, parallelTeam, solventMask, neutron,
                SolventModel.POLYNOMIAL);
    }

    /**
     * Crystal Reciprocal Space constructor, all parameters provided
     *
     * @param reflectionlist the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventMask true if this is a bulk solvent mask
     * @param neutron true if this is a neutron structure
     * @param solventModel bulk solvent model type
     * @see CrystalReciprocalSpace.SolventModel
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventMask, boolean neutron, int solventModel) {
        this.reflectionList = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.resolution = reflectionlist.resolution;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;
        this.solvent = solventMask;
        this.solventModel = solventModel;
        this.neutron = neutron;
        this.nSymm = 1;
        threadCount = parallelTeam.getThreadCount();
        // necssary for the bulksolvent expansion!
        bulkNSymm = crystal.spaceGroup.symOps.size();
        double density = 2.0 / resolution.sampling_limit();
        double res = resolution.res_limit();

        // Set default FFT grid size from unit cell dimensions.
        int nX = (int) Math.floor(crystal.a * density / res) + 1;
        int nY = (int) Math.floor(crystal.b * density / res) + 1;
        int nZ = (int) Math.floor(crystal.c * density / res) + 1;
        if (nX % 2 != 0) {
            nX += 1;
        }
        // Use preferred dimensions.
        while (!Complex.preferredDimension(nX)) {
            nX += 2;
        }
        while (!Complex.preferredDimension(nY)) {
            nY += 1;
        }
        while (!Complex.preferredDimension(nZ)) {
            nZ += 1;
        }

        fftX = nX;
        fftY = nY;
        fftZ = nZ;
        halfFFTX = nX / 2;
        halfFFTY = nY / 2;
        halfFFTZ = nZ / 2;
        fftScale = crystal.volume;
        complexFFT3DSpace = fftX * fftY * fftZ * 2;
        densityGrid = new double[complexFFT3DSpace];

        String solventName = "none";
        if (solvent) {
            bAdd = 0.0;
            atomFormFactors = null;
            double vdwr;
            switch (solventModel) {
                case SolventModel.BINARY:
                    solventName = "binary";
                    solventA = 1.0;
                    solventB = 1.0;
                    solventGrid = new double[complexFFT3DSpace];
                    solventFormFactors = new SolventBinaryFormFactor[bulkNSymm][nAtoms];
                    for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                        for (int i = 0; i < nAtoms; i++) {
                            vdwr = atoms[i].getVDWType().radius * 0.5;
                            solventFormFactors[iSymm][i] = new SolventBinaryFormFactor(atoms[i], vdwr + solventA);
                        }
                    }
                    break;
                case SolventModel.GAUSSIAN:
                    solventName = "Gaussian";
                    solventA = 11.5;
                    solventB = 0.55;
                    solventGrid = new double[complexFFT3DSpace];
                    solventFormFactors = new SolventGaussFormFactor[bulkNSymm][nAtoms];
                    for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                        for (int i = 0; i < nAtoms; i++) {
                            vdwr = atoms[i].getVDWType().radius * 0.5;
                            solventFormFactors[iSymm][i] = new SolventGaussFormFactor(atoms[i], vdwr * solventB);
                        }
                    }
                    break;
                case SolventModel.POLYNOMIAL:
                    solventName = "polynomial switch";
                    solventA = 0.0;
                    solventB = 0.8;
                    solventGrid = new double[complexFFT3DSpace];
                    solventFormFactors = new SolventPolyFormFactor[bulkNSymm][nAtoms];
                    for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                        for (int i = 0; i < nAtoms; i++) {
                            vdwr = atoms[i].getVDWType().radius * 0.5;
                            solventFormFactors[iSymm][i] = new SolventPolyFormFactor(atoms[i], vdwr + solventA, solventB);
                        }
                    }
                    break;
                default:
                    solventGrid = null;
                    solventFormFactors = null;
                    break;
            }
        } else {
            bAdd = 2.0;
            if (neutron) {
                atomFormFactors = new NeutronFormFactor[bulkNSymm][nAtoms];
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        atomFormFactors[iSymm][i] = new NeutronFormFactor(atoms[i], bAdd);
                    }
                }
            } else {
                atomFormFactors = new XRayFormFactor[bulkNSymm][nAtoms];
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        atomFormFactors[iSymm][i] = new XRayFormFactor(atoms[i], useThreeGaussians, bAdd);
                    }
                }
            }
            solventGrid = null;
            solventFormFactors = null;
        }

        // determine number of grid points to sample density on
        aRad = -1.0;
        for (Atom a : atoms) {
            double vdwr = a.getVDWType().radius * 0.5;
            if (!solvent) {
                aRad = Math.max(aRad, a.getFormFactorWidth());
            } else {
                switch (solventModel) {
                    case SolventModel.BINARY:
                        aRad = Math.max(aRad, vdwr + solventA + 0.2);
                        break;
                    case SolventModel.GAUSSIAN:
                        aRad = Math.max(aRad, vdwr * solventB + 2.0);
                        break;
                    case SolventModel.POLYNOMIAL:
                        aRad = Math.max(aRad, vdwr + solventB + 0.2);
                        break;
                }
            }
        }
        aRadGrid = (int) Math.floor(aRad * nX / crystal.a) + 1;

        // local copy of coordinates - note use of bulknsym
        coordinates = new double[bulkNSymm][3][nAtoms];
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        for (int i = 0; i < bulkNSymm; i++) {
            double x[] = coordinates[i][0];
            double y[] = coordinates[i][1];
            double z[] = coordinates[i][2];
            for (int j = 0; j < nAtoms; j++) {
                Atom aj = atoms[j];
                crystal.applySymOp(aj.getXYZ(), xyz, symops.get(i));
                x[j] = xyz[0];
                y[j] = xyz[1];
                z[j] = xyz[2];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            if (solvent) {
                sb.append(String.format("  Bulk Solvent Grid\n"));
                sb.append(String.format("  Bulk solvent model type:       %s\n",
                        solventName));
            } else {
                sb.append(String.format("  Atomic Grid\n"));
            }
            sb.append(String.format("  Form factor grid radius (radius):  %d (%8.3f)\n",
                    aRadGrid, aRad));
            sb.append(String.format("  Grid density:               %8.3f\n",
                    density));
            sb.append(String.format("  Grid dimensions:           (%d,%d,%d)\n",
                    fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        String gridString = System.getProperty("grid-method", "SLICE").toUpperCase();
        GridMethod tempGrid;
        try {
            tempGrid = GridMethod.valueOf(gridString);
        } catch (Exception e) {
            tempGrid = GridMethod.SLICE;
        }
        gridMethod = tempGrid;

        logger.log( Level.INFO, " X-ray Refinement Parallelization Method: {0}", gridMethod.toString());

        if (solvent) {
            int minWork = nSymm;
            if (solventModel != SolventModel.NONE) {
                solventGradientRegion = new SolventGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES);
                /**
                 * Null out scattering instances used without bulk solvent.
                 */
                atomicGradientRegion = null;
                atomicDensityLoops = null;
                atomicSliceLoops = null;
                switch (gridMethod) {
                    case SPATIAL:
                        atomicDensityRegion
                                = new SpatialDensityRegion(fftX, fftY, fftZ,
                                        densityGrid, (aRadGrid + 2) * 2, nSymm, minWork,
                                        threadCount, crystal, atoms, coordinates);
                        solventDensityRegion
                                = new BulkSolventDensityRegion(fftX, fftY, fftZ,
                                        solventGrid, (aRadGrid + 2) * 2, bulkNSymm, minWork,
                                        threadCount, crystal, atoms, coordinates, 4.0, parallelTeam);
                        if (solventModel == SolventModel.GAUSSIAN) {
                            atomicDensityRegion.setInitValue(0.0);
                        } else {
                            atomicDensityRegion.setInitValue(1.0);
                        }
                        solventDensityLoops = new SolventDensityLoop[threadCount];
                        bulkSolventDensityLoops = new SolventDensityLoop[threadCount];
                        for (int i = 0; i < threadCount; i++) {
                            solventDensityLoops[i] = new SolventDensityLoop(atomicDensityRegion);
                            bulkSolventDensityLoops[i] = new SolventDensityLoop(solventDensityRegion);
                        }
                        atomicDensityRegion.setDensityLoop(solventDensityLoops);
                        solventDensityRegion.setDensityLoop(bulkSolventDensityLoops);

                        atomicSliceRegion = null;
                        solventSliceRegion = null;
                        bulkSolventSliceLoops = null;
                        solventSliceLoops = null;
                        break;
                    case SLICE:
                    default:
                        atomicSliceRegion
                                = new SliceRegion(fftX, fftY, fftZ,
                                        densityGrid, (aRadGrid + 2) * 2, nSymm,
                                        threadCount, crystal, atoms, coordinates);
                        solventSliceRegion
                                = new BulkSolventSliceRegion(fftX, fftY, fftZ,
                                        solventGrid, (aRadGrid + 2) * 2, bulkNSymm,
                                        threadCount, crystal, atoms, coordinates, 4.0, parallelTeam);
                        if (solventModel == SolventModel.GAUSSIAN) {
                            atomicSliceRegion.setInitValue(0.0);
                        } else {
                            atomicSliceRegion.setInitValue(1.0);
                        }
                        solventSliceLoops = new SolventSliceLoop[threadCount];
                        bulkSolventSliceLoops = new SolventSliceLoop[threadCount];
                        for (int i = 0; i < threadCount; i++) {
                            solventSliceLoops[i] = new SolventSliceLoop(atomicSliceRegion);
                            bulkSolventSliceLoops[i] = new SolventSliceLoop(solventSliceRegion);
                        }
                        atomicSliceRegion.setDensityLoop(solventSliceLoops);
                        solventSliceRegion.setDensityLoop(bulkSolventSliceLoops);

                        atomicDensityRegion = null;
                        solventDensityRegion = null;
                        bulkSolventDensityLoops = null;
                        solventDensityLoops = null;
                }
            } else {
                atomicDensityRegion = null;
                solventDensityRegion = null;
                atomicDensityLoops = null;
                solventDensityLoops = null;
                bulkSolventDensityLoops = null;
                atomicGradientRegion = null;
                solventGradientRegion = null;
                atomicSliceRegion = null;
                atomicSliceLoops = null;
                solventSliceRegion = null;
                bulkSolventSliceLoops = null;
                solventSliceLoops = null;
            }
        } else {
            atomicGradientRegion = new AtomicGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES);
            /**
             * Null out bulk solvent instances.
             */
            solventDensityRegion = null;
            solventDensityLoops = null;
            bulkSolventDensityLoops = null;
            solventSliceRegion = null;
            solventSliceLoops = null;
            bulkSolventSliceLoops = null;
            solventGradientRegion = null;
            /**
             * Create nSymm pieces of work per thread; the empty pieces will be
             * removed leaving 1 piece of work per thread.
             */
            switch (gridMethod) {
                case SPATIAL:
                    int minWork = nSymm;
                    atomicDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ,
                            densityGrid, (aRadGrid + 2) * 2, nSymm, minWork,
                            threadCount, crystal, atoms, coordinates);
                    atomicDensityRegion.setInitValue(0.0);
                    atomicDensityLoops = new AtomicDensityLoop[threadCount];
                    for (int i = 0; i < threadCount; i++) {
                        atomicDensityLoops[i]
                                = new AtomicDensityLoop(atomicDensityRegion);
                    }
                    atomicDensityRegion.setDensityLoop(atomicDensityLoops);
                    atomicSliceRegion = null;
                    atomicSliceLoops = null;
                    break;
                case SLICE:
                default:
                    atomicSliceRegion = new SliceRegion(fftX, fftY, fftZ,
                            densityGrid, (aRadGrid + 2) * 2, nSymm,
                            threadCount, crystal, atoms, coordinates);
                    atomicSliceRegion.setInitValue(0.0);
                    atomicSliceLoops = new AtomicSliceLoop[threadCount];
                    for (int i = 0; i < threadCount; i++) {
                        atomicSliceLoops[i] = new AtomicSliceLoop(atomicSliceRegion);
                    }
                    atomicSliceRegion.setDensityLoop(atomicSliceLoops);
                    atomicDensityRegion = null;
                    atomicDensityLoops = null;
            }
        }

        extractRegion = new ExtractRegion(threadCount);
        babinetRegion = new BabinetRegion(threadCount);
        initRegion = new InitRegion(threadCount);
        solventScaleRegion = new SolventScaleRegion(threadCount);
        atomicScaleRegion = new AtomicScaleRegion(threadCount);
        complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam);
    }

    /**
     * Set the current value of the state variable.
     *
     * @param lambda a double.
     */
    protected void setLambda(double lambda) {
        this.lambda = lambda;
    }

    /**
     * set the bulk solvent parameters
     *
     * @param a added atom width (binary, polynomial only)
     * @param b falloff of the switch (Gaussian, polynomial only)
     */
    public void setSolventAB(double a, double b) {
        double vdwr;
        this.solventA = a;
        this.solventB = b;
        if (!solvent) {
            return;
        }
        switch (solventModel) {
            case SolventModel.BINARY:
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventFormFactors[iSymm][i] = new SolventBinaryFormFactor(atoms[i], vdwr + solventA);
                    }
                }
                break;
            case SolventModel.GAUSSIAN:
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventFormFactors[iSymm][i] = new SolventGaussFormFactor(atoms[i], vdwr * solventB);
                    }
                }
                break;
            case SolventModel.POLYNOMIAL:
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventFormFactors[iSymm][i] = new SolventPolyFormFactor(atoms[i], vdwr + solventA, solventB);
                    }
                }
                break;
        }
    }

    /**
     * should the structure factor computation use 3 Gaussians or 6 for atoms?
     *
     * @param useThreeGaussians if true, use 3 Gaussian
     */
    public void setUse3G(boolean useThreeGaussians) {
        this.useThreeGaussians = useThreeGaussians;
        if (solvent || neutron) {
            return;
        }
        for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
            for (int i = 0; i < nAtoms; i++) {
                atomFormFactors[iSymm][i] = new XRayFormFactor(atoms[i], useThreeGaussians, bAdd);
            }
        }
    }

    /**
     * return dataset weight
     *
     * @return the weight wA
     */
    public double getWeight() {
        return weight;
    }

    /**
     * set the dataset weight
     *
     * @param weight desired weight wA
     */
    public void setWeight(double weight) {
        this.weight = weight;
    }

    /**
     * offset X coordinates (mostly for finite difference checks)
     *
     * @param n {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    public void deltaX(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[0] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][0][n] = symxyz[0];
        }
    }

    /**
     * offset Y coordinates (mostly for finite difference checks)
     *
     * @param n {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    public void deltaY(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[1] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][1][n] = symxyz[1];
        }
    }

    /**
     * offset Z coordinates (mostly for finite difference checks)
     *
     * @param n {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    public void deltaZ(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[2] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][2][n] = symxyz[2];
        }
    }

    /**
     * set atom coordinates
     *
     * @param coords new coordinate positions (3 params per atom)
     */
    public void setCoordinates(double coords[]) {
        assert (coords != null);
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];

        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            xyz[0] = coords[index++];
            xyz[1] = coords[index++];
            xyz[2] = coords[index++];
            for (int j = 0; j < nSymm; j++) {
                crystal.applySymOp(xyz, symxyz, symops.get(j));
                coordinates[j][0][i] = symxyz[0];
                coordinates[j][1][i] = symxyz[1];
                coordinates[j][2][i] = symxyz[2];
            }
        }
    }

    /**
     * parallelized computation of structure factors
     *
     * @param hklData structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeDensity(double hklData[][]) {
        computeDensity(hklData, false);
    }

    /**
     * parallelized computation of structure factors
     *
     * @param hklData structure factor list to fill in
     * @param print if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    public void computeDensity(double hklData[][], boolean print) {
        if (solvent) {
            if (solventModel != SolventModel.NONE) {
                computeSolventDensity(hklData, print);
            }
        } else {
            computeAtomicDensity(hklData, print);
        }
    }

    /**
     * compute inverse FFT to determine atomic gradients
     *
     * @param hklData structure factors to apply inverse FFT
     * @param freer array of free r flags corresponding to hkldata
     * @param flag Rfree flag value
     * @param refinementMode
     * {@link RefinementMinimize.RefinementMode refinement mode}
     * @see RefinementMinimize.RefinementMode
     * @see DiffractionRefinementData
     */
    public void computeAtomicGradients(double hklData[][],
            int freer[], int flag, RefinementMode refinementMode) {
        computeAtomicGradients(hklData, freer, flag, refinementMode, false);
    }

    /**
     * compute inverse FFT to determine atomic gradients
     *
     * @param hklData structure factors to apply inverse FFT
     * @param freer array of free r flags corresponding to hkldata
     * @param flag Rfree flag value
     * @param refinementMode
     * {@link RefinementMinimize.RefinementMode refinement mode}
     * @param print if true, print information on timings during the calculation
     * @see RefinementMinimize.RefinementMode
     * @see DiffractionRefinementData
     */
    public void computeAtomicGradients(double hklData[][],
            int freer[], int flag, RefinementMode refinementMode,
            boolean print) {

        if (solvent && solventModel == SolventModel.NONE) {
            return;
        }

        // zero out the density
        for (int i = 0; i < complexFFT3DSpace; i++) {
            densityGrid[i] = 0.0;
        }

        int nfree = 0;
        StringBuilder sb = new StringBuilder();
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        // int nsym = 1;
        List<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        ComplexNumber cj = new ComplexNumber();
        HKL ij = new HKL();
        for (HKL ih : reflectionList.hkllist) {
            double fc[] = hklData[ih.index()];
            if (Double.isNaN(fc[0])) {
                continue;
            }

            // cross validation check!!!
            if (freer != null) {
                if (freer[ih.index()] == flag) {
                    nfree++;
                    continue;
                }
            }

            c.re(fc[0]);
            c.im(fc[1]);
            // scale
            c.times_ip(2.0 / fftScale);

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                cj.copy(c);
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    cj.phase_shift_ip(shift);
                    densityGrid[ii] += cj.re();
                    densityGrid[ii + 1] += -cj.im();
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    cj.phase_shift_ip(shift);
                    densityGrid[ii] += cj.re();
                    densityGrid[ii + 1] += cj.im();
                }
            }
        }
        symtime += System.nanoTime();

        long startTime = System.nanoTime();
        complexFFT3D.ifft(densityGrid);
        long fftTime = System.nanoTime() - startTime;


        /*
         CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
         mapout.write(densityGrid);
         */
        startTime = System.nanoTime();
        long permanentDensityTime = 0;
        try {
            if (solvent) {
                solventGradientRegion.setRefinementMode(refinementMode);
                parallelTeam.execute(solventGradientRegion);
            } else {
                atomicGradientRegion.setRefinementMode(refinementMode);
                parallelTeam.execute(atomicGradientRegion);
            }
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Exception computing atomic gradients.";
            logger.log(Level.SEVERE, message, e);
        }

        if (solvent) {
            sb.append(String.format(" Solvent symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(String.format(" Solvent inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(String.format(" Grid solvent gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)\n", nfree));
        } else {
            sb.append(String.format(" Atomic symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(String.format(" Atomic inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(String.format(" Grid atomic gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)\n", nfree));
        }

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    /**
     * parallelized computation of structure factors (identical to
     * compuateDensity)
     *
     * @param hklData structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeAtomicDensity(double hklData[][]) {
        computeAtomicDensity(hklData, false);
    }

    /**
     * parallelized computation of structure factors (identical to
     * computeDensity)
     *
     * @param hklData structure factor list to fill in
     * @param print if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    public void computeAtomicDensity(double hklData[][], boolean print) {
        /**
         * Zero out reflection data.
         */
        long initTime = -System.nanoTime();
        try {
            initRegion.setHKL(hklData);
            parallelTeam.execute(initRegion);
        } catch (Exception e) {
            String message = "Fatal exception zeroing out reflection data.";
            logger.log(Level.SEVERE, message, e);
        }
        initTime += System.nanoTime();

        /**
         * Assign atomic electron density to FFT grid.
         */
        long atomicGridTime = -System.nanoTime();
        try {
            switch (gridMethod) {
                case SPATIAL:
                    atomicDensityRegion.assignAtomsToCells();
                    parallelTeam.execute(atomicDensityRegion);
                    break;
                case SLICE:
                default:
                    parallelTeam.execute(atomicSliceRegion);
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating atomic electron density.";
            logger.log(Level.SEVERE, message, e);
        }
        atomicGridTime += System.nanoTime();

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(densityGrid);
        /**
         * Compute model structure factors via an FFT of the electron density.
         */
        long startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        /**
         * Extract and scale structure factors.
         */
        long symTime = -System.nanoTime();
        try {
            extractRegion.setHKL(hklData);
            parallelTeam.execute(extractRegion);
            double scale = (fftScale) / (fftX * fftY * fftZ);
            atomicScaleRegion.setHKL(hklData);
            atomicScaleRegion.setScale(scale);
            parallelTeam.execute(atomicScaleRegion);
        } catch (Exception e) {
            String message = "Fatal exception extracting structure factors.";
            logger.log(Level.SEVERE, message, e);
        }
        symTime += System.nanoTime();

        if (logger.isLoggable(Level.INFO) && print) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\n Fc Initialization:         %8.4f\n", initTime * toSeconds));
            sb.append(String.format(" Atomic grid density:       %8.4f\n", atomicGridTime * toSeconds));
            sb.append(String.format(" Atomic FFT:                %8.4f\n", fftTime * toSeconds));
            sb.append(String.format(" Atomic symmetry & scaling: %8.4f\n", symTime * toSeconds));
            logger.info(sb.toString());
        }
    }

    /**
     * parallelized computation of bulk solvent structure factors
     *
     * @param hklData structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeSolventDensity(double hklData[][]) {
        computeSolventDensity(hklData, false);
    }

    /**
     * parallelized computation of bulk solvent structure factors
     *
     * @param hklData structure factor list to fill in
     * @param print if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    public void computeSolventDensity(double hklData[][], boolean print) {
        /**
         * Zero out the reflection data.
         */
        long initTime = -System.nanoTime();
        try {
            initRegion.setHKL(hklData);
            parallelTeam.execute(initRegion);
        } catch (Exception e) {
            String message = "Fatal exception initializing structure factors.";
            logger.log(Level.SEVERE, message, e);
        }
        initTime += System.nanoTime();

        long solventGridTime = -System.nanoTime();
        try {
            switch (gridMethod) {
                case SPATIAL:
                    atomicDensityRegion.assignAtomsToCells();
                    parallelTeam.execute(atomicDensityRegion);
                    break;
                case SLICE:
                default:
                    parallelTeam.execute(atomicSliceRegion);
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }
        solventGridTime += System.nanoTime();

        long expTime = -System.nanoTime();
        if (solventModel == SolventModel.BINARY) {
            /*
             * need to shrink mask
             * first, copy densitygrid to df_map
             * temporarily to store original map
             */
            int nmap = densityGrid.length;
            System.arraycopy(densityGrid, 0, solventGrid, 0, nmap);
            // Logic to loop within the cutoff box.
            double xyz[] = {solventB, solventB, solventB};
            double uvw[] = new double[3];
            crystal.toFractionalCoordinates(xyz, uvw);
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        if (densityGrid[ii] == 1.0) {
                            continue;
                        }
                        for (int iz = k - ifrz; iz <= k + ifrz; iz++) {
                            int giz = Crystal.mod(iz, fftZ);
                            for (int iy = j - ifry; iy <= j + ifry; iy++) {
                                int giy = Crystal.mod(iy, fftY);
                                for (int ix = i - ifrx; ix <= i + ifrx; ix++) {
                                    int gix = Crystal.mod(ix, fftX);
                                    final int jj = iComplex3D(gix, giy, giz, fftX, fftY);
                                    if (densityGrid[jj] == 1.0) {
                                        solventGrid[ii] = 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Copy the completed solvent grid back.
            System.arraycopy(solventGrid, 0, densityGrid, 0, nmap);
        }

        // Copy the grid over for derivatives.
        System.arraycopy(densityGrid, 0, solventGrid, 0, densityGrid.length);
        /**
         * Babinet Principle.
         */
        try {
            parallelTeam.execute(babinetRegion);
        } catch (Exception e) {
            String message = "Fatal exception Babinet Principle.";
            logger.log(Level.SEVERE, message, e);
        }
        expTime += System.nanoTime();

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(densityGrid);
        long fftTime = -System.nanoTime();
        complexFFT3D.fft(densityGrid);
        fftTime += System.nanoTime();

        /**
         * Extract and scale structure factors.
         */
        long symTime = -System.nanoTime();
        try {
            extractRegion.setHKL(hklData);
            parallelTeam.execute(extractRegion);
            final double scale = (fftScale) / (fftX * fftY * fftZ);
            solventScaleRegion.setHKL(hklData);
            solventScaleRegion.setScale(scale);
            parallelTeam.execute(solventScaleRegion);
        } catch (Exception e) {
            String message = "Fatal exception extracting structure factors.";
            logger.log(Level.SEVERE, message, e);
        }
        symTime += System.nanoTime();

        /**
         * Expand derivative mask so it includes nearby local density.
         */
        long solventDensityTime = -System.nanoTime();
        try {
            switch (gridMethod) {
                case SPATIAL:
                    solventDensityRegion.assignAtomsToCells();
                    //parallelTeam.execute(solventDensityRegion);
                    break;
                case SLICE:
                default:
                //parallelTeam.execute(solventSliceRegion);
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }
        solventDensityTime += System.nanoTime();

        long solventExpTime = -System.nanoTime();
        if (solventModel == SolventModel.GAUSSIAN) {
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        solventGrid[ii] = exp(-solventA * solventGrid[ii]);
                    }
                }
            }
        }
        solventExpTime += System.nanoTime();

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(solventGrid);
        if (logger.isLoggable(Level.INFO) && print) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("\n Fc Initialization:                     %8.4f\n", initTime * toSeconds));
            sb.append(String.format(" Solvent grid density:                  %8.4f\n", solventGridTime * toSeconds));
            sb.append(String.format(" Bulk solvent exponentiation:           %8.4f\n", expTime * toSeconds));
            sb.append(String.format(" Solvent FFT:                           %8.4f\n", fftTime * toSeconds));
            sb.append(String.format(" Solvent symmetry & scaling:            %8.4f\n", symTime * toSeconds));
            sb.append(String.format(" Solvent grid expansion:                %8.4f\n", solventDensityTime * toSeconds));
            sb.append(String.format(" Solvent grid expansion exponentiation: %8.4f\n", solventExpTime * toSeconds));
            logger.info(sb.toString());
        }
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
     * <p>
     * densityNorm</p>
     *
     * @param data an array of double.
     * @param meansd an array of double.
     * @param norm a boolean.
     */
    public void densityNorm(double data[], double meansd[], boolean norm) {
        double mean, sd;

        mean = sd = 0.0;
        int n = 0;
        for (int k = 0; k < fftZ; k++) {
            for (int j = 0; j < fftY; j++) {
                for (int i = 0; i < fftX; i++) {
                    int index = iComplex3D(i, j, k, fftX, fftY);
                    n++;
                    mean += (data[index] - mean) / n;
                }
            }
        }

        n = 0;
        for (int k = 0; k < fftZ; k++) {
            for (int j = 0; j < fftY; j++) {
                for (int i = 0; i < fftX; i++) {
                    int index = iComplex3D(i, j, k, fftX, fftY);
                    sd += Math.pow(data[index] - mean, 2.0);
                    n++;
                }
            }
        }
        sd = Math.sqrt(sd / n);

        if (meansd != null) {
            meansd[0] = mean;
            meansd[1] = sd;
        }

        if (norm) {
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        int index = iComplex3D(i, j, k, fftX, fftY);
                        data[index] = (data[index] - mean) / sd;
                    }
                }
            }
        }
    }

    private class AtomicDensityLoop extends SpatialDensityLoop {

        final double xyz[] = new double[3];
        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public AtomicDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            if (!atoms[n].isActive()) {
                return;
            }

            final double lambdai = atoms[n].applyLambda() ? lambda : 1.0;
            xyz[0] = coordinates[iSymm][0][n];
            xyz[1] = coordinates[iSymm][1][n];
            xyz[2] = coordinates[iSymm][2][n];
            FormFactor atomff = atomFormFactors[iSymm][n];
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = Math.min(aRadGrid,
                    (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                xf[2] = iz / (double) fftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = atomff.rho(grid[ii], lambdai, xc);
                    }
                }
            }
        }
    }

    private class SolventDensityLoop extends SpatialDensityLoop {

        final double xyz[] = new double[3];
        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public SolventDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            if (!atoms[n].isActive()) {
                return;
            }
            final double lambdai = atoms[n].applyLambda() ? lambda : 1.0;
            xyz[0] = coordinates[iSymm][0][n];
            xyz[1] = coordinates[iSymm][1][n];
            xyz[2] = coordinates[iSymm][2][n];
            FormFactor solventff = solventFormFactors[iSymm][n];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[n].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case SolventModel.BINARY:
                    frad = Math.min(aRadGrid,
                            (int) Math.floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.GAUSSIAN:
                    frad = Math.min(aRadGrid,
                            (int) Math.floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.POLYNOMIAL:
                    frad = Math.min(aRadGrid,
                            (int) Math.floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
            }

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                xf[2] = iz / (double) fftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = solventff.rho(grid[ii], lambdai, xc);
                    }
                }
            }
        }
    }

    private class AtomicSliceLoop extends SliceLoop {

        final double xyz[] = new double[3];
        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public AtomicSliceLoop(SliceRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].isActive()) {
                return;
            }

            final double lambdai = atoms[iAtom].applyLambda() ? lambda : 1.0;

            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            FormFactor atomff = atomFormFactors[iSymm][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = min(aRadGrid, (int) floor(atoms[iAtom].getFormFactorWidth() * fftX / crystal.a) + 1);

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                if (lb <= giz && giz <= ub) {
                    xf[2] = iz / (double) fftZ;
                    for (int iy = ifry - frad; iy <= ifryu; iy++) {
                        int giy = Crystal.mod(iy, fftY);
                        xf[1] = iy / (double) fftY;
                        for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                            int gix = Crystal.mod(ix, fftX);
                            xf[0] = ix / (double) fftX;
                            crystal.toCartesianCoordinates(xf, xc);
                            final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                            grid[ii] = atomff.rho(grid[ii], lambdai, xc);
                        }
                    }
                }
            }
        }
    }

    private class SolventSliceLoop extends SliceLoop {

        final double xyz[] = new double[3];
        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public SolventSliceLoop(SliceRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].isActive()) {
                return;
            }
            final double lambdai = atoms[iAtom].applyLambda() ? lambda : 1.0;
            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            FormFactor formFactor = solventFormFactors[iSymm][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[iAtom].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case SolventModel.BINARY:
                    frad = min(aRadGrid, (int) floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.GAUSSIAN:
                    frad = min(aRadGrid, (int) floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.POLYNOMIAL:
                    frad = min(aRadGrid, (int) floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
            }

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                if (lb > giz || giz > ub) {
                    continue;
                }
                xf[2] = iz / (double) fftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = formFactor.rho(grid[ii], lambdai, xc);
                    }
                }
            }
        }
    }

    private class AtomicGradientRegion extends ParallelRegion {

        private final AtomicGradientLoop atomicGradientLoop[];
        private RefinementMode refinementmode;

        public AtomicGradientRegion(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
            atomicGradientLoop = new AtomicGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                atomicGradientLoop[i] = new AtomicGradientLoop();
            }
        }

        public RefinementMode getRefinementMode() {
            return refinementmode;
        }

        public void setRefinementMode(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, atomicGradientLoop[getThreadIndex()]);
            } catch (Exception e) {
                e.printStackTrace();
                logger.severe(e.toString());
            }
        }

        private class AtomicGradientLoop extends IntegerForLoop {

            final double xyz[] = new double[3];
            final double uvw[] = new double[3];
            final double xc[] = new double[3];
            final double xf[] = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    if (!atoms[n].isActive()) {
                        continue;
                    }
                    if (lambdaTerm) {
                        if (!atoms[n].applyLambda()) {
                            continue;
                        }
                    }
                    xyz[0] = coordinates[0][0][n];
                    xyz[1] = coordinates[0][1][n];
                    xyz[2] = coordinates[0][2][n];
                    FormFactor atomff = atomFormFactors[0][n];
                    atomff.update(xyz, 0.0);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    final int dfrad = Math.min(aRadGrid,
                            (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;

                    for (int ix = ifrx - dfrad; ix <= ifrx + dfrad; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        for (int iy = ifry - dfrad; iy <= ifry + dfrad; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy / (double) fftY;
                            for (int iz = ifrz - dfrad; iz <= ifrz + dfrad; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                atomff.rho_grad(xc, weight * densityGrid[ii], refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    private class SolventGradientRegion extends ParallelRegion {

        private final SolventGradientLoop solventGradientLoop[];
        private RefinementMode refinementmode;

        public SolventGradientRegion(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
            solventGradientLoop = new SolventGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                solventGradientLoop[i] = new SolventGradientLoop();
            }
        }

        public RefinementMode getRefinementMode() {
            return refinementmode;
        }

        public void setRefinementMode(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, solventGradientLoop[getThreadIndex()]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class SolventGradientLoop extends IntegerForLoop {

            double xyz[] = new double[3];
            double uvw[] = new double[3];
            double xc[] = new double[3];
            double xf[] = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    if (!atoms[n].isActive()) {
                        continue;
                    }
                    if (lambdaTerm) {
                        if (!atoms[n].applyLambda()) {
                            continue;
                        }
                    }
                    xyz[0] = coordinates[0][0][n];
                    xyz[1] = coordinates[0][1][n];
                    xyz[2] = coordinates[0][2][n];
                    FormFactor solventff = solventFormFactors[0][n];
                    solventff.update(xyz);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    double vdwr = atoms[n].getVDWType().radius * 0.5;
                    double dfcmult = 1.0;
                    int dfrad = aRadGrid;
                    switch (solventModel) {
                        case SolventModel.BINARY:
                            dfrad = Math.min(aRadGrid,
                                    (int) Math.floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                            break;
                        case SolventModel.GAUSSIAN:
                            dfrad = Math.min(aRadGrid,
                                    (int) Math.floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                            dfcmult = solventA;
                            break;
                        case SolventModel.POLYNOMIAL:
                            dfrad = Math.min(aRadGrid,
                                    (int) Math.floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                            break;
                    }

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;
                    final int ifrxu = ifrx + dfrad;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;
                    final int ifryu = ifry + dfrad;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;
                    final int ifrzu = ifrz + dfrad;

                    for (int ix = ifrx - dfrad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        for (int iy = ifry - dfrad; iy <= ifryu; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy / (double) fftY;
                            for (int iz = ifrz - dfrad; iz <= ifrzu; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                solventff.rho_grad(xc,
                                        weight * densityGrid[ii] * dfcmult * solventGrid[ii],
                                        refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    private class BabinetRegion extends ParallelRegion {

        BabinetLoop babinetLoops[];

        public BabinetRegion(int nThreads) {
            babinetLoops = new BabinetLoop[nThreads];

        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (babinetLoops[ti] == null) {
                babinetLoops[ti] = new BabinetLoop();
            }

            try {
                execute(0, fftZ - 1, babinetLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class BabinetLoop extends IntegerForLoop {

            @Override
            public void run(int lb, int ub) throws Exception {
                switch (solventModel) {
                    case BINARY:
                    case POLYNOMIAL:
                        for (int k = lb; k <= ub; k++) {
                            for (int j = 0; j < fftY; j++) {
                                for (int i = 0; i < fftX; i++) {
                                    final int ii = iComplex3D(i, j, k, fftX, fftY);
                                    densityGrid[ii] = 1.0 - densityGrid[ii];
                                }
                            }
                        }
                        break;
                    case GAUSSIAN:
                        for (int k = lb; k <= ub; k++) {
                            for (int j = 0; j < fftY; j++) {
                                for (int i = 0; i < fftX; i++) {
                                    final int ii = iComplex3D(i, j, k, fftX, fftY);
                                    densityGrid[ii] = 1.0 - exp(-solventA * densityGrid[ii]);
                                }
                            }
                        }
                        break;
                    default:
                }
            }
        }
    }

    private class InitRegion extends ParallelRegion {

        InitLoop initLoops[];
        FormFactorUpdateLoop formFactorUpdateLoops[];
        int nHKL = reflectionList.hkllist.size();
        double hkldata[][] = null;

        public InitRegion(int nThreads) {
            initLoops = new InitLoop[nThreads];
            formFactorUpdateLoops = new FormFactorUpdateLoop[nThreads];
        }

        public void setHKL(double hkldata[][]) {
            this.hkldata = hkldata;
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (initLoops[ti] == null) {
                initLoops[ti] = new InitLoop();
                formFactorUpdateLoops[ti] = new FormFactorUpdateLoop();
            }

            try {
                execute(0, nHKL - 1, initLoops[ti]);
                execute(0, nAtoms - 1, formFactorUpdateLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class InitLoop extends IntegerForLoop {

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int i = lb; i <= ub; i++) {
                    hkldata[i][0] = 0.0;
                    hkldata[i][1] = 0.0;
                }
            }
        }

        private class FormFactorUpdateLoop extends IntegerForLoop {

            private final double xyz[] = new double[3];

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = lb; i <= ub; i++) {
                        xyz[0] = coordinates[iSymm][0][i];
                        xyz[1] = coordinates[iSymm][1][i];
                        xyz[2] = coordinates[iSymm][2][i];
                        if (solventFormFactors != null) {
                            solventFormFactors[iSymm][i].update(xyz);
                        }
                        if (atomFormFactors != null) {
                            atomFormFactors[iSymm][i].update(xyz, bAdd);
                        }
                    }
                }
            }
        }
    }

    private class AtomicScaleRegion extends ParallelRegion {

        AtomicScaleLoop atomicScaleLoops[];
        int nHKL = reflectionList.hkllist.size();
        double hklData[][] = null;
        double scale;

        public AtomicScaleRegion(int nThreads) {
            atomicScaleLoops = new AtomicScaleLoop[nThreads];
        }

        public void setHKL(double hkldata[][]) {
            this.hklData = hkldata;
        }

        public void setScale(double scale) {
            this.scale = scale;
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (atomicScaleLoops[ti] == null) {
                atomicScaleLoops[ti] = new AtomicScaleLoop();
            }

            try {
                execute(0, nHKL - 1, atomicScaleLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class AtomicScaleLoop extends IntegerForLoop {

            ComplexNumber c;

            public AtomicScaleLoop() {
                c = new ComplexNumber();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int i = lb; i <= ub; i++) {
                    HKL ih = reflectionList.hkllist.get(i);
                    double fc[] = hklData[ih.index()];
                    c.re(fc[0]);
                    c.im(fc[1]);
                    // Remove Badd
                    double s = Crystal.invressq(crystal, ih);
                    c.times_ip(scale * exp(0.25 * bAdd * s));
                    c.conjugate_ip();
                    fc[0] = c.re();
                    fc[1] = c.im();
                }
            }
        }
    }

    private class SolventScaleRegion extends ParallelRegion {

        SolventScaleLoop solventScaleLoops[];
        int nHKL = reflectionList.hkllist.size();
        double hkldata[][] = null;
        double scale;

        public SolventScaleRegion(int nThreads) {
            solventScaleLoops = new SolventScaleLoop[nThreads];
        }

        public void setHKL(double hkldata[][]) {
            this.hkldata = hkldata;
        }

        public void setScale(double scale) {
            this.scale = scale;
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (solventScaleLoops[ti] == null) {
                solventScaleLoops[ti] = new SolventScaleLoop();
            }

            try {
                execute(0, nHKL - 1, solventScaleLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class SolventScaleLoop extends IntegerForLoop {

            ComplexNumber c;

            public SolventScaleLoop() {
                c = new ComplexNumber();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int i = lb; i <= ub; i++) {
                    HKL ih = reflectionList.hkllist.get(i);
                    double fc[] = hkldata[ih.index()];
                    c.re(fc[0]);
                    c.im(fc[1]);
                    c.times_ip(scale);
                    c.conjugate_ip();
                    // negative: babinet
                    fc[0] = -c.re();
                    fc[1] = -c.im();
                }
            }
        }
    }

    private class ExtractRegion extends ParallelRegion {

        ExtractLoop extractLoops[];
        int nHKL = reflectionList.hkllist.size();
        double hkldata[][] = null;

        public ExtractRegion(int nThreads) {
            extractLoops = new ExtractLoop[nThreads];
        }

        public void setHKL(double hkldata[][]) {
            this.hkldata = hkldata;
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (extractLoops[ti] == null) {
                extractLoops[ti] = new ExtractLoop();
            }

            try {
                execute(0, nHKL - 1, extractLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class ExtractLoop extends IntegerForLoop {

            ComplexNumber c;
            int nsym;
            List<SymOp> symops;
            HKL ij;

            public ExtractLoop() {
                c = new ComplexNumber();
                nsym = crystal.spaceGroup.symOps.size();
                symops = crystal.spaceGroup.symOps;
                ij = new HKL();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int i = lb; i <= ub; i++) {
                    HKL ih = reflectionList.hkllist.get(i);
                    double fc[] = hkldata[ih.index()];
                    // Apply symmetry
                    for (int j = 0; j < nsym; j++) {
                        crystal.applyTransSymRot(ih, ij, symops.get(j));
                        double shift = Crystal.sym_phase_shift(ih, symops.get(j));
                        int h = Crystal.mod(ij.h(), fftX);
                        int k = Crystal.mod(ij.k(), fftY);
                        int l = Crystal.mod(ij.l(), fftZ);
                        if (h < halfFFTX + 1) {
                            final int ii = iComplex3D(h, k, l, fftX, fftY);
                            c.re(densityGrid[ii]);
                            c.im(densityGrid[ii + 1]);
                            c.phase_shift_ip(shift);
                            fc[0] += c.re();
                            fc[1] += c.im();
                        } else {
                            h = (fftX - h) % fftX;
                            k = (fftY - k) % fftY;
                            l = (fftZ - l) % fftZ;
                            final int ii = iComplex3D(h, k, l, fftX, fftY);
                            c.re(densityGrid[ii]);
                            c.im(-densityGrid[ii + 1]);
                            c.phase_shift_ip(shift);
                            fc[0] += c.re();
                            fc[1] += c.im();
                        }
                    }
                }
            }
        }
    }

}
