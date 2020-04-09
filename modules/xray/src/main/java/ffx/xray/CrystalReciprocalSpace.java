//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.xray;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedIntegerArray;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.crystal.SymOp;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Complex3DParallel;
import ffx.numerics.math.ComplexNumber;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.RowLoop;
import ffx.potential.nonbonded.RowRegion;
import ffx.potential.nonbonded.SliceLoop;
import ffx.potential.nonbonded.SliceRegion;
import ffx.potential.nonbonded.SpatialDensityLoop;
import ffx.potential.nonbonded.SpatialDensityRegion;
import ffx.xray.RefinementMinimize.RefinementMode;
import static ffx.numerics.fft.Complex3D.iComplex3D;

/**
 * Structure factor calculation (including bulk solvent structure factors)
 *
 * @author Timothy D. Fenn
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @since 1.0
 */
public class CrystalReciprocalSpace {

    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());
    private static double toSeconds = 1.0e-9;
    private final boolean neutron;
    private final double bAdd;
    private final double fftScale;
    private final double[] densityGrid;
    private final double[] solventGrid;
    private final int nSymm;
    private final int bulkNSymm;
    private final int nAtoms;
    private final int fftX, fftY, fftZ;
    private final double ifftX, ifftY, ifftZ;
    private final int halfFFTX;
    private final int complexFFT3DSpace;
    private final int threadCount;
    private final int bufferSize = 3;
    // Slice Scheduling
    private final SharedIntegerArray optSliceWeightAtomic;
    private final SharedIntegerArray optSliceWeightSolvent;
    private final SharedIntegerArray optSliceWeightBulkSolvent;
    private final int[] previousOptSliceWeightAtomic;
    private final int[] previousOptSliceWeightSolvent;
    private final int[] previousOptSliceWeightBulkSolvent;
    private final SliceSchedule atomicSliceSchedule;
    private final SliceSchedule solventSliceSchedule;
    private final SliceSchedule bulkSolventSliceSchedule;
    // Row Scheduling
    private final SharedIntegerArray optRowWeightAtomic;
    private final SharedIntegerArray optRowWeightSolvent;
    private final SharedIntegerArray optRowWeightBulkSolvent;
    private final int[] previousOptRowWeightAtomic;
    private final int[] previousOptRowWeightSolvent;
    private final int[] previousOptRowWeightBulkSolvent;
    private final RowSchedule atomicRowSchedule;
    private final RowSchedule solventRowSchedule;
    private final RowSchedule bulkSolventRowSchedule;

    private final ParallelTeam parallelTeam;
    private final Atom[] atoms;
    private final Crystal crystal;
    private final ReflectionList reflectionList;
    private final FormFactor[][] atomFormFactors;
    private final FormFactor[][] solventFormFactors;
    private final GridMethod gridMethod;
    private final SharedIntegerArray optAtomicGradientWeight;
    private final int[] previousOptAtomicGradientWeight;
    private final GradientSchedule atomicGradientSchedule;
    /**
     * Parallelization of putting atomic form factors onto the 3D grid using a
     * 3D spatial decomposition.
     */
    private final SpatialDensityRegion atomicDensityRegion;
    private final SpatialDensityRegion solventDensityRegion;
    /**
     * Parallelization of putting atomic form factors onto the 3D grid using a
     * slice-based decomposition.
     */
    private final SliceRegion atomicSliceRegion;
    private final AtomicSliceLoop[] atomicSliceLoops;
    private final SliceRegion solventSliceRegion;
    private final SolventSliceLoop[] solventSliceLoops;
    private final BulkSolventSliceLoop[] bulkSolventSliceLoops;
    /**
     * Parallelization of putting atomic form factors onto the 3D grid using a
     * row-based decomposition.
     */
    private final RowRegion atomicRowRegion;
    private final AtomicRowLoop[] atomicRowLoops;
    private final RowRegion solventRowRegion;
    private final SolventRowLoop[] solventRowLoops;
    private final BulkSolventRowLoop[] bulkSolventRowLoops;
    /**
     * Parallelization over atomic structure factor loops.
     */
    private final InitRegion initRegion;
    private final ExtractRegion extractRegion;
    private final AtomicScaleRegion atomicScaleRegion;
    private final AtomicGradientRegion atomicGradientRegion;
    /**
     * Parallelization over solvent structure factor loops.
     */
    private final BabinetRegion babinetRegion;
    private final SolventGridRegion solventGridRegion;
    private final SolventScaleRegion solventScaleRegion;
    private final SolventGradientRegion solventGradientRegion;
    private final Complex3DParallel complexFFT3D;
    boolean lambdaTerm = false;
    double solventA;
    double solventB;
    protected SolventModel solventModel;
    private boolean solvent;
    private boolean useThreeGaussians = true;
    // not final for purposes of finite differences
    private double[][][] coordinates;
    private double weight = 1.0;
    private int aRadGrid;
    /**
     * If the "Native Environment Approximation" is true, the "use" flag is ignored.
     */
    private boolean nativeEnvironmentApproximation = false;

    /**
     * Crystal Reciprocal Space constructor, assumes this is not a bulk solvent
     * mask and is not a neutron data set
     *
     * @param reflectionList the {@link ffx.crystal.ReflectionList} to fill with structure factors
     * @param atoms          array of {@link ffx.potential.bonded.Atom atoms} for structure factor computation
     * @param fftTeam        {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam   {@link edu.rit.pj.ParallelTeam} for parallelization
     */
    public CrystalReciprocalSpace(ReflectionList reflectionList,
                                  Atom[] atoms, ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this(reflectionList, atoms, fftTeam, parallelTeam, false, false,
                SolventModel.POLYNOMIAL, GridMethod.SLICE);
    }

    /**
     * Crystal Reciprocal Space constructor, assumes this is not a neutron data
     * set and implements a polynomial bulk solvent mask if needed
     *
     * @param reflectionList the {@link ffx.crystal.ReflectionList} to fill with structure factors
     * @param atoms          array of {@link ffx.potential.bonded.Atom atoms} for structure factor computation
     * @param fftTeam        {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam   {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventMask    true if this is a bulk solvent mask
     */
    public CrystalReciprocalSpace(ReflectionList reflectionList, Atom[] atoms,
                                  ParallelTeam fftTeam, ParallelTeam parallelTeam,
                                  boolean solventMask) {
        this(reflectionList, atoms, fftTeam, parallelTeam, solventMask, false,
                SolventModel.POLYNOMIAL, GridMethod.SLICE);
    }

    /**
     * Crystal Reciprocal Space constructor, assumes a polynomial bulk solvent
     * mask if needed
     *
     * @param reflectionList the {@link ffx.crystal.ReflectionList} to fill with structure factors
     * @param atoms          array of {@link ffx.potential.bonded.Atom atoms} for structure factor computation
     * @param fftTeam        {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam   {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventMask    true if this is a bulk solvent mask
     * @param neutron        true if this is a neutron structure
     */
    public CrystalReciprocalSpace(ReflectionList reflectionList,
                                  Atom[] atoms, ParallelTeam fftTeam, ParallelTeam parallelTeam,
                                  boolean solventMask, boolean neutron) {
        this(reflectionList, atoms, fftTeam, parallelTeam, solventMask, neutron,
                SolventModel.POLYNOMIAL, GridMethod.SLICE);
    }

    /**
     * Crystal Reciprocal Space constructor, all parameters provided
     *
     * @param reflectionlist the {@link ffx.crystal.ReflectionList} to fill with structure factors
     * @param atoms          array of {@link ffx.potential.bonded.Atom atoms} for structure factor computation
     * @param fftTeam        {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam   {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventMask    true if this is a bulk solvent mask
     * @param neutron        true if this is a neutron structure
     * @param solventModel   bulk solvent model type
     * @param gridMethod     parallel method for putting density onto the FFT grid.
     * @see CrystalReciprocalSpace.SolventModel
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist, Atom[] atoms,
                                  ParallelTeam fftTeam, ParallelTeam parallelTeam,
                                  boolean solventMask, boolean neutron,
                                  SolventModel solventModel, GridMethod gridMethod) {
        this.reflectionList = reflectionlist;
        this.atoms = atoms;
        this.parallelTeam = parallelTeam;
        this.solvent = solventMask;
        this.neutron = neutron;
        this.solventModel = solventModel;
        this.gridMethod = gridMethod;

        nAtoms = atoms.length;
        crystal = reflectionlist.crystal;
        nSymm = 1;
        threadCount = parallelTeam.getThreadCount();

        // Necssary for the bulk-solvent expansion!
        bulkNSymm = crystal.spaceGroup.symOps.size();
        Resolution resolution = reflectionlist.resolution;
        double gridFactor = resolution.samplingLimit() / 2.0;
        double gridStep = resolution.resolutionLimit() * gridFactor;
        // Set default FFT grid size from unit cell dimensions.
        int nX = (int) Math.floor(crystal.a / gridStep) + 1;
        int nY = (int) Math.floor(crystal.b / gridStep) + 1;
        int nZ = (int) Math.floor(crystal.c / gridStep) + 1;
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
        ifftX = 1.0 / (double) fftX;
        ifftY = 1.0 / (double) fftY;
        ifftZ = 1.0 / (double) fftZ;
        fftScale = crystal.volume;
        complexFFT3DSpace = fftX * fftY * fftZ * 2;
        densityGrid = new double[complexFFT3DSpace];

        // Slice Method Variables
        optSliceWeightAtomic = new SharedIntegerArray(fftZ);
        optSliceWeightSolvent = new SharedIntegerArray(fftZ);
        optSliceWeightBulkSolvent = new SharedIntegerArray(fftZ);
        previousOptSliceWeightAtomic = new int[fftZ];
        previousOptSliceWeightSolvent = new int[fftZ];
        previousOptSliceWeightBulkSolvent = new int[fftZ];

        atomicSliceSchedule = new SliceSchedule(threadCount, fftZ);
        solventSliceSchedule = new SliceSchedule(threadCount, fftZ);
        bulkSolventSliceSchedule = new SliceSchedule(threadCount, fftZ);

        // Row Method Variables
        optRowWeightAtomic = new SharedIntegerArray(fftZ * fftY);
        optRowWeightSolvent = new SharedIntegerArray(fftZ * fftY);
        optRowWeightBulkSolvent = new SharedIntegerArray(fftZ * fftY);
        previousOptRowWeightAtomic = new int[fftZ * fftY];
        previousOptRowWeightSolvent = new int[fftZ * fftY];
        previousOptRowWeightBulkSolvent = new int[fftZ * fftY];
        atomicRowSchedule = new RowSchedule(threadCount, fftZ, fftY);
        solventRowSchedule = new RowSchedule(threadCount, fftZ, fftY);
        bulkSolventRowSchedule = new RowSchedule(threadCount, fftZ, fftY);

        optAtomicGradientWeight = new SharedIntegerArray(nAtoms);
        previousOptAtomicGradientWeight = new int[nAtoms];
        atomicGradientSchedule = new GradientSchedule(threadCount, nAtoms);

        if (solvent) {
            bAdd = 0.0;
            atomFormFactors = null;
            double vdwr;
            switch (solventModel) {
                case BINARY:
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
                case GAUSSIAN:
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
                case POLYNOMIAL:
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
                case NONE:
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

        // Determine number of grid points to sample density on
        double aRad = -1.0;
        for (Atom a : atoms) {
            double vdwr = a.getVDWType().radius * 0.5;
            if (!solvent) {
                aRad = Math.max(aRad, a.getFormFactorWidth());
            } else {
                switch (solventModel) {
                    case BINARY:
                        aRad = Math.max(aRad, vdwr + solventA + 0.2);
                        break;
                    case GAUSSIAN:
                        aRad = Math.max(aRad, vdwr * solventB + 2.0);
                        break;
                    case POLYNOMIAL:
                        aRad = Math.max(aRad, vdwr + solventB + 0.2);
                        break;
                    case NONE:
                    default:
                        aRad = 0.0;
                }
            }
        }
        aRadGrid = (int) Math.floor(aRad * nX / crystal.a) + 1;

        // local copy of coordinates - note use of bulknsym
        coordinates = new double[bulkNSymm][3][nAtoms];
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double[] xyz = new double[3];
        for (int i = 0; i < bulkNSymm; i++) {
            double[] x = coordinates[i][0];
            double[] y = coordinates[i][1];
            double[] z = coordinates[i][2];
            for (int j = 0; j < nAtoms; j++) {
                Atom aj = atoms[j];
                crystal.applySymOp(aj.getXYZ(null), xyz, symops.get(i));
                x[j] = xyz[0];
                y[j] = xyz[1];
                z[j] = xyz[2];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            if (solvent) {
                sb.append("  Bulk Solvent Grid\n");
                sb.append(format("  Bulk solvent model type:           %8s\n", solventModel.toString()));
            } else {
                sb.append("  Atomic Scattering Grid\n");
            }
            sb.append(format("  Form factor grid points:             %8d\n", aRadGrid * 2));
            sb.append(format("  Form factor grid diameter:           %8.3f\n", aRad * 2));
            sb.append(format("  Grid density:                        %8.3f\n", 1.0 / gridStep));
            sb.append(format("  Grid dimensions:                (%3d,%3d,%3d)\n", fftX, fftY, fftZ));
            sb.append(format("  Grid method:                         %8s\n", gridMethod.toString()));
            logger.info(sb.toString());
        }

        if (solvent) {
            int minWork = nSymm;
            if (solventModel != SolventModel.NONE) {
                solventGradientRegion = new SolventGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES);
                // Null out scattering instances used without bulk solvent.
                atomicGradientRegion = null;
                atomicSliceLoops = null;
                atomicRowLoops = null;
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
                        SolventDensityLoop[] solventDensityLoops = new SolventDensityLoop[threadCount];
                        SolventDensityLoop[] bulkSolventDensityLoops = new SolventDensityLoop[threadCount];
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

                        atomicRowRegion = null;
                        solventRowRegion = null;
                        bulkSolventRowLoops = null;
                        solventRowLoops = null;
                        break;
                    case ROW:
                        atomicRowRegion
                                = new RowRegion(fftX, fftY, fftZ, densityGrid, nSymm, threadCount, atoms, coordinates);
                        solventRowRegion
                                = new BulkSolventRowRegion(fftX, fftY, fftZ, solventGrid, bulkNSymm,
                                threadCount, crystal, atoms, coordinates, 4.0, parallelTeam);
                        if (solventModel == SolventModel.GAUSSIAN) {
                            atomicRowRegion.setInitValue(0.0);
                        } else {
                            atomicRowRegion.setInitValue(1.0);
                        }
                        solventRowLoops = new SolventRowLoop[threadCount];
                        bulkSolventRowLoops = new BulkSolventRowLoop[threadCount];
                        for (int i = 0; i < threadCount; i++) {
                            solventRowLoops[i] = new SolventRowLoop(atomicRowRegion);
                            bulkSolventRowLoops[i] = new BulkSolventRowLoop(solventRowRegion);
                        }
                        atomicRowRegion.setDensityLoop(solventRowLoops);
                        solventRowRegion.setDensityLoop(bulkSolventRowLoops);
                        atomicDensityRegion = null;
                        solventDensityRegion = null;
                        atomicSliceRegion = null;
                        solventSliceRegion = null;
                        bulkSolventSliceLoops = null;
                        solventSliceLoops = null;
                        break;
                    case SLICE:
                    default:
                        atomicSliceRegion
                                = new SliceRegion(fftX, fftY, fftZ, densityGrid, nSymm, threadCount, atoms, coordinates);
                        solventSliceRegion
                                = new BulkSolventSliceRegion(fftX, fftY, fftZ, solventGrid, bulkNSymm,
                                threadCount, crystal, atoms, coordinates, 4.0, parallelTeam);
                        if (solventModel == SolventModel.GAUSSIAN) {
                            atomicSliceRegion.setInitValue(0.0);
                        } else {
                            atomicSliceRegion.setInitValue(1.0);
                        }
                        solventSliceLoops = new SolventSliceLoop[threadCount];
                        bulkSolventSliceLoops = new BulkSolventSliceLoop[threadCount];
                        for (int i = 0; i < threadCount; i++) {
                            solventSliceLoops[i] = new SolventSliceLoop(atomicSliceRegion);
                            bulkSolventSliceLoops[i] = new BulkSolventSliceLoop(solventSliceRegion);
                        }
                        atomicSliceRegion.setDensityLoop(solventSliceLoops);
                        solventSliceRegion.setDensityLoop(bulkSolventSliceLoops);

                        atomicDensityRegion = null;
                        solventDensityRegion = null;
                        atomicRowRegion = null;
                        solventRowRegion = null;
                        bulkSolventRowLoops = null;
                        solventRowLoops = null;

                }
            } else {
                atomicDensityRegion = null;
                solventDensityRegion = null;
                atomicGradientRegion = null;
                solventGradientRegion = null;
                atomicSliceRegion = null;
                atomicSliceLoops = null;
                solventSliceRegion = null;
                bulkSolventSliceLoops = null;
                solventSliceLoops = null;
                atomicRowRegion = null;
                atomicRowLoops = null;
                solventRowRegion = null;
                bulkSolventRowLoops = null;
                solventRowLoops = null;
            }
        } else {
            atomicGradientRegion = new AtomicGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES);

            // Null out bulk solvent instances.
            solventDensityRegion = null;
            solventSliceRegion = null;
            solventSliceLoops = null;
            bulkSolventSliceLoops = null;
            solventRowRegion = null;
            solventRowLoops = null;
            bulkSolventRowLoops = null;
            solventGradientRegion = null;

            /*
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
                    AtomicDensityLoop[] atomicDensityLoops = new AtomicDensityLoop[threadCount];
                    for (int i = 0; i < threadCount; i++) {
                        atomicDensityLoops[i]
                                = new AtomicDensityLoop(atomicDensityRegion);
                    }
                    atomicDensityRegion.setDensityLoop(atomicDensityLoops);
                    atomicSliceRegion = null;
                    atomicSliceLoops = null;
                    atomicRowRegion = null;
                    atomicRowLoops = null;
                    break;

                case ROW:
                    atomicRowRegion = new RowRegion(fftX, fftY, fftZ, densityGrid, nSymm, threadCount, atoms, coordinates);
                    atomicRowRegion.setInitValue(0.0);
                    atomicRowLoops = new AtomicRowLoop[threadCount];
                    for (int i = 0; i < threadCount; i++) {
                        atomicRowLoops[i] = new AtomicRowLoop(atomicRowRegion);
                    }
                    atomicRowRegion.setDensityLoop(atomicRowLoops);
                    atomicDensityRegion = null;
                    atomicSliceRegion = null;
                    atomicSliceLoops = null;
                    break;

                case SLICE:
                default:
                    atomicSliceRegion = new SliceRegion(fftX, fftY, fftZ, densityGrid, nSymm, threadCount, atoms, coordinates);
                    atomicSliceRegion.setInitValue(0.0);
                    atomicSliceLoops = new AtomicSliceLoop[threadCount];
                    for (int i = 0; i < threadCount; i++) {
                        atomicSliceLoops[i] = new AtomicSliceLoop(atomicSliceRegion);
                    }
                    atomicSliceRegion.setDensityLoop(atomicSliceLoops);
                    atomicDensityRegion = null;
                    atomicRowRegion = null;
                    atomicRowLoops = null;
            }
        }

        extractRegion = new ExtractRegion(threadCount);
        babinetRegion = new BabinetRegion(threadCount);
        solventGridRegion = new SolventGridRegion(threadCount);
        initRegion = new InitRegion(threadCount);
        solventScaleRegion = new SolventScaleRegion(threadCount);
        atomicScaleRegion = new AtomicScaleRegion(threadCount);
        complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam);
    }

    /**
     * Parse a solvent model string and return a SolventModel enum.
     *
     * @param str The solvent model String.
     * @return The SolventModel to use.
     */
    public static SolventModel parseSolventModel(String str) {
        try {
            return SolventModel.valueOf(str.toUpperCase().replaceAll("\\s+", ""));
        } catch (Exception e) {
            logger.info(format(" Could not parse %s as a Solvent Model; defaulting to Polynomial.", str));
            return SolventModel.POLYNOMIAL;
        }
    }

    /**
     * <p>Setter for the field <code>nativeEnvironmentApproximation</code>.</p>
     *
     * @param nativeEnvironmentApproximation a boolean.
     */
    void setNativeEnvironmentApproximation(boolean nativeEnvironmentApproximation) {
        this.nativeEnvironmentApproximation = nativeEnvironmentApproximation;
    }

    /**
     * <p>Getter for the field <code>densityGrid</code>.</p>
     *
     * @return the densityGrid
     */
    public double[] getDensityGrid() {
        return densityGrid;
    }

    /**
     * <p>Getter for the field <code>solventGrid</code>.</p>
     *
     * @return the solventGrid
     */
    double[] getSolventGrid() {
        return solventGrid;
    }

    /**
     * <p>Setter for the field <code>lambdaTerm</code>.</p>
     *
     * @param lambdaTerm a boolean.
     */
    void setLambdaTerm(boolean lambdaTerm) {
        this.lambdaTerm = lambdaTerm;
    }

    /**
     * set the bulk solvent parameters
     *
     * @param a added atom width (binary, polynomial only)
     * @param b falloff of the switch (Gaussian, polynomial only)
     */
    void setSolventAB(double a, double b) {
        double vdwr;
        this.solventA = a;
        this.solventB = b;
        if (!solvent) {
            return;
        }
        switch (solventModel) {
            case BINARY:
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventFormFactors[iSymm][i] = new SolventBinaryFormFactor(atoms[i], vdwr + solventA);
                    }
                }
                break;
            case GAUSSIAN:
                for (int iSymm = 0; iSymm < bulkNSymm; iSymm++) {
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventFormFactors[iSymm][i] = new SolventGaussFormFactor(atoms[i], vdwr * solventB);
                    }
                }
                break;
            case POLYNOMIAL:
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
    void setUse3G(boolean useThreeGaussians) {
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
     * @param n     {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    void deltaX(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double[] xyz = new double[3];
        double[] symxyz = new double[3];
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
     * @param n     {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    void deltaY(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double[] xyz = new double[3];
        double[] symxyz = new double[3];
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
     * @param n     {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    void deltaZ(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double[] xyz = new double[3];
        double[] symxyz = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[2] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][2][n] = symxyz[2];
        }
    }

    /**
     * Set atomic coordinates.
     *
     * @param coords new coordinate positions (3 params per atom)
     */
    public void setCoordinates(double[] coords) {
        assert (coords != null);
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double[] xyz = new double[3];
        double[] symXYZ = new double[3];
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            if (!atoms[i].isActive()) {
                // Coordinates for inactive atoms do not change.
                continue;
            }
            xyz[0] = coords[index++];
            xyz[1] = coords[index++];
            xyz[2] = coords[index++];
            for (int j = 0; j < nSymm; j++) {
                crystal.applySymOp(xyz, symXYZ, symops.get(j));
                coordinates[j][0][i] = symXYZ[0];
                coordinates[j][1][i] = symXYZ[1];
                coordinates[j][2][i] = symXYZ[2];
            }
        }
    }

    /**
     * parallelized computation of structure factors
     *
     * @param hklData structure factor list to fill in
     * @see DiffractionRefinementData
     */
    void computeDensity(double[][] hklData) {
        computeDensity(hklData, false);
    }

    /**
     * parallelized computation of structure factors
     *
     * @param hklData structure factor list to fill in
     * @param print   if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    void computeDensity(double[][] hklData, boolean print) {
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
     * @param hklData        structure factors to apply inverse FFT
     * @param freer          array of free r flags corresponding to hkldata
     * @param flag           Rfree flag value
     * @param refinementMode {@link RefinementMinimize.RefinementMode refinement mode}
     * @see RefinementMinimize.RefinementMode
     * @see DiffractionRefinementData
     */
    public void computeAtomicGradients(double[][] hklData, int[] freer, int flag,
                                       RefinementMode refinementMode) {
        computeAtomicGradients(hklData, freer, flag, refinementMode, false);
    }

    /**
     * compute inverse FFT to determine atomic gradients
     *
     * @param hklData        structure factors to apply inverse FFT
     * @param freer          array of free r flags corresponding to hkldata
     * @param flag           Rfree flag value
     * @param refinementMode {@link RefinementMinimize.RefinementMode refinement mode}
     * @param print          if true, print information on timings during the calculation
     * @see RefinementMinimize.RefinementMode
     * @see DiffractionRefinementData
     */
    void computeAtomicGradients(double[][] hklData, int[] freer, int flag,
                                RefinementMode refinementMode, boolean print) {

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
            double[] fc = hklData[ih.index()];
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
            c.timesIP(2.0 / fftScale);

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
                    cj.phaseShiftIP(shift);
                    densityGrid[ii] += cj.re();
                    densityGrid[ii + 1] += -cj.im();
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    cj.phaseShiftIP(shift);
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
                for (int i = 0; i < nAtoms; i++) {
                    optAtomicGradientWeight.set(i, 0);
                }
                atomicGradientSchedule.updateWeights(previousOptAtomicGradientWeight);
                atomicGradientRegion.setRefinementMode(refinementMode);
                parallelTeam.execute(atomicGradientRegion);
                for (int i = 0; i < nAtoms; i++) {
                    previousOptAtomicGradientWeight[i] = optAtomicGradientWeight.get(i);
                }
            }
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Exception computing atomic gradients.";
            logger.log(Level.SEVERE, message, e);
        }

        if (solvent) {
            sb.append(format(" Solvent symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(format(" Solvent inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(format(" Grid solvent gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(format(" %d reflections ignored (cross validation set)\n", nfree));
        } else {
            sb.append(format(" Atomic symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(format(" Atomic inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(format(" Grid atomic gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(format(" %d reflections ignored (cross validation set)\n", nfree));
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
    void computeAtomicDensity(double[][] hklData) {
        computeAtomicDensity(hklData, false);
    }

    /**
     * parallelized computation of structure factors (identical to
     * computeDensity)
     *
     * @param hklData structure factor list to fill in
     * @param print   if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    private void computeAtomicDensity(double[][] hklData, boolean print) {
        // Zero out reflection data.
        long initTime = -System.nanoTime();
        try {
            initRegion.setHKL(hklData);
            parallelTeam.execute(initRegion);
        } catch (Exception e) {
            String message = "Fatal exception zeroing out reflection data.";
            logger.log(Level.SEVERE, message, e);
        }
        initTime += System.nanoTime();

        // Assign atomic electron density to FFT grid.
        long atomicGridTime = -System.nanoTime();
        try {
            switch (gridMethod) {
                case SPATIAL:
                    atomicDensityRegion.assignAtomsToCells();
                    parallelTeam.execute(atomicDensityRegion);
                    break;

                case ROW:
                    for (int i = 0; i < fftZ * fftY; i++) {
                        optRowWeightAtomic.set(i, 0);
                    }
                    atomicRowSchedule.updateWeights(previousOptRowWeightAtomic);
                    parallelTeam.execute(atomicRowRegion);
                    for (int i = 0; i < fftZ * fftY; i++) {
                        previousOptRowWeightAtomic[i] = optRowWeightAtomic.get(i);
                    }
                    break;
                case SLICE:
                default:
                    for (int i = 0; i < fftZ; i++) {
                        optSliceWeightAtomic.set(i, 0);
                    }
                    atomicSliceSchedule.updateWeights(previousOptSliceWeightAtomic);
                    parallelTeam.execute(atomicSliceRegion);
                    for (int i = 0; i < fftZ; i++) {
                        previousOptSliceWeightAtomic[i] = optSliceWeightAtomic.get(i);
                    }
                    break;
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating atomic electron density.";
            logger.log(Level.SEVERE, message, e);
        }
        atomicGridTime += System.nanoTime();

        // Compute model structure factors via an FFT of the electron density.
        long startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        // Extract and scale structure factors.
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
            sb.append(format("\n Fc Initialization:         %8.4f\n", initTime * toSeconds));
            sb.append(format(" Atomic grid density:       %8.4f\n", atomicGridTime * toSeconds));
            sb.append(format(" Atomic FFT:                %8.4f\n", fftTime * toSeconds));
            sb.append(format(" Atomic symmetry & scaling: %8.4f\n", symTime * toSeconds));
            logger.info(sb.toString());

            StringBuilder ASB = new StringBuilder();
            switch (gridMethod) {

                case ROW:
                    double atomicRowTotal = 0;
                    double atomicRowMin = atomicRowLoops[0].getThreadTime();
                    double atomicRowMax = 0;
                    double atomicRowWeightTotal = 0;

                    for (int i = 0; i < threadCount; i++) {
                        atomicRowTotal += atomicRowLoops[i].getThreadTime();
                        atomicRowMax = max(atomicRowLoops[i].getThreadTime(), atomicRowMax);
                        atomicRowMin = min(atomicRowLoops[i].getThreadTime(), atomicRowMin);
                        atomicRowWeightTotal += atomicRowLoops[i].getThreadWeight();
                    }

                    //Atomic timing and balance analysis
                    ASB.append(format("\n RowLoop (Atomic): %7.5f (sec)                 | Total Weight: %7.0f\n",
                            atomicGridTime * toSeconds, atomicRowWeightTotal));
                    ASB.append(" Thread     LoopTime    Balance(%)  Normalized   |      Rows       Weight    Balance(%)  Normalized\n");

                    //check for no weights issued, then set to 1 so a 0 is printed instead of NaN
                    if (atomicRowWeightTotal == 0) {
                        atomicRowWeightTotal = 1;
                    }

                    for (int i = 0; i < threadCount; i++) {
                        ASB.append(format("     %3d     %7.5f     %7.1f     %7.1f     |   %7d     %7d     %7.1f     %7.1f\n", i,
                                (double) (atomicRowLoops[i].getThreadTime() * toSeconds),
                                ((atomicRowLoops[i].getThreadTime()) / (atomicRowTotal)) * 100.00,
                                ((atomicRowLoops[i].getThreadTime()) / (atomicRowTotal)) * (100.00 * threadCount),
                                atomicRowLoops[i].getNumberofSlices(), atomicRowLoops[i].getThreadWeight(),
                                100.00 * (atomicRowLoops[i].getThreadWeight() / atomicRowWeightTotal),
                                (100.00 * threadCount) * (atomicRowLoops[i].getThreadWeight() / atomicRowWeightTotal)));
                    }
                    ASB.append(format("    Min      %7.5f\n", atomicRowMin * toSeconds));
                    ASB.append(format("    Max      %7.5f\n", atomicRowMax * toSeconds));
                    ASB.append(format("    Delta    %7.5f\n", (atomicRowMax - atomicRowMin) * toSeconds));
                    logger.info(ASB.toString());

                    break;
                case SPATIAL:
                    break;
                case SLICE:
                    double atomicSliceTotal = 0;
                    double atomicSliceMin = atomicSliceLoops[0].getThreadTime();
                    double atomicSliceMax = 0;
                    double atomicSliceWeightTotal = 0;

                    for (int i = 0; i < threadCount; i++) {
                        atomicSliceTotal += atomicSliceLoops[i].getThreadTime();
                        atomicSliceMax = max(atomicSliceLoops[i].getThreadTime(), atomicSliceMax);
                        atomicSliceMin = min(atomicSliceLoops[i].getThreadTime(), atomicSliceMin);
                        atomicSliceWeightTotal += atomicSliceLoops[i].getThreadWeight();
                    }

                    //Atomic timing and balance analysis
                    ASB.append(format("\n SliceLoop (Atomic): %7.5f (sec)               | Total Weight: %7.0f\n",
                            atomicGridTime * toSeconds, atomicSliceWeightTotal));
                    ASB.append(" Thread     LoopTime    Balance(%)  Normalized   |      Slices    Weight    Balance(%)  Normalized\n");

                    //check for no weights issued, then set to 1 so a 0 is printed instead of NaN
                    if (atomicSliceWeightTotal == 0) {
                        atomicSliceWeightTotal = 1;
                    }

                    for (int i = 0; i < threadCount; i++) {
                        ASB.append(format("     %3d     %7.5f     %7.1f     %7.1f     |   %7d     %7d     %7.1f     %7.1f\n", i,
                                (double) (atomicSliceLoops[i].getThreadTime() * toSeconds),
                                ((atomicSliceLoops[i].getThreadTime()) / (atomicSliceTotal)) * 100.00,
                                ((atomicSliceLoops[i].getThreadTime()) / (atomicSliceTotal)) * (100.00 * threadCount),
                                atomicSliceLoops[i].getNumberofSlices(),
                                atomicSliceLoops[i].getThreadWeight(),
                                100.00 * (atomicSliceLoops[i].getThreadWeight() / atomicSliceWeightTotal),
                                (100.00 * threadCount) * (atomicSliceLoops[i].getThreadWeight() / atomicSliceWeightTotal)));
                    }
                    ASB.append(format("    Min      %7.5f\n", atomicSliceMin * toSeconds));
                    ASB.append(format("    Max      %7.5f\n", atomicSliceMax * toSeconds));
                    ASB.append(format("    Delta    %7.5f\n", (atomicSliceMax - atomicSliceMin) * toSeconds));
                    logger.info(ASB.toString());
                default:
            }
        }

    }

    /**
     * parallelized computation of bulk solvent structure factors
     *
     * @param hklData structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeSolventDensity(double[][] hklData) {
        computeSolventDensity(hklData, false);
    }

    /**
     * parallelized computation of bulk solvent structure factors
     *
     * @param hklData structure factor list to fill in
     * @param print   if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    private void computeSolventDensity(double[][] hklData, boolean print) {
        // Zero out the reflection data.
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
                case ROW:
                    for (int i = 0; i < fftZ * fftY; i++) {
                        optRowWeightSolvent.set(i, 0);
                    }
                    solventRowSchedule.updateWeights(previousOptRowWeightSolvent);
                    parallelTeam.execute(atomicRowRegion);
                    for (int i = 0; i < fftZ * fftY; i++) {
                        previousOptRowWeightSolvent[i] = optRowWeightSolvent.get(i);
                    }
                    break;
                case SLICE:
                default:
                    for (int i = 0; i < fftZ; i++) {
                        optSliceWeightSolvent.set(i, 0);
                    }
                    solventSliceSchedule.updateWeights(previousOptSliceWeightSolvent);
                    parallelTeam.execute(atomicSliceRegion);
                    for (int i = 0; i < fftZ; i++) {
                        previousOptSliceWeightSolvent[i] = optSliceWeightSolvent.get(i);
                    }
                    break;
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }
        solventGridTime += System.nanoTime();

        long expTime = -System.nanoTime();
        if (solventModel == SolventModel.BINARY) {
            /*
             * Need to shrink mask. First, copy densityGrid to solventGrid
             * temporarily to store original map.
             */
            int nmap = densityGrid.length;
            arraycopy(densityGrid, 0, solventGrid, 0, nmap);
            // Logic to loop within the cutoff box.
            double[] xyz = {solventB, solventB, solventB};
            double[] uvw = new double[3];
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
            arraycopy(solventGrid, 0, densityGrid, 0, nmap);
        }

        // Copy the grid over for derivatives.
        ArrayCopyRegion arrayCopyRegion = new ArrayCopyRegion();
        try {
            parallelTeam.execute(arrayCopyRegion);
        } catch (Exception e) {
            logger.info(e.toString());
        }

        // Babinet Principle.
        try {
            parallelTeam.execute(babinetRegion);
        } catch (Exception e) {
            String message = "Fatal exception Babinet Principle.";
            logger.log(Level.SEVERE, message, e);
        }
        expTime += System.nanoTime();

        long fftTime = -System.nanoTime();
        complexFFT3D.fft(densityGrid);
        fftTime += System.nanoTime();

        // Extract and scale structure factors.
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

        // Expand derivative mask so it includes nearby local density.
        long solventDensityTime = -System.nanoTime();
        try {
            switch (gridMethod) {
                case SPATIAL:
                    solventDensityRegion.assignAtomsToCells();
                    parallelTeam.execute(solventDensityRegion);
                    break;

                case ROW:
                    for (int i = 0; i < (fftZ * fftY); i++) {
                        optRowWeightBulkSolvent.set(i, 0);
                    }
                    bulkSolventRowSchedule.updateWeights(previousOptRowWeightBulkSolvent);
                    parallelTeam.execute(solventRowRegion);
                    for (int i = 0; i < fftZ * fftY; i++) {
                        previousOptRowWeightBulkSolvent[i] = optRowWeightBulkSolvent.get(i);
                    }
                    break;
                case SLICE:
                default:
                    for (int i = 0; i < (fftZ); i++) {
                        optSliceWeightBulkSolvent.set(i, 0);
                    }

                    bulkSolventSliceSchedule.updateWeights(previousOptSliceWeightBulkSolvent);
                    parallelTeam.execute(solventSliceRegion);
                    for (int i = 0; i < fftZ; i++) {
                        previousOptSliceWeightBulkSolvent[i] = optSliceWeightBulkSolvent.get(i);
                    }
                    break;
            }
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }
        solventDensityTime += System.nanoTime();

        long solventExpTime = 0;
        if (solventModel == SolventModel.GAUSSIAN) {
            solventExpTime = -System.nanoTime();
            try {
                parallelTeam.execute(solventGridRegion);
            } catch (Exception e) {
                String message = "Fatal exception solvent grid region.";
                logger.log(Level.SEVERE, message, e);
            }
            solventExpTime += System.nanoTime();
        }

        if (logger.isLoggable(Level.INFO) && print) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n Fc Initialization:                     %8.4f\n", initTime * toSeconds));
            sb.append(format(" Solvent grid density:                  %8.4f\n", solventGridTime * toSeconds));
            sb.append(format(" Bulk solvent exponentiation:           %8.4f\n", expTime * toSeconds));
            sb.append(format(" Solvent FFT:                           %8.4f\n", fftTime * toSeconds));
            sb.append(format(" Solvent symmetry & scaling:            %8.4f\n", symTime * toSeconds));
            sb.append(format(" Solvent grid expansion:                %8.4f\n", solventDensityTime * toSeconds));
            if (solventModel == SolventModel.GAUSSIAN) {
                sb.append(format(" Gaussian grid exponentiation: %8.4f\n",
                        solventExpTime * toSeconds));
            }
            logger.info(sb.toString());

            StringBuilder SSB = new StringBuilder();
            StringBuilder BSSB = new StringBuilder();
            switch (gridMethod) {
                case ROW:
                    double solventRowTotal = 0;
                    double solventRowMin = solventRowLoops[0].getThreadTime();
                    double solventRowMax = 0;
                    double solventRowWeightTotal = 0;

                    double bulkSolventRowTotal = 0;
                    double bulkSolventRowMin = bulkSolventRowLoops[0].getTimePerThread();
                    double bulkSolventRowMax = 0;
                    double bulkSolventRowWeightTotal = 0;

                    for (int i = 0; i < threadCount; i++) {
                        solventRowTotal += solventRowLoops[i].getThreadTime();
                        solventRowMax = max(solventRowLoops[i].getThreadTime(), solventRowMax);
                        solventRowMin = min(solventRowLoops[i].getThreadTime(), solventRowMin);
                        solventRowWeightTotal += solventRowLoops[i].getThreadWeight();

                        bulkSolventRowTotal += bulkSolventRowLoops[i].getTimePerThread();
                        bulkSolventRowMax = max(bulkSolventRowLoops[i].getTimePerThread(), bulkSolventRowMax);
                        bulkSolventRowMin = min(bulkSolventRowLoops[i].getTimePerThread(), bulkSolventRowMin);
                        bulkSolventRowWeightTotal += bulkSolventRowLoops[i].getWeightPerThread()[i];

                    }

                    // Solvent timing and balance analysis
                    SSB.append(format("\n RowLoop (Solvent): %7.5f (sec)                | Total Weight: %7.0f\n", solventGridTime * toSeconds, solventRowWeightTotal));
                    SSB.append(" Thread     LoopTime    Balance(%)  Normalized   |      Rows       Weight     Balance(%)  Normalized\n");

                    // Check for no weights issued, then set to 1 so a 0 is printed instead of NaN
                    if (solventRowWeightTotal == 0) {
                        solventRowWeightTotal = 1;
                    }

                    for (int i = 0; i < threadCount; i++) {
                        SSB.append(format("     %3d     %7.5f     %7.1f     %7.1f     |   %7d     %7d     %7.1f     %7.1f\n", i, (double) (solventRowLoops[i].getThreadTime() * toSeconds),
                                ((solventRowLoops[i].getThreadTime()) / (solventRowTotal)) * 100.00, ((solventRowLoops[i].getThreadTime()) / (solventRowTotal)) * (100.00 * threadCount),
                                solventRowLoops[i].getNumberofSlices(), solventRowLoops[i].getThreadWeight(), 100.00 * (solventRowLoops[i].getThreadWeight() / solventRowWeightTotal),
                                (100.00 * threadCount) * (solventRowLoops[i].getThreadWeight() / solventRowWeightTotal)));
                    }
                    SSB.append(format("    Min      %7.5f\n", solventRowMin * toSeconds));
                    SSB.append(format("    Max      %7.5f\n", solventRowMax * toSeconds));
                    SSB.append(format("    Delta    %7.5f\n", (solventRowMax - solventRowMin) * toSeconds));
                    logger.info(SSB.toString());

                    // Bulk solvent timing and balance analysis
                    BSSB.append(format("\n RowLoop (Bulk Solvent): %7.5f (sec)           | Total Weight: %7.0f\n", bulkSolventRowTotal * toSeconds, bulkSolventRowWeightTotal));
                    BSSB.append(" Thread     LoopTime    Balance(%)  Normalized   |      Rows       Weight     Balance(%)  Normalized\n");

                    //check for no weights issued, then set to 1 so a 0 is printed instead of NaN
                    if (bulkSolventRowWeightTotal == 0) {
                        bulkSolventRowWeightTotal = 1;
                    }

                    for (int i = 0; i < threadCount; i++) {
                        BSSB.append(format("     %3d     %7.5f     %7.1f     %7.1f     |   %7d     %7d     %7.1f     %7.1f\n", i, (double) (bulkSolventRowLoops[i].getTimePerThread() * toSeconds),
                                ((bulkSolventRowLoops[i].getTimePerThread()) / (bulkSolventRowTotal)) * 100.00, ((bulkSolventRowLoops[i].getTimePerThread()) / (bulkSolventRowTotal)) * (100.00 * threadCount),
                                bulkSolventRowLoops[i].getBoundsPerThread()[i], bulkSolventRowLoops[i].getWeightPerThread()[i], 100.00 * (bulkSolventRowLoops[i].getWeightPerThread()[i] / bulkSolventRowWeightTotal),
                                (100.00 * threadCount) * (bulkSolventRowLoops[i].getWeightPerThread()[i] / bulkSolventRowWeightTotal)));
                    }
                    BSSB.append(format("    Min      %7.5f\n", bulkSolventRowMin * toSeconds));
                    BSSB.append(format("    Max      %7.5f\n", bulkSolventRowMax * toSeconds));
                    BSSB.append(format("    Delta    %7.5f\n", (bulkSolventRowMax - bulkSolventRowMin) * toSeconds));
                    logger.info(BSSB.toString());
                    break;
                case SPATIAL:
                    break;
                case SLICE:
                    double solventSliceTotal = 0;
                    double solventSliceMin = solventSliceLoops[0].getThreadTime();
                    double solventSliceMax = 0;
                    double solventSliceWeightTotal = 0;

                    double bulkSolventSliceTotal = 0;
                    double bulkSolventSliceMin = bulkSolventSliceLoops[0].getTimePerThread();
                    double bulkSolventSliceMax = 0;
                    double bulkSolventSliceWeightTotal = 0;

                    for (int i = 0; i < threadCount; i++) {
                        solventSliceTotal += solventSliceLoops[i].getThreadTime();
                        solventSliceMax = max(solventSliceLoops[i].getThreadTime(), solventSliceMax);
                        solventSliceMin = min(solventSliceLoops[i].getThreadTime(), solventSliceMin);
                        solventSliceWeightTotal += solventSliceLoops[i].getThreadWeight();

                        bulkSolventSliceTotal += bulkSolventSliceLoops[i].getTimePerThread();
                        bulkSolventSliceMax = max(bulkSolventSliceLoops[i].getTimePerThread(), bulkSolventSliceMax);
                        bulkSolventSliceMin = min(bulkSolventSliceLoops[i].getTimePerThread(), bulkSolventSliceMin);
                        bulkSolventSliceWeightTotal += bulkSolventSliceLoops[i].getWeightPerThread()[i];

                    }

                    // Solvent timing and balance analysis
                    SSB.append(format("\n SliceLoop (Solvent): %7.5f (sec)               | Total Weight: %7.0f\n", solventGridTime * toSeconds, solventSliceWeightTotal));
                    SSB.append(" Thread     LoopTime    Balance(%)  Normalized   |      Slices    Weight     Balance(%)  Normalized\n");

                    // Check for no weights issued, then set to 1 so a 0 is printed instead of NaN
                    if (solventSliceWeightTotal == 0) {
                        solventSliceWeightTotal = 1;
                    }

                    for (int i = 0; i < threadCount; i++) {
                        SSB.append(format("     %3d     %7.5f     %7.1f     %7.1f     |   %7d     %7d     %7.1f     %7.1f\n", i,
                                (double) (solventSliceLoops[i].getThreadTime() * toSeconds),
                                ((solventSliceLoops[i].getThreadTime()) / (solventSliceTotal)) * 100.00,
                                ((solventSliceLoops[i].getThreadTime()) / (solventSliceTotal)) * (100.00 * threadCount),
                                solventSliceLoops[i].getNumberofSlices(), solventSliceLoops[i].getThreadWeight(),
                                100.00 * (solventSliceLoops[i].getThreadWeight() / solventSliceWeightTotal),
                                (100.00 * threadCount) * (solventSliceLoops[i].getThreadWeight() / solventSliceWeightTotal)));
                    }
                    SSB.append(format("    Min      %7.5f\n", solventSliceMin * toSeconds));
                    SSB.append(format("    Max      %7.5f\n", solventSliceMax * toSeconds));
                    SSB.append(format("    Delta    %7.5f\n", (solventSliceMax - solventSliceMin) * toSeconds));
                    logger.info(SSB.toString());

                    // Bulk solvent timing and balance analysis
                    BSSB.append(format("\n SliceLoop (Bulk Solvent): %7.5f (sec)          | Total Weight: %7.0f\n", bulkSolventSliceTotal * toSeconds, bulkSolventSliceWeightTotal));
                    BSSB.append(" Thread     LoopTime    Balance(%)  Normalized   |      Slices     Weight     Balance(%)  Normalized\n");

                    // Check for no weights issued, then set to 1 so a 0 is printed instead of NaN
                    if (bulkSolventSliceWeightTotal == 0) {
                        bulkSolventSliceWeightTotal = 1;
                    }

                    for (int i = 0; i < threadCount; i++) {
                        BSSB.append(format("     %3d     %7.5f     %7.1f     %7.1f     |   %7d     %7d     %7.1f     %7.1f\n", i, (double) (bulkSolventSliceLoops[i].getTimePerThread() * toSeconds),
                                ((bulkSolventSliceLoops[i].getTimePerThread()) / (bulkSolventSliceTotal)) * 100.00, ((bulkSolventSliceLoops[i].getTimePerThread()) / (bulkSolventSliceTotal)) * (100.00 * threadCount),
                                bulkSolventSliceLoops[i].getBoundsPerThread()[i], bulkSolventSliceLoops[i].getWeightPerThread()[i], 100.00 * (bulkSolventSliceLoops[i].getWeightPerThread()[i] / bulkSolventSliceWeightTotal),
                                (100.00 * threadCount) * (bulkSolventSliceLoops[i].getWeightPerThread()[i] / bulkSolventSliceWeightTotal)));
                    }
                    BSSB.append(format("    Min      %7.5f\n", bulkSolventSliceMin * toSeconds));
                    BSSB.append(format("    Max      %7.5f\n", bulkSolventSliceMax * toSeconds));
                    BSSB.append(format("    Delta    %7.5f\n", (bulkSolventSliceMax - bulkSolventSliceMin) * toSeconds));
                    logger.info(BSSB.toString());
                    break;
                default:
            }
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
     * @param data   an array of double.
     * @param meansd an array of double.
     * @param norm   a boolean.
     */
    public void densityNorm(double[] data, double[] meansd, boolean norm) {
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
                    sd += pow(data[index] - mean, 2);
                    n++;
                }
            }
        }
        sd = sqrt(sd / n);

        if (meansd != null) {
            meansd[0] = mean;
            meansd[1] = sd;
        }

        if (norm) {
            double isd = 1.0 / sd;
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        int index = iComplex3D(i, j, k, fftX, fftY);
                        data[index] = (data[index] - mean) * isd;
                    }
                }
            }
        }
    }

    /**
     * The possible solvent model methods
     */
    public enum SolventModel {

        /**
         * Do not model solvent scattering.
         */
        NONE,
        /**
         * The classic binary (0, 1) model.
         */
        BINARY,
        /**
         * Smooth the boundary of the classic model using Gaussians.
         */
        GAUSSIAN,
        /**
         * Smooth the boundar of the classic model using a cubic polynomial
         * switch (default).
         */
        POLYNOMIAL
    }

    public enum GridMethod {

        SPATIAL, SLICE, ROW
    }

    private class AtomicDensityLoop extends SpatialDensityLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;

        AtomicDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            if (!atoms[n].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[n].applyLambda()) {
                return;
            }

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
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = atomff.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class SolventDensityLoop extends SpatialDensityLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;

        SolventDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            if (!atoms[n].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[n].applyLambda()) {
                return;
            }

            xyz[0] = coordinates[iSymm][0][n];
            xyz[1] = coordinates[iSymm][1][n];
            xyz[2] = coordinates[iSymm][2][n];
            FormFactor solventff = solventFormFactors[iSymm][n];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[n].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case BINARY:
                    frad = Math.min(aRadGrid,
                            (int) Math.floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case GAUSSIAN:
                    frad = Math.min(aRadGrid,
                            (int) Math.floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case POLYNOMIAL:
                    frad = Math.min(aRadGrid,
                            (int) Math.floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
                case NONE:
                default:
                    return;
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
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = solventff.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class AtomicRowLoop extends RowLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;
        final int[] optLocal;
        long timer;
        int previousUB, previousLB;
        int actualWeight;

        AtomicRowLoop(RowRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
            optLocal = new int[fftZ * fftY];
        }

        double getThreadTime() {
            return timer;
        }

        int getNumberofSlices() {
            return (previousUB - previousLB + 1);
        }

        int getThreadWeight() {
            return actualWeight;
        }

        @Override
        public IntegerSchedule schedule() {
            return atomicRowSchedule;
        }

        @Override
        public void start() {
            fill(optLocal, 0);
            actualWeight = 0;
            timer = -System.nanoTime();
        }

        @Override
        public void finish() {
            for (int i = 0; i < fftZ * fftY; i++) {
                optRowWeightAtomic.addAndGet(i, optLocal[i]);
            }
            timer += System.nanoTime();
        }

        @Override
        public boolean checkList(int[][][] zyAtListBuild, int buff) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    if (rowRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;
                        final int previousZ = zyAtListBuild[iSymm][iAtom][0];

                        final double fry = fftY * uvw[1];
                        final int ifry = (int) fry;
                        final int previousY = zyAtListBuild[iSymm][iAtom][1];
                        if (abs(ifrz - previousZ) >= buff || abs(ifry - previousY) >= buff) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        @Override
        public void saveZYValues(int[][][] zyAtListBuild) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    zyAtListBuild[iSymm][iAtom][0] = 0;
                    zyAtListBuild[iSymm][iAtom][1] = 0;
                    if (rowRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;

                        final double fry = fftY * uvw[1];
                        final int ifry = (int) fry;
                        zyAtListBuild[iSymm][iAtom][0] = ifrz;
                        zyAtListBuild[iSymm][iAtom][1] = ifry;
                    }
                }
            }
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            boolean boundsChange = false;
            if (previousLB != lb || previousUB != ub) {
                boundsChange = true;
            }
            previousLB = lb;
            previousUB = ub;
            if (rebuildList || boundsChange) {
                buildListA = new ArrayList<>();
                buildListS = new ArrayList<>();
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                        if (rowRegion.select[iSymm][iAtom]) {
                            buildList(iSymm, iAtom, lb, ub);
                        }
                    }
                }
            }
            for (int i = 0; i < buildListA.size(); i++) {
                if (rowRegion.select[buildListS.get(i)][buildListA.get(i)]) {
                    gridDensity(buildListS.get(i), buildListA.get(i), lb, ub);
                }
            }
        }

        public void buildList(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }
            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = min(aRadGrid, (int) floor(atoms[iAtom].getFormFactorWidth() * fftX / crystal.a) + 1);

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;
            final int ifrzl = ifrz - frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;
            final int ifryl = ifry - frad;

            // Loop over allowed z coordinates for this Loop
            // Check if the current atom is close enough
            // If so, add to list.
            int buff = bufferSize;

            int lbZ = rowRegion.zFromRowIndex(lb);
            int ubZ = rowRegion.zFromRowIndex(ub);

            for (int iz = ifrzl - buff; iz <= ifrzu + buff; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                if (lbZ > giz || giz > ubZ) {
                    continue;
                }
                int rowLB = rowRegion.rowIndexForYZ(Crystal.mod(ifryl - buff, fftY), giz);
                int rowUB = rowRegion.rowIndexForYZ(Crystal.mod(ifryu + buff, fftY), giz);
                if (lb >= rowLB || rowUB <= ub) {
                    buildListA.add(iAtom);
                    buildListS.add(iSymm);
                    break;
                }
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[iAtom].applyLambda()) {
                return;
            }

            int lbZ = rowRegion.zFromRowIndex(lb);
            int ubZ = rowRegion.zFromRowIndex(ub);

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
                if (lbZ > giz || giz > ubZ) {
                    continue;
                }
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    int rowIndex = rowRegion.rowIndexForYZ(giy, giz);
                    if (lb > rowIndex || rowIndex > ub) {
                        continue;
                    }
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        optLocal[rowIndex]++;
                        actualWeight++;
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = atomff.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class SolventRowLoop extends RowLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;
        final int[] optLocal;
        long timer;
        int previousUB, previousLB;
        int actualWeight;

        public SolventRowLoop(RowRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
            optLocal = new int[fftZ * fftY];
        }

        public double getThreadTime() {
            return timer;
        }

        public int getNumberofSlices() {
            return (previousUB - previousLB + 1);
        }

        public int getThreadWeight() {
            return actualWeight;
        }

        @Override
        public boolean checkList(int[][][] zyAtListBuild, int buff) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    if (rowRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;
                        final int previousZ = zyAtListBuild[iSymm][iAtom][0];

                        final double fry = fftY * uvw[1];
                        final int ifry = (int) fry;
                        final int previousY = zyAtListBuild[iSymm][iAtom][1];
                        if (abs(ifrz - previousZ) >= buff || abs(ifry - previousY) >= buff) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        @Override
        public void saveZYValues(int[][][] zyAtListBuild) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    zyAtListBuild[iSymm][iAtom][0] = 0;
                    zyAtListBuild[iSymm][iAtom][1] = 0;
                    if (rowRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;

                        final double fry = fftY * uvw[1];
                        final int ifry = (int) fry;
                        zyAtListBuild[iSymm][iAtom][0] = ifrz;
                        zyAtListBuild[iSymm][iAtom][1] = ifry;
                    }
                }
            }
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            boolean boundsChange = false;
            if (previousLB != lb || previousUB != ub) {
                boundsChange = true;
            }
            previousLB = lb;
            previousUB = ub;
            if (rebuildList || boundsChange) {
                buildListA = new ArrayList<>();
                buildListS = new ArrayList<>();
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                        if (rowRegion.select[iSymm][iAtom]) {
                            buildList(iSymm, iAtom, lb, ub);
                        }
                    }
                }
            }
            for (int i = 0; i < buildListA.size(); i++) {
                if (rowRegion.select[buildListS.get(i)][buildListA.get(i)]) {
                    gridDensity(buildListS.get(i), buildListA.get(i), lb, ub);
                }
            }
        }

        public void buildList(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }
            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = min(aRadGrid, (int) floor(atoms[iAtom].getFormFactorWidth() * fftX / crystal.a) + 1);

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;
            final int ifrzl = ifrz - frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;
            final int ifryl = ifry - frad;

            // Loop over allowed z coordinates for this Loop
            // Check if the current atom is close enough
            // If so, add to list.
            int buff = bufferSize;

            int lbZ = rowRegion.zFromRowIndex(lb);
            int ubZ = rowRegion.zFromRowIndex(ub);

            for (int iz = ifrzl - buff; iz <= ifrzu + buff; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                if (lbZ > giz || giz > ubZ) {
                    continue;
                }
                int rowLB = rowRegion.rowIndexForYZ(Crystal.mod(ifryl - buff, fftY), giz);
                int rowUB = rowRegion.rowIndexForYZ(Crystal.mod(ifryu + buff, fftY), giz);
                if (lb >= rowLB || rowUB <= ub) {
                    buildListA.add(iAtom);
                    buildListS.add(iSymm);
                    break;
                }
            }
        }

        @Override
        public IntegerSchedule schedule() {
            return solventRowSchedule;
        }

        @Override
        public void start() {
            fill(optLocal, 0);
            timer = -System.nanoTime();
        }

        @Override
        public void finish() {
            timer += System.nanoTime();
            for (int i = 0; i < fftZ * fftY; i++) {
                optRowWeightSolvent.addAndGet(i, optLocal[i]);
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {

            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[iAtom].applyLambda()) {
                return;
            }

            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            FormFactor formFactor = solventFormFactors[iSymm][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[iAtom].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case BINARY:
                    frad = min(aRadGrid, (int) floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case GAUSSIAN:
                    frad = min(aRadGrid, (int) floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case POLYNOMIAL:
                    frad = min(aRadGrid, (int) floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
                case NONE:
                default:
                    return;
            }

            int ubZ = rowRegion.zFromRowIndex(ub);
            int lbZ = rowRegion.zFromRowIndex(lb);

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
                if (lbZ > giz || giz > ubZ) {
                    continue;
                }
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    int rowIndex = rowRegion.rowIndexForYZ(giy, giz);
                    if (lb > rowIndex || rowIndex > ub) {
                        continue;
                    }
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        optLocal[rowIndex]++;
                        actualWeight++;
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = formFactor.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class BulkSolventRowLoop extends RowLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;
        final int[] optLocal;
        long timer;
        int[] threadWeights;
        int[] threadBounds;

        public BulkSolventRowLoop(RowRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
            optLocal = new int[fftZ * fftY];
        }

        public double getTimePerThread() {
            return timer;
        }

        public int[] getBoundsPerThread() {
            return threadBounds;
        }

        public int[] getWeightPerThread() {
            return optLocal;
        }

        @Override
        public IntegerSchedule schedule() {
            return bulkSolventRowSchedule;
        }

        @Override
        public void start() {
            timer = -System.nanoTime();
            fill(optLocal, 0);
        }

        @Override
        public void finish() {
            for (int i = 0; i < fftZ * fftY; i++) {
                optRowWeightBulkSolvent.addAndGet(i, optLocal[i]);
            }
            timer += System.nanoTime();
            threadBounds = bulkSolventRowSchedule.getLowerBounds().clone();
            threadWeights = bulkSolventRowSchedule.getThreadWeights().clone();
            for (int i = threadCount - 1; i > 0; i--) {
                threadWeights[i] -= threadWeights[i - 1];
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[iAtom].applyLambda()) {
                return;
            }

            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            FormFactor formFactor = solventFormFactors[iSymm][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[iAtom].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case BINARY:
                    frad = min(aRadGrid, (int) floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case GAUSSIAN:
                    frad = min(aRadGrid, (int) floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case POLYNOMIAL:
                    frad = min(aRadGrid, (int) floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
                case NONE:
                default:
                    return;
            }

            int ubZ = rowRegion.zFromRowIndex(ub);
            int lbZ = rowRegion.zFromRowIndex(lb);

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
                if (lbZ > giz || giz > ubZ) {
                    continue;
                }
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    int rowIndex = rowRegion.rowIndexForYZ(giy, giz);
                    if (lb > rowIndex || rowIndex > ub) {
                        continue;
                    }
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        optLocal[rowIndex]++;
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = formFactor.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class AtomicSliceLoop extends SliceLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;
        final int[] optLocal;
        long timer;
        int previousUB, previousLB;
        int actualWeight;

        public AtomicSliceLoop(SliceRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
            optLocal = new int[fftZ];
        }

        public double getThreadTime() {
            return timer;
        }

        public int getNumberofSlices() {
            return (previousUB - previousLB + 1);
        }

        public int getThreadWeight() {
            return actualWeight;
        }

        @Override
        public IntegerSchedule schedule() {
            return atomicSliceSchedule;
        }

        @Override
        public void start() {
            fill(optLocal, 0);
            actualWeight = 0;
            timer = -System.nanoTime();
        }

        @Override
        public void finish() {
            for (int i = 0; i < fftZ; i++) {
                optSliceWeightAtomic.addAndGet(i, optLocal[i]);
            }
            timer += System.nanoTime();
        }

        @Override
        public boolean checkList(int[][] zAtListBuild, int buff) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    if (sliceRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;
                        final int previousZ = zAtListBuild[iSymm][iAtom];
                        if (abs(ifrz - previousZ) >= buff) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        @Override
        public void saveZValues(int[][] zAtListBuild) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    zAtListBuild[iSymm][iAtom] = 0;
                    if (sliceRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;
                        zAtListBuild[iSymm][iAtom] = ifrz;
                    }
                }
            }
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            boolean boundsChange = false;
            if (previousLB != lb || previousUB != ub) {
                boundsChange = true;
            }
            previousLB = lb;
            previousUB = ub;
            if (rebuildList || boundsChange) {
                buildListA = new ArrayList<>();
                buildListS = new ArrayList<>();
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                        if (sliceRegion.select[iSymm][iAtom]) {
                            buildList(iSymm, iAtom, lb, ub);
                        }
                    }
                }
            }
            for (int i = 0; i < buildListA.size(); i++) {
                if (sliceRegion.select[buildListS.get(i)][buildListA.get(i)]) {
                    gridDensity(buildListS.get(i), buildListA.get(i), lb, ub);
                }
            }
        }

        public void buildList(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }
            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = min(aRadGrid, (int) floor(atoms[iAtom].getFormFactorWidth() * fftX / crystal.a) + 1);
            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;
            final int ifrzl = ifrz - frad;
            // Loop over allowed z coordinates for this Loop
            // Check if the current atom is close enough
            // If so, add to list.
            int buff = bufferSize;
            for (int iz = ifrzl - buff; iz <= ifrzu + buff; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                if (giz >= lb && giz <= ub) {
                    buildListA.add(iAtom);
                    buildListS.add(iSymm);
                    break;
                }
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[iAtom].applyLambda()) {
                return;
            }

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
                if (lb > giz || giz > ub) {
                    continue;
                }
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        optLocal[giz]++;
                        actualWeight++;
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = atomff.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class SolventSliceLoop extends SliceLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;
        final int[] optLocal;
        long timer;
        int previousUB, previousLB;
        int actualWeight;

        public SolventSliceLoop(SliceRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
            optLocal = new int[fftZ];
        }

        public double getThreadTime() {
            return timer;
        }

        public int getNumberofSlices() {
            return (previousUB - previousLB + 1);
        }

        public int getThreadWeight() {
            return actualWeight;
        }

        @Override
        public IntegerSchedule schedule() {
            return solventSliceSchedule;
        }

        public void buildList(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }
            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = min(aRadGrid, (int) floor(atoms[iAtom].getFormFactorWidth() * fftX / crystal.a) + 1);
            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;
            final int ifrzl = ifrz - frad;
            // Loop over allowed z coordinates for this Loop
            // Check if the current atom is close enough
            // If so, add to list.
            int buff = bufferSize;
            for (int iz = ifrzl - buff; iz <= ifrzu + buff; iz++) {
                int giz = Crystal.mod(iz, fftZ);
                if (giz >= lb && giz <= ub) {
                    buildListA.add(iAtom);
                    buildListS.add(iSymm);
                    break;
                }
            }
        }

        @Override
        public boolean checkList(int[][] zAtListBuild, int buff) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    if (sliceRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;
                        final int previousZ = zAtListBuild[iSymm][iAtom];
                        if (abs(ifrz - previousZ) >= buff) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        @Override
        public void saveZValues(int[][] zAtListBuild) {
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                    zAtListBuild[iSymm][iAtom] = 0;
                    if (sliceRegion.select[iSymm][iAtom]) {
                        if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                            continue;
                        }
                        xyz[0] = coordinates[iSymm][0][iAtom];
                        xyz[1] = coordinates[iSymm][1][iAtom];
                        xyz[2] = coordinates[iSymm][2][iAtom];
                        crystal.toFractionalCoordinates(xyz, uvw);
                        final double frz = fftZ * uvw[2];
                        final int ifrz = (int) frz;
                        zAtListBuild[iSymm][iAtom] = ifrz;
                    }
                }
            }
        }

        @Override
        public void start() {
            fill(optLocal, 0);
            actualWeight = 0;
            timer = -System.nanoTime();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            boolean boundsChange = false;
            if (previousLB != lb || previousUB != ub) {
                boundsChange = true;
            }
            previousLB = lb;
            previousUB = ub;
            if (rebuildList || boundsChange) {
                buildListA = new ArrayList<>();
                buildListS = new ArrayList<>();
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                        if (sliceRegion.select[iSymm][iAtom]) {
                            buildList(iSymm, iAtom, lb, ub);
                        }
                    }
                }
            }
            for (int i = 0; i < buildListA.size(); i++) {
                if (sliceRegion.select[buildListS.get(i)][buildListA.get(i)]) {
                    gridDensity(buildListS.get(i), buildListA.get(i), lb, ub);
                }
            }
        }

        @Override
        public void finish() {
            timer += System.nanoTime();
            for (int i = 0; i < fftZ; i++) {
                optSliceWeightSolvent.addAndGet(i, optLocal[i]);
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[iAtom].applyLambda()) {
                return;
            }

            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            FormFactor formFactor = solventFormFactors[iSymm][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[iAtom].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case BINARY:
                    frad = min(aRadGrid, (int) floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case GAUSSIAN:
                    frad = min(aRadGrid, (int) floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case POLYNOMIAL:
                    frad = min(aRadGrid, (int) floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
                case NONE:
                default:
                    return;
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
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        optLocal[giz]++;
                        actualWeight++;
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = formFactor.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class BulkSolventSliceLoop extends SliceLoop {

        final double[] xyz = new double[3];
        final double[] uvw = new double[3];
        final double[] xc = new double[3];
        final double[] xf = new double[3];
        final double[] grid;
        final int[] optLocal;
        long timer;
        int[] threadBounds;
        int[] threadWeights;

        public BulkSolventSliceLoop(SliceRegion region) {
            super(region.getNatoms(), region.getNsymm(), region);
            grid = region.getGrid();
            optLocal = new int[fftZ];
        }

        public double getTimePerThread() {
            return timer;
        }

        public int[] getBoundsPerThread() {
            return threadBounds;
        }

        public int[] getWeightPerThread() {
            return threadWeights;
        }

        @Override
        public IntegerSchedule schedule() {
            return bulkSolventSliceSchedule;
        }

        @Override
        public void start() {
            fill(optLocal, 0);
            timer = -System.nanoTime();
        }

        @Override
        public void finish() {
            timer += System.nanoTime();
            for (int i = 0; i < fftZ; i++) {
                optSliceWeightBulkSolvent.addAndGet(i, optLocal[i]);
            }
            threadBounds = bulkSolventSliceSchedule.getLowerBounds().clone();
            threadWeights = bulkSolventSliceSchedule.getThreadWeights().clone();
            for (int i = threadBounds.length - 1; i > 0; i--) {
                threadBounds[i] -= threadBounds[i - 1];
            }
        }

        @Override
        public void gridDensity(int iSymm, int iAtom, int lb, int ub) {
            if (!atoms[iAtom].getUse() && !nativeEnvironmentApproximation) {
                return;
            }

            if (lambdaTerm && atoms[iAtom].applyLambda()) {
                return;
            }

            xyz[0] = coordinates[iSymm][0][iAtom];
            xyz[1] = coordinates[iSymm][1][iAtom];
            xyz[2] = coordinates[iSymm][2][iAtom];
            FormFactor formFactor = solventFormFactors[iSymm][iAtom];
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[iAtom].getVDWType().radius * 0.5;
            int frad = aRadGrid;
            switch (solventModel) {
                case BINARY:
                    frad = min(aRadGrid, (int) floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                    break;
                case GAUSSIAN:
                    frad = min(aRadGrid, (int) floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                    break;
                case POLYNOMIAL:
                    frad = min(aRadGrid, (int) floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                    break;
                case NONE:
                default:
                    return;
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
                xf[2] = iz * ifftZ;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy * ifftY;
                    for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix * ifftX;
                        crystal.toCartesianCoordinates(xf, xc);
                        optLocal[giz]++;
                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = formFactor.rho(grid[ii], 1.0, xc);
                    }
                }
            }
        }
    }

    private class AtomicGradientRegion extends ParallelRegion {

        private final AtomicGradientLoop[] atomicGradientLoop;
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

            final double[] xyz = new double[3];
            final double[] uvw = new double[3];
            final double[] xc = new double[3];
            final double[] xf = new double[3];
            final int[] optLocal = new int[nAtoms];
            long timer;

            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void start() {
                fill(optLocal, 0);
                timer = -System.nanoTime();
            }

            @Override
            public void finish() {
                timer += System.nanoTime();
                for (int i = 0; i < nAtoms; i++) {
                    optAtomicGradientWeight.addAndGet(i, optLocal[i]);
                }
            }

            @Override
            public IntegerSchedule schedule() {
                return atomicGradientSchedule;
            }

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    if (!atoms[n].getUse() && !nativeEnvironmentApproximation) {
                        continue;
                    }
                    if (lambdaTerm && atoms[n].applyLambda()) {
                        continue;
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
                        xf[0] = ix * ifftX;
                        for (int iy = ifry - dfrad; iy <= ifry + dfrad; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy * ifftY;
                            for (int iz = ifrz - dfrad; iz <= ifrz + dfrad; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz * ifftZ;
                                crystal.toCartesianCoordinates(xf, xc);
                                optLocal[n]++;
                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                atomff.rhoGrad(xc, weight * getDensityGrid()[ii], refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    private class SolventGradientRegion extends ParallelRegion {

        private final SolventGradientLoop[] solventGradientLoop;
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

            double[] xyz = new double[3];
            double[] uvw = new double[3];
            double[] xc = new double[3];
            double[] xf = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    if (!atoms[n].getUse() && !nativeEnvironmentApproximation) {
                        continue;
                    }
                    if (lambdaTerm && atoms[n].applyLambda()) {
                        continue;
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
                        case BINARY:
                            dfrad = Math.min(aRadGrid,
                                    (int) Math.floor((vdwr + solventA + 0.2) * fftX / crystal.a) + 1);
                            break;
                        case GAUSSIAN:
                            dfrad = Math.min(aRadGrid,
                                    (int) Math.floor((vdwr * solventB + 2.0) * fftX / crystal.a) + 1);
                            dfcmult = solventA;
                            break;
                        case POLYNOMIAL:
                            dfrad = Math.min(aRadGrid,
                                    (int) Math.floor((vdwr + solventB + 0.2) * fftX / crystal.a) + 1);
                            break;
                        case NONE:
                            return;
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
                        xf[0] = ix * ifftX;
                        for (int iy = ifry - dfrad; iy <= ifryu; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy * ifftY;
                            for (int iz = ifrz - dfrad; iz <= ifrzu; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz * ifftZ;
                                crystal.toCartesianCoordinates(xf, xc);
                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                solventff.rhoGrad(xc,
                                        weight * getDensityGrid()[ii] * dfcmult * getSolventGrid()[ii],
                                        refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    private class BabinetRegion extends ParallelRegion {

        BabinetLoop[] babinetLoops;

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
                                    densityGrid[ii] = 1.0 - getDensityGrid()[ii];
                                }
                            }
                        }
                        break;
                    case GAUSSIAN:
                        for (int k = lb; k <= ub; k++) {
                            for (int j = 0; j < fftY; j++) {
                                for (int i = 0; i < fftX; i++) {
                                    final int ii = iComplex3D(i, j, k, fftX, fftY);
                                    densityGrid[ii] = 1.0 - exp(-solventA * getDensityGrid()[ii]);
                                }
                            }
                        }
                        break;
                    default:
                }
            }
        }
    }

    private class SolventGridRegion extends ParallelRegion {

        SolventGridLoop[] solventGridLoops;

        public SolventGridRegion(int nThreads) {
            solventGridLoops = new SolventGridLoop[nThreads];

        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (solventGridLoops[ti] == null) {
                solventGridLoops[ti] = new SolventGridLoop();
            }

            try {
                execute(0, fftZ - 1, solventGridLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class SolventGridLoop extends IntegerForLoop {

            @Override
            public void run(int lb, int ub) throws Exception {
                // case GAUSSIAN:
                for (int k = lb; k <= ub; k++) {
                    for (int j = 0; j < fftY; j++) {
                        for (int i = 0; i < fftX; i++) {
                            final int ii = iComplex3D(i, j, k, fftX, fftY);
                            solventGrid[ii] = exp(-solventA * getSolventGrid()[ii]);
                        }
                    }
                }
            }
        }
    }

    private class InitRegion extends ParallelRegion {

        InitLoop[] initLoops;
        FormFactorUpdateLoop[] formFactorUpdateLoops;
        int nHKL = reflectionList.hkllist.size();
        double[][] hkldata = null;

        public InitRegion(int nThreads) {
            initLoops = new InitLoop[nThreads];
            formFactorUpdateLoops = new FormFactorUpdateLoop[nThreads];
        }

        public void setHKL(double[][] hkldata) {
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

            private final double[] xyz = new double[3];

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

        AtomicScaleLoop[] atomicScaleLoops;
        int nHKL = reflectionList.hkllist.size();
        double[][] hklData = null;
        double scale;

        public AtomicScaleRegion(int nThreads) {
            atomicScaleLoops = new AtomicScaleLoop[nThreads];
        }

        public void setHKL(double[][] hkldata) {
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
                    double[] fc = hklData[ih.index()];
                    c.re(fc[0]);
                    c.im(fc[1]);
                    // Remove Badd
                    double s = Crystal.invressq(crystal, ih);
                    c.timesIP(scale * exp(0.25 * bAdd * s));
                    c.conjugateIP();
                    fc[0] = c.re();
                    fc[1] = c.im();
                }
            }
        }
    }

    private class SolventScaleRegion extends ParallelRegion {

        SolventScaleLoop[] solventScaleLoops;
        int nHKL = reflectionList.hkllist.size();
        double[][] hkldata = null;
        double scale;

        public SolventScaleRegion(int nThreads) {
            solventScaleLoops = new SolventScaleLoop[nThreads];
        }

        public void setHKL(double[][] hkldata) {
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
                    double[] fc = hkldata[ih.index()];
                    c.re(fc[0]);
                    c.im(fc[1]);
                    c.timesIP(scale);
                    c.conjugateIP();
                    // negative: babinet
                    fc[0] = -c.re();
                    fc[1] = -c.im();
                }
            }
        }
    }

    private class ExtractRegion extends ParallelRegion {

        ExtractLoop[] extractLoops;
        int nHKL = reflectionList.hkllist.size();
        double[][] hkldata = null;

        public ExtractRegion(int nThreads) {
            extractLoops = new ExtractLoop[nThreads];
        }

        public void setHKL(double[][] hkldata) {
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
                    double[] fc = hkldata[ih.index()];
                    // Apply symmetry
                    for (int j = 0; j < nsym; j++) {
                        crystal.applyTransSymRot(ih, ij, symops.get(j));
                        double shift = Crystal.sym_phase_shift(ih, symops.get(j));
                        int h = Crystal.mod(ij.h(), fftX);
                        int k = Crystal.mod(ij.k(), fftY);
                        int l = Crystal.mod(ij.l(), fftZ);
                        if (h < halfFFTX + 1) {
                            final int ii = iComplex3D(h, k, l, fftX, fftY);
                            c.re(getDensityGrid()[ii]);
                            c.im(getDensityGrid()[ii + 1]);
                            c.phaseShiftIP(shift);
                            fc[0] += c.re();
                            fc[1] += c.im();
                        } else {
                            h = (fftX - h) % fftX;
                            k = (fftY - k) % fftY;
                            l = (fftZ - l) % fftZ;
                            final int ii = iComplex3D(h, k, l, fftX, fftY);
                            c.re(getDensityGrid()[ii]);
                            c.im(-getDensityGrid()[ii + 1]);
                            c.phaseShiftIP(shift);
                            fc[0] += c.re();
                            fc[1] += c.im();
                        }
                    }
                }
            }
        }
    }

    private class ArrayCopyRegion extends ParallelRegion {

        ArrayCopyLoop[] arrayCopyLoops;

        public ArrayCopyRegion() {
            arrayCopyLoops = new ArrayCopyLoop[threadCount];
        }

        @Override
        public void run() throws Exception {
            int ti = getThreadIndex();

            if (arrayCopyLoops[ti] == null) {
                arrayCopyLoops[ti] = new ArrayCopyLoop();
            }

            try {
                execute(0, complexFFT3DSpace - 1, arrayCopyLoops[ti]);
            } catch (Exception e) {
                logger.info(e.toString());
            }
        }

        private class ArrayCopyLoop extends IntegerForLoop {

            @Override
            public void run(int lb, int ub) throws Exception {
                arraycopy(getDensityGrid(), lb, getSolventGrid(), lb, ub - lb);
            }
        }
    }
}
