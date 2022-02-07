// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.potential.utils;

import static ffx.crystal.Crystal.mod;
import static ffx.numerics.math.DoubleMath.dist;
import static ffx.potential.parsers.DistanceMatrixFilter.writeDistanceMatrixRow;
import static ffx.potential.utils.Superpose.applyRotation;
import static ffx.potential.utils.Superpose.applyTranslation;
import static ffx.potential.utils.Superpose.calculateRotation;
import static ffx.potential.utils.Superpose.calculateTranslation;
import static ffx.potential.utils.Superpose.rmsd;
import static ffx.potential.utils.Superpose.rotate;
import static ffx.potential.utils.Superpose.translate;
import static ffx.potential.utils.Gyrate.radiusOfGyration;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static java.util.Arrays.sort;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cbrt;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.numerics.math.RunningStatistics;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.DistanceMatrixFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.DoubleIndexPair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

/**
 * Class ProgressiveAlignmentOfCrystals holds the majority of the functionality necessary to quantify
 * crystal similarity following the PAC method.
 *
 * @author Okimasa OKADA, Aaron J. Nessler,  and Michael J. Schnieders
 * @since 1.0
 */
public class ProgressiveAlignmentOfCrystals {

    /**
     * Logger for the ProgressiveAlignmentOfCrystals Class.
     */
    private static final Logger logger = Logger
            .getLogger(ProgressiveAlignmentOfCrystals.class.getName());
    /**
     * SystemFilter containing structures for crystal 1.
     */
    SystemFilter baseFilter;
    /**
     * Number of structures stored in SystemFilter for crystal 1.
     */
    private final int baseSize;
    /**
     * Label for the first crystal.
     */
    private final String baseLabel;
    /**
     * SystemFilter containing structures for crystal 2.
     */
    SystemFilter targetFilter;
    /**
     * Number of structures stored in SystemFilter for crystal 2.
     */
    private final int targetSize;
    /**
     * The amount of work per row for each process.
     */
    private final int numWorkItems;
    /**
     * If this flag is true, then the RMSD matrix is symmetric (e.g., comparing an archive of
     * structures to itself).
     */
    private final boolean isSymmetric;
    /**
     * Label for the second crystal.
     */
    private final String targetLabel;
    /**
     * Label to use for the RSMD logging
     */
    private String rmsdLabel;
    /**
     * Row of RMSD values (length = targetSize).
     */
    public final double[] distRow;
    /**
     * The default restart row is 0. A larger value may be set by the "readMatrix" method if a restart
     * is requested.
     */
    private int restartRow = 0;
    /**
     * The default restart column is 0. Larger values being set by the "readMatrix" method during
     * restarts are not currently supported. A restart must begin from the beginning of a row for
     * simplicity.
     */
    private int restartColumn = 0;
    /**
     * Parallel Java world communicator.
     */
    private final Comm world;
    /**
     * If false, do not use MPI communication.
     */
    private final boolean useMPI;
    /**
     * Number of processes.
     */
    private final int numProc;
    /**
     * Rank of this process.
     */
    private final int rank;
    /**
     * The distances matrix stores a single RSMD value from each process. The array is of size
     * [numProc][1].
     */
    private final double[][] distances;
    /**
     * Each distance is wrapped inside a DoubleBuf for MPI communication.
     */
    private final DoubleBuf[] buffers;
    /**
     * Convenience reference to the RMSD value for this process.
     */
    private final double[] myDistances;
    /**
     * Convenience reference for the DoubleBuf of this process.
     */
    private final DoubleBuf myBuffer;
    /**
     * Mass of N asymmetric units.
     */
    private double[] massN;
    /**
     * Sum of masses for one asymmetric unit.
     */
    private double massSum = 0;
    /**
     * Indices of atoms from first crystal to be used in comparison (Z' = 1).
     */
    private int[] comparisonAtoms;
    /**
     * List of each index's distance to center of expanded crystal for first crystal.
     */
    private DoubleIndexPair[] molDist1;
    /**
     * Working copy of molDist1, updated when treating a new AU as center.
     */
    private DoubleIndexPair[] molDist1_2;
    /**
     * List of each index's distance to center of expanded crystal for second crystal.
     */
    private DoubleIndexPair[] molDist2;
    /**
     * Working copy of molDist2, updated when treating a new AU as center.
     */
    private DoubleIndexPair[] molDist2_2;
    /**
     * List containing coordinates for each unique AU in first crystal.
     */
    private final ArrayList<double[]> tempListXYZ = new ArrayList<>();
    /**
     * Array containing coordinates for each unique AU in first crystal.
     */
    private double[][] baseXYZs;
    /**
     * Array containing coordinates for each unique AU in second crystal.
     */
    private double[][] targetXYZs;
    /**
     * List containing indices for each unique AU in first crystal.
     */
    private final List<Integer> tempListIndices = new ArrayList<>();
    /**
     * Array containing indices for each unique AU in first crystal.
     */
    private Integer[] baseIndices;
    /**
     * Array containing indices for each unique AU in second crystal.
     */
    private Integer[] targetIndices;
    /**
     * Array containing XYZ coordinates for first crystal.
     */
    private double[] baseXYZ;
    /**
     * Array containing XYZ coordinates for second crystal.
     */
    private double[] targetXYZ;
    /**
     * Coordinates of N AUs from first crystal of the closest match (lowest RMSD).
     */
    private double[] bestBaseNAUs;
    /**
     * Coordinates of N AUs from second crystal of the closest match (lowest RMSD).
     */
    private double[] bestTargetNAUs;
    /**
     * Center of masses (or geometric center) for every AU in first replicates crystal
     */
    private double[][] baseCoM;
    /**
     * Coordinates for central AU of first crystal.
     */
    private double[] baseAU;
    /**
     * Coordinates for central 3 AUs of first crystal.
     */
    private double[] base3AUs;
    /**
     * Coordinates for central N AUs of first crystal.
     */
    private double[] baseNAUs;
    /**
     * Center of masses (or geometric center) for every AU in second replicates crystal
     */
    private double[][] targetCoM;
    /**
     * Coordinates for central AU of second crystal.
     */
    private double[] targetAU;
    /**
     * Coordinates for central 3 AUs of second crystal.
     */
    private double[] target3AUs;
    /**
     * Coordinates for central N AUs of second crystal.
     */
    private double[] targetNAUs;
    /**
     * Indices and distances of matching AUs between first and second crystal.
     */
    private DoubleIndexPair[] matchAUs;
    /**
     * String builder to attempt to off load work from logger.
     */
    private final StringBuilder stringBuilder = new StringBuilder();
    /**
     * If molecules between two crystals differ below this tolerance, it is assumed they are equivalent.
     */
    private static final double MATCH_TOLERANCE = 1.0E-12;

    /**
     * Constructor for the ProgressiveAlignmentOfCrystals class.
     *
     * @param baseFilter   SystemFilter containing a set of crystal structures to compare.
     * @param targetFilter SystemFilter containing the other set of crystals to compare.
     */
    public ProgressiveAlignmentOfCrystals(SystemFilter baseFilter,
                                          SystemFilter targetFilter, boolean isSymmetric) {
        this.baseFilter = baseFilter;
        this.targetFilter = targetFilter;
        this.isSymmetric = isSymmetric;

        // Number of models to be evaluated.
        baseSize = baseFilter.countNumModels();
        baseLabel = getName(baseFilter.getFile().getAbsolutePath());
        targetSize = targetFilter.countNumModels();
        targetLabel = getName(targetFilter.getFile().getAbsolutePath());

        assert !isSymmetric || (baseSize == targetSize);

        logger.info(format(" %s conformations: %d", baseLabel, baseSize));
        logger.info(format(" %s conformations: %d", targetLabel, targetSize));

        // Distance matrix to store compared values (dimensions are "human readable" [m x n]).
        distRow = new double[targetSize];

        CompositeConfiguration properties = baseFilter.getActiveMolecularSystem().getProperties();
        useMPI = properties.getBoolean("pj.use.mpi", true);

        if (useMPI) {
            world = Comm.world();
            // Number of processes is equal to world size (often called size).
            numProc = world.size();
            // Each processor gets its own rank (ID of sorts).
            rank = world.rank();
        } else {
            world = null;
            numProc = 1;
            rank = 0;
        }

        // Padding of the target array size (inner loop limit) is for parallelization.
        // Target conformations are parallelized over available nodes.
        // For example, if numProc = 8 and targetSize = 12, then paddedTargetSize = 16.
        int extra = targetSize % numProc;
        int paddedTargetSize = targetSize;
        if (extra != 0) {
            paddedTargetSize = targetSize - extra + numProc;
        }
        numWorkItems = paddedTargetSize / numProc;

        if (numProc > 1) {
            logger.info(format(" Number of MPI Processes:  %d", numProc));
            logger.info(format(" Rank of this MPI Process: %d", rank));
            logger.info(format(" Work per process per row: %d", numWorkItems));
        }

        // Initialize array as -1.0 as -1.0 is not a viable RMSD.
        fill(distRow, -1.0);

        // Each process will complete the following amount of work per row.
        distances = new double[numProc][numWorkItems];

        // Initialize each distance as -1.0.
        for (int i = 0; i < numProc; i++) {
            fill(distances[i], -1.0);
        }

        // DoubleBuf is a wrapper used by MPI Comm methods to transfer data between processors.
        buffers = new DoubleBuf[numProc];
        for (int i = 0; i < numProc; i++) {
            buffers[i] = DoubleBuf.buffer(distances[i]);
        }

        // Convenience reference to the storage for each process.
        myDistances = distances[rank];
        myBuffer = buffers[rank];
    }

    /**
     * Perform default comparison
     *
     * @return RunningStatistics of results.
     */
    public RunningStatistics comparisons() {
        return comparisons(15, 500, 0.1, -1, -1,
                true, false, false, 0, false, 0, false,
                false, false, 0, "default");
    }

    /**
     * Compare the crystals within the SystemFilters that were inputted into the constructor of this
     * class.
     *
     * @param nAU             Number of asymmetric units to compare.
     * @param inflatedAU      Minimum number of asymmetric units in inflated crystal
     * @param matchTol        Tolerance to determine whether two AUs are the same.
     * @param zPrime          Number of asymmetric units in first crystal.
     * @param zPrime2         Number of asymmetric units in second crystal.
     * @param alphaCarbons    Perform comparisons on only alpha carbons.
     * @param noHydrogen      Perform comparisons without hydrogen atoms.
     * @param permute         Compare all unique AUs between crystals.
     * @param save            Save out files of the resulting superposition.
     * @param restart         Try to restart from a previous job.
     * @param write           Save out a PAC RMSD file.
     * @param machineLearning Save out CSV files for machine learning input (saves PDBs as well).
     * @param pacFileName     The filename to use.
     * @return RunningStatistics Statistics for comparisons performed.
     */
    public RunningStatistics comparisons(int nAU, int inflatedAU, double matchTol, int zPrime, int zPrime2,
                                         boolean alphaCarbons, boolean noHydrogen, boolean massWeighted,
                                         int crystalPriority, boolean permute, int save, boolean restart, boolean write,
                                         boolean machineLearning, int linkage, String pacFileName) {
        //TODO: Incorporate graphic user interface (gui: ffx)
        //TODO: Save systems out as original molecule regardless of selection
        //TODO: Handle ring flipping or atoms in equivalent positions (mislabeled atoms)
        //TODO: Handle Z' > 1 for heterogenous/co-crystals.
        RunningStatistics runningStatistics;
        if (restart) {
            runningStatistics = readMatrix(pacFileName, isSymmetric, baseSize, targetSize);
            if (runningStatistics == null) {
                runningStatistics = new RunningStatistics();
            }
        } else {
            runningStatistics = new RunningStatistics();
            File file = new File(pacFileName);
            if (file.exists() && file.delete()) {
                logger.info(format(" PAC RMSD file (%s) was deleted.", pacFileName));
                logger.info(" To restart from a previous run, use the '-r' flag.");
            }
        }

        // restartRow and restartColumn are initialized to zero when this class was constructed.
        // They are updated by the "readMatrix" method if a restart is requested.

        // Read ahead to the base starting conformation.
        for (int row = 0; row < restartRow; row++) {
            baseFilter.readNext(false, false);
        }

        MolecularAssembly baseAssembly = baseFilter.getActiveMolecularSystem();

        // Atom arrays from the 1st assembly.
        Atom[] atoms = baseAssembly.getAtomArray();
        int nAtoms = atoms.length;

        // Collect selected atoms.
        ArrayList<Integer> activeIndices = new ArrayList<>();
        determineActiveAtoms(baseAssembly, activeIndices, alphaCarbons, noHydrogen);

        if (activeIndices.size() < 1) {
            logger.info("\n No atoms were selected for the PAC RMSD in first crystal.");
            return null;
        }

        int[] comparisonAtoms = activeIndices.stream().mapToInt(i -> i).toArray();

        activeIndices.clear();
        // Atom arrays from the 2nd assembly.
        MolecularAssembly targetAssembly = targetFilter.getActiveMolecularSystem();
        Atom[] atoms2 = targetAssembly.getAtomArray();
        int nAtoms2 = atoms2.length;

        // Collect selected atoms.
        determineActiveAtoms(targetAssembly, activeIndices, alphaCarbons, noHydrogen);

        if (activeIndices.size() < 1) {
            logger.info("\n No atoms were selected for the PAC RMSD in second crystal.");
            return null;
        }
        int[] comparisonAtoms2 = activeIndices.stream().mapToInt(i -> i).toArray();

        int compareAtomsSize = comparisonAtoms.length;
        int compareAtomsSize2 = comparisonAtoms2.length;

        //Determine number of species within asymmetric unit.
        //TODO: Handle Z' > 1 for heterogenous/co-crystals.
        int z1 = guessZPrime(zPrime, baseAssembly.getMolecules().size());
        int z2 = guessZPrime(zPrime2, targetAssembly.getMolecules().size());
        // Each ASU contains z * comparisonAtoms species so treat each species individually.
        if (z1 > 1) {
            compareAtomsSize /= z1;
        }
        if (z2 > 1) {
            compareAtomsSize2 /= z2;
        }

        if (compareAtomsSize != compareAtomsSize2) {
            logger.warning(" Selected atom sizes differ between crystals.");
        }

        Crystal baseCrystal = baseAssembly.getCrystal().getUnitCell();
        Crystal targetCrystal = targetFilter.getActiveMolecularSystem().getCrystal().getUnitCell();
        int baseSearchValue = (baseCrystal.spaceGroup.respectsChirality()) ? z1 : 2 * z1;
        int targetSearchValue = (targetCrystal.spaceGroup.respectsChirality()) ? z2 : 2 * z2;

        // Number of used coordinates for atoms in one AU.
        int nCoords = compareAtomsSize * 3;
        // Remove duplicated atoms from Z' > 1.
        if (this.comparisonAtoms == null) {
            this.comparisonAtoms = new int[compareAtomsSize];
            arraycopy(comparisonAtoms, 0, this.comparisonAtoms, 0, compareAtomsSize);
        }
        if (bestBaseNAUs == null) {
            bestBaseNAUs = new double[nAU * nCoords];
        }
        if (bestTargetNAUs == null) {
            bestTargetNAUs = new double[nAU * nCoords];
        }
        if (baseAU == null) {
            baseAU = new double[nCoords];
        }
        if (targetAU == null) {
            targetAU = new double[nCoords];
        }
        if (matchAUs == null || matchAUs.length != Math.max(3, nAU)) {
            matchAUs = new DoubleIndexPair[Math.max(3, nAU)];
        }
        if (base3AUs == null) {
            base3AUs = new double[nCoords * 3];
        }
        if (target3AUs == null) {
            target3AUs = new double[nCoords * 3];
        }
        if (baseNAUs == null) {
            baseNAUs = new double[nAU * nCoords];
        }
        if (targetNAUs == null) {
            targetNAUs = new double[nAU * nCoords];
        }

        int massIndex = 0;
        double[] mass = new double[compareAtomsSize];
        for (Integer value : this.comparisonAtoms) {
            Atom atom = atoms[value];
            double m = atom.getMass();
            mass[massIndex++] = (massWeighted) ? m : 1.0;
        }

        massIndex = 0;
        double[] mass2 = new double[compareAtomsSize2];
        for (Integer value : comparisonAtoms2) {
            if(massIndex == compareAtomsSize2){
                //ComparisonAtoms2 contains all indices for the unit cell.
                break;
            }
            Atom atom = atoms2[value];
            double m = atom.getMass();
            mass2[massIndex++] = (massWeighted) ? m : 1.0;
        }

        if(!Arrays.equals(mass, mass2)){
            logger.warning(" Mass arrays do not match. Check atom ordering.");
            if(logger.isLoggable(Level.FINER)){
                for(int i = 0; i < compareAtomsSize; i++){
                    if(mass[i] != mass2[i]){
                        logger.finer(format(" Index: %d Mass 1: %4.4f Mass 2: %4.4f", i, mass[i], mass2[i]));
                    }
                }
            }
        }

        //Mass of N species
        if (massN == null) {
            massN = (nAU > 3) ? new double[compareAtomsSize * nAU] : new double[compareAtomsSize * 3];
            for (int i = 0; i < nAU; i++) {
                arraycopy(mass, 0, massN, i * compareAtomsSize, compareAtomsSize);
            }
        }
        // Sum of all masses for one species.
        if (massSum == 0) {
            for (int i = 0; i < compareAtomsSize; i++) {
                massSum += massN[i];
            }
        }

        if (machineLearning) {
            save = 1;
        }
        if (save > 2) {
            save = 0;
            logger.info(" Save flag specified incorrectly (1:PDB; 2:XYZ). Not saving files.");
        }

        if (linkage == 0) {
            logger.finer(" Single linkage will be used.");
        } else if (linkage == 2) {
            logger.finer(" Complete linkage will be used.");
        } else if (linkage == 1) {
            logger.finer(" Average linkage will be used.");
        } else {
            logger.warning(
                    "Prioritization method specified incorrectly (--pm {0, 1, 2}). Using default of average linkage.");
            linkage = 1;
        }

        // Number of atoms included in the PAC RMSD.
        logger.info(format("\n %d atoms will be used for the PAC RMSD out of %d in first crystal.", compareAtomsSize * z1, nAtoms));
        logger.info(format(" %d atoms will be used for the PAC RMSD out of %d in second crystal.\n", compareAtomsSize2 * z2, nAtoms2));

        // Label for logging.
        rmsdLabel = format("RMSD_%d", nAU);

        // Minimum amount of time for a single comparison.
        double minTime = Double.MAX_VALUE;
        // Loop over conformations in the base assembly.
        for (int row = restartRow; row < baseSize; row++) {
            // Initialize the distance this rank is responsible for to zero.
            fill(myDistances, -1.0);
            int myIndex = 0;
            // Base unit cell for logging.
            baseAssembly = baseFilter.getActiveMolecularSystem();
            baseCrystal = baseAssembly.getCrystal().getUnitCell();
            double baseDensity = baseCrystal.getDensity(baseAssembly.getMass());
            stringBuilder.setLength(0);
            if (baseCrystal == null || baseCrystal.aperiodic()) {
                stringBuilder.append(" WARNING: Base structure does not have a crystal.\n");
                continue;
            }
            for (int column = restartColumn; column < targetSize; column++) {
                int targetRank = column % numProc;
                if (targetRank == rank) {
                    long time = -System.nanoTime();
                    targetAssembly = targetFilter.getActiveMolecularSystem();
                    targetCrystal = targetAssembly.getCrystal().getUnitCell();
                    if (targetCrystal == null || targetCrystal.aperiodic()) {
                        stringBuilder.append(" WARNING: Target structure does not have a crystal.\n");
                        continue;
                    }
                    double rmsd = -2.0;
                    if (isSymmetric && row == column) {
                        stringBuilder.append(format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s\n",
                                row + 1, baseCrystal.toShortString(), baseLabel,
                                column + 1, targetCrystal.toShortString(), targetLabel));
                        // Fill the diagonal.
                        rmsd = 0.0;
                        // Log the final result.
                        stringBuilder.append(format("\n PAC %s: %12s %7.4f A\n", rmsdLabel, "", rmsd));
                    } else if (isSymmetric && row > column) {
                        // Do not compute lower triangle values.
                        rmsd = -3.0;
                    } else {
                        stringBuilder.append(format("\n Comparing Model %d (%s) of %s\n with      Model %d (%s) of %s\n",
                                row + 1, baseCrystal.toShortString(), baseLabel,
                                column + 1, targetCrystal.toShortString(), targetLabel));
                        double[] gyrations = new double[2];

                        //Setup for comparison with crystal specific information.
                        // Prioritize crystal order based on user specification (High/low density or file order).
                        double targetDensity = targetCrystal.getDensity(targetAssembly.getMass());
                        if (logger.isLoggable(Level.FINER)) {
                            stringBuilder.append(format(" Base Density: %4.4f Target Density: %4.4f\n", baseDensity, targetDensity));
                        }
                        boolean densityCheck = (crystalPriority == 1) ? baseDensity < targetDensity : baseDensity > targetDensity;
                        MolecularAssembly staticAssembly;
                        MolecularAssembly mobileAssembly;
                        // Flip system order if needed.
                        if (densityCheck || crystalPriority == 2) {
                            staticAssembly = baseAssembly;
                            mobileAssembly = targetAssembly;
                        } else {
                            staticAssembly = targetAssembly;
                            mobileAssembly = baseAssembly;
                            int[] tempAtoms = comparisonAtoms.clone();
                            comparisonAtoms = comparisonAtoms2.clone();
                            comparisonAtoms2 = tempAtoms.clone();
                            int temp = z1;
                            z1 = z2;
                            z2 = temp;
                        }

                        //Remove atoms not used in comparisons from the original molecular assembly (crystal 1).
                        staticAssembly.moveAllIntoUnitCell();

                        double[] reducedBaseCoords = reduceSystem(staticAssembly, comparisonAtoms);
                        baseCrystal = staticAssembly.getCrystal();
                        baseXYZ = generateInflatedSphere(baseCrystal, reducedBaseCoords, mass, inflatedAU);

                        //Remove atoms not used in comparisons from the original molecular assembly (crystal 2).
                        mobileAssembly.moveAllIntoUnitCell();
                        double[] reducedTargetCoords = reduceSystem(mobileAssembly, comparisonAtoms2);
                        targetCrystal = mobileAssembly.getCrystal();
                        targetXYZ = generateInflatedSphere(targetCrystal, reducedTargetCoords, mass2,
                                inflatedAU);

                        // Compute the PAC RMSD.
                        rmsd = compare(staticAssembly, mobileAssembly, compareAtomsSize, nAU, baseSearchValue, targetSearchValue,
                                matchTol, row * targetSize + column, gyrations, permute, save,
                                machineLearning, linkage, stringBuilder);
                        time += System.nanoTime();
                        double timeSec = time * 1.0e-9;
                        // Record the fastest comparison.
                        if (timeSec < minTime) {
                            minTime = timeSec;
                        }
                        // Log the final result.
                        double avgGyration = (gyrations[0] + gyrations[1]) / 2;
                        double diff = Math.max(gyrations[0], gyrations[1]) - avgGyration;
                        stringBuilder.append(format("\n PAC %s: %12s %7.4f A (%5.3f sec) G(r) %7.4f A +/- %7.4f\n", rmsdLabel, "", rmsd,
                                timeSec, avgGyration, diff));
                        if(logger.isLoggable(Level.FINER)){
                            stringBuilder.append(format(" Gyration Crystal 1 (%s): %7.4f Crystal 2 (%s): %7.4f.\n", staticAssembly.getName(),
                                    gyrations[0], mobileAssembly.getName(), gyrations[1]));
                        }
                    }
                    myDistances[myIndex] = rmsd;
                    myIndex++;
                }
                targetFilter.readNext(false, false);
            }
            restartColumn = 0;
            targetFilter.readNext(true, false);
            baseFilter.readNext(false, false);

            // Gather RMSDs for this row.
            gatherRMSDs(row, runningStatistics);

            // Write out this row.
            if (rank == 0 && write) {
                int firstColumn = 0;
                if (isSymmetric) {
                    firstColumn = row;
                }
                writeDistanceMatrixRow(pacFileName, distRow, firstColumn);
            }
            logger.info(stringBuilder.toString());
        }

        if (minTime < Double.MAX_VALUE) {
            logger.info(format("\n Minimum PAC time: %7.4f", minTime));
        }

        baseFilter.closeReader();
        targetFilter.closeReader();

        logger.info(format(" RMSD Minimum:  %8.6f", runningStatistics.getMin()));
        logger.info(format(" RMSD Maximum:  %8.6f", runningStatistics.getMax()));
        logger.info(format(" RMSD Mean:     %8.6f", runningStatistics.getMean()));
        double variance = runningStatistics.getVariance();
        if (!Double.isNaN(variance)) {
            logger.info(format(" RMSD Variance: %8.6f", variance));
        }

        // Return distMatrix for validation if this is for the test script
        return runningStatistics;
    }

    /**
     * Perform single comparison between two crystals.
     *
     * @param compareAtomsSize      Number of active atoms in asymmetric unit of first crystal.
     * @param nAU                   Number of asymmetric units to compare.
     * @param matchTol              Tolerance to determine whether two AUs are the same.
     * @param compNum               Comparison number based on all file submitted (logging).
     * @param permute               Compare all unique AUs between crystals.
     * @param save                  Save out files of compared crystals.
     * @param machineLearning       Save out PDBs and CSVs of compared crystals.
     * @param linkage               Criteria to select nearest AUs (0=single, 1=average, 2=complete linkage).
     * @return the computed RMSD.
     */
    private double compare(MolecularAssembly staticAssembly, MolecularAssembly mobileAssembly, int compareAtomsSize, int nAU,
                           int baseSearchValue, int targetSearchValue, double matchTol, int compNum,
                           double[] gyrations, boolean permute, int save, boolean machineLearning,
                           int linkage, StringBuilder stringBuilder) {
        // TODO: Does PAC work for a combination of molecules and polymers?
        // Yes and no. It does not compare them on an individual basis, but can compare AUs as a whole.

        int nCoords = compareAtomsSize * 3;
        //Number of species in expanded crystals.
        int nBaseMols = baseXYZ.length / nCoords;
        int nTargetMols = targetXYZ.length / nCoords;

        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(format(" Number of copies to compare:    %4d\n" +
                    " Number entities in base sphere: %4d\n" +
                    " Number entities in target sphere: %d\n", nAU, nBaseMols, nTargetMols));
        }

        // Translate asymmetric unit of 0th index (closest to all atom center) to the origin.
        translateAUtoOrigin(baseXYZ, baseAU, massN, nCoords, 0);
        molDist1 = new DoubleIndexPair[nBaseMols];
        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(" Prioritize Base System.\n");
        }
        // Center of Masses for crystal 1.
        baseCoM = new double[nBaseMols][3];
        //Update coordinates for each new comparison.
        centerOfMass(baseCoM, baseXYZ, massN, massSum, compareAtomsSize);
        prioritizeReplicates(baseXYZ, massN, massSum, baseCoM, compareAtomsSize, molDist1, 0, linkage);

        //Used for debugging. can be removed.
        if (logger.isLoggable(Level.FINEST)) {
            int printSize = 20;
            stringBuilder.append(" System 1 distances to center of sphere:\n");
            for (int i = 0; i < printSize; i++) {
                stringBuilder.append(format(" %d\t%16.8f\n", molDist1[i].getIndex(), molDist1[i].getDoubleValue()));
            }
        }

        // Translate system to the origin.
        translateAUtoOrigin(targetXYZ, targetAU, massN, nCoords, 0);
        molDist2 = new DoubleIndexPair[nTargetMols];

        // Reorder molDist2 as we shift a different molecule (m) to the center each loop.
        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(" Prioritize target system.\n");
        }
        targetCoM = new double[nTargetMols][3];
        centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
        prioritizeReplicates(targetXYZ, massN, massSum, targetCoM, compareAtomsSize, molDist2, 0, linkage);

        if (logger.isLoggable(Level.FINEST)) {
            int printSize = 20;
            stringBuilder.append(" System 2 distances to center of sphere:\n");
            for (int i = 0; i < printSize; i++) {
                stringBuilder.append(format(" %d\t%16.8f\n", molDist2[i].getIndex(), molDist2[i].getDoubleValue()));
            }
        }

        //Determine if AUs in first system are same hand as center most in first (stereoisomer handling).
        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(" Search Conformations of Base Crystal:\n");
        }
        tempListXYZ.clear();
        tempListIndices.clear();
        numberUniqueAUs(baseXYZ, molDist1, baseAU, nCoords, baseSearchValue, permute,
                nBaseMols, massN, matchTol);
        if(baseXYZs == null || baseXYZs.length != tempListXYZ.size()) {
            baseXYZs = new double[tempListXYZ.size()][nCoords];
        }
        tempListXYZ.toArray(baseXYZs);
        if(baseIndices == null || baseIndices.length != tempListIndices.size()) {
            baseIndices = new Integer[tempListIndices.size()];
        }
        tempListIndices.toArray(baseIndices);
        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(format(" %d conformations detected out of %d in base crystal.\n" +
                    " Search Conformations of Target Crystal:\n", baseIndices.length, baseSearchValue));
        }
        tempListXYZ.clear();
        tempListIndices.clear();
        numberUniqueAUs(targetXYZ, molDist2, targetAU, nCoords, targetSearchValue, permute,
                nTargetMols, massN, matchTol);
        if(targetXYZs == null || targetXYZs.length != tempListXYZ.size()) {
            targetXYZs = new double[tempListXYZ.size()][nCoords];
        }
        tempListXYZ.toArray(targetXYZs);
        int targetLength = tempListIndices.size();
        if(targetIndices == null || targetIndices.length != targetLength) {
            targetIndices = new Integer[tempListIndices.size()];
        }
        tempListIndices.toArray(targetIndices);
        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(format(" %d conformations detected out of %d in target crystal.\n", targetIndices.length, targetSearchValue));
        }

        // Determine which unique AUs are most similar between crystals
        //Minimum difference between each unique target AU and closest matching base AU.
        double[] targetBaseDiff = new double[targetLength];
        // Index of the closest matching base AU to each target AU.
        int[] targetBaseIndices = new int[targetLength];
        int baseLength = baseXYZs.length;
        targetLength = targetXYZs.length;
        for (int i = 0; i < targetLength; i++) {
            int minIndex = -1;
            double minDiff = Double.MAX_VALUE;
            for (int j = 0; j < baseLength; j++) {
                double value = RMSD_1(targetXYZs[i], baseXYZs[j], massN);
                if (value < minDiff) {
                    minDiff = value;
                    minIndex = j;
                }
            }
            targetBaseDiff[i] = minDiff;
            targetBaseIndices[i] = baseIndices[minIndex];
        }
        if (logger.isLoggable(Level.FINER)) {
            stringBuilder.append(" Minimum RMSD_1 Between Unique Target and Base AUs:\n i tInd RMSD_1    bInd\n");
            int tbIndicesLength = targetBaseIndices.length;
            for (int i = 0; i < tbIndicesLength; i++) {
                stringBuilder.append(format(" %d %d %4.4f %d", i, targetIndices[i], targetBaseDiff[i], targetBaseIndices[i]));
            }
        }

        // Coordinate arrays to save out structures at the end.
        double bestRMSD = Double.MAX_VALUE;
        if (logger.isLoggable(Level.FINE)) {
            stringBuilder.append(format("\n  Trial     RMSD_1 (%7s)  RMSD_3 (%7s)  %7s  G(r1)   G(r2)\n",
                    rmsdLabel, rmsdLabel, rmsdLabel));
        }
        molDist1_2 = new DoubleIndexPair[nBaseMols];
        molDist2_2 = new DoubleIndexPair[nTargetMols];
        // Begin comparison
        // Integer used only for user display logging.
        int currentComparison = 1;
        for (Integer index : baseIndices) {
            int center = molDist1[index].getIndex();

            //Re-prioritize based on center-most molecule if different from first linkage.
            if (center != 0) {
                // Reprioritize replicates crystal based on new center.
                prioritizeReplicates(baseXYZ, massN, massSum, baseCoM, compareAtomsSize, molDist1_2, center, linkage);
            } else {
                // Reuse original reprioritization.
                arraycopy(molDist1, 0, molDist1_2, 0, nBaseMols);
            }

            //Translate base system based on center-most molecule
            int baseAUIndex = molDist1_2[0].getIndex() * nCoords;
            arraycopy(baseXYZ, baseAUIndex, baseAU, 0, nCoords);
            double[] translation = calculateTranslation(baseAU, massN);
            applyTranslation(baseAU, translation);
            applyTranslation(baseXYZ, translation);

            //Update CoMs with translation
            centerOfMass(baseCoM, baseXYZ, massN, massSum, compareAtomsSize);

            // Acquire coordinates based on center 3 molecules
            if (logger.isLoggable(Level.FINER)) {
                stringBuilder.append(" Base 3 Conformations\n");
            }
            for (int i = 0; i < 3; i++) {
                int baseIndex = molDist1_2[i].getIndex() * nCoords;
                arraycopy(baseXYZ, baseIndex, base3AUs, i * nCoords, nCoords);
            }

            // Acquire coordinates for final comparison
            for (int i = 0; i < nAU; i++) {
                int molIndex = molDist1_2[i].getIndex() * nCoords;
                arraycopy(baseXYZ, molIndex, baseNAUs, i * nCoords, nCoords);
            }

            int targetConformations = targetXYZs.length;
            for (int m = 0; m < targetConformations; m++) {
                if (permute || Objects.equals(targetBaseIndices[m], index)) {
                    // Switch m center most molecules (looking for stereoisomers)
                    center = molDist2[targetIndices[m]].getIndex();
                    //Re-prioritize based on central AU if different from first prioritization.
                    if (logger.isLoggable(Level.FINER)) {
                        stringBuilder.append(" Re-prioritize target system.\n");
                    }
                    if (center != 0) {
                        // Reprioritize replicates crystal based on new center.
                        prioritizeReplicates(targetXYZ, massN, massSum, targetCoM, compareAtomsSize, molDist2_2, center, linkage);
                    } else {
                        // Reuse original reprioritization.
                        arraycopy(molDist2, 0, molDist2_2, 0, nTargetMols);
                    }

                    if (logger.isLoggable(Level.FINER)) {
                        stringBuilder.append(" Rotation 1:\n");
                    }

                    firstRotation(nCoords, molDist2_2[0].getIndex());

                    // At this point both systems have completed first rotation/translation
                    //  Therefore both center-most molecules should be overlapped.
                    if (logger.isLoggable(Level.FINER)) {
                        stringBuilder.append(" Match molecules between systems.\n");
                    }

                    //Update center of masses with the first trans/rot
                    centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
                    matchMolecules();

                    double checkRMSD1 = -1.0;
                    double n1RMSD = -1.0;
                    if (logger.isLoggable(Level.FINE)) {
                        int targetCenterMol = molDist2_2[0].getIndex() * nCoords;
                        arraycopy(targetXYZ, targetCenterMol, targetAU, 0, nCoords);
                        checkRMSD1 = rmsd(baseAU, targetAU, massN);
                        if (logger.isLoggable(Level.FINER)) {
                            stringBuilder.append(format(" Center Molecule RMSD after rot 1: %16.8f\n", checkRMSD1));
                        }
                        for (int i = 0; i < nAU; i++) {
                            int offset = i * nCoords;
                            int molIndex = matchAUs[i].getIndex() * nCoords;
                            arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
                        }
                        n1RMSD = rmsd(baseNAUs, targetNAUs, massN);
                        if (logger.isLoggable(Level.FINEST) && save > 0) {
                            saveAssembly(staticAssembly, baseNAUs, comparisonAtoms, "_c1", compNum, save);
                            saveAssembly(mobileAssembly, targetNAUs, comparisonAtoms, "_c2", compNum, save);
                        }
                    }

                    if (logger.isLoggable(Level.FINER)) {
                        stringBuilder.append(" Rotation 2:\n");
                    }

                    secondRotation(nCoords);

                    double checkRMSD2 = -1.0;
                    double n3RMSD = -1.0;
                    if (logger.isLoggable(Level.FINE)) {
                        checkRMSD2 = rmsd(base3AUs, target3AUs, massN);
                        for (int i = 0; i < nAU; i++) {
                            int offset = i * nCoords;
                            int molIndex = matchAUs[i].getIndex() * nCoords;
                            arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
                        }
                        n3RMSD = rmsd(baseNAUs, targetNAUs, massN);
                        if (logger.isLoggable(Level.FINEST) && save > 0) {
                            saveAssembly(staticAssembly, baseNAUs, comparisonAtoms, "_c1", compNum, save);
                            saveAssembly(mobileAssembly, targetNAUs, comparisonAtoms, "_c2", compNum, save);
                        }
                    }

                    //Update center of masses with the second rot (only one crystal moves).
                    centerOfMass(targetCoM, targetXYZ, massN, massSum, compareAtomsSize);
                    // Rotations 1 and 2 have been completed and both systems should be overlapped
                    //  Isolate center most nAU from System 1 and matching molecules from System 2
                    if (logger.isLoggable(Level.FINER)) {
                        stringBuilder.append(" Match Molecules:\n");
                    }
                    if (matchAUs.length != nAU) {
                        matchAUs = new DoubleIndexPair[nAU];
                    }
                    matchMolecules();

                    for (int i = 0; i < nAU; i++) {
                        int offset = i * nCoords;
                        int molIndex = matchAUs[i].getIndex() * nCoords;
                        arraycopy(targetXYZ, molIndex, targetNAUs, offset, nCoords);
                    }

                    translate(baseNAUs, massN, targetNAUs, massN);
                    rotate(baseNAUs, targetNAUs, massN);
                    double rmsdSymOp = rmsd(baseNAUs, targetNAUs, massN);

                    double baseGyration = radiusOfGyration(baseNAUs);
                    double targetGyration = radiusOfGyration(targetNAUs);

                    if (logger.isLoggable(Level.FINE)) {
                        int totalComparisons = (permute) ? baseXYZs.length * targetConformations : Math.min(baseXYZs.length, targetConformations);
                        String output = format(" %2d of %2d: %7.4f (%7.4f) %7.4f (%7.4f) %7.4f A  %7.4f A %7.4f A",
                                currentComparison, totalComparisons, checkRMSD1, n1RMSD, checkRMSD2, n3RMSD, rmsdSymOp,
                                baseGyration, targetGyration);

                        if (logger.isLoggable(Level.FINER)) {
                            output += format(" %d=%d", index, targetBaseIndices[m]);
                        }
                        stringBuilder.append(output).append("\n");
                        if (logger.isLoggable(Level.FINEST) && save > 0) {
                            saveAssembly(staticAssembly, baseNAUs, comparisonAtoms, "_c1", compNum, save);
                            saveAssembly(mobileAssembly, targetNAUs, comparisonAtoms, "_c2", compNum, save);
                        }
                    }

                    if (rmsdSymOp < bestRMSD) {
                        gyrations[0] = baseGyration;
                        gyrations[1] = targetGyration;
                        bestRMSD = rmsdSymOp;
                        bestBaseNAUs = baseNAUs;
                        bestTargetNAUs = targetNAUs;
                    }
                    currentComparison++;
                }
            }
        }

        double finalRMSD;
        if (bestRMSD < Double.MAX_VALUE) {
            finalRMSD = bestRMSD;
        } else {
            stringBuilder.append(" This RMSD was filtered out! Try the --ex flag." +
                    "\nAlternatively increase --ns and/or --ns2.\n");
            // TODO: Double.NaN causes an error in RunningStatistics... Set to -4.0 for now...
            finalRMSD = -4.0;
        }
        if (save > 0) {
            if (machineLearning) {
                saveAssembly(staticAssembly, bestBaseNAUs, comparisonAtoms, "_c1", 0.000, compNum, save);
                saveAssembly(mobileAssembly, bestTargetNAUs, comparisonAtoms, "_c2", finalRMSD, compNum, save);
            }else {
                saveAssembly(staticAssembly, bestBaseNAUs, comparisonAtoms, "_c1", compNum, save);
                saveAssembly(mobileAssembly, bestTargetNAUs, comparisonAtoms, "_c2", compNum, save);
            }
        }

        return finalRMSD;
    }

    /**
     * Read in the distance matrix.
     *
     * @param filename        The PAC RMSD matrix file to read from.
     * @param isSymmetric     Is the distance matrix symmetric.
     * @param expectedRows    The expected number of rows.
     * @param expectedColumns The expected number of columns.
     * @return Stats for all read in distance matrix values.
     */
    private RunningStatistics readMatrix(String filename, boolean isSymmetric, int expectedRows,
                                         int expectedColumns) {
        restartRow = 0;
        restartColumn = 0;

        DistanceMatrixFilter distanceMatrixFilter = new DistanceMatrixFilter();
        RunningStatistics runningStatistics = distanceMatrixFilter.readDistanceMatrix(
                filename, expectedRows, expectedColumns);

        if (runningStatistics != null && runningStatistics.getCount() > 0) {
            restartRow = distanceMatrixFilter.getRestartRow();
            restartColumn = distanceMatrixFilter.getRestartColumn();

            if (isSymmetric) {
                // Only the diagonal entry (0.0) is on the last row for a symmetric matrix.
                if (restartRow == expectedRows && restartColumn == 1) {
                    logger.info(format(" Complete symmetric distance matrix found (%d x %d).", restartRow,
                            restartRow));
                } else {
                    restartColumn = 0;
                    logger.info(format(
                            " Incomplete symmetric distance matrix found.\n Restarting at row %d, column %d.",
                            restartRow + 1, restartColumn + 1));
                }
            } else if (restartRow == expectedRows && restartColumn == expectedColumns) {
                logger.info(format(" Complete distance matrix found (%d x %d).", restartRow, restartColumn));
            } else {
                restartColumn = 0;
                logger.info(format(" Incomplete distance matrix found.\n Restarting at row %d, column %d.",
                        restartRow + 1, restartColumn + 1));
            }
        }

        return runningStatistics;
    }

    /**
     * This method calls <code>world.gather</code> to collect numProc PAC RMSD values.
     *
     * @param row               Current row of the PAC RMSD matrix.
     * @param runningStatistics Stats for the RMSDs.
     */
    private void gatherRMSDs(int row, RunningStatistics runningStatistics) {
        if (useMPI) {
            try {
                if (logger.isLoggable(Level.FINER)) {
                    logger.finer(" Receiving results.");
                }
                // TODO: Node 0 is the only process that updates. Remove from other nodes.
                world.gather(0, myBuffer, buffers);
                for (int workItem = 0; workItem < numWorkItems; workItem++) {
                    for (int proc = 0; proc < numProc; proc++) {
                        int column = numProc * workItem + proc;
                        // Do not include padded results.
                        if (column < targetSize) {
                            distRow[column] = distances[proc][workItem];
                            if (!isSymmetric) {
                                runningStatistics.addValue(distRow[column]);
                            } else if (column >= row) {
                                // Only collect stats for the upper triangle.
                                runningStatistics.addValue(distRow[column]);
                            }
                            if (logger.isLoggable(Level.FINER)) {
                                logger.finer(format(" %d %d %16.8f", row, column, distances[proc][workItem]));
                            }
                        }
                    }
                }
            } catch (Exception e) {
                logger.severe(" Exception collecting distance values." + e);
            }
        } else {
            for (int i = 0; i < targetSize; i++) {
                distRow[i] = myDistances[i];
                if (!isSymmetric) {
                    runningStatistics.addValue(distRow[i]);
                } else if (i >= row) {
                    // Only collect stats for the upper triangle.
                    runningStatistics.addValue(distRow[i]);
                }
            }
        }
    }

    /**
     * Try to automatically determine number of species in asymmetric unit (only works for molecules).
     *
     * @param zPrime     User input overrides detection method.
     * @param numEntities Number of species detected.
     * @return Number of expected species in asymmetric unit.
     */
    private static int guessZPrime(int zPrime, int numEntities) {
        int z = (zPrime > 0) ? zPrime : Math.max(numEntities, 0);
        if (z < 1) {
            z = 1;
        }
        if(z!=zPrime){
            logger.info(" Z' automatically generated. Does not work for heterogenous/co-crystals (Z' can be set manually with the --zp and --zp2 flags).");
        }
        if (logger.isLoggable(Level.FINER)) {
            logger.finer(format(" Number of species in asymmetric unit (Z'): %d", z));
        }
        return z;
    }

    /**
     * Determine the number of unique AUs within the replicates crystal to a tolerance.
     *
     * @param allCoords       Coordinates for every atom in replicates crystal ([x1, y1, z1, x2, y2, z2...].
     * @param molIndices      Prioritization of molecules.
     * @param auCoords        Coordinates for single AU.
     * @param nCoords         Number of coordinates in an AU (number of atoms * 3).
     * @param upperLimit      The largest number of unique AUs (0 for no upper limit).
     * @param permute         Search entire replicates crystal if true, otherwise only the expected.
     * @param nAUinReplicates Number of AUs in replicates crystal.
     * @param massN            Array containing masses for each atom in AU.
     * @param matchTol        Tolerance to determine whether two AUs are the same.
     */
    private void numberUniqueAUs(double[] allCoords, DoubleIndexPair[] molIndices, double[] auCoords, int nCoords, int upperLimit, boolean permute,
                                        int nAUinReplicates, double[] massN, double matchTol) {
        // uniqueDiffs is only recorded for logging... could remove later.
        // List of differences (RMSD_1) for AUs in replicates crystal.
        List<Double> uniqueDiffs = new ArrayList<>();
        tempListIndices.add(0);
        uniqueDiffs.add(rmsd(auCoords, auCoords, massN));
        tempListXYZ.add(auCoords);
        // Start from 1 as zero is automatically added.
        int index = 1;
        //Determine number of conformations in first crystal
        int numConfCheck = (permute || upperLimit <= 0) ? nAUinReplicates : upperLimit;
        if (logger.isLoggable(Level.FINEST)) {
            logger.finest("RMSD Differences in Replicates Crystal:");
        }
        while (uniqueDiffs.size() < numConfCheck) {
            double[] baseCheckMol = new double[nCoords];
            arraycopy(allCoords, molIndices[index].getIndex() * nCoords, baseCheckMol, 0,
                    nCoords);
            double value = RMSD_1(tempListXYZ.get(0), baseCheckMol, massN);
            if (logger.isLoggable(Level.FINEST)) {
                logger.finest(format("%d %4.4f", molIndices[index].getIndex(), value));
            }
            if (addLooseUnequal(uniqueDiffs, value, matchTol)) {
                tempListIndices.add(index);
                tempListXYZ.add(baseCheckMol);
            }
            index++;
            if (index >= molIndices.length) {
                break;
            }
        }
        if (logger.isLoggable(Level.FINER)) {
            logger.finer(" RMSD_1 from 1st AU:\n i RMSD     AUIndex");
            int numUniqueIndices = tempListIndices.size();
            for (int i = 0; i < numUniqueIndices; i++) {
                logger.finer(format(" %d %4.4f %d", i, uniqueDiffs.get(i), tempListIndices.get(i)));
            }
        }
    }

    /**
     * Perform first rotation to match center most AUs.
     *
     * @param nCoords   Number of coordinates in an entity.
     * @param index     Index of target entity to move to base entity.
     */
    private void firstRotation(int nCoords, int index) {
        translateAUtoOrigin(targetXYZ, targetAU, massN, nCoords, index);

        // Rotate the target molecule onto the base molecule.
        double[][] rotation = calculateRotation(baseAU, targetAU, massN);
        applyRotation(targetAU, rotation);
        applyRotation(targetXYZ, rotation);
    }

    /**
     * Translate an asymmetric unit to the origin.
     *
     * @param systemXYZ Coordinates of the system.
     * @param auCoords  Coordinates for an AU.
     * @param mass      Masses of atoms in one asymmetric unit.
     * @param nCoords   Number of coordinates in an entity.
     * @param index     Index of the asymmetric unit to move.
     */
    private static void translateAUtoOrigin(double[] systemXYZ, double[] auCoords, double[] mass, int nCoords, int index) {

        // Load the coordinates.
        int auIndex = index * nCoords;
        arraycopy(systemXYZ, auIndex, auCoords, 0, nCoords);

        double[] translation = calculateTranslation(auCoords, mass);
        applyTranslation(auCoords, translation);
        applyTranslation(systemXYZ, translation);
    }

    /**
     * Perform second rotation to better match crystal systems.
     *
     * @param nCoords   Number of coordinates in an entity.
     */
    private void secondRotation(int nCoords) {

        // Load coordinates for 3 molecules for the base and target systems
        for (int i = 0; i < 3; i++) {
            int targetIndex = matchAUs[i].getIndex() * nCoords;
            int index = i * nCoords;
            arraycopy(targetXYZ, targetIndex, target3AUs, index, nCoords);
        }

        // Calculate the rotation matrix and apply it to the target system.
        double[][] rotation = calculateRotation(base3AUs, target3AUs, massN);
        applyRotation(target3AUs, rotation);
        applyRotation(targetXYZ, rotation);
    }

    /**
     * Pair species between two crystals based on center of mass distances.
     *
     */
    private void matchMolecules() {
        int desiredMols = matchAUs.length;
        tempListIndices.clear();
        // List of indexes for second system.
        for (DoubleIndexPair doubleIndexPair : molDist2_2) {
            // Only search molecules within range of the desired number of molecules.
            // Must have enough molecules for matching (using exhaustive till better heuristic is determined)
            tempListIndices.add(doubleIndexPair.getIndex());
        }

        // Compare distances between center of masses from system 1 and 2.
        for (int i = 0; i < desiredMols; i++) {
            double minDist = Double.MAX_VALUE;
            Integer minIndex = -1;
            double[] baseCoMCurr = baseCoM[molDist1_2[i].getIndex()];
            for (Integer target : tempListIndices) {
                double dist = dist(baseCoMCurr, targetCoM[target]);
                if (dist < minDist) {
                    minIndex = target;
                    minDist = dist;
                    if(logger.isLoggable(Level.FINEST)) {
                        stringBuilder.append(format(" \t Index: %d Dist: %7.4f\n", minIndex, minDist));
                    }
                }
                if (abs(minDist) < MATCH_TOLERANCE) {
                    // Distance between center of masses is ~0 is the best scenario assuming no coordinate overlaps.
                    if(logger.isLoggable(Level.FINEST)){
                        stringBuilder.append(" \tExit out of match loop.\n");
                    }
                    break;
                }
            }
            matchAUs[i] = new DoubleIndexPair(minIndex, minDist);
            if (!tempListIndices.remove(minIndex)) {
                logger.warning(format(" Index value of %d was not found (%4.4f).", minIndex, minDist));
            }
            if (logger.isLoggable(Level.FINER)) {
                stringBuilder.append(
                        format("\n Base position:   %d: %8.4f %8.4f %8.4f\n", i, baseCoM[molDist1_2[i].getIndex()][0],
                                baseCoM[molDist1_2[i].getIndex()][1],
                                baseCoM[molDist1_2[i].getIndex()][2]));
                stringBuilder.append(format(" Match Distance:  %d: %8.4f\n", i, matchAUs[i].getDoubleValue()));
                stringBuilder.append(format(" Target position: %d: %8.4f %8.4f %8.4f\n", i,
                        targetCoM[matchAUs[i].getIndex()][0],
                        targetCoM[matchAUs[i].getIndex()][1], targetCoM[matchAUs[i].getIndex()][2]));
            }
        }
        if (logger.isLoggable(Level.FINEST)) {
            stringBuilder.append("  Distance between pairs after rot 2:\n" +
                    " Index  Base Index  Target Index    Match Index    Distance\n");
            for (int i = 0; i < matchAUs.length; i++) {
                stringBuilder.append(format(" %5d %10d %14d %14d %10.4f\n", i,
                        molDist1_2[i].getIndex(), molDist2_2[i].getIndex(), matchAUs[i].getIndex(), matchAUs[i].getDoubleValue()));
            }
        }
    }

    /**
     * Determine the RMSD between two AUs (should be same AUs).
     *
     * @param baseAU   Coordinates for first AU.
     * @param targetAU Coordinates for second AU.
     * @param mass      Mass of atoms within AU.
     * @return RMSD between AUs.
     */
    private static double RMSD_1(double[] baseAU, double[] targetAU, double[] mass) {
        translate(baseAU, mass, targetAU, mass);
        rotate(baseAU, targetAU, mass);
        return rmsd(baseAU, targetAU, mass);
    }

    /**
     * Calculate the center of mass for a given set of masses for the asymmetric unit and coordinates
     * (xyz)
     *
     * @param centersOfMass Returned center of mass for each asymmetric unit
     * @param coords        Coordinates of every atom in system.
     * @param mass          Masses of each atom in asymmetric unit.
     * @param massSum       Sum of masses within asymmetric unit.
     * @param nAtoms        Number of coordinates in an entity.
     */
    private static void centerOfMass(double[][] centersOfMass, double[] coords, double[] mass,
                                     double massSum, int nAtoms) {
        int size = centersOfMass.length;
        for (int i = 0; i < size; i++) {
            int molIndex = i * nAtoms * 3;
            double[] comI = new double[3];
            for (int j = 0; j < nAtoms; j++) {
                int atomIndex = j * 3;
                double m = mass[j];
                for (int k = 0; k < 3; k++) {
                    comI[k] += coords[molIndex + atomIndex + k] * m;
                }
            }
            for (int j = 0; j < 3; j++) {
                comI[j] /= massSum;
            }
            centersOfMass[i] = comI;
        }
    }

    /**
     * Reduce asymmetric unit to atoms that are going to be used in final RMSD.
     *
     * @param assembly        Asymmetric unit we wish to reduce.
     * @param comparisonAtoms Atoms of interest within asymmetric unit.
     * @return Linear coordinates for only atoms of interest.
     */
    private static double[] reduceSystem(MolecularAssembly assembly, int[] comparisonAtoms) {
        Atom[] atoms = assembly.getAtomArray();
        // Collect asymmetric unit atomic coordinates.
        double[] reducedCoords = new double[comparisonAtoms.length * 3];
        int coordIndex = 0;
        for (Integer value : comparisonAtoms) {
            Atom atom = atoms[value];
            reducedCoords[coordIndex++] = atom.getX();
            reducedCoords[coordIndex++] = atom.getY();
            reducedCoords[coordIndex++] = atom.getZ();
        }
        return reducedCoords;
    }

    /**
     * Determine the indices of the atoms from the assembly that are active for this comparison.
     *
     * @param assembly Assembly of interest to compare.
     * @param indices  Array list containing atom indices that will be used for this comparison.
     * @param alphaCarbons Boolean whether to include only alpha carbons/nitrogens.
     * @param noHydrogen Boolean whether to include hydrogens.
     */
    private static void determineActiveAtoms(MolecularAssembly assembly, ArrayList<Integer> indices, boolean alphaCarbons,
                                             boolean noHydrogen) {
        Atom[] atoms = assembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            if (atom.isActive()) {
                if (alphaCarbons) {
                    String atomName = atom.getName();
                    int atomAtNum = atom.getAtomicNumber();
                    boolean proteinCheck = atomName.equals("CA") && atomAtNum == 6;
                    boolean aminoCheck = (atomName.equals("N1") || atomName.equals("N9")) && atomAtNum == 7;
                    if (proteinCheck || aminoCheck) {
                        indices.add(i);
                    }
                } else if (!noHydrogen || !atom.isHydrogen()) {
                    indices.add(i);
                }
            }
            // Reset all atoms to active once the selection is recorded.
            atom.setActive(true);
        }
    }

    /**
     * Generate and expanded sphere of asymmetric unit with the intention of observing a crystals
     * distribution of replicates rather to facilitate comparisons that go beyond lattice parameters.
     *
     * @param crystal       Crystal to define expansion.
     * @param reducedCoords Coordinates of asymmetric unit we wish to expand.
     * @param mass          Masses for atoms within reduced asymmetric unit.
     * @param inflatedAU    Number of asymmetric units after inflation.
     * @return double[] containing the coordinates for the expanded crystal.
     */
    private static double[] generateInflatedSphere(Crystal crystal, double[] reducedCoords,
                                                   double[] mass, int inflatedAU) {
        int nAtoms = reducedCoords.length / 3;
        int zAtoms = mass.length;
        int zPrime = nAtoms/mass.length;
        // Collect asymmetric unit atomic coordinates.
        double[] x = new double[nAtoms];
        double[] y = new double[nAtoms];
        double[] z = new double[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            int atomIndex = i * 3;
            x[i] = reducedCoords[atomIndex];
            y[i] = reducedCoords[atomIndex + 1];
            z[i] = reducedCoords[atomIndex + 2];
        }

        // When the system was read in, a replicates crystal may have been created to satisfy the cutoff.
        // Retrieve a reference to the unit cell (not the replicates crystal).
        // Here we will use the unit cell, to create a new replicates crystal that may be
        // a different size (i.e. larger).
        Crystal unitCell = crystal.getUnitCell();

        double asymmetricUnitVolume = unitCell.volume / unitCell.getNumSymOps() / zPrime;

        // Estimate a radius that will include desired number of asymmetric units (inflatedAU).
        double radius = cbrt(inflatedAU * asymmetricUnitVolume) / 2;
        Crystal replicatesCrystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, radius * 2.0);

        // Symmetry coordinates for each molecule in replicates crystal
        int nSymm = replicatesCrystal.getNumSymOps();

        if (logger.isLoggable(Level.FINER)) {
            logger.finer(format(" Desired copies in target sphere:     %3d", inflatedAU));
            logger.finer(format(" Asymmetric Unit Volume:  %4.2f", asymmetricUnitVolume));
            logger.finer(format(" Estimated spherical radius:  %4.2f", radius));
            logger.finer(" Replicates crystal " + replicatesCrystal);
            logger.finer(format(" Number of entities in replicates: %3d\n\n", nSymm * zPrime));
        }

        double[][] xS = new double[nSymm][nAtoms];
        double[][] yS = new double[nSymm][nAtoms];
        double[][] zS = new double[nSymm][nAtoms];
        // Cartesian center of each molecule
        int numEntities = nSymm * zPrime;
        double[][] centerMolsCart = new double[numEntities][3];

        // Loop over replicate crystal SymOps
        List<SymOp> inflatedSymOps = replicatesCrystal.spaceGroup.symOps;
        for (int iSym = 0; iSym < nSymm; iSym++) {
            SymOp symOp = inflatedSymOps.get(iSym);
            // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
            replicatesCrystal.applySymOp(nAtoms, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
            for(int zp = 0; zp < zPrime; zp++) {
                // Compute center-of-mass (CoM) for Cartesian coordinates
                double[] centerOfMass = new double[3];
                double totalMass = 0.0;
                for (int i = 0; i < zAtoms; i++) {
                    double m = mass[i];
                    int symIndex = zp * zAtoms;
                    centerOfMass[0] += xS[iSym][i + symIndex] * m;
                    centerOfMass[1] += yS[iSym][i + symIndex] * m;
                    centerOfMass[2] += zS[iSym][i + symIndex] * m;
                    totalMass += m;
                }
                centerOfMass[0] /= totalMass;
                centerOfMass[1] /= totalMass;
                centerOfMass[2] /= totalMass;

                double[] translate = moveIntoCrystal(replicatesCrystal, centerOfMass);
                for (int i = 0; i < zAtoms; i++) {
                    int symIndex = zp * zAtoms;
                    xS[iSym][i + symIndex] += translate[0];
                    yS[iSym][i + symIndex] += translate[1];
                    zS[iSym][i + symIndex] += translate[2];
                }

                // Save CoM cartesian coordinates
                centerMolsCart[iSym * zPrime + zp] = centerOfMass;
            }
        }

        //Determine molecular distances to "center" of sphere.
        //  In PACCOM the center is the geometric average of coordinates.
        //  In FFX the center is the center of the replicates crystal.
        DoubleIndexPair[] molsDists = new DoubleIndexPair[numEntities];
        double[] cartCenter = new double[3];

        // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
        // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
        double[] fracCenter = {0.5, 0.5, 0.5};
        replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

        if (logger.isLoggable(Level.FINER)) {
            logger.finer(format(" Expanded Crystal Center: %16.8f %16.8f %16.8f",
                    cartCenter[0], cartCenter[1], cartCenter[2]));
        }

        for (int i = 0; i < numEntities; i++) {
            // Then compute Euclidean distance from Cartesian center of the replicates cell
            molsDists[i] = new DoubleIndexPair(i, dist(cartCenter, centerMolsCart[i]));
        }

        // Sort the molecules by their distance from the center.
        // Note that the smallest distances are first in the array after the sort.
        sort(molsDists);

        if (logger.isLoggable(Level.FINEST)) {
            logger.finest("\n Copy  SymOp        Distance");
        }
        double[] systemCoords = new double[numEntities * zAtoms * 3];
        for (int n = 0; n < nSymm; n++) {
            for(int zp = 0; zp < zPrime; zp++) {
                int index = n * zPrime + zp;
                // Current molecule
                int iSym = molsDists[index].getIndex();
                double distance = molsDists[index].getDoubleValue();
                if (logger.isLoggable(Level.FINEST)) {
                    logger.finest(format(" %4d  %5d  %16.8f", index, iSym, distance));
                }

                // Create a new set of Atoms for each SymOp molecule
                for (int i = 0; i < zAtoms; i++) {
                    int symIndex = index * zAtoms * 3;
                    int atomIndex = i * 3;
                    int conversion = iSym/zPrime;
                    int shift = (iSym % zPrime) * zAtoms;
                    systemCoords[symIndex + atomIndex] = xS[conversion][i + shift];
                    systemCoords[symIndex + atomIndex + 1] = yS[conversion][i + shift];
                    systemCoords[symIndex + atomIndex + 2] = zS[conversion][i + shift];
                }
            }
        }
        return systemCoords;
    }

    /**
     * Produce a translation vector necessary to move an object with the current center of mass (com)
     * into the provided crystal.
     *
     * @param crystal Replicates crystal within whom coordinates should be moved.
     * @param com     Center of mass (x, y, z) for the object of concern
     * @return double[] translation vector to move the object to within the provided crystal.
     */
    private static double[] moveIntoCrystal(Crystal crystal, double[] com) {

        double[] translate = new double[3];
        double[] currentCoM = new double[3];
        currentCoM[0] = com[0];
        currentCoM[1] = com[1];
        currentCoM[2] = com[2];

        // Move the COM to the Replicates Crystal.
        crystal.toFractionalCoordinates(com, translate);
        translate[0] = mod(translate[0], 1.0);
        translate[1] = mod(translate[1], 1.0);
        translate[2] = mod(translate[2], 1.0);
        crystal.toCartesianCoordinates(translate, translate);

        // Correct center of mass.
        com[0] = translate[0];
        com[1] = translate[1];
        com[2] = translate[2];

        // The translation vector is difference between the new location and the current COM.
        translate[0] -= currentCoM[0];
        translate[1] -= currentCoM[1];
        translate[2] -= currentCoM[2];

        if (logger.isLoggable(Level.FINEST)) {
            logger.finest(format(" Center of Mass Prior: %16.8f %16.8f %16.8f",
                    currentCoM[0], currentCoM[1], currentCoM[2]));
            logger.finest(format(" Center of Mass Post: %16.8f %16.8f %16.8f", com[0], com[1], com[2]));
        }

        return translate;
    }

    /**
     * Prioritize asymmetric units within the system based on distance to specified index.
     *
     * @param coordsXYZ      Coordinates for expanded crystal (should contain 3 * nAtoms * nMols entries).
     * @param mass           Mass of atoms within asymmetric unit (should contain one mass per atom in asym unit).
     * @param massSum        Sum of atomic masses within asymmetric unit.
     * @param centerOfMasses Center of masses for each replicate within inflated crystal.
     * @param nAtoms         Number of coordinates in an entity.
     * @param molDists       Prioritization of AUs in expanded system based on linkage criteria
     * @param index          Index of molecules to be center.
     * @param linkage        User specified criteria to determine prioritization.
     */
    private static void prioritizeReplicates(double[] coordsXYZ, double[] mass,
                                             double massSum, double[][] centerOfMasses, int nAtoms,
                                             DoubleIndexPair[] molDists, int index, int linkage) {
        // Find AU to be treated as the new center.
        // AUs added to system based on distance to center of all atoms. Index = 0 AU should be closest to all atom center.
        int nCoords = nAtoms * 3;
        int nMols = coordsXYZ.length / (nCoords);
        if (linkage == 0) {
            // Prioritize based on closest atomic distances.
            int centerIndex = index * nCoords;
            for (int i = 0; i < nMols; i++) {
                double tempDist = Double.MAX_VALUE;
                int molIndex = i * nCoords;
                for (int j = 0; j < nAtoms; j++) {
                    int centerAtomIndex = j * 3;
                    double[] centerXYZ = {coordsXYZ[centerIndex + centerAtomIndex],
                            coordsXYZ[centerIndex + centerAtomIndex + 1],
                            coordsXYZ[centerIndex + centerAtomIndex + 2]};
                    for (int k = 0; k < nAtoms; k++) {
                        int atomIndex = k * 3;
                        double[] xyz = {coordsXYZ[molIndex + atomIndex],
                                coordsXYZ[molIndex + atomIndex + 1],
                                coordsXYZ[molIndex + atomIndex + 2]};
                        double currDist = dist(centerXYZ, xyz);
                        if (currDist < tempDist) {
                            tempDist = currDist;
                        }
                        if (abs(tempDist) < MATCH_TOLERANCE) {
                            break;
                        }
                    }
                    if (abs(tempDist) < MATCH_TOLERANCE) {
                        break;
                    }
                }
                molDists[i] = new DoubleIndexPair(i, tempDist);
            }
            // Sort so the smallest distance is at position 0.
            Arrays.sort(molDists);

            if (logger.isLoggable(Level.FINEST)) {
                int numCheck = Math.min(5, molDists.length);
                double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
                centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
                for (int i = 0; i < numCheck; i++) {
                    logger.finest(format(" 1AU value %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                            molDists[i].getIndex(),
                            molDists[i].getDoubleValue(), targetMol[molDists[i].getIndex()][0],
                            targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
                }
            }
        } else if (linkage == 1) {
            // Prioritize based on geometric center/center of mass
            double[] coordCenter = centerOfMasses[index];
            for (int i = 0; i < nMols; i++) {
                double[] moleculeCenter = centerOfMasses[i];
                molDists[i] = new DoubleIndexPair(i, dist(coordCenter, moleculeCenter));
            }
            // Reorder based on distance to AU closest to Index.
            sort(molDists);

            if (logger.isLoggable(Level.FINEST)) {
                int numCheck = Math.min(5, molDists.length);
                double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
                centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
                for (int i = 0; i < numCheck; i++) {
                    logger.finest(format(" 1AU Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                            molDists[i].getIndex(),
                            molDists[i].getDoubleValue(), targetMol[molDists[i].getIndex()][0],
                            targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
                }
            }

            // Molecules in crystal sorted based on distance to center of mass of center most molecule
            // Want the first two molecules chosen in this manner, but third molecule to be closest to both
            // Assigning distance from molDists ensures correct ordering of center most molecule.
            DoubleIndexPair[] molDists2 = new DoubleIndexPair[nMols];
            molDists2[0] = new DoubleIndexPair(molDists[0].getIndex(), molDists[0].getDoubleValue());

            double[] avgCenter = new double[3];
            avgCenter[0] =
                    (centerOfMasses[molDists[0].getIndex()][0] + centerOfMasses[molDists[1].getIndex()][0]) / 2;
            avgCenter[1] =
                    (centerOfMasses[molDists[0].getIndex()][1] + centerOfMasses[molDists[1].getIndex()][1]) / 2;
            avgCenter[2] =
                    (centerOfMasses[molDists[0].getIndex()][2] + centerOfMasses[molDists[1].getIndex()][2]) / 2;

            for (int i = 1; i < nMols; i++) {
                double[] moleculeCenter = centerOfMasses[molDists[i].getIndex()];
                molDists2[i] = new DoubleIndexPair(molDists[i].getIndex(), dist(avgCenter, moleculeCenter));
            }
            sort(molDists2);
            if (logger.isLoggable(Level.FINEST)) {
                int numCheck = Math.min(5, molDists2.length);
                double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
                centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
                for (int i = 0; i < numCheck; i++) {
                    logger.finest(format(" 2AU Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                            molDists2[i].getIndex(),
                            molDists2[i].getDoubleValue(), targetMol[molDists2[i].getIndex()][0],
                            targetMol[molDists2[i].getIndex()][1], targetMol[molDists2[i].getIndex()][2]));
                }
            }
            // Molecules in crystal sorted based on distance to center of mass of center 2 most molecule
            // Want the first three molecules chosen in this manner, but rest to be closest to all three
            // Assigning distance from molDists2 ensures correct ordering of center most molecule.
            DoubleIndexPair[] molDists3 = new DoubleIndexPair[nMols];
            molDists3[0] = new DoubleIndexPair(molDists2[0].getIndex(), molDists2[0].getDoubleValue());
            molDists3[1] = new DoubleIndexPair(molDists2[1].getIndex(), molDists2[1].getDoubleValue());
            avgCenter = new double[3];
            for (int i = 0; i < 3; i++) {
                avgCenter[0] += centerOfMasses[molDists2[i].getIndex()][0];
                avgCenter[1] += centerOfMasses[molDists2[i].getIndex()][1];
                avgCenter[2] += centerOfMasses[molDists2[i].getIndex()][2];
            }
            for (int i = 0; i < 3; i++) {
                avgCenter[i] /= 3;
            }

            for (int i = 2; i < nMols; i++) {
                double[] moleculeCenter = centerOfMasses[molDists2[i].getIndex()];
                molDists3[i] = new DoubleIndexPair(molDists2[i].getIndex(), dist(avgCenter, moleculeCenter));
            }
            //Reorder based on center point between center-most AU to all atom center and closest AU to center-most AU.
            Arrays.sort(molDists3);
            arraycopy(molDists3, 0, molDists, 0, nMols);
            if (logger.isLoggable(Level.FINEST)) {
                int numCheck = Math.min(5, molDists3.length);
                double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
                centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
                for (int i = 0; i < numCheck; i++) {
                    logger.finest(format(" 3AU Rank %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                            molDists3[i].getIndex(),
                            molDists3[i].getDoubleValue(), targetMol[molDists3[i].getIndex()][0],
                            targetMol[molDists3[i].getIndex()][1], targetMol[molDists3[i].getIndex()][2]));
                }
            }
        } else if (linkage == 2) {
            // Prioritize based on minimum distance between the farthest atom.
            int centerIndex = index * nCoords;
            for (int i = 0; i < nMols; i++) {
                double tempDist = 0.0;
                int molIndex = i * nCoords;
                for (int j = 0; j < nAtoms; j++) {
                    int centerAtomIndex = j * 3;
                    double[] centerXYZ = {coordsXYZ[centerIndex + centerAtomIndex],
                            coordsXYZ[centerIndex + centerAtomIndex + 1],
                            coordsXYZ[centerIndex + centerAtomIndex + 2]};
                    for (int k = 0; k < nAtoms; k++) {
                        int atomIndex = k * 3;
                        double[] xyz = {coordsXYZ[molIndex + atomIndex],
                                coordsXYZ[molIndex + atomIndex + 1],
                                coordsXYZ[molIndex + atomIndex + 2]};
                        double currDist = dist(centerXYZ, xyz);
                        if (currDist > tempDist) {
                            tempDist = currDist;
                        }
                        if (abs(tempDist) < MATCH_TOLERANCE) {
                            break;
                        }
                    }
                    if (abs(tempDist) < MATCH_TOLERANCE) {
                        break;
                    }
                }
                molDists[i] = new DoubleIndexPair(i, tempDist);
            }
            // Sort so the smallest distance is at position 0.
            Arrays.sort(molDists);

            if (logger.isLoggable(Level.FINEST)) {
                int numCheck = Math.min(5, molDists.length);
                double[][] targetMol = new double[coordsXYZ.length / (nCoords)][3];
                centerOfMass(targetMol, coordsXYZ, mass, massSum, nAtoms);
                for (int i = 0; i < numCheck; i++) {
                    logger.finest(format(" 1AU value %d Target: %d Index: %d Dist %4.4f (%4.4f %4.4f %4.4f)", i, index,
                            molDists[i].getIndex(),
                            molDists[i].getDoubleValue(), targetMol[molDists[i].getIndex()][0],
                            targetMol[molDists[i].getIndex()][1], targetMol[molDists[i].getIndex()][2]));
                }
            }
        }
    }

    /**
     * Add a value to a list of doubles if its difference to all listed values is greater than the
     * tolerance.
     *
     * @param values List of values already found.
     * @param value  Potential new value if it is not already in list.
     * @param tol   Tolerance that determine whether values are equal.
     */
    private static boolean addLooseUnequal(List<Double> values, double value, double tol) {
        boolean found = false;
        for (Double dbl : values) {
            if (abs(dbl - value) < tol) {
                found = true;
                break;
            }
        }
        if (!found) {
            values.add(value);
        }
        // If value is not found it was added. Otherwise, not added.
        return !found;
    }

    /**
     * Save the provided coordinates as a PDB file.
     *
     * @param molecularAssembly Asymmetric unit that forms the crystal of interest.
     * @param coords            Coordinates to be saved within the PDB.
     * @param comparisonAtoms   Atoms of interest within the initial asymmetric unit.
     * @param description       Unique identifier that will be added to PDB file name.
     * @param compNum           Unique number for the current comparison
     * @param save              Type of file to save (0=none, 1=PDB, 2=XYZ)
     */
    private static void saveAssembly(MolecularAssembly molecularAssembly, double[] coords,
                                     int[] comparisonAtoms, String description, int compNum, int save) {
        //TODO: Save systems out as original molecule regardless of selection
        String fileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getName());
        File saveLocation;
        if (save == 2) {
            saveLocation = new File(fileName + description + ".xyz");
        } else {
            saveLocation = new File(fileName + description + ".pdb");
        }
        // Save aperiodic system of n_mol closest atoms for visualization
        MolecularAssembly currentAssembly = new MolecularAssembly(molecularAssembly.getName());
        List<Bond> bondList = molecularAssembly.getBondList();
        ArrayList<Atom> newAtomList = new ArrayList<>();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int atomIndex = 0;
        int compareAtomsSize = comparisonAtoms.length;
        int nCoords = compareAtomsSize * 3;
        int nAtoms = atoms.length;
        // Reset atom Index for indexing new atoms
        atomIndex = 1;
        int numMols = coords.length / (3 * compareAtomsSize);
        for (int n = 0; n < numMols; n++) {
            // Obtain atoms from moved AU (create copy to move to origin)
            // move original and copy to origin
            // rotate original to match copy
            // translate back to moved location
            // Add atoms from moved original to atom list.
            ArrayList<Atom> atomList = new ArrayList<>();
            // Create a new set of Atoms for each SymOp molecule
            int atomValue = 0;
            //Add atoms from comparison to output assembly.
            for (Integer i : comparisonAtoms) {
                Atom a = atoms[i];
                double[] xyz = new double[3];
                xyz[0] = coords[n * nCoords + atomValue];
                xyz[1] = coords[n * nCoords + atomValue + 1];
                xyz[2] = coords[n * nCoords + atomValue + 2];
                Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
                atomList.add(atom);
                atomValue += 3;
            }
            // Setup bonds for AUs in comparison.
            if(save == 2) {
                //TODO make more robust. Currently only handles all atom selection. Not subsets.
                for (Bond bond : bondList) {
                    Atom a1 = bond.getAtom(0);
                    Atom a2 = bond.getAtom(1);
                    //Indices stored as human-readable.
                    int a1Ind = a1.getIndex() - 1;
                    int a2Ind = a2.getIndex() - 1;
                    if (IntStream.of(comparisonAtoms).anyMatch(x -> x == a1Ind) && IntStream.of(comparisonAtoms).anyMatch(x -> x == a2Ind)) {
                            Atom newA1 = atomList.get(a1Ind);
                            Atom newA2 = atomList.get(a2Ind);
                            Bond b = new Bond(newA1, newA2);
                            b.setBondType(bond.getBondType());
                    }
                }
            }
            newAtomList.addAll(atomList);
        }

        // Construct the force field for the expanded set of molecules
        ForceField forceField = molecularAssembly.getForceField();

        // Clear all periodic boundary keywords to achieve an aperiodic system.
        forceField.clearProperty("a-axis");
        forceField.clearProperty("b-axis");
        forceField.clearProperty("c-axis");
        forceField.clearProperty("alpha");
        forceField.clearProperty("beta");
        forceField.clearProperty("gamma");
        forceField.clearProperty("spacegroup");
        currentAssembly.setForceField(forceField);

        // The biochemistry method is designed to load chemical entities into the
        // Polymer, Molecule, Water and Ion data structure.
        Utilities.biochemistry(currentAssembly, newAtomList);

        currentAssembly.setFile(molecularAssembly.getFile());
        if (save == 2) {
            File key = new File(fileName + ".key");
            File properties = new File(fileName + ".properties");
            if (key.exists()) {
                File keyComparison = new File(fileName + description + ".key");
                try {
                    if (keyComparison.createNewFile()) {
                        FileUtils.copyFile(key, keyComparison);
                    } else {
                        logger.info(" Could not create properties file.");
                    }
                } catch (Exception ex) {
                    // Likely using properties file.
                    logger.finest(ex.toString());
                }
            } else if (properties.exists()) {
                File propertiesComparison = new File(fileName + description + ".properties");
                try {
                    if (propertiesComparison.createNewFile()) {
                        FileUtils.copyFile(properties, propertiesComparison);
                    } else {
                        logger.info(" Could not create properties file.");
                    }
                } catch (Exception ex) {
                    // Likely not using a key/properties file... so PDB?
                    logger.info(" Neither key nor properties file detected therefore only creating XYZ.");
                    logger.finest(ex.toString());
                }
            }
            XYZFilter xyzFilter = new XYZFilter(saveLocation, currentAssembly, forceField, currentAssembly.getProperties());
            xyzFilter.writeFile(saveLocation, true);
        } else {
            PDBFilter pdbfilter = new PDBFilter(saveLocation, currentAssembly, forceField,
                    currentAssembly.getProperties());
            pdbfilter.setModelNumbering(compNum);
            pdbfilter.writeFile(saveLocation, true, false, false);
        }
    }

    /**
     * Save the provided coordinates as a PDB file with accompanying CSV containing RMSD.
     *
     * @param molecularAssembly Asymmetric unit that forms the crystal of interest.
     * @param coords            Coordinates to be saved within the PDB.
     * @param comparisonAtoms   Atoms of interest within the initial asymmetric unit.
     * @param description       Unique identifier that will be added to PDB file name.
     * @param finalRMSD         RMSD to be saved to CSV file.
     * @param compNum           Unique number for the current comparison
     * @param save              Type of file to save (0=none, 1=PDB, 2=XYZ)
     */
    private static void saveAssembly(MolecularAssembly molecularAssembly, double[] coords,
                                     int[] comparisonAtoms, String description, double finalRMSD, int compNum, int save) {
        saveAssembly(molecularAssembly, coords, comparisonAtoms, description, compNum, save);
        String fileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getName());
        try {
            // Needs same name as PDB so follow PDB format.
            File csv = PDBFilter.version(new File(fileName + description + ".csv"));
            if (csv.createNewFile()) {
                BufferedWriter bw = new BufferedWriter(new FileWriter(csv, false));
                bw.append("rms\n");
                bw.append(format("%4.4f\n", finalRMSD));
                bw.close();
            } else {
                logger.warning(" Could not create CSV file \"" + csv.getName() + "\"");
            }
        } catch (Exception ex) {
            logger.info(ex.toString());
        }
    }
}

//        // Used in QEtoXYZ.groovy which is not ready for git which is why this method appears unused.
//    /**
//     * Orient coordinates so that the second index is on the x axis, and the third index is on the X-Y
//     * plane. First index should be at the origin (0, 0, 0).
//     *
//     * @param coordsXYZ   An array of XYZ positions (e.g. [x0, y0, z0, x1, y1, z1, x2, y2, z2]
//     * @param atomIndices Indices for three desired sets from the XYZ list (e.g. [0, 1, 2]). Index
//     *                    0 should be at origin.
//     */
//    public static void standardOrientation(double[] coordsXYZ, int[] atomIndices) {
//        int numCoords = coordsXYZ.length / 3;
//        double[] atomCoords = new double[3 * 3];
//        atomCoords[0] = coordsXYZ[atomIndices[0]];
//        atomCoords[1] = coordsXYZ[atomIndices[0] + 1];
//        atomCoords[2] = coordsXYZ[atomIndices[0] + 2];
//        atomCoords[3] = coordsXYZ[atomIndices[1]];
//        atomCoords[4] = coordsXYZ[atomIndices[1] + 1];
//        atomCoords[5] = coordsXYZ[atomIndices[1] + 2];
//        atomCoords[6] = coordsXYZ[atomIndices[2]];
//        atomCoords[7] = coordsXYZ[atomIndices[2] + 1];
//        atomCoords[8] = coordsXYZ[atomIndices[2] + 2];
//
//        // TODO: Delete coordsXYZOrig?
//        double[] coordsXYZOrig = new double[numCoords * 3];
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            coordsXYZOrig[atomIndex] = coordsXYZ[atomIndex];
//            coordsXYZOrig[atomIndex + 1] = coordsXYZ[atomIndex + 1];
//            coordsXYZOrig[atomIndex + 2] = coordsXYZ[atomIndex + 2];
//        }
//
//        // TODO: Delete atomsCoordsOrig?
//        double[] atomsCoordsOrig = new double[3 * 3];
//        arraycopy(atomCoords, 0, atomsCoordsOrig, 0, 9);
//
//        logger.fine(
//                format(" START: N1:\t%16.15f %16.15f %16.15f", atomCoords[0], atomCoords[1], atomCoords[2]));
//        logger.fine(
//                format(" START: N2:\t%16.15f %16.15f %16.15f", atomCoords[3], atomCoords[4], atomCoords[5]));
//        logger.fine(
//                format(" START: N3:\t%16.15f %16.15f %16.15f", atomCoords[6], atomCoords[7], atomCoords[8]));
//
//        double p1n2 = coordsXYZ[atomIndices[1]];
//        double q1n2 = coordsXYZ[atomIndices[1] + 1];
//        double r1n2 = coordsXYZ[atomIndices[1] + 2];
//
//        // Calculation of sigma, phai, and cita angles needed to get specified atoms to desired loci
//        double cita0 = acos(p1n2 / sqrt(p1n2 * p1n2 + q1n2 * q1n2));
//        double phai0 = acos(sqrt(p1n2 * p1n2 + q1n2 * q1n2) /
//                sqrt(p1n2 * p1n2 + q1n2 * q1n2 + r1n2 * r1n2));
//        if (q1n2 < 0.0) {
//            cita0 = -cita0;
//        }
//
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            double ptmp = coordsXYZ[atomIndex] * cos(cita0) + coordsXYZ[atomIndex + 1] * sin(cita0);
//            double qtmp = -coordsXYZ[atomIndex] * sin(cita0) + coordsXYZ[atomIndex + 1] * cos(cita0);
//            coordsXYZ[atomIndex] = ptmp;
//            coordsXYZ[atomIndex + 1] = qtmp;
//        }
//
//        p1n2 = coordsXYZ[atomIndices[1]];
//        q1n2 = coordsXYZ[atomIndices[1] + 1];
//
//        if (r1n2 > 0.0) {
//            phai0 = -phai0;
//        }
//
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            double ptmp = coordsXYZ[atomIndex] * cos(phai0) - coordsXYZ[atomIndex + 2] * sin(phai0);
//            double rtmp = coordsXYZ[atomIndex] * sin(phai0) + coordsXYZ[atomIndex + 2] * cos(phai0);
//            coordsXYZ[atomIndex] = ptmp;
//            coordsXYZ[atomIndex + 2] = rtmp;
//        }
//
//        p1n2 = coordsXYZ[atomIndices[1]];
//        r1n2 = coordsXYZ[atomIndices[1] + 2];
//        double p1n3 = coordsXYZ[atomIndices[2]];
//        double q1n3 = coordsXYZ[atomIndices[2] + 1];
//        double r1n3 = coordsXYZ[atomIndices[2] + 2];
//
//        double sigma0 = acos(q1n3 / sqrt(q1n3 * q1n3 + r1n3 * r1n3));
//        if (r1n3 < 0.0) {
//            sigma0 = -sigma0;
//        }
//
//        for (int i = 0; i < numCoords; i++) {
//            int atomIndex = i * 3;
//            double qtmp = coordsXYZ[atomIndex + 1] * cos(sigma0) + coordsXYZ[atomIndex + 2] * sin(sigma0);
//            double rtmp = -coordsXYZ[atomIndex + 1] * sin(sigma0) + coordsXYZ[atomIndex + 2] * cos(sigma0);
//            coordsXYZ[atomIndex + 1] = qtmp;
//            coordsXYZ[atomIndex + 2] = rtmp;
//        }
//
//        q1n2 = coordsXYZ[atomIndices[1] + 1];
//        r1n2 = coordsXYZ[atomIndices[1] + 2];
//        q1n3 = coordsXYZ[atomIndices[2] + 1];
//        r1n3 = coordsXYZ[atomIndices[2] + 2];
//
//        if (logger.isLoggable(Level.FINER)) {
//            logger.finer(
//                    format(" DONE N1: %16.15f %16.15f %16.15f", atomCoords[0], atomCoords[1], atomCoords[2]));
//            logger.finer(format(" DONE N2: %16.15f %16.15f %16.15f", p1n2, q1n2, r1n2));
//            logger.finer(format(" DONE N3: %16.15f %16.15f %16.15f", p1n3, q1n3, r1n3));
//        }
//    }

//  Frank-Kasper phases of metallic ions can reach a coordination number of 16...