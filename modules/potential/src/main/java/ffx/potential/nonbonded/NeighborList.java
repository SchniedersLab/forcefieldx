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
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.util.Range;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;

/**
 * The NeighborList class builds Verlet lists in parallel via a spatial
 * decomposition.
 * <br>
 * <ol>
 * <li>
 * The unit cell is partitioned into <code>nA * nB * nC</code> smaller
 * axis-aligned cells, where {nA, nB, nC} are chosen as large as possible
 * subject to the criteria that the length of each side of a sub-volume (rCellA,
 * rCellB, rCellC) multiplied by (nEdgeA, nEdgeB, nEdgeC), respectively, must be
 * greater than the cutoff distance <code>Rcut</code> plus a buffer distance
 * <code>delta</code>:
 * <br>
 * <code>rCellA * nEdgeA .GE. (Rcut + delta)</code>
 * <br>
 * <code>rCellB * nEdgeB .GE. (Rcut + delta)</code>
 * <br>
 * <code>rCellC * nEdgeC .GE. (Rcut + delta)</code>
 * <br>
 * All neighbors of an atom are in a block of
 * (2*nEdgeA+1)(2*nEdgeB+1)(2*nEdgeC+1) neighborCells.
 * </li>
 * <li>
 * Interactions between an atom and neighbors in the asymmetric unit require
 * only half the neighboring cells to be searched to avoid double counting.
 * However, enumeration of interactions between an atom in the asymmetric unit
 * and its neighbors in a symmetry mate require all cells to be searched.
 * </li>
 * <li>
 * Verlet lists from the search are stored, which reduces the number of
 * neighbors whose distances must be calculated by a factor of approximately:
 * <br>
 * <code>(4/3*Pi*Rcut^3)/(neighborCells*Vcell)</code>
 * About 1/3 as many interactions are contained in the Verlet lists as in the neighboring
 * cells.
 * </li>
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class NeighborList extends ParallelRegion {

    /**
     * The logger.
     */
    private static final Logger logger = Logger.getLogger(NeighborList.class.getName());
    /**
     * The crystal object defines the unit cell dimensions and space group.
     */
    private Crystal crystal;
    /**
     * The number of asymmetric units in the unit cell.
     */
    private int nSymm;
    /**
     * The number of atoms in the asymmetric unit.
     */
    private int nAtoms;
    /**
     * The masking rules to apply when building the neighbor list.
     */
    private final MaskingInterface maskingRules;
    /**
     * Reduced coordinates for each symmetry copy. [nSymm][3*nAtoms]
     */
    private double[][] coordinates;
    /**
     * The reduced coordinates of the asymmetric unit when the list was last
     * rebuilt.
     */
    private double[] previous;
    /**
     * The Verlet lists. [nSymm][nAtoms][nNeighbors]
     */
    private int[][][] lists;
    /**
     * Number of interactions per atom.
     */
    private int[] listCount;
    /**
     * Total number of interactions.
     */
    private final SharedInteger sharedCount;
    /**
     * Optimal pairwise ranges.
     */
    private final Range[] ranges;
    /**
     * Pairwise ranges for load balancing.
     */
    private PairwiseSchedule pairwiseSchedule;
    /**
     * Number of interactions between atoms in the asymmetric unit.
     */
    private int asymmetricUnitCount;
    /**
     * Number of interactions between atoms in the asymmetric unit and atoms in
     * symmetry mates.
     */
    private int symmetryMateCount;
    /**
     * Number of atoms in the asymmetric unit with interactions lists greater
     * than 0.
     */
    private int atomsWithIteractions;
    /**
     * The number of subcells that must be searched along the a-axis to find all
     * neighbors within the cutoff + buffer distance.
     * <p>
     * If each nEdgeX == 1 for (X=A,B,C} then all neighbors will be found in
     * 3x3x3 = 27 cells. If each nEdgeX == 2, then all neighbors will be found
     * in 5x5x5 = 125 cells (in this case the cells are smaller).
     */
    private int nEdgeA;
    /**
     * The number of subcells that must be searched along the b-axis to find all
     * neighbors within the cutoff + buffer distance.
     */
    private int nEdgeB;
    /**
     * The number of subcells that must be searched along the c-axis to find all
     * neighbors within the cutoff + buffer distance.
     */
    private int nEdgeC;
    /**
     * The number of divisions along the A-axis.
     */
    private int nA;
    /**
     * The number of divisions along the B-axis.
     */
    private int nB;
    /**
     * The number of divisions along the C-Axis.
     */
    private int nC;
    /**
     * The number of cells in one plane (nDivisions^2).
     */
    private int nAB;
    /**
     * The number of cells (nDivisions^3).
     */
    private int nCells;
    /**
     * A temporary array that holds the index of the cell each atom is assigned
     * to. [nSymm][nAtoms]
     */
    private int[][] cellIndex;
    /**
     * The cell indices of each atom along a A-axis.
     */
    private int[] cellA;
    /**
     * The cell indices of each atom along a B-axis.
     */
    private int[] cellB;
    /**
     * The cell indices of each atom along a C-axis.
     */
    private int[] cellC;
    /**
     * The list of atoms in each cell. [nSymm][nAtoms] = atom index
     */
    private int[][] cellList;
    /**
     * The offset of each atom from the start of the cell. The first atom atom
     * in the cell has 0 offset. [nSymm][nAtom] = offset of the atom
     */
    private int[][] cellOffset;
    /**
     * The number of atoms in each cell. [nSymm][nCell]
     */
    private int[][] cellCount;
    /**
     * The index of the first atom in each cell. [nSymm][nCell]
     */
    private int[][] cellStart;
    /**
     * The cutoff beyond which the pairwise energy is zero.
     */
    private final double cutoff;
    /**
     * A buffer, which is added to the cutoff distance, such that the Verlet
     * lists do not need to be calculated for all coordinate changes.
     */
    private final double buffer;
    /**
     * The maximum squared displacement allowed before list rebuild.
     */
    private final double motion2;
    /**
     * The sum of the cutoff + buffer.
     */
    private final double cutoffPlusBuffer;
    /**
     * Total^2 for distance comparisons without taking a sqrt.
     */
    private final double cutoffPlusBuffer2;
    /**
     * The array of fractional "a", "b", and "c" coordinates.
     */
    private double[] frac;
    /**
     * Atoms being used list.
     */
    private boolean[] use;
    /**
     * List of atoms.
     */
    private Atom[] atoms;

    // *************************************************************************
    // Parallel variables.
    /**
     * The ParallelTeam coordinates use of threads and their schedules.
     */
    private final ParallelTeam parallelTeam;
    /**
     * Number of threads used by the parallelTeam.
     */
    private final int threadCount;
    /**
     * A Verlet list loop for each thread.
     */
    private final NeighborListLoop[] verletListLoop;
    private long time;
    /**
     * Include intermolecular interactions.
     */
    private boolean intermolecular = true;
    /**
     * Molecule number for each atom.
     */
    private int[] molecules = null;
    /**
     * If true, interactions between two inactive atoms are included in the Neighborlist.
     * Set to true to match OpenMM behavior (i.e. interactions between two atoms with zero mass are included).
     * Set to false for more efficient "pure Java" optimizations.
     */
    private boolean inactiveInteractions = true;
    /**
     * Disable updates to the NeighborList; use with caution.
     */
    private boolean disableUpdates = false;

    /**
     * Constructor for the NeighborList class.
     *
     * @param maskingRules This parameter may be null.
     * @param crystal      Definition of the unit cell and space group.
     * @param atoms        The atoms to generate Verlet lists for.
     * @param cutoff       The cutoff distance.
     * @param buffer       The buffer distance.
     * @param parallelTeam Specifies the parallel environment.
     * @since 1.0
     */
    public NeighborList(MaskingInterface maskingRules, Crystal crystal,
                        Atom[] atoms, double cutoff, double buffer,
                        ParallelTeam parallelTeam) {
        this.maskingRules = maskingRules;
        this.crystal = crystal;
        this.cutoff = cutoff;
        this.buffer = buffer;
        this.parallelTeam = new ParallelTeam(parallelTeam.getThreadCount());
        this.atoms = atoms;
        nAtoms = atoms.length;

        // Configure the neighbor cutoff and list rebuilding criteria.
        cutoffPlusBuffer = cutoff + buffer;
        cutoffPlusBuffer2 = cutoffPlusBuffer * cutoffPlusBuffer;
        motion2 = (buffer / 2.0) * (buffer / 2.0);

        // Initialize parallel constructs.
        threadCount = parallelTeam.getThreadCount();
        sharedCount = new SharedInteger();
        ranges = new Range[threadCount];

        verletListLoop = new NeighborListLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            verletListLoop[i] = new NeighborListLoop();
        }

        // Initialize the neighbor list builder subcells.
        boolean print = logger.isLoggable(Level.FINE);
        initNeighborList(print);
    }

    /**
     * <p>Setter for the field <code>intermolecular</code>.</p>
     *
     * @param intermolecular a boolean.
     * @param molecules      an array of {@link int} objects.
     */
    public void setIntermolecular(boolean intermolecular, int[] molecules) {
        this.intermolecular = intermolecular;
        this.molecules = molecules;
    }

    /**
     * If disableUpdates true, disable updating the neighbor list upon motion.
     * Use with caution; best recommendation is to only use if all atoms have a
     * coordinate restraint.
     *
     * @param disableUpdate Disable updating the neighbor list
     */
    public void setDisableUpdates(boolean disableUpdate) {
        this.disableUpdates = disableUpdate;
    }

    /**
     * Returns the cutoff distance used internally by NeighborList.
     *
     * @return Cutoff distance in Angstroms.
     */
    public double getCutoff() {
        return cutoff;
    }

    /**
     * Getter for the disableUpdates field.
     *
     * @return If neighbor list updates disabled
     */
    public boolean getDisableUpdates() {
        return disableUpdates;
    }

    private void initNeighborList(boolean print) {

        // Allocate memory for fractional coordinates and subcell pointers for each atom.
        if (cellA == null || cellA.length < nAtoms) {
            cellA = new int[nAtoms];
            cellB = new int[nAtoms];
            cellC = new int[nAtoms];
            frac = new double[3 * nAtoms];
            previous = new double[3 * nAtoms];
            listCount = new int[nAtoms];
            pairwiseSchedule = new PairwiseSchedule(threadCount, nAtoms, ranges);
        } else {
            pairwiseSchedule.setAtoms(nAtoms);
        }

        // Set the number of symmetry operators.
        int newNSymm = crystal.spaceGroup.symOps.size();
        if (nSymm != newNSymm) {
            nSymm = newNSymm;
        }

        // Find the largest sphere that is enclosed by the unit cell.
        final double sphere = min(min(crystal.interfacialRadiusA,
                crystal.interfacialRadiusB), crystal.interfacialRadiusC);

        /*
          Assert that the boundary conditions defined by the crystal allow use
          of the minimum image condition.
         */
        if (!crystal.aperiodic()) {
            assert (sphere >= cutoffPlusBuffer);
        }

        // nEdgeA, nEdgeB and nEdgeC must be >= 1.
        nEdgeA = 2;
        nEdgeB = nEdgeA;
        nEdgeC = nEdgeA;

        /*
          totalA is the smallest a-axis for a subcell, given nEdgeA subcells
          along the a-axis, that assures all neighbors will be found.
         */
        double totalA = cutoffPlusBuffer / (double) nEdgeA;
        double totalB = cutoffPlusBuffer / (double) nEdgeB;
        double totalC = cutoffPlusBuffer / (double) nEdgeC;
        nA = (int) floor(crystal.interfacialRadiusA / totalA);
        nB = (int) floor(crystal.interfacialRadiusB / totalB);
        nC = (int) floor(crystal.interfacialRadiusC / totalC);
        if (nA < nEdgeA * 2 + 1) {
            nA = 1;
        }
        if (nB < nEdgeB * 2 + 1) {
            nB = 1;
        }
        if (nC < nEdgeC * 2 + 1) {
            nC = 1;
        }
        nAB = nA * nB;
        nCells = nAB * nC;

        if (print) {
            StringBuilder sb = new StringBuilder("  Neighbor List Builder\n");
            sb.append(format("   Sub-volumes:                        %8d\n", nCells));
            sb.append(format("   Average Atoms Per Sub-Volume:       %8d",
                    nAtoms * nSymm / nCells));
            logger.info(sb.toString());
        }

        // Allocate memory, if necessary.
        if (cellList == null) {
            cellList = new int[nSymm][nAtoms];
            cellIndex = new int[nSymm][nAtoms];
            cellOffset = new int[nSymm][nAtoms];
        } else if (cellList.length < nSymm) {
            logger.info(String.format("  Neighbor-List: Increasing memory for nSymm (%d -> %d)",
                    cellList.length, nSymm));
            int nAtomMax = cellList[0].length;
            cellList = new int[nSymm][nAtomMax];
            cellIndex = new int[nSymm][nAtomMax];
            cellOffset = new int[nSymm][nAtomMax];
        } else if (cellList[0].length < nAtoms) {
            int nSymmMax = cellList.length;
            cellList = new int[nSymmMax][nAtoms];
            cellIndex = new int[nSymmMax][nAtoms];
            cellOffset = new int[nSymmMax][nAtoms];
        }

        if (cellStart == null) {
            cellStart = new int[nSymm][nCells];
            cellCount = new int[nSymm][nCells];
        } else if (cellStart.length < nSymm) {
            int maxCells = Math.max(nCells, cellStart[0].length);
            cellStart = new int[nSymm][maxCells];
            cellCount = new int[nSymm][maxCells];
        } else if (cellStart[0].length < nCells) {
            logger.info(String.format("  Neighbor-List: Increasing memory for nCell (%d -> %d)",
                    cellStart[0].length, nCells));
            for (int i = 0; i < cellStart.length; i++) {
                cellStart[i] = new int[nCells];
                cellCount[i] = new int[nCells];
            }
        }
    }

    /**
     * The NeighborList will be re-configured, if necessary, for the supplied
     * Crystal. Changes to both unit cell parameters and number of symmetry
     * operators are also acceptable.
     *
     * @param crystal A crystal defining boundary conditions and symmetry.
     */
    public void setCrystal(Crystal crystal) {
        this.crystal = crystal;
        initNeighborList(false);
    }

    /**
     * The NeighborList will be re-configured, if necessary, for the supplied
     * atom list.
     *
     * @param atoms A new list of atoms to operate on.
     */
    public void setAtoms(Atom[] atoms) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        initNeighborList(false);
    }

    /**
     * This method can be called as necessary to build/rebuild the neighbor
     * lists.
     *
     * @param coordinates  The coordinates of each atom [nSymm][nAtoms*3].
     * @param lists        The neighbor lists [nSymm][nAtoms][nPairs].
     * @param forceRebuild If true, the list is rebuilt even if no atom has
     *                     moved half the buffer size.
     * @param use          an array of boolean.
     * @param print        a boolean.
     * @since 1.0
     */
    public void buildList(final double[][] coordinates, final int[][][] lists,
                          boolean[] use, boolean forceRebuild, boolean print) {
        if (disableUpdates) {
            return;
        }
        this.coordinates = coordinates;
        this.lists = lists;
        this.use = use;
        if (forceRebuild || motion()) {

            // Save the current coordinates.
            double[] current = coordinates[0];
            for (int i = 0; i < nAtoms; i++) {
                int i3 = i * 3;
                int iX = i3 + XX;
                int iY = i3 + YY;
                int iZ = i3 + ZZ;
                previous[iX] = current[iX];
                previous[iY] = current[iY];
                previous[iZ] = current[iZ];
            }

            assignAtomsToCells();
            createNeighborList();
            if (print) {
                print();
            }

            pairwiseSchedule.updateRanges(sharedCount.get(), atomsWithIteractions, listCount);
        }
    }

    /**
     * Return the Verlet list.
     *
     * @return The Verlet list of size [nSymm][nAtoms][nNeighbors].
     */
    public int[][][] getNeighborList() {
        return lists;
    }

    /**
     * Takes a list of Atoms and obtains the local indices for these atoms.
     *
     * @param atomsToIndex Atoms to index
     * @return Indices thereof
     */
    private List<Integer> atomsToIndices(List<Atom> atomsToIndex) {
        Set<Atom> indexedAtoms = new HashSet<>(atomsToIndex);
        List<Integer> indices = new ArrayList<>(atomsToIndex.size());
        for (int i = 0; i < atoms.length; i++) {
            if (indexedAtoms.contains(atoms[i])) {
                indices.add(i);
            }
        }
        return indices;
    }

    /**
     * Returns the indices of atoms neighboring the passed set of atom indices.
     * Returned Set is exclusive of the passed-in indices.
     *
     * @param atomIndices Atom indices
     * @param maxDist     Maximum distance to consider
     * @return Set of neighboring atoms indices
     */
    private Set<Integer> getNeighborIndices(List<Integer> atomIndices, double maxDist) {
        double md2 = maxDist * maxDist;
        Set<Integer> neighbors = new HashSet<>();
        for (int i : atomIndices) {
            double[] xyzI = new double[3];
            atoms[i].getXYZ(xyzI);
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                int[] listi = lists[iSymm][i];
                for (int indexJ : listi) {
                    if (atomIndices.contains(indexJ) || neighbors.contains(indexJ)) {
                        continue;
                    }
                    double[] xyzJ = new double[3];
                    atoms[indexJ].getXYZ(xyzJ);
                    crystal.applySymOp(xyzJ, xyzJ, symOp);
                    for (int k = 0; k < 3; k++) {
                        xyzJ[k] -= xyzI[k];
                    }
                    if (crystal.image(xyzJ) <= md2) {
                        neighbors.add(indexJ);
                    }
                }
            }
        }

        return neighbors;
    }

    /**
     * Returns a set of Atoms neighboring those passed in. If their indices are
     * already available, preferentially use getNeighborIndices instead. Is
     * exclusive of passed atoms.
     *
     * @param atomList Atoms to find neighbors of
     * @param maxDist  Maximum distance to consider a neighbor
     * @return Set of neighboring Atoms
     */
    public Set<Atom> getNeighborAtoms(List<Atom> atomList, double maxDist) {
        Set<Integer> atomIndices = getNeighborIndices(atomsToIndices(atomList), maxDist);
        Set<Atom> atomSet = new HashSet<>();
        for (Integer ai : atomIndices) {
            atomSet.add(atoms[ai]);
        }
        return atomSet;
    }

    /**
     * <p>
     * Getter for the field <code>pairwiseSchedule</code>.</p>
     *
     * @return a {@link ffx.potential.nonbonded.PairwiseSchedule} object.
     */
    public PairwiseSchedule getPairwiseSchedule() {
        return pairwiseSchedule;
    }

    private void print() {
        StringBuilder sb = new StringBuilder(format("   Buffer:                                %5.2f (A)\n", buffer));
        sb.append(format("   Cut-off:                               %5.2f (A)\n", cutoff));
        sb.append(format("   Total:                                 %5.2f (A)\n", cutoffPlusBuffer));
        sb.append(format("   Neighbors in the asymmetric unit:%11d\n", asymmetricUnitCount));
        if (nSymm > 1) {
            int num = (asymmetricUnitCount + symmetryMateCount) * nSymm;
            sb.append(format("   Neighbors in symmetry mates:    %12d\n", symmetryMateCount));
            sb.append(format("   Neighbors in the unit cell:     %12d\n", num));
        }
        logger.info(sb.toString());
    }

    /**
     * Assign asymmetric and symmetry mate atoms to cells.
     * <p>
     * This is very fast; there is little to be gained from parallelizing it at
     * this point.
     *
     * @since 1.0
     */
    private void assignAtomsToCells() {
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int[] cellIndexs = cellIndex[iSymm];
            final int[] cellCounts = cellCount[iSymm];
            final int[] cellStarts = cellStart[iSymm];
            final int[] cellLists = cellList[iSymm];
            final int[] cellOffsets = cellOffset[iSymm];
            // Zero out the cell counts.
            for (int i = 0; i < nCells; i++) {
                cellCounts[i] = 0;
            }
            // Convert to fractional coordinates.
            final double[] xyz = coordinates[iSymm];
            crystal.toFractionalCoordinates(nAtoms, xyz, frac);
            // Assign each atom to a cell using fractional coordinates.
            for (int i = 0; i < nAtoms; i++) {
                int i3 = i * 3;
                double xu = frac[i3 + XX];
                double yu = frac[i3 + YY];
                double zu = frac[i3 + ZZ];
                // Move the atom into the range 0.0 <= x < 1.0
                while (xu < 0.0) {
                    xu += 1.0;
                }
                while (xu >= 1.0) {
                    xu -= 1.0;
                }
                while (yu < 0.0) {
                    yu += 1.0;
                }
                while (yu >= 1.0) {
                    yu -= 1.0;
                }
                while (zu < 0.0) {
                    zu += 1.0;
                }
                while (zu >= 1.0) {
                    zu -= 1.0;
                }
                // The cell indices of this atom.
                final int a = (int) floor(xu * nA);
                final int b = (int) floor(yu * nB);
                final int c = (int) floor(zu * nC);
                if (iSymm == 0) {
                    cellA[i] = a;
                    cellB[i] = b;
                    cellC[i] = c;
                }
                // The cell index of this atom.
                final int index = a + b * nA + c * nAB;
                cellIndexs[i] = index;
                // The offset of this atom from the beginning of the cell.
                cellOffsets[i] = cellCounts[index]++;
            }
            // Define the starting indices.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }
            // Move atom locations into a list ordered by cell.
            for (int i = 0; i < nAtoms; i++) {
                final int index = cellIndexs[i];
                cellLists[cellStarts[index]++] = i;
            }
            // Define the starting indices again.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }
        }
    }

    /**
     * Execute the parallel Verlet list builder.
     *
     * @since 1.0
     */
    private void createNeighborList() {
        asymmetricUnitCount = 0;
        symmetryMateCount = 0;
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            if (lists[iSymm] == null || lists[iSymm].length < nAtoms) {
                lists[iSymm] = new int[nAtoms][];
            }
        }
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = "Fatal exception building neighbor list.\n";
            logger.log(Level.SEVERE, message, e);
        }
        int[][] list = lists[0];
        atomsWithIteractions = 0;
        for (int i = 0; i < nAtoms; i++) {
            asymmetricUnitCount += list[i].length;
            if (listCount[i] > 0) {
                atomsWithIteractions++;
            }
        }
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
            list = lists[iSymm];
            for (int i = 0; i < nAtoms; i++) {
                symmetryMateCount += list[i].length;
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 1.0
     */
    @Override
    public void start() {
        time = System.nanoTime();
        sharedCount.set(0);
    }

    /**
     * {@inheritDoc}
     * <p>
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 1.0
     */
    @Override
    public void run() {
        try {
            execute(0, nAtoms - 1, verletListLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception building neighbor list in thread: " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * This is method should not be called; it is invoked by Parallel Java.
     *
     * @since 1.0
     */
    @Override
    public void finish() {
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format("   Parallel Neighbor List:           %10.3f sec",
                    (System.nanoTime() - time) * 1e-9));
        }
    }

    /**
     * Detect if any atom has moved 1/2 the buffer size.
     *
     * @return True if an atom has moved 1/2 the buffer size.
     * @since 1.0
     */
    private boolean motion() {
        double[] current = coordinates[0];
        for (int i = 0; i < nAtoms; i++) {
            int i3 = i * 3;
            int iX = i3 + XX;
            int iY = i3 + YY;
            int iZ = i3 + ZZ;
            double dx = previous[iX] - current[iX];
            double dy = previous[iY] - current[iY];
            double dz = previous[iZ] - current[iZ];
            double dr2 = crystal.image(dx, dy, dz);
            if (dr2 > motion2) {
                return true;
            }
        }
        return false;
    }

    /**
     * The VerletListLoop class encapsulates thread local variables and methods
     * for building Verlet lists based on a spatial decomposition of the unit
     * cell.
     *
     * @author Michael J. Schnieders
     * @since 1.0
     */
    private class NeighborListLoop extends IntegerForLoop {

        private int n;
        private int iSymm;
        private int atomIndex;
        private boolean iactive = true;
        private int count;
        private int[] asymmetricIndex;
        private double[] xyz;
        private int[] pairs;
        private double[] mask;
        private boolean[] vdw14;
        private final IntegerSchedule schedule;
        // Extra padding to avert cache interference.
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        NeighborListLoop() {
            int len = 1000;
            pairs = new int[len];
            schedule = IntegerSchedule.dynamic(10);
        }

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void start() {
            xyz = coordinates[0];
            count = 0;
            if (mask == null || mask.length < nAtoms) {
                mask = new double[nAtoms];
                fill(mask, 1.0);
                vdw14 = new boolean[nAtoms];
                fill(vdw14, false);
            }
        }

        @Override
        public void finish() {
            sharedCount.addAndGet(count);
        }

        @Override
        public void run(final int lb, final int ub) {
            asymmetricIndex = cellIndex[0];
            for (iSymm = 0; iSymm < nSymm; iSymm++) {
                int[][] list = lists[iSymm];
                // Loop over all atoms.
                for (atomIndex = lb; atomIndex <= ub; atomIndex++) {
                    n = 0;

                    if (iSymm == 0) {
                        listCount[atomIndex] = 0;
                    }

                    if (use == null || use[atomIndex]) {

                        iactive = atoms[atomIndex].isActive();

                        final int a = cellA[atomIndex];
                        final int b = cellB[atomIndex];
                        final int c = cellC[atomIndex];

                        final int index = a + b * nA + c * nAB;

                        int a1 = a + 1;
                        int aStart = a - nEdgeA;
                        int aStop = a + nEdgeA;
                        int b1 = b + 1;
                        int bStart = b - nEdgeB;
                        int bStop = b + nEdgeB;
                        int c1 = c + 1;
                        int cStart = c - nEdgeC;
                        int cStop = c + nEdgeC;

                        /*
                          If the number of divisions is 1 in any direction then
                          set the loop limits to the current cell value.
                         */
                        if (nA == 1) {
                            aStart = a;
                            aStop = a;
                        }
                        if (nB == 1) {
                            bStart = b;
                            bStop = b;
                        }
                        if (nC == 1) {
                            cStart = c;
                            cStop = c;
                        }

                        if (iSymm == 0) {
                            // Interactions within the "self-volume".
                            atomCellPairs(index);

                            // Half of the neighboring volumes are searched to avoid double counting.

                            // (a, b+1..b+nE, c)
                            for (int bi = b1; bi <= bStop; bi++) {
                                atomCellPairs(image(a, bi, c));
                            }
                            // (a, b-nE..b+nE, c+1..c+nE)
                            for (int bi = bStart; bi <= bStop; bi++) {
                                for (int ci = c1; ci <= cStop; ci++) {
                                    atomCellPairs(image(a, bi, ci));
                                }
                            }
                            // (a+1..a+nE, b-nE..b+nE, c-nE..c+nE)
                            for (int bi = bStart; bi <= bStop; bi++) {
                                for (int ci = cStart; ci <= cStop; ci++) {
                                    for (int ai = a1; ai <= aStop; ai++) {
                                        atomCellPairs(image(ai, bi, ci));
                                    }
                                }
                            }
                        } else {

                            // Interactions with all adjacent symmetry mate cells.
                            for (int ai = aStart; ai <= aStop; ai++) {
                                for (int bi = bStart; bi <= bStop; bi++) {
                                    for (int ci = cStart; ci <= cStop; ci++) {
                                        atomCellPairs(image(ai, bi, ci));
                                    }
                                }
                            }
                        }
                    }

                    list[atomIndex] = new int[n];
                    listCount[atomIndex] += n;
                    count += n;
                    arraycopy(pairs, 0, list[atomIndex], 0, n);
                }
            }
        }

        /**
         * If the index is >= to nX, it is mapped back into the periodic unit
         * cell by subtracting nX. If the index is less than 0, it is mapped
         * into the periodic unit cell by adding nX. The Neighbor list algorithm
         * never requires multiple additions or subtractions of nX.
         *
         * @param i The index along the a-axis.
         * @param j The index along the b-axis.
         * @param k The index along the c-axis.
         * @return The pointer into the 1D cell array.
         */
        private int image(int i, int j, int k) {
            if (i >= nA) {
                i -= nA;
            } else if (i < 0) {
                i += nA;
            }
            if (j >= nB) {
                j -= nB;
            } else if (j < 0) {
                j += nB;
            }
            if (k >= nC) {
                k -= nC;
            } else if (k < 0) {
                k += nC;
            }
            return i + j * nA + k * nAB;
        }

        private void atomCellPairs(final int pairCellIndex) {
            final int atomCellIndex = asymmetricIndex[atomIndex];
            final int i3 = atomIndex * 3;
            final double xi = xyz[i3 + XX];
            final double yi = xyz[i3 + YY];
            final double zi = xyz[i3 + ZZ];
            final int[] pairCellAtoms = cellList[iSymm];
            int start = cellStart[iSymm][pairCellIndex];
            final int pairStop = start + cellCount[iSymm][pairCellIndex];
            final double[] pair = coordinates[iSymm];

            // Check if this pair search is over atoms in the asymmetric unit.
            if (iSymm == 0) {

                // Interactions between atoms in the asymmetric unit may be masked.
                if (maskingRules != null) {
                    maskingRules.applyMask(atomIndex, vdw14, mask);
                }

                // If the self-volume is being searched for pairs, we must avoid double counting.
                if (atomCellIndex == pairCellIndex) {
                    /*
                      The cellOffset is the index of the current atom in the
                      cell that is being searched for neighbors.
                     */
                    start += cellOffset[0][atomIndex] + 1;
                }
            }

            // Loop over atoms in the "pair" cell.
            for (int j = start; j < pairStop; j++) {
                final int aj = pairCellAtoms[j];
                if (use != null && !use[aj]) {
                    continue;
                }
                boolean jactive = atoms[aj].isActive();
                if (!iactive && !jactive && !inactiveInteractions) {
                    continue;
                }
                if (!intermolecular
                        && (molecules[atomIndex] != molecules[aj])) {
                    continue;
                }
                if (mask[aj] > 0 && (iSymm == 0 || aj >= atomIndex)) {
                    int aj3 = aj * 3;
                    final double xj = pair[aj3 + XX];
                    final double yj = pair[aj3 + YY];
                    final double zj = pair[aj3 + ZZ];
                    final double xr = xi - xj;
                    final double yr = yi - yj;
                    final double zr = zi - zj;
                    final double d2 = crystal.image(xr, yr, zr);
                    if (d2 <= cutoffPlusBuffer2) {

                        // Warn about close overlaps.
                        if (d2 < crystal.specialPositionCutoff2 && logger.isLoggable(Level.FINE)) {
                            logger.fine(format(" Close overlap (%6.3f) between atoms (iSymm = %d):\n %s\n %s\n",
                                    sqrt(d2), iSymm, atoms[atomIndex].toString(), atoms[aj].toString()));
                        }

                        // Add the pair to the list, reallocating the array size if necessary.
                        try {
                            pairs[n++] = aj;
                        } catch (Exception e) {
                            n = pairs.length;
                            pairs = copyOf(pairs, n + 100);
                            pairs[n++] = aj;
                        }
                    }
                }
            }

            // Interactions between atoms in the asymmetric unit may be masked.
            if (iSymm == 0 && maskingRules != null) {
                maskingRules.removeMask(atomIndex, vdw14, mask);
            }
        }
    }

    /**
     * <p>destroy.</p>
     *
     * @throws java.lang.Exception if any.
     */
    public void destroy() throws Exception {
        parallelTeam.shutdown();
    }

    /**
     * Moves an array of doubles to be within 0.0 and 1.0 by addition or
     * subtraction of a multiple of 1.0. Typical use is moving an atom placed
     * outside crystal boundaries from the symmetry mate back into the crystal.
     *
     * @param valuesToMove Doubles to be moved between 0 and 1.
     */
    public static void moveValuesBetweenZeroAndOne(double[] valuesToMove) {
        for (int i = 0; i < valuesToMove.length; i++) {
            valuesToMove[i] = moveBetweenZeroAndOne(valuesToMove[i]);
        }
    }

    /**
     * Moves a double to be within 0.0 and 1.0 by addition or subtraction of a
     * multiple of 1.0. Typical use is moving an atom place outside crystal
     * boundaries from the symmetry mate back into the crystal.
     *
     * @param value Double to be moved between 0 and 1.
     * @return Shifted double.
     */
    private static double moveBetweenZeroAndOne(double value) {
        if (value < 0.0) {
            int belowZero = (int) (value / 1.0);
            belowZero = 1 + (-1 * belowZero);
            value = value + belowZero;
        } else {
            value = value % 1.0;
        }
        return value;
    }

    private final static int XX = 0;
    private final static int YY = 1;
    private final static int ZZ = 2;
}
