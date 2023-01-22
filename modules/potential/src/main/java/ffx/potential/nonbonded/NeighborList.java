// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.nonbonded;

import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.BarrierAction;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedBoolean;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.util.Range;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The NeighborList class builds Verlet lists in parallel via a spatial decomposition. <br>
 *
 * <ol>
 *   <li>The unit cell is partitioned into <code>nA * nB * nC</code> smaller axis-aligned cells,
 *       where {nA, nB, nC} are chosen as large as possible subject to the criteria that the length
 *       of each side of a sub-volume (rCellA, rCellB, rCellC) multiplied by (nEdgeA, nEdgeB,
 *       nEdgeC), respectively, must be greater than the cutoff distance <code>Rcut</code> plus a
 *       buffer distance <code>delta</code>: <br>
 *       <code>rCellA * nEdgeA .GE. (Rcut + delta)</code> <br>
 *       <code>rCellB * nEdgeB .GE. (Rcut + delta)</code> <br>
 *       <code>rCellC * nEdgeC .GE. (Rcut + delta)</code> <br>
 *       All neighbors of an atom are in a block of (2*nEdgeA+1)(2*nEdgeB+1)(2*nEdgeC+1)
 *       neighborCells.
 *   <li>Interactions between an atom and neighbors in the asymmetric unit require only half the
 *       neighboring cells to be searched to avoid double counting. However, enumeration of
 *       interactions between an atom in the asymmetric unit and its neighbors in a symmetry mate
 *       require all cells to be searched.
 *   <li>Verlet lists from the search are stored, which reduces the number of neighbors whose
 *       distances must be calculated by a factor of approximately: <br>
 *       <code>(4/3*Pi*Rcut^3)/(neighborCells*Vcell)</code> About 1/3 as many interactions are
 *       contained in the Verlet lists as in the neighboring cells.
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class NeighborList extends ParallelRegion {

  /** The logger. */
  private static final Logger logger = Logger.getLogger(NeighborList.class.getName());

  private static final int XX = 0;
  private static final int YY = 1;
  private static final int ZZ = 2;
  /** The masking rules to apply when building the neighbor list. */
  private final MaskingInterface maskingRules;
  /** The cutoff beyond which the pairwise energy is zero. */
  private final double cutoff;
  /**
   * A buffer, which is added to the cutoff distance, such that the Verlet lists do not need to be
   * calculated for all coordinate changes.
   */
  private final double buffer;
  /** The maximum squared displacement allowed before list rebuild. */
  private final double motion2;
  /** The sum of the cutoff + buffer. */
  private final double cutoffPlusBuffer;
  /** Total^2 for distance comparisons without taking a sqrt. */
  private final double cutoffPlusBuffer2;
  /** The ParallelTeam coordinates use of threads and their schedules. */
  private final ParallelTeam parallelTeam;
  /** Number of threads used by the parallelTeam. */
  private final int threadCount;
  /**
   * If this flag is set, then the Verlet lists will be rebuilt even in the absence of motion.
   */
  private boolean forceRebuild = false;
  /**
   * If this flag is set then the Verlet lists will be printed to the logger.
   */
  private boolean print = false;
  /**
   * Shared Boolean used to detect list rebuild. If any thread detects motion in the MotionLoop, the
   * value is set to true.
   */
  private final SharedBoolean sharedMotion;
  /**
   * Detect motion in the asymmetric unit.
   */
  private final MotionLoop[] motionLoops;
  /**
   * Perform some initialization tasks prior to rebuilding the list.
   */
  private final ListInitBarrierAction listInitBarrierAction;
  /**
   * Assign atoms to cells in parallel. This loop is executed once per list rebuild.
   */
  private final AssignAtomsToCellsLoop[] assignAtomsToCellsLoops;
  /** A Verlet list loop for each thread. */
  private final NeighborListLoop[] verletListLoop;
  /** Total number of interactions. */
  private final SharedInteger sharedCount;
  /** Optimal pairwise ranges. */
  private final Range[] ranges;
  /**
   * The crystal object defines the unit cell dimensions and space group. If replicates of the unit
   * cell are being used, this should be the overall replicates crystal and not the small unit cell.
   */
  private Crystal crystal;
  /** The number of asymmetric units in the unit cell. */
  private int nSymm;
  /** The number of atoms in the asymmetric unit. */
  private int nAtoms;
  /** Reduced coordinates for each symmetry copy. [nSymm][3*nAtoms] */
  private double[][] coordinates;
  /** The reduced coordinates of the asymmetric unit when the list was last rebuilt. */
  private double[] previous;
  /** The Verlet lists. [nSymm][nAtoms][nNeighbors] */
  private int[][][] lists;
  /** Number of interactions per atom. */
  private int[] listCount;
  /** Pairwise ranges for load balancing. */
  private PairwiseSchedule pairwiseSchedule;
  /** Number of interactions between atoms in the asymmetric unit. */
  private int asymmetricUnitCount;
  /** Number of interactions between atoms in the asymmetric unit and atoms in symmetry mates. */
  private int symmetryMateCount;
  /** Number of atoms in the asymmetric unit with interactions lists greater than 0. */
  private int atomsWithIteractions;
  /**
   * The number of subcells that must be searched along the an axis to find all neighbors within the
   * cutoff + buffer distance.
   *
   * <p>If the nEdge == 1 for {X=A,B,C} then all neighbors will be found in 3x3x3 = 27 cells.
   * If each nEdge == 2, then all neighbors will be found in 5x5x5 = 125 cells (in this case the
   * cells are smaller).
   */
  private int nEdge;
  /** The number of divisions along the A-axis. */
  private int nA;
  /** The number of divisions along the B-axis. */
  private int nB;
  /**
   * The number of divisions along the C-Axis.
   */
  private int nC;
  /** The cell indices of each atom along a A-axis. */
  private int[] cellA;
  /** The cell indices of each atom along a B-axis. */
  private int[] cellB;
  /** The cell indices of each atom along a C-axis. */
  private int[] cellC;
  /**
   * The fractional sub-volumes of the unit cell.
   */
  private Cell[][][] cells;
  /** Atoms being used list. */
  private boolean[] use;
  /** List of atoms. */
  private Atom[] atoms;
  /**
   * Time building the neighbor list.
   */
  private long time;
  /** Include intermolecular interactions. */
  private boolean intermolecular = true;
  /** Molecule number for each atom. */
  private int[] molecules = null;
  /**
   * If true, interactions between two inactive atoms are included in the Neighborlist. Set to true
   * to match OpenMM behavior (i.e. interactions between two atoms with zero mass are included). Set
   * to false for more efficient "pure Java" optimizations.
   */
  private final boolean inactiveInteractions = true;
  /** Disable updates to the NeighborList; use with caution. */
  private boolean disableUpdates = false;

  /**
   * Constructor for the NeighborList class.
   *
   * @param maskingRules This parameter may be null.
   * @param crystal Definition of the unit cell and space group.
   * @param atoms The atoms to generate Verlet lists for.
   * @param cutoff The cutoff distance.
   * @param buffer The buffer distance.
   * @param parallelTeam Specifies the parallel environment.
   * @since 1.0
   */
  public NeighborList(MaskingInterface maskingRules, Crystal crystal,
      Atom[] atoms, double cutoff, double buffer, ParallelTeam parallelTeam) {
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

    sharedMotion = new SharedBoolean();
    motionLoops = new MotionLoop[threadCount];
    listInitBarrierAction = new ListInitBarrierAction();
    verletListLoop = new NeighborListLoop[threadCount];
    assignAtomsToCellsLoops = new AssignAtomsToCellsLoop[threadCount];
    for (int i = 0; i < threadCount; i++) {
      motionLoops[i] = new MotionLoop();
      verletListLoop[i] = new NeighborListLoop();
      assignAtomsToCellsLoops[i] = new AssignAtomsToCellsLoop();
    }

    // Initialize the neighbor list builder subcells.
    boolean print = logger.isLoggable(Level.FINE);
    initNeighborList(print);
  }

  /**
   * This method can be called as necessary to build/rebuild the neighbor lists.
   *
   * @param coordinates The coordinates of each atom [nSymm][nAtoms*3].
   * @param lists The neighbor lists [nSymm][nAtoms][nPairs].
   * @param forceRebuild If true, the list is rebuilt even if no atom has moved half the buffer
   *     size.
   * @param use an array of boolean.
   * @param print a boolean.
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
    this.forceRebuild = forceRebuild;
    this.print = print;

    try {
      parallelTeam.execute(this);
    } catch (Exception e) {
      String message = "Fatal exception building neighbor list.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * destroy.
   *
   * @throws java.lang.Exception if any.
   */
  public void destroy() throws Exception {
    parallelTeam.shutdown();
  }

  /**
   * {@inheritDoc}
   *
   * <p>This is method should not be called; it is invoked by Parallel Java.
   *
   * @since 1.0
   */
  @Override
  public void finish() {
    if (disableUpdates) {
      return;
    }
    if (!forceRebuild && !sharedMotion.get()) {
      return;
    }

    // Collect interactions.
    atomsWithIteractions = 0;
    asymmetricUnitCount = 0;
    symmetryMateCount = 0;
    int[][] list = lists[0];
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

    // Print the interactions counts.
    if (print) {
      print();
    }

    // Update the pairwise schedule.
    pairwiseSchedule.updateRanges(sharedCount.get(), atomsWithIteractions, listCount);

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format("   Parallel Neighbor List:           %10.3f sec",
          (System.nanoTime() - time) * 1e-9));
    }
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
   * If disableUpdates true, disable updating the neighbor list upon motion. Use with caution; best
   * recommendation is to only use if all atoms have a coordinate restraint.
   *
   * @param disableUpdate Disable updating the neighbor list
   */
  public void setDisableUpdates(boolean disableUpdate) {
    this.disableUpdates = disableUpdate;
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
   * Getter for the field <code>pairwiseSchedule</code>.
   *
   * @return a {@link ffx.potential.nonbonded.PairwiseSchedule} object.
   */
  public PairwiseSchedule getPairwiseSchedule() {
    return pairwiseSchedule;
  }

  /**
   * {@inheritDoc}
   *
   * <p>This is method should not be called; it is invoked by Parallel Java.
   *
   * @since 1.0
   */
  @Override
  public void run() {
    if (disableUpdates) {
      return;
    }

    int threadIndex = getThreadIndex();
    try {
      if (!forceRebuild) {
        // Check for motion.
        execute(0, nAtoms - 1, motionLoops[threadIndex]);
      }
      if (!forceRebuild && !sharedMotion.get()) {
        return;
      }
      barrier(listInitBarrierAction);
      execute(0, nAtoms - 1, assignAtomsToCellsLoops[threadIndex]);
      execute(0, nAtoms - 1, verletListLoop[threadIndex]);
    } catch (Exception e) {
      String message =
          "Fatal exception building neighbor list in thread: " + getThreadIndex() + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * The NeighborList will be re-configured, if necessary, for the supplied atom list.
   *
   * @param atoms A new list of atoms to operate on.
   */
  public void setAtoms(Atom[] atoms) {
    this.atoms = atoms;
    this.nAtoms = atoms.length;
    initNeighborList(false);
  }

  /**
   * The NeighborList will be re-configured, if necessary, for the supplied Crystal. Changes to both
   * unit cell parameters and number of symmetry operators are also acceptable.
   *
   * @param crystal A crystal defining boundary conditions and symmetry.
   */
  public void setCrystal(Crystal crystal) {
    this.crystal = crystal;
    initNeighborList(false);
  }

  /**
   * Setter for the field <code>intermolecular</code>.
   *
   * @param intermolecular a boolean.
   * @param molecules an array of {@link int} objects.
   */
  public void setIntermolecular(boolean intermolecular, int[] molecules) {
    this.intermolecular = intermolecular;
    this.molecules = molecules;
  }

  /**
   * {@inheritDoc}
   *
   * <p>This is method should not be called; it is invoked by Parallel Java.
   *
   * @since 1.0
   */
  @Override
  public void start() {
    if (disableUpdates) {
      return;
    }
    time = System.nanoTime();
    sharedMotion.set(false);
  }

  private void initNeighborList(boolean print) {

    // Allocate memory for fractional coordinates and subcell pointers for each atom.
    if (cellA == null || cellA.length < nAtoms) {
      cellA = new int[nAtoms];
      cellB = new int[nAtoms];
      cellC = new int[nAtoms];
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
    final double sphere = min(min(crystal.interfacialRadiusA, crystal.interfacialRadiusB),
        crystal.interfacialRadiusC);

    /*
     Assert that the boundary conditions defined by the crystal allow use
     of the minimum image condition.
    */
    if (!crystal.aperiodic()) {
      assert (sphere >= cutoffPlusBuffer);
    }

    // nEdge >= 1.
    // For nEdge == 1, then 3x3x3 subcells are searched for neighbors (mean of ~440 water atoms per cell).
    // For nEdge == 2, then 5x5x5 subcells are searched for neighbors (mean of ~40 water atoms per cell).
    // For nEdge == 3, then 7x7x7 subcells are searched for neighbors (mean of ~12 water atoms per cell).
    nEdge = 2;

    /*
     The targetInterfacialRadius is the smallest interfacial radius for a subcell,
     given nEdge subcells along an axis, that assures all neighbors will be found.
    */
    double targetInterfacialRadius = cutoffPlusBuffer / (double) nEdge;
    nA = (int) floor(2 * crystal.interfacialRadiusA / targetInterfacialRadius);
    nB = (int) floor(2 * crystal.interfacialRadiusB / targetInterfacialRadius);
    nC = (int) floor(2 * crystal.interfacialRadiusC / targetInterfacialRadius);
    if (nA < 2 * nEdge + 1) {
      nA = 1;
    }
    if (nB < 2 * nEdge + 1) {
      nB = 1;
    }
    if (nC < 2 * nEdge + 1) {
      nC = 1;
    }

    cells = new Cell[nA][nB][nC];
    for (int i = 0; i < nA; i++) {
      for (int j = 0; j < nB; j++) {
        for (int k = 0; k < nC; k++) {
          cells[i][j][k] = new Cell(i, j, k);
        }
      }
    }

    if (print) {
      int nCells = nA * nB * nC;
      StringBuilder sb = new StringBuilder("  Neighbor List Builder\n");
      sb.append(format("   Total Cells:          %8d\n", nCells));
      if (nCells > 1) {
        sb.append(format("   Interfacial radius    %8.3f A\n", targetInterfacialRadius));
        int searchA = nA == 1 ? 1 : 2 * nEdge + 1;
        int searchB = nB == 1 ? 1 : 2 * nEdge + 1;
        int searchC = nC == 1 ? 1 : 2 * nEdge + 1;
        sb.append(format("   Domain Decomposition: %8d %4d %4d\n", nA, nB, nC));
        sb.append(format("   Neighbor Search:      %8d x%3d x%3d\n", searchA, searchB, searchC));
        sb.append(format("   Mean Atoms per Cell:  %8d", nAtoms * nSymm / nCells));
      }
      logger.info(sb.toString());
    }
  }

  private void print() {
    StringBuilder sb =
        new StringBuilder(format("   Buffer:                                %5.2f (A)\n", buffer));
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
   * Detect if any atom has moved 1/2 the buffer size.
   *
   * @since 1.0
   */
  private class MotionLoop extends IntegerForLoop {

    @Override
    public void run(int lb, int ub) {
      double[] current = coordinates[0];
      for (int i = lb; i <= ub; i++) {
        int i3 = i * 3;
        int iX = i3 + XX;
        int iY = i3 + YY;
        int iZ = i3 + ZZ;
        double dx = previous[iX] - current[iX];
        double dy = previous[iY] - current[iY];
        double dz = previous[iZ] - current[iZ];
        double dr2 = crystal.image(dx, dy, dz);
        if (dr2 > motion2) {
          sharedMotion.set(true);
          return;
        }
      }
    }
  }

  /**
   * Perform some initialization tasks prior to list rebuild.
   */
  private class ListInitBarrierAction extends BarrierAction {

    @Override
    public void run() {
      sharedCount.set(0);

      // Clear cell contents.
      for (int i = 0; i < nA; i++) {
        for (int j = 0; j < nB; j++) {
          for (int k = 0; k < nC; k++) {
            cells[i][j][k].clear();
          }
        }
      }

      // Allocate memory for neighbor lists.
      for (int iSymm = 0; iSymm < nSymm; iSymm++) {
        if (lists[iSymm] == null || lists[iSymm].length < nAtoms) {
          lists[iSymm] = new int[nAtoms][];
        }
      }
    }
  }

  /**
   * Assign atoms to cells.
   *
   * @since 1.0
   */
  private class AssignAtomsToCellsLoop extends IntegerForLoop {

    private final double[] cart = new double[3];
    private final double[] frac = new double[3];

    @Override
    public void run(int lb, int ub) throws Exception {

      // Save the current coordinates.
      double[] current = coordinates[0];
      for (int i = lb; i <= ub; i++) {
        int i3 = i * 3;
        int iX = i3 + XX;
        int iY = i3 + YY;
        int iZ = i3 + ZZ;
        previous[iX] = current[iX];
        previous[iY] = current[iY];
        previous[iZ] = current[iZ];
      }

      for (int iSymm = 0; iSymm < nSymm; iSymm++) {
        // Convert to fractional coordinates.
        final double[] xyz = coordinates[iSymm];
        // Assign each atom to a cell using fractional coordinates.
        for (int i = lb; i <= ub; i++) {
          int i3 = i * 3;
          cart[0] = xyz[i3 + XX];
          cart[1] = xyz[i3 + YY];
          cart[2] = xyz[i3 + ZZ];
          crystal.toFractionalCoordinates(cart, frac);
          double xu = frac[0];
          double yu = frac[1];
          double zu = frac[2];
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
          cells[a][b][c].add(i, iSymm);
        }
      }
    }
  }

  /**
   * The VerletListLoop class encapsulates thread local variables and methods for building Verlet
   * lists based on a spatial decomposition of the unit cell.
   *
   * @author Michael J. Schnieders
   * @since 1.0
   */
  private class NeighborListLoop extends IntegerForLoop {

    private final IntegerSchedule schedule;
    private int n;
    private int iSymm;
    private int atomIndex;
    private boolean iactive = true;
    private int count;
    private double[] xyz;
    private int[] pairs;
    private Cell cellForCurrentAtom;
    private int[] pairCellAtoms;
    private double[] mask;
    private boolean[] vdw14;
    // Extra padding to avert cache interference.
    private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
    private long pad8, pad9, pada, padb, padc, padd, pade, padf;

    NeighborListLoop() {
      int len = 1000;
      pairs = new int[len];
      schedule = IntegerSchedule.dynamic(10);
      pairCellAtoms = new int[len];
    }

    @Override
    public void finish() {
      sharedCount.addAndGet(count);
    }

    @Override
    public void run(final int lb, final int ub) {
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

            cellForCurrentAtom = cells[a][b][c];

            int a1 = a + 1;
            int aStart = a - nEdge;
            int aStop = a + nEdge;
            int b1 = b + 1;
            int bStart = b - nEdge;
            int bStop = b + nEdge;
            int c1 = c + 1;
            int cStart = c - nEdge;
            int cStop = c + nEdge;

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
              atomCellPairs(cellForCurrentAtom);
              // Search half of the neighboring volumes to avoid double counting.
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

    /**
     * If the index is >= to nX, it is mapped back into the periodic unit cell by subtracting nX. If
     * the index is less than 0, it is mapped into the periodic unit cell by adding nX. The Neighbor
     * list algorithm never requires multiple additions or subtractions of nX.
     *
     * @param i The index along the a-axis.
     * @param j The index along the b-axis.
     * @param k The index along the c-axis.
     * @return The requested Cell.
     */
    private Cell image(int i, int j, int k) {
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
      return cells[i][j][k];
    }

    private void atomCellPairs(final Cell cell) {
      final int i3 = atomIndex * 3;
      final double xi = xyz[i3 + XX];
      final double yi = xyz[i3 + YY];
      final double zi = xyz[i3 + ZZ];

      if (pairCellAtoms.length < cell.getCount()) {
        pairCellAtoms = new int[cell.getCount()];
      }
      // Load the cell's atom indices into an array for the specified SymOp.
      int num = cell.getSymOpAtoms(iSymm, pairCellAtoms);
      final double[] pair = coordinates[iSymm];

      // Check if this pair search is over atoms in the asymmetric unit.
      if (iSymm == 0 && maskingRules != null) {
        // Interactions between atoms in the asymmetric unit may be masked.
        maskingRules.applyMask(atomIndex, vdw14, mask);
      }

      // Loop over atoms in the "pair" cell.
      for (int j = 0; j < num; j++) {
        final int aj = pairCellAtoms[j];
        // If the self-volume is being searched for pairs, avoid double counting.
        if (iSymm == 0 && cell == cellForCurrentAtom) {
          // Only consider atoms with a higher index than the current atom.
          if (aj <= atomIndex) {
            continue;
          }
        }

        if (use != null && !use[aj]) {
          continue;
        }
        boolean jactive = atoms[aj].isActive();
        if (!iactive && !jactive && !inactiveInteractions) {
          continue;
        }
        if (!intermolecular && (molecules[atomIndex] != molecules[aj])) {
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
            if (d2 < crystal.getSpecialPositionCutoff2() && logger.isLoggable(Level.FINE)) {
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
   * Hold the atom index and its symmetry operator.
   */
  private static class AtomIndex {

    public final int iSymm;
    public final int i;

    public AtomIndex(int iSymm, int i) {
      this.iSymm = iSymm;
      this.i = i;
    }
  }

  /**
   * Hold the atoms in each cell.
   */
  private static class Cell {

    /**
     * The cell index along the A-axis.
     */
    final int a;
    /**
     * The cell index along the B-axis.
     */
    final int b;
    /**
     * The cell index along the C-axis.
     */
    final int c;

    /**
     * The list of atoms in the cell, together with their symmetry operator.
     */
    final List<AtomIndex> list;

    public Cell(int a, int b, int c) {
      this.a = a;
      this.b = b;
      this.c = c;
      list = Collections.synchronizedList(new ArrayList<>());
    }

    /**
     * Add an atom to the cell.
     *
     * @param atomIndex The atom index.
     * @param symOpIndex The symmetry operator index.
     */
    public void add(int atomIndex, int symOpIndex) {
      list.add(new AtomIndex(symOpIndex, atomIndex));
    }

    public AtomIndex get(int index) {
      if (index >= list.size()) {
        return null;
      }
      return list.get(index);
    }

    public int getCount() {
      return list.size();
    }

    /**
     * Return the number of atoms in the cell for a given symmetry operator.
     *
     * @param symOpIndex The symmetry operator index.
     * @param index The list of indexes for the given symmetry operator.
     * @return The number of atoms in the cell for the symmetry operator.
     */
    public int getSymOpAtoms(int symOpIndex, int[] index) {
      int count = 0;
      for (AtomIndex atomIndex : list) {
        if (atomIndex.iSymm == symOpIndex) {
          if (count >= index.length) {
            throw new IndexOutOfBoundsException("Index out of bounds: count "
                + count + " " + index.length + " " + list.size());
          }
          index[count++] = atomIndex.i;
        }
      }
      return count;
    }

    /**
     * Clear the list of atoms in the cell.
     */
    public void clear() {
      list.clear();
    }

  }
}
