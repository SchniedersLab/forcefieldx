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
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import com.github.quickhull3d.QuickHull3D;
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
import ffx.potential.utils.ConvexHullOps;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The NeighborList class builds Verlet lists in parallel via a spatial decomposition. <br>
 *
 * <ol>
 *   <li>The unit cell is partitioned into <code>nA * nB * nC</code> smaller axis-aligned cells,
 *       where {nA, nB, nC} are chosen as large as possible subject to the criteria that the length
 *       of each side of a sub-volume (rCellA, rCellB, rCellC) multiplied by nEdge
 *       must be greater than the cutoff distance <code>Rcut</code> plus a
 *       buffer distance <code>delta</code>: <br>
 *       <code>rCellA * nEdge .GE. (Rcut + delta)</code> <br>
 *       <code>rCellB * nEdge .GE. (Rcut + delta)</code> <br>
 *       <code>rCellC * nEdge .GE. (Rcut + delta)</code> <br>
 *       All neighbors of an atom are in a block of (2*nEdge+1)^3 neighborCells.
 *   <li>Interactions between an atom and neighbors in the asymmetric unit require only half the
 *       neighboring cells to be searched to avoid double counting. However, enumeration of
 *       interactions between an atom in the asymmetric unit and its neighbors in a symmetry mate
 *       require all cells to be searched.
 *   <li>Verlet lists from the search are stored, which reduces the number of neighbors whose
 *       distances must be calculated by a factor of approximately: <br>
 *       <code>(4/3*Pi*Rcut^3)/(neighborCells*Vcell)</code>
 *       About 1/3 as many interactions are contained in the Verlet lists compared to
 *       the total amount in the neighboring cells.
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class NeighborList extends ParallelRegion {

  /** The logger. */
  private static final Logger logger = Logger.getLogger(NeighborList.class.getName());

  private final DomainDecomposition domainDecomposition;
  private DomainDecompositionOctree domainDecompositionOctree;
  private boolean useOctree = true;
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
  /**
   * Time to check for motion.
   */
  private long motionTime;
  /**
   * Time for some initializations.
   */
  private long initTime;
  /**
   * Time to assign atoms to cells.
   */
  private long assignAtomsToCellsTime;
  /**
   * Time to build the Verlet lists.
   */
  private long verletListTime;
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
  /** The Cartesian coordinates of the asymmetric unit when the list was last rebuilt. */
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
  /**
   * If true, interactions between two inactive atoms are included in the Neighborlist. Set to true
   * to match OpenMM behavior (i.e. interactions between two atoms with zero mass are included). Set
   * to false for more efficient "pure Java" optimizations.
   */
  private final boolean includeInactivePairs = true;
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
  public NeighborList(MaskingInterface maskingRules, Crystal crystal, Atom[] atoms, double cutoff,
      double buffer, ParallelTeam parallelTeam) {
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
    domainDecomposition = new DomainDecomposition(nAtoms, crystal, cutoffPlusBuffer);

    if(useOctree) {
      QuickHull3D quickHull3D = ConvexHullOps.constructHull(atoms);
      double maxDist = ConvexHullOps.maxDist(quickHull3D);
      domainDecompositionOctree = new DomainDecompositionOctree(nAtoms, crystal, maxDist, 8);
    }
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
  public void buildList(final double[][] coordinates, final int[][][] lists, boolean[] use,
      boolean forceRebuild, boolean print) {
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

  public void buildOctree(int nAtoms, Crystal crystal, double cutoffPlusBuffer, int dimension){

//    domainDecompositionOctree.initDomainDecomposition(nAtoms, crystal);
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
//    logger.info("NL dDO log in finish method");
//    domainDecompositionOctree.log();

    // Update the pair-wise schedule.
    long scheduleTime = -System.nanoTime();
    pairwiseSchedule.updateRanges(sharedCount.get(), atomsWithIteractions, listCount);
    scheduleTime += System.nanoTime();

    if (logger.isLoggable(Level.FINE)) {
      time = System.nanoTime() - time;
      StringBuilder sb = new StringBuilder();
      sb.append(format("   Motion Check:           %6.4f sec\n", motionTime * 1e-9));
      sb.append(format("   List Initialization:    %6.4f sec\n", initTime * 1e-9));
      sb.append(format("   Assign Atoms to Cells:  %6.4f sec\n", assignAtomsToCellsTime * 1e-9));
      sb.append(format("   Create Vertlet Lists:   %6.4f sec\n", verletListTime * 1e-9));
      sb.append(format("   Parallel Schedule:      %6.4f sec\n", scheduleTime * 1e-9));
      sb.append(format("   Neighbor List Total:    %6.4f sec\n", time * 1e-9));
      logger.fine(sb.toString());
    }
  }

  public Cell[][][] getCells() {
    return domainDecompositionOctree.getCells();
  }

  /**
   * Returns the cutoff distance used internally by NeighborList.
   *
   * @return Cutoff distance in Angstroms.
   */
  public double getCutoff() {
    return cutoff;
  }

  public int[] getNumDivisions() {
    return domainDecompositionOctree.getNumDivisions();
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
        if (threadIndex == 0) {
          motionTime = -System.nanoTime();
        }
        execute(0, nAtoms - 1, motionLoops[threadIndex]);
        if (threadIndex == 0) {
          motionTime += System.nanoTime();
        }
      }
      if (forceRebuild || sharedMotion.get()) {
        // Complete some initializations.
        barrier(listInitBarrierAction);
        // Assign atoms to cells.
        if (threadIndex == 0) {
          assignAtomsToCellsTime = -System.nanoTime();
        }
        execute(0, nAtoms - 1, assignAtomsToCellsLoops[threadIndex]);
        if (threadIndex == 0) {
          assignAtomsToCellsTime += System.nanoTime();
        }
        // Build the Verlet list.
        if (threadIndex == 0) {
          verletListTime = -System.nanoTime();
        }
        execute(0, nAtoms - 1, verletListLoop[threadIndex]);
        if (threadIndex == 0) {
          verletListTime += System.nanoTime();
        }
      }
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
   * @param intermolecular If true, the NeighborList will include intermolecular interactions.
   */
  public void setIntermolecular(boolean intermolecular) {
    this.intermolecular = intermolecular;
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
    motionTime = 0;
    initTime = 0;
    assignAtomsToCellsTime = 0;
    verletListTime = 0;
    sharedMotion.set(false);
  }

  private void initNeighborList(boolean print) {

    // Set the number of symmetry operators.
    int newNSymm = crystal.spaceGroup.symOps.size();
    if (nSymm != newNSymm) {
      nSymm = newNSymm;
    }

    // Allocate memory for fractional coordinates and subcell pointers for each atom.
    if (previous == null || previous.length < 3 * nAtoms) {
      previous = new double[3 * nAtoms];
      listCount = new int[nAtoms];
      pairwiseSchedule = new PairwiseSchedule(threadCount, nAtoms, ranges);
    } else {
      pairwiseSchedule.setAtoms(nAtoms);
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

    domainDecomposition.initDomainDecomposition(nAtoms, crystal);
    if (print) {
      domainDecomposition.log();
    }

    // Initialize DomainDecompositionOctree
    if(useOctree) {
    domainDecompositionOctree.initDomainDecomposition(nAtoms, crystal);
//    logger.info("NL dDO log in initNeighborList method");
    domainDecompositionOctree.log();
    }
  }

  private void print() {
    StringBuilder sb = new StringBuilder(
        format("   Buffer:                                %5.2f (A)\n", buffer));
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

  private static final int XX = 0;
  private static final int YY = 1;
  private static final int ZZ = 2;

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
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Motion detected for atom %d (%8.6f A).", i, sqrt(dr2)));
          }
          sharedMotion.set(true);
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
      initTime = -System.nanoTime();
      // Clear the list count.
      sharedCount.set(0);

      domainDecomposition.clear();

      // Allocate memory for neighbor lists.
      for (int iSymm = 0; iSymm < nSymm; iSymm++) {
        if (lists[iSymm] == null || lists[iSymm].length < nAtoms) {
          lists[iSymm] = new int[nAtoms][];
        }
      }
      initTime += System.nanoTime();
    }
  }

  /**
   * Assign atoms to cells.
   *
   * @since 1.0
   */
  private class AssignAtomsToCellsLoop extends IntegerForLoop {

    /**
     * The cartesian coordinates of the atom in the asymmetric unit.
     */
    private final double[] auCart = new double[3];
    /**
     * The cartesian coordinates of the atom after symmetry application.
     */
    private final double[] cart = new double[3];
    /**
     * The fractional coordinates of the atom after symmetry application.
     */
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

      final double[] auXYZ = coordinates[0];
      final double spCutoff2 = crystal.getSpecialPositionCutoff2();
      for (int iSymm = 0; iSymm < nSymm; iSymm++) {
        final double[] xyz = coordinates[iSymm];
        // Assign each atom to a cell using fractional coordinates.
        for (int i = lb; i <= ub; i++) {
          int i3 = i * 3;
          cart[0] = xyz[i3 + XX];
          cart[1] = xyz[i3 + YY];
          cart[2] = xyz[i3 + ZZ];
          if (iSymm != 0) {
            auCart[0] = auXYZ[i3 + XX];
            auCart[1] = auXYZ[i3 + YY];
            auCart[2] = auXYZ[i3 + ZZ];
            double dx = auCart[0] - cart[0];
            double dy = auCart[1] - cart[1];
            double dz = auCart[2] - cart[2];
            double dr2 = crystal.image(dx, dy, dz);
            if (dr2 < spCutoff2) {
              int molecule = atoms[i].getMoleculeNumber();
              if (atoms[i].isActive()) {
                logger.info(format(
                    "   Active Atom %s at Special Position in Molecule %d with SymOp %d (%8.6f A).",
                    atoms[i], molecule, iSymm, sqrt(dr2)));
              } else if (logger.isLoggable(Level.FINE)) {
                logger.fine(
                    format("   Atom %s at Special Position in Molecule %d with SymOp %d (%8.6f A).",
                        atoms[i], molecule, iSymm, sqrt(dr2)));
              }
              // Exclude molecule-molecule interactions between the asymmetric unit and the iSymm symmetry mate.
              domainDecomposition.addSpecialPositionExclusion(molecule, iSymm);
            }
          }
          // Convert to fractional coordinates.
          crystal.toFractionalCoordinates(cart, frac);
          domainDecomposition.addAtomToCell(i, iSymm, frac);
          domainDecompositionOctree.addAtomToCell(i, iSymm, frac);
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
      int nEdge = domainDecomposition.nEdge;
      int nA = domainDecomposition.nA;
      int nB = domainDecomposition.nB;
      int nC = domainDecomposition.nC;
      for (iSymm = 0; iSymm < nSymm; iSymm++) {
        int[][] list = lists[iSymm];
        // Loop over all atoms.
        for (atomIndex = lb; atomIndex <= ub; atomIndex++) {
          n = 0;

          if (iSymm == 0) {
            listCount[atomIndex] = 0;
            int molecule = atoms[atomIndex].getMoleculeNumber();
            List<Integer> symOpList = domainDecomposition.getSpecialPositionSymOps(molecule);
            atoms[atomIndex].setSpecialPositionSymOps(symOpList);
          }

          if (use == null || use[atomIndex]) {

            iactive = atoms[atomIndex].isActive();

            cellForCurrentAtom = domainDecomposition.getCellForAtom(atomIndex);
            int a = cellForCurrentAtom.a;
            int b = cellForCurrentAtom.b;
            int c = cellForCurrentAtom.c;
            int a1 = a + 1;
            int b1 = b + 1;
            int c1 = c + 1;

            /*
             If the number of divisions is 1 in any direction then
             set the loop limits to the current cell value.
             Otherwise search from a - nEdge to a + nEdge.
            */
            int aStart = nA == 1 ? a : a - nEdge;
            int aStop = nA == 1 ? a : a + nEdge;
            int bStart = nB == 1 ? b : b - nEdge;
            int bStop = nB == 1 ? b : b + nEdge;
            int cStart = nC == 1 ? c : c - nEdge;
            int cStop = nC == 1 ? c : c + nEdge;

            if (iSymm == 0) {
              // Interactions within the "self-volume".
              atomCellPairs(cellForCurrentAtom);
              // Search half of the neighboring volumes to avoid double counting.
              // (a, b+1..b+nE, c)
              for (int bi = b1; bi <= bStop; bi++) {
                atomCellPairs(domainDecomposition.image(a, bi, c));
              }
              // (a, b-nE..b+nE, c+1..c+nE)
              for (int bi = bStart; bi <= bStop; bi++) {
                for (int ci = c1; ci <= cStop; ci++) {
                  atomCellPairs(domainDecomposition.image(a, bi, ci));
                }
              }
              // (a+1..a+nE, b-nE..b+nE, c-nE..c+nE)
              for (int bi = bStart; bi <= bStop; bi++) {
                for (int ci = cStart; ci <= cStop; ci++) {
                  for (int ai = a1; ai <= aStop; ai++) {
                    atomCellPairs(domainDecomposition.image(ai, bi, ci));
                  }
                }
              }
            } else {
              // Interactions with all adjacent symmetry mate cells.
              for (int ai = aStart; ai <= aStop; ai++) {
                for (int bi = bStart; bi <= bStop; bi++) {
                  for (int ci = cStart; ci <= cStop; ci++) {
                    atomCellPairs(domainDecomposition.image(ai, bi, ci));
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

    private void atomCellPairs(final Cell cell) {
      if (pairCellAtoms.length < cell.getCount()) {
        pairCellAtoms = new int[cell.getCount()];
      }

      // Load the cell's atom indices into an array for the specified SymOp.
      int num = cell.getSymOpAtoms(iSymm, pairCellAtoms);
      if (num == 0) {
        return;
      }

      // Check if this pair search is over atoms in the asymmetric unit.
      if (iSymm == 0 && maskingRules != null) {
        // Interactions between atoms in the asymmetric unit may be masked.
        maskingRules.applyMask(atomIndex, vdw14, mask);
      }

      final int moleculeIndex = atoms[atomIndex].getMoleculeNumber();

      boolean logClosePairs = logger.isLoggable(Level.FINE);

      final int i3 = atomIndex * 3;
      final double xi = xyz[i3 + XX];
      final double yi = xyz[i3 + YY];
      final double zi = xyz[i3 + ZZ];
      final double[] pair = coordinates[iSymm];

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

        if (!includeInactivePairs && !iactive) {
          if (!atoms[aj].isActive()) {
            continue;
          }
        }

        int moleculeIndexJ = atoms[aj].getMoleculeNumber();
        // Check if this pair search is only intra-molecular.
        if (!intermolecular && (moleculeIndex != moleculeIndexJ)) {
          continue;
        }

        // Check for a special position exclusion between two atoms in the same molecule.
        if (iSymm != 0 && (moleculeIndex == moleculeIndexJ)) {
          if (domainDecomposition.isSpecialPositionExclusion(moleculeIndex, iSymm)) {
            if (logger.isLoggable(Level.FINEST)) {
              logger.finest(
                  format("   Excluding Interaction for Atoms %d and %d in Molecule %d for SymOp %d.",
                      atomIndex, aj, moleculeIndex, iSymm));
            }
            continue;
          }
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
            if (logClosePairs && d2 < crystal.getSpecialPositionCutoff2()) {
              // Warn about close interactions.
              logger.fine(
                  format(" Close interaction (%6.3f) between atoms (iSymm = %d):\n %s\n %s\n",
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
   * Break down the simulation domain into a 3D grid of cells.
   */
  private static class DomainDecomposition {

    /**
     * The crystal that defines the simulation domain.
     */
    private Crystal crystal;
    /**
     * The targetInterfacialRadius is the smallest interfacial radius for a subcell, given a search
     * of nEdge subcells along an axis, that assures all neighbors will be found.
     */
    private final double targetInterfacialRadius;
    /**
     * The number of subcells that must be searched along the an axis to find all neighbors within
     * the cutoff + buffer distance.
     *
     * <p>If the nEdge == 1 for {X=A,B,C} then all neighbors will be found in 3x3x3 = 27 cells.
     * If each nEdge == 2, then all neighbors will be found in 5x5x5 = 125 cells (in this case the
     * cells are smaller).
     */
    private final int nEdge;
    /**
     * The number of subcells that must be searched along an axis to find all neighbors within the
     * cutoff + buffer distance. nSearch = 2 * nEdge + 1
     */
    private final int nSearch;
    /** The number of divisions along the A-axis. */
    private int nA;
    /** The number of divisions along the B-axis. */
    private int nB;
    /**
     * The number of divisions along the C-Axis.
     */
    private int nC;
    /**
     * The number of atoms.
     */
    private int nAtoms;
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
    /**
     * Special position exclusions. Map<Molecule, List<SymOpIndex>>
     */
    private Map<Integer, List<Integer>> specialPositionExclusionMap;


    /**
     * DomainDecomposition Constructor.
     *
     * @param nAtoms Number of atoms.
     * @param crystal The crystal.
     * @param cutoffPlusBuffer The cutoff plus buffer distance.
     */
    public DomainDecomposition(int nAtoms, Crystal crystal, double cutoffPlusBuffer) {
      this.nEdge = 2;
      this.nSearch = 2 * nEdge + 1;
      targetInterfacialRadius = cutoffPlusBuffer / (double) nEdge;
      initDomainDecomposition(nAtoms, crystal);
    }

    /**
     * Initialize the domain decomposition.
     */
    public void initDomainDecomposition(int nAtoms, Crystal crystal) {
      this.nAtoms = nAtoms;
      this.crystal = crystal;

      // Allocate memory for fractional coordinates and subcell pointers for each atom.
      if (cellA == null || cellA.length < nAtoms) {
        cellA = new int[nAtoms];
        cellB = new int[nAtoms];
        cellC = new int[nAtoms];
      }

      // Determine the number of subcells along each axis.
      int maxA = (int) floor(2 * crystal.interfacialRadiusA / targetInterfacialRadius);
      int maxB = (int) floor(2 * crystal.interfacialRadiusB / targetInterfacialRadius);
      int maxC = (int) floor(2 * crystal.interfacialRadiusC / targetInterfacialRadius);

      /*
        To satisfy the cutoff + buffer distance, the number of subcells (maxA, maxB, maxC)
        is guaranteed to be at least nSearch - 1 for periodic crystals.
        */
      if (!crystal.aperiodic()) {
        assert (maxA >= nSearch - 1);
        assert (maxB >= nSearch - 1);
        assert (maxC >= nSearch - 1);
      }

      /*
        If the number of subcells less than nSearch, then turn off the overhead of the
        subcell search by setting the number of cells along that axis to 1.
       */
      nA = maxA < nSearch ? 1 : maxA;
      nB = maxB < nSearch ? 1 : maxB;
      nC = maxC < nSearch ? 1 : maxC;

      // Allocate memory for the subcells.
      if (cells == null || cells.length != nA || cells[0].length != nB || cells[0][0].length != nC) {
        cells = new Cell[nA][nB][nC];
        for (int i = 0; i < nA; i++) {
          for (int j = 0; j < nB; j++) {
            for (int k = 0; k < nC; k++) {
              cells[i][j][k] = new Cell(i, j, k, targetInterfacialRadius);
            }
          }
        }
      } else {
        clear();
      }

    }

    /**
     * Clear the sub-cells of the atom list.
     */
    public void clear() {
      // Clear cell contents.
      for (int i = 0; i < nA; i++) {
        for (int j = 0; j < nB; j++) {
          for (int k = 0; k < nC; k++) {
            cells[i][j][k].clear();
          }
        }
      }

      // Clear special position exclusions.
      // if (specialPositionExclusionMap != null) {
      // Clear the lists, but don't remove the keys.
      //   for (List<Integer> list : specialPositionExclusionMap.values()) {
      //     list.clear();
      //   }
      //   specialPositionExclusionMap.clear();
      // }
    }

    /**
     * Add a special position exclusion.
     *
     * @param molecule The molecule.
     * @param symOpIndex The symmetry operation index.
     */
    public void addSpecialPositionExclusion(int molecule, int symOpIndex) {
      // Initialize the map.
      if (specialPositionExclusionMap == null) {
        specialPositionExclusionMap = new HashMap<>();
      }
      // Initialize the list for the molecule.
      List<Integer> list = specialPositionExclusionMap.get(molecule);
      if (list == null) {
        list = new ArrayList<>();
        specialPositionExclusionMap.put(molecule, list);
      }
      // Add the symOpIndex to the list.
      if (!list.contains(symOpIndex)) {
        list.add(symOpIndex);
      }
    }

    /**
     * Check if a special position exclusion exists.
     *
     * @param molecule The molecule.
     * @param symOpIndex The symmetry operation index.
     * @return True if the special position exclusion exists.
     */
    public boolean isSpecialPositionExclusion(int molecule, int symOpIndex) {
      if (specialPositionExclusionMap == null) {
        return false;
      }
      List<Integer> list = specialPositionExclusionMap.get(molecule);
      if (list == null) {
        return false;
      }
      return list.contains(symOpIndex);
    }

    /**
     * Get the special position SymOps for a molecule.
     *
     * @param molecule The molecule.
     * @return The special position SymOps.
     */
    public List<Integer> getSpecialPositionSymOps(int molecule) {
      if (specialPositionExclusionMap == null) {
        return null;
      }
      if (specialPositionExclusionMap.containsKey(molecule)) {
        List<Integer> list = specialPositionExclusionMap.get(molecule);
        if (list.isEmpty()) {
          return null;
        } else {
          return list;
        }
      }
      return null;
    }

    /**
     * Log information about the domain decomposition.
     */
    public void log() {
      int nCells = nA * nB * nC;
      int nSymm = crystal.spaceGroup.getNumberOfSymOps();
      StringBuilder sb = new StringBuilder("  Neighbor List Builder\n");
      sb.append(format("   Total Cells:          %8d\n", nCells));
      if (nCells > 1) {
        sb.append(format("   Interfacial radius    %8.3f A\n", targetInterfacialRadius));
        int searchA = nA == 1 ? 1 : nSearch;
        int searchB = nB == 1 ? 1 : nSearch;
        int searchC = nC == 1 ? 1 : nSearch;
        sb.append(format("   Domain Decomposition: %8d %4d %4d\n", nA, nB, nC));
        sb.append(format("   Neighbor Search:      %8d x%3d x%3d\n", searchA, searchB, searchC));
        sb.append(format("   Mean Atoms per Cell:  %8d", nAtoms * nSymm / nCells));
      }
      logger.info(sb.toString());
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
    public Cell image(int i, int j, int k) {
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

    /**
     * Add an atom to a sub-cell.
     *
     * @param i The index of the atom.
     * @param iSymm The index of the symmetry operator.
     * @param frac The fractional coordinates of the atom.
     */
    public void addAtomToCell(int i, int iSymm, double[] frac) {
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
        // Set the sub-cell indices for an asymmetric unit atom.
        cellA[i] = a;
        cellB[i] = b;
        cellC[i] = c;
      }
      cells[a][b][c].add(i, iSymm);
    }

    /**
     * Get the sub-cell for an asymmetric unit atom.
     *
     * @param i The index of the atom.
     * @return The sub-cell for the atom.
     */
    public Cell getCellForAtom(int i) {
      return cells[cellA[i]][cellB[i]][cellC[i]];
    }
  }

  /**
   * Break down the simulation domain into a 3D grid of cells.
   */
  private static class DomainDecompositionOctree {

    /**
     * The crystal that defines the simulation domain.
     */
    private Crystal crystal;
    /**
     * The targetInterfacialRadius is the smallest interfacial radius for a subcell, given a search
     * of nEdge subcells along an axis, that assures all neighbors will be found.
     */
    private final double targetInterfacialRadius;
    /**
     * The number of subcells that must be searched along an axis to find all neighbors within
     * the cutoff + buffer distance.
     *
     * <p>If the nEdge == 1 for {X=A,B,C} then all neighbors will be found in 3x3x3 = 27 cells.
     * If each nEdge == 2, then all neighbors will be found in 5x5x5 = 125 cells (in this case the
     * cells are smaller).
     */
    private final int nEdge;
    /**
     * The number of subcells that must be searched along an axis to find all neighbors within the
     * cutoff + buffer distance. nSearch = 2 * nEdge + 1
     */
    private final int nSearch;
    /** The number of divisions along the A-axis. */
    private int nA;
    /** The number of divisions along the B-axis. */
    private int nB;
    /**
     * The number of divisions along the C-Axis.
     */
    private int nC;
    /**
     * The number of atoms.
     */
    private int nAtoms;
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
    /**
     * Special position exclusions. Map<Molecule, List<SymOpIndex>>
     */
    private Map<Integer, List<Integer>> specialPositionExclusionMap;


    /**
     * DomainDecomposition Constructor.
     *
     * @param nAtoms Number of atoms.
     * @param crystal The crystal.
     * @param cutoffPlusBuffer The cutoff plus buffer distance.
     */
    public DomainDecompositionOctree(int nAtoms, Crystal crystal, double cutoffPlusBuffer, int dimension) {
//      this.nEdge = 2;
      this.nEdge = dimension;
//      this.nSearch = 2 * nEdge + 1;
      this.nSearch = dimension;
      targetInterfacialRadius = cutoffPlusBuffer / (double) nEdge;
      logger.info(format("** nEdge %d and nSearch %d", nEdge, nSearch));
      initDomainDecomposition(nAtoms, crystal);
    }

    /**
     * Initialize the domain decomposition.
     */
    public void initDomainDecomposition(int nAtoms, Crystal crystal) {
      this.nAtoms = nAtoms;
      this.crystal = crystal;

      // Allocate memory for fractional coordinates and subcell pointers for each atom.
      if (cellA == null || cellA.length < nAtoms) {
        cellA = new int[nAtoms];
        cellB = new int[nAtoms];
        cellC = new int[nAtoms];
      }

      // Determine the number of subcells along each axis.
      int maxA = nSearch;
      int maxB = nSearch;
      int maxC = nSearch;
      logger.info(format("** Max A %d Max B %d Max C %d\n",maxA, maxB, maxC));
      logger.info(format("** Target interfacial radius %4.4f\n",targetInterfacialRadius));
      logger.info(format("** Crystal interfacial radii %4.4f %4.4f %4.4f\n",crystal.interfacialRadiusA, crystal.interfacialRadiusB, crystal.interfacialRadiusC));

      /*
        To satisfy the cutoff + buffer distance, the number of subcells (maxA, maxB, maxC)
        is guaranteed to be at least nSearch - 1 for periodic crystals.
        */
//      if (!crystal.aperiodic()) {
//        assert (maxA >= nSearch - 1);
//        assert (maxB >= nSearch - 1);
//        assert (maxC >= nSearch - 1);
//      }

      /*
        If the number of subcells less than nSearch, then turn off the overhead of the
        subcell search by setting the number of cells along that axis to 1.
       */
      if (maxA < nSearch) nA = 1;
      else nA = maxA;
      nB = maxB < nSearch ? 1 : maxB;
      nC = maxC < nSearch ? 1 : maxC;

//      logger.info(format("** number of subcells %d %d %d",nA,nB,nC));

      // Allocate memory for the subcells.
      if (cells == null || cells.length != nA || cells[0].length != nB || cells[0][0].length != nC) {
        cells = new Cell[nA][nB][nC];
        for (int i = 0; i < nA; i++) {
          for (int j = 0; j < nB; j++) {
            for (int k = 0; k < nC; k++) {
              cells[i][j][k] = new Cell(i, j, k, targetInterfacialRadius);
            }
          }
        }
      } else {
        clear();
      }

    }

    /**
     * Clear the sub-cells of the atom list.
     */
    public void clear() {
      // Clear cell contents.
      for (int i = 0; i < nA; i++) {
        for (int j = 0; j < nB; j++) {
          for (int k = 0; k < nC; k++) {
            cells[i][j][k].clear();
          }
        }
      }

      // Clear special position exclusions.
      // if (specialPositionExclusionMap != null) {
      // Clear the lists, but don't remove the keys.
      //   for (List<Integer> list : specialPositionExclusionMap.values()) {
      //     list.clear();
      //   }
      //   specialPositionExclusionMap.clear();
      // }
    }

    /**
     * Add a special position exclusion.
     *
     * @param molecule The molecule.
     * @param symOpIndex The symmetry operation index.
     */
    public void addSpecialPositionExclusion(int molecule, int symOpIndex) {
      // Initialize the map.
      if (specialPositionExclusionMap == null) {
        specialPositionExclusionMap = new HashMap<>();
      }
      // Initialize the list for the molecule.
      List<Integer> list = specialPositionExclusionMap.get(molecule);
      if (list == null) {
        list = new ArrayList<>();
        specialPositionExclusionMap.put(molecule, list);
      }
      // Add the symOpIndex to the list.
      if (!list.contains(symOpIndex)) {
        list.add(symOpIndex);
      }
    }

    /**
     * Check if a special position exclusion exists.
     *
     * @param molecule The molecule.
     * @param symOpIndex The symmetry operation index.
     * @return True if the special position exclusion exists.
     */
    public boolean isSpecialPositionExclusion(int molecule, int symOpIndex) {
      if (specialPositionExclusionMap == null) {
        return false;
      }
      List<Integer> list = specialPositionExclusionMap.get(molecule);
      if (list == null) {
        return false;
      }
      return list.contains(symOpIndex);
    }

    public Cell[][][] getCells() {
      return cells;
    }

    public int[] getNumDivisions() {
      return new int[]{nA, nB, nC};
    }

    /**
     * Get the special position SymOps for a molecule.
     *
     * @param molecule The molecule.
     * @return The special position SymOps.
     */
    public List<Integer> getSpecialPositionSymOps(int molecule) {
      if (specialPositionExclusionMap == null) {
        return null;
      }
      if (specialPositionExclusionMap.containsKey(molecule)) {
        List<Integer> list = specialPositionExclusionMap.get(molecule);
        if (list.isEmpty()) {
          return null;
        } else {
          return list;
        }
      }
      return null;
    }

    /**
     * Log information about the domain decomposition.
     */
    public void log() {
      int nCells = nA * nB * nC;
      int nSymm = crystal.spaceGroup.getNumberOfSymOps();
      StringBuilder sb = new StringBuilder("  Neighbor List Builder\n");
      sb.append(format("   Total Cells:          %8d\n", nCells));
      if (nCells > 1) {
        sb.append(format("   Interfacial radius    %8.3f A\n", targetInterfacialRadius));
        int searchA;
        if (nA == 1) searchA = 1;
        else searchA = nSearch;
        int searchB = nB == 1 ? 1 : nSearch;
        int searchC = nC == 1 ? 1 : nSearch;
        sb.append(format("   Domain Decomposition: %8d %4d %4d\n", nA, nB, nC));
        sb.append(format("   Neighbor Search:      %8d x%3d x%3d\n", searchA, searchB, searchC));
        sb.append(format("   Total Atoms in System %8d\n", nAtoms));
        sb.append(format("   Mean Atoms per Cell:  %8d\n", nAtoms * nSymm / nCells));
      }
      int sum = 0;
      for (int i = 0; i < nA; i++){
        for (int j = 0; j < nB; j++){
          for (int k = 0; k < nC; k++){
            Cell cell = cells[i][j][k];
            sb.append(format("nLL Number of atoms in cell %d %d %d : %d\n",i,j,k,cell.getCount()));
            sum += cell.getCount();
          }
        }
      }
      sb.append(format("Total number of atoms in cells %d\n", sum));
      logger.info(sb.toString());
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
    public Cell image(int i, int j, int k) {
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

    /**
     * Add an atom to a sub-cell.
     *
     * @param i The index of the atom.
     * @param iSymm The index of the symmetry operator.
     * @param frac The fractional coordinates of the atom.
     */
    public void addAtomToCell(int i, int iSymm, double[] frac) {
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
        // Set the sub-cell indices for an asymmetric unit atom.
        cellA[i] = a;
        cellB[i] = b;
        cellC[i] = c;
      }
      cells[a][b][c].add(i, iSymm);
    }

    /**
     * Get the sub-cell for an asymmetric unit atom.
     *
     * @param i The index of the atom.
     * @return The sub-cell for the atom.
     */
    public Cell getCellForAtom(int i) {
      return cells[cellA[i]][cellB[i]][cellC[i]];
    }
  }

  /**
   * Hold the atom index and its symmetry operator.
   */
  public static class AtomIndex {

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
  public static class Cell {

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

    final double sideLength;

    /**
     * The list of atoms in the cell, together with their symmetry operator.
     */
    final List<AtomIndex> list;

    public Cell(int a, int b, int c, double sideLength) {
      this.a = a;
      this.b = b;
      this.c = c;
      this.sideLength = sideLength;
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

    /**
     * Returns the AtomIndex list of atoms in the cell
     *
     * @return AtomIndex list
     */
    public List<AtomIndex> getAtomIndexList() {
      return list;
    }

    /**
     * Return integer array of atom indices
     * @return int indices
     */
    public int[] getAtomIndices(){
      int[] indices = new int[list.size()];
      for (int i = 0; i < list.size(); i++){
        indices[i] = list.get(i).i;
      }
      return indices;
    }

    public int getCount() {
      return list.size();
    }

    /**
     * Returns the cell indices along the a, b, and c axes
     *
     * @return a, b, and c indices
     */
    public int[] getIndices() {
      return new int[]{a, b, c};
    }

    public double getSideLength() {return sideLength;}

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
            throw new IndexOutOfBoundsException(
                "Index out of bounds: count " + count + " " + index.length + " " + list.size());
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
