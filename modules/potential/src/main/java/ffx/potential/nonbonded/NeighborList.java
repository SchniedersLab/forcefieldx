// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.max;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

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
  /** Lists of groups */
  private ArrayList<int[]> groups;
  /** Span of each group across zCell indices (each at most length M ) */
  private ArrayList<HashSet<Integer>> zSets;
  /** Lists of interactions between groups of atoms. nGroups is not known until after groups are built. */
  private int[][][] groupLists;
  /** Size of each group of atoms */
  private int M = -1;
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

  public void buildMxNList(int M, final double[][] coordinates, final int[][][] lists, boolean[] use, boolean forceRebuild, boolean print){
    if (disableUpdates) {
      return;
    }
    this.coordinates = coordinates;
    this.lists = lists;
    this.M = M;
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
        // Create groups of M and allocate mem
        if (M != -1){
          groups();
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
   * Build groups of size M by walking up Z columns. Must be called
   * after domain decomposition & assignments of atoms into cells.
   * <p></p>
   * Allocates memory for groupLists and cell maps in domain decomposition
   * object. Builds groups and zSets array.
   */
  private void groups() {
    // Use cells to build groups
    int nA = domainDecomposition.nA;
    int nB = domainDecomposition.nB;
    int nC = domainDecomposition.nC;
    zSets = new ArrayList<>();
    for(int i = 0; i < nA; i++){
      for(int j = 0; j < nB; j++){
        HashSet<Integer> zSet = new HashSet<>();
        int[] group = new int[M];
        // Slow in-order add, but typically fast due to small # zCol atoms
        PriorityQueue<AtomIndex> atomQueue = new PriorityQueue<>();
        for(int k = 0; k < nC; k++){
          // Walk up z-columns adding to groups of size M
          Cell cell = domainDecomposition.getCell(i, j, k);
          assert(cell != null); // after domain decomposition
          atomQueue.addAll(cell.list);
          while(atomQueue.size() >= M){
            for(int l = 0; l < M; l++){
              AtomIndex atom = atomQueue.poll();
              assert(atom != null); // with while-loop cond.
              // Unique index for every atom
              group[l] = nSymm > nAtoms ? atom.iSymm*nSymm + atom.i : atom.i*nAtoms + atom.iSymm;
              double[] xyz = new double[]{ coordinates[atom.iSymm][atom.i + XX],
                      coordinates[atom.iSymm][atom.i + YY],
                      coordinates[atom.iSymm][atom.i + ZZ]};
              double[] frac = new double[3];
              crystal.toFractionalCoordinates(xyz, frac);
              // Keep track of groups across multiple z cells
              zSet.add((int) floor(frac[ZZ]*nC));
            }
            groups.add(group);
            group = new int[4];
            zSets.add(zSet);
            zSet.clear();
          }
        }
        if(atomQueue.isEmpty()){ continue; }
        // Fill in last group with -1 filler
        for(int k = 0; k < M; k++){
          if(!atomQueue.isEmpty()){
            AtomIndex atom = atomQueue.poll();
            assert(atom != null);
            group[k] = nSymm > nAtoms ? atom.iSymm*nSymm + atom.i : atom.i*nAtoms + atom.iSymm;
            double[] xyz = new double[]{ coordinates[atom.iSymm][atom.i + XX],
                    coordinates[atom.iSymm][atom.i + YY],
                    coordinates[atom.iSymm][atom.i + ZZ]};
            double[] frac = new double[3];
            crystal.toFractionalCoordinates(xyz, frac);
            zSet.add((int) floor(frac[ZZ]*nC));
          } else{
            group[k] = -1;
          }
          groups.add(group);
          zSets.add(zSet);
          zSet.clear();
        }
      }
    }
    // Old implementation that doesn't use pre-built cells
    /*
    final SortedMap<Double, Integer>[][] map = new SortedMap[nA][nB];
    int[] groupZ = new int[nAtoms*nSymm];
    for(int i = 0; i < nA; i++) {
      for(int j = 0; j < nB; j++) {
        map[i][j] = new TreeMap<>(); // Sorted map
        for(int m = 0; m < nSymm; m++) {
          final double[] xyz = coordinates[m];
          for (int k = 0; k < nAtoms; k++) { // For every atom (including nSymm atoms)
            final double[] pos = new double[3];
            final double[] posFrac = new double[3];
            pos[0] = xyz[k + XX];
            pos[1] = xyz[k + YY];
            pos[2] = xyz[k + ZZ];
            crystal.toFractionalCoordinates(pos, posFrac);
            // Shift into 0.0 <= x,y,z <= 1.0 fractional coords
            for(int ii = 0; ii < posFrac.length; ii++)
              while(posFrac[ii] < 0.0 || posFrac[ii] > 1.0)
                posFrac[ii] += posFrac[ii] < 0.0 ? 1.0 : -1.0;

            int xIndex = (int) floor(posFrac[0]*nA);
            int yIndex = (int) floor(posFrac[1]*nB);
            // Unique index for every atom including symm atoms
            int atomIndex = nSymm > nAtoms ? m*nSymm + k : k*nAtoms + m;
            groupZ[atomIndex] = (int) floor(posFrac[2]*nC);
            // log((n/(nA*nB)) put onto map complexity -> still very fast if uniform density
            map[xIndex][yIndex].put(posFrac[2], atomIndex);
          }
        }
      }
    }
    // Place most into groups of M (others < M with -1 as filler)
    groups = new ArrayList<>((int) ((double)(nSymm*nAtoms)/M*1.01)); // 1% room for extra groups
    zSets = new HashSet[(int) ((double)(nSymm*nAtoms)/M*1.01)]; // Dynamic reallocation if required
    for (int i = 0; i < nA; i++) {
      for (int j = 0; j < nB; j++) {
        if (map[i][j].isEmpty()) continue;
        // This is a list of atomIndex as defined 14 lines above
        List<Integer> columnAtoms = map[i][j].values().stream().toList();
        map[i][j].clear();
        int[] group = new int[M];
        int count = 0;
        while (columnAtoms.size() - count >= M) {
          // Keep span of group across z
          HashSet<Integer> zSet = new HashSet<>();
          for (int k = 0; k < M; k++) {
            group[k] = columnAtoms.get(count);
            zSet.add(groupZ[group[k]]);
            count++;
          }
          // Store & reset
          groups.add(group);
          if (groups.size() - 1 == zSets.length){ // Dynamic allocation
            HashSet<Integer>[] temp = new HashSet[(int) (zSets.length*1.25)];
            System.arraycopy(zSets, 0, temp, 0, zSets.length);
            zSets = temp;
          }
          else {
            zSets[groups.size() - 1] = zSet;
          }
          group = new int[M];
        }
        if (columnAtoms.size() - count == 0) {
          continue;
        }
        // Handle any left behind with -1 as filler
        HashSet<Integer> zSet = new HashSet<>();
        for (int k = 0; k < M; k++) {
          if (k <= (columnAtoms.size() - count)) {
            group[k] = columnAtoms.get(count);
            zSet.add(groupZ[group[k]]);
          } else {
            group[k] = -1;
          }
        }
        groups.add(group);
        // Place group into respective cells
        if (groups.size() - 1 == zSets.length){ // Dynamic allocation
          HashSet<Integer>[] temp = new HashSet[(int) (zSets.length*1.25)];
          System.arraycopy(zSets, 0, temp, 0, zSets.length);
          zSets = temp;
        }
        else {
          zSets[groups.size() - 1] = zSet;
        }
      }
    }
    groups.trimToSize();
     */
    groupLists = new int[nSymm][groups.size()][];
    domainDecomposition.allocateGroupCellMap(groups.size());
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
        }

        if(M != -1) {
          // TODO: Independent call for this = better thread saturation
          int bottom = groups.size() > lb ? lb : 0;
          int top = groups.size() > ub ? ub : 0;
          assert (bottom <= top);
          for (int i = bottom; i <= top; i++) {
            int i3 = i * 3;
            cart[0] = xyz[i3 + XX];
            cart[1] = xyz[i3 + YY];
            cart[2] = xyz[i3 + ZZ];
            crystal.toFractionalCoordinates(cart, frac);
            // This index will be used later
            int groupIndex = nSymm > groups.size() ? iSymm * nSymm + i : i * groups.size() + iSymm;
            domainDecomposition.addGroupToCell(groupIndex, frac);
          }
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
              // (a, b+1..b+nE, c) -> radial line in y-dir
              for (int bi = b1; bi <= bStop; bi++) {
                atomCellPairs(domainDecomposition.image(a, bi, c));
              }
              // (a, b-nE..b+nE, c+1..c+nE) -> semicircle in in z,y plane
              for (int bi = bStart; bi <= bStop; bi++) {
                for (int ci = c1; ci <= cStop; ci++) {
                  atomCellPairs(domainDecomposition.image(a, bi, ci));
                }
              }
              // (a+1..a+nE, b-nE..b+nE, c-nE..c+nE) -> hemisphere pushed in x-dir
              for (int bi = bStart; bi <= bStop; bi++) {
                for (int ci = cStart; ci <= cStop; ci++) {
                  for (int ai = a1; ai <= aStop; ai++) {
                    atomCellPairs(domainDecomposition.image(ai, bi, ci));
                  }
                }
              }
            } else { // Full circle because????
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

        // TODO: Separate call for this -> better thread saturation
        // Either way we need to do atom based to get masking rules
        if( M != -1) {
          int[][] groupList = groupLists[iSymm];
          int bottom = groups.size() > lb ? lb : 0;
          int top = groups.size() > ub ? ub : 0;
          assert (lb <= ub);
          int cubesA = nA == 1 ? 0 : domainDecomposition.nEdge;
          int cubesB = nB == 1 ? 0 : domainDecomposition.nEdge;
          int cubesC = nC == 1 ? 0 : domainDecomposition.nEdge;
          for (int i = bottom; i <= top; i++) {
            int groupIndex = nSymm > groups.size() ? iSymm * nSymm + i : i * groups.size() + iSymm;
            Cell groupCell = domainDecomposition.getCellForGroup(groupIndex);
            int a = groupCell.a;
            int b = groupCell.b;
            int c = groupCell.c;
            boolean[] assigned = new boolean[groups.size()];
            ArrayList<Integer> neighbors = new ArrayList<>();
            for (int j = -cubesA + a; j <= cubesA + a; j++) {
              for (int k = -cubesB + b; k <= cubesB + b; k++) {
                // Bump search in z direction by number of cells this group spans in z direction
                cubesC += zSets.get(i).size() == 1 ? 0 : zSets.get(i).size() - 1;
                cubesC = nC == 1 ? 0 : cubesC;
                for (int l = -cubesC + c; l <= cubesC + c; l++) {
                  Cell target = domainDecomposition.image(j, k, l);
                  for (Integer group : target.groupList) {
                    if (!assigned[group]) { // No duplicates due to z cell uncertainty
                      assigned[group] = true;
                      neighbors.add(group);
                    }
                  }
                }
                cubesC -= zSets.get(i).size() == 1 ? 0 : zSets.get(i).size() - 1;
              }
            }
            neighbors.trimToSize();
            groupList[i] = new int[neighbors.size()];
            for (int j = 0; j < neighbors.size(); j++) {
              groupList[j][j] = neighbors.get(j);
            }
          }
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
    /** The cell indices of each atom along the A-axis. */
    private int[] cellA;
    /** The cell indices of each atom along the B-axis. */
    private int[] cellB;
    /** The cell indices of each atom along the C-axis. */
    private int[] cellC;
    private int[] groupCellA;
    private int[] groupCellB;
    private int[] groupCellC;
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
              cells[i][j][k] = new Cell(i, j, k);
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
      cells[a][b][c].add(i, iSymm, zu);
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

    public Cell getCell(int a, int b, int c) {
      if (a < 0 || a >= nA || b < 0 || b >= nB || c < 0 || c >= nC) {
        return null;
      }
      return cells[a][b][c];
    }

    public void addGroupToCell(int groupID, double[] frac) {
      for(int j = 0; j < frac.length; j++)
        while(frac[j] < 0.0 || frac[j] > 1.0)
          frac[j] = frac[j] < 0.0 ? frac[j]+1.0 : frac[j]-1.0;
      final int a = (int) floor(frac[XX] * nA);
      final int b = (int) floor(frac[YY] * nB);
      final int c = (int) floor(frac[ZZ] * nC);
      groupCellA[groupID] = a;
      groupCellB[groupID] = b;
      groupCellC[groupID] = c;
      cells[a][b][c].groupAdd(groupID);
    }

    public void allocateGroupCellMap(int size) {
      this.groupCellA = new int[crystal.getNumSymOps()*size];
      this.groupCellB = new int[crystal.getNumSymOps()*size];
      this.groupCellC = new int[crystal.getNumSymOps()*size];
    }

    public Cell getCellForGroup(int i) {
      return cells[groupCellA[i]][groupCellB[i]][groupCellC[i]];
    }
  }

  /**
   * Hold the atom index and its symmetry operator.
   */
  public static class AtomIndex implements Comparator<Double> {

    public final int iSymm;
    public final int i;
    public final double z;

    public AtomIndex(int iSymm, int i, double z) {
      this.iSymm = iSymm;
      this.i = i;
      this.z = z;
    }

    @Override
    public int compare(Double o1, Double o2) {
      return Double.compare(o1, o2);
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

    /**
     * The list of atoms in the cell, together with their symmetry operator.
     */
    final List<AtomIndex> list;
    /**
     * List of groups of atoms (only contains indices) that contain at least one atom in this cell.
     */
    final List<Integer> groupList;

    public Cell(int a, int b, int c) {
      this.a = a;
      this.b = b;
      this.c = c;
      list = Collections.synchronizedList(new ArrayList<>());
      groupList = Collections.synchronizedList(new ArrayList<>());
    }

    /**
     * Add an atom to the cell.
     *
     * @param atomIndex The atom index.
     * @param symOpIndex The symmetry operator index.
     */
    public void add(int atomIndex, int symOpIndex, double z) {
      list.add(new AtomIndex(symOpIndex, atomIndex, z));
    }

    public void groupAdd(int groupIndex){
      groupList.add(groupIndex);
    }

    public AtomIndex get(int index) {
      if (index >= list.size()) {
        return null;
      }
      return list.get(index);
    }

    public int getGroup(int groupIndex){
      if(groupIndex >= groupList.size()){
        return -1;
      }
      return groupList.get(groupIndex);
    }

    public int getCount() {
      return list.size();
    }

    public int getGroupCount(){
      return groupList.size();
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
