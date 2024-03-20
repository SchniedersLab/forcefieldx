//******************************************************************************
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
//******************************************************************************

package ffx.algorithms.optimize;

import com.google.common.collect.Lists;
import com.google.common.collect.MinMaxPriorityQueue;
import ffx.numerics.math.HilbertCurveTransforms;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Molecule;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.FastMath;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static ffx.potential.utils.Superpose.applyRotation;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * TorsionSearch class for performing a torsion scan on a molecule in a molecular assembly. It
 * provides options to run a static comparison of each bond on it's own, or run the full
 * exponential number of configurations. Options to run in parallel with multiple workers
 * exists. It is also set up so that multiple workers can run on different molecules if runWorker
 * is not used.
 *
 * @author Matthew J. Speranza
 * @author Arron J. Nessler
 */
public class TorsionSearch {
  private static final Logger logger = Logger.getLogger(TorsionSearch.class.getName());

  private final MolecularAssembly molecularAssembly;
  private final ForceFieldEnergy forceFieldEnergy;
  private final double[] x;
  private final boolean run;

  private final List<Bond> torsionalBonds;
  private final List<Atom[]> atomGroups;

  private final int nTorsionsPerBond;
  private final int nBits;
  private int nTorsionalBonds;
  private final int returnedStates;
  private long numConfigs;
  private long numIndices;
  private long end;
  private double minEnergy = Double.MAX_VALUE;
  private final List<Double> energies = new ArrayList<>();
  private final List<Long> hilbertIndices = new ArrayList<>();
  private final List<AssemblyState> states = new ArrayList<>();
  private MinMaxPriorityQueue<StateContainer> queue = null;

  private long[] workerAssignments;
  private long indicesPerAssignment;
  private boolean runWorker = false;

  public TorsionSearch(MolecularAssembly molecularAssembly,
                       Molecule molecule,
                       int nTorsionsPerBond,
                       int returnedStates
  ) {
    this.molecularAssembly = molecularAssembly;
    if (ArrayUtils.contains(molecularAssembly.getMoleculeArray(), molecule)) {
      run = true;
    } else {
      logger.warning("Molecule is not part of the assembly. Torsion scan will not run.");
      run = false;
    }
    this.forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    this.x = new double[forceFieldEnergy.getNumberOfVariables()];
    this.nTorsionsPerBond = nTorsionsPerBond; // Powers of two recommended
    this.nBits = (int) FastMath.ceil(FastMath.log(2, nTorsionsPerBond));
    this.torsionalBonds = getTorsionalBonds(molecule);
    this.atomGroups = getRotationGroups(torsionalBonds);
    this.nTorsionalBonds = torsionalBonds.size();
    this.numConfigs = (long) FastMath.pow(this.nTorsionsPerBond, torsionalBonds.size());
    this.numIndices = (long) FastMath.pow(2, this.nBits * torsionalBonds.size());
    if (returnedStates != -1) {
      this.returnedStates = returnedStates;
    } else {
      this.returnedStates = (int) numConfigs;
    }
    this.end = numIndices;
  }

  /**
   * Scanning through nTorsionsPerBond^nTorsionsPerBond (nTorsionsPerBond x nTorsionsPerBond x nTorsionsPerBond x ...)
   * array. Uses a hilbert curve to get to each point in the discrete space with minimal rotations
   * between each state without recursion. Hilbert curve implementation based on OpenMM's
   * c++ implementation. Runs through ALL indices with this worker. Use other methods to run just
   * one workers indices assigned by the buildWorker() method.
   */
  public void spinTorsions() {
    logger.info("\n ----------------- Starting torsion scan -----------------");
    logger.info(format(" Number of configurations: %d", numConfigs));
    logger.info(format(" Number of indices: %d", numIndices));
    this.spinTorsions(0, this.end);
  }

  /**
   * Similar to above, but with a limited number of hilbert indices to run through with
   * this worker alone.
   */
  public void spinTorsions(long start, long end) {
    logger.info("\n ----------------- Starting torsion scan -----------------");
    logger.info(format(" Number of configurations: %d", numConfigs));
    logger.info(format(" Number of indices: %d", numIndices));
    this.spinTorsions(start, end, true);
  }

  /**
   * Similar to above, but with a limited number of hilbert indices to run through with
   * this worker alone. Can also choose whether to update the lists of energies, indices, and states.
   * Private to prevent user from running with false, since this will not update the lists, which
   * will waste time & resources if the lists are never updated.
   */
  private void spinTorsions(long start, long end, boolean updateLists) {
    if (!run) {
      logger.warning("Torsion spin returning early since molecule not part of molecular assembly.");
      return;
    }

    // Initialize the coordinates to all zeros --> torsions depend on initial configuration
    long[] currentState = new long[nTorsionalBonds];

    // Create a priority queue to store the lowest energy states
    if (queue == null) {
      queue = MinMaxPriorityQueue
          .maximumSize(returnedStates)
          .create();
    }

    int progress = 0; // Worker jobs and specified start/end jobs cannot finish early -> can't know how many will exist
    while (start <= end && progress < numConfigs) {
      long[] newState = HilbertCurveTransforms.hilbertIndexToCoordinates(nTorsionalBonds, nBits, start);
      if (ArrayUtils.contains(newState, nTorsionsPerBond)) {
        start++;
        continue;
      }
      //logger.info(format(" Hilbert Index: %d; Coordinate State: " + Arrays.toString(newState), start));
      // Permute from currentState to newState
      changeState(currentState, newState, nTorsionsPerBond, torsionalBonds, atomGroups);
      // Update coordinates
      forceFieldEnergy.getCoordinates(x);
      // Calculate energy
      double energy = forceFieldEnergy.energy(x);
      // Log minimum energy
      if (energy < minEnergy) {
        minEnergy = energy;
        logger.info(format("\n New minimum energy: %12.5f", minEnergy));
        logger.info(format(" Hilbert Index: %d; Coordinate State: " + Arrays.toString(newState), start));
      }
      // Add to queue
      queue.add(new StateContainer(new AssemblyState(molecularAssembly), energy, start));
      currentState = newState;
      start++;
      progress++;
    }
    // Revert back to original state
    changeState(currentState, new long[nTorsionalBonds], nTorsionsPerBond, torsionalBonds, atomGroups);
    if (progress == numConfigs) {
      logger.info("\n Completed all configurations before end index.");
    }
    if (updateLists) {
      this.updateInfoLists();
    }
  }

  /**
   * Static analysis of torsional bonds. Optionally remove bonds that cause large energy changes
   * outside of a specified threshold. Only one worker should run this method on the same structure.
   *
   * @param numRemove            Number of bonds to remove if they exceed the threshold.
   * @param eliminationThreshold Threshold for removing bonds and logging as high energy.
   */
  public void staticAnalysis(int numRemove, double eliminationThreshold) {
    if (!run) {
      logger.warning("Static analysis returning early since molecule not part of molecular assembly.");
      return;
    }
    if (queue == null) {
      queue = MinMaxPriorityQueue
          .maximumSize(returnedStates)
          .create();
    }
    eliminationThreshold = FastMath.abs(eliminationThreshold);
    // Update coordinates
    forceFieldEnergy.getCoordinates(x);
    // Calculate energy as a reference
    double initialE = forceFieldEnergy.energy(x);
    ArrayList<Integer> remove = new ArrayList<>();
    long[] state = new long[nTorsionalBonds];
    long[] oldState = new long[nTorsionalBonds];
    // Loop through all torsional bonds up to the number of torsions per bond finding energies and resetting
    for (int i = 0; i < nTorsionalBonds; i++) {
      for (int j = i == 0 ? 0 : 1; j < nTorsionsPerBond; j++) {
        state[i] = j;
        changeState(oldState, state, nTorsionsPerBond, torsionalBonds, atomGroups);
        // Update coordinates
        forceFieldEnergy.getCoordinates(x);
        // Calculate energy
        double newEnergy = forceFieldEnergy.energy(x);
        if (newEnergy - initialE > eliminationThreshold && !remove.contains(i)) {
          remove.add(i);
        } else {
          long index = HilbertCurveTransforms.coordinatesToHilbertIndex(nTorsionalBonds, nBits, state);
          queue.add(new StateContainer(new AssemblyState(molecularAssembly), newEnergy, index));
        }
        changeState(state, oldState, nTorsionsPerBond, torsionalBonds, atomGroups);
      }
      state[i] = 0;
    }
    remove.sort(Collections.reverseOrder());
    logger.info("\n " + remove.size() + " bonds that cause large energy increase: " + remove);
    if (remove.size() > numRemove) {
      remove = Lists.newArrayList(remove.subList(0, numRemove));
    }
    for (int i = remove.size() - 1; i >= 0; i--) {
      logger.info(" Removing bond: " + torsionalBonds.get(remove.get(i)));
      logger.info(" Bond index: " + remove.get(i));
      torsionalBonds.set(remove.get(i), null);
      atomGroups.set(remove.get(i), null);
    }
    torsionalBonds.removeAll(Collections.singleton(null));
    atomGroups.removeAll(Collections.singleton(null));
    this.nTorsionalBonds = torsionalBonds.size();
    // Update parameters based on new set of bonds
    this.end = (long) FastMath.pow(2, nBits * nTorsionalBonds);
    this.numConfigs = (long) FastMath.pow(nTorsionsPerBond, nTorsionalBonds);
    this.numIndices = this.end;
    logger.info(" Finished static analysis.");
    logger.info(format(" Number of configurations after elimination: %d", numConfigs));
    logger.info(format(" Number of indices after elimination: %d", numIndices));
    this.updateInfoLists();
  }

  /**
   * Builds the worker assignments for each rank. Should be run before spinTorsions.
   *
   * @param rank
   * @param worldSize
   * @return whether the worker was built successfully
   */
  public boolean buildWorker(int rank, int worldSize) {
    if (!run) {
      logger.warning("Build worker returning early since molecule not part of molecular assembly.");
      return false;
    }
    if (rank >= worldSize) {
      logger.warning(" Rank is greater than world size.");
      return false;
    }
    runWorker = true;
    // Assign each worker the same number of indices
    this.workerAssignments = new long[worldSize];
    long jobsPerWorker = numIndices / worldSize;
    this.indicesPerAssignment = jobsPerWorker / worldSize;
    logger.info("\n Jobs per worker: " + jobsPerWorker);
    logger.info(" Jobs per worker split: " + indicesPerAssignment);
    // Assign starting indexes to each worker
    for (int i = 0; i < worldSize; i++) {
      workerAssignments[i] = i * jobsPerWorker + rank * indicesPerAssignment;
    }
    logger.info(" Worker " + rank + " assigned indices: " + Arrays.toString(workerAssignments));
    return true;
  }

  /**
   * Runs this worker with the indices assigned to it by the buildWorker() method.
   */
  public void runWorker() {
    if (!run) {
      logger.warning("Worker returning early since molecule not part of molecular assembly.");
      return;
    } else if (!runWorker) {
      logger.warning("Worker returning early since worker not built or is invalid.");
      return;
    }
    logger.info("\n ----------------- Starting torsion scan -----------------");
    logger.info(format("\n Number of configurations before worker starts: %d", numConfigs));
    logger.info(format(" Number of indices before worker starts: %d", numIndices));
    for (int i = 0; i < workerAssignments.length; i++) {
      logger.info(format(" Worker torsion assignment %3d of %3d: %12d to %12d.", i + 1, workerAssignments.length,
          workerAssignments[i], workerAssignments[i] + indicesPerAssignment - 1));
      // Only update lists on last assignment
      this.spinTorsions(workerAssignments[i], workerAssignments[i] + indicesPerAssignment - 1,
          i == workerAssignments.length - 1);
    }
  }

  /**
   * List of energies for each state in order of lowest to highest energy.
   */
  public List<Double> getEnergies() {
    return energies;
  }

  /**
   * List of hilbert indices for each state in order of lowest to highest energy.
   */
  public List<Long> getHilbertIndices() {
    return hilbertIndices;
  }

  /**
   * List of states in order of lowest to highest energy.
   */
  public List<AssemblyState> getStates() {
    return states;
  }

  /**
   * @return the end index of the system
   */
  public long getEnd() {
    return end;
  }

  /**
   * Updates state, energy, and hilbert index lists after a run of one of the public methods.
   */
  private void updateInfoLists() {
    if (!states.isEmpty()) {
      // Refill limited-size priority queue with previous information
      for (int i = 0; i < states.size(); i++) {
        queue.add(new StateContainer(states.get(i), energies.get(i), hilbertIndices.get(i)));
      }
      // Clear lists
      states.clear();
      energies.clear();
      hilbertIndices.clear();
    }
    // Refill lists
    while (!queue.isEmpty()) {
      StateContainer toBeSaved = queue.removeFirst();
      states.add(toBeSaved.getState());
      energies.add(toBeSaved.getEnergy());
      hilbertIndices.add(toBeSaved.getIndex());
    }
  }

  /**
   * Finds torisonal bonds in a molecule ignoring methyl groups, hydrogens, and up to 6-membered rings.
   *
   * @param molecule molecule of interest
   * @return list of torsional bonds
   */
  private static List<Bond> getTorsionalBonds(Molecule molecule) {
    List<Bond> torsionalBonds = new ArrayList<>();
    for (Bond bond : molecule.getBondList()) {
      Atom a1 = bond.getAtom(0);
      Atom a2 = bond.getAtom(1);
      List<Bond> bond1 = a1.getBonds();
      int b1 = bond1.size();
      List<Bond> bond2 = a2.getBonds();
      int b2 = bond2.size();
      // Ignore single bonds (ex. alcohol groups)
      if (b1 > 1 && b2 > 1) {
        // Ignore methyl groups
        if (a1.getAtomicNumber() == 6 && a1.getNumberOfBondedHydrogen() == 3 || a2.getAtomicNumber() == 6 && a2.getNumberOfBondedHydrogen() == 3) {
          continue;
        }
        // Ignore ring atoms
        if (a1.isRing(a2)) {
          continue;
        }
        torsionalBonds.add(bond);
        logger.info(" Bond " + bond + " is a torsional bond.");
      }
    }
    return torsionalBonds;
  }

  /**
   * Get the smallest group of atoms to rotate about a bond.
   *
   * @param bonds Torsional bonds
   * @return List of atoms to rotate about each bond in the given list
   */
  private static List<Atom[]> getRotationGroups(List<Bond> bonds) {
    List<Atom[]> rotationGroups = new ArrayList<>();
    for (Bond bond : bonds) {
      Atom a1 = bond.getAtom(0);
      Atom a2 = bond.getAtom(1);
      // We should have two atoms with a spinnable bond
      List<Atom> a1List = new ArrayList<>();
      List<Atom> a2List = new ArrayList<>();
      // Search for rotation groups
      searchTorsions(a1, a1List, a2);
      searchTorsions(a2, a2List, a1);
      // Convert List to array
      Atom[] a1Array = new Atom[a1List.size()];
      Atom[] a2Array = new Atom[a2List.size()];
      a1List.toArray(a1Array);
      a2List.toArray(a2Array);
      // Add smaller list to rotation groups
      if (a1List.size() > a2List.size()) {
        rotationGroups.add(a2Array);
      } else {
        rotationGroups.add(a1Array);
      }
    }
    return rotationGroups;
  }

  /**
   * Rotate an atom about another.
   *
   * @param u     a unit vector to rotate {@link ffx.potential.bonded.Atom} object about.
   * @param a2    a {@link ffx.potential.bonded.Atom} object to create axis of rotation.
   * @param theta Amount to rotate by in degrees.
   */
  private static void rotateAbout(double[] u, Atom a2, double theta) {

    theta = toRadians(theta);

    // Define quaternion from axis-angle
    double[] quaternion = new double[4];
    quaternion[0] = cos(theta / 2);
    quaternion[1] = u[0] * sin(theta / 2);
    quaternion[2] = u[1] * sin(theta / 2);
    quaternion[3] = u[2] * sin(theta / 2);
    // Normalize quaternion
    double quaternionNorm = 1 / Math.sqrt(quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1]
        + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3]);
    for (int i = 0; i < 4; i++) {
      quaternion[i] *= quaternionNorm;
    }
    // Useful storage
    double q1q1 = quaternion[1] * quaternion[1];
    double q2q2 = quaternion[2] * quaternion[2];
    double q3q3 = quaternion[3] * quaternion[3];
    double q0q1 = quaternion[0] * quaternion[1];
    double q0q2 = quaternion[0] * quaternion[2];
    double q0q3 = quaternion[0] * quaternion[3];
    double q1q2 = quaternion[1] * quaternion[2];
    double q1q3 = quaternion[1] * quaternion[3];
    double q2q3 = quaternion[2] * quaternion[3];
    // Quaternion rotation matrix
    double[][] rotation2 = new double[3][3];
    rotation2[0][0] = 1 - 2 * (q2q2 + q3q3);
    rotation2[0][1] = 2 * (q1q2 - q0q3);
    rotation2[0][2] = 2 * (q1q3 + q0q2);
    rotation2[1][0] = 2 * (q1q2 + q0q3);
    rotation2[1][1] = 1 - 2 * (q1q1 + q3q3);
    rotation2[1][2] = 2 * (q2q3 - q0q1);
    rotation2[2][0] = 2 * (q1q3 - q0q2);
    rotation2[2][1] = 2 * (q2q3 + q0q1);
    rotation2[2][2] = 1 - 2 * (q1q1 + q2q2);

    double[] a2XYZ = new double[3];
    a2.getXYZ(a2XYZ);
    applyRotation(a2XYZ, rotation2);
    a2.setXYZ(a2XYZ);
  }

  /**
   * Identify atoms that should be rotated.
   *
   * @param seed    an {@link ffx.potential.bonded.Atom} object to rotate about.
   * @param atoms   a list of {@link ffx.potential.bonded.Atom} objects to rotate.
   * @param notAtom Avoid this atom (wrong side of bond).
   */
  private static void searchTorsions(@Nullable Atom seed, List<Atom> atoms, Atom notAtom) {
    if (seed == null) {
      return;
    }
    atoms.add(seed);
    for (Bond b : seed.getBonds()) {
      Atom nextAtom = b.get1_2(seed);
      if (nextAtom == notAtom || atoms.contains(nextAtom)) { //nextAtom.getParent() != null ||
        continue;
      }
      searchTorsions(nextAtom, atoms, notAtom);
    }
  }

  /**
   * Change the state of the molecule (rotations) to match newState based on changes from the oldState.
   *
   * @param oldState  current positions
   * @param newState  new positions
   * @param nTorsions number of torsions per bond
   * @param bonds     list of bonds that the states define
   * @param atoms     list of atoms to rotate with each bond
   */
  private static void changeState(long[] oldState, long[] newState, int nTorsions,
                                  List<Bond> bonds, List<Atom[]> atoms) {
    for (int i = 0; i < oldState.length; i++) {
      if (oldState[i] != newState[i]) {
        // Add change required to go from oldState to newState to change array
        int change = (int) (newState[i] - oldState[i]);
        // Apply change at this index to atoms
        int turnDegrees = (change * (360 / nTorsions));
        // Get vector from atom to atom of bond i
        double[] u = new double[3];
        double[] translation = new double[3];
        double[] a1 = bonds.get(i).getAtom(0).getXYZ(new double[3]);
        double[] a2 = bonds.get(i).getAtom(1).getXYZ(new double[3]);
        if (atoms.get(i)[0] == bonds.get(i).getAtom(0)) { // if a1 had the smaller rotation group
          for (int j = 0; j < 3; j++) {
            // Rotation about this vector is same as the same rotation about the opposite vector?
            u[j] = a1[j] - a2[j]; // Vector defined as line segment from a2 to a1
            translation[j] = -a1[j];
          }
          // Unit vector
          double unit = 1 / Math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
          for (int j = 0; j < 3; j++) {
            u[j] *= unit;
          }
          for (int j = 0; j < atoms.get(i).length; j++) {
            atoms.get(i)[j].move(translation);
            rotateAbout(u, atoms.get(i)[j], turnDegrees);
            atoms.get(i)[j].move(a1);
          }
        } else {
          for (int j = 0; j < 3; j++) {
            u[j] = a2[j] - a1[j]; // Vector defined as line segment from a1 to a2
            translation[j] = -a2[j];
          }
          double unit = 1 / Math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
          for (int j = 0; j < 3; j++) {
            u[j] *= unit;
          }
          for (int j = 0; j < atoms.get(i).length; j++) {
            atoms.get(i)[j].move(translation);
            rotateAbout(u, atoms.get(i)[j], turnDegrees);
            atoms.get(i)[j].move(a2); // This line different then above
          }
        }
      }
    }
  }

  /**
   * Implements StateContainer to store the coordinates of a state and its energy
   */
  private record StateContainer(AssemblyState state, double e,
                                long hilbertIndex) implements Comparable<StateContainer> {
    AssemblyState getState() {
      return state;
    }

    double getEnergy() {
      return e;
    }

    long getIndex() {
      return hilbertIndex;
    }

    @Override
    public int compareTo(StateContainer o) {
      return Double.compare(e, o.getEnergy());
    }
  }
}
